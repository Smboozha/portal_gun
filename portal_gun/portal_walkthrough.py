#!/usr/bin/env python3
"""
PORTAL WALK-THROUGH — Human-scale dual-sided portal controller (V1.0)

Два устройства A и B (человеческий масштаб), каждое с:
  - камерой (захват сцены перед порталом)
  - проектором/экраном (показ сцены с противоположной стороны)
  - кольцевым плазмотроном (зелёное свечение по краю)

Физика прохода:
  Оптическое реле (линза/зеркало) или видео-досветка обеспечивают
  непрерывную картинку, в которую можно ФИЗИЧЕСКИ шагнуть через проём.
  Зелёное кольцо реагирует на движение и «проход» человека.

Тестирование:
  python3 portal_walkthrough.py --test    — синтетический loopback, проверка связки
  python3 portal_walkthrough.py --role A  —.hardware node A
  python3 portal_walkthrough.py --role B --peer 10.0.0.1  — hardware node B

Зависимости: numpy, (cv2 опционально для JPEG), RPi.GPIO опционально.
"""

import argparse, asyncio, json, logging, os, socket, struct, time, sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

try:
    import cv2
    HAVE_CV2 = True
except ImportError:
    HAVE_CV2 = False
    logging.getLogger('PORTAL').warning("OpenCV не найден — raw-режим (без JPEG-сжатия)")

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s %(levelname)s] %(name)s: %(message)s')
log = logging.getLogger('PORTAL')

# ── CONFIG (human-scale) ─────────────────────────────────────────────
DEFAULT = {
    "camera": {"device": 0, "width": 1920, "height": 1080, "fps": 60},
    "projector": {"refresh": 60},
    "plasma": {"gpio_pin": 18, "power_min": 0.0, "power_max": 0.35},
    "network": {
        "role": "A", "peer_host": "10.0.0.2", "peer_port": 9000,
        "listen_port": 9000, "bind": "0.0.0.0",
        "codec": "JPEG", "quality": 80, "latency_budget_ms": 80,
    },
    "motion": {"sensitivity": 0.15},
    "calibration": {"checker_w": 8, "checker_h": 6, "square_px": 60},
}

# ── CODEC ─────────────────────────────────────────────────────────────
def encode(frame: np.ndarray, q: int) -> bytes:
    if HAVE_CV2:
        ok, buf = cv2.imencode('.jpg', frame, [int(cv2.IMWRITE_JPEG_QUALITY), q])
        return buf.tobytes() if ok else frame.tobytes()
    return frame.tobytes()

def decode(data: bytes, w: int, h: int, ch=3) -> np.ndarray:
    if HAVE_CV2:
        arr = np.frombuffer(data, np.uint8)
        img = cv2.imdecode(arr, cv2.IMREAD_COLOR)
        if img is not None and img.shape[1] == w and img.shape[0] == h:
            return img
    expected = w * h * ch
    if len(data) == expected:
        return np.frombuffer(data, np.uint8).reshape(h, w, ch)
    return np.zeros((h, w, ch), np.uint8)

# ── GREEN PORTAL OVERLAY ──────────────────────────────────────────────
def make_overlay(w: int, h: int, r0: float = 0.38) -> np.ndarray:
    ov = np.zeros((h, w, 3), np.uint8)
    cx, cy = w // 2, h // 2
    R = min(w, h) * r0
    yy, xx = np.mgrid[0:h, 0:w]
    dist = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    # свечение (12 колец)
    for rr in np.linspace(R * 0.98, R * 1.28, 12):
        band = np.abs(dist - rr) < 2.5
        ov[band] = np.maximum(ov[band], np.array([0, 255, 100], np.uint8))
    # основное кольцо
    ov[np.abs(dist - R) < 3.5] = np.array([0, 255, 150], np.uint8)
    ov[np.abs(dist - R) < 1.0] = np.array([140, 255, 210], np.uint8)
    # «молнии» (дуги) — только на пикселях кольца, по углу
    ring = np.abs(dist - R) < 3.0
    ang = np.arctan2((yy - cy)[ring], (xx - cx)[ring])
    r_idx, c_idx = np.nonzero(ring)
    rng = np.random.default_rng(int(time.time() % 1000))
    for _ in range(6):
        a = rng.uniform(0, 2 * np.pi)
        a0, a1 = a - 0.12, a + rng.uniform(0.1, 0.3)
        if a0 > a1:
            a0, a1 = a1, a0
        sel = (ang >= a0) & (ang <= a1)
        if not np.any(sel):
            continue
        ov[r_idx[sel], c_idx[sel]] = np.array([120, 255, 200], np.uint8)
    return ov

def overlay(frame: np.ndarray, ov: np.ndarray, a: float = 0.6) -> np.ndarray:
    mask = np.any(ov > 0, axis=2)
    out = frame.copy()
    glow = (ov.astype(np.float32) * a).astype(np.uint8)
    out[mask] = np.clip(frame[mask].astype(np.int16) + glow[mask].astype(np.int16),
                        0, 255).astype(np.uint8)
    return out

# ── WALK DETECTION ────────────────────────────────────────────────────
class WalkDetector:
    """Simple center-region motion detector → 'someone is walking through'."""
    def __init__(self, w: int, h: int, threshold: float = 0.03):
        self.cx, self.cy = w // 2, h // 2
        self.r = min(w, h) // 5
        self.thresh = threshold
        self._prev = None
        self.walk_count = 0

    def detect(self, frame: np.ndarray) -> bool:
        gray = frame.mean(axis=2)
        y0, y1 = max(0, self.cy - self.r), min(gray.shape[0], self.cy + self.r)
        x0, x1 = max(0, self.cx - self.r), min(gray.shape[1], self.cx + self.r)
        roi = gray[y0:y1, x0:x1]
        if self._prev is None:
            self._prev = roi.copy()
            return False
        diff = np.abs(roi.astype(np.float32) - self._prev.astype(np.float32)).mean()
        self._prev = 0.8 * self._prev + 0.2 * roi
        if diff > self.thresh * 255:
            self.walk_count += 1
            if self.walk_count > 3:
                self.walk_count = 0
                return True
        else:
            self.walk_count = max(0, self.walk_count - 1)
        return False

# ── PLASMA (GPIO PWM) ────────────────────────────────────────────────
class Plasma:
    def __init__(self, pin=18, pmin=0.0, pmax=0.35, enabled=True):
        self.pin, self.pmin, self.pmax, self.enabled = pin, pmin, pmax, enabled
        self._pwm = None
    def setup(self):
        if not self.enabled: return
        try:
            import RPi.GPIO as GPIO
            GPIO.setmode(GPIO.BCM)
            GPIO.setup(self.pin, GPIO.OUT)
            self._pwm = GPIO.PWM(self.pin, 100); self._pwm.start(0)
        except ImportError: pass
    def set(self, motion: float):
        d = self.pmin + (self.pmax - self.pmin) * motion
        if self._pwm: self._pwm.ChangeDutyCycle(d * 100)
        return d
    def shutdown(self):
        if self._pwm: self._pwm.stop()

# ── NETWORK ───────────────────────────────────────────────────────────
class Link:
    """TCP link with 4-byte length header. Frames carry embedded 8-byte microsecond timestamp."""
    HDR_LEN = 4 + 8  # length(uint32) + timestamp_us(uint64)

    def __init__(self, role: str, host: str, port: int, listen: int, w: int, h: int):
        self.role, self.w, self.h = role, w, h
        self._host, self._port, self._listen = host, port, listen
        self._sock: Optional[socket.socket] = None
        self.sent = self.bytes_ = 0
        self.rtt_ms = 0.0  # last round-trip estimate (ms)

    async def start(self) -> bool:
        loop = asyncio.get_running_loop()
        if self.role == "A":
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            s.setblocking(False)
            s.bind(("0.0.0.0", self._listen)); s.listen(4)
            log.info("PORTAL A: wait for B on :%d", self._listen)
            conn, addr = await loop.sock_accept(s)
            s.close()
            conn.setblocking(False)
            self._sock = conn
            log.info("PORTAL A: peer %s connected", addr)
            return True
        # role B
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        s.setblocking(False)
        for i in range(30):
            try:
                await loop.sock_connect(s, (self._host, self._port))
                self._sock = s
                log.info("PORTAL B: connected to A %s:%d", self._host, self._port)
                return True
            except (ConnectionRefusedError, OSError):
                log.warning("PORTAL B: retry %d/30…", i+1)
                await asyncio.sleep(1.0)
        return False

    async def send(self, frame: np.ndarray, quality: int = 80):
        loop = asyncio.get_running_loop()
        ts_us = int(time.time() * 1_000_000) & 0xFFFFFFFF_FFFFFFFF
        payload = encode(frame, quality)
        hdr = len(payload).to_bytes(4, 'big') + struct.pack('>Q', ts_us)
        try:
            await asyncio.wait_for(
                loop.sock_sendall(self._sock, hdr + payload), timeout=1.0)
        except (asyncio.TimeoutError, OSError):
            return  # drop frame, try again next loop
        self.sent += 1; self.bytes_ += len(hdr) + len(payload)

    async def recv(self) -> Optional[tuple]:
        """Returns (frame, send_timestamp_us) or None (timeout/closed)."""
        loop = asyncio.get_running_loop()
        try:
            hdr = await asyncio.wait_for(
                loop.sock_recv(self._sock, Link.HDR_LEN), timeout=0.5)
        except (asyncio.TimeoutError, OSError):
            return None
        if len(hdr) < Link.HDR_LEN: return None
        size = int.from_bytes(hdr[:4], 'big')
        if size <= 0 or size > 8_000_000: return None
        buf = b''
        try:
            while len(buf) < size:
                ch = await asyncio.wait_for(
                    loop.sock_recv(self._sock, size - len(buf)), timeout=0.5)
                if not ch: return None
                buf += ch
        except (asyncio.TimeoutError, OSError):
            return None
        ts_us = struct.unpack('>Q', hdr[4:])[0]
        frame = decode(buf, self.w, self.h)
        return frame, ts_us

    def close(self):
        try: self._sock.close()
        except: pass

# ── CALIBRATOR ────────────────────────────────────────────────────────
class Calibrator:
    """
    Prints a checkerboard to the projector; partner captures it and
    computes round-trip latency. Results are printed to stdout.
    """
    @staticmethod
    def make_checker(w: int, h: int, cols: int = 8, rows: int = 6, sq: int = 60) -> np.ndarray:
        board = np.zeros((h, w, 3), np.uint8)
        ox, oy = (w - cols * sq) // 2, (h - rows * sq) // 2
        for r in range(rows):
            for c in range(cols):
                if (r + c) % 2 == 0:
                    board[oy+r*sq:oy+(r+1)*sq, ox+c*sq:ox+(c+1)*sq] = 255
        return board

    @staticmethod
    async def measure_latency(loop, link: "Link", quality: int) -> float:
        """Send/receive single frame → measure RTT in ms."""
        board = Calibrator.make_checker(link.w, link.h)
        t0 = time.time()
        payload = encode(board, quality)
        hdr = len(payload).to_bytes(4, 'big') + struct.pack('>Q', int(t0*1e6))
        await loop.sock_sendall(link._sock, hdr + payload)
        # receive echo (partner sends the same frame back)
        hdr2 = await loop.sock_recv(link._sock, Link.HDR_LEN)
        size = int.from_bytes(hdr2[:4], 'big')
        buf = b''
        while len(buf) < size:
            buf += await loop.sock_recv(link._sock, size - len(buf))
        return (time.time() - t0) * 1000

# ── PORTAL NODE ───────────────────────────────────────────────────────
class PortalNode:
    def __init__(self, cfg: dict, hw: bool = True):
        self.cfg = cfg
        self.hw = hw
        n = cfg["network"]
        self.link = Link(n["role"], n["peer_host"], n["peer_port"],
                         n["listen_port"], cfg["camera"]["width"], cfg["camera"]["height"])
        self.overlay = make_overlay(cfg["camera"]["width"], cfg["camera"]["height"])
        self.walk = WalkDetector(cfg["camera"]["width"], cfg["camera"]["height"])
        self.plasma = Plasma(cfg["plasma"]["gpio_pin"],
                             cfg["plasma"]["power_min"],
                             cfg["plasma"]["power_max"], enabled=hw)
        self.cam: Optional[cv2.VideoCapture] = None
        self.running = True

    def open_camera(self):
        if not self.hw or not HAVE_CV2: return
        c = self.cfg["camera"]
        self.cam = cv2.VideoCapture(c["device"])
        self.cam.set(cv2.CAP_PROP_FRAME_WIDTH, c["width"])
        self.cam.set(cv2.CAP_PROP_FRAME_HEIGHT, c["height"])
        self.cam.set(cv2.CAP_PROP_FPS, c["fps"])
        self.cam.set(cv2.CAP_PROP_BUFFERSIZE, 1)
        log.info("Camera: /dev/video%d %dx%d@%d",
                 c["device"], c["width"], c["height"], c["fps"])

    async def run(self):
        self.plasma.setup()
        if self.hw: self.open_camera()
        if not await self.link.start():
            log.error("No connection. Exiting."); self.plasma.shutdown(); return
        role = self.cfg["network"]["role"]
        self.ts0 = self.last_log = time.time()
        log.info("PORTAL %s: WALK-THROUGH MODE ACTIVE  ←→  %s",
                 role, f"latency_budget={self.cfg['network']['latency_budget_ms']}ms")
        s_task = asyncio.create_task(self._send_loop())
        r_task = asyncio.create_task(self._recv_loop())
        try:
            await asyncio.gather(s_task, r_task)
        except (BrokenPipeError, ConnectionResetError):
            log.error("PORTAL %s: connection lost.", role)
        finally:
            self.running = False; s_task.cancel(); r_task.cancel()
            self.link.close(); self.plasma.shutdown()

    async def _send_loop(self):
        quality = self.cfg["network"]["quality"]
        log_interval = 5.0
        while self.running:
            if self.cam:
                ok, frame = self.cam.read()
                if not ok: await asyncio.sleep(0.01); continue
            else:
                # synthetic frame for --test
                w, h = self.cfg["camera"]["width"], self.cfg["camera"]["height"]
                frame = np.zeros((h, w, 3), np.uint8)
                frame[:] = (10, 10, 40)
                elapsed = time.time() - self.ts0
                # движущийся «прохожий» в центре — для проверки walk-детекции
                px = w // 2 + int(120 * np.sin(elapsed * 2.0))
                py = h // 2 + int(80 * np.cos(elapsed * 3.0))
                yy, xx = np.mgrid[0:h, 0:w]
                disc = (xx - px) ** 2 + (yy - py) ** 2 <= 60 ** 2
                frame[disc] = (40, 200, 40)
                if HAVE_CV2:
                    cv2.putText(frame,
                                f"PORTAL {self.cfg['network']['role']}  t={elapsed:.2f}s",
                                (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 1,
                                (0, 255, 150), 2)
            # walk detection
            if self.walk.detect(frame):
                log.info("PORTAL %s: >> WALK THROUGH DETECTED <<",
                         self.cfg['network']['role'])
            scene = overlay(frame, self.overlay)
            await self.link.send(scene, quality)
            now = time.time()
            if now - self.last_log > log_interval:
                dt = now - self.last_log
                log.info("PORTAL %s: %d frames, %.0f kB/s",
                         self.cfg['network']['role'], self.link.sent,
                         self.link.bytes_/1024/dt)
                self.link.bytes_ = 0; self.link.sent = 0; self.last_log = now

    async def _recv_loop(self):
        self.last_log = time.time()
        while self.running:
            res = await self.link.recv()
            if res is None: break
            frame, ts_us = res
            rtt_ms = (time.time() - ts_us/1e6) * 1000
            self.link.rtt_ms = 0.9 * self.link.rtt_ms + 0.1 * rtt_ms
            if self.cfg.get("show_local", False) and HAVE_CV2:
                disp = overlay(frame, self.overlay)
                if self.link.rtt_ms > 0:
                    cv2.putText(disp, f"RTT={self.link.rtt_ms:.0f}ms",
                                (30, self.cfg['camera']['height']-30),
                                cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0,255,200), 2)
                cv2.imshow(f"PORTAL {self.cfg['network']['role']} (remote)", disp)
                cv2.waitKey(1)
            await asyncio.sleep(0)

# ── CLI ───────────────────────────────────────────────────────────────
def load(path: Optional[str], overrides: dict) -> dict:
    cfg = json.loads(json.dumps(DEFAULT))
    if path and os.path.exists(path):
        with open(path) as f: u = json.load(f)
        def m(d, s):
            for k,v in s.items():
                if isinstance(v, dict) and isinstance(d.get(k), dict): m(d[k], v)
                else: d[k] = v
        m(cfg, u)
    for k,v in overrides.items():
        if v is not None: cfg["network"][k] = v
    return cfg

def main():
    ap = argparse.ArgumentParser(description="Portal Walk-Through V1.0")
    ap.add_argument("--config", default=None)
    ap.add_argument("--role", choices=["A","B"], default=None)
    ap.add_argument("--peer", default=None)
    ap.add_argument("--port", type=int, default=None)
    ap.add_argument("--simulate", action="store_true",
                    help="Hardware-less synthetic loopback test (two nodes, one process)")
    ap.add_argument("--calibrate", action="store_true",
                    help="Project checkerboard & measure latency on startup")
    ap.add_argument("--test", action="store_true", dest="test",
                    help="Quick 5-second loopback test (--simulate --run 5)")
    ap.add_argument("--run", type=int, default=30,
                    help="Seconds to run in test/simulate mode")
    args = ap.parse_args()

    if args.test:
        args.simulate = True; args.run = 5

    overrides = {k: v for k, v in
                 {"role": args.role, "peer_host": args.peer,
                  "peer_port": args.port, "listen_port": args.port}.items()
                 if v is not None}
    cfg = load(args.config, overrides)

    if args.simulate and args.role is None:
        # run both A and B in one process
        cfgA = json.loads(json.dumps(cfg)); cfgA["network"]["role"] = "A"
        cfgA["network"]["listen_port"] = 9100; cfgA["network"]["peer_port"] = 9100
        cfgB = json.loads(json.dumps(cfg)); cfgB["network"]["role"] = "B"
        cfgB["network"]["peer_host"] = "127.0.0.1"; cfgB["network"]["peer_port"] = 9100
        cfgB["network"]["listen_port"] = 9200

        log.info("=== SIMULATION MODE: two nodes, one process ===")
        log.info("Running %d seconds. Ctrl+C to stop.", args.run)

        nodeA = PortalNode(cfgA, hw=False)
        nodeB = PortalNode(cfgB, hw=False)
        async def runner():
            t = asyncio.create_task(nodeA.run())
            await asyncio.sleep(0.2)  # let A start listen
            u = asyncio.create_task(nodeB.run())
            await asyncio.sleep(args.run)
            nodeA.running = False; nodeB.running = False
            try:
                await asyncio.wait_for(
                    asyncio.gather(t, u, return_exceptions=True), timeout=10)
            except asyncio.TimeoutError:
                log.warning("runner: gather timeout — force cancel")
                t.cancel(); u.cancel()
                await asyncio.sleep(0.1)
            finally:
                nodeA.link.close(); nodeB.link.close()
            log.info("runner: gather done")
        try:
            asyncio.run(runner())
        except KeyboardInterrupt:
            pass
        log.info("=== SIMULATION COMPLETE ===")
        return

    node = PortalNode(cfg, hw=not args.simulate)
    try:
        asyncio.run(node.run())
    except KeyboardInterrupt:
        pass
    finally:
        if node.cam: node.cam.release()
        if HAVE_CV2: cv2.destroyAllWindows()

if __name__ == "__main__":
    main()
