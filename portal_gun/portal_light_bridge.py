#!/usr/bin/env python3
"""
PORTAL LIGHT BRIDGE — контроллер двустороннего портала (V1.0).

Два устройства A и B, каждое с:
  - камерой (захват сцены перед порталом)
  - проектором (показ сцены с противоположной стороны)
  - кольцевым плазмотроном (зелёное свечение по краю)

Топология:
   [ПОРТАЛ A] <──→ МОСТ (P2P, WiFi/5G/оптоволокно) <──→ [ПОРТАЛ B]

  Камера A ──► Проектор B   (видно, что «за» порталом A)
  Камера B ──► Проектор A   (видно, что «за» порталом B)

Работает на Raspberry Pi 5 / Jetson Orin + камера CSI/USB.
"""

import argparse
import asyncio
import json
import logging
import os
import socket
import struct
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np

try:
    import cv2
    HAVE_CV2 = True
except ImportError:
    HAVE_CV2 = False
    logging.getLogger('PORTAL').warning("OpenCV не найден — работаю в режиме raw (без сжатия)")

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
)
log = logging.getLogger('PORTAL')

# ============================================================
# CONFIG
# ============================================================
DEFAULT_CONFIG = {
    "camera": {
        "device": 0,          # /dev/video0 или CSI
        "width": 1280,
        "height": 720,
        "fps": 60,            # необходимо для низкой задержки
    },
    "projector": {
        "device": "EDID-дисплей / HDMI 1.4+",
        "refresh": 60,
    },
    "plasma": {
        "gpio_pin": 18,       # для плазмотрона (ШИМ)
        "power_min": 0.0,     # 0..1
        "power_max": 0.35,
    },
    "network": {
        "role": "A",          # 'A' или 'B'
        "peer_host": "10.0.0.2",
        "peer_port": 9000,
        "listen_port": 9000,
        "bind": "0.0.0.0",
        "codec": "MJPEG",     # JPEG-кадры (низкая задержка)
        "quality": 80,        # 0..100
        "keyframe_interval": 30,
        "delay_budget_ms": 80,  # целевая задержка видео
    },
    "motion": {
        "sensitivity": 0.15,  # порог детекта движения
    },
}

# ============================================================
# CODEC: JPEG (cv2) или raw bytes (фолбек без OpenCV)
# ============================================================
def encode_frame(frame: np.ndarray, quality: int) -> bytes:
    """Кодирование кадра. С cv2 -> JPEG; без cv2 -> raw bytes."""
    if HAVE_CV2:
        ok, jpg = cv2.imencode('.jpg', frame,
                               [int(cv2.IMWRITE_JPEG_QUALITY), quality])
        return jpg.tobytes() if ok else frame.tobytes()
    return frame.tobytes()


def decode_frame(data: bytes, w: int, h: int, channels: int = 3) -> np.ndarray:
    """Декодирование кадра."""
    if HAVE_CV2:
        arr = np.frombuffer(data, dtype=np.uint8)
        dec = cv2.imdecode(arr, cv2.IMREAD_COLOR)
        if dec is not None:
            return dec
    expected = w * h * channels
    if len(data) == expected:
        return np.frombuffer(data, dtype=np.uint8).reshape(h, w, channels)
    # Данные повреждены или другого размера — пустой кадр
    return np.zeros((h, w, channels), dtype=np.uint8)


# ============================================================
# GREEN PORTAL FRAME OVERLAY
# ============================================================
def make_portal_overlay(width: int, height: int, r0: float = 0.38) -> np.ndarray:
    """Зелёное кольцо портала + свечение. Только край — центр прозрачен."""
    overlay = np.zeros((height, width, 3), dtype=np.uint8)
    cx, cy = width // 2, height // 2
    R = min(width, height) * r0

    # Матрица расстояний до центра
    yy, xx = np.mgrid[0:height, 0:width]
    dist = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)

    # Внешнее свечение (несколько колец)
    for rr in np.linspace(R * 0.98, R * 1.28, 12):
        band = np.abs(dist - rr) < 2.0
        overlay[band] = np.maximum(overlay[band],
                                   np.array([0, 255, 100], dtype=np.uint8))

    # Основное кольцо (толстое + тонкое яркое)
    band = np.abs(dist - R) < 3.0
    overlay[band] = np.array([0, 255, 150], dtype=np.uint8)
    band2 = np.abs(dist - R) < 0.9
    overlay[band2] = np.array([140, 255, 210], dtype=np.uint8)

    # Маленькие «молнии» на кольце (короткие дуги)
    rng = np.random.default_rng(int(time.time() % 1000))
    for _ in range(6):
        ang = rng.uniform(0, 2 * np.pi)
        a0 = ang - 0.12
        a1 = ang + rng.uniform(0.1, 0.3)
        seg = np.linspace(a0, a1, 10)
        sx = cx + R * np.cos(seg)
        sy = cy + R * np.sin(seg)
        for i in range(len(seg) - 1):
            band = np.abs(dist - np.hypot(sx[i], sy[i])) < 1.0
            band &= (np.abs(xx - sx[i]) < 0.9) & (np.abs(yy - sy[i]) < 0.9)
            overlay[band] = np.array([120, 255, 200], dtype=np.uint8)

    return overlay


def apply_overlay(frame: np.ndarray, overlay: np.ndarray, alpha_glow: float = 0.6) -> np.ndarray:
    """Накладывает зелёное кольцо поверх кадра."""
    mask = np.any(overlay > 0, axis=2)
    out = frame.copy()
    glow = (overlay.astype(np.float32) * alpha_glow).astype(np.uint8)
    out[mask] = np.clip(frame[mask].astype(np.int16)
                        + glow[mask].astype(np.int16), 0, 255).astype(np.uint8)
    return out


# ============================================================
# MOTION / PROXIMITY (для «яркости» при подходе)
# ============================================================
class MotionDetector:
    def __init__(self, sensitivity: float):
        self.sensitivity = sensitivity
        self._bg: Optional[np.ndarray] = None
        if HAVE_CV2:
            self._sub = cv2.createBackgroundSubtractorMOG2(history=20,
                                                           varThreshold=25,
                                                           detectShadows=False)

    def update(self, frame: np.ndarray) -> float:
        """Возвращает 0..1 — насколько сильное движение в кадре."""
        if HAVE_CV2:
            gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            fg = self._sub.apply(gray, learningRate=0.1)
            fg = fg > 200
            ratio = np.count_nonzero(fg) / fg.size
            return min(1.0, ratio / self.sensitivity)
        # Без cv2 — грубая оценка изменения яркости кадра за время
        gray = frame.mean(axis=2)
        if self._bg is None:
            self._bg = gray.copy()
            return 0.0
        diff = np.abs(gray.astype(np.int16) - self._bg.astype(np.int16)).mean()
        self._bg = 0.9 * self._bg + 0.1 * gray
        return min(1.0, diff / (self.sensitivity * 255.0 * 40.0))


# ============================================================
# NETWORK: P2P битовый мост
# ============================================================
class PortalLink:
    """
    P2P соединение между порталами A и B.
    Один узел слушает (listen), второй подключается (peer).
    Кадры сжимаются в JPEG и шлются без обратной связи (UDP-ish).
    """

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.role = cfg["network"]["role"]
        self.quality = cfg["network"]["quality"]
        self._sock: Optional[socket.socket] = None
        self._timer = 0.0
        self.frames_sent = 0
        self.bytes_sent = 0

    async def start(self) -> bool:
        cfg = self.cfg["network"]

        if self.role == "A":
            # A — сервер-слушатель
            self._sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            self._sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self._sock.bind((cfg["bind"], cfg["listen_port"]))
            self._sock.listen(1)
            log.info("Портал A: ожидаю подключение B на :%d", cfg["listen_port"])
            conn, addr = self._sock.accept()
            self._sock = conn
            log.info("Портал A: подключён партнёр %s", addr)
            return True
        else:
            # B — клиент
            self._sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            for attempt in range(30):
                try:
                    self._sock.connect((cfg["peer_host"], cfg["peer_port"]))
                    log.info("Портал B: подключён к A %s:%d",
                             cfg["peer_host"], cfg["peer_port"])
                    return True
                except (ConnectionRefusedError, OSError):
                    log.warning("Портал B: попытка %d/30, повтор...", attempt + 1)
                    await asyncio.sleep(1.0)
            return False

    async def send_frame(self, frame: np.ndarray) -> None:
        """Кодирует кадр и отправляет партнёру."""
        data = encode_frame(frame, self.quality)
        header = len(data).to_bytes(4, 'big')
        try:
            self._sock.sendall(header)
            self._sock.sendall(data)
            self.frames_sent += 1
            self.bytes_sent += len(data)
        except (BrokenPipeError, ConnectionResetError) as e:
            log.error("Соединение потеряно: %s", e)
            raise

    async def recv_frame(self) -> Optional[np.ndarray]:
        """Принимает кадр от партнёра."""
        try:
            header = await asyncio.to_thread(self._sock.recv, 4)
            if not header or len(header) < 4:
                return None
            size = int.from_bytes(header, 'big')
            if size <= 0 or size > 8_000_000:
                return None
            buf = b''
            while len(buf) < size:
                chunk = await asyncio.to_thread(self._sock.recv, size - len(buf))
                if not chunk:
                    return None
                buf += chunk
            return decode_frame(buf,
                                self.cfg["camera"]["width"],
                                self.cfg["camera"]["height"])
        except (BrokenPipeError, ConnectionResetError) as e:
            log.error("Чтение прервано: %s", e)
            raise

    def close(self):
        if self._sock:
            try:
                self._sock.close()
            except OSError:
                pass

# ============================================================
# PLASMA RING CONTROL (GPIO PWM)
# ============================================================
class PlasmaRing:
    """Кольцевой плазмотрон: ШИМ по GPIO. Яркость зависит от движения."""

    def __init__(self, cfg: dict, enabled: bool = True):
        self.enabled = enabled
        self.pin = cfg.get("gpio_pin", 18)
        self.power_min = cfg.get("power_min", 0.0)
        self.power_max = cfg.get("power_max", 0.35)
        self.power = 0.0
        self._pwm = None

    def setup(self):
        if not self.enabled:
            return
        try:
            import RPi.GPIO as GPIO
            GPIO.setmode(GPIO.BCM)
            GPIO.setup(self.pin, GPIO.OUT)
            self._pwm = GPIO.PWM(self.pin, 100)   # 100 Гц
            self._pwm.start(0)
            log.info("Плазмотрон: инициализирован (GPIO %d)", self.pin)
        except ImportError:
            log.warning("Плазмотрон: RPi.GPIO не найден — симуляция (quiet)")

    def update(self, motion: float):
        """motion = 0..1; плазмотрон разгорается при активности."""
        duty = self.power_min + (self.power_max - self.power_min) * motion
        self.power = duty
        if self._pwm:
            self._pwm.ChangeDutyCycle(duty * 100)
        else:
            # Тихо (лог каждые ~2с)
            now = time.time()
            if now - getattr(self, '_last_log', 0) > 2.0:
                log.info("Плазмотрон: power=%.2f (motion=%.2f)", duty, motion)
                self._last_log = now

    def shutdown(self):
        if self._pwm:
            self._pwm.stop()

# ============================================================
# MAIN LOOP (портал-узел)
# ============================================================
class PortalNode:
    def __init__(self, cfg: dict, enable_hw: bool = True):
        self.cfg = cfg
        self.enable_hw = enable_hw
        self.link = PortalLink(cfg)
        self.plasma = PlasmaRing(cfg.get("plasma", {}), enabled=enable_hw)
        self.motion = MotionDetector(cfg["motion"]["sensitivity"])
        self.cam: Optional[cv2.VideoCapture] = None
        self.overlay = make_portal_overlay(cfg["camera"]["width"],
                                           cfg["camera"]["height"])
        self.on_proximity = None   # колбэк: при «проходе» человека
        self.running = True

    def open_camera(self):
        c = self.cfg["camera"]
        if HAVE_CV2:
            self.cam = cv2.VideoCapture(c["device"])
            self.cam.set(cv2.CAP_PROP_FRAME_WIDTH, c["width"])
            self.cam.set(cv2.CAP_PROP_FRAME_HEIGHT, c["height"])
            self.cam.set(cv2.CAP_PROP_FPS, c["fps"])
            self.cam.set(cv2.CAP_PROP_BUFFERSIZE, 1)
            log.info("Камера: открыта /dev/video%d %dx%d@%d",
                     c["device"], c["width"], c["height"], c["fps"])
        else:
            log.warning("Камера: OpenCV не установлен — используйте --simulate")

    async def run(self):
        self.plasma.setup()
        if self.enable_hw:
            self.open_camera()
        if not await self.link.start():
            log.error("Портал: не удалось установить связь. Выход.")
            self.plasma.shutdown()
            return

        log.info("Портал: МОСТ АКТИВЕН")
        self.last_stats = time.time()

        send_task = asyncio.create_task(self._send_loop())
        recv_task = asyncio.create_task(self._recv_loop())

        try:
            await asyncio.gather(send_task, recv_task)
        except (BrokenPipeError, ConnectionResetError):
            log.error("Портал: соединение оборвано.")
        finally:
            self.running = False
            send_task.cancel()
            recv_task.cancel()
            self.link.close()
            self.plasma.shutdown()

    async def _send_loop(self):
        """Захват, оверлей, движение → отправка."""
        while self.running:
            if self.cam:
                ok, frame = self.cam.read()
                if not ok:
                    await asyncio.sleep(0.01)
                    continue
            else:
                # Режим без железа: синтетический кадр (тест сети)
                frame = np.zeros((self.cfg["camera"]["height"],
                                  self.cfg["camera"]["width"], 3),
                                 dtype=np.uint8)
                frame[:] = (10, 10, 40)
                ts = time.time()
                if HAVE_CV2:
                    cv2.putText(frame,
                                f"PORTAL {self.cfg['network']['role']} t={ts:.2f}",
                                (50, 50), cv2.FONT_HERSHEY_SIMPLEX, 1,
                                (0, 255, 150), 2)

            motion = self.motion.update(frame)
            self.plasma.update(motion)
            scene = apply_overlay(frame, self.overlay)
            await self.link.send_frame(scene)

            now = time.time()
            if now - self.last_stats > 5.0:
                log.info("Портал: отправлено %d кадров, %.1f кБ/с",
                         self.link.frames_sent,
                         self.link.bytes_sent / 5.0 / 1024)
                self.last_stats = now

    async def _recv_loop(self):
        """Приём + отображение (у головы) + обратный оверлей для отрисовки."""
        while self.running:
            frame = await self.link.recv_frame()
            if frame is None:
                break
            # Показывает партнёрскую сцену (можно развернуть в окно)
            if self.cfg.get("show_local", False) and HAVE_CV2:
                display = apply_overlay(frame, self.overlay)
                cv2.imshow("PORTAL (обратная сцена)", display)
                cv2.waitKey(1)
            await asyncio.sleep(0)


# ============================================================
# CLI + CONFIG LOAD
# ============================================================
def load_config(path: Optional[str], overrides: dict) -> dict:
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))  # deep copy
    if path and os.path.exists(path):
        with open(path) as f:
            user_cfg = json.load(f)
        # рекурсивный merge
        def merge(dst, src):
            for k, v in src.items():
                if isinstance(v, dict) and isinstance(dst.get(k), dict):
                    merge(dst[k], v)
                else:
                    dst[k] = v
        merge(cfg, user_cfg)
    for k, v in overrides.items():
        if v is not None:
            cfg["network"][k] = v
    return cfg


def main():
    ap = argparse.ArgumentParser(description="Portal Light Bridge (V1.0)")
    ap.add_argument("--config", default=None, help="JSON файл конфигурации")
    ap.add_argument("--role", choices=["A", "B"], default=None,
                    help="Роль узла: A (слушатель) / B (клиент)")
    ap.add_argument("--peer", default=None, help="Адрес партнёра (для B)")
    ap.add_argument("--port", type=int, default=None, help="Порт")
    ap.add_argument("--simulate", action="store_true",
                    help="Запуск без железа (тест сети + видео)")
    args = ap.parse_args()

    overrides = {"role": args.role, "peer_host": args.peer}
    if args.port is not None:
        overrides["peer_port"] = args.port
        overrides["listen_port"] = args.port
    overrides = {k: v for k, v in overrides.items() if v is not None}
    cfg = load_config(args.config, overrides)

    node = PortalNode(cfg, enable_hw=not args.simulate)
    try:
        asyncio.run(node.run())
    except KeyboardInterrupt:
        log.info("Портал: остановлен оператором.")
    finally:
        if node.cam:
            node.cam.release()
        cv2.destroyAllWindows()


if __name__ == "__main__":
    main()