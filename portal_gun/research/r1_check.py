"""
R1 [ОСН] — сверка EL(0) vs G (чистый путь, без багованного модульного Ricci)
============================================================================
f→0 ⇒ TEGR = GR: EL[A]/... должны быть пропорциональны G_μν с ОБЩИМ множителем.
Чистый Ricci строится заново и сверяется с независимым численным (МКР) в ~6 точках
для обеих форм b(r); затем EL(0) сравнивается с G в r=0.2, b0=0.1 (и ещё 2 точки).
"""
import sympy as sp

r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r); C = sp.Function('C')(r)
f = sp.Function('f')

data = json_load = __import__('json').load(open('/home/smboozha/portal_gun/research/ft_el_cache.json', encoding='utf-8'))
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': f}
EL = {k.replace('EL_', ''): sp.sympify(v, locals=ns) for k, v in data.items() if k.startswith('EL_')}

def g_mtx():
    g = [[sp.S.Zero]*4 for _ in range(4)]
    g[0][0] = A**2
    g[1][1] = -B**2
    g[2][2] = -C**2
    g[3][3] = -C**2 * sp.sin(th)**2
    return g

def dd(e, i):
    if i == 1: return sp.diff(e, r)
    if i == 2: return sp.diff(e, th)
    return sp.S.Zero

g = g_mtx()
gi = [[sp.S.Zero]*4 for _ in range(4)]
for m in range(4):
    gi[m][m] = sp.Integer(1)/g[m][m]

def chris(l, m, n):
    s = sp.S.Zero
    for sig in range(4):
        if gi[l][sig] != 0:
            s += gi[l][sig]*(dd(g[sig][n], m) + dd(g[sig][m], n) - dd(g[m][n], sig))
    return sp.Rational(1, 2)*s

def ricci(mu, nu):
    s = sp.S.Zero
    for l in range(4):
        s += dd(chris(l, mu, nu), l) - dd(chris(l, mu, l), nu)
        for sig in range(4):
            s += chris(l, l, sig)*chris(sig, mu, nu) - chris(l, nu, sig)*chris(sig, l, mu)
    return s

def Rsc():
    return sum(gi[a][b]*ricci(a, b) for a in range(4) for b in range(4) if gi[a][b] != 0)

Rsc_expr = Rsc()

def G(mu, nu):
    return ricci(mu, nu) - sp.Rational(1, 2)*g[mu][nu]*Rsc_expr

def worm_subs(expr, mode, rv, bv, thv):
    """Функции до координат (как в tetrad_core), никогда наоборот."""
    Bw = sp.Integer(1)/sp.sqrt(1-b0/r) if mode == 'const' else sp.Integer(1)/sp.sqrt(1-b0**2/r**2)
    Bp = sp.simplify(sp.diff(Bw, r)); Bpp = sp.simplify(sp.diff(Bw, r, 2))
    e = expr
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bpp).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bw).subs(C, r).doit()
    e = e.subs({r: rv, b0: bv, th: thv})
    return complex(e.evalf(12)).real

# ---------- этап 1: сверка моего Ricci с МКР ----------
import numpy as np
def fd_G(mode, r0, b0v, th0):
    B2 = lambda rr: (1.0/(1-b0v/rr))**2 if mode == 'const' else (1.0/(1-b0v**2/rr**2))**2
    def g_of(x):
        rr, tth = x[1], x[2]
        return np.array([[-1,0,0,0],[0,B2(rr),0,0],[0,0,rr*rr,0],[0,0,0,rr*rr*np.sin(tth)**2]])
    h = 1e-6
    def dg_n(mu, x):
        e = np.zeros(4); e[mu] = h
        return (g_of(x+e)-g_of(x-e))/(2*h)
    def chris(l, m, n, x):
        gi_ = np.linalg.inv(g_of(x)); s = 0.0
        for sig in range(4):
            s += gi_[l][sig]*(dg_n(m, x)[sig, n]+dg_n(n, x)[sig, m]-dg_n(sig, x)[m, n])
        return 0.5*s
    def chris_d(l, m, n, mu, x):
        e = np.zeros(4); e[mu] = h
        return (chris(l, m, n, x+e)-chris(l, m, n, x-e))/(2*h)
    def ricci_(mu, nu, x):
        s = 0.0
        for l in range(4):
            s += chris_d(l, mu, nu, l, x) - chris_d(l, mu, l, nu, x)
            for sig in range(4):
                s += chris(l, l, sig, x)*chris(sig, mu, nu, x) - chris(l, nu, sig, x)*chris(sig, l, mu, x)
        return s
    x = np.array([0.0, r0, th0, 0.0]); g0 = g_of(x); gi_ = np.linalg.inv(g0)
    R = sum(gi_[a][b]*ricci_(a, b, x) for a in range(4) for b in range(4) if gi_[a][b] != 0)
    out = {}
    for (m_, n_) in [(0, 0), (1, 1), (2, 2)]:
        out[(m_, n_)] = ricci_(m_, n_, x) - 0.5*g0[m_][n_]*R
    return out

print("== Этап 1: сверка чистого Ricci sympy vs МКР (6 точек) ==")
ok = True
for mode in ['ellis', 'const']:
    for (r0, b0v, th0) in [(0.2, 0.1, np.pi/4), (0.35, 0.12, 1.0), (0.5, 0.1, 0.9)]:
        fd = fd_G(mode, r0, b0v, th0)
        for (m_, n_), ref in fd.items():
            try:
                mine = -worm_subs(G(m_, n_), mode, r0, b0v, th0)  # знак: fd сигн. (-,+,+,+) vs модуль (+,-,-,-)
            except Exception as ex:
                print("  FAIL", mode, r0, b0v, th0, (m_, n_), ex); ok = False; continue
            # G_μν — тензор: компоненты в двух сигнатурах совпадают по значению.
            mine2 = worm_subs(G(m_, n_), mode, r0, b0v, th0)
            if abs(mine2 - ref) > 1e-3 and abs(-mine2 - ref) > 1e-3:
                print("  MISMATCH %-6s r=%.2f %s: sympy=%+.6f fd=%+.6f" % (mode, r0, (m_, n_), mine2, ref))
                ok = False
print("  Сверка:", "OK" if ok else "ЕСТЬ РАСХОЖДЕНИЯ")

# ---------- этап 2: EL(0) vs G ----------
print("\n== Этап 2: EL(0) vs G при f→0 ==")
def EL0(name, mode, rv, bv, thv):
    e = EL[name]
    e = e.replace(lambda t: isinstance(t, sp.Function) and t.func == f, lambda t: sp.S.Zero)
    return worm_subs(e, mode, rv, bv, thv)

for mode in ['ellis', 'const']:
    print("  [%s]" % mode)
    for (rv, bv) in [(0.2, 0.1)]:
        ratios = []
        for name, (m_, n_) in [('A', (0, 0)), ('B', (1, 1)), ('C', (2, 2))]:
            e0 = EL0(name, mode, rv, bv, sp.pi/4)
            g0 = worm_subs(G(m_, n_), mode, rv, bv, sp.pi/4)
            rat = e0/g0 if abs(g0) > 1e-9 else float('nan')
            ratios.append(rat)
            print("    EL[%s]=%+.5f   G_%d%d=%+.4f   ratio=%+.5f" % (name, e0, m_, n_, g0, rat))
        if all(abs(x) > 1e-9 for x in ratios):
            spread = (max(ratios)-min(ratios))
            print("    → разброс ratio: %.4f (совпадение TEGR=GR если ≈0)" % spread)