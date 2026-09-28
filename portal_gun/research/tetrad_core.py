"""
tetrad_core.py — телепараллельное ядро для R3/R4/R1
====================================================
Диагональный войеров для статично-сферической метрики:
  ds² = A²dt² − B²dr² − C²dΩ²,  A=e^Φ, C=r, B=1/√(1−b(r)/r).

Даёт (символьно, затем численно):
  · T(r) — торсионный скаляр (мы проверяем: T = S^{μν}_ρ·T^ρ_{μν})
  · S^{μν}_ρ — суперпотенциал (нужен для ADM-энергии Малуфа)
  · выбор формы: 'ellis'  b=b0²/r  и  'const' b=b0
"""
import sympy as sp

r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r)
B = sp.Function('B')(r)
C = sp.Function('C')(r)

def dd(expr, i):
    if i == 1:
        return sp.diff(expr, r)
    if i == 2:
        return sp.diff(expr, th)
    return sp.Integer(0)

Etet = [
    [A, 0, 0, 0],
    [0, B, 0, 0],
    [0, 0, C, 0],
    [0, 0, 0, C * sp.sin(th)],
]
g = [[0]*4 for _ in range(4)]
ginv = [[0]*4 for _ in range(4)]
for m in range(4):
    g[m][m] = sp.simplify(sum(Etet[a][m]**2 * (1 if a in (0, 2) else -1) for a in range(4)))
    ginv[m][m] = sp.Integer(1) / g[m][m]
Einv = [[0]*4 for _ in range(4)]
for a in range(4):
    Einv[a][a] = sp.Integer(1) / Etet[a][a]

def Tcomp(rho, mu, nu):
    s = sp.S.Zero
    for a in range(4):
        e_nu = Etet[a][nu]
        e_mu = Etet[a][mu]
        s += Einv[rho][a] * (dd(e_nu, mu) - dd(e_mu, nu))
    return s

def T_low(sigma, mu, nu):
    return g[sigma][sigma] * Tcomp(sigma, mu, nu)

def K(mu, nu, rho):
    return sp.Rational(1, 2) * (T_low(mu, nu, rho) + T_low(rho, mu, nu) - T_low(nu, mu, rho))

def Superpotential(muu, nuu, rhol):
    """S^{μν}_ρ (μ,ν сверху охоты, ρ — нижний, индекс торсионного скаляра)."""
    Kup = sp.S.Zero
    for a in range(4):
        for b in range(4):
            Kup += ginv[muu][a] * ginv[nuu][b] * K(a, b, rhol)
    Tsvs = lambda nuz: sum(Tcomp(s, s, nuz) for s in range(4))
    Tsus = lambda muz: sum(Tcomp(s, muz, s) for s in range(4))
    return sp.Rational(1, 2) * (Kup
                                + sp.KroneckerDelta(muu, rhol) * Tsvs(nuu)
                                - sp.KroneckerDelta(nuu, rhol) * Tsus(muu))

def T_scalar_expr():
    s = sp.S.Zero
    for rho_ in range(4):
        for mu_ in range(4):
            for nu_ in range(4):
                s += Superpotential(mu_, nu_, rho_) * Tcomp(rho_, mu_, nu_)
    return s

def _etp(a, mu):
    return sp.Integer(1)/Etet[a][mu] if a == mu else sp.Integer(0)

# Бильдинг форм червоточины
def worm_exprs(mode):
    """Возвращает подстановки {B: Bexpr, ...} и производные для формы mode."""
    Bexpr = 1/sp.sqrt(1 - b0/r) if mode == 'const' else 1/sp.sqrt(1 - b0**2/r**2)
    return Bexpr

def S_component(muu, nuu, rhol):
    """Символьная компонента суперпотенциала в терминах A,B,C."""
    return Superpotential(muu, nuu, rhol)

def S_worm_numeric(muu, nuu, rhol, mode, b0v, rv, thetav=None, extra_sub=None):
    """Численное значение S^{μν}_ρ в червоточине mode."""
    ve = S_component(muu, nuu, rhol)
    Bexpr = worm_exprs(mode)
    Bp = sp.simplify(sp.diff(Bexpr, r))
    Bpp = sp.simplify(sp.diff(Bexpr, r, 2))
    e = ve
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bpp).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bexpr).subs(C, r)
    e = e.doit()
    if thetav is not None:
        e = e.subs(th, thetav)
    e = e.subs({r: rv, b0: b0v})
    return complex(e.evalf(30))

def T_worm_numeric(mode, b0v, rv, thetav=None):
    e = T_scalar_expr()
    Bexpr = worm_exprs(mode)
    Bp = sp.simplify(sp.diff(Bexpr, r))
    Bpp = sp.simplify(sp.diff(Bexpr, r, 2))
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bpp).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bexpr).subs(C, r).doit()
    if thetav is not None:
        e = e.subs(th, thetav)
    e = e.subs({r: rv, b0: b0v})
    return complex(e.evalf(30))

if __name__ == '__main__':
    print("self-test (GR-limit torsion scalar at r=2 b0, th=pi/4):")
    for mode in ['ellis', 'const']:
        print("  %-6s  T = %s" % (mode, T_worm_numeric(mode, 0.1, 0.2, sp.pi/4)))
    print("  S^{00r} (ADM-kernel) at r=2 b0:")
    for mode in ['ellis', 'const']:
        v = S_worm_numeric(0, 0, 1, mode, 0.1, 0.2, sp.pi/4)
        print("  %-6s  S00r = %s" % (mode, v))