"""
RESEARCH 1 — ПОЛЕВЫЕ УРАВНЕНИЯ f(T)-ГРАВИТАЦИИ ДЛЯ ЧЕРВОТОЧИНЫ (Symbolic EL)
=============================================================================
Действие:  S = ∫ d⁴x e [ T + f(T) ] / 2κ   (телепараллельная f(T)-гравитация)

Диагональный войеров (anisotropic fluid source, статично-сферически-симметрично):
  ds² = A(r)²dt² − B(r)²dr² − C(r)²(dθ² + sin²θ dφ²)
  A=e^Φ,  B=1/√(1−b/r),  C=r   (Моррис–Торн)

Отсюда по кручению Тейценбёка строится скаляр T и вариация ЭЛ:
  δL/δq − d/dr(δL/δq') = 0,  L = A·B·C²·[T + f(T)].

Проверка GR-предела: при f→0 уравнения TEGR обязаны совпасть с уравнениями GR
(эйнштейновский тензор, посчитанный Риччи) с точностью до одной нормировки κ.
"""
import time
import json
import sympy as sp

t_start = time.time()
print("=" * 78)
print("  RESEARCH 1 — f(T) WORMHOLE FIELD EQUATIONS (pure sympy EL)")
print("=" * 78)

# Символы и функции
r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r)
B = sp.Function('B')(r)
C = sp.Function('C')(r)

def dd(expr, i):
    """Частная производная по координате i: 0=t, 1=r, 2=θ, 3=φ."""
    if i == 1:
        return sp.diff(expr, r)
    if i == 2:
        return sp.diff(expr, th)
    return sp.Integer(0)

# ---- Vierbein (диагональный) и метрика ----
Etet = [
    [A, sp.S.Zero, sp.S.Zero, sp.S.Zero],
    [sp.S.Zero, B, sp.S.Zero, sp.S.Zero],
    [sp.S.Zero, sp.S.Zero, C, sp.S.Zero],
    [sp.S.Zero, sp.S.Zero, sp.S.Zero, C * sp.sin(th)],
]

g = [[sp.S.Zero] * 4 for _ in range(4)]
ginv = [[sp.S.Zero] * 4 for _ in range(4)]
for m in range(4):
    g[m][m] = sp.simplify(sum(Etet[a][m] ** 2 * (1 if a in (0, 2) else -1)
                              for a in range(4) if Etet[a][m] != 0))
    ginv[m][m] = sp.Integer(1) / g[m][m]

Einv = [[sp.S.Zero] * 4 for _ in range(4)]
for a in range(4):
    for m in range(4):
        if a == m:
            Einv[a][m] = sp.Integer(1) / Etet[a][a]

# ---- Кручение: T^ρ_{μν} = e^ρ_a (∂_μ e^a_ν − ∂_ν e^a_μ) ----
def Tcomp(rho, mu, nu):
    s = sp.S.Zero
    for a in range(4):
        ea_nu = Etet[a][nu]
        ea_mu = Etet[a][mu]
        if ea_nu != 0 and ea_mu != 0:
            s += Einv[rho][a] * (dd(ea_nu, mu) - dd(ea_mu, nu))
        elif ea_nu != 0:
            s += Einv[rho][a] * dd(ea_nu, mu)
        elif ea_mu != 0:
            s += Einv[rho][a] * (-dd(ea_mu, nu))
    return s

def T_low(sigma, mu, nu):
    """T_{σ μ ν} = g_{σ ρ} T^ρ_{μ ν}."""
    if g[sigma][sigma] == 0:
        return sp.S.Zero
    return g[sigma][sigma] * Tcomp(sigma, mu, nu)

# ---- Конторсия ----
def K(mu, nu, rho):
    """K_{μνρ} = ½(T_{μνρ} + T_{ρμν} − T_{νμρ})"""
    return sp.Rational(1, 2) * (T_low(mu, nu, rho) + T_low(rho, mu, nu)
                                - T_low(nu, mu, rho))

def T_sigma_nu_sigma(nuu):
    """T^{σν}{}_σ = Σ_σ Σ_β g^{νβ} T^σ_{βσ}"""
    s = sp.S.Zero
    for si in range(4):
        for be in range(4):
            if ginv[nuu][be] != 0:
                s += ginv[nuu][be] * Tcomp(si, be, si)
    return s

def T_sigma_mu_sigma(muu):
    """T^{σμ}{}_σ = Σ_σ Σ_β g^{μβ} T^σ_{βσ}"""
    s = sp.S.Zero
    for si in range(4):
        for be in range(4):
            if ginv[muu][be] != 0:
                s += ginv[muu][be] * Tcomp(si, be, si)
    return s

# ---- Суперпотенциал и скаляр T ----
def Superpotential(muu, nuu, rhol):
    """S^{μν}_ρ = ½( K^{μν}_ρ + δ_ρ^μ T^{σν}{}_σ − δ_ρ^ν T^{σμ}{}_σ )"""
    Kup = sp.S.Zero
    for a in range(4):
        for b in range(4):
            if ginv[muu][a] != 0 and ginv[nuu][b] != 0:
                Kup += ginv[muu][a] * ginv[nuu][b] * K(a, b, rhol)
    return sp.Rational(1, 2) * (Kup
                                + sp.KroneckerDelta(muu, rhol) * T_sigma_nu_sigma(nuu)
                                - sp.KroneckerDelta(nuu, rhol) * T_sigma_mu_sigma(muu))

# Скаляр T = S^{μν}_ρ · T^ρ_{μν}
T_scalar = sp.S.Zero
for rho_ in range(4):
    for mu_ in range(4):
        for nu_ in range(4):
            T_scalar += Superpotential(mu_, nu_, rho_) * Tcomp(rho_, mu_, nu_)

print("  T построен, размер выражения: %d тыс. символов" % (len(str(T_scalar)) // 1000))
T_scalar = sp.simplify(T_scalar)
print("  После simplify:")
print(T_scalar)
print()

# ============================================================
# ВАРИАЦИЯ ЭЙЛЕРА–ЛАГРАНЖА
# ============================================================
F = sp.Function('f')
L = A * B * C ** 2 * (T_scalar + F(T_scalar))

EL = {}
for name, q in [('A', A), ('B', B), ('C', C)]:
    qp = q.diff(r)
    dLdq = sp.diff(L, q)
    dLdqp = sp.diff(L, qp)
    EL[name] = sp.simplify(dLdq - sp.diff(dLdqp, r))
    print("  EL[%s]: %d тыс. символов" % (name, len(str(EL[name])) // 1000))

# ============================================================
# GR-ПРОВЕРКА (f→0) — через эйнштейновский тензор G_μν
# ============================================================
print()
print("  [ПРОВЕРКА] f→0 ⇒ TEGR = GR")
gmt = [[A**2, 0, 0, 0],
       [0, -B**2, 0, 0],
       [0, 0, -C**2, 0],
       [0, 0, 0, -C**2 * sp.sin(th)**2]]
gmt2 = [[sp.S.Zero]*4 for _ in range(4)]
for m in range(4):
    gmt2[m][m] = sp.Integer(1) / gmt[m][m]

def chris(l, m, n):
    s = sp.S.Zero
    for sig in range(4):
        if gmt2[l][sig] != 0:
            s += gmt2[l][sig] * (dd(gmt[sig][n], m) + dd(gmt[sig][m], n) - dd(gmt[m][n], sig))
    return sp.Rational(1, 2) * s

def ricci(mu, nu):
    """R_{μν} = Σ_l [∂_l Γ^l_{μν} − ∂_ν Γ^l_{μl}]
                  + Σ_{l,sig} [Γ^l_{l,sig}Γ^sig_{μν} − Γ^l_{ν,sig}Γ^sig_{l,μ}]
    ВАЖНО: второй член (∂_ν…) ОДИН на каждый l, НЕ внутри цикла sig."""
    s = sp.S.Zero
    for l in range(4):
        s += dd(chris(l, mu, nu), l) - dd(chris(l, mu, l), nu)
        for sig in range(4):
            s += chris(l, l, sig) * chris(sig, mu, nu)
            s -= chris(l, nu, sig) * chris(sig, l, mu)
    return s

R_scalar_ = sp.S.Zero
for a in range(4):
    for b in range(4):
        if gmt2[a][b] != 0:
            R_scalar_ += gmt2[a][b] * ricci(a, b)

def G(mu, nu):
    return ricci(mu, nu) - sp.Rational(1, 2) * gmt[mu][nu] * R_scalar_

print("  (Риччи построен)")

# ============================================================
# ЧИСЛЕННОЕ СРАВНЕНИЕ ПРИ f=0
# ============================================================
# EL при f→0, f'→0 — ЧИСЛОВОЙ путь (без тяжёлого simplify)
# ============================================================
def _worm_num(expr, rv, bv, thv=None):
    """MT-ответ полностью численно: A=1, B=1/sqrt(1-b0^2/r^2), C=r.
    Порядок ОБЯЗАТЕЛЕН: производные-узлы → функции → координаты (иначе
    A(0.2)-узлы не матчатся). Без тяжёлого simplify.
    [ТР3-ФИКС] B' = -½(1-b0^2/r^2)^{-3/2}·(2b0^2/r^3); у знаков Bpn и Bppn
    были перепутаны ЗНАКИ (было: Bpn=+1.92, правда −1.92).
    """
    thn = thv if thv is not None else sp.pi/4
    db = 1 - bv**2/rv**2
    dbp = 2*bv**2/rv**3
    dbpp = -6*bv**2/rv**4
    Bn = db**-0.5
    Bpn = -0.5*db**-1.5*dbp
    Bppn = 0.75*db**-2.5*dbp**2 - 0.5*db**-1.5*dbpp
    e = expr
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bppn).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bpn).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bn).subs(C, rv)
    e = e.doit()
    e = e.subs({r: rv, b0: bv, th: thn})
    return complex(e.evalf(30))

def EL0(name, rv, bv):
    e = EL[name]
    e = e.replace(lambda t: isinstance(t, sp.Function)
                  and getattr(t, 'func', None) == sp.Function('f'),
                  lambda t: sp.S.Zero)
    return _worm_num(e, rv, bv)

def Gval(mu, nu, rv, bv):
    return _worm_num(G(mu, nu), rv, bv)

# ---- кэш символьных EL (для быстрой отладки без перевывода) ----
out = {'T_scalar': sp.srepr(T_scalar),
       'EL_A': sp.srepr(EL['A']), 'EL_B': sp.srepr(EL['B']), 'EL_C': sp.srepr(EL['C'])}
with open('/home/smboozha/portal_gun/research/ft_el_cache.json', 'w', encoding='utf-8') as f:
    json.dump(out, f)
print("  Кэш: research/ft_el_cache.json сохранён (%d щелчков)" % sum(len(v) for v in out.values()))
print("  Время вывода: %.1f с" % (time.time() - t_start))

rv, bv = 0.2, 0.1
print("\n  Сравнение EL(0) vs G_μν при r=%s, b0=%s:" % (rv, bv))
try:
    for comp in [(0,0),(1,1),(2,2),(3,3)]:
        print("    debug R_%d%d = %+.6f" % (comp[0], comp[1], _worm_num(ricci(*comp), rv, bv).real))
    vR = _worm_num(R_scalar_, rv, bv).real
    print("    debug R = %+.6f" % vR)
    for comp in [(0,0),(1,1),(2,2)]:
        vG = _worm_num(G(*comp), rv, bv).real
        print("    debug G_%d%d = %+.6f" % (comp[0], comp[1], vG))
    for name, (m_, n_) in [('A', (0, 0)), ('B', (1, 1)), ('C', (2, 2))]:
        e0 = EL0(name, rv, bv).real
        g0 = Gval(m_, n_, rv, bv).real
        ratio = e0 / g0 if abs(g0) > 1e-12 else float('nan')
        print("    EL[%s]=%.6e   G_%d%d=%.6e   ratio=%.4f" % (name, e0, m_, n_, g0, ratio))
    print("\n  (если все три ratio ~ совпадают и ≠0 — TEGR совпал с GR, нормировка κ общая)")
except Exception as exc:
    print("  РАССЛЕДОВАНИЕ R1 (float-баг):", type(exc).__name__, exc)
print("  Время полное: %.1f с" % (time.time() - t_start))
print("=" * 78)