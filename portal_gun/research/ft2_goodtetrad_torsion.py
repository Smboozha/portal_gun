"""
R2 [ПОМ] — независимая проверка: «хорошая» vs «диагональная» тетрада в f(T)
============================================================================
Мотивация (Tamanini & Boehmer, arXiv:1204.4593): для статической сферически-
симметричной метрики ДИАГОНАЛЬНАЯ тетрада Вайтценбёка в f(T)-гравитации является
"bad tetrad": уравнения переопределены / зажимают f(T) до линейной (TEGR).
Главная модель [ОСН] строит R1 именно на диагональной тетраде — проверяем,
не испортила ли это зависимость поля от f(T).

Здесь: свежая независимая реализация (другой путь, чем R1)
  1) «good» (symmetric, O-rotated) тетрада  →  T_good(r)
  2) диагональная тетрада (как в R1)        →  T_bad(r)
  3) сравнение + численная сверка на Morris–Thorne (A=1, B=1/sqrt(1-b0^2/r^2))
"""
import time, json
import sympy as sp

t0 = time.time()
r, th, ph, b0 = sp.symbols('r th ph b0', positive=True)

def dd(expr, var):
    return sp.diff(expr, var)

# ---------------------------------------------------------------------------
# Метрика Morris-Thorne: A=1, B = 1/sqrt(1 - b0^2/r^2),  (т.е. b(r)=b0^2/r)
# ---------------------------------------------------------------------------
Aval = sp.Integer(1)
Bval = 1 / sp.sqrt(1 - b0**2 / r**2)

def metric_from_tetrad(e):
    g = [[sp.S.Zero]*4 for _ in range(4)]
    for m in range(4):
        for n in range(4):
            s = sp.S.Zero
            for a in range(4):
                s += sp.S(1 if a == 0 else -1) * e[a][m] * e[a][n]
            g[m][n] = sp.simplify(s)
    return g

# ---------------------------------------------------------------------------
# 1) ХОРОШАЯ тетрада (симметричная, экваториальный поворот)
# ---------------------------------------------------------------------------
def good_tetrad():
    A, B = Aval, Bval
    return [
        [A, 0, 0, 0],
        [0, B*sp.sin(th)*sp.cos(ph), r*sp.cos(th)*sp.cos(ph), -r*sp.sin(th)*sp.sin(ph)],
        [0, B*sp.sin(th)*sp.sin(ph), r*sp.cos(th)*sp.sin(ph),  r*sp.sin(th)*sp.cos(ph)],
        [0, B*sp.cos(th),           -r*sp.sin(th),              0],
    ]

# ---------------------------------------------------------------------------
# 2) ДИАГОНАЛЬНАЯ тетрада (как в R1 главной модели)
# ---------------------------------------------------------------------------
def diag_tetrad():
    A, B = Aval, Bval
    return [
        [A, 0, 0, 0],
        [0, B, 0, 0],
        [0, 0, r, 0],
        [0, 0, 0, r*sp.sin(th)],
    ]

# ---------------------------------------------------------------------------
# Кручение / конторсия / суперпотенциал / скаляр T  (независимая реализация)
# ---------------------------------------------------------------------------
def torsion_scalar(e):
    E = [[sp.Integer(1)/e[a][a] if e[a][a] != 0 else sp.S.Zero for a in range(4)]]
    Einv = [[sp.Integer(1)/e[a][a] if e[a][a] != 0 else sp.S.Zero for _ in range(4)] for a in range(4)]
    g = metric_from_tetrad(e)
    gin = [[sp.S.Zero]*4 for _ in range(4)]
    for m in range(4):
        for n in range(4):
            if g[m][n] != 0:
                gin[m][n] = 1/g[m][n]
    # частные производные e^a_mu
    def de(a, mu, x):
        return sp.diff(e[a][mu], x)
    xs = [r, th, ph]
    # T^rho_{mu nu} = e_a^rho (d_mu e^a_nu - d_nu e^a_mu)
    def Tcomp(rho, mu, nu):
        s = sp.S.Zero
        for a in range(4):
            if e[a][mu] == 0 and e[a][nu] == 0:
                continue
            d1 = sp.S.Zero if e[a][nu] == 0 else de(a, nu, r if mu == 1 else th if mu == 2 else ph)
            d2 = sp.S.Zero if e[a][mu] == 0 else de(a, mu, r if nu == 1 else th if nu == 2 else ph)
            s += Einv[rho][a] * (d1 - d2)
        return s
    def T_low(sg, mu, nu):
        if g[sg][sg] == 0:
            return sp.S.Zero
        return g[sg][sg] * Tcomp(sg, mu, nu)
    def K(mu, nu, rho):
        return sp.Rational(1,2)*(T_low(mu,nu,rho)+T_low(rho,mu,nu)-T_low(nu,mu,rho))
    def Tsv(nuu):
        s = sp.S.Zero
        for si in range(4):
            for be in range(4):
                if gin[nuu][be] != 0:
                    s += gin[nuu][be]*Tcomp(si, be, si)
        return s
    def Tsm(muu):
        s = sp.S.Zero
        for si in range(4):
            for be in range(4):
                if gin[muu][be] != 0:
                    s += gin[muu][be]*Tcomp(si, be, si)
        return s
    def S(muu, nuu, rhol):
        Kup = sp.S.Zero
        for a in range(4):
            for b in range(4):
                if gin[muu][a] != 0 and gin[nuu][b] != 0:
                    Kup += gin[muu][a]*gin[nuu][b]*K(a, b, rhol)
        return sp.Rational(1,2)*(Kup
                                 + sp.KroneckerDelta(muu, rhol)*Tsv(nuu)
                                 - sp.KroneckerDelta(nuu, rhol)*Tsm(muu))
    Ts = sp.S.Zero
    for rho_ in range(4):
        for mu_ in range(4):
            for nu_ in range(4):
                Ts += S(mu_, nu_, rho_)*Tcomp(rho_, mu_, nu_)
    return sp.trigsimp(sp.simplify(Ts))

print("="*72)
print("  R2 [ПОМ] — torsion scalar: good vs diagonal tetrad (fresh code)")
print("="*72)

for name, build in [("GOOD (symmetric)", good_tetrad), ("DIAG (bad? R1)", diag_tetrad)]:
    t1 = time.time()
    e = build()
    g = metric_from_tetrad(e)
    diag_ok = all(g[m][n] == 0 for m in range(4) for n in range(4) if m != n)
    print(f"\n  --- {name}: {time.time()-t1:.1f}s ---")
    print("  metric diag? ", diag_ok)
    if diag_ok:
        print("    g_tt=%s  g_rr=%s" % (sp.simplify(g[0][0]), sp.simplify(g[1][1])))
        print("    g_thth=%s  g_phph=%s" % (sp.simplify(g[2][2]), sp.simplify(g[3][3])))
    t2 = time.time()
    T = torsion_scalar(e)
    Ts = sp.simplify(T)
    print("  T(r) size=%d chars, time %.1fs" % (len(str(Ts)), time.time()-t2))
    print("  T(r) =", Ts)
    # быстрая численная сверка: значение в двух точках r при th=pi/4, ph=pi/4
    import mpmath
    for rv in [0.25, 0.5, 1.0]:
        val = sp.N(Ts.subs({th: sp.pi/4, ph: sp.pi/4, r: rv, b0: sp.Rational(1,10)}), 20)
        print("    r=%.2f: T = %s" % (rv, val))

print("\n  TIME %.1f s" % (time.time()-t0))