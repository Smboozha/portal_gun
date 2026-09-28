"""
r3_debug_worm_num.py [ТР3] — пошаговый трейс _worm_num(): самодостаточный,
копирует построение EL из wormhole_ft_derive.py и печатает free_symbols
после КАЖДОГО шага подстановки для каждого EL[name].
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
    if g[sigma][sigma] == 0:
        return sp.S.Zero
    return g[sigma][sigma] * Tcomp(sigma, mu, nu)

def K(mu, nu, rho):
    return sp.Rational(1, 2) * (T_low(mu, nu, rho) + T_low(rho, mu, nu)
                                - T_low(nu, mu, rho))

def T_sigma_nu_sigma(nuu):
    s = sp.S.Zero
    for si in range(4):
        for be in range(4):
            if ginv[nuu][be] != 0:
                s += ginv[nuu][be] * Tcomp(si, be, si)
    return s

def T_sigma_mu_sigma(muu):
    s = sp.S.Zero
    for si in range(4):
        for be in range(4):
            if ginv[muu][be] != 0:
                s += ginv[muu][be] * Tcomp(si, be, si)
    return s

def Superpotential(muu, nuu, rhol):
    Kup = sp.S.Zero
    for a in range(4):
        for b in range(4):
            if ginv[muu][a] != 0 and ginv[nuu][b] != 0:
                Kup += ginv[muu][a] * ginv[nuu][b] * K(a, b, rhol)
    return sp.Rational(1, 2) * (Kup
                                + sp.KroneckerDelta(muu, rhol) * T_sigma_nu_sigma(nuu)
                                - sp.KroneckerDelta(nuu, rhol) * T_sigma_mu_sigma(muu))

T_scalar = sp.S.Zero
for rho_ in range(4):
    for mu_ in range(4):
        for nu_ in range(4):
            T_scalar += Superpotential(mu_, nu_, rho_) * Tcomp(rho_, mu_, nu_)
T_scalar = sp.simplify(T_scalar)
print("T_scalar =", T_scalar)

F = sp.Function('f')
L = A * B * C ** 2 * (T_scalar + F(T_scalar))
EL = {}
for name, q in [('A', A), ('B', B), ('C', C)]:
    qp = q.diff(r)
    EL[name] = sp.simplify(sp.diff(L, q) - sp.diff(sp.diff(L, qp), r))
    print("EL[%s] size=%d тыс" % (name, len(str(EL[name])) // 1000))

# G_μν  — те же формулы, что в модуле
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
    s = sp.S.Zero
    for l in range(4):
        s += dd(chris(l, mu, nu), l)
        for sig in range(4):
            s -= dd(chris(l, mu, l), nu)
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

def show(label, e):
    print("   %-42s free=%s" % (label, sorted(e.free_symbols, key=str)))
    print("   %-42s funcs=%s" % ("", sorted({f.func.__name__ for f in e.atoms(sp.Function)})))

rv, bv, thn = 0.2, 0.1, sp.pi/4
db = 1 - bv**2/rv**2
dbp = 2*bv**2/rv**3
dbpp = -6*bv**2/rv**4
Bn = db**-0.5
Bpn = 0.5*db**-1.5*dbp
Bppn = -0.75*db**-2.5*dbp**2 + 0.5*db**-1.5*dbpp

for name in ['A', 'B', 'C']:
    print("#" * 78)
    print("# EL[%s]  (f->0)  пошаговый трейс" % name)
    print("#" * 78)
    e = EL[name]
    e = e.replace(lambda t: isinstance(t, sp.Function)
                  and getattr(t, 'func', None) == sp.Function('f'),
                  lambda t: sp.S.Zero)
    show("исходное", e)
    e = e.subs({r: rv, b0: bv, th: thn})
    show("после r,b0,th", e)
    orig = e
    e = e.subs(A.diff(r, 2), Bppn*0).subs(B.diff(r, 2), Bppn).subs(C.diff(r, 2), 0)
    show("после 2-х производных", e)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bpn).subs(C.diff(r), 1)
    show("после 1-х производных", e)
    e = e.subs(A, 1).subs(B, Bn).subs(C, rv).doit()
    show("после A,B,C + doit", e)
    try:
        print("   EVALF =", e.evalf(30))
    except Exception as ex:
        print("   EVALF FAILED:", type(ex).__name__, ex)

# то же для G(0,0), G(1,1), G(2,2)
print("#" * 78)
print("# G_μν пошаговый трейс (те же подстановки)")
print("#" * 78)
for (m_, n_) in [(0, 0), (1, 1), (2, 2)]:
    print("# G_%d%d" % (m_, n_))
    e = G(m_, n_)
    show("исходное", e)
    e = e.subs({r: rv, b0: bv, th: thn})
    e = e.subs(A.diff(r, 2), Bppn*0).subs(B.diff(r, 2), Bppn).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bpn).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bn).subs(C, rv).doit()
    show("после всех subs", e)
    try:
        print("   EVALF =", e.evalf(30))
    except Exception as ex:
        print("   EVALF FAILED:", type(ex).__name__, ex)