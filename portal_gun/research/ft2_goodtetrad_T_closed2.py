"""
R2 [ПОМ] — закрытая форма T_good(r) для симметричной тетрады, ЧЕСТНО:
дифференцирование по символическим theta, phi, ЗАТЕМ числовые углы, потом simplify.
"""
import time
import sympy as sp

t0 = time.time()
r, th, ph, b0 = sp.symbols('r th ph b0', positive=True)

A = sp.Integer(1)
B = 1/sp.sqrt(1 - b0**2/r**2)
e = [
    [A, 0, 0, 0],
    [0, B*sp.sin(th)*sp.cos(ph), r*sp.cos(th)*sp.cos(ph), -r*sp.sin(th)*sp.sin(ph)],
    [0, B*sp.sin(th)*sp.sin(ph), r*sp.cos(th)*sp.sin(ph),  r*sp.sin(th)*sp.cos(ph)],
    [0, B*sp.cos(th),           -r*sp.sin(th),              0],
]

def de(a, mu):
    return [sp.diff(e[a][mu], r), sp.diff(e[a][mu], th), sp.diff(e[a][mu], ph)]

def D(a, name, coord):
    if coord == 1 or coord == 0 and name == 0:
        pass
    if coord == 0:
        return sp.Integer(0)
    idx = 0 if coord == 1 else 1 if coord == 2 else 2
    return de(a, name)[idx]

Einv = sp.Matrix(e).inv()
g = sp.zeros(4)
for m in range(4):
    for n in range(4):
        s = 0
        for a in range(4):
            s += (1 if a == 0 else -1)*e[a][m]*e[a][n]
        g[m, n] = s
gin = g.inv()

def Tcomp(rho, mu, nu):
    s = 0
    for a in range(4):
        s += Einv[rho, a]*(D(a, nu, mu) - D(a, mu, nu))
    return s

def Tl(sg, mu, nu):
    return g[sg, sg]*Tcomp(sg, mu, nu)

def K(mu, nu, rho):
    return sp.Rational(1, 2)*(Tl(mu, nu, rho) + Tl(rho, mu, nu) - Tl(nu, mu, rho))

def Tsv(nuu):
    s = 0
    for si in range(4):
        for be in range(4):
            s += gin[nuu, be]*Tcomp(si, be, si)
    return s

def Tsm(muu):
    s = 0
    for si in range(4):
        for be in range(4):
            s += gin[muu, be]*Tcomp(si, be, si)
    return s

def Sup(muu, nuu, rhol):
    s = 0
    for a in range(4):
        for b in range(4):
            s += gin[muu, a]*gin[nuu, b]*K(a, b, rhol)
    return sp.Rational(1, 2)*(s
                              + sp.KroneckerDelta(muu, rhol)*Tsv(nuu)
                              - sp.KroneckerDelta(nuu, rhol)*Tsm(muu))

T = 0
for rh in range(4):
    for mu in range(4):
        for nu in range(4):
            T += Sup(mu, nu, rh)*Tcomp(rh, mu, nu)
print("T построен (символьный, th/ph живы): %.1f s" % (time.time()-t0))

Th = sp.expand(T.subs({th: sp.pi/4, ph: sp.pi/4}))
print("после подстановки углов: %.1f s, size %d" % (time.time()-t0, len(str(Th))))
Ts = sp.simplify(Th)
print("simplify: %.1f s" % (time.time()-t0))
Ts = sp.factor(Ts)
print("factor: %.1f s" % (time.time()-t0))
print("T_good(r) =", Ts)
print("TIME %.1f s" % (time.time()-t0))