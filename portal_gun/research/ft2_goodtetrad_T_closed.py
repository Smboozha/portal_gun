"""
R2 [ПОМ] — закрытая форма T_good(r) для симметричной (good) тетрады.
Подставляем углы ЧИСЛАМИ (th=pi/4, ph=pi/4) до алгебры -> тригонометрия исчезает,
упрощение символьное быстрое. Проверка: результат не должен зависеть от выбора углов.
"""
import time
import sympy as sp

t0 = time.time()
r, b0 = sp.symbols('r b0', positive=True)
r2 = sp.Integer(1)/sp.sqrt(2)  # sin(pi/4)=cos(pi/4)

def T_symbolic_for_angles(th, ph):
    S, C, sp_, cp = sp.sin(th), sp.cos(th), sp.sin(ph), sp.cos(ph)
    A = sp.Integer(1)
    B = 1/sp.sqrt(1 - b0**2/r**2)
    e = [
        [A, 0, 0, 0],
        [0, B*S*cp, r*C*cp, -r*S*sp_],
        [0, B*S*sp_, r*C*sp_, r*S*cp],
        [0, B*C, -r*S, 0],
    ]
    Einv = sp.Matrix(e).inv()
    g = sp.zeros(4)
    for m in range(4):
        for n in range(4):
            s = 0
            for a in range(4):
                s += (1 if a == 0 else -1)*e[a][m]*e[a][n]
            g[m, n] = sp.simplify(s)
    gin = g.inv()

    def de(a, mu, x):
        if x == 1:
            return sp.diff(e[a][mu], r)
        return sp.Integer(0)

    def Tcomp(rho, mu, nu):
        s = 0
        for a in range(4):
            s += Einv[rho, a]*(de(a, nu, mu) - de(a, mu, nu))
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
    return sp.simplify(T)

print("=" * 72)
print("  R2 [ПОМ] — closed form T_good(r),  th=pi/4, ph=pi/4")
print("=" * 72)

T1 = T_symbolic_for_angles(sp.pi/4, sp.pi/4)
print("  time %.1f s" % (time.time()-t0))
print("  T_good(r) =", T1)
try:
    Tf = sp.factor(T1)
    print("  factor    =", Tf)
except Exception as ex:
    print("  factor fail:", ex)

# Проверка угловой независимости численно (разные углы -> те же числа)
import mpmath as mp
for rv in [mp.mpf('0.5'), mp.mpf('1.0')]:
    b = mp.mpf('0.1')
    v1 = sp.N(T1.subs({r: rv, b0: b}), 25)
    # пересчёт чистым численным кодом с другими углами
    exec(open('ft2_goodtetrad_torsion_num.py').read().split('b0 = mp.mpf')[0])
    v2 = T_value(rv, mp.pi/3, mp.pi/6, b, good=True)
    print("  r=%.2f : symbolic=%s  numeric(th=pi/3) =%.15e" % (rv, v1, v2))
print("  TIME %.1f s" % (time.time()-t0))