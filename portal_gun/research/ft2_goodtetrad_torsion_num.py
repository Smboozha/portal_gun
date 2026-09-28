"""
R2 [ПОМ] — torsion scalar T(r): good vs diagonal tetrad (NUMERIC, fast)
========================================================================
Считаем T = S^{nu rho}_mu * T^mu_{nu rho} численно в точке (r, th, ph).
Независимая реализация (без тяжелого sympy-упрощения тригонометрии).

Morris-Thorne: A=1, B=1/sqrt(1-b0^2/r^2), т.е. b(r)=b0^2/r.
"""
import mpmath as mp

mp.mp.dps = 15
pi = mp.pi

# --- координаты: 0=t,1=r,2=th,3=ph ---
# аналитическая производная по координате
def dB_dr(r, b0):
    return (b0**2/r**3) * (1 - b0**2/r**2) ** mp.mpf('-1.5')

def T_value(rv, thv, phv, b0, good=True):
    if good:
        S, C, sp_, cp = mp.sin(thv), mp.cos(thv), mp.sin(phv), mp.cos(phv)
        A = mp.mpf(1); B = 1/mp.sqrt(1 - b0**2/rv**2)
        # e[a][mu]
        e = [
            [A, mp.mpf(0), mp.mpf(0), mp.mpf(0)],
            [0, B*S*cp, rv*C*cp, -rv*S*sp_],
            [0, B*S*sp_, rv*C*sp_, rv*S*cp],
            [0, B*C, -rv*S, mp.mpf(0)],
        ]
    else:
        A = mp.mpf(1); B = 1/mp.sqrt(1 - b0**2/rv**2)
        S = mp.sin(thv)
        e = [
            [A, mp.mpf(0), mp.mpf(0), mp.mpf(0)],
            [0, B, mp.mpf(0), mp.mpf(0)],
            [0, mp.mpf(0), rv, mp.mpf(0)],
            [0, mp.mpf(0), mp.mpf(0), rv*S],
        ]

    # de[a][mu][x] по x: 0=t,1=r,2=th,3=ph
    de = [[[mp.mpf(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    if good:
        S, C, sp_, cp = mp.sin(thv), mp.cos(thv), mp.sin(phv), mp.cos(phv)
        B = 1/mp.sqrt(1 - b0**2/rv**2); dB = dB_dr(rv, b0)
        # производные по r: только там, где B/r
        de[1][1][1] = dB*S*cp;   de[2][1][1] = dB*S*sp_;   de[3][1][1] = dB*C
        de[1][2][1] = C*cp;      de[2][2][1] = C*sp_;      de[3][2][1] = -S
        de[1][3][1] = -S*sp_;    de[2][3][1] = S*cp
        # производные по theta
        de[1][1][2] = B*C*cp;    de[2][1][2] = B*C*sp_;    de[3][1][2] = -B*S
        de[1][2][2] = -rv*S*cp;  de[2][2][2] = -rv*S*sp_;  de[3][2][2] = -rv*C
        de[1][3][2] = -rv*C*sp_; de[2][3][2] = rv*C*cp
        # производные по phi
        de[1][1][3] = -B*S*sp_;  de[2][1][3] = B*S*cp
        de[1][2][3] = -rv*C*sp_; de[2][2][3] = rv*C*cp
        de[1][3][3] = -rv*S*cp;  de[2][3][3] = -rv*S*sp_
    else:
        A = mp.mpf(1); B = 1/mp.sqrt(1 - b0**2/rv**2); dB = dB_dr(rv, b0)
        S, C = mp.sin(thv), mp.cos(thv)
        de[1][1][1] = dB
        de[3][3][2] = rv*C

    # обратная тетрада E = inverse(e)  (e — матрица 4x4, строки = Lorentz-индекс a)
    Efull = mp.matrix([[e[a][m] for m in range(4)] for a in range(4)])
    Einv_full = mp.inverse(Efull)
    Einv = [[Einv_full[a, m] for m in range(4)] for a in range(4)]

    # метрика
    g = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    for m in range(4):
        for n in range(4):
            s = mp.mpf(0)
            for a in range(4):
                s += (mp.mpf(1) if a == 0 else mp.mpf(-1)) * e[a][m]*e[a][n]
            g[m][n] = s
    gin = [[mp.mpf(0) if i != j else (1/g[i][i] if g[i][i] != 0 else mp.mpf(0)) for j in range(4)] for i in range(4)]

    def Tcomp(rho, mu, nu):
        s = mp.mpf(0)
        for a in range(4):
            s += Einv[rho][a] * (de[a][nu][mu] - de[a][mu][nu])
        return s

    nz = lambda x: (x if not isinstance(x, int) else mp.mpf(x))
    def T_low(sg, mu, nu):
        return nz(g[sg][sg]) * Tcomp(sg, mu, nu)

    def K(mu, nu, rho):
        return mp.mpf('0.5')*(T_low(mu,nu,rho)+T_low(rho,mu,nu)-T_low(nu,mu,rho))

    def Tsv(nuu):
        s = mp.mpf(0)
        for si in range(4):
            for be in range(4):
                s += nz(gin[nuu][be]) * Tcomp(si, be, si)
        return s
    def Tsm(muu):
        s = mp.mpf(0)
        for si in range(4):
            for be in range(4):
                s += nz(gin[muu][be]) * Tcomp(si, be, si)
        return s

    def S(muu, nuu, rhol):
        Kup = mp.mpf(0)
        for a in range(4):
            for b in range(4):
                Kup += nz(gin[muu][a])*nz(gin[nuu][b])*K(a, b, rhol)
        return mp.mpf('0.5')*(Kup + (mp.mpf(1) if muu == rhol else mp.mpf(0))*Tsv(nuu)
                              - (mp.mpf(1) if nuu == rhol else mp.mpf(0))*Tsm(muu))

    Ts = mp.mpf(0)
    for rho_ in range(4):
        for mu_ in range(4):
            for nu_ in range(4):
                Ts += S(mu_, nu_, rho_)*Tcomp(rho_, mu_, nu_)
    return Ts

b0 = mp.mpf('0.1')
print("="*72)
print("  R2 [ПОМ] — T(r) численно: GOOD vs DIAGONAL (Morris-Thorne, b0=0.1)")
print("="*72)

for rv in [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0')]:
    row = []
    for (thv, phv) in [(pi/4, pi/4), (pi/3, pi/6), (pi/2, 0)]:
        tg = T_value(rv, thv, phv, b0, good=True)
        row.append(tg)
    same = all(abs(row[0]-x) < 1e-9 for x in row)
    print("  GOOD  r=%.2f: T=%.12e (th/ph-variation same? %s)" % (rv, row[0], same))

print()
for rv in [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0')]:
    td = T_value(rv, pi/4, pi/4, b0, good=False)
    print("  DIAG  r=%.2f: T=%.12e" % (rv, td))

# сравнение: DIAG должна совпадать с символьным T_scalar из R1 (там b0^2/r^2)
print()
print("  Ожидаемая форма GOOD (по литературе T = 2/r^2*(1 - 1/B^2) и т.п.):")
for rv in [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0')]:
    tg = T_value(rv, pi/4, pi/4, b0, good=True)
    B = 1/mp.sqrt(1 - b0**2/rv**2)
    cand1 = (mp.mpf(2)/rv**2)*(1 - 1/B**2)
    cand2 = -(mp.mpf(2)/rv**2)*(1 - 1/B**2)
    print("    r=%.2f: T_good=%.10e | 2/r^2*(1-1/B^2)=%.6e | -то же=%.6e" % (rv, tg, cand1, cand2))