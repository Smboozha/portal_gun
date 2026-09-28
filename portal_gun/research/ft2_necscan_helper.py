"""
R2 [ПОМ] — NEC/WEC scan для f(T)-червоточины НА ХОРОШЕЙ ТЕТРАДЕ.
============================================================================
Контекст:
  - R1 главной модели построена на ДИАГОНАЛЬНОЙ (bad) тетраде: T=1/r^2,
    не зависит от формы горловины -> f(T)-динамика там не видит червоточину.
  - Здесь: GOOD (symmetric) тетрада, T_good(r)=4(r-sqrt(r^2-b0^2))^2/r^4.
  Численные функциональные производные действия S=∫ A·B·C²·(T+f(T)) dr дают
  уравнения поля; калибровка по GR при f->0; читаем (rho, p_r, p_t) и NEC/WEC.
  - Morris-Thorne: A=1, B=1/sqrt(1-b0^2/r^2), C=r. b0=10 см.
"""
import json, time
import mpmath as mp

mp.mp.dps = 25
pi = mp.pi

# ---------------- параметризованный торсионный скаляр (good tetrad) ----------
def torsion_at(r, b0, A, dA, B, dB, C, dC=1):
    """T в точке (r, pi/4, pi/4) для good-тетрады, профиль (A,dA,B,dB,C)."""
    thv = pi/4; phv = pi/4
    S, Cp, sp_, cp = mp.sin(thv), mp.cos(thv), mp.sin(phv), mp.cos(phv)
    e = [
        [A, 0, 0, 0],
        [0, B*S*cp, C*Cp*cp, -C*S*sp_],
        [0, B*S*sp_, C*Cp*sp_, C*S*cp],
        [0, B*Cp, -C*S, 0],
    ]
    de = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    # по r
    de[0][0][1] = dA
    de[1][1][1] = dB*S*cp;  de[2][1][1] = dB*S*sp_;  de[3][1][1] = dB*Cp
    de[1][2][1] = dC*cp*Cp; de[2][2][1] = dC*sp_*Cp; de[3][2][1] = -dC*S
    de[1][3][1] = -dC*S*sp_; de[2][3][1] = dC*S*cp
    # по theta
    de[1][1][2] = B*Cp*cp;  de[2][1][2] = B*Cp*sp_;  de[3][1][2] = -B*S
    de[1][2][2] = -C*S*cp;  de[2][2][2] = -C*S*sp_;  de[3][2][2] = -C*Cp
    de[1][3][2] = -C*Cp*sp_; de[2][3][2] = C*Cp*cp
    # по phi
    de[1][1][3] = -B*S*sp_; de[2][1][3] = B*S*cp
    de[1][2][3] = -C*Cp*sp_; de[2][2][3] = C*Cp*cp
    de[1][3][3] = -C*S*cp; de[2][3][3] = -C*S*sp_

    Efull = mp.matrix([[e[a][m] for m in range(4)] for a in range(4)])
    Einv = mp.inverse(Efull)
    g = [[0]*4 for _ in range(4)]
    for m in range(4):
        for n in range(4):
            s = 0
            for a in range(4):
                s += (1 if a == 0 else -1)*e[a][m]*e[a][n]
            g[m][n] = s
    gin = [[1/g[i][i] if i == j else 0 for j in range(4)] for i in range(4)]

    def Tcomp(rho, mu, nu):
        s = 0
        for a in range(4):
            s += Einv[rho, a]*(de[a][nu][mu] - de[a][mu][nu])
        return s

    def Tl(sg, mu, nu):
        return g[sg][sg]*Tcomp(sg, mu, nu)

    def Tsv(nuu):
        s = 0
        for si in range(4):
            for be in range(4):
                s += gin[nuu][be]*Tcomp(si, be, si)
        return s

    Ts = 0
    for rh in range(4):
        for mu in range(4):
            for nu in range(4):
                # ко нторция K_{mu nu rho} = 1/2(T_munu_rho + T_rhonu_mu? - T_nu_mu_rho)
                Kup = 0
                for a in range(4):
                    for b in range(4):
                        Kdown = mp.mpf('0.5')*(Tl(mu,nu,a)+Tl(a,mu,nu)-Tl(nu,mu,a))
                        Kup += gin[mu][a]*gin[nu][b]*Kdown
                Sup = mp.mpf('0.5')*(Kup + (1 if mu == rh else 0)*Tsv(nu)
                                     - (1 if nu == rh else 0)*Tsv(mu))
                Ts += Sup*Tcomp(rh, mu, nu)
    return Ts

# ---------------- Morris-Thorne ----------------------------------------------
def Bmet(r, b0): return 1/mp.sqrt(1 - b0*b0/(r*r))
def dBmet(r, b0): return (b0*b0/(r*r*r))*(1 - b0*b0/(r*r))**mp.mpf('-1.5')
def T_MT(r, b0):
    u = mp.sqrt(r*r - b0*b0)
    return 4*(r - u)**2/(r**4)

b0 = mp.mpf('0.10')
print("check torsion_at vs closed form T_MT:")
for rv in [mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0')]:
    ta = torsion_at(rv, b0, 1, 0, Bmet(rv, b0), dBmet(rv, b0), rv)
    tc = T_MT(rv, b0)
    print("  r=%s  numeric=%.14e  closed=%.14e  reldev=%+.2e" % (rv, ta, tc, (ta-tc)/tc))

# ---------------- EL (численные функциональные производные) -------------------
EPS = mp.mpf('1e-8'); H = mp.mpf('1e-6')

def Lag(r, b0, A, dA, B, dB, C, f, dC=1):
    T = torsion_at(r, b0, A, dA, B, dB, C, dC)
    return A*B*C**2*(T + f(T))

def Lq(q):  # local Lagrangian for the MT-like profile given by dict q
    return Lag(q['r'], q['b0'], q['A'], q['dA'], q['B'], q['dB'], q['C'], q['f'])

def dL_dq(r, b0, which, f, q):
    eps = EPS
    q0 = dict(q)
    def stim(w, sh):
        qq = dict(q0); qq[w] = qq[w] + sh
        return Lq(qq)
    lo = stim(which, -eps); hi = stim(which, +eps)
    return (hi - lo)/(2*eps)

def dL_dqprime(r, b0, which, f, q):
    eps = EPS
    q0 = dict(q)
    def stim(w, sh):
        qq = dict(q0)
        key = 'dC' if w == 'C' else 'd'+w
        qq[key] = qq[key] + sh
        return Lq(qq)
    lo = stim(which, -eps); hi = stim(which, +eps)
    return (hi - lo)/(2*eps)

def EL_at(r, b0, which, f):
    q = {'r': r, 'b0': b0, 'f': f,
         'A': 1, 'dA': 0, 'B': Bmet(r, b0), 'dB': dBmet(r, b0), 'C': r, 'dC': 1}
    h = H
    dLdq = dL_dq(r, b0, which, f, q)
    def dp(x):
        qx = {'r': x, 'b0': b0, 'f': f,
              'A': 1, 'dA': 0, 'B': Bmet(x, b0), 'dB': dBmet(x, b0), 'C': x, 'dC': 1}
        return dL_dqprime(x, b0, which, f, qx)
    der = (dp(r+h) - dp(r-h))/(2*h)
    return dLdq - der

# --------- калибровка: GR-предел f=0 -> сопоставление EL и (rho,p_r,p_t) ------
# В GR (TEGR), для MT: rho = b'(r)/(2 r^2)? , p_r = -b'(r)/(2 r^2)? Прямая проверка
# по определению: G_tt должен дать rho. Сверяем EL_A с аналитикой.
def rho_gr(r, b0):
    # b(r)=b0^2/r, b' = -2 b0^2/r^3
    return -2*b0*b0/(r**3)/(2*r*r)   # -b' / (2 r^2)? знак уточнить
def p_r_gr(r, b0):
    return -b0*b0/(r**3)              # -b/(2 r^3)? знак уточнить

print("\nCalibration f=0:")
f0 = lambda _: mp.mpf(0)
for rv in [mp.mpf('0.5'), mp.mpf('1.0')]:
    ea = EL_at(rv, b0, 'A', f0)
    eb = EL_at(rv, b0, 'B', f0)
    ec = EL_at(rv, b0, 'C', f0)
    print("  r=%s  EL_A=%+.8e  EL_B=%+.8e  EL_C=%+.8e" % (rv, ea, eb, ec))

if __name__ == "__main__":
    pass