"""
r3_numeric_ricci.py [ТР3] — полностью численный независимый расчёт Риччи/Эйнштейна
для метрики Морриса–Торн ds² = dt² − B(r)²dr² − r²dΩ², B=(1−b0²/r²)^(−1/2),
метод: конечные разности символов Кристоффеля, mpmath 50 знаков.
Версия 2: Γ^l_{mn} = ½ g^{ll}(∂_m g_{ln} + ∂_n g_{lm} − ∂_l g_{mn}) для
диагональной метрики (в сумме только p=l); ковариантная Ricci собрана явно.
"""

import mpmath as mp
mp.mp.dps = 50

b0 = mp.mpf('0.1')
r0 = mp.mpf('0.2')
th0 = mp.pi / 4
H = mp.mpf('1e-4')

def Bv(r):
    return 1 / mp.sqrt(1 - b0**2 / r**2)

def Bpv(r):
    # B' = d/dr (1-b0^2/r^2)^(-1/2) = -b0^2/r^3 * (1-b0^2/r^2)^(-3/2)
    return -b0**2 / r**3 * (1 - b0**2 / r**2) ** mp.mpf('-1.5')

def make_ch(r, th):
    """ch(l,m,n) = Γ^l_{mn} в точке (r,th)."""
    B_ = Bv(r)
    gd = [mp.mpf(1), -B_**2, -r**2, -r**2 * mp.sin(th)**2]
    gi = [1 / gd[0], 1 / gd[1], 1 / gd[2], 1 / gd[3]]

    def dg(m, n, k):
        if m != n:
            return 0
        if k == 1:
            return {0: 0, 1: -2 * B_ * Bpv(r), 2: -2 * r,
                    3: -2 * r * mp.sin(th)**2}[m]
        if k == 2:
            return {0: 0, 1: 0, 2: 0,
                    3: -2 * r**2 * mp.sin(th) * mp.cos(th)}[m]
        return 0

    def ch(l, m, n):
        # диагональная метрика: вклад только p=l
        s = gi[l] * (dg(l, n, m) + dg(l, m, n) - dg(m, n, l))
        return s / 2
    return ch

def D(f, k):
    """Численная ∂_k функции f(r,th) в (r0,th0). k=1->r, k=2->th."""
    if k == 1:
        return (f(r0 + H, th0) - f(r0 - H, th0)) / (2 * H)
    if k == 2:
        return (f(r0, th0 + H) - f(r0, th0 - H)) / (2 * H)
    return mp.mpf(0)

def ricci(mu, nu):
    ch = make_ch(r0, th0)
    S1 = 0; S2 = 0; S3 = 0; S4 = 0
    for l in range(4):
        S1 += D(lambda r_, t_, l_=l: make_ch(r_, t_)(l_, mu, nu), l)
        for p in range(4):
            S3 += ch(l, l, p) * ch(p, mu, nu)
            S4 -= ch(p, mu, l) * ch(l, nu, p)
    # член −∂_ν Γ^λ_{λμ}: ОДИН раз (вне цикла по l!)
    S2 -= D(lambda r_, t_, mu_=mu: sum(make_ch(r_, t_)(l, mu_, l)
                                       for l in range(4)), nu)
    return S1 + S2 + S3 + S4

B = Bv(r0)
R = (ricci(0, 0) - ricci(1, 1) / B**2
     - ricci(2, 2) / r0**2 - ricci(3, 3) / (r0**2 * mp.sin(th0)**2))

print("mpmath  Ric_00 =", ricci(0, 0))
print("mpmath  Ric_11 =", ricci(1, 1))
print("mpmath  Ric_22 =", ricci(2, 2))
print("mpmath  Ric_33 =", ricci(3, 3))
print("mpmath  R      =", R)
for (m, n) in [(0, 0), (1, 1), (2, 2)]:
    gv = {0: 1, 1: -B**2, 2: -r0**2, 3: -r0**2 * mp.sin(th0)**2}[m]
    print("mpmath  G_%d%d   = %s" % (m, n, ricci(m, n) - mp.mpf(1) / 2 * gv * R))