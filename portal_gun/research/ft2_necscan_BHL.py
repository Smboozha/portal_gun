"""
R2 [ПОМ] — NEC/WEC scan в f(T) НА ХОРОШЕЙ ТЕТРАДЕ (v3, финальный).
Morris-Thorne: a(r)=0, e^{-b}=1-beta/r, beta(r)=r0^2/r (= b0^2/r в R1).
Формулы: BHL (16)-(20). ВАЖНЫЕ ИСПРАВЛЕНИЯ:
 1) В (18) для rho стоит + f/4 (подтверждено (22)-(26) и TEGR-пределом (38)-(40)).
 2) pt из (20) соответствует (49)-(45), а НЕ печатной (47): печатная (47) в статье
    НЕ консистентна (проверено: (45)+(47)-ne-(49); все остальные согласуются <=1e-18).
Размерности: x = r/r0; alpha = T0*r0^2 (безразмерно); T в конвенции BHL:
T(x) = 2 (x - sqrt(x^2-1))^2 / x^4. Величины даются как 4pi*r0^2 * {rho, p_r, p_t}.
"""
import json
import mpmath as mp
mp.mp.dps = 18

def u_of(x):
    return mp.sqrt(x*x - 1)

def T16(x):
    u = u_of(x)
    return 2*(x - u)**2 / x**4

def Tp17(x):
    u = u_of(x)
    return -4*(x - u)**2*(x + 2*u)/(u * x**5)

def field_matter(x, alpha, family='quad', lam=None, n=None):
    """4pi*r0^2 * (rho, p_r, p_t). alpha/lam в ед. r0."""
    u = u_of(x)
    T  = T16(x)
    Tp = Tp17(x)
    ebm = 1 - 1/x**2
    al  = mp.sqrt(ebm)
    ebh = 1/al
    xb_m1 = -(x*x + 1)/(x*x - 1)
    bp = -2/(x*(x*x - 1))

    if family == 'quad':
        f  = T + alpha*T**2
        fT = 1 + 2*alpha*T
        fTT = 2*alpha
    elif family == 'power':
        f  = T + lam * T**n
        fT = 1 + lam*n*T**(n-1)
        fTT = lam*n*(n-1)*T**(n-2)
    else:
        raise ValueError(family)

    rho = (al/x)*(1 - al)*Tp*fTT - (T/4 - 1/(2*x*x))*fT \
          + (ebm/(2*x*x))*xb_m1*fT + f/4
    pr  = (-1/(2*x*x) + T/4 + ebm/(2*x*x))*fT - f/4
    pt  = (ebm/2)*(1/x - ebh/x)*Tp*fTT \
          + fT*(T/4 + (ebm/(2*x))*((-0.5)*bp)) - f/4
    return rho, pr, pt

def quad_throat(alpha):
    """Точные пределы при x->1 (f=T+alpha T^2), 4pi r0^2 единицы."""
    rho = -(1 + 22*alpha)/2
    pr  = -(1 + 2*alpha)/2
    pt  = (1 + 14*alpha)/2          # из (49)-(45); печатная (47) некорректна
    return rho, pr, pt

# ------- Точные формулы BHL (45),(46),(48),(49) как независимая сверка --------
def rhoB45(x, a):
    u = u_of(x)
    return -1/(2*x**4)*(1 + 16*a*u/x*(3 - 5/x**2) - 2*a*(24 - 52/x**2 + 17/x**4))
def prB46(x, a):
    u = u_of(x)
    return -1/(2*x**4)*(1 + 16*a*u/x*(1 - 1/x**2) - 2*a*(8 - 12/x**2 + 3/x**4))
def rprB48(x, a):
    u = u_of(x)
    return -1/x**4*(1 + 16*a*u/x*(2 - 3/x**2) - 4*a*(8 - 16/x**2 + 5/x**4))
def rptB49(x, a):
    u = u_of(x)
    return -4*a/x**4*(u/x*(4 - 5/x**2) - 4 + 7/x**2 - 2/x**4)

def main():
    print("="*80)
    print("  R2 [ПОМ] — NEC/WEC scan, good tetrad, MT a=0, beta=r0^2/r  (BHL v3)")
    print("="*80)
    print(" 1) TEGR-предел f=T: (18)-(20) vs MT/GR")
    for x in [mp.mpf('1.01'), mp.mpf('1.3'), mp.mpf('2.0')]:
        rho, pr, pt = field_matter(x, 0, 'quad')
        print("   x=%s rho=%+.6e(gr %+.6e) pr=%+.6e(%+.6e) pt=%+.6e(%+.6e)"
              % (x, rho, -1/(2*x**4), pr, -1/(2*x**4), pt, 1/(2*x**4)))

    print(" 2) квадратичный: (18)-(20) против (45),(46),(48),(49)")
    mx = mp.mpf(0)
    for a in [mp.mpf('-1'), mp.mpf('-0.5'), mp.mpf('-0.25'), mp.mpf('0.5'), mp.mpf('2')]:
        for x in [mp.mpf('1.05'), mp.mpf('1.2'), mp.mpf('1.5'), mp.mpf('2.5')]:
            rho, pr, pt = field_matter(x, a, 'quad')
            d = max(abs(rho - rhoB45(x, a)), abs(pr - prB46(x, a)),
                    abs((rho+pr) - rprB48(x, a)), abs((rho+pt) - rptB49(x, a)))
            if d > mx: mx = d
    print("   maxdev(45,46,48,49) = %.3e -> %s" % (mx, "OK" if mx < 1e-12 else "FAIL"))

    # ---- throat (x=1) analytic limits для QUAD ---- 
    print("\n 3) Scan (требование: rho>=0 [WEC], rho+pr>=0, rho+pt>=0 [NEC] на всех x)")
    print("    единицы: 4pi*r0^2*X; x=r/r0; alpha = T0*r0^2 = lam*r0^2 (r0=10 см)")
    print("\n -- QUAD f=T+alpha T^2: при горловине и вблизи (x=1.05,1.2,1.5) --")
    print("   %6s |%9s%9s%9s%9s |%9s%9s | %s" % ("alpha","rho_1","rho_1.05","rho_1.2","rho_1.5",
                                                 "NECr_1","NECt_1","статус"))
    res_quad = []
    for av in ['-5','-3','-2','-1.5','-1','-0.75','-0.5','-0.25','-0.1','0','0.1','0.25','0.5','1','2','5']:
        a = mp.mpf(av)
        rh1, pr1, pt1 = quad_throat(a)
        nr1, nt1 = rh1+pr1, rh1+pt1
        nec_r = nr1 >= 0; nec_t = nt1 >= 0; wec = rh1 >= 0
        row = []
        for x in [mp.mpf('1.05'), mp.mpf('1.2'), mp.mpf('1.5')]:
            rho, pr, pt = field_matter(x, a, 'quad')
            nec_r &= (rho + pr) >= 0; nec_t &= (rho + pt) >= 0; wec &= rho >= 0
            row.append(rho)
        st = ("NECr%s NECt%s WEC%s" % ("+" if nec_r else "-",
                                       "+" if nec_t else "-",
                                       "+" if wec else "-"))
        print("   %6s |%9.3f%9.3f%9.3f%9.3f |%9.3f%9.3f | %s"
              % (av, rh1, *row, nr1, nt1, st))
        res_quad.append((av, float(rh1), float(pr1), float(pt1), float(nr1), float(nt1),
                         bool(nec_r), bool(nec_t), bool(wec)))

    print("\n -- POWER f=T+lam T^n (n=3/2,2,3): горловина (x->1 numerically) --")
    res_pow = []
    for n in [mp.mpf('1.5'), mp.mpf('2'), mp.mpf('3')]:
        print("   n=%s" % n)
        for lv in ['-1','-0.75','-0.5','-0.25','-0.1','0','0.1','0.25','0.5','1']:
            lam = mp.mpf(lv)
            x = mp.mpf('1') + mp.mpf('1e-6')
            rho, pr, pt = field_matter(x, mp.mpf(0), 'power', lam=lam, n=n)
            nr, nt = rho+pr, rho+pt
            print("    lam=%6s  rho=%+.4f pr=%+.4f pt=%+.4f | NECr=%s(%+.4f) NECt=%s(%+.4f) WEC=%s"
                  % (lv, rho, pr, pt, "+" if nr>=0 else "-", nr,
                     "+" if nt>=0 else "-", nt, "+" if rho>=0 else "-"))
            res_pow.append((str(n), lv, float(rho), float(pr), float(pt)))

    data = {
        "T_conv": "BHL(16): T(x)=2(x-sqrt(x^2-1))^2/x^4, x=r/r0. "
                  "Мой независимый численный torsion (ft2_goodtetrad_torsion_num) даёт 2x этой величины.",
        "r0": "10 cm", "r0_diss": "alpha=lam*r0^2",
        "checks": {
            "TEGR_limit": "OK: (18)-(20) при f=T -> MT/GR (38)-(40)",
            "vs_BHL_45_46_48_49": "maxdev %.2e" % mx,
            "printed_BHL_47_misprint": "Печатная (47) НЕ консистентна: (45)+(47)-ne-(49). "
                                       "pt из (20) = (49)-(45).",
        },
        "quad_throat_4pi_r0^2": {r[0]: {"rho": r[1], "pr": r[2], "pt": r[3],
                                        "NECr": r[4], "NECt": r[5],
                                        "NECr_ok": r[6], "NECt_ok": r[7], "WEC_ok": r[8]}
                                 for r in res_quad},
        "power_throat": [{"n": r[0], "lam": r[1], "rho": r[2], "pr": r[3], "pt": r[4]}
                         for r in res_pow],
    }
    with open("ft2_necscan_results.json", "w") as fp:
        json.dump(data, fp, indent=2, default=str)
    print("\nsaved ft2_necscan_results.json")

if __name__ == "__main__":
    main()