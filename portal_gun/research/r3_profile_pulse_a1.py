"""
R3 — Профили T(r), ρ(r), p_r(r), p_t(r) + импульсный режим (1 мкс)
=====================================================================
f(T) = T + αT² + γT⁻¹,  α = −0.1,  γ = −0.01
MT/Ellis червоточина: a(r)=0, β(r)=r0²/r, good tetrad, r0 = 0.1 м.

Единицы:
  Безразмерные: 4πr0² × {ρ, p_r, p_t};  x = r/r0.
  СИ: ρ_СИ = ρ_dim × c⁴/(8πG·r0²),  p_СИ = p_dim × c⁴/(8πG·r0²)
      c⁴/G = 1.209×10⁴⁴ Н → c⁴/(8πG·r0²) = 1.209e44/(8π·0.01) ≈ 4.812e44 Па
"""
import json
import mpmath as mp

mp.mp.dps = 40

# ── Физические константы ──────────────────────────────────────────────────────
r0 = mp.mpf('0.1')                         # м
c4_over_G = mp.mpf('1.209e44')             # Н = кг·м/с²
c4_over_8piGr0sq = c4_over_G / (8 * mp.pi * r0**2)   # Па
alpha = mp.mpf('-0.1')
gamma = mp.mpf('-0.01')

# ── Скаляр кручения и производные ─────────────────────────────────────────────

def u_of(x):
    return mp.sqrt(x*x - 1)

def T_of(x):
    """T(x) = 2(x − √(x²−1))² / x⁴  — BHL (16)"""
    u = u_of(x)
    return 2*(x - u)**2 / x**4

def Tp_of(x):
    """dT/dx = −4(x−u)²(x+2u)/(u·x⁵)  — BHL (17)"""
    u = u_of(x)
    return -4*(x - u)**2*(x + 2*u)/(u * x**5)

# ── f(T) и производные ────────────────────────────────────────────────────────

def f_funcs(T, a, g):
    """f = T + αT² + γ/T,  fT = 1+2αT−γ/T²,  fTT = 2α+2γ/T³"""
    f  = T + a*T**2 + g/T
    fT = 1 + 2*a*T - g/T**2
    fTT = 2*a + 2*g/T**3
    return f, fT, fTT

# ── Аналитика на горловине x=1 ───────────────────────────────────────────────

def throat_analytic(a, g):
    """eq.(11)–(13):  ρ, p_r, p_t  на горловине."""
    rho = -(1 + 22*a + 24*g) / 2
    pr  = -(1 + 2*a + 4*g) / 2
    pt  = (1 + 14*a + 22*g) / 2
    return rho, pr, pt

# ── Численные формулы (eq.1–eq.3) ───────────────────────────────────────────

def field_matter(x, a, g):
    """4πr0² × (ρ, p_r, p_t)  для f=T+αT²+γT⁻¹ на MT/Ellis, good tetrad."""
    if abs(x - 1) < 1e-15:
        return throat_analytic(a, g)

    u  = u_of(x)
    T  = T_of(x)
    Tp = Tp_of(x)
    ebm = 1 - 1/(x*x)
    al  = mp.sqrt(ebm)
    xb_m1 = -(x*x + 1)/(x*x - 1)
    bp = -2/(x*(x*x - 1))
    ebh = 1/al

    f, fT, fTT = f_funcs(T, a, g)

    rho = (al/x)*(1 - al)*Tp*fTT \
          - (T/4 - 1/(2*x*x))*fT \
          + (ebm/(2*x*x))*xb_m1*fT \
          + f/4

    pr  = (-1/(2*x*x) + T/4 + ebm/(2*x*x))*fT - f/4

    ebm_bracket = ebm*(1/x - ebh/x)
    pt  = (ebm_bracket/2)*Tp*fTT \
          + fT*(T/4 + (ebm/(2*x))*((-0.5)*bp)) \
          - f/4

    return rho, pr, pt

# ── Пульс: модуляция T(r,t)=T0(r)·Θ_pulse(t), τ=1μs ────────────────────────

def pulse_energy_shell(rho_dim, pr_dim, pt_dim, x, r0_val, tau_pulse_s):
    """
    Энергия-импульс для сферического слоя толщиной Δr = r0·dx (профиль),
    при пиковой плотности ρ_peak(=ρ_stat).
    w = −ρ (для «пустоты») или +ρ (для «вещества»), берём |ρ| как
    плотность энергии в единицах 4πr0².
    T_μν ~ diag(ρ, −p_r, −p_t, −p_t).
    Энергия на единицу площади поверхности сферы: e = ∫|ρ|·dr  (безразм.)
    Полная энергия слоя: E = 4πr²·e·(c⁴/(8πG r0²)).
    """
    dr_dim = r0_val  # толщина слоя ~ r0 (один e-фолд профиля)
    rho_cgs = mp.fabs(rho_dim) * c4_over_8piGr0sq
    E = 4 * mp.pi * (x * r0_val)**2 * rho_cgs * dr_dim   # Дж
    return float(E)

# ── Сеть точек ────────────────────────────────────────────────────────────────

x_points = [1.0, 1.05, 1.2, 1.5, 10.0, 100.0]
tau_pulse = mp.mpf('1e-6')   # 1 мкс

def main():
    print("="*80)
    print("  R3 PROFILE+PULSE: α=%.2f, γ=%.3f, r0=%.1f м" % (float(alpha), float(gamma), float(r0)))
    print("="*80)

    # ── Горловина (аналитика) ──
    rho_th, pr_th, pt_th = throat_analytic(alpha, gamma)
    T_th = T_of(mp.mpf('1.0'))
    print("\n(1) Горловина x=1 (аналитика):")
    print("    T(1)    = %.6f" % float(T_th))
    print("    ρ_th    = %.6f  (4πr0²)" % float(rho_th))
    print("    p_r,th  = %.6f  (4πr0²)" % float(pr_th))
    print("    p_t,th  = %.6f  (4πr0²)" % float(pt_th))
    print("    ρ_СИ    = %.6e Па" % float(rho_th * c4_over_8piGr0sq))
    print("    p_r,СИ  = %.6e Па" % float(pr_th * c4_over_8piGr0sq))
    print("    p_t,СИ  = %.6e Па" % float(pt_th * c4_over_8piGr0sq))

    # ── (2) Профили по x ──
    print("\n(2) Профили T(x), ρ(x), p_r(x), p_t(x):")
    hdr = "%6s | %14s %14s %14s %14s | %14s %14s %14s" % (
        "x", "T(x)", "ρ(4πr0²)", "p_r(4πr0²)", "p_t(4πr0²)",
        "ρ_СИ(Па)", "p_r,СИ(Па)", "p_t,СИ(Па)")
    print("  " + hdr)
    print("  " + "-"*len(hdr))

    rows = []
    for x in x_points:
        xv = mp.mpf(str(x))
        T_val = T_of(xv)
        rho, pr, pt = field_matter(xv, alpha, gamma)
        rho_si = float(rho * c4_over_8piGr0sq)
        pr_si  = float(pr  * c4_over_8piGr0sq)
        pt_si  = float(pt  * c4_over_8piGr0sq)
        rows.append({
            'x': x,
            'r_m': float(xv * r0),
            'T': float(T_val),
            'rho_dim': float(rho),
            'pr_dim': float(pr),
            'pt_dim': float(pt),
            'rho_SI': rho_si,
            'pr_SI': pr_si,
            'pt_SI': pt_si,
        })
        print("  %6.2f | %14.6e %14.6e %14.6e %14.6e | %14.6e %14.6e %14.6e" % (
            x, float(T_val), float(rho), float(pr), float(pt),
            rho_si, pr_si, pt_si))

    # ── (3) Импульсный режим ──
    print("\n(3) Импульсный режим τ_pulse = 1 мкс:")
    print("    T(r,t) = T₀(r)·Θ_pulse(t);  пик = статический профиль.")
    print()
    hdr2 = "%6s | %14s %14s %14s | %14s %14s %14s | %12s" % (
        "x", "ρ_peak,СИ", "p_r,peak,СИ", "p_t,peak,СИ",
        "F_ρ(N/m²)", "F_pr(N/m²)", "F_pt(N/m²)", "E_pulse(Дж)")
    print("  " + hdr2)
    print("  " + "-"*len(hdr2))

    pulse_rows = []
    for row in rows:
        x = row['x']
        xv = mp.mpf(str(x))
        rho_dim = mp.mpf(str(row['rho_dim']))
        pr_dim  = mp.mpf(str(row['pr_dim']))
        pt_dim  = mp.mpf(str(row['pt_dim']))

        # Пиковые = статические (модуляция Θ_pulse не усиливает)
        rho_peak = float(rho_dim * c4_over_8piGr0sq)
        pr_peak  = float(pr_dim  * c4_over_8piGr0sq)
        pt_peak  = float(pt_dim  * c4_over_8piGr0sq)

        # Импульс: F = ∫ P dt = P·τ  (при прямоугольном импульсе)
        F_rho = rho_peak * float(tau_pulse)
        F_pr  = pr_peak  * float(tau_pulse)
        F_pt  = pt_peak  * float(tau_pulse)

        # Энергия слоя (при |ρ|>0)
        E_layer = pulse_energy_shell(rho_dim, pr_dim, pt_dim, xv, r0, tau_pulse)

        pulse_rows.append({
            'x': x,
            'rho_peak_SI': rho_peak,
            'pr_peak_SI': pr_peak,
            'pt_peak_SI': pt_peak,
            'F_rho_Nm2': F_rho,
            'F_pr_Nm2': F_pr,
            'F_pt_Nm2': F_pt,
            'E_pulse_J': E_layer,
        })
        print("  %6.2f | %14.6e %14.6e %14.6e | %14.6e %14.6e %14.6e | %12.4e" % (
            x, rho_peak, pr_peak, pt_peak, F_rho, F_pr, F_pt, E_layer))

    # ── (4) Асимптотика x→∞: γT⁻¹ растёт ──
    print("\n(4) Асимптотика x→∞ (проблема γT⁻¹):")
    print("    T ~ 2/x⁴ → 0,  γ/T ~ γx⁴/2 → −∞ (γ<0).")
    print("    ρ_dim ~ γx⁴/2 → −∞, WEC нарушен far-zone.")
    asy_rows = []
    for x_val in [10.0, 100.0, 1000.0]:
        xv = mp.mpf(str(x_val))
        T_v = T_of(xv)
        rho_v, pr_v, pt_v = field_matter(xv, alpha, gamma)
        print("    x=%6.0f  T=%.4e  ρ=%.4e  p_r=%.4e  p_t=%.4e  ρ_СИ=%.4e Па" %
              (x_val, float(T_v), float(rho_v), float(pr_v), float(pt_v),
               float(rho_v * c4_over_8piGr0sq)))
        asy_rows.append({'x': x_val, 'T': float(T_v),
                         'rho_dim': float(rho_v), 'pr_dim': float(pr_v), 'pt_dim': float(pt_v),
                         'rho_SI': float(rho_v * c4_over_8piGr0sq)})

    # ── (5) Сводная таблица (markdown) ──
    print("\n" + "="*80)
    print("  СВОДНАЯ ТАБЛИЦА (markdown)")
    print("="*80)
    print()
    print("| x | T(x) | ρ_dim | p_r,dim | p_t,dim | ρ_СИ, Па | p_r,СИ | p_t,СИ |")
    print("|---|------|-------|---------|---------|----------|--------|--------|")
    for row in rows:
        print("| %.2f | %.4e | %+.6f | %+.6f | %+.6f | %.4e | %.4e | %.4e |" % (
            row['x'], row['T'],
            row['rho_dim'], row['pr_dim'], row['pt_dim'],
            row['rho_SI'], row['pr_SI'], row['pt_SI']))
    print()

    print("### Импульсный режим (τ=1 мкс)")
    print()
    print("| x | ρ_peak,СИ | p_r,peak | p_t,peak | E_pulse, Дж |")
    print("|---|-----------|----------|----------|-------------|")
    for pr_ in pulse_rows:
        print("| %.2f | %.4e | %.4e | %.4e | %.4e |" % (
            pr_['x'], pr_['rho_peak_SI'], pr_['pr_peak_SI'],
            pr_['pt_peak_SI'], pr_['E_pulse_J']))
    print()

    # ── Ключевые числа ──
    T_peak = float(T_th)                         # T(1) в безразмерных
    T_peak_SI = T_peak / (r0**2)                 # T в м⁻²  (T(x)=T_bezr/r0²)
    rho_peak_SI = float(rho_th * c4_over_8piGr0sq)
    p_star_SI   = float(pr_th * c4_over_8piGr0sq)  # p_r на горловине

    print("### Ключевые числа:")
    print("  T_peak (x=1)  = %.6f  (безразм.)  = %.6f  м⁻²" % (T_peak, T_peak_SI))
    print("  ρ_peak (x=1)  = %.6f  (4πr0²)    = %.6e  Па" % (float(rho_th), rho_peak_SI))
    print("  p_r,peak      = %.6f  (4πr0²)    = %.6e  Па" % (float(pr_th), p_star_SI))
    print("  p_t,peak      = %.6f  (4πr0²)    = %.6e  Па" % (
        float(pt_th), float(pt_th * c4_over_8piGr0sq)))
    print("  c⁴/(8πG·r0²) = %.6e  Па" % float(c4_over_8piGr0sq))
    print()

    # ── JSON ──
    output = {
        "model": "f(T) = T + αT² + γT⁻¹",
        "params": {"alpha": float(alpha), "gamma": float(gamma), "r0_m": float(r0)},
        "units": {
            "rho_dim": "4πr0² (безразм.)",
            "rho_SI": "Па  (ρ_dim × c⁴/(8πG·r0²))",
            "c4_over_8piGr0sq_Pa": float(c4_over_8piGr0sq),
        },
        "throat_analytic": {
            "T_throat": float(T_th),
            "rho_dim": float(rho_th),
            "pr_dim": float(pr_th),
            "pt_dim": float(pt_th),
            "rho_SI": float(rho_th * c4_over_8piGr0sq),
            "pr_SI": float(pr_th * c4_over_8piGr0sq),
            "pt_SI": float(pt_th * c4_over_8piGr0sq),
        },
        "profile": rows,
        "pulse": {
            "tau_s": float(tau_pulse),
            "mode": "T(r,t)=T0(r)·Θ_pulse(t), peak = static",
            "data": pulse_rows,
        },
        "asymptotic": asy_rows,
        "key_numbers": {
            "T_peak_bezr": T_peak,
            "T_peak_m_inv2": T_peak_SI,
            "rho_peak_SI_Pa": rho_peak_SI,
            "p_star_SI_Pa": p_star_SI,
        },
        "gammaT_inv_problem": {
            "note": "T→0 при x→∞, γ/T→−∞ → ρ→−∞, WEC далёк от горловины.",
            "onset_x": float(1/abs(gamma)**0.25),   # x ~ |γ|^{-1/4}
            "onset_note": "При x > %.1f отклонение от Minkowski O(10%%)" % float(1/abs(gamma)**0.25),
        },
    }

    outpath = "/home/smboozha/portal_gun/research/r3_profile_pulse_a1_results.json"
    with open(outpath, "w") as fp:
        json.dump(output, fp, indent=2, default=str)
    print("  → %s сохранён" % outpath)

if __name__ == "__main__":
    main()
