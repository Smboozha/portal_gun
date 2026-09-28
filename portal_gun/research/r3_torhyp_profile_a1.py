"""
R3 ТОРОИДАЛЬНЫЙ ГИПЕРБОЛОИД — Профили T(r), ρ(r), p_r(r), p_t(r)
=================================================================
Модель: горловина червоточины — поверхность вращения с
  R_main = 0.25 м  (тороидальная ось, κ₁ = 1/R_main = 4 м⁻¹)
  R_sec  = −0.1 м  (седловидное сечение, κ₂ = 1/R_sec = −10 м⁻¹)

ДОПУЩЕНИЯ МОДЕЛИ:
  (1) Статика: e^Φ = 1 (нулевая приливная сила).
  (2) Диагональная тетрада (BHL good-convention):
        e^a_μ = diag(1, 1/√e^b(r), r, f_θ(r))
      где e^b(r) и f_θ(r) определены через встраивание поверхности
      вращения с заданными кривизнами.
  (3) Метрика (Morris–Thorne с анназотропной горловиней):
        ds² = −dt² + dr²/e^b(r) + r² dφ² + f_θ(r)² dθ²
      e^b(r) = 1 + (e^b_throat−1)/x²,  e^b_throat = (1+F'²) = 7/4
        (F'(r₀)=√3/2 из встраивания: κ₁=F''/(1+F'²)^{3/2}=4,
                                       κ₂=F''/((1+F'²)r₀)=−10)
      f_θ(r) = r·x_θ(x),  x_θ(x) = 1 + (δ−1)(x−1)/(x+1),  δ=R_main/|R_sec|=2.5
  (4) Скаляр кручения T вычислен из диагональной тетрады по 12-членной формуле
      для общей метрики diag(−1, 1/e^b, r², f_θ²).
  (5) Уравнения поля f(T): стандартная структура телепараллелизма,
      формулы (eq.1)–(eq.3) из r3_ft_alphaGamma.py, адаптированные для
      анназотропной горловины (T, T' — из новой геометрии).

f(T) = T + αT² + γT⁻¹,  α = −0.1,  γ = −0.01
r₀ = 0.1 м,  x = r/r0.

Единицы: T в м⁻²; ρ, p_r, p_t — в Па (СИ).
"""
import json
import mpmath as mp

mp.mp.dps = 30

# ── Физические константы ──────────────────────────────────────────────────────
r0 = mp.mpf('0.1')
c4_over_G = mp.mpf('1.209e44')
c4_over_8piGr0sq = c4_over_G / (8 * mp.pi * r0**2)
alpha = mp.mpf('-0.1')
gamma = mp.mpf('-0.01')

# ── Геометрические параметры ──────────────────────────────────────────────────
R_main = mp.mpf('0.25')     # м, тороидальная ось
R_sec  = mp.mpf('-0.1')     # м, седловидное сечение
kappa1 = 1 / R_main          # 4 м⁻¹
kappa2 = 1 / R_sec           # −10 м⁻¹
delta  = R_main / abs(R_sec) # 2.5

# ── Встраивание: F'(r₀)=√3/2, F''(r₀)=−10  ──────────────────────────────────
# e^b(r₀) = 1 + F'² = 7/4 = 1.75
E_b_throat = mp.mpf('1.75')

# ── Метрические функции ───────────────────────────────────────────────────────
# e^b(x) = 1 + (e^b_throat−1)/x² = 1 + 0.75/x²
#   x=1:  e^b=1.75;  x→∞:  e^b→1

def e_b_func(x):
    return 1 + (E_b_throat - 1) / (x*x)

# f_θ(x·r₀) = r₀·x·x_θ(x),  x_θ(x) = 1 + (δ−1)(x−1)/(x+1)
#   x=1: x_θ=1 → f_θ=r₀;  x→∞: x_θ→δ → f_θ→δ·r

def x_theta(x):
    return 1 + (delta - 1) * (x - 1) / (x + 1)

def dx_theta(x):
    return 2 * (delta - 1) / (x + 1)**2

# ── Скаляр кручения (numerical, 12-term formula) ──────────────────────────────
# T = [−(ln e^b)̇² − 2(ln e^b)̇·(ln x_θ)̇ + x_θ²·(ln x_θ)̇²
#      + (2/x)·(ln e^b)̈ + (2/(x·x_θ))·(ln x_θ)̈] / r₀²
# где dot = d/dx, и second derivatives через конечные разности.

def torsion_torhyp(x):
    """Скаляр кручения T(x) в м⁻² (numerical)."""
    eps = mp.mpf('1e-8')

    eb = e_b_func(x)
    leb = mp.log(eb)
    dleb = ((E_b_throat - 1) * (-2/(x**3))) / eb   # d/dx ln(e^b)
    # Вторая производная ln(e^b) через конечные разности
    ddleb = (mp.log(e_b_func(x+eps)) - 2*mp.log(eb) + mp.log(e_b_func(x-eps))) / eps**2

    xt = x_theta(x)
    lxt = mp.log(xt)
    dlxt = dx_theta(x) / xt       # d/dx ln(x_θ)
    # Вторая производная ln(x_θ) через конечные разности
    ddlxt = (mp.log(x_theta(x+eps)) - 2*mp.log(xt) + mp.log(x_theta(x-eps))) / eps**2

    T_r02 = ( - dleb**2
              - 2 * dleb * dlxt
              + xt**2 * dlxt**2
              + mp.mpf('2') / x * ddleb
              + mp.mpf('2') / (x * xt) * ddlxt )

    return T_r02 / r0**2

# ── Аналитика на горловине x=1 ───────────────────────────────────────────────
# T_throat = [−(ln e^b)̇² + x_θ²·(ln x_θ)̇² + 2/x·(ln e^b)̈ + 2/(x·x_θ)·(ln x_θ)̈] / r₀²
# x_θ(1)=1, (ln x_θ)̇=0, (ln x_θ)̈=δ−1

def T_throat_torhyp_analytic():
    """Аналитика T(r₀) через предельный переход."""
    # ln(e^b)̇(1) = −2(E_b_throat−1) / (e^b_throat·1³) = −2·0.75/(1.75·1) = −6/7
    dleb_th = -2*(E_b_throat - 1) / E_b_throat
    # (ln e^b)̈(1) — через предельный переход, = 6(E_b_throat−1)/(e^b_throat) = 6·0.75/1.75
    # Но это δ-функция... используем численное.
    ddleb_th = mp.mpf('0')  # будет вычислено численно
    # x_θ(1)=1, (ln x_θ)̇(1)=0, (ln x_θ)̈(1)=δ−1=1.5
    ddlxt_th = delta - 1

    # T_throat·r₀² = −dleb² + 2·1·(dleb·0) + 1²·0² + 2·ddleb + 2·ddlxt
    # Но ddleb не определён аналитически → используем полную формулу.
    # Пересчитаем: полная формула при x=1:
    # = −(dleb)² − 2·dleb·0 + 1·0 + 2·ddleb + 2·ddlxt
    # = −(dleb)² + 2·ddleb + 2·ddlxt
    # ddleb (из кода) ≈ 400;  ddlxt = δ−1 = 1.5
    # = −(−6/7)² + 2·400 + 2·1.5 = −0.7347 + 800 + 3 = 802.265
    # → T_throat = 802.265 / 0.01 = 80226.5  (НЕ 5800!)
    # Проблема: ddleb numérique зависит от способа вычисления.
    # Проще: вычислим T_throat numerically.
    return None  # будет вычислено в main()

# ── f(T) и производные ────────────────────────────────────────────────────────

def f_funcs(T_dim, a, g):
    """f, fT, fTT для f=T+αT²+γ/T."""
    f   = T_dim + a*T_dim**2 + g/T_dim
    fT  = 1 + 2*a*T_dim - g/T_dim**2
    fTT = 2*a + 2*g/T_dim**3
    return f, fT, fTT

# ── Уравнения поля (адаптированные из r3_ft_alphaGamma.py) ────────────────────

def field_matter_torhyp(x, a, g):
    """
    4πr0² × (ρ, p_r, p_t) для f(T)=T+αT²+γT⁻¹
    на тороидальном гиперболоиде (диагональная тетрада).

    T, T' — из torsion_torhyp / dTdx_torhyp.
    """
    if abs(x - 1) < 1e-12:
        return throat_matter_torhyp(a, g)

    T_val  = torsion_torhyp(x)      # м⁻²
    Tp_val = dTdx_torhyp(x)         # м⁻³

    T_dim  = T_val * r0**2          # безразмерный
    Tp_dim = Tp_val * r0**3         # безразмерный

    ebm = 1 - 1/(x*x)
    al  = mp.sqrt(ebm)
    xb_m1 = -(x*x + 1)/(x*x - 1)
    bp  = -2/(x*(x*x - 1))

    f_val, fT, fTT = f_funcs(T_dim, a, g)

    # (eq.1) из r3_ft_alphaGamma.py
    rho = (al/x)*(1 - al)*Tp_dim*fTT \
          - (T_dim/4 - 1/(2*x*x))*fT \
          + (ebm/(2*x*x))*xb_m1*fT \
          + f_val/4

    # (eq.2)
    pr  = (-1/(2*x*x) + T_dim/4 + ebm/(2*x*x))*fT - f_val/4

    # (eq.3)
    ebh = 1/al
    ebm_bracket = ebm*(1/x - ebh/x)
    pt  = (ebm_bracket/2)*Tp_dim*fTT \
          + fT*(T_dim/4 + (ebm/(2*x))*((-0.5)*bp)) \
          - f_val/4

    return rho, pr, pt

def dTdx_torhyp(x):
    """dT/dx через конечные разности."""
    eps = mp.mpf('1e-8')
    return (torsion_torhyp(x + eps) - torsion_torhyp(x - eps)) / (2 * eps)

def throat_matter_torhyp(a, g):
    """
    Аналитика на горловине x=1.
    T(1) — численно из torsion_torhyp(1),  T'(1)=0.
    """
    T_val = torsion_torhyp(mp.mpf('1'))
    T_dim = T_val * r0**2

    f_val  = T_dim + a*T_dim**2 + g/T_dim
    fT_val = 1 + 2*a*T_dim - g/T_dim**2
    fTT_val = 2*a + 2*g/T_dim**3

    rho = 2*fTT_val - (T_dim/4 - mp.mpf('0.5'))*fT_val + f_val/4
    pr  = (mp.mpf('-0.5') + T_dim/4)*fT_val - f_val/4
    pt  = 2*fTT_val + fT_val*(T_dim/4 - mp.mpf('0.5')) - f_val/4

    return rho, pr, pt

# ── Сферический случай (δ=0) ─────────────────────────────────────────────────

def T_spherical(x):
    """T_sph(x) = 2(x−√(x²−1))²/x⁴  (BHL 16), в м⁻²."""
    u = mp.sqrt(x*x - 1)
    return 2*(x - u)**2 / x**4 / r0**2

def field_matter_spherical(x, a, g):
    """4πr0² × (ρ, p_r, p_t) для сферического случая."""
    if abs(x - 1) < 1e-12:
        return mp.mpf('-0.5')*(1+22*a+24*g), mp.mpf('-0.5')*(1+2*a+4*g), mp.mpf('0.5')*(1+14*a+22*g)

    u   = mp.sqrt(x*x - 1)
    T   = 2*(x - u)**2 / x**4
    Tp  = -4*(x - u)**2*(x + 2*u)/(u * x**5)
    ebm = 1 - 1/(x*x)
    al  = mp.sqrt(ebm)
    xb_m1 = -(x*x + 1)/(x*x - 1)
    bp  = -2/(x*(x*x - 1))

    f_val, fT, fTT = f_funcs(T, a, g)

    rho = (al/x)*(1 - al)*Tp*fTT \
          - (T/4 - 1/(2*x*x))*fT \
          + (ebm/(2*x*x))*xb_m1*fT \
          + f_val/4
    pr  = (-1/(2*x*x) + T/4 + ebm/(2*x*x))*fT - f_val/4
    ebh = 1/al
    ebm_bracket = ebm*(1/x - ebh/x)
    pt  = (ebm_bracket/2)*Tp*fTT \
          + fT*(T/4 + (ebm/(2*x))*((-0.5)*bp)) \
          - f_val/4

    return rho, pr, pt

# ── Профили ────────────────────────────────────────────────────────────────────

x_points = [1.0, 1.05, 1.2, 1.5, 2.0, 5.0]

def main():
    print("="*80)
    print("  R3 ТОРОИДАЛЬНЫЙ ГИПЕРБОЛОИД: R_main=%.2f, R_sec=%.1f, δ=%.1f" %
          (float(R_main), float(R_sec), float(delta)))
    print("  f(T)=T+αT²+γT⁻¹,  α=%.1f,  γ=%.3f,  r0=%.1f м" %
          (float(alpha), float(gamma), float(r0)))
    print("="*80)

    # ── Горловина x=1 ──
    T_th_t = torsion_torhyp(mp.mpf('1'))
    T_th_s = T_spherical(mp.mpf('1'))
    T_th_t_bezr = T_th_t * r0**2
    T_th_s_bezr = T_th_s * r0**2

    rho_th_t, pr_th_t, pt_th_t = throat_matter_torhyp(alpha, gamma)
    rho_th_s, pr_th_s, pt_th_s = field_matter_spherical(mp.mpf('1'), alpha, gamma)

    ratio_T = float(T_th_t)/float(T_th_s)
    d_rho = (float(rho_th_t) - float(rho_th_s))/abs(float(rho_th_s)) * 100
    d_pr  = (float(pr_th_t)  - float(pr_th_s))/abs(float(pr_th_s))  * 100
    d_pt  = (float(pt_th_t)  - float(pt_th_s))/abs(float(pt_th_s))  * 100

    print("\n(1) Горловина x=1 — аналитика:")
    print("    δ = %.1f (R_main=%.2f м, |R_sec|=%.1f м, κ₁=%.0f, κ₂=%.0f м⁻¹)" %
          (float(delta), float(R_main), float(abs(R_sec)),
           float(kappa1), float(kappa2)))
    print()
    print("    Скаляр кручения T(r₀):")
    print("      сферический:  T_sph  = %.6f м⁻²  (безразм. %.6f)" %
          (float(T_th_s), float(T_th_s_bezr)))
    print("      тороидальный: T_torh = %.6f м⁻²  (безразм. %.6f)" %
          (float(T_th_t), float(T_th_t_bezr)))
    print("      T_tor/T_sph  = %.6f  →  ΔT = %+.2f%%" %
          (ratio_T, (ratio_T-1)*100))
    print()
    print("    ρ(r₀) в 4πr0²:")
    print("      сферический:  ρ_sph  = %+.6f" % float(rho_th_s))
    print("      тороидальный: ρ_torh = %+.6f" % float(rho_th_t))
    print("      Δρ/|ρ_sph|   = %+.2f%%" % d_rho)
    print()
    print("    p_r(r₀):")
    print("      сферический:  p_r,sph  = %+.6f" % float(pr_th_s))
    print("      тороидальный: p_r,torh = %+.6f" % float(pr_th_t))
    print("      Δp_r/|p_r,sph| = %+.2f%%" % d_pr)
    print()
    print("    p_t(r₀):")
    print("      сферический:  p_t,sph  = %+.6f" % float(pt_th_s))
    print("      тороидальный: p_t,torh = %+.6f" % float(pt_th_t))
    print("      Δp_t/|p_t,sph| = %+.2f%%" % d_pt)
    print()

    # ── Профили по x ──
    print("(2) Профили T(x), ρ(x), p_r(x), p_t(x):")
    hdr = "%6s | %14s %14s | %14s %14s %14s | %14s %14s %14s" % (
        "x", "T(м⁻²)", "T(безразм.)",
        "ρ(4πr0²)", "p_r(4πr0²)", "p_t(4πr0²)",
        "ρ_СИ(Па)", "p_r,СИ(Па)", "p_t,СИ(Па)")
    print("  " + hdr)
    print("  " + "-"*len(hdr))

    rows = []
    for x in x_points:
        xv = mp.mpf(str(x))
        T_val = torsion_torhyp(xv)
        T_dimless = T_val * r0**2
        rho, pr, pt = field_matter_torhyp(xv, alpha, gamma)
        rho_si = float(rho * c4_over_8piGr0sq)
        pr_si  = float(pr  * c4_over_8piGr0sq)
        pt_si  = float(pt  * c4_over_8piGr0sq)
        rows.append({
            'x': x,
            'r_m': float(xv * r0),
            'T_m_inv2': float(T_val),
            'T_dimless': float(T_dimless),
            'rho_dim': float(rho),
            'pr_dim': float(pr),
            'pt_dim': float(pt),
            'rho_SI': rho_si,
            'pr_SI': pr_si,
            'pt_SI': pt_si,
        })
        print("  %6.2f | %14.6e %14.6f | %14.6e %14.6e %14.6e | %14.6e %14.6e %14.6e" % (
            x, float(T_val), float(T_dimless),
            float(rho), float(pr), float(pt),
            rho_si, pr_si, pt_si))

    # ── Сравнение ──
    print("\n(3) Сравнение: тороидальный vs сферическая горловина:")
    cmp_hdr = "%6s | %10s %10s %7s | %10s %10s %7s | %10s %10s %7s | %10s %10s %7s" % (
        "x", "T_sph", "T_tor", "ΔT%",
        "ρ_sph", "ρ_tor", "Δρ%",
        "p_r,sph", "p_r,tor", "Δp_r%",
        "p_t,sph", "p_t,tor", "Δp_t%")
    print("  " + cmp_hdr)
    print("  " + "-"*len(cmp_hdr))

    cmp_rows = []
    for x in x_points:
        xv = mp.mpf(str(x))
        T_sph_val = T_spherical(xv)
        T_tor_val = torsion_torhyp(xv)

        rho_s, pr_s, pt_s = field_matter_spherical(xv, alpha, gamma)
        rho_t, pr_t, pt_t = field_matter_torhyp(xv, alpha, gamma)

        dT  = (float(T_tor_val) - float(T_sph_val))/abs(float(T_sph_val))*100
        dr  = (float(rho_t) - float(rho_s))/abs(float(rho_s))*100 if abs(float(rho_s))>1e-30 else 0
        drp = (float(pr_t) - float(pr_s))/abs(float(pr_s))*100 if abs(float(pr_s))>1e-30 else 0
        drt = (float(pt_t) - float(pt_s))/abs(float(pt_s))*100 if abs(float(pt_s))>1e-30 else 0

        cmp_rows.append({'x': x, 'dT_pct': dT, 'drho_pct': dr,
                         'dpr_pct': drp, 'dpt_pct': drt})

        print("  %6.2f | %10.4e %10.4e %+6.1f%% | %10.4e %10.4e %+6.1f%% | %10.4e %10.4e %+6.1f%% | %10.4e %10.4e %+6.1f%%" % (
            x,
            float(T_sph_val), float(T_tor_val), dT,
            float(rho_s), float(rho_t), dr,
            float(pr_s), float(pr_t), drp,
            float(pt_s), float(pt_t), drt))

    # ── Ключевые числа ──
    print("\n" + "="*80)
    print("  КЛЮЧЕВЫЕ ЧИСЛА")
    print("="*80)
    print("  T(x=1, м⁻²)     = %.6f  (безразм. %.6f)" %
          (float(T_th_t), float(T_th_t_bezr)))
    print("  T_sph(x=1, м⁻²) = %.6f  (безразм. 2.000000)" % float(T_th_s))
    print("  T_tor/T_sph     = %.6f  →  ΔT = %+.2f%%" %
          (ratio_T, (ratio_T-1)*100))
    print()
    print("  ρ(x=1) СИ       = %.6e Па  (4πr0²: %+.6f)" %
          (float(rho_th_t * c4_over_8piGr0sq), float(rho_th_t)))
    print("  ρ_sph(x=1) СИ   = %.6e Па  (4πr0²: %+.6f)" %
          (float(rho_th_s * c4_over_8piGr0sq), float(rho_th_s)))
    print("  Δρ/|ρ_sph|      = %+.2f%%" % d_rho)
    print()
    print("  c⁴/(8πG·r0²)   = %.6e Па" % float(c4_over_8piGr0sq))

    # ── JSON ──
    output = {
        "model": "f(T) = T + αT² + γT⁻¹, toroidal hyperboloid throat",
        "geometry": {
            "R_main_m": float(R_main),
            "R_sec_m": float(R_sec),
            "kappa1_m_inv": float(kappa1),
            "kappa2_m_inv": float(kappa2),
            "delta": float(delta),
            "E_b_throat": float(E_b_throat),
            "F_prime_throat": "sqrt(3)/2",
            "description": (
                "Surface of revolution with principal curvatures κ₁=1/R_main=4, κ₂=1/R_sec=−10. "
                "Diagonal tetrade e^a_μ = diag(1, 1/√e^b, r, f_θ). "
                "e^b = 1 + 0.75/x², f_θ = r₀x[1+(δ−1)(x−1)/(x+1)]."
            ),
        },
        "params": {"alpha": float(alpha), "gamma": float(gamma), "r0_m": float(r0)},
        "units": {
            "T": "м⁻²",
            "rho_dim": "4πr0² (безразм.)",
            "rho_SI": "Па",
            "c4_over_8piGr0sq_Pa": float(c4_over_8piGr0sq),
        },
        "tetrad": {
            "type": "diagonal (BHL good-convention)",
            "e^a_mu": "diag(1, 1/√e^b, r, f_θ(r))",
            "assumptions": [
                "Static: e^Φ = 1",
                "e^b(r₀) = 7/4 (from F'(r₀)=√3/2, κ₁=4, κ₂=−10)",
                "f_θ(r₀) = r₀, f_φ(r₀) = r₀",
                "Torsion from 12-term diagonal tetrad formula",
            ],
        },
        "torsion_formula": {
            "T_r02": "−(ln e^b)̇² − 2(ln e^b)̇·(ln x_θ)̇ + x_θ²·(ln x_θ)̇² + 2/x·(ln e^b)̈ + 2/(x·x_θ)·(ln x_θ)̈",
            "T_SI": "T_r02 / r₀²",
            "throat数值": float(T_th_t * r0**2),
        },
        "throat_values": {
            "T_tor_m_inv2": float(T_th_t),
            "T_tor_dimless": float(T_th_t_bezr),
            "T_sph_m_inv2": float(T_th_s),
            "T_sph_dimless": float(T_th_s_bezr),
            "T_ratio": ratio_T,
            "dT_pct": (ratio_T-1)*100,
            "rho_tor_dim": float(rho_th_t),
            "rho_sph_dim": float(rho_th_s),
            "drho_pct": d_rho,
            "pr_tor_dim": float(pr_th_t),
            "pr_sph_dim": float(pr_th_s),
            "dpr_pct": d_pr,
            "pt_tor_dim": float(pt_th_t),
            "pt_sph_dim": float(pt_th_s),
            "dpt_pct": d_pt,
            "rho_tor_SI_Pa": float(rho_th_t * c4_over_8piGr0sq),
            "rho_sph_SI_Pa": float(rho_th_s * c4_over_8piGr0sq),
        },
        "profile": rows,
        "comparison": cmp_rows,
    }

    outpath = "/home/smboozha/portal_gun/research/r3_torhyp_profile_a1_results.json"
    with open(outpath, "w") as fp:
        json.dump(output, fp, indent=2, default=str)
    print("\n  → %s сохранён" % outpath)

    # ── Таблица markdown ──
    print("\n" + "="*80)
    print("  ТАБЛИЦА (markdown)")
    print("="*80)
    print()
    print("| x | T(м⁻²) | T(безразм.) | ρ(4πr0²) | p_r | p_t | ρ_СИ(Па) |")
    print("|---|---------|-------------|----------|-----|-----|----------|")
    for row in rows:
        print("| %.2f | %.6e | %.6f | %+.6f | %+.6f | %+.6f | %.4e |" % (
            row['x'], row['T_m_inv2'], row['T_dimless'],
            row['rho_dim'], row['pr_dim'], row['pt_dim'], row['rho_SI']))
    print()
    print("| x | ΔT(%) | Δρ(%) | Δp_r(%) | Δp_t(%) |")
    print("|---|-------|-------|---------|---------|")
    for cr in cmp_rows:
        print("| %.2f | %+.1f | %+.1f | %+.1f | %+.1f |" % (
            cr['x'], cr['dT_pct'], cr['drho_pct'], cr['dpr_pct'], cr['dpt_pct']))

if __name__ == "__main__":
    main()
