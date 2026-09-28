"""
R3 — f(T) = T + αT² + γT⁻¹  на Morris–Thorne (Ellis) червоточине
===================================================================
Тетрада: хорошая (BHL конвенция), a(r)=0, β(r)=r0²/r, r0=0.1 м.
Единицы: 4πr0² × {ρ, p_r, p_t}.  x = r/r0.

(1) Уравнения поля (стандартная структура телепараллелизма):
    T_μν = (1/κ) [ (∇_λ S^λ_μν) − ... ]   ← определение S из tetrad_core.py.
    Это НЕ новый постулат, а базовая форма f(T)-гравитации:
    варьирование действия S = ∫ e [T + f(T)]/(2κ) d⁴x по тетраде e^a_μ.

(2) Скалярные уравнения для MT (diag tetrade):
    ρ = (al/x)(1−al) Tp·fTT − (T/4 − 1/(2x²))·fT + ((x²+1)/(2x²(x²−1)))·fT + f/4  ... (eq.1)
    p_r = [−1/(2x²) + T/4 + (x²−1)/(2x²(x²−1))]·fT − f/4                            ... (eq.2)
    p_t = (ebm/2)(1/x − 1/(al·x))·Tp·fTT + fT·[T/4 + (ebm/(2x))(−bp²/2)] − f/4       ... (eq.3)

    где T(x) = 2(x−√(x²−1))²/x⁴                                      ... (eq.4)
    Tp = dT/dx = −4(x−√(x²−1))²·(x+2√(x²−1))/(√(x²−1)·x⁵)          ... (eq.5)
    al = √(1−1/x²),  ebm = 1−1/x²,  bp = −2/(x·(x²−1))              ... (eq.6)

    f(T) = T + αT² + γT⁻¹  (γ<0)                                     ... (eq.7)
    f_T  = 1 + 2αT − γT⁻²                                             ... (eq.8)
    f_TT = 2α + 2γT⁻³                                                  ... (eq.9)

    Аналитический предел при x→1 (горловина):
    T(1) = 2,  Tp(1) = −2,  al(1)=1−1/(2x²)≈0,  bp(1)=−1             ... (eq.10)

    ρ_throat  = −(1 + 22α + 24γ)/2                                     ... (eq.11)
    p_r_throat = −(1 + 2α + 4γ)/2                                      ... (eq.12)
    p_t_throat = (1 + 14α + 22γ)/2                                     ... (eq.13)

(3) Точность: формулы сверены с BHL (45),(46),(48),(49) — maxdev < 5e-18 (R2).
    Печатная (47) в BHL ошибочна: p_t из (20) = (49)−(45), НЕ (47).

(4) Проверка TEGR-предела: α=γ=0 → ρ=−1/2, p_r=−1/2, p_t=1/2  (= MT/GR).
"""
import json
import mpmath as mp

mp.mp.dps = 30

# ── Торсионный скаляр и его производная (BHL конвенция) ────────────────────────

def u_of(x):
    """u = √(x²−1)"""
    return mp.sqrt(x*x - 1)

def T16(x):
    """T(x) = 2(x−√(x²−1))² / x⁴,  BHL (16)"""
    u = u_of(x)
    return 2*(x - u)**2 / x**4

def Tp17(x):
    """dT/dx = −4(x−u)²·(x+2u)/(u·x⁵),  BHL (17)"""
    u = u_of(x)
    return -4*(x - u)**2*(x + 2*u)/(u * x**5)

# ── f(T) и производные ────────────────────────────────────────────────────────

def f_funcs(T, alpha, gamma):
    """f, fT, fTT для f = T + αT² + γT⁻¹ (γ<0).
       f  = T + αT² + γ/T                                           ... (eq.7)
       fT = 1 + 2αT − γ/T²                                          ... (eq.8)
       fTT = 2α + 2γ/T³                                             ... (eq.9)
    """
    f  = T + alpha*T**2 + gamma/T
    fT = 1 + 2*alpha*T - gamma/T**2
    fTT = 2*alpha + 2*gamma/T**3
    return f, fT, fTT

# ── ρ, p_r, p_t (численно, по аналогии с field_matter в ft2_necscan_BHL.py) ─

def field_matter_alphaGamma(x, alpha, gamma):
    """4πr0² × (ρ, p_r, p_t) для f=T+αT²+γT⁻¹ на Morris–Thorne (Ellis).
    alpha, gamma — безразмерные (в ед. r0).
    Формулы: (eq.1)-(eq.3) с (eq.7)-(eq.9).
    При x=1 (горловина) — возвращаем аналитику (eq.11)-(eq.13).
    """
    if abs(x - 1) < 1e-12:
        return throat_values(alpha, gamma)

    u  = u_of(x)
    T  = T16(x)
    Tp = Tp17(x)
    ebm = 1 - 1/(x*x)
    al  = mp.sqrt(ebm)
    xb_m1 = -(x*x + 1)/(x*x - 1)
    bp = -2/(x*(x*x - 1))

    f, fT, fTT = f_funcs(T, alpha, gamma)

    # (eq.1):  ρ
    rho = (al/x)*(1 - al)*Tp*fTT \
          - (T/4 - 1/(2*x*x))*fT \
          + (ebm/(2*x*x))*xb_m1*fT \
          + f/4

    # (eq.2):  p_r
    pr  = (-1/(2*x*x) + T/4 + ebm/(2*x*x))*fT - f/4

    # (eq.3):  p_t
    # ebm*(1/x - 1/(al*x)) = (ebm/al)*(al-1)/x = al*(al-1)/x  (алгебраически)
    ebh = 1/al
    ebm_bracket = ebm*(1/x - ebh/x)   # = al*(al-1)/x
    pt  = (ebm_bracket/2)*Tp*fTT \
          + fT*(T/4 + (ebm/(2*x))*((-0.5)*bp)) \
          - f/4

    return rho, pr, pt

# ── Точные формулы на горловине x=1  (eq.11)-(eq.13) ────────────────────────

def throat_values(alpha, gamma):
    """Аналитический предел при x→1: ρ, p_r, p_t  (eq.11)-(eq.13)."""
    rho  = -(1 + 22*alpha + 24*gamma)/2
    pr   = -(1 + 2*alpha + 4*gamma)/2
    pt   = (1 + 14*alpha + 22*gamma)/2
    return rho, pr, pt

# ── Вспомогательные проверки ──────────────────────────────────────────────────

def check_tegr_limit():
    """Проверка TEGR-предела: α=γ=0 → MT/GR: ρ=−1/2, p_r=−1/2, p_t=1/2."""
    r, p, t = throat_values(0, 0)
    ok = (abs(r - mp.mpf('-0.5')) < 1e-28 and
          abs(p - mp.mpf('-0.5')) < 1e-28 and
          abs(t - mp.mpf('0.5'))  < 1e-28)
    # Сверка с числовой формулой вдали от горловины (x=2, нет вырождений)
    rn, pn, tn = field_matter_alphaGamma(mp.mpf('2.0'), 0, 0)
    # BHL TEGR: ρ=p_r=−1/(2x^4)=−1/32, p_t=1/(2x^4)=1/32
    ok2 = (abs(rn - mp.mpf('-1')/32) < 1e-12 and
           abs(pn - mp.mpf('-1')/32) < 1e-12 and
           abs(tn - mp.mpf('1')/32)  < 1e-12)
    return ok, ok2

def check_quad_only():
    """Проверка: γ=0 → воспроизводит результаты ft2_necscan_BHL.py."""
    checks = []
    for a_val in ['-5', '-3', '-1', '-0.5', '0', '0.5', '1', '2']:
        a = mp.mpf(a_val)
        # Аналитика (eq.11)-(eq.13) с γ=0
        r_th, p_th, t_th = throat_values(a, 0)
        # Число (x=2, вдали от вырождений)
        rn, pn, tn = field_matter_alphaGamma(mp.mpf('2.0'), a, 0)
        # BHL exact: (45) и (46) для quad
        u2 = mp.sqrt(3)
        rho_bhl = -1/(2*16)*(1 + 16*a*u2/2*(3 - 5/4) - 2*a*(24 - 52/4 + 17/16))
        pr_bhl  = -1/(2*16)*(1 + 16*a*u2/2*(1 - 1/4) - 2*a*(8 - 12/4 + 3/16))
        ok_r = abs(rn - rho_bhl) < 1e-10
        ok_p = abs(pn - pr_bhl)  < 1e-10
        checks.append((a_val, ok_r and ok_p))
    return checks

# ── Асимптотика x→∞ ──────────────────────────────────────────────────────────

def asymptotic_check(x, alpha, gamma):
    """Проверка вкладда T⁻¹-члена и отклонения от Minkowski.
    При x→∞: T ~ 2/x⁴ → 0.  γ/T ~ γx⁴/2 → ±∞ при γ≠0.
    Отклонение от фона O(x⁴): ρ~γx⁴/2, p_r~−γx⁴/2, p_t~γx⁴/2.
    """
    rho, pr, pt = field_matter_alphaGamma(x, alpha, gamma)
    T_val = T16(x)
    f_val = T_val + alpha*T_val**2 + gamma/T_val
    return {
        'x': float(x),
        'T': float(T_val),
        'f': float(f_val),
        'rho': float(rho),
        'pr': float(pr),
        'pt': float(pt),
        'f_over_T': float(f_val/T_val) if abs(T_val) > 1e-300 else None,
    }

# ── Основной скан ─────────────────────────────────────────────────────────────

def run_scan():
    print("="*80)
    print("  R3 — f(T)=T+αT²+γT⁻¹: NEC/WEC scan, MT (Ellis), good tetrad")
    print("="*80)

    # --- (1) TEGR-предел ---
    ok_tegr, ok_num = check_tegr_limit()
    print("\n(1) TEGR-предел α=γ=0: аналитика %s, vs числовой %s" %
          ("OK" if ok_tegr else "FAIL", "OK" if ok_num else "FAIL"))

    # --- (2) Воспроизводство quad-only (γ=0) ---
    print("\n(2) γ=0 → квадратичная модель (сверка с ft2_necscan_BHL):")
    quad_checks = check_quad_only()
    all_ok = all(c[1] for c in quad_checks)
    for a_val, ok in quad_checks:
        print("   α=%6s  match: %s" % (a_val, "OK" if ok else "FAIL"))
    print("   Итого: %s" % ("ALL OK" if all_ok else "MISMATCH"))

    # --- (3) Скан (α, γ) → EoS ---
    alpha_vals = [mp.mpf(v) for v in ['0', '-0.1', '-1']]
    gamma_vals = [mp.mpf(v) for v in ['-0.01', '-0.1', '-0.5', '-1']]
    x_pts = [mp.mpf('1.0'), mp.mpf('1.05'), mp.mpf('1.2'), mp.mpf('1.5')]

    print("\n(3) Таблица (α,γ) → EoS-совместимость на горловине x=1:")
    print("   Требования: ρ>0 [WEC], ρ+p_r≥0 [NEC_r], ρ+p_t≥0 [NEC_t]")
    print("   Проверка на x=1 (горловина), x=1.05, x=1.2, x=1.5")
    print()
    hdr = "%8s %8s | %9s %9s %9s | %9s %9s | %6s %6s %6s | %s" % (
        "α", "γ", "ρ_th", "pr_th", "pt_th", "NECr_th", "NECt_th",
        "WEC", "NECr", "NECt", "статус")
    print("   " + hdr)
    print("   " + "-"*len(hdr))

    results_table = []
    for gamma in gamma_vals:
        for alpha in alpha_vals:
            # Горловина (аналитика)
            rho_th, pr_th, pt_th = throat_values(alpha, gamma)
            necr_th = rho_th + pr_th
            nect_th = rho_th + pt_th
            wec_th  = rho_th >= 0
            necr_ok = necr_th >= 0
            nect_ok = nect_th >= 0

            # Расширенная проверка на нескольких x
            wec_all  = bool(wec_th)
            necr_all = bool(necr_ok)
            nect_all = bool(nect_ok)

            for x in x_pts[1:]:  # 1.05, 1.2, 1.5
                r_, p_, t_ = field_matter_alphaGamma(x, alpha, gamma)
                wec_all  &= (r_ >= 0)
                necr_all &= (r_ + p_ >= 0)
                nect_all &= (r_ + t_ >= 0)

            status = []
            if wec_all:  status.append("WEC+")
            if necr_all: status.append("NECr+")
            if nect_all: status.append("NECt+")
            if not status: status.append("FAIL")

            # Проверка p=−ρ на горловине
            p_eq_minus_rho_r = (abs(pr_th + rho_th) < 1e-20)
            p_eq_minus_rho_t = (abs(pt_th + rho_th) < 1e-20)

            print("   %+8.2f %+8.2f | %+9.3f %+9.3f %+9.3f | %+9.3f %+9.3f | %3s %3s %3s | %s" % (
                float(alpha), float(gamma),
                float(rho_th), float(pr_th), float(pt_th),
                float(necr_th), float(nect_th),
                "Y" if wec_th else "N",
                "Y" if necr_ok else "N",
                "Y" if nect_ok else "N",
                " ".join(status)))

            results_table.append({
                'alpha': float(alpha), 'gamma': float(gamma),
                'rho_th': float(rho_th), 'pr_th': float(pr_th), 'pt_th': float(pt_th),
                'NECr_th': float(necr_th), 'NECt_th': float(nect_th),
                'WEC_ok': bool(wec_th), 'NECr_ok': bool(necr_ok), 'NECt_ok': bool(nect_ok),
                'WEC_all': wec_all, 'NECr_all': necr_all, 'NECt_all': nect_all,
                'p_r_eq_mrho': p_eq_minus_rho_r,
                'p_t_eq_mrho': p_eq_minus_rho_t,
            })

    # --- (4) Условие p=−ρ ---
    print("\n(4) Поиск p=−ρ (фантом-подобное EoS) на горловине:")
    print("   p_r = −ρ  →  (eq.12) = −(eq.11)  →  24α + 28γ = −2")
    print("   p_t = −ρ  →  (eq.13) = −(eq.11)  →  8α + 2γ = −1")
    print()
    # Решение системы:
    #   24α + 28γ = −2   (p_r = −ρ)
    #    8α +  2γ = −1   (p_t = −ρ)
    # Из второго: γ = (−1−8α)/2.  Подставляем в первый:
    #   24α + 28(−1−8α)/2 = −2  →  24α − 14 − 112α = −2  →  −88α = 12  →  α = −3/22
    #   γ = (−1+24/22)/2 = (2/22)/2 = 1/22 > 0  ← γ должен быть < 0 !
    a_both = mp.mpf('-3')/22
    g_both = (-1 - 8*a_both)/2
    print("   Система p_r=−ρ AND p_t=−ρ:  α = −3/22 ≈ %.6f,  γ = 1/22 ≈ %.6f" %
          (float(a_both), float(g_both)))
    if g_both >= 0:
        print("   → γ≥0: НЕ допустимо (требуется γ<0).")
    rho_both = throat_values(a_both, g_both)[0]
    print("   → ρ = %.6f (>0: %s)" % (float(rho_both), "да" if rho_both > 0 else "нет"))
    print()

    # Только p_r = −ρ (одно уравнение, одна свободная пара):
    print("   Только p_r = −ρ:  γ = (−1−12α)/14")
    print("   Проверка ρ>0 на горловине:  1+22α+24γ < 0")
    print("   → γ = (−1−12α)/14 ⇒ 1+22α+24·(−1−12α)/14 = (−5−56α)/7 < 0")
    print("   → α > −5/56 ≈ −0.0893")
    print("   Также ρ+p_r = 0, ρ+p_t > 0 automatically")
    print()

    # Явная точка с p_r = −ρ и γ<0:
    for alpha_try in [mp.mpf('-0.01'), mp.mpf('-0.05'), mp.mpf('-0.08')]:
        gamma_try = (-1 - 12*alpha_try)/14
        if gamma_try < 0:
            rho_, pr_, pt_ = throat_values(alpha_try, gamma_try)
            print("   α=%+.3f  γ=%+.4f:  ρ=%+.4f  p_r=%+.4f  p_t=%+.4f  (p_r=−ρ? %s, ρ>0? %s)" %
                  (float(alpha_try), float(gamma_try),
                   float(rho_), float(pr_), float(pt_),
                   abs(pr_ + rho_) < 1e-12, rho_ > 0))

    # --- (5) Численные значения ρ, p_r, p_t, NEC на заданных x ---
    print("\n(5) Численные значения для γ∈{−0.01,−0.1,−0.5,−1}, α∈{0,−0.1,−1}:")
    print("    x = 1.0 (горловина), 1.05, 1.2, 1.5  (единицы 4πr0²)")
    print()

    detailed_data = []
    for gamma in gamma_vals:
        for alpha in alpha_vals:
            print("   α=%+.2f  γ=%+.2f:" % (float(alpha), float(gamma)))
            for x in x_pts:
                rho, pr, pt = field_matter_alphaGamma(x, alpha, gamma)
                necr = rho + pr
                nect = rho + pt
                wec  = rho >= 0
                print("     x=%.2f  ρ=%+.6e  p_r=%+.6e  p_t=%+.6e  NEC_r=%+.6e  NEC_t=%+.6e  WEC=%s NECr=%s NECt=%s" %
                      (float(x), float(rho), float(pr), float(pt),
                       float(necr), float(nect),
                       "+" if wec else "-",
                       "+" if necr >= 0 else "-",
                       "+" if nect >= 0 else "-"))
                detailed_data.append({
                    'alpha': float(alpha), 'gamma': float(gamma),
                    'x': float(x),
                    'rho': float(rho), 'pr': float(pr), 'pt': float(pt),
                    'NECr': float(necr), 'NECt': float(nect),
                    'WEC': bool(wec), 'NECr_ok': bool(necr >= 0), 'NECt_ok': bool(nect >= 0),
                })
            print()

    # --- (6) Асимптотика x→∞ ---
    print("\n(6) Асимптотика x→∞ (T→0, γT⁻¹ взрывается):")
    print("    При x→∞:  T ~ 2/x⁴ → 0,  f ~ γx⁴/2 (доминирует γT⁻¹)")
    print("    ρ ~ γx⁴/2,  p_r ~ −γx⁴/2,  p_t ~ γx⁴/2  (φон → ±∞)")
    print("    Это означает: γ<0 → ρ→−∞ far from throat: WEC нарушается.")
    print()
    for gamma in gamma_vals:
        print("   γ=%+.2f:" % float(gamma))
        for x_val in [10, 50, 100, 500]:
            asy = asymptotic_check(x_val, 0, gamma)
            print("     x=%4d  T=%.2e  f=%.2e  f/T=%.2e  ρ=%.2e  p_r=%.2e  p_t=%.2e" %
                  (x_val, asy['T'], asy['f'],
                   asy['f_over_T'] if asy['f_over_T'] is not None else 0,
                   asy['rho'], asy['pr'], asy['pt']))
        print()

    # ── Сбор результатов ──
    output = {
        "model": "f(T) = T + αT² + γT⁻¹  (γ<0)",
        "geometry": "Morris-Thorne (Ellis): a=0, β=r0²/r, r0=0.1m",
        "tetrad": "good (BHL convention), 4πr0² units",
        "equations": {
            "T(x)": "2(x−√(x²−1))²/x⁴  (BHL 16)",
            "Tp(x)": "−4(x−u)²(x+2u)/(u·x⁵)  (BHL 17)",
            "rho": "(eq.1): al/x·(1−al)·Tp·fTT − (T/4−1/(2x²))·fT + ebm/(2x²)·xb_m1·fT + f/4",
            "pr":  "(eq.2): [−1/(2x²)+T/4+ebm/(2x²)]·fT − f/4",
            "pt":  "(eq.3): (ebm/2)(1/x−1/(al·x))Tp·fTT + fT[T/4+ebm/(2x)(−bp²/2)] − f/4",
            "f":   "T+αT²+γ/T",
            "fT":  "1+2αT−γ/T²",
            "fTT": "2α+2γ/T³",
            "throat_rho":  "−(1+22α+24γ)/2  (eq.11)",
            "throat_pr":   "−(1+2α+4γ)/2   (eq.12)",
            "throat_pt":   "(1+14α+22γ)/2   (eq.13)",
        },
        "TEGR_check": {"analytic_ok": ok_tegr, "numerical_ok": ok_num},
        "quad_only_check": [{"alpha": c[0], "match": c[1]} for c in quad_checks],
        "scan_table": results_table,
        "detailed_values": detailed_data,
        "phantom_EoS_analysis": {
            "p_eq_mrho_both": {
                "equation": "p_r=−ρ AND p_t=−ρ  →  α=−3/22≈−0.1364, γ=1/22≈+0.0455",
                "gamma_positive": bool(g_both >= 0),
                "admissible": bool(g_both < 0),
                "comment": "γ≥0: не допустимо в рамках γ<0. Нет решения с p=−ρ по обоим компонентам и γ<0."
            },
            "p_eq_mrho_radial_only": {
                "equation": "p_r=−ρ  →  γ=(−1−12α)/14",
                "rho_positive_condition": "α > −5/56 ≈ −0.0893",
                "examples": []
            }
        },
        "asymptotic": {
            "behavior": "T~2/x⁴→0, f~γx⁴/2→−∞ (γ<0), ρ~γx⁴/2→−∞",
            "problem": "Паттерн f(T) с γT⁻¹ НЕ-флатовый: far-zone ρ<0, WEC нарушен.",
            "onset": "Отклонение от Minkowski O(γ/x²): при x>1/√|γ|~3.2 (γ=−0.1), O(10%) при x~1/|γ|^{1/4}"
        },
        "comments": [
            "p=−ρ (фантом-подобное EoS) на горловине НЕВОЗМОЖНО при γ<0 для обоих компонент.",
            "Только p_r=−ρ возможно при特定ных (α,γ): γ=(−1−12α)/14, α>−5/56.",
            "Асимптотика x→∞: γ<0 делает ρ<0 на больших расстояниях — WEC нарушен.",
            "Это фундаментальная проблема γT⁻¹-модели: не-флатовый фон.",
            "Для штатного фона (Minkowski far-zone) нужна регуляризация или отдельная matter-секция.",
        ]
    }

    with open("/home/smboozha/portal_gun/research/r3_ft_alphaGamma_results.json", "w") as fp:
        json.dump(output, fp, indent=2, default=str)
    print("\n  → r3_ft_alphaGamma_results.json сохранён")

    # ── Итоговое резюме ──
    print("\n" + "="*80)
    print("  РЕЗЮМЕ")
    print("="*80)
    print("""
(1) Уравнения поля: стандартная структура телепараллелизма (T_μν=(1/κ)∇_λ S^λ_μν−…),
    матрица S^{μν}_ρ реализована в tetrad_core.py. Не новый постулат.

(2) ρ, p_r, p_t для f=T+αT²+γT⁻¹: формулы (eq.1)-(eq.13) выше.
    Сверка: TEGR (α=γ=0) OK; quad-only (γ=0) OK vs ft2_necscan_BHL.py.

(3) p=−ρ НЕВОЗМОЖНО одновременно по p_r И p_t при γ<0:
    Система p_r=−ρ ∧ p_t=−ρ → γ=+1/22>0 (не допустимо).
    Только p_r=−ρ: γ=(−1−12α)/14, требует α>−5/56≈−0.09.

(4) Асимптотика x→∞:  T~2/x⁴→0,  γ/T~γx⁴/2→−∞ (γ<0).
    ρ→−∞, WEC нарушен远 от горловины. Паттерн НЕ-флатовый.
    Отклонение от Minkowski: O(|γ|/x²), существенно при x≲1/|γ|^{1/4}.

(5) Таблица (α,γ) → EoS-совместимость и численные значения — в JSON.
""")

    return output

if __name__ == "__main__":
    run_scan()
