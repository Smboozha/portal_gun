"""
R4 [АГЕНТ-4, ИНЖЕНЕР-ФИЗИК] — «Портальный картридж»: сверхпроводящее кольцо
================================================================================
R=0.05 м, I=1000 А (персистентный ток, YBCO-ВТСП).

ЧАСТЬ 1: РЕАЛЬНЫЕ параметры кольца (магнитное поле, энергия, давление,
         самоиндукция, критический ток YBCO, поток в квантах, фаза AB).
ЧАСТЬ 2: ГЛАВНЫЙ ВОПРОС — может ли ток кольца создать кручение (torsion)
         пространства-времени?
ЧАСТЬ 3: Таблица + JSON.

Конвенции проекта: r0=0.1 м (ft2_necscan_BHL.py), α=T0·r0² (r4_qi_ft.py).
Порог снятия QI: α≤−1/22 → T0_crit=−1/(22·r0²) ≈ −4.545 м⁻².

[РЕАЛЬНО] — всё основано на стандартной физике (SI-единицы).
[ПОСТУЛАТ] — явно помечены гипотезы о связи tor↔EM.
"""
import json
import numpy as np

# ================================ КОНСТАНТЫ (SI) =============================
mu0 = 4.0 * np.pi * 1e-7          # Тл·м/А
hbar = 1.054571817e-34             # Дж·с
e_charge = 1.602176634e-19          # Кл
c = 2.99792458e8                    # м/с
G = 6.67430e-11                     # м³/(кг·с²)
pi = np.pi
r0 = 0.10                           # м — конвенция проекта

# ================================ ПАРАМЕТРЫ КОЛЬЦА ===========================
R = 0.05                            # м — радиус кольца
I = 1000.0                          # А — персистентный ток
a_w = 1.0e-3                        # м — радиус проволоки (1 мм)

# ======================== 1. МАГНИТНОЕ ПОЛЕ (РЕАЛЬНО) ========================

def B_center(mu0, I, R):
    """[1] Магнитное поле в центре кольца: B_0 = μ0·I/(2R) [Тл]."""
    return mu0 * I / (2.0 * R)

def B_axis(mu0, I, R, z):
    """[2] Магнитное поле на оси: B(z) = μ0·I·R² / (2·(R²+z²)^{3/2}) [Тл]."""
    return mu0 * I * R**2 / (2.0 * (R**2 + z**2)**1.5)

def self_inductance_ring(mu0, R, a_w):
    """[3] Самоиндукция кольца: L ≈ μ0·R·[ln(8R/a_w) − 2] [Гн].
    Аппроксимация для тонкого кольца (a_w << R)."""
    return mu0 * R * (np.log(8.0 * R / a_w) - 2.0)

def magnetic_energy(L, I):
    """[4] Запасённая энергия: E_mag = ½·L·I² [Дж]."""
    return 0.5 * L * I**2

def magnetic_pressure(B):
    """[5] Магнитное давление: p_B = B²/(2μ0) [Па]."""
    return B**2 / (2.0 * mu0)

# ======================== 2. ВТСП: YBCO (РЕАЛЬНО) ===========================
# YBCO (YBa₂Cu₃O₇): T_c ≈ 92 К, при T=77 К (жидкий азот):
#   B_c2(77 К) ≈ 1.5–2 Тл (верхнее критическое поле)
#   J_c(77 К, B=0) ≈ 10⁴–10⁵ А/см² для монокристаллов
#   J_c(77 К, B≈0.01 Тл) ≈ 10³–10⁴ А/см² (типично для ленточных проводников)

T_c_YBCO = 92.0                     # К — критическая температура
T_op = 77.0                          # К — рабочая температура (жидкий азот)
B_c2_77 = 1.5                        # Тл — верхнее критическое поле при 77 К (типичное)
J_c_typical = 3.0e3                  # А/см² — типичное J_c при 77 К, B~0.01 Тл (ленточный проводник)
J_c_optimistic = 1.0e4               # А/см² — хороший J_c при 77 К, B→0 ( SuAM/SuperPower )

def ybco_cross_section(I, J_c):
    """[6] Минимальное сечение проводника: A = I / J_c [м²]."""
    return I / J_c

# ======================== 3. МАГНИТНЫЙ ПОТОК И AB (РЕАЛЬНО) ==================

Phi_0 = hbar / (2.0 * e_charge)     # ≈ 3.29×10⁻¹⁶ Вб — ħ/(2e)
Phi_0_h = 6.62607015e-34 / (2.0 * e_charge)   # ≈ 2.07×10⁻¹⁵ Вб — h/(2e)
# ПРИМЕЧАНИЕ (честно): стандартный сверхпроводниковый квант потока (fluxoid)
# = h/2e = 2.07e-15 Вб. Задание просит «Φ0=ħ/2e≈2.07e-15» — это НЕСОГЛАСОВАННО:
# ħ/2e = 3.29e-16, а 2.07e-15 — это именно h/2e. Ниже даны ОБА числа: физический
# fluxoid N(Φ0=h/2e) и вариант с ħ/2e из задания. Разница — множитель 2π.

def flux_ring(B, R):
    """[7] Полный магнитный поток через кольцо: Φ = B·πR² [Вб]."""
    return B * pi * R**2

def flux_quanta_ratio(Phi):
    """[8] Число квантов: N = Φ/Φ0."""
    return Phi / Phi_0

def ab_phase(Phi):
    """[9] Фаза Ааронова–Бома: φ = eΦ/ħ = 2π·N [рад]."""
    return e_charge * Phi / hbar

# ======================== 4. TORSION: ГИПОТЕЗА [ПОСТУЛАТ] ===================
# [ПОСТУЛАТ]: если бы связь «ток → кручение» была ξ·T0, с T0 ∝ B·(e/ħ),
# то T0 = B·(e/ħ) (в ед. м⁻²). Это ЧИСТАЯ ОЦЕНКА ПОРЯДКА.

def T0_hypothesis(B):
    """[10] [ПОСТУЛАТ] T0 = B·(e/ħ) [м⁻²] — если бы EM-поле давало torsion.
    Ниже сравнивается с T0_crit из r4_qi_ft."""
    return B * e_charge / hbar

# Конвенция проекта: α = T0·r0² (безразмерно), см. r4_qi_ft.py, r4_bubound.py
# Порог снятия QI: q = −22α; q≥1 → ρ_req≥0, экзотика снята.
# → α_crit = −1/22 ≈ −0.04545
# → T0_crit = α_crit / r0² = −1/(22·r0²) = −4.545 м⁻²

alpha_crit = -1.0 / 22.0            # ≈ −0.04545
T0_crit = alpha_crit / r0**2         # ≈ −4.545 м⁻²

# ======================== 5. ГРАВИТАЦИЯ: Schwarzschild-поправки ==============

def mass_equivalent(E_mag):
    """[11] Гравитационный эквивалент энергии: m = E/c² [кг].
    Это НЕ масса кольца, а масса-эквивалент энергии поля."""
    return E_mag / c**2

def schwarzschild_radius(m):
    """[12] Радиус Шварцшильда: r_s = 2Gm/c² [м]."""
    return 2.0 * G * m / c**2

def tidal_accel(m, r, dl):
    """[13] Приливное ускорение на расстоянии r: δa ≈ 2Gm·dl/r³ [м/с²].
    dl — размер тестового тела."""
    return 2.0 * G * m * dl / r**3

# ============================= ГЛАВНЫЙ РАСЧЁТ ================================

def main():
    print("=" * 110)
    print("  R4 [АГЕНТ-4]  «ПОРТАЛЬНЫЙ КАРТРИДЖ» — сверхпроводящее кольцо R=0.05 м, I=1000 А")
    print("  Конвенция: r0=0.1 м (ft2_necscan_BHL.py); α=T0·r0² (r4_qi_ft.py)")
    print("=" * 110)

    # ---- 1. Магнитное поле ----
    B0 = B_center(mu0, I, R)
    print("\n[1] МАГНИТНОЕ ПОЛЕ В ЦЕНТРЕ")
    print("    B_0 = μ0·I/(2R) = %.4e Тл = %.2f мТл" % (B0, B0 * 1e3))

    print("\n[2] ПОЛЕ НА ОСИ B(z):")
    for z_m in [0.0, 0.05, 0.10, 0.20, 0.50]:
        Bz = B_axis(mu0, I, R, z_m)
        print("    z = %.2f м: B = %.4e Тл = %.4f мТл" % (z_m, Bz, Bz * 1e3))

    # ---- 2. Самоиндукция и энергия ----
    L = self_inductance_ring(mu0, R, a_w)
    E = magnetic_energy(L, I)
    p_B = magnetic_pressure(B0)
    print("\n[3] САМОИНДУКЦИЯ")
    print("    L ≈ μ0·R·[ln(8R/a_w)−2] = %.4e Гн" % L)
    print("    ln(8R/a_w) = %.4f" % np.log(8.0 * R / a_w))

    print("\n[4] ЗАПАСЁННАЯ ЭНЕРГИЯ")
    print("    E_mag = ½·L·I² = %.4e Дж = %.4e МДж" % (E, E / 1e6))

    m_eq = mass_equivalent(E)
    r_s = schwarzschild_radius(m_eq)
    print("    m_eq = E/c² = %.4e кг" % m_eq)
    print("    r_s = 2Gm/c² = %.4e м (≈ %.2e R)" % (r_s, r_s / R))

    print("\n[5] МАГНИТНОЕ ДАВЛЕНИЕ")
    print("    p_B = B²/(2μ0) = %.4e Па = %.2f атм" % (p_B, p_B / 101325.0))

    # ---- 3. YBCO: критический ток ----
    print("\n[6] ВТСП: YBCO ПРИ 77 К (жидкий азот)")
    print("    T_c = %.0f К,  T_раб = %.0f К,  B_c2(77 К) ≈ %.1f Тл" %
          (T_c_YBCO, T_op, B_c2_77))
    print("    B_0 = %.4f мТл ≪ B_c2 = %.1f Тл → кольцо в пределе критического тока" %
          (B0 * 1e3, B_c2_77))

    A_typ = ybco_cross_section(I, J_c_typical)
    A_opt = ybco_cross_section(I, J_c_optimistic)
    d_typ = np.sqrt(A_typ / pi)      # радиус эквивалентного проводника
    d_opt = np.sqrt(A_opt / pi)
    print("    J_c(типичное, 77К, B~0.01Тл) ≈ %g А/см² → A_min = %.0f мм² (∅~%.0f мм)"
          % (J_c_typical, A_typ * 1e4, 2 * d_typ * 1e3))
    print("    J_c(хорошее, 77К, B→0) ≈ %g А/см² → A_min = %.0f мм² (∅~%.0f мм)"
          % (J_c_optimistic, A_opt * 1e4, 2 * d_opt * 1e3))
    print("    Вывод: для I=1000 А нужно суммарное сечение сверхпроводника ~%.0f–%.0f мм²." %
          (A_typ * 1e4, A_opt * 1e4))
    print("    Это РЕАЛИЗУЕМО как СТОПКА (не один) параллельных YBCO-лент.")
    print("    Термостабильность: B_0=%.1f мТл ≪ B_c2=%.1f Тл — глубокий сверхпроводящий"
          % (B0 * 1e3, B_c2_77))
    print("    режим; p_B≈%.0f Па (%.1e атм) — механическая нагрузка невелика." % (p_B, p_B / 101325.0))

    # ---- 4. Магнитный поток и фаза AB ----
    Phi = flux_ring(B0, R)
    phi_AB = ab_phase(Phi)
    N_hbar2e = Phi / Phi_0          # с ħ/2e (из задания)
    N_h2e = Phi / Phi_0_h           # с физическим fluxoid h/2e
    print("\n[7] МАГНИТНЫЙ ПОТОК")
    print("    Φ = B_0·πR² = %.4e Вб" % Phi)
    print("    Φ₀,fluxoid = h/2e    = %.4e Вб  (стандартный сверхпроводниковый квант)" % Phi_0_h)
    print("    Φ₀,ħ/2e     = ħ/2e   = %.4e Вб  (вариант из задания)" % Phi_0)
    print("    N = Φ/(h/2e)        = %.4e fluxoid квантов   <-- физический ответ" % N_h2e)
    print("    N = Φ/(ħ/2e)        = %.4e  (согласованно с φ_AB=eΦ/ħ)" % N_hbar2e)

    print("\n[8] ФАЗА ААРОНОВА–БОМА")
    print("    φ_AB = eΦ/ħ = 2πN(ħ/2e) = %.4e рад = %.4e·2π" % (phi_AB, N_hbar2e))
    print("    Полное число fluxoid-квантов: ~%.0f (h/2e); фаза AB = %d·2π (ħ/2e)." %
          (N_h2e, N_hbar2e))
    print("    Физика: кольцо с персистентным током НЕ создаёт кручение;")
    print("    фаза AB — квантовомеханический эффект на заряженные частицы,")
    print("    а не геометрический torsion-ток.")

    # ---- 5. Torsion-гипотеза [ПОСТУЛАТ] ----
    T0_hyp = T0_hypothesis(B0)
    alpha_hyp = T0_hyp * r0**2
    print("\n" + "=" * 110)
    print("[ПОСТУЛАТ] СВЯЗЬ «EM → TORSION»: T0 = B·(e/ħ)")
    print("=" * 110)
    print("    T0(гипотеза) = B_0·(e/ħ) = %.4e м⁻²" % T0_hyp)
    print("    α(гипотеза) = T0·r0² = %.4e" % alpha_hyp)
    print()
    print("    T0_crit (из r4_qi_ft: α≤−1/22) = %.4e м⁻²" % T0_crit)
    print("    α_crit = −1/22 = %.6f" % alpha_crit)
    print()
    print("    СРАВНЕНИЕ:")
    print("    |T0(гипотеза)| / |T0_crit| = %.4e" % abs(T0_hyp / T0_crit))
    print("    α(гипотеза) / α_crit = %.4e" % abs(alpha_hyp / alpha_crit))
    print()
    if T0_hyp > 0 and T0_crit < 0:
        print("    знак: T0(гипотеза) > 0, T0_crit < 0 — даже знак НЕ совпадает!")
        print("    (T0_crit < 0 required для q=−22α>0; T0_hyp>0 даёт α>0 → q<0 → ужесточение)")
    print("    Порядок: гипотеза даёт T0 ~ %.1e м⁻², порог ~ %.1f м⁻²." %
          (T0_hyp, abs(T0_crit)))
    print("    → ГИПОТЕЗА «B·e/ħ» ЗАПРЕЩЁНА (дает α>0, ужесточает QI вместо снятия).")

    # ---- 6. Schwarzschild-поправки (.Gravity, НЕ torsion) ----
    print("\n" + "-" * 110)
    print("[РЕАЛЬНО] ГРАВИТАЦИОННЫЙ ЭФФЕКТ ЭНЕРГИИ ПОЛЯ (Schwarzschild):")
    print("    E_mag/c² = %.4e кг → r_s = %.4e м" % (m_eq, r_s))
    print("    a_тилд(R) = Gm/R² = %.4e м/с² (attractive)" % (G * m_eq / R**2))
    print("    Приливное ускорение на δl=0.05 м, r=0.1 м: δa = %.4e м/с²" %
          tidal_accel(m_eq, 0.1, 0.05))
    print("    → Гравитационный эффект энергии поля: r_s/R ~ %.1e (≈0)." % (r_s / R))
    print("    Schwarzschild-поправки ЕСТЬ, но ОЧЕНЬ малы; это КРИВИЗНА (g_μν),")
    print("    а НЕ телепораллельное кручение (не торсион).")

    # ---- 7. ИТОГОВАЯ ТАБЛИЦА ----
    print("\n" + "=" * 110)
    print("  ИТОГОВАЯ ТАБЛИЦА ПАРАМЕТРОВ КОЛЬЦА")
    print("=" * 110)
    hdr = ("  Параметр".ljust(45) + "Значение".rjust(25) + "  Примечание")
    print(hdr)
    print("  " + "-" * 106)
    rows = [
        ("B_0 (центр)",                      "%.4e Тл" % B0,               "[РЕАЛЬНО] μ0I/(2R)"),
        ("B(z=0.05м) = B(R)",                "%.4e Тл" % B_axis(mu0,I,R,R),"[РЕАЛЬНО]"),
        ("B(z=0.20м)",                       "%.4e Тл" % B_axis(mu0,I,R,0.2),"[РЕАЛЬНО]"),
        ("L (самоиндукция)",                  "%.4e Гн" % L,                "[РЕАЛЬНО] ln(8R/a_w)−2"),
        ("E_mag (энергия)",                   "%.4e Дж" % E,                "[РЕАЛЬНО] ½LI²"),
        ("p_B (магн. давление)",              "%.4e Па" % p_B,              "[РЕАЛЬНО] B²/(2μ0)"),
        ("Φ (поток)",                         "%.4e Вб" % Phi,              "[РЕАЛЬНО] B·πR²"),
        ("Φ/(h/2e) fluxoid квантов",          "%.4e" % N_h2e,               "[РЕАЛЬНО] физич. N"),
        ("Φ/(ħ/2e) (вариант задания)",        "%.4e" % N_hbar2e,            "[РЕАЛЬНО] =φ_AB/2π"),
        ("φ_AB (фаза AB)",                    "%.4e рад" % phi_AB,          "[РЕАЛЬНО] eΦ/ħ"),
        ("YBCO: сечение min (типичн.)",       "%.0f мм²" % (A_typ*1e4),    "[РЕАЛЬНО] J_c=3kA/cm²"),
        ("YBCO: сечение min (хорошее)",       "%.0f мм²" % (A_opt*1e4),    "[РЕАЛЬНО] J_c=10kA/cm²"),
        ("m_eq (гравит. эквивалент)",        "%.4e кг" % m_eq,            "[РЕАЛЬНО] E/c²"),
        ("r_s (Шварцшильд)",                 "%.4e м" % r_s,              "[РЕАЛЬНО] 2Gm/c²"),
        ("r_s / R",                            "%.1e" % (r_s/R),           "[РЕАЛЬНО] ~10⁻⁴⁴"),
        ("T0(гипотеза B·e/ħ)",               "%.4e м⁻²" % T0_hyp,         "[ПОСТУЛАТ]"),
        ("α(гипотеза)",                       "%.4e" % alpha_hyp,          "[ПОСТУЛАТ] T0·r0²"),
        ("T0_crit (порог QI)",                "%.4e м⁻²" % T0_crit,        "[РЕАЛЬНО] −1/(22r0²)"),
        ("α_crit (порог QI)",                 "%.6f" % alpha_crit,         "[РЕАЛЬНО] −1/22"),
        ("|T0_hyp| / |T0_crit|",             "%.1e" % abs(T0_hyp/T0_crit),"[СРАВНЕНИЕ]"),
    ]
    for name, val, note in rows:
        print("  %-44s %25s  %s" % (name, val, note))

    print("\n" + "-" * 110)
    print("  ЧЕСТНЫЙ ВЕРДИКТ ПО МЕХАНИЗМУ «ТОК → КРУЧЕНИЕ → УДЕРЖАНИЕ ГОРЛОВИНЫ»")
    print("-" * 110)
    print("""
  ФИЗИКА МАГНИТНОГО ПОЛЯ:
    ЭМ-поле кольца создаёт T_μν ≠ 0 (тензор энергии-импульса).
    По Эйнштейну, T_μν ИСТОЧНИК КРИВИЗНЫ (g_μν): Schwarzschild-подобные поправки.
    r_s/R ~ 10⁻⁴² → кривизна ЕСТЬ, но практически нулевая.

  ФИЗИКА TORSION (TЕЛЕПАРАЛЛЕЛЬНaЯ):
    В телепараллельной гравитации (TEGR/f(T)) кручение — это фоновая
    неголономная реперная структура (tetrad e^a_μ), а НЕ результат
    электромагнитного поля как такового. torsion T^ρ_μν связан с геометрией
    реперов, не с T_μν материи напрямую.

  ВЫВОД (ЧЕСТНЫЙ):
    [РЕАЛЬНО] Магнитное поле кольца создаёт КРИВИЗНУ (Schwarzschild),
    а НЕ torsion. Нет установленного механизма «ток I → torsion T0».
    [ПОСТУЛАТ] Даже если T0 = B·(e/ħ) (квантовый масштаб):
      — T0_hyp ~ %.1e м⁻², T0_crit = −%.2f м⁻².
      — |T0_hyp|/|T0_crit| ~ %.1e: гипотетический torsion ~10¹²·(2π) раз
        ВЫШЕ критического (но с НЕВЕРНЫМ ЗНАКОМ: α>0 → q<0 → ужесточает QI).
      — даже если бы знак совпал, необходимая величина кручения не связана
        с физикой кольца: нужен отдельный (неконвенциональный) источник.

    УДЕРЖАНИЕ ГОРЛОВИНЫ МАГНИТНЫМ КОЛЬЦОМ (R=0.05 м, I=1000 А):
    — НЕ ДОКАЗАНО и НЕ РЕАЛИЗУЕМО текущей физикой.
    — Schwarzschild-поправка: r_s/R ~ 4×10⁻⁴⁴ (нерелевантна).
    — Torsion из EM: НЕТ устоявшегося механизма; гипотеза запрещена знаком.
    — Для поддержки горловины в f(T) нужен α<0 (T0<0), что требует
      ОТДЕЛЬНОГО источника кручения (не ЭМ-поля).
""" % (T0_hyp, abs(T0_crit), abs(T0_hyp / T0_crit)))

    # ---- 8. JSON ----
    Phi_json = float(Phi)
    phi_AB_json = float(phi_AB)
    out = {
        "dataset": "r3_ring_cartridge",
        "ring": {"R_m": R, "I_A": I, "a_w_m": a_w},
        "B0_T": float(B0),
        "B_on_axis": {str(z): float(B_axis(mu0, I, R, z))
                      for z in [0.0, 0.05, 0.1, 0.2, 0.5]},
        "L_H": float(L),
        "E_mag_J": float(E),
        "p_B_Pa": float(p_B),
        "p_B_atm": float(p_B / 101325.0),
        "YBCO": {
            "T_c_K": T_c_YBCO, "T_op_K": T_op, "B_c2_77K_T": B_c2_77,
            "J_c_typical_A_cm2": J_c_typical,
            "J_c_optimistic_A_cm2": J_c_optimistic,
            "A_min_typical_mm2": float(A_typ * 1e4),
            "A_min_optimistic_mm2": float(A_opt * 1e4),
            "verdict": "realizuem kak stаck parallellnykh YBCO-lent; B<<B_c2, glubokii SC-rezhim"
        },
        "flux": {
            "Phi_Wb": Phi_json,
            "Phi0_h2e_Wb": float(Phi_0_h),
            "Phi0_hbar2e_Wb": float(Phi_0),
            "N_fluxoid_h2e": float(N_h2e),
            "N_hbar2e": float(N_hbar2e),
            "phase_AB_rad": phi_AB_json,
            "phase_AB_times_2pi": float(N_hbar2e),
            "note": "fizicheskii fluxoid = h/2e = 2.07e-15 Wb; variant zadaniya 'hbar/2e' "
                    "est' fakticheski h/2e (oshibka 2pi); N(h2e)=4.77e10, N(hbar2e)=2.999e11"
        },
        "gravitational": {
            "m_eq_kg": float(m_eq),
            "r_s_m": float(r_s),
            "r_s_over_R": float(r_s / R),
            "verdict": "Schwarzschild-korrektciia EST' no ~10^-42 ot R, nerelevantna"
        },
        "torsion_hypothesis": {
            "label": "[POSTULAT] T0 = B*(e/hbar)",
            "T0_hypothesis_m2": float(T0_hyp),
            "alpha_hypothesis": float(alpha_hyp),
            "T0_crit_m2": float(T0_crit),
            "alpha_crit": alpha_crit,
            "ratio_T0_hyp_over_T0_crit": float(abs(T0_hyp / T0_crit)),
            "sign_mismatch": True,
            "verdict": "ZAPRESHCHENO: alpha>0 → q<0 → ushchestochet QI vmesto sniatiia; "
                       "|T0_hyp|/|T0_crit| ~ 4.2e12 (nevernyi znak + nevernyi masshtab)"
        },
        "verdict_final": (
            "UDERZHANIE GORLOVINY MAGNITNYM KOLTSOM: NE DOKAZANO. "
            "Mag. pole sozdaet kriviznu (Schwarzschild ~ 4e-44 ot R), a NE torsion. "
            "Mekhanizm 'tok → torsion' neustanovlen; gipoteza T0=B*(e/hbar) "
            "zapreshchena znakom (alpha>0) i masshtabom (~4.2e12). "
            "Dlya f(T)-garpunirovaniia nuzhen OTDELNYI istochnik T0<0."
        ),
    }
    with open("/home/smboozha/portal_gun/research/r3_ring_cartridge_values.json", "w",
              encoding="utf-8") as fp:
        json.dump(out, fp, ensure_ascii=False, indent=2, default=str)

    print("  JSON: research/r3_ring_cartridge_values.json сохранён")
    print("=" * 110)

if __name__ == "__main__":
    main()
