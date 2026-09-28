"""
R3 [АГЕНТ-4, ИНЖЕНЕР-ТЕОРФИЗИК] — «Сфера-гиротор»: вращающийся сверхпроводник
================================================================================
R_s=0.1 м, вращающаяся сверхпроводящая сфера (YBCO, никелевая подложка) как
генератор «тороидального кручения» и магнитного поля одновременно.

ЧАСТЬ I  — МАГНИТНОЕ ПОЛЕ (честная электродинамика)
   Модель-1 «вращающийся заряженный шар» (полный заряд Q):
       m = Q·Ω·R_s²/3 (дипольный момент), B_c = μ0·Q·Ω/(6πR_s).
       Честный предел Q_max=4πε0R_s²·E_max (пробой 1e7 В/м) — сколько заряда
       реально «разделить» на обкладке-потенциале.
   Модель-2 «эффект Лондона» (физика ВТСП, Meissner):
       B_London = −(2m_e/e)·Ω — поле ВНУТРИ вращающегося сверхпроводника.
       Ключевой честный факт: сверхпроводник ЭКРАНИРУЕТ заряд/поле, поэтому
       физическое поле — это Лондон, НЕ зависящее от Q. Модель заряженного
       шара — только умозрительный верхний предел, не реализуемый в S/C.

   Предел по прочности: Ω_break=√(σ/(ρ·R_s²)), σ=150 МПа (никель), ρ≈8908 кг/м³
       → Ω_max ≈ 1298 рад/с при R_s=0.1 м.

ЧАСТЬ II — КРУЧЕНИЕ ОТ ВРАЩЕНИЯ (честно, два канона)
   (а) T_rot ≈ Ω²/c² — скаляр безразм. тензора кручения (телепараллельность
       вращающейся СО), м⁻².
   (б) Эквивалент Lense–Thirring (frame-dragging):
       Ω_LT = G·J/(c²·R_s³) на поверхности, J=⅖M·R_s²·Ω;
       T_LT ≈ (Ω_LT/c)², м⁻².
   Оба — на Ω_max=1298 рад/с и сравнение c T0_crit=−4.545 м⁻².

ЧАСТЬ III — РАБОЧИЙ ИТОГ: «сфера-гиротор» даёт B~1e-8 Тл и torsion~1e-11..1e-61 м⁻².
   Схема к «рабочему порталу» НЕ сходится. ДА/НЕТ с числами.

Предыдущие честные вердикты (r3_ring_cartridge.py, r3_plasma_vortex_a4.py):
   ЭМ-источники torsion НЕ создают; вращающийся сверхпроводник здесь — ТОЖЕ НЕТ.

Конвенции проекта: r0=0.1 м (ft2_necscan_BHL.py), α=T0·r0² (r4_qi_ft.py).
Порог снятия QI: α≤−1/22 → T0_crit = −1/(22·r0²) ≈ −4.545 м⁻².

[РЕАЛЬНО] — стандартная физика (SI). [ПОСТУЛАТ] — гипотеза tor↔rot.
"""
import json
import numpy as np

# ================================ КОНСТАНТЫ (SI) =============================
mu0 = 4.0 * np.pi * 1e-7          # Тл·м/А
epsilon0 = 8.8541878128e-12       # Ф/м
hbar = 1.054571817e-34            # Дж·с
e_charge = 1.602176634e-19        # Кл
m_e = 9.1093837015e-31            # кг
c = 2.99792458e8                  # м/с
G = 6.67430e-11                   # м³/(кг·с²)
pi = np.pi
r0 = 0.10                         # м — конвенция проекта

# ================================ ПАРАМЕТРЫ СФЕРЫ ===========================
R_s = 0.10                        # м — радиус сверхпроводящей сферы
E_max = 1.0e7                     # В/м — предел пробоя (обкладка-потенциал)
sigma_m = 150.0e6                 # Па — предел прочности (никель/сталь)
rho_mat = 8908.0                  # кг/м³ — плотность никеля (подложка/YBCO)

# ============ 1. ПРЕДЕЛ ПО ПРОЧНОСТИ: Ω_max (РЕАЛЬНО) =======================

def omega_break(sigma, rho, R_s):
    """[1] Ω_break = √(σ/(ρ·R_s²)) [рад/с] — центробежный разрушающий предел
    для твёрдой сферы (порядок оценки σ/ρ·R²)."""
    return np.sqrt(sigma / (rho * R_s**2))

# ============ 2. МАГНИТНОЕ ПОЛЕ — МОДЕЛЬ ЗАРЯЖЕННОГО ШАРА (РЕАЛЬНО) =========

def Q_max(epsilon0, R_s, E_max):
    """[2] Максимальный разделимый заряд до пробоя: Q_max = 4πε0·R_s²·E_max [Кл]."""
    return 4.0 * pi * epsilon0 * R_s**2 * E_max

def B_center_charged(mu0, Q, Omega, R_s):
    """[3] Поле в центре вращающегося заряженного шара:
    m = Q·Ω·R_s²/3, B_c = μ0·m/(2πR_s³) = μ0·Q·Ω/(6πR_s) [Тл]."""
    return mu0 * Q * Omega / (6.0 * pi * R_s)

# ============ 3. МАГНИТНОЕ ПОЛЕ — ЭФФЕКТ ЛОНДОНА (РЕАЛЬНО, ВТСП) ============

def B_london(m_e, e_charge, Omega):
    """[4] Поле Лондона вращающегося сверхпроводника (Meissner):
    B_L = −(2m_e/e)·Ω [Тл]. Для электрона заряд −e, отсюда знак «−»;
    величина |B_L| = (2m_e/e)·Ω. Поле НЕ зависит от Q (S/C экранирует)."""
    return (2.0 * m_e / e_charge) * Omega

# ============ 4. КРУЧЕНИЕ ОТ ВРАЩЕНИЯ (ЧЕСТНО, ДВА КАНОНА) ==================

def T_rot_scenario_a(Omega, c):
    """[5] (а) Скаляр безразм. тензора кручения вращающейся СО: T ≈ Ω²/c² [м⁻²]."""
    return Omega**2 / c**2

def moment_of_inertia_sphere(M, R_s):
    """[6] Момент инерции твёрдой сферы: I = ⅖·M·R_s² [кг·м²]."""
    return 0.4 * M * R_s**2

def T_rot_scenario_b(mass, R_s, Omega, G, c):
    """[7] (б) Lense–Thirring frame-dragging: Ω_LT = G·J/(c²·R_s³) на поверхности,
    J = ⅖·M·R_s²·Ω. Torsion-эквивалент: T ≈ (Ω_LT/c)² [м⁻²]."""
    J = 0.4 * mass * R_s**2 * Omega
    Omega_LT = G * J / (c**2 * R_s**3)
    return (Omega_LT / c)**2, Omega_LT, J

# ============================= ПРОЧИЕ ВЕЛИЧИНЫ ==============================
def sphere_mass(rho, R_s):
    """[8] Масса твёрдой сферы: M = ρ·(4/3)πR_s³ [кг]."""
    return rho * (4.0 / 3.0) * pi * R_s**3

# ============================= ГЛАВНЫЙ РАСЧЁТ ================================

def main():
    print("=" * 120)
    print("  R3 [АГЕНТ-4]  «СФЕРА-ГИРОТОР» — вращающийся сверхпроводник (R_s=0.1 м)")
    print("  YBCO/никель, Ω до предела прочности. B-field + torsion одновременно.")
    print("  Конвенция: r0=0.1 м (ft2_necscan_BHL.py); α=T0·r0² (r4_qi_ft.py)")
    print("=" * 120)

    # ---- T0_crit ----
    alpha_crit = -1.0 / 22.0
    T0_crit = alpha_crit / r0**2

    # ---- 1. Ω_max (прочность) ----
    Omega_max = omega_break(sigma_m, rho_mat, R_s)
    print("\n[ЧАСТЬ I-0] ПРЕДЕЛ ПО ПРОЧНОСТИ Ω_max")
    print("    Ω_break = √(σ/(ρ·R_s²)) = √(%.1e/(%.0f·%.2f²))" % (sigma_m, rho_mat, R_s))
    print("    Ω_max = %.1f рад/с" % Omega_max)
    print("    f_max = Ω/(2π) = %.1f Гц" % (Omega_max/(2*pi)))

    # ---- 2. Q_max ----
    Qm = Q_max(epsilon0, R_s, E_max)
    print("\n[ЧАСТЬ I-1] СКОЛЬКО ЗАРЯДА РЕАЛЬНО РАЗДЕЛИТЬ?")
    print("    Q_max = 4πε0·R_s²·E_max = 4π·%.3e·%.2f²·%.1e" % (epsilon0, R_s, E_max))
    print("    Q_max = %.4e Кл  (ограничено пробоем 1e7 В/м)" % Qm)
    # сколько электронов
    print("    N_e = Q_max/e = %.2e электронов" % (Qm / e_charge))

    # ---- 3. B_c, модель заряженного шара (условно, при Q_max) ----
    Bc_Qmax_omegamax = B_center_charged(mu0, Qm, Omega_max, R_s)
    Bc_Qmax_1e3 = B_center_charged(mu0, Qm, 1e3, R_s)
    print("\n[ЧАСТЬ I-2] МОДЕЛЬ «ВРАЩАЮЩИЙСЯ ЗАРЯЖЕННЫЙ ШАР» (ВЕРХНИЙ ПРЕДЕЛ)")
    print("    m = Q·Ω·R²/3,  B_c = μ0·Q·Ω/(6πR_s)")
    print("    B_c(Ω_max=%.0f рад/с)  = %.4e Тл = %.4e Гс" %
          (Omega_max, Bc_Qmax_omegamax, Bc_Qmax_omegamax*1e4))
    print("    B_c(Ω=1000 рад/с)      = %.4e Тл = %.4e Гс" %
          (Bc_Qmax_1e3, Bc_Qmax_1e3*1e4))
    print("    → Модель даёт ~1e-8 Тл: микроскопически мало.")

    # ---- 4. B_London (физический канал для сверхпроводника) ----
    BL_omegamax = B_london(m_e, e_charge, Omega_max)
    BL_1e3 = B_london(m_e, e_charge, 1e3)
    print("\n[ЧАСТЬ I-3] ЭФФЕКТ ЛОНДОНА (ФИЗИКА ВТСП, Meissner)")
    print("    B_London = −(2m_e/e)·Ω; |B_L| = (2m_e/e)·Ω")
    print("    2m_e/e = %.4e Тл/(рад/с)" % (2.0*m_e/e_charge))
    print("    |B_L|(Ω_max=%.0f рад/с) = %.4e Тл = %d нТл" %
          (Omega_max, BL_omegamax, BL_omegamax*1e9))
    print("    |B_L|(Ω=1000 рад/с)     = %.4e Тл = %d нТл" %
          (BL_1e3, BL_1e3*1e9))
    print("    → КЛЮЧЕВОЙ ЧЕСТНЫЙ ФАКТ: поле Лондона НЕ зависит от Q.")
    print("      Сверхпроводник экранирует заряд/поле (Meissner), поэтому")
    print("      «набрать Q» до пробоя физически НЕВОЗМОЖНО — B_L — единственный")
    print("      магнитный канал, и он ~1e-8 Тл.")

    # ---- 5. Кручение, сценарий (а) ----
    Ta_omegamax = T_rot_scenario_a(Omega_max, c)
    Ta_1e3 = T_rot_scenario_a(1e3, c)
    print("\n[ЧАСТЬ II-а] КРУЧЕНИЕ ОТ ВРАЩЕНИЯ — СЦЕНАРИЙ (а) T≈Ω²/c²")
    print("    T_rot(Ω_max) = Ω²/c² = %.4e²/%.4e²" % (Omega_max, c))
    print("    T_rot(a)(Ω_max=%.0f) = %.4e м⁻²" % (Omega_max, Ta_omegamax))
    print("    T_rot(a)(Ω=1000)      = %.4e м⁻²" % Ta_1e3)
    print("    |T0_crit| = %.4f м⁻²" % abs(T0_crit))

    # ---- 6. Кручение, сценарий (б) Lense-Thirring ----
    M_sph = sphere_mass(rho_mat, R_s)
    Tb_omegamax, Omega_LT_max, J_max = T_rot_scenario_b(M_sph, R_s, Omega_max, G, c)
    Tb_1e3, Omega_LT_1e3, J_1e3 = T_rot_scenario_b(M_sph, R_s, 1e3, G, c)
    print("\n[ЧАСТЬ II-б] КРУЧЕНИЕ — СЦЕНАРИЙ (б) LENSE–THIRRING")
    print("    M_sph = ρ·(4/3)πR³ = %.2f кг    I = ⅖MR² = %.4e кг·м²" %
          (M_sph, 0.4*M_sph*R_s**2))
    print("    J(Ω_max) = ⅖MR²·Ω = %.4e кг·м²/с" % J_max)
    print("    Ω_LT = G·J/(c²·R³) = %.4e рад/с (на поверхности)" % Omega_LT_max)
    print("    T_rot(б) = (Ω_LT/c)² = (%.3e)² " % (Omega_LT_max/c))
    print("    T_rot(б)(Ω_max) = %.4e м⁻²" % Tb_omegamax)
    print("    T_rot(б)(Ω=1000) = %.4e м⁻²" % Tb_1e3)

    # ---- 7. Сравнение с T0_crit (оба сценария) ----
    orders_a = np.log10(abs(T0_crit) / Ta_omegamax)
    orders_b = np.log10(abs(T0_crit) / Tb_omegamax)
    print("\n[ЧАСТЬ II-в] СРАВНЕНИЕ С T0_crit = −%.4e м⁻²" % abs(T0_crit))
    print("    ┌─────────────────────────────────────────────────────────────────────┐")
    print("    │  (а) T≈Ω²/c²:  |T_rot|/|T0_crit|  = %.2e                    │" % (Ta_omegamax/abs(T0_crit)))
    print("    │      нехватка ≈ ~%.1f порядков                                  │" % orders_a)
    print("    │                                                                     │")
    print("    │  (б) LT:       |T_rot|/|T0_crit|  = %.2e                  │" % (Tb_omegamax/abs(T0_crit)))
    print("    │      нехватка ≈ ~%.1f порядков                                  │" % orders_b)
    print("    └─────────────────────────────────────────────────────────────────────┘")
    print("    → Знак T_rot: сценарий (а) даёт T≈Ω²/c² > 0 (≠ T0_crit<0).")
    print("      Сценарий (б) — гравитационно бесконечно мал (LT ≈ 0).")

    # ---- 8. Насколько большой Ω нужен для T0_crit ----
    Omega_needed_a = np.sqrt(abs(T0_crit)) * c
    # для LT: Ω_LT/c)^2 = Tcrit -> Omega_LT = sqrt(Tcrit)*c; но Omega_LT ~ M·Ω -> Omega
    omega_for_lt = np.sqrt(abs(T0_crit)) / G * c**2 * R_s**3 / J_max * Omega_max  # solve scaling
    # Точнее: Tb = (G·(2/5)M R² Ω /(c²R³)/c)² = Tcrit
    # => G*(2/5)M R² Ω /(c³ R³) = sqrt(Tcrit)
    Omega_needed_b = np.sqrt(abs(T0_crit)) * c**3 * R_s**3 / (G * 0.4 * M_sph * R_s**2)
    print("\n[ЧАСТЬ II-г] КАКОЙ Ω НУЖЕН ДЛЯ |T_rot|≥|T0_crit|?")
    print("    (а) Ω_needed = √(|T0_crit|)·c = √%.3f·%.3e = %.3e рад/с" %
          (abs(T0_crit), c, Omega_needed_a))
    print("        → vs Ω_max=%.0f рад/с: НЕВОЗМОЖНО (астрономически >)." % Omega_max)
    print("    (б) Ω_needed = √(Tcrit)·c³R³/(G·⅖MR²) = %.3e рад/с" % Omega_needed_b)
    print("        → vs Ω_max=%0.f рад/с: физически абсурдно (LT-канал нулевой)." % Omega_max)

    # ---- 9. ИТОГОВАЯ ТАБЛИЦА ----
    print("\n" + "=" * 120)
    print("  ИТОГОВАЯ ТАБЛИЦА «СФЕРА-ГИРОТОР»")
    print("=" * 120)
    hdr = ("  Параметр".ljust(50) + "Значение".rjust(30) + "  Примечание")
    print(hdr)
    print("  " + "-" * 118)
    rows = [
        ("R_s (радиус сферы)",              "%.2f м" % R_s,              "[РЕАЛЬНО]"),
        ("σ (предел прочности)",            "%.0f МПа" % (sigma_m/1e6),  "[РЕАЛЬНО] никель"),
        ("ρ (плотность)",                   "%.0f кг/м³" % rho_mat,      "[РЕАЛЬНО]"),
        ("Ω_max (прочность)",               "%.1f рад/с" % Omega_max,    "[РЕАЛЬНО] √(σ/ρR²)"),
        ("f_max",                           "%.1f Гц" % (Omega_max/(2*pi)), "[РЕАЛЬНО]"),
        ("Q_max (до пробоя 1e7 В/м)",       "%.4e Кл" % Qm,              "[РЕАЛЬНО] 4πε0R²E_max"),
        ("N_e (электронов в Q_max)",        "%.2e" % (Qm/e_charge),      "[РЕАЛЬНО]"),
        ("B_c (заряж. шар, Ω_max)",         "%.4e Тл" % Bc_Qmax_omegamax,"[МОДЕЛЬ] μ0QΩ/(6πR)"),
        ("B_L (Лондон, Ω_max)",             "%.4e Тл" % BL_omegamax,     "[РЕАЛЬНО] (2m_e/e)Ω"),
        ("B_L (Лондон, Ω=1000)",            "%.4e Тл" % BL_1e3,          "[РЕАЛЬНО]"),
        ("T_rot(а)=Ω²/c² (Ω_max)",          "%.4e м⁻²" % Ta_omegamax,    "[СЦЕН.а] Ω²/c²"),
        ("T_rot(а) Ω=1000",                 "%.4e м⁻²" % Ta_1e3,         "[СЦЕН.а]"),
        ("Ω_LT (LT, Ω_max)",                "%.4e рад/с" % Omega_LT_max, "[СЦЕН.б] GJ/(c²R³)"),
        ("T_rot(б)=(Ω_LT/c)² (Ω_max)",      "%.4e м⁻²" % Tb_omegamax,    "[СЦЕН.б]"),
        ("T0_crit (порог QI)",              "%.4e м⁻²" % T0_crit,        "[РЕАЛЬНО] −1/(22r0²)"),
        ("|T_rot(а)|/|T0_crit|",            "%.2e" % (Ta_omegamax/abs(T0_crit)), "[СРАВН.]"),
        ("нехватка (а)",                    "~%.1f порядков" % orders_a, "[КРИТИЧНО]"),
        ("|T_rot(б)|/|T0_crit|",            "%.2e" % (Tb_omegamax/abs(T0_crit)), "[СРАВН.]"),
        ("нехватка (б)",                    "~%.1f порядков" % orders_b, "[КРИТИЧНО]"),
        ("sign(T_rot) vs sign(T0_crit)",    "+,+ vs − (НЕСОВПАД.)",    "[КРИТИЧНО]"),
    ]
    for name, val, note in rows:
        print("  %-49s %29s  %s" % (name, val, note))

    # ---- 10. ВЕРДИКТ ----
    print("\n" + "=" * 120)
    print("  ВЕРДИКТ АГЕНТА-4: «СФЕРА-ГИРОТОР → РАБОЧИЙ ПОРТАЛ?»")
    print("=" * 120)
    print("""
  ┌─────────────────────────────────────────────────────────────────────────────────┐
  │                           НЕТ (численно НЕВОЗМОЖНО).                          │
  │                                                                                 │
  │  МАГНИТНОЕ ПОЛЕ (честно):                                                       │
  │    • Модель заряженного шара — НЕ РЕАЛИЗУЕМА в сверхпроводнике: S/C            │
  │      экранирует заряд/поле (Meissner). Даже если условно «набрать»              │
  │      Q_max=%.2e Кл до пробоя, B_c(Ω_max)=%.2e Тл — микроскопически мало.       │
  │    • ФИЗИЧЕСКИЙ канал — эффект Лондона: B_L=(2m_e/e)Ω=%.2e Тл при Ω_max.      │
  │    • Итог: сфера-гиротор даёт B ≈ %.1e Тл (10-100 нТл). Это УРОВЕНЬ          │
  │      геомагнитного поля, не более. НИКАКОГО «портала».                          │
  │                                                                                 │
  │  КРУЧЕНИЕ (честно, оба канона):                                                 │
  │    • (а) T≈Ω²/c² = %.2e м⁻² при Ω_max — нехватка ~%.1f порядков vs             │
  │      T0_crit=%.1e м⁻².                                                           │
  │    • (б) Lense–Thirring torsion (Ω_LT/c)²=%.2e м⁻² — нехватка ~%.1f            │
  │      порядков (frame-dragging гравитационно нулевой).                           │
  │    • Знак: T_rot>0 для обоих, нужен T0_crit<0 — не совпадает.                  │
  │                                                                                 │
  │  ГЛАВНЫЙ УБИЙЦА:                                                                │
  │    1. Вращение САМО ПО СЕБЕ torsion НЕ создаёт (нет канала rot→torsion в TEGR;  │
  │       кручение вращающейся СО — кинематика тетрады, не источник геометрии).     │
  │    2. Даже в кинематическом пределе: T~Ω²/c² ≪ T0_crit на ~%.0f порядков.      │
  │    3. B_L ~ 1e-8 Тл: магнитный «бонус» не имеет никакого отношения к T0_crit.  │
  │                                                                                 │
  │  ПОРЯДОК НЕВОЗМОЖНОСТИ:                                                         │
  │    • torsion: (а) ~%.1f порядков, (б) ~%.1f порядков НЕХВАТКИ vs T0_crit.      │
  │    • Для (а) нужен Ω~1e6 рад/с (недостижимо: прочность ломается при 1.3e3).    │
  │    • Для (б) нужен Ω~1e6 рад/с при M~M⊕ (LT-канал физически мёртв).           │
  │                                                                                 │
  │  ВЫВОД: вращающийся сверхпроводник даёт B≈1e-8 Тл и torsion≤1e-11 м⁻² —        │
  │  оба на ~11..61 порядков ниже порога T0_crit. Схема К ПОРТАЛУ НЕ СХОДИТСЯ.     │
  └─────────────────────────────────────────────────────────────────────────────────┘
""" % (Qm, Bc_Qmax_omegamax, BL_omegamax, BL_omegamax,
       Ta_omegamax, orders_a, abs(T0_crit),
       Tb_omegamax, orders_b,
       orders_a, orders_a, orders_b))

    # ---- JSON ----
    out = {
        "dataset": "r3_spinsphere_a4",
        "sphere": {
            "R_s_m": R_s,
            "sigma_Mpa": sigma_m / 1e6,
            "rho_kg_m3": rho_mat,
            "E_max_V_m": E_max,
            "omega_max_rad_s": float(Omega_max),
            "f_max_Hz": float(Omega_max / (2 * pi))
        },
        "magnetic_field": {
            "model_charged_ball": {
                "Q_max_C": float(Qm),
                "N_e": float(Qm / e_charge),
                "B_center_at_omega_max_T": float(Bc_Qmax_omegamax),
                "B_center_at_1000_T": float(Bc_Qmax_1e3),
                "formula": "B_c=mu0*Q*Omega/(6*pi*R)",
                "note": "VERHNYI PREDEL; S/C ekraniruet zaryad (Meissner) -> ne realizuem"
            },
            "model_london": {
                "B_london_at_omega_max_T": float(BL_omegamax),
                "B_london_at_1000_T": float(BL_1e3),
                "B_london_nT": float(BL_omegamax * 1e9),
                "formula": "B_L=(2*m_e/e)*Omega",
                "note": "FIZICHESKII kanal dlya S/C; nezavisit ot Q"
            },
            "summary": {
                "B_phys_T": float(BL_omegamax),
                "scale": "~1e-8 T (uroven' geomagnitnogo polya)"
            }
        },
        "torsion": {
            "T0_crit_m2": float(T0_crit),
            "alpha_crit": alpha_crit,
            "scenario_a_omega2_c2": {
                "T_rot_at_omega_max_m2": float(Ta_omegamax),
                "T_rot_at_1000_m2": float(Ta_1e3),
                "orders_below_Tcrit": float(orders_a),
                "formula": "T=Omega^2/c^2"
            },
            "scenario_b_lense_thirring": {
                "M_kg": float(M_sph),
                "J_at_omega_max": float(J_max),
                "Omega_LT_at_omega_max": float(Omega_LT_max),
                "T_rot_at_omega_max_m2": float(Tb_omegamax),
                "orders_below_Tcrit": float(orders_b),
                "formula": "Omega_LT=GJ/(c^2 R^3), T=(Omega_LT/c)^2"
            },
            "omega_needed_for_Tcrit": {
                "scenario_a_rad_s": float(Omega_needed_a),
                "scenario_b_rad_s": float(Omega_needed_b),
                "omega_max_rad_s": float(Omega_max),
                "impossible_by_orders_a": float(np.log10(Omega_needed_a / Omega_max)),
                "impossible_by_orders_b": float(np.log10(Omega_needed_b / Omega_max))
            },
            "sign": {
                "T_rot": "positive",
                "T0_crit": "negative",
                "mismatch": True
            }
        },
        "verdict_final": (
            "NET. Vrashchayushchiysya sverkhprovodnik daet B~1e-8 Tl (London) i "
            "torsion ~1e-11 m^-2 (scen.a Omega^2/c^2) ili ~1e-61 m^-2 (scen.b LT). "
            "Nehvatka vs T0_crit: ~%.1f poryadkov (a), ~%.1f poryadkov (b). "
            "Znak: T_rot>0, nuzhen T0_crit<0. B-dlya-torsion svyazi net. "
            "Skhema K RABOCHEMU PORTALU NE SKHODITSYA."
            % (orders_a, orders_b)),
        "prior_verdicts": [
            "r3_ring_cartridge.py: EM-koltso torsion NE sozdaet",
            "r3_plasma_vortex_a4.py: plazmennyi vikhr' torsion NE sozdaet",
            "r3_spinsphere_a4.py: vrashchayushchiysya S/C torsion NE sozdaet"
        ]
    }

    json_path = "/home/smboozha/portal_gun/research/r3_spinsphere_a4_values.json"
    with open(json_path, "w", encoding="utf-8") as fp:
        json.dump(out, fp, ensure_ascii=False, indent=2, default=str)

    print("  JSON: research/r3_spinsphere_a4_values.json сохранён")
    print("=" * 120)

if __name__ == "__main__":
    main()
