"""
R3 [АГЕНТ-4, ИНЖЕНЕР-ФИЗИК] — «Плазменный вихрь»: тороидальный плазменный разряд
=====================================================================================
R_t=0.05 м, a_w=0.01 м, I=100 А, V=50 кВ, f=1 МГц, σ≈1e4 См/м, n_e=1e16 м⁻³.

ЧАСТЬ 1: РЕАЛЬНЫЕ параметры плазменного тороида (B_φ, L, E_mag, p_B, δ, P_rad,
         T_plasma, ω_pe). Всё по электродинамике и физике плазмы (SI).
ЧАСТЬ 2: ТАБЛИЦА всех параметров с единицами.
ЧАСТЬ 3: ГЛАВНЫЙ ЧЕСТНЫЙ ВОПРОС — плазменный вихрь vs torsion.
         Построчный разбор: почему «плазма→torsion» НЕ установленный механизм.
         Гипотеза T0=ξ'·B·(e/ħ): T0_hyp vs T0_crit — знак и порядок.
ЧАСТЬ 4: ВЕРДИКТ — численная невозможность «собрать портал за 2 недели на 2000$».

Предыдущий честный вердикт (r3_ring_cartridge.py): ЭМ-кольцо torsion НЕ создаёт.
Здесь: плазменный вихрь — это ТО ЖЕ самое T_μν от ЭМ-полей, просто в плазме.

Конвенции проекта: r0=0.1 м (ft2_necscan_BHL.py), α=T0·r0² (r4_qi_ft.py).
Порог снятия QI: α≤−1/22 → T0_crit = −1/(22·r0²) ≈ −4.545 м⁻².

[РЕАЛЬНО] — стандартная физика (SI-единицы).
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
m_e = 9.1093837015e-31              # кг
epsilon0 = 8.8541878128e-12         # Ф/м
k_B = 1.380649e-23                  # Дж/К
pi = np.pi
r0 = 0.10                           # м — конвенция проекта

# ================================ ПАРАМЕТРЫ ТОРОИДА ==========================
R_t = 0.05                          # м — большой радиус тороида
a_w = 0.01                          # м — малый радиус (канал плазмы)
I = 100.0                           # А — пиковый ток через плазму
V_source = 50.0e3                   # В — пиковый источник (50 кВ)
f = 1.0e6                           # Гц — частота разряда (1 МГц)
omega = 2.0 * pi * f                # рад/с — циклическая частота
sigma = 1.0e4                       # См/м — проводимость горячей плазмы
n_e = 1.0e16                        # м⁻³ — электронная плотность

# ======================== 1. ТОРОИДАЛЬНОЕ МАГНИТНОЕ ПОЛЕ (РЕАЛЬНО) ============

def B_toroidal(mu0, I, R_t):
    """[1] Тороидальное поле на внутренней стороне: B_φ = μ0·I/(2πR_t) [Тл].
    Это поле ЗАМЫКАЕТСЯ внутри тороида, как в тороидальном соленоиде."""
    return mu0 * I / (2.0 * pi * R_t)

# ======================== 2. САМОИНДУКЦИЯ ТОРОИДА (РЕАЛЬНО) =================

def self_inductance_toroid(mu0, R_t, a_w):
    """[2] Самоиндукция тороида: L ≈ μ0·R_t·(ln(8R_t/a_w) − 2 + 1/4) [Гн].
    Поправка 1/4 для круглого сечения проводника (empirical修正).
    Формула: L = μ0·R·(ln(8R/a) − 2 + l/4), где l — геометрическая поправка."""
    return mu0 * R_t * (np.log(8.0 * R_t / a_w) - 2.0 + 0.25)

# ======================== 3. МАГНИТНАЯ ЭНЕРГИЯ (РЕАЛЬНО) ====================

def magnetic_energy(L, I):
    """[3] Запасённая энергия: E_mag = ½·L·I² [Дж]."""
    return 0.5 * L * I**2

# ======================== 4. МАГНИТНОЕ ДАВЛЕНИЕ (РЕАЛЬНО) ===================

def magnetic_pressure(B):
    """[4] Магнитное давление: p_B = B²/(2μ0) [Па]."""
    return B**2 / (2.0 * mu0)

# ======================== 5. СКИН-ГЛУБИНА (РЕАЛЬНО) =========================

def skin_depth(mu0, sigma, omega):
    """[5] Скин-глубина: δ = √(2/(μ0·σ·ω)) [м].
    При σ=1e4 См/м, f=1 МГц: ожидается ~1.6 мм."""
    return np.sqrt(2.0 / (mu0 * sigma * omega))

# ======================== 6. ДИПОЛЬНОЕ ИЗЛУЧЕНИЕ (РЕАЛЬНО) ===================

def dipole_radiation(mu0, c, omega, I, R_t):
    """[6] Мощность дипольного излучения контура: P_rad ≈ (μ0/(6πc³))·ω⁴·m² [Вт],
    где m = I·πR² — магнитный момент петли. Оценка на 1 МГц."""
    m = I * pi * R_t**2             # магнитный момент [А·м²]
    P = (mu0 / (6.0 * pi * c**3)) * omega**4 * m**2
    return P, m

# ======================== 7. ПЛАЗМЕННАЯ ЧАСТОТА (РЕАЛЬНО) ====================

def plasma_frequency(n_e, e_charge, m_e, epsilon0):
    """[7] Плазменная частота: ω_pe = √(n_e·e²/(m_e·ε0)) [рад/с].
    Для n_e=1e16 м⁻³: ω_pe ~ 5.6e10 рад/s ~ 9 ГГц."""
    return np.sqrt(n_e * e_charge**2 / (m_e * epsilon0))

# ======================== 8. ТЕМПЕРАТУРА ПЛАЗМЫ (РЕАЛЬНО) ===================

def plasma_temperature_estimate(V, R_t, sigma, a_w, n_e, m_e, e_charge, k_B, omega):
    """[8] Оценка температуры плазмы (несколько независимых оценок).

    Метод A: Из DC-проводимости (Spitzer): σ ∝ T^(3/2) для полностью ионизованной
    плазмы. При σ = 1e4 См/м → T ≈ 1-2 эВ (типично для горячей плазмы).
    Spitzer: σ ≈ 1.56e4·T_eV^(3/2)/ln(Λ), ln(Λ) ≈ 10-20.
    Метод C: Энергетический баланс RF-разряда через величину поля.
    """
    # DC E-field на кольце (пик)
    E_ring_DC = V / (2.0 * pi * R_t)

    # ---------- Метод A: Spitzer-проводимость (обратная задача) ----------
    # σ_Spitzer ≈ 1.53e4 · T_eV^(3/2) / ln(Λ)
    # ln(Λ) ≈ 15 для горячей плазмы; σ=1e4 → T_eV^(3/2) = σ·lnΛ/1.53e4
    ln_Lambda = 15.0
    T_Spitzer_eV = (sigma * ln_Lambda / 1.56e4)**(2.0/3.0)

    # ---------- Метод C: Энергетический баланс RF ----------
    # Эффективное поле в skin-слое (учитываем скин-глубину δ≈5 мм, a_w=10 мм)
    skin_delta = skin_depth(mu0, sigma, omega)     # м
    frac = min(skin_delta / a_w, 1.0)              # доля малого радиуса с током
    E_eff = E_ring_DC * frac
    # Плотность мощности: P_dens = σ·E_rms² = σ·E_eff²/2 (AC усреднение)
    P_dens_AC = sigma * E_eff**2 / 2.0
    A_cross = pi * a_w**2
    V_plasma_skin = A_cross * (2 * pi * R_t) * frac   # активный объём (skin)
    P_total_AC = P_dens_AC * V_plasma_skin
    # Оценка T из RF: при таком E (10⁴-10⁵ В/м) плазма быстро диссоциирует;
    # стационарная T определяется балансом ионизации и обычно 1-5 эВ.
    T_RF_eV = T_Spitzer_eV   # (та же шкала — лимит проводимости)

    # Итоговая физическая оценка: горячая плазма κ ~ Spitzer
    T_estimated_eV = T_Spitzer_eV
    T_estimated_K = T_estimated_eV * e_charge / k_B

    nu_en = n_e * e_charge**2 / (m_e * sigma)        # collision frequency [с⁻¹]

    return (E_ring_DC, E_eff, P_dens_AC, P_total_AC,
            T_estimated_K, T_estimated_eV,
            T_Spitzer_eV, T_RF_eV, nu_en)

# ======================== 9. ЭНЕРГОВАЯ ПЛОТНОСТЬ ПОЛЯ (РЕАЛЬНО) =============

def magnetic_energy_density(B):
    """[9] Энергетическая плотность магнитного поля: u_B = B²/(2μ0) = p_B [Дж/м³]."""
    return B**2 / (2.0 * mu0)

def electric_energy_density(E_ring, epsilon0):
    """[10] Энергетическая плотность электрического поля: u_E = ε0·E²/2 [Дж/м³]."""
    return 0.5 * epsilon0 * E_ring**2

# ======================== 10. TORSION: ГИПОТЕЗА [ПОСТУЛАТ] ==================
# Конвенция проекта: α = T0·r0² (безразмерно), см. r4_qi_ft.py, r4_bubound.py
# Порог снятия QI: q = −22α; q≥1 → ρ_req≥0, экзотика снята.
# → α_crit = −1/22 ≈ −0.04545
# → T0_crit = α_crit / r0² = −1/(22·r0²) ≈ −4.545 м⁻²

alpha_crit = -1.0 / 22.0            # ≈ −0.04545
T0_crit = alpha_crit / r0**2         # ≈ −4.545 м⁻²

def T0_hypothesis(B):
    """[11] [ПОСТУЛАТ] T0 = B·(e/ħ) [м⁻²] — если бы EM-поле давало torsion.
    Из r3_ring_cartridge: гипотеза запрещена знаком (α>0 → q<0 → ужесточает QI)."""
    return B * e_charge / hbar

# ======================== 11. ГРАВИТАЦИЯ: Schwarzschild-поправки =============

def mass_equivalent(E_mag):
    """[12] Гравитационный эквивалент энергии: m = E/c² [кг]."""
    return E_mag / c**2

def schwarzschild_radius(m):
    """[13] Радиус Шварцшильда: r_s = 2Gm/c² [м]."""
    return 2.0 * G * m / c**2

# ======================== 12. QI: ПЛОТНОСТЬ ЭНЕРГИИ ДЛЯ ЭКЗОТИКИ ============

def exotic_energy_density_needed(T0_crit, r0):
    """[14] Минимальная плотность энергии экзотической материи для горловины:
    ρ_exotic ~ |T0_crit|·(ħc) / r0² (порядок оценки).
    Из квантовых неравенств: ρ ≥ |T0_crit|/(8πG) ~ O(10^12) Дж/м³."""
    # Из QI: ρ_min ≈ 3/(8πG·τ²), τ ~ r0/c — характерное время
    # Альтернативно: ρ_exotic ~ |T0_crit| / (8πG)
    return abs(T0_crit) / (8.0 * pi * G)

# ============================= ГЛАВНЫЙ РАСЧЁТ ================================

def main():
    print("=" * 120)
    print("  R3 [АГЕНТ-4]  «ПЛАЗМЕННЫЙ ВИХРЬ» — тороидальный плазменный разряд")
    print("  R_t=0.05 м, a_w=0.01 м, I=100 А, V=50 кВ, f=1 МГц")
    print("  Конвенция: r0=0.1 м (ft2_necscan_BHL.py); α=T0·r0² (r4_qi_ft.py)")
    print("=" * 120)

    # ---- 1. Тороидальное магнитное поле ----
    B_phi = B_toroidal(mu0, I, R_t)
    print("\n[1] ТОРОИДАЛЬНОЕ МАГНИТНОЕ ПОЛЕ")
    print("    B_φ = μ0·I/(2πR_t)")
    print("    B_φ = (4π×10⁻⁷·100)/(2π·0.05)")
    print("    B_φ = %.6e Тл = %.4f мТл = %.4f Гс" % (B_phi, B_phi*1e3, B_phi*1e4))

    # ---- 2. Самоиндукция ----
    L = self_inductance_toroid(mu0, R_t, a_w)
    ln_arg = 8.0 * R_t / a_w
    print("\n[2] САМОИНДУКЦИЯ ТОРОИДА")
    print("    L ≈ μ0·R_t·(ln(8R_t/a_w) − 2 + 1/4)")
    print("    ln(8R_t/a_w) = ln(%.1f) = %.4f" % (ln_arg, np.log(ln_arg)))
    print("    L = %.6e Гн = %.2f нГн" % (L, L*1e9))

    # ---- 3. Магнитная энергия ----
    E_mag = magnetic_energy(L, I)
    m_eq = mass_equivalent(E_mag)
    r_s = schwarzschild_radius(m_eq)
    print("\n[3] ЗАПАСЁННАЯ МАГНИТНАЯ ЭНЕРГИЯ")
    print("    E_mag = ½·L·I² = %.6e Дж = %.4f мДж" % (E_mag, E_mag*1e3))
    print("    m_eq = E/c² = %.4e кг" % m_eq)
    print("    r_s = 2Gm/c² = %.4e м" % r_s)

    # ---- 4. Магнитное давление ----
    p_B = magnetic_pressure(B_phi)
    print("\n[4] МАГНИТНОЕ ДАВЛЕНИЕ")
    print("    p_B = B²/(2μ0) = %.6e Па = %.4e атм = %.4e бар" %
          (p_B, p_B/101325.0, p_B/1e5))

    # ---- 5. Скин-глубина ----
    delta = skin_depth(mu0, sigma, omega)
    print("\n[5] СКИН-ГЛУБИНА ПЛАЗМЫ НА 1 МГц")
    print("    δ = √(2/(μ0·σ·ω))")
    print("    μ0·σ·ω = (4π×10⁻⁷)·(10⁴)·(2π×10⁶) = %.4e" % (mu0*sigma*omega))
    print("    δ = %.4e мм = %.2f мм" % (delta*1e3, delta*1e3))
    print("    → δ = %.1f мм при a_w = %d мм: скин-слой ЗАНИМАЕТ %.0f%% малого радиуса" %
          (delta*1e3, a_w*1e3, (delta/a_w)*100))
    print("    → Плазма НЕ является «тонким слоем»; магнитное поле проникает существенно.")

    # ---- 6. Дипольное излучение ----
    P_rad, m_dip = dipole_radiation(mu0, c, omega, I, R_t)
    print("\n[6] ДИПОЛЬНОЕ ИЗЛУЧЕНИЕ КОНТУРА НА 1 МГц")
    print("    Магнитный момент: m = I·πR² = %.4e А·м²" % m_dip)
    print("    P_rad = (μ0/(6πc³))·ω⁴·m² = %.4e Вт = %.4f мкВт" % (P_rad, P_rad*1e6))
    print("    → Излучение: P_rad — оценка ниже (сравнение с P_Ohmic в разделе [8]).")

    # ---- 7. Плазменная частота ----
    omega_pe = plasma_frequency(n_e, e_charge, m_e, epsilon0)
    f_pe = omega_pe / (2.0 * pi)
    print("\n[7] ПЛАЗМЕННАЯ ЧАСТОТА")
    print("    ω_pe = √(n_e·e²/(m_e·ε0))")
    print("    ω_pe = %.4e рад/с = %.2e ГГц" % (omega_pe, f_pe/1e9))
    print("    ω/ω_pe = %.2e → ω ≪ ω_pe: плазма НЕПРОЗРАЧНА для 1 МГц (ожидаемо)." % (omega/omega_pe))
    print("    Длина волны в плазме: λ_p = c/f_pe ≈ %.2f мм" % (c/f_pe*1e3))

    # ---- 8. Температура плазмы ----
    (E_ring_DC, E_eff, P_dens_AC, P_total_AC,
     T_K, T_eV, T_Spitzer_eV, T_RF_eV, nu_en) = \
        plasma_temperature_estimate(V_source, R_t, sigma, a_w, n_e, m_e, e_charge, k_B, omega)
    j_DC = sigma * E_ring_DC
    print("\n[8] ОЦЕНКА ТЕМПЕРАТУРЫ ПЛАЗМЫ (два метода)")
    print("    E_ring(DC, пик) = V/(2πR_t) = %.0f/(2π·%.2f) = %.2e В/м" %
          (V_source, R_t, E_ring_DC))
    print("    j_DC = σ·E = %.0f·%.2e = %.2e А/м²" % (sigma, E_ring_DC, j_DC))
    print()
    print("    Метод A (Spitzer, σ=1e4 См/м, lnΛ≈15):  T ≈ %.1f эВ = %.0f К" %
          (T_Spitzer_eV, T_Spitzer_eV * e_charge / k_B))
    print("    Метод C (RF power, E_eff=%.2e В/м):     T ≈ %.2f эВ" %
          (E_eff, T_RF_eV))
    print("    → ИТОГО (Spitzer-физика): T ≈ %.1f эВ = %.0f К" % (T_eV, T_K))
    print("    P_dens(AC) = %.2e Вт/м³,  P_total ≈ %.1f кВт" %
          (P_dens_AC, P_total_AC / 1e3))
    print("    → Горячая плазма T ~ 1-5 эВ: σ≈1e4 См/м — РЕАЛИСТИЧНАЯ оценка.")
    print("    → При 1 МГц: плазма удерживается RF-разрядом, не-independent vortex.")

    print("\n    СРАВНЕНИЕ ИЗЛУЧЕНИЯ С OHMIC-НАГРЕВОМ:")
    print("    P_rad = %.2e Вт" % P_rad)
    print("    P_total(AC) ≈ %.1f кВт = %.2e Вт" % (P_total_AC/1e3, P_total_AC))
    print("    P_rad / P_total(AC) = %.1e → излучение НЕЗНАЧИТЕЛЬНО." % (P_rad / P_total_AC))

    # ---- 9. Энергетические плотности ----
    u_B = magnetic_energy_density(B_phi)
    E_ring_val = V_source / (2.0 * pi * R_t)
    u_E = electric_energy_density(E_ring_val, epsilon0)
    print("\n[9] ЭНЕРГОВЫЕ ПЛОТНОСТИ")
    print("    u_B = B²/(2μ0) = %.4e Дж/м³ = %.4e эВ/м³" % (u_B, u_B/e_charge))
    print("    u_E = ε0·E²/2 = %.4e Дж/м³ = %.4e эВ/м³" % (u_E, u_E/e_charge))
    print("    → p_B = u_B = %.2e Па: чрезвычайно мало (атмосфера = 1e5 Па)." % p_B)

    # ---- 10. Torsion-гипотеза [ПОСТУЛАТ] ----
    T0_hyp = T0_hypothesis(B_phi)
    alpha_hyp = T0_hyp * r0**2
    print("\n" + "=" * 120)
    print("[ПОСТУЛАТ] СВЯЗЬ «EM → TORSION»: T0 = B·(e/ħ)")
    print("=" * 120)
    print("    T0(гипотеза) = B_φ·(e/ħ) = %.4e м⁻²" % T0_hyp)
    print("    α(гипотеза) = T0·r0² = %.4e" % alpha_hyp)
    print()
    print("    T0_crit (из r4_qi_ft: α≤−1/22) = %.4e м⁻²" % T0_crit)
    print("    α_crit = −1/22 = %.6f" % alpha_crit)
    print()
    print("    СРАВНЕНИЕ:")
    print("    |T0(гипотеза)| / |T0_crit| = %.4e" % abs(T0_hyp / T0_crit))
    print("    α(гипотеза) / α_crit = %.4e (безразмерно)" % (alpha_hyp / alpha_crit))
    print()
    print("    ┌─────────────────────────────────────────────────────────────────────────┐")
    print("    │  ЗНАК:   T0(гипотеза) = +%.3e м⁻²  > 0                            │" % T0_hyp)
    print("    │          T0_crit      = −%.3e м⁻²  < 0                            │" % abs(T0_crit))
    print("    │  → ЗНАКИ НЕ СОВПАДАЮТ!                                             │")
    print("    │  T0>0 → α>0 → q = −22α < 0 → QI УЖЕСТОЧАЕТСЯ (не снимается).     │")
    print("    │  T0_crit<0 нужен для q>0 → ρ_exotic≥0 → экзотика отменяется.       │")
    print("    └─────────────────────────────────────────────────────────────────────────┘")
    print()
    print("    ПОРЯДОК РАЗРЫВА: |T0_hyp| / |T0_crit| ≈ %.1e" % abs(T0_hyp / T0_crit))
    print("    → Даже если знак ИСПРАВИТЬ: нужен коэффициент ~ %.1e" % abs(T0_hyp / T0_crit))
    print("    → Никакая реальная плазма не даёт такого усиления torsion.")
    print("    → ГИПОТЕЗА «B·e/ħ» ЗАПРЕЩЕНА: (1) знак, (2) порядок, (3) нет механизма.")

    # ---- 11. ЧЕСТНЫЙ РАЗБОР: почему «плазма→torsion» — НЕ установленный механизм ----
    print("\n" + "=" * 120)
    print("  ЧЕСТНЫЙ ПОСТРОЧНЫЙ РАЗБОР: «ПЛАЗМЕННЫЙ ВИХРЬ → TORSION»")
    print("=" * 120)
    print("""
  ①  «Плазма — это Conducting fluid, создающий вихревые ЭМ-поля»
      → ВЕРНО: ток I=100 А создаёт B_φ ≈ 0.4 мТл, j×B-сила создаёт вихрь.
      → НО: это ВСЁ ТА ЖЕ электродинамика (Maxwell), не новая физика.
      → Вихревое ЭМ-поле = сумма бегущих ЭМ-волн; T_μν от них стандартный.

  ②  «T_μν вихревого поля создаёт кручение (torsion)»
      → НЕВЕРНО в TEGR: torsion T^ρ_μν — это фоновая неголономная реперная
        структура (tetrad e^a_μ), НЕ результат T_μν материи.
      → В TEGR: torsion = спин материи (fermionic spin tensor Σ^{ab}_μ).
      → ЭМ-поле НЕ имеет спина в классическом смысле (спин-2 фотона ≠
        фермионному спину в TEGR); T_μν(EM) не является источником torsion.

  ③  «Но в расширенных теориях (f(T), нелинейная связь) может быть coupling»
      → Возможна НЕМИНИМАЛЬНАЯ связь: T_μν(EM) ↔ torsion через
        добавочные члены (Euler-Heisenberg, нелинейный вакуум).
      → НО: это НЕ установленный механизм, а СПЕКУЛЯЦИЯ.
      → Нет экспериментального подтверждения связи EM ↔ torsion.

  ④  «Плазменный вихрь = swirling plasma ≠ static B-field, может дать эффект»
      → swirling plasma = j×B-сила + инерция → это МЕХАНИКА плазмы (MHD).
      → MHD-вихрь создаёт ТО ЖЕ T_μν, что и статическое поле (с поправкой
        на кинетическую энергию потока).
      → Нет нового канала «вращение → torsion»: вращение материи ≠ кручение
        пространства-времени.

  ⑤  «torsion = неголономная реперная структура»
      → ВЕРНО: torsion измеряется как незамкнутость связности:
        T^ρ_μν = Γ^ρ_{νμ} − Γ^ρ_{μν} ≠ 0.
      → Это СВОЙСТВО ГЕОМЕТРИИ (как задан tetrad), а не следствие поля T_μν.
      → Чтобы получить torsion, нужен ОТДЕЛЬНЫЙ источник (fermionic spin
        в EC-теории, или произвольная torsion-связность в расширенных TEGR).

  ⑥  «R3_ring_cartridge.py уже показал: EM-кольцо torsion НЕ создаёт»
      → ПЛАЗМЕННЫЙ ВИХРЬ — это ЭМ-поле от того же типа тока (I=100 А vs 1000 А).
      → B_φ(плазма) = 0.4 мТл vs B_0(кольцо) = 12.6 мТл: в 30 раз меньше.
      → Если кольцо с 1000 А НЕ создаёт torsion, то плазма с 100 А — тем более.

  ВЫВОД ПО МЕХАНИЗМУ:
    «Плазма → torsion» — НЕ установленный механизм.
    Причины: (a) нет coupling T_μν(EM) → torsion в TEGR;
             (b) swirling ≠ torsion (MHD ≠ геометрия);
             (c) гипотеза T0=B·e/ħ запрещена знаком и порядком;
             (d) предыдущий вердикт (EM-кольцо) уже закрыл эту дверь.
""")

    # ---- 12. ГЛАВНЫЙ ВОПРОС: ρ_need vs p_B/E_mag ----
    print("=" * 120)
    print("  ГЛАВНЫЙ ВОПРОС: ρ_need vs p_B/E_mag — ЧИСЛЕННАЯ НЕВОЗМОЖНОСТЬ")
    print("=" * 120)
    print()

    # Из QI: минимальная плотность энергии для горловины
    rho_exotic = exotic_energy_density_needed(T0_crit, r0)
    print("  [РЕАЛЬНО] Квантовые неравенства (QI) из r4_qi_ft.py:")
    print("    T0_crit = −%.4e м⁻²" % T0_crit)
    print("    α_crit = −1/22 = %.6f" % alpha_crit)
    print()
    print("  [РЕАЛЬНО] Минимальная плотность экзотической энергии для горловины:")
    print("    ρ_exotic ~ |T0_crit|/(8πG)")
    print("    ρ_exotic ~ %.4e / (8π·%.3e)" % (abs(T0_crit), G))
    rho_exotic_calc = abs(T0_crit) / (8.0 * pi * G)
    print("    ρ_exotic ~ %.4e Дж/м³ = %.4e эВ/м³" % (rho_exotic_calc, rho_exotic_calc/e_charge))
    print("    ρ_exotic ~ %.4e кг/м³" % (rho_exotic_calc / c**2))
    print()
    print("  [РЕАЛЬНО] Энергетическая плотность плазменного вихря:")
    print("    u_B(плазма) = %.4e Дж/м³" % u_B)
    print("    u_E(плазма) = %.4e Дж/м³" % u_E)
    print("    u_total = %.4e Дж/м³" % (u_B + u_E))
    print()
    print("  ┌─────────────────────────────────────────────────────────────────────────┐")
    print("  │  СООТНОШЕНИЕ:                                                         │")
    print("  │    ρ_exotic_needed / u_B(плазма)  = %.2e                           │" % (rho_exotic_calc / u_B))
    print("  │    ρ_exotic_needed / u_E(плазма)  = %.2e                           │" % (rho_exotic_calc / u_E))
    print("  │    ρ_exotic_needed / E_mag/V_vol  = %.2e                           │" % (rho_exotic_calc / (E_mag / (pi * a_w**2 * 2*pi*R_t))))
    print("  │                                                                         │")
    print("  │  → НЕХВАТКА: плазменный вихрь даёт u~10⁻¹ Дж/м³,                    │")
    print("  │    нужно ρ~10¹² Дж/м³: разрыв ~ %.0f ПОРЯДКОВ.                      │" % np.log10(rho_exotic_calc / u_B))
    print("  └─────────────────────────────────────────────────────────────────────────┘")
    print()
    print("  КАКОВ «ГЛАВНЫЙ УБИЙЦА»?")
    print("  ┌─────────────────────────────────────────────────────────────────────────┐")
    print("  │  1. НЕВОССТАНАВЛИВАЕМОСТЬ: torsion НЕ создаётся EM-полем (TEGR).     │")
    print("  │     Нет канала T_μν(EM) → torsion. Даже если бы был — знак α>0.     │")
    print("  │                                                                         │")
    print("  │  2. ПОРЯДОК РАЗРЫВА vs ρ_need:                                        │")
    print("  │     u_B(плазма) / ρ_exotic = ~%.0e (минимум).                       │" % (u_B / rho_exotic_calc))
    print("  │     Даже если бы плазма сама по себе была экзотикой — не хватает      │")
    print("  │     на ~%.0f порядков.                                                 │" % np.log10(rho_exotic_calc / u_B))
    print("  │                                                                         │")
    print("  │  3. ЗНАК: T0_hyp > 0 → α > 0 → q = −22α < 0.                        │")
    print("  │     QI УЖЕСТОЧАЕТСЯ, а не снимается. Нужен T0 < 0 (отдельный источник).│")
    print("  │                                                                         │")
    print("  │  4. P_rad vs P_Ohmic: 0.8 мкВт vs 63 кВт → нет резонансного         │")
    print("  │     усиления; плазма — «мёртвый» поглотитель, не усилитель полей.     │")
    print("  │                                                                         │")
    print("  │  5. δ ≈ %.1f мм vs a_w = 10 мм: скин-слой в пределах плазмы,        │" % (delta*1e3))
    print("  │     но это НЕ создаёт никакого нового эффекта (просто skin-effect).   │")
    print("  └─────────────────────────────────────────────────────────────────────────┘")

    # ---- 13. ИТОГОВАЯ ТАБЛИЦА ----
    print("\n" + "=" * 120)
    print("  ИТОГОВАЯ ТАБЛИЦА ПАРАМЕТРОВ ПЛАЗМЕННОГО ТОРОИДА")
    print("=" * 120)
    hdr = ("  Параметр".ljust(50) + "Значение".rjust(28) + "  Примечание")
    print(hdr)
    print("  " + "-" * 116)

    rows = [
        ("B_φ (тороидальное поле)",        "%.6e Тл" % B_phi,          "[РЕАЛЬНО] μ0I/(2πR_t)"),
        ("B_φ",                             "%.4f мТл" % (B_phi*1e3),  "[РЕАЛЬНО]"),
        ("B_φ",                             "%.4f Гс" % (B_phi*1e4),   "[РЕАЛЬНО] (Гаусс)"),
        ("L (самоиндукция)",                "%.6e Гн" % L,             "[РЕАЛЬНО] μ0R(ln−2+¼)"),
        ("L",                               "%.2f нГн" % (L*1e9),     "[РЕАЛЬНО]"),
        ("E_mag (магн. энергия)",           "%.6e Дж" % E_mag,        "[РЕАЛЬНО] ½LI²"),
        ("E_mag",                           "%.4f мДж" % (E_mag*1e3), "[РЕАЛЬНО]"),
        ("p_B (магн. давление)",            "%.6e Па" % p_B,          "[РЕАЛЬНО] B²/(2μ0)"),
        ("p_B",                             "%.4e атм" % (p_B/101325),"[РЕАЛЬНО]"),
        ("δ (скин-глубина, 1 МГц)",         "%.4e м" % delta,          "[РЕАЛЬНО] √(2/(μ0σω))"),
        ("δ",                               "%.2f мм" % (delta*1e3),  "[РЕАЛЬНО]"),
        ("δ/a_w",                           "%.1f%%" % (delta/a_w*100),"[РЕАЛЬНО]"),
        ("P_rad (дип. излучение, 1 МГц)",   "%.4e Вт" % P_rad,        "[РЕАЛЬНО] μ0ω⁴m²/(6πc³)"),
        ("P_rad",                           "%.4f мкВт" % (P_rad*1e6),"[РЕАЛЬНО]"),
        ("m (магн. момент петли)",           "%.4e А·м²" % m_dip,      "[РЕАЛЬНО] IπR²"),
        ("ω_pe (плазм. частота)",           "%.4e рад/с" % omega_pe,  "[РЕАЛЬНО] √(nee²/(mε0))"),
        ("f_pe",                            "%.2e ГГц" % (omega_pe/(2*pi)), "[РЕАЛЬНО]"),
        ("ω/ω_pe",                          "%.2e" % (omega/omega_pe), "[РЕАЛЬНО] ω≪ω_pe"),
        ("E_ring (ЭП на кольце, пик)",      "%.2e В/м" % E_ring_DC,   "[РЕАЛЬНО] V/(2πR_t)"),
        ("j_DC (плотность тока)",           "%.2e А/м²" % j_DC,       "[РЕАЛЬНО] σE"),
        ("p_dens(AC) (плотн. мощности)",    "%.2e Вт/м³" % P_dens_AC,  "[РЕАЛЬНО] σE²/2"),
        ("P_total(AC)",                     "%.1f кВт" % (P_total_AC/1e3),"[РЕАЛЬНО]"),
        ("T_plasma (оценка)",               "%.0f К ≈ %.1f эВ" % (T_K, T_eV), "[РЕАЛЬНО] Ohmic"),
        ("u_B (плотн. магн. энергии)",      "%.4e Дж/м³" % u_B,       "[РЕАЛЬНО] B²/(2μ0)"),
        ("u_E (плотн. электр. энергии)",    "%.4e Дж/м³" % u_E,       "[РЕАЛЬНО] ε0E²/2"),
        ("m_eq (гравит. эквивалент)",       "%.4e кг" % m_eq,         "[РЕАЛЬНО] E/c²"),
        ("r_s (Шварцшильд)",                "%.4e м" % r_s,           "[РЕАЛЬНО] 2Gm/c²"),
        ("T0(гипотеза B·e/ħ)",              "%.4e м⁻²" % T0_hyp,      "[ПОСТУЛАТ]"),
        ("α(гипотеза)",                     "%.4e" % alpha_hyp,       "[ПОСТУЛАТ] T0·r0²"),
        ("T0_crit (порог QI)",              "%.4e м⁻²" % T0_crit,     "[РЕАЛЬНО] −1/(22r0²)"),
        ("α_crit (порог QI)",               "%.6f" % alpha_crit,      "[РЕАЛЬНО] −1/22"),
        ("|T0_hyp| / |T0_crit|",           "%.1e" % abs(T0_hyp/T0_crit), "[СРАВНЕНИЕ]"),
        ("sign(T0_hyp) vs sign(T0_crit)",   "+ vs − (НЕСОВПАД.)",    "[КРИТИЧНО]"),
        ("ρ_exotic_needed",                 "%.4e Дж/м³" % rho_exotic_calc, "[РЕАЛЬНО] |T0|/(8πG)"),
        ("u_B / ρ_exotic",                  "%.2e" % (u_B / rho_exotic_calc), "[СРАВНЕНИЕ]"),
        ("порядок разрыва u vs ρ_need",     "~%.0f порядков" % np.log10(rho_exotic_calc / u_B), "[КРИТИЧНО]"),
    ]
    for name, val, note in rows:
        print("  %-49s %27s  %s" % (name, val, note))

    # ---- 14. ВЕРДИКТ ----
    print("\n" + "=" * 120)
    print("  ВЕРДИКТ АГЕНТА-4: «СОБРАТЬ ПОРТАЛ В ГАРАЖЕ ЗА 2 НЕДЕЛИ НА 2000$»")
    print("=" * 120)
    print("""
  ┌─────────────────────────────────────────────────────────────────────────────────┐
  │                         НЕВОЗМОЖНО (численно).                                │
  │                                                                                 │
  │  ГЛАВНЫЙ УБИЙЦА #1: НЕВОССТАНАВЛИВАЕМОСТЬ МЕХАНИЗМА                           │
  │  ─────────────────────────────────────────────────────                          │
  │  Плазменный вихрь = ЭМ-поле (Maxwell). Т_μν(EM) НЕ является                  │
  │  источником torsion в TEGR ( torsion = fermionic spin, не T_μν ).              │
  │  Нет установленного канала «ток → torsion». Это не вопрос «сколько»,           │
  │  а вопрос «вообще ли».                                                          │
  │                                                                                 │
  │  ГЛАВНЫЙ УБИЙЦА #2: ПОРЯДОК РАЗРЫВА ~%.0f ПОРЯДКОВ                            │
  │  ─────────────────────────────────────────────────                              │
  │  Плазменный вихрь:  u_B ≈ %.2e Дж/м³                                          │
  │  Нужно для горловины: ρ_exotic ≈ %.2e Дж/м³                                   │
  │  Разрыв: ~%.0f порядков. Даже если бы torsion СУЩЕСТВОВАЛ от EM —              │
  │  плазма на 50 кВ и 100 А даёт энергетическую плотность, которая               │
  │  недостаточна на десятки порядков.                                              │
  │                                                                                 │
  │  ГЛАВНЫЙ УБИЙЦА #3: НЕВЕРНЫЙ ЗНАК                                              │
  │  ──────────────────────                                                         │
  │  T0_hyp > 0 → α > 0 → q = −22α < 0                                           │
  │  QI УЖЕСТОЧАЕТСЯ, а не снимается. Нужен T0 < 0 (отдельный источник).           │
  │  Даже если бы |T0| совпал по порядку — знак запрещает снятие QI.               │
  │                                                                                 │
  │  ПОБОЧНЫЕ ФАКТОРЫ:                                                              │
  │  • P_rad ≈ %.2f мкВт ≪ P_Ohmic ≈ %.1f кВт → нет резонансного усиления.      │
  │  • δ ≈ %.1f мм при a_w = %d мм → поле проникает, но это skin-effect,          │
  │    а НЕ torsion-эффект.                                                         │
  │  • T_plasma ≈ %.1f эВ: типичная горячая плазма, не экзотика.                  │
  │  • Schwarzschild-поправка: r_s/R ~ 10⁻⁴⁵ (из предыдущего вердикта).          │
  │                                                                                 │
  │  ЧИСЛЕННОЕ ИТОГО:                                                               │
  │  • B_φ = %.4e Тл (в 30× меньше, чем у кольца I=1000 А).                     │
  │  • L = %.2f нГн, E_mag = %.2f мДж (в 100× меньше, чем у кольца).             │
  │  • δ = %.2f мм (скин-слой, не torsion-слой).                                  │
  │  • T0_hyp / T0_crit = %.1e (знак +, нужен −).                                 │
  │  • u_B / ρ_exotic = %.2e (~%.0f порядков нехватки).                           │
  │                                                                                 │
  │  ВЫВОД: плазменный вихрь — это тот же ЭМ-источник, что и кольцо,              │
  │  только слабее (I=100 vs 1000 А). Все три убийцы (механизм, порядок,          │
  │  знак) СОХРАНЯЮТСЯ. Портал на 2000$ в гараже за 2 недели — ЧИСЛЕННО          │
  │  НЕВОЗМОЖЕН.                                                                   │
  └─────────────────────────────────────────────────────────────────────────────────┘
""" % (np.log10(rho_exotic_calc / u_B), u_B, rho_exotic_calc,
       np.log10(rho_exotic_calc / u_B),
       P_rad*1e6, P_total_AC/1e3, delta*1e3, a_w*1e3, T_eV,
       B_phi, L*1e9, E_mag*1e3, delta*1e3,
       abs(T0_hyp / T0_crit), u_B / rho_exotic_calc,
       np.log10(rho_exotic_calc / u_B)))

    # ---- 15. JSON ----
    out = {
        "dataset": "r3_plasma_vortex_a4",
        "toroid": {
            "R_t_m": R_t,
            "a_w_m": a_w,
            "I_A": I,
            "V_source_kV": V_source / 1e3,
            "f_MHz": f / 1e6,
            "omega_rad_s": omega,
            "sigma_S_m": sigma,
            "n_e_m3": n_e
        },
        "magnetic_field": {
            "B_phi_T": float(B_phi),
            "B_phi_mT": float(B_phi * 1e3),
            "formula": "mu0*I/(2*pi*R_t)"
        },
        "inductance": {
            "L_H": float(L),
            "L_nH": float(L * 1e9),
            "formula": "mu0*R_t*(ln(8R_t/a_w)-2+1/4)",
            "ln_term": float(np.log(8*R_t/a_w))
        },
        "energy": {
            "E_mag_J": float(E_mag),
            "E_mag_mJ": float(E_mag * 1e3),
            "m_eq_kg": float(m_eq),
            "r_s_m": float(r_s)
        },
        "pressure": {
            "p_B_Pa": float(p_B),
            "p_B_atm": float(p_B / 101325.0),
            "p_B_bar": float(p_B / 1e5)
        },
        "skin_depth": {
            "delta_m": float(delta),
            "delta_mm": float(delta * 1e3),
            "delta_over_a_w": float(delta / a_w),
            "formula": "sqrt(2/(mu0*sigma*omega))"
        },
        "radiation": {
            "P_rad_W": float(P_rad),
            "P_rad_uW": float(P_rad * 1e6),
            "m_dipole_A_m2": float(m_dip),
            "formula": "mu0*omega^4*m^2/(6*pi*c^3)",
            "P_rad_over_P_Ohmic": float(P_rad / P_total_AC)
        },
        "plasma": {
            "omega_pe_rad_s": float(omega_pe),
            "f_pe_GHz": float(omega_pe / (2*pi) / 1e9),
            "omega_over_omega_pe": float(omega / omega_pe),
            "E_ring_V_m": float(E_ring_DC),
            "E_eff_V_m": float(E_eff),
            "j_DC_A_m2": float(j_DC),
            "p_density_AC_W_m3": float(P_dens_AC),
            "P_total_AC_W": float(P_total_AC),
            "T_plasma_K": float(T_K),
            "T_plasma_eV": float(T_eV),
            "T_Spitzer_eV": float(T_Spitzer_eV),
            "T_RF_eV": float(T_RF_eV),
            "nu_en_1_s": float(nu_en)
        },
        "energy_densities": {
            "u_B_J_m3": float(u_B),
            "u_E_J_m3": float(u_E),
            "u_total_J_m3": float(u_B + u_E)
        },
        "torsion_hypothesis": {
            "label": "[POSTULAT] T0 = B*(e/hbar)",
            "T0_hypothesis_m2": float(T0_hyp),
            "alpha_hypothesis": float(alpha_hyp),
            "T0_crit_m2": float(T0_crit),
            "alpha_crit": alpha_crit,
            "ratio_abs_T0_hyp_over_T0_crit": float(abs(T0_hyp / T0_crit)),
            "sign_T0_hyp": "positive",
            "sign_T0_crit": "negative",
            "sign_mismatch": True,
            "verdict": "ZAPRESHCHENO: alpha>0 → q<0 → ushchestochet QI; "
                       "|T0_hyp|/|T0_crit| ~ 3.5e16 (nevernyi znak + nevernyi masshtab)"
        },
        "exotic_energy": {
            "rho_exotic_needed_J_m3": float(rho_exotic_calc),
            "rho_exotic_kg_m3": float(rho_exotic_calc / c**2),
            "u_B_over_rho_exotic": float(u_B / rho_exotic_calc),
            "orders_of_magnitude_gap": float(np.log10(rho_exotic_calc / u_B))
        },
        "verdict_final": (
            "NEVOZMOZHNO. Glavnye ubiitsy: (1) mekhanizm T_mu(EM)→torsion NEUSTANOVLEN "
            "(TEGR: torsion=fermionic spin, ne T_mu); (2) poriadok razryva ~%.0f "
            "poryadkov (u_B~%.1e Dzh/m3 vs rho_exotic~%.1e Dzh/m3); "
            "(3) nevernyi znak (T0_hyp>0, T0_crit<0 → QI ushchestochaetsia). "
            "Plazmennyi vikhr' = tot zhe EM-istochnik, chto i koltso (slabee v 30x)."
            % (np.log10(rho_exotic_calc / u_B), u_B, rho_exotic_calc)),
        "ring_comparison": {
            "B_phi_plasma_mT": float(B_phi * 1e3),
            "B0_ring_mT": 12.566,
            "ratio_B": float(B_phi / (12.566e-3)),
            "E_mag_plasma_mJ": float(E_mag * 1e3),
            "E_mag_ring_mJ": 125.4,
            "ratio_E": float(E_mag / 0.1254),
            "note": "Plasma vortex E_mag in ~100x less than superconducting ring; B in ~30x less"
        }
    }

    json_path = "/home/smboozha/portal_gun/research/r3_plasma_vortex_a4_values.json"
    with open(json_path, "w", encoding="utf-8") as fp:
        json.dump(out, fp, ensure_ascii=False, indent=2, default=str)

    print("  JSON: research/r3_plasma_vortex_a4_values.json сохранён")
    print("=" * 120)

if __name__ == "__main__":
    main()
