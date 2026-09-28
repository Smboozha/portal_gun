"""
R3 [АГЕНТ-2 ФИЗИКА ОГРАНИЧЕНИЙ] — QI Форд-Роман × МОДУЛИРОВАННЫЙ РЕЖИМ
==========================================================================
Сценарий: отрицательная плотность включается импульсами длиной τ_on с
дuty-cycle f = τ_on / T_cycle.

Ключевой факт: временная модуляция НЕ обходит FR-QI для мгновенной
локальной меры плотности. Bound применяется к трубке, покрывающей
случайный момент включения.

Два честных случая сэмплинга:
(а) τ₀ = τ_on  — трубка перекрывает отдельный импульс
(б) τ₀ = T_cycle — трубка усредняет полный цикл

Источник FR-QI: Ford & Roman, PRD 53 (1996) 5496
"""
import json
import numpy as np

# ========================= физические константы (SI) ========================
hbar = 1.054571817e-34
c    = 2.99792458e8
G    = 6.67430e-11
pi   = np.pi

R0   = 0.10
C4_G = c**4 / G

# ========================= FR-предел ========================================
def eps_FR_max(tau_s):
    """|ε|_max = 3ħ/(32π²c³τ⁴) [Дж/м³]"""
    return 3.0 * hbar / (32.0 * pi**2 * c**3 * tau_s**4)

# ========================= ρ_need из f(T) ===================================
def rho_need_SI(alpha, gamma, r0=R0):
    """ρ_th = −(1+22α+24γ)/2  [4πr0²-единицы]
       ρ_SI = ρ_dim × c⁴/(8πG·r0²)  [Дж/м³]"""
    rho_dim = -(1.0 + 22.0*alpha + 24.0*gamma) / 2.0
    rho_si  = rho_dim * C4_G / (8.0 * pi * r0**2)
    return rho_dim, rho_si

# ========================= main ==============================================
def main():
    alpha = -0.1
    gamma = -0.01

    rho_dim, rho_SI = rho_need_SI(alpha, gamma)

    # Параметры модуляции
    tau_on    = 1e-6       # 1 мкс — длительность импульса
    T_cycle   = 1e-3       # 1 мс  — полный цикл
    f_duty    = tau_on / T_cycle  # = 1/1000

    print("=" * 100)
    print("  R3 [АГЕНТ-2]  FR-QI × МОДУЛИРОВАННЫЙ РЕЖИМ: duty-cycle")
    print("=" * 100)

    # ---- Параметры модуляции ----
    print(f"\n  ПАРАМЕТРЫ МОДУЛЯЦИИ:")
    print(f"    τ_on     = {tau_on:.1e} с  (1 мкс)")
    print(f"    T_cycle  = {T_cycle:.1e} с  (1 мс)")
    print(f"    f = τ_on/T_cycle = {f_duty:.4f}  (1/1000)")
    print(f"    Частота повторения = {1/T_cycle:.0f} Гц = {1/T_cycle/1e3:.1f} кГц")

    # ---- ρ_need ----
    print(f"\n  ρ_need ИЗ f(T) (α={alpha}, γ={gamma}, r0={R0} м)")
    print(f"    ρ_dim = {rho_dim:.4f}  [4πr0²-единицы]")
    print(f"    |ρ_need| = {abs(rho_SI):.6e} Дж/м³")

    # ---- Случай (а): τ₀ = τ_on ----
    print(f"\n{'─'*100}")
    print(f"  СЛУЧАЙ (а): τ₀ = τ_on = {tau_on:.1e} с  (трубка перекрывает отдельный импульс)")
    print(f"{'─'*100}")

    eps_max_a = eps_FR_max(tau_on)
    ratio_a   = abs(rho_SI) / eps_max_a
    log_a     = np.log10(ratio_a)
    viol_a    = "ДА" if ratio_a > 1 else "НЕТ"

    print(f"    |ε|_max(τ_on) = 3ħ/(32π²c³τ_on⁴)")
    print(f"                 = {eps_max_a:.6e} Дж/м³")
    print(f"    │ρ_need│/│ε│_max = {ratio_a:.4e}")
    print(f"    log₁₀(нарушение) = {log_a:.1f}")
    print(f"    Нарушение QI: {viol_a}  (~{abs(log_a):.0f} порядков)")

    # ---- Случай (б): τ₀ = T_cycle ----
    print(f"\n{'─'*100}")
    print(f"  СЛУЧАЙ (б): τ₀ = T_cycle = {T_cycle:.1e} с  (трубка усредняет полный цикл)")
    print(f"{'─'*100}")

    eps_max_b = eps_FR_max(T_cycle)
    rho_avg   = f_duty * abs(rho_SI)       # средняя плотность за цикл
    ratio_b   = rho_avg / eps_max_b
    log_b     = np.log10(ratio_b) if ratio_b > 0 else -np.inf
    viol_b    = "ДА" if ratio_b > 1 else "НЕТ"

    print(f"    |ε|_max(T_cycle) = 3ħ/(32π²c³T_cycle⁴)")
    print(f"                    = {eps_max_b:.6e} Дж/м³")
    print(f"    ⟨ρ⟩ = f × |ρ_need| = {f_duty:.6f} × {abs(rho_SI):.6e}")
    print(f"        = {rho_avg:.6e} Дж/м³")
    print(f"    ⟨ρ⟩/│ε│_max(T_cycle) = {ratio_b:.4e}")
    print(f"    log₁₀(нарушение)     = {log_b:.1f}")
    print(f"    Нарушение QI: {viol_b}  (~{abs(log_b):.0f} порядков)")

    # ---- Ключевое сравнение ----
    print(f"\n{'═'*100}")
    print(f"  КЛЮЧЕВОЕ СРАВНЕНИЕ: СНИМАЕТ ЛИ DUTY-CYCLE НАРУШЕНИЕ?")
    print(f"{'═'*100}")
    print(f"    Случай (а)  τ₀=1 мкс:  ratio = {ratio_a:.4e}  (log₁₀ = {log_a:.1f})")
    print(f"    Случай (б)  τ₀=1 мс:   ratio = {ratio_b:.4e}  (log₁₀ = {log_b:.1f})")
    print()

    # Анализ масштабирования
    scale_bound = (tau_on / T_cycle)**4   # = f⁴ = 10⁻¹²
    scale_avg   = f_duty                  # = 10⁻³
    ratio_scale = scale_avg / scale_bound  # = f⁴/f = f³ = 10⁹

    print(f"    Масштабирование:")
    print(f"      Bound: |ε|_max ~ τ⁻⁴ → |ε|_max(1мс)/|ε|_max(1мкс) = (1μs/1ms)⁴ = {scale_bound:.4e}")
    print(f"      Средняя плотность: ⟨ρ⟩/|ρ| = f = {f_duty:.4e}")
    print(f"      Отношение (средняя/ bound) scales as: f / f⁴ = f⁻³ = {1/scale_bound * f_duty:.4e}")
    print(f"      → Duty-cycle УХУДШАЕТ нарушение в {ratio_scale:.0e} раз!")
    print()
    print(f"    Почему: bound падает как τ⁻⁴ (быстро),")
    print(f"      а средняя плотность падает только как f (медленно).")
    print(f"      f = 10⁻³: bound падает на 10¹², плотность — на 10³.")
    print(f"      Итого: ratio растёт в 10⁹ раз при усреднении по циклу.")

    # ---- Итоговая таблица ----
    print(f"\n{'═'*100}")
    print(f"  ИТОГОВАЯ ТАБЛИЦА: 4 строки")
    print(f"{'═'*100}")

    table_data = [
        {"τ₀": "1 мкс (τ_on)",   "|ε|_max [Дж/м³]": eps_max_a, "плотность": abs(rho_SI), "type": "пиковая", "ratio": ratio_a, "log10": log_a, "violation": viol_a},
        {"τ₀": "1 мс (T_cycle)", "|ε|_max [Дж/м³]": eps_max_b, "плотность": rho_avg,     "type": "средняя", "ratio": ratio_b, "log10": log_b, "violation": viol_b},
        {"τ₀": "1 мкс (τ_on)",   "|ε|_max [Дж/м³]": eps_max_a, "плотность": abs(rho_SI), "type": "пиковая", "ratio": ratio_a, "log10": log_a, "violation": viol_a},
        {"τ₀": "1 мс (T_cycle)", "|ε|_max [Дж/м³]": eps_max_b, "плотность": abs(rho_SI), "type": "пиковая", "ratio": abs(rho_SI)/eps_max_b, "log10": np.log10(abs(rho_SI)/eps_max_b), "violation": "ДА"},
    ]

    # Красивая таблица
    header = f"  {'τ₀':>14s}  {'|ε|_max [Дж/м³]':>18s}  {'плотность [Дж/м³]':>18s}  {'тип':>8s}  {'ratio':>14s}  {'log₁₀':>8s}  {'наруш.':>6s}"
    print(header)
    print("  " + "─" * 90)
    for row in table_data:
        print(f"  {row['τ₀']:>14s}  {row['|ε|_max [Дж/м³]']:18.6e}  {row['плотность']:18.6e}  {row['type']:>8s}  {row['ratio']:14.4e}  {row['log10']:8.1f}  {row['violation']:>6s}")

    print()
    print("  Примечания:")
    print("    1) Строка 1: случай (а) — пиковая плотность vs bound при τ₀=τ_on")
    print("    2) Строка 2: случай (б) — средняя плотность vs bound при τ₀=T_cycle")
    print("    3) Строка 3: та же что 1 (для наглядности)")
    print("    4) Строка 4: если «обмануть» и взять пик при большом τ₀ — всё равно нарушение")

    # ---- Вердикт ----
    print(f"\n{'═'*100}")
    print(f"  ВЕРДИКТ")
    print(f"{'═'*100}")
    print(f"""
    ┌──────────────────────────────────────────────────────────────────────────────────┐
    │  МОДУЛЯЦИЯ НЕ ОБХОДИТ FR-QI                                                    │
    │                                                                                  │
    │  Причина: FR-QI — это неравенство для мгновенной меры плотности энергии,        │
    │  интегрированной по световой трубке (worldline of an inertial observer).        │
    │  Bound применяется к трубке, покрывающей ЛЮБОЙ момент включения импульса.      │
    │                                                                                  │
    │  Случай (а): τ₀ = τ_on = 1 мкс                                                 │
    │    → |ρ|/|ε|_max ≈ {ratio_a:.1e}  (~{abs(log_a):.0f} порядков нарушения)                   │
    │    → Duty-cycle НЕ влияет: в момент включения импульс длительностью τ_on        │
    │      полностью покрыт трубкой и bound такой же как в статическом случае.        │
    │                                                                                  │
    │  Случай (б): τ₀ = T_cycle = 1 мс                                               │
    │    → ⟨ρ⟩/|ε|_max ≈ {ratio_b:.1e}  (~{abs(log_b):.0f} порядков нарушения)                   │
    │    → Усреднение НЕ помогает: bound падает как τ⁻⁴ = {scale_bound:.1e},             │
    │      а средняя плотность падает только как f = {f_duty:.1e}.                         │
    │      Итого: ratio растёт в ~{ratio_scale:.0e} раз по сравнению со случаем (а).          │
    │                                                                                  │
    │  Физический факт из литературы:                                                  │
    │    Временная модуляция НЕ обходит QI для мгновенной локальной меры.              │
    │    Если ρ(t₀) < 0 в момент t₀, tube с τ₀ = время жизни импульса                │
    │    фиксирует нарушение. Duty-cycle только добавляет «пустые» периоды,            │
    │    но не снижает пиковую плотность в момент включения.                           │
    │                                                                                  │
    │  ВЕРДИКТ: duty-cycle НЕ снимает QI.                                             │
    │    Случай (а): ~{abs(log_a):.0f} порядков нарушения (как в статическом режиме)      │
    │    Случай (б): ~{abs(log_b):.0f} порядков нарушения (ЕЩЁ ХУЖЕ из-за τ⁻⁴)          │
    │                                                                                  │
    │  Общий факт: модуляция НЕ обходит QI.                                           │
    └──────────────────────────────────────────────────────────────────────────────────┘
""")

    # ---- JSON ----
    out = {
        "file": "r3_qi_modulated_a2",
        "agent": "R3-АГЕНТ-2-ФИЗИКА-ОГРАНИЧЕНИЙ",
        "scenario": "modulated_duty_cycle",
        "FR_QI_formula": {
            "source": "Ford & Roman, PRD 53 (1996) 5496",
            "integral_form": "(tau0/pi) * int <T_mu nu u^mu u^nu>/(t^2+tau0^2) dt >= -3*hbar/(32*pi^2*c^3*tau0^4)",
            "local_limit": "|epsilon|_max = 3*hbar/(32*pi^2*c^3*tau^4) [J/m^3]",
        },
        "modulation_params": {
            "tau_on_s": tau_on,
            "T_cycle_s": T_cycle,
            "f_duty": f_duty,
            "frequency_Hz": 1.0 / T_cycle,
            "tau_on_label": "1 us",
            "T_cycle_label": "1 ms",
        },
        "rho_need_from_ft": {
            "model": "f(T) = T + alpha*T^2 + gamma*T^-1",
            "alpha": alpha,
            "gamma": gamma,
            "r0_m": R0,
            "rho_th_dim_4piR0sq": rho_dim,
            "rho_SI_J_per_m3": rho_SI,
            "abs_rho_SI": abs(rho_SI),
        },
        "case_a_tau0_eq_tau_on": {
            "tau0_s": tau_on,
            "eps_FR_max_J_per_m3": eps_max_a,
            "rho_used": abs(rho_SI),
            "rho_type": "peak (full amplitude during ON time)",
            "ratio": ratio_a,
            "log10_violation": log_a,
            "violated": ratio_a > 1,
            "orders_of_violation": abs(log_a),
        },
        "case_b_tau0_eq_T_cycle": {
            "tau0_s": T_cycle,
            "eps_FR_max_J_per_m3": eps_max_b,
            "rho_avg_J_per_m3": rho_avg,
            "rho_type": "average over full cycle",
            "ratio": ratio_b,
            "log10_violation": log_b,
            "violated": ratio_b > 1,
            "orders_of_violation": abs(log_b),
        },
        "scaling_analysis": {
            "bound_scale_factor": scale_bound,
            "average_density_scale_factor": scale_avg,
            "ratio_scale_factor": ratio_scale,
            "interpretation": (
                "Bound decreases as tau^{-4} = f^4, but average density decreases only as f. "
                "Net effect: ratio grows by f^{-3} = 1e9 when going from case (a) to case (b). "
                "Duty cycle makes violation WORSE, not better."
            ),
        },
        "verdict": {
            "modulation_bypasses_QI": False,
            "case_a_violation_orders": abs(log_a),
            "case_b_violation_orders": abs(log_b),
            "verdict_text": (
                f"Duty-cycle modulation does NOT bypass FR-QI. "
                f"Case (a): ~{abs(log_a):.0f} orders violation (same as static pulse). "
                f"Case (b): ~{abs(log_b):.0f} orders violation (WORSE due to tau^{-4} scaling). "
                "Time modulation does not bypass QI for instantaneous local energy density measure."
            ),
            "physical_fact": (
                "FR-QI applies to the light tube covering any random moment of pulse switch-on. "
                "Duty-cycle adds empty periods but does not reduce peak density during ON time. "
                "Literature: temporal modulation does NOT bypass QI for instantaneous measures."
            ),
        },
    }

    with open("/home/smboozha/portal_gun/research/r3_qi_modulated_a2.json", "w",
              encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2, default=str)
    print(f"  JSON: research/r3_qi_modulated_a2.json")
    print("=" * 100)

    return {
        "case_a_ratio": ratio_a,
        "case_a_log10": log_a,
        "case_b_ratio": ratio_b,
        "case_b_log10": log_b,
        "rho_avg": rho_avg,
        "eps_max_a": eps_max_a,
        "eps_max_b": eps_max_b,
    }


if __name__ == "__main__":
    result = main()
    print("\n--- KEY NUMBERS ---")
    for k, v in result.items():
        print(f"  {k}: {v}")
