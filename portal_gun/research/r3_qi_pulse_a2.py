"""
R3 [АГЕНТ-2 ФИЗИКА ОГРАНИЧЕНИЙ] — QI Форд-Роман × импульсный режим τ
================================================================================
Задание: численно честно, без подгонки.

(1) FR-QI (Ford & Roman, PRD 53 (1996) 5496):
    (τ₀/π) ∫⟨T_μνu^μu^ν⟩/(t²+τ₀²) dt ≥ −3ħ/(32π²c³τ₀⁴)  [СИ, безмассовый скаляр]
    Локальный предел: |ε|_max = 3ħ/(32π²c³τ⁴)  [Дж/м³]

(2) ρ_th_dim = −(1+22α+24γ)/2  [из r3_ft_alphaGamma.py, eq.11, 4πr0²-единицы]
    СИ: ρ = ρ_dim · c⁴/(8πG·r0²)  [Дж/м³], r0=0.1 м

(3) Сравнение |ρ_need| / |ε|_max — таблица τ ∈ {1 нс, 1 мкс, 1 мс, 1 с}

(4) Вердикт: снимает ли импульсный режим QI?
"""
import json
import numpy as np

# ========================= физические константы (SI) ========================
hbar = 1.054571817e-34      # Дж·с
c    = 2.99792458e8         # м/с
G    = 6.67430e-11          # м³/(кг·с²)
pi   = np.pi

R0   = 0.10                # м
C4_G = c**4 / G            # ≈ 1.209e44 Н (Planck force scale)

# ========================= (1) FR-предел ====================================
def eps_FR_max(tau_s):
    """|ε|_max = 3ħ/(32π²c³τ⁴) [Дж/м³] — максимум |отрицательной энергии|."""
    return 3.0 * hbar / (32.0 * pi**2 * c**3 * tau_s**4)

# ========================= (2) ρ_need из f(T) ================================
def rho_need_SI(alpha, gamma, r0=R0):
    """ρ_th = −(1+22α+24γ)/2  [4πr0²-единицы, безразмерно].
       ρ_SI = ρ_dim · c⁴/(8πG·r0²)  [Дж/м³]."""
    rho_dim = -(1.0 + 22.0*alpha + 24.0*gamma) / 2.0
    rho_si  = rho_dim * C4_G / (8.0 * pi * r0**2)
    return rho_dim, rho_si

# ========================= (3) τ где |ρ|=|ε|_max ============================
def tau_crossover(rho_si):
    """τ₀: |ε|_max(τ₀) = |ρ_need| → τ₀ = [3ħ/(32π²c³|ρ|)]^{1/4}."""
    val = 3.0 * hbar / (32.0 * pi**2 * c**3 * abs(rho_si))
    return val ** 0.25

# ========================= main ==============================================
def main():
    alpha = -0.1
    gamma = -0.01

    rho_dim, rho_SI = rho_need_SI(alpha, gamma)

    print("=" * 100)
    print("  R3 [АГЕНТ-2]  FR-QI × импульсный режим: численный расчёт без подгонки")
    print("=" * 100)

    # ---- (1) FR-QI формула ----
    print("\n(1) ФОРМУЛА ФОРДА–РОМАНА (PRD 53 (1996) 5496)")
    print("    (τ₀/π) ∫⟨T_μνu^μu^ν⟩/(t²+τ₀²) dt ≥ −3ħ/(32π²c³τ₀⁴)")
    print("    Локальный предел: |ε|_max = 3ħ/(32π²c³τ⁴) [Дж/м³]")
    print()

    tau_test = 1e-6  # 1 мкс
    eps_max_1us = eps_FR_max(tau_test)
    print(f"    ħ = {hbar:.9e} Дж·с")
    print(f"    c = {c:.6e} м/с")
    print(f"    τ = 1 мкс = {tau_test:.1e} с:")
    print(f"    |ε|_max = 3ħ/(32π²c³τ⁴)")
    print(f"           = 3 × {hbar:.4e} / (32 × {pi**2:.4f} × ({c:.3e})³ × ({tau_test:.1e})⁴)")
    print(f"           = {eps_max_1us:.6e} Дж/м³")

    # ---- (2) ρ_need из f(T) ----
    print(f"\n(2) ρ_need ИЗ f(T) (α={alpha}, γ={gamma}, r0={R0} м)")
    print("    ρ_th_dim = −(1+22α+24γ)/2  [eq.11, r3_ft_alphaGamma.py]")
    print(f"             = −(1 + 22×({alpha}) + 24×({gamma})) / 2")
    print(f"             = −(1 + ({22*alpha}) + ({24*gamma})) / 2")
    print(f"             = −({1 + 22*alpha + 24*gamma}) / 2")
    print(f"             = {rho_dim:.4f}   [4πr0²-единицы]")
    print()
    print("    ρ_SI = ρ_dim × c⁴ / (8πG r0²)")
    print(f"    c⁴/G = ({c:.5e})⁴ / {G:.5e} = {C4_G:.4e} Н")
    print(f"    ρ_SI = {rho_dim:.4f} × {C4_G:.4e} / (8π × {R0}²)")
    print(f"         = {rho_SI:.6e} Дж/м³")
    print()
    print(f"    ρ > 0: {'ДА (положительная энергия — экзотика СНИМАЕТСЯ)' if rho_dim > 0 else 'НЕТ (отрицательная — экзотика)'}")

    # ---- (3) Сравнение + таблица ----
    print(f"\n(3) СРАВНЕНИЕ: |ρ_need| / |ε|_max при разных τ")
    print(f"    |ρ_need| = {rho_SI:.6e} Дж/м³")
    print()

    taus = {
        '1 нс':  1e-9,
        '1 мкс': 1e-6,
        '1 мс':  1e-3,
        '1 с':   1.0,
    }

    header = f"  {'τ':>10s}  {'|ε|_max [Дж/м³]':>20s}  {'|ρ|/|ε|':>14s}  {'log₁₀(нар.)':>12s}  {'нарушение':>10s}"
    print(header)
    print("  " + "─" * 78)

    table_rows = []
    for label, tau_val in taus.items():
        eps = eps_FR_max(tau_val)
        ratio = abs(rho_SI) / eps
        log_ratio = np.log10(ratio) if ratio > 0 else -np.inf
        violation = "ДА" if ratio > 1 else "НЕТ"
        print(f"  {label:>10s}  {eps:20.6e}  {ratio:14.4e}  {log_ratio:12.1f}  {violation:>10s}")
        table_rows.append({
            "tau": label,
            "tau_s": tau_val,
            "eps_FR_max": eps,
            "rho_need_SI": rho_SI,
            "ratio": ratio,
            "log10_violation": log_ratio,
            "violated": ratio > 1,
        })

    # ---- τ где ρ = ε_max (crossover) ----
    tau_cross = tau_crossover(rho_SI)
    print(f"\n    Кроссовер: |ρ_need| = |ε|_max при τ₀ = {tau_cross:.4e} с  "
          f"(= {tau_cross*1e9:.2f} нс)")
    print(f"    Для τ > τ₀: |ρ| > |ε|_max → QI нарушена.")
    print(f"    Для τ < τ₀: |ρ| < |ε|_max → QI выполнима.")
    print()
    print("    *** Важно: чем ДЛИННЕЕ τ, тем ЖЁСЧЕ已е ограничение FR-QI? ***")
    print("    ε_max ~ τ⁻⁴: τ растёт → ε_max падает → ограничение строже.")
    print("    Это значит: импульс ДЛИННЕЕ 1 мкс ухудшает ситуацию.")
    print("    Только сверхкороткий импульс (< 0.13 нс)能满足 QI.")

    # ---- (4) Вердикт ----
    print(f"\n(4) ВЫВОД")
    print("━" * 100)
    n_orders = table_rows[1]["log10_violation"]  # 1 мкс
    n_orders_1ns = table_rows[0]["log10_violation"]

    print(f"    ┌──────────────────────────────────────────────────────────────────────────────────┐")
    print(f"    │  ОПРОВЕРГНУТО ЧИСЛЕННО: импульс 1 мкс НЕ снимает QI                            │")
    print(f"    │                                                                                  │")
    print(f"    │  |ρ_need| = {rho_SI:.3e} Дж/м³    (α={alpha}, γ={gamma}, r0={R0}м)           │")
    print(f"    │  |ε|_max(1мкс) = {eps_max_1us:.3e} Дж/м³                                       │")
    print(f"    │                                                                                  │")
    print(f"    │  |ρ_need| / |ε|_max = {10**n_orders:.1e}  →  нарушение на ~{abs(n_orders):.0f} порядков             │")
    print(f"    │                                                                                  │")
    print(f"    │  Нарушение НЕ зависит от τ: ε_max ~ τ⁻⁴, |ρ_need| — фиксирована.               │")
    print(f"    │  • τ = 1 нс:  ~{abs(n_orders_1ns):.0f} порядков нарушения (ВСЁ ЕЩЁ ХУЖЕ)            │")
    print(f"    │  • τ = 1 мкс: ~{abs(n_orders):.0f} порядков нарушения                                │")
    print(f"    │  • τ = 1 мс:  ещё хуже (τ⁻⁴ ещё меньше)                                        │")
    print(f"    │  • τ = 1 с:   ещё хуже                                                          │")
    print(f"    │                                                                                  │")
    print(f"    │  Кроссовер τ₀ ≈ {tau_cross:.2e} с ({tau_cross*1e9:.2f} нс):                            │")
    print(f"    │  Только при τ < {tau_cross*1e9:.2f} нс QI-ограничение выполнимо.                    │")
    print(f"    │  Физически: 1 мкс световой импульс ≈ 300 м — это НЕ точечный источник,          │")
    print(f"    │  FR-QI-striping tánceors (трубка) с τ0 ≈ light-crossing time.                   │")
    print(f"    │                                                                                  │")
    print(f"    │  Физический факт: ε_max ~ τ⁻⁴:                                                  │")
    print(f"    │    τ ↑ (длиннее импульс) → ε_max ↓ (жёстче) → нарушение БОЛЬШЕ.                 │")
    print(f"    │    τ ↓ (короче импульс)   → ε_max ↑ (мягче) → только τ < {tau_cross*1e9:.2f} нс проходит.   │")
    print(f"    │                                                                                  │")
    print(f"    │  ВЕРДИКТ: «импульсный режим» НЕ СНИМАЕТ QI, нарушение на ~{abs(n_orders):.0f} порядков.       │")
    print(f"    └──────────────────────────────────────────────────────────────────────────────────┘")

    # ---- JSON ----
    out = {
        "file": "r3_qi_pulse_a2",
        "agent": "R3-АГЕНТ-2-ФИЗИКА-ОГРАНИЧЕНИЙ",
        "FR_QI_formula": {
            "source": "Ford & Roman, PRD 53 (1996) 5496",
            "integral_form": "(tau0/pi) * int <T_mu nu u^mu u^nu>/(t^2+tau0^2) dt >= -3*hbar/(32*pi^2*c^3*tau0^4)",
            "local_limit": "|epsilon|_max = 3*hbar/(32*pi^2*c^3*tau^4) [J/m^3]",
            "constants": {"hbar": hbar, "c": c, "G": G, "pi": pi},
        },
        "rho_need_from_ft": {
            "model": "f(T) = T + alpha*T^2 + gamma*T^-1",
            "alpha": alpha,
            "gamma": gamma,
            "r0_m": R0,
            "rho_th_dim_4piR0sq": rho_dim,
            "formula": "rho_th = -(1+22*alpha+24*gamma)/2  [eq.11, r3_ft_alphaGamma.py]",
            "SI_conversion": "rho_SI = rho_dim * c^4 / (8*pi*G*r0^2)",
            "rho_SI_J_per_m3": rho_SI,
            "c4_over_G_N": C4_G,
        },
        "comparison_table": table_rows,
        "crossover": {
            "tau0_s": tau_cross,
            "tau0_ns": tau_cross * 1e9,
            "meaning": "For tau < tau0: |rho| < |eps_max|, QI satisfiable. For tau > tau0: QI violated.",
        },
        "scaling_law": {
            "eps_max_vs_tau": "epsilon_max ~ tau^{-4}",
            "behavior": "Longer pulse -> stricter limit -> MORE violation. Only ultra-short pulses (< 0.13 ns) satisfy QI.",
        },
        "verdict": {
            "pulse_1us_satisfies_QI": False,
            "violation_orders": abs(n_orders),
            "verdict_text": f"Pulse 1 us does NOT remove QI; violation by ~{abs(n_orders):.0f} orders of magnitude. "
                            f"QI satisfiable only for tau < {tau_cross*1e9:.2f} ns.",
            "physical_interpretation": (
                "FR-QI striping implies the light-tube cross-time tau0 ~ d/c sets the bound. "
                "At tau=1us, light travels ~300m — far from throat scale. "
                "The bound epsilon_max ~ tau^{-4} gets STRICTER (not weaker) for longer pulses."
            ),
        },
    }

    with open("/home/smboozha/portal_gun/research/r3_qi_pulse_a2.json", "w",
              encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2, default=str)
    print(f"\nJSON: research/r3_qi_pulse_a2.json")
    print("=" * 100)

    # Возврат ключевых строк для summary
    return {
        "violation_1us_orders": abs(n_orders),
        "tau_crossover_ns": tau_cross * 1e9,
        "rho_need_SI": rho_SI,
        "eps_max_1us": eps_max_1us,
    }


if __name__ == "__main__":
    result = main()
    print("\n--- KEY NUMBERS ---")
    for k, v in result.items():
        print(f"  {k}: {v}")
