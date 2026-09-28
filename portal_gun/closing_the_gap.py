#!/usr/bin/env python3
"""
CLOSING THE GAP — движок закрытия энергетического разрыва портала.

Честно перебирает ВСЕ известные масштабы длины a (зазор Казимировой полости)
и все реальные/обсуждаемые усиливающие механизмы, чтобы найти конфигурацию,
где разрыв между доступной отрицательной энергией и нужной ρ_need(r₀=0.1 м)
падает до ~1 порядка (≤10×).

Формулы (SI):
  ρ_Casimir(a) = π²ħc / (720 a⁴)          — плотность энергии Казимира
  ρ_need(r₀)   = c⁴ / (8πG r₀²)           — требование горловины Morris-Thorne
  QI-предел    = ħc³ / (8π² t₀⁴), t₀ = a/c — квантовое неравенство (микро)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================== КОНСТАНТЫ (SI) =================
hbar = 1.054571817e-34     # Дж·с
c = 2.99792458e8           # м/с
G = 6.67430e-11            # м³/(кг·с²)
l_P = 1.616255e-35         # планковская длина
pi2 = np.pi**2

# ============================== БАЗОВЫЕ ФОРМУЛЫ =================
def rho_casimir(a: float) -> float:
    """Плотность энергии Казимира при зазоре a (Дж/м³)."""
    return pi2 * hbar * c / (720.0 * a**4)

def rho_need(r0: float = 0.1) -> float:
    """Нужная плотность отрицательной энергии горловины r0 (Дж/м³)."""
    return c**4 / (8.0 * np.pi * G * r0**2)

def qi_bound(t0: float) -> float:
    """Квантовое неравенство: |ρ|_max при времени выборки t0."""
    return hbar * c**3 / (8.0 * pi2 * t0**4)

def gap(a: float, r0: float = 0.1) -> float:
    """Разрыв = нужное/доступное (≤1 — закрыто, 10 = один порядок)."""
    return rho_need(r0) / max(rho_casimir(a), 1e-300)

# ============================== ШКАЛЫ a =========================
SCALES = [
    # (a, название, маркировка)
    (3.0e-9,  "Литография EUV 2025",          "[РЕАЛЬНО]"),
    (1.0e-9,  "Литография 1 нм (передовая)",  "[РЕАЛЬНО]"),
    (3.0e-10, "Одиночные атомы (STM)",        "[РЕАЛЬНО]"),
    (2.4e-10, "Графен: межслойный зазор",     "[РЕАЛЬНО]"),
    (1.0e-10, "Ангстрем (хим. связь)",        "[РЕАЛЬНО]"),
    (5.3e-11, "Боровский радиус ЯБP:",         "[ФИЗИКА]"),
    (3.86e-13,"Комптоновская длина электрона","[ФИЗИКА]"),
    (1.0e-15, "Ядро (ферми)",                 "[ФИЗИКА]"),
    (2.0e-16, "КХД / 1 ГэВ",                  "[ФИЗИКА]"),
    (2.0e-18, "Электро-слабая / 100 ГэВ",     "[ФИЗИКА/ПОСТУЛАТ]"),
    (8.0e-19, "Вакуум Хиггса v=246 ГэВ",      "[ФИЗИКА/ПОСТУЛАТ]"),
    (1.6e-35, "Планковская длина",            "[ПОСТУЛАТ]"),
]

# ============================== РЫЧАГИ УСИЛЕНИЯ =================
# Мультипликативный вклад в ρ_Casimir. Только честные физические
# механизмы с подтверждёнными/разумными величинами.
LEVERS = {
    "ENZ-метаматериал n=1e-6 (заполнение)":     1e6,
    "Гиперболич. метаматериал (плотность мод)": 1e4,
    "Сжатый вакуум e^(2r), r=5":                2.0e4,
    "Резонатор Q=1e9 (добротность)":            1e9,
}

def order_zero(a: float) -> float:
    """Плотность Казимира без усиления."""
    return rho_casimir(a)

# отдельная обработка: QI запрещает БЕСКОНЕЧНО большое усиление —
# проверяем после каждого умножения.

# ============================== ГЛАВНАЯ ТАБЛИЦА =================
def main():
    rn = rho_need()
    print("="*100)
    print("CLOSING THE GAP — полный перебор масштабов зазора a")
    print("  ρ_need(0.1 м) = %.3e Дж/м³" % rn)
    print("="*100)
    print()
    print(f"{'зазор a, м':>12} | {'ρ_Casimir':>10} | {'разрыв':>8} | "
          f"{'QI-макс':>10} | {'OK?':>4} | {'масштаб':<28} | статус")
    print("-"*100)

    results = []
    for a, name, status in SCALES:
        rC = rho_casimir(a)
        gp = gap(a)
        t0 = a / c
        qm = qi_bound(t0)
        ok = "✓" if rC <= qm else "✗!"
        results.append((a, rC, gp, qm, name, status))
        print(f"{a:12.2e} | {rC:10.2e} | {gp:8.1e} | {qm:10.1e} | {ok:>4} | "
              f"{name:<28} | {status}")

    # ---------- точная критическая точка (разрыв = 1) ----------
    a_crit = (pi2 * hbar * c / (720.0 * rn)) ** 0.25
    print()
    print("="*100)
    print("КРИТИЧЕСКАЯ ТОЧКА: зазор, при котором разрыв = 1")
    print(f"  a_крит = {a_crit:.3e} м = {a_crit*1e18:.3f} зептометров (=аттометр)")
    print(f"  ρ_Casimir(a_крит) = {rho_casimir(a_crit):.3e} Дж/м³ ≈ ρ_need")
    print(f"  разрыв = {gap(a_crit):.4f}")
    a_EW = 2.0e-18
    a_H  = 8.0e-19
    print()
    print("  Сравнение с известными шкалами:")
    print(f"    Электро-слабая (100 ГэВ):  a={a_EW:.0e} м → ρ_C={rho_casimir(a_EW):.3e}, "
          f"разрыв={gap(a_EW):.2f}")
    print(f"    Вакуум Хиггса (246 ГэВ):   a={a_H:.0e} м → ρ_C={rho_casimir(a_H):.3e}, "
          f"разрыв={gap(a_H):.2f}")
    print(f"    Геометр. среднее √(l_P·r₀)= {np.sqrt(l_P*0.1):.3e} м")
    print()

    # ---------- рычаги усиления от атомного пола ----------
    print("="*100)
    print("РЫЧАГИ УСИЛЕНИЯ от атомного масштаба (a₀=2.4e-10 м):")
    a0 = 2.4e-10
    print(f"  базовый ρ_Casimir(2.4e-10 м) = {rho_casimir(a0):.2e} Дж/м³")
    print(f"  разрыв без рычагов = {gap(a0):.1e}")
    combo = order_zero(a0)
    combos = []
    levers_total = 1.0
    for i, (name, mult) in enumerate(LEVERS.items()):
        combo *= mult
        levers_total *= mult
        combos.append((name, mult, combo))
        print(f"    + {name:<45} ×{mult:8.0e} → ρ={combo:.2e}, разрыв={rn/max(combo,1e-300):.1e}")

    # ---------- закрытие: что даёт сумма всех рычагов ----------
    print()
    print("  Сума всех рычагов сверху от a₀=2.4e-10 м:")
    combo_top = order_zero(a0) * levers_total
    print(f"    ρ_max = {combo_top:.2e} Дж/м³ → разрыв = {rn/combo_top:.1e}")

    # ---------- вывод: единственное закрытие ----------
    print()
    print("="*100)
    print("ВЫВОД (честно):")
    print("  1) Все рычаги от атомного масштаба (ENZ+метаматериал+сжатие+Q)")
    print(f"     дают ρ_max ≈ {combo_top:.2e} Дж/м³ → разрыв остаётся ~{rn/combo_top:.0e} (10 порядков).")
    print("  2) Единственная точка пересечения разрыва = 1 — зазор a ≈ 1e-18 м")
    print("     (электро-слабая / вакуум-Хиггса шкала, 100–250 ГэВ).")
    print("  3) ЗАМКНУТАЯ ФОРМУЛА: a_крит = (π³/90)^(1/4) · √(l_P·r₀) = 0.766·√(l_P·r₀)")
    print("     для r₀=0.1 м → a_крит = 9.74e-19 м ≈ 1 аттометр.")
    print("  4) Квантовое неравенство этому НЕ мешает: на уровне полости")
    print("     t₀=a/c, QI-макс ~1e94 ≫ 1e44 (запас 50 порядков).")
    print("  5) Настроечный рычаг: менять r₀ горловины — разрыв масштабируется")
    print("     как (r₀)⁰.⁵: большая горловина требует бóльший зазор при том же ρ.")
    print("="*100)

    # ============================== ГРАФИК =================
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_facecolor('#0a0a2e')
    for spine in ax.spines.values():
        spine.set_color('#333366')
    ax.tick_params(colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')

    a_range = np.logspace(-35, -9, 2000)
    gp_range = np.array([gap(a) for a in a_range])
    rC_range = np.array([rho_casimir(a) for a in a_range])

    ax.loglog(a_range*1e9, gp_range, color='#00ff41', linewidth=2,
              label='разрыв ρ_need/ρ_Casimir')
    ax.axhline(1, color='white', linestyle=':', alpha=0.6)
    ax.axhline(10, color='#ffaa00', linestyle='--', linewidth=1.5,
               label='разрыв = 1 порядок (10×)')

    # Точки шкал
    for a, name, _ in SCALES:
        g = gap(a)
        color = '#ff5555' if g > 10 else '#00ffaa'
        ax.scatter([a*1e9], [g], color=color, s=60, zorder=5,
                   edgecolors='white', linewidths=0.5)
        if a >= 1e-18:
            ax.annotate(name.split('(')[0].strip(),
                        xy=(a*1e9, g), xytext=(a*2.5e-9+1e-10, g*2),
                        color='white', fontsize=8,
                        arrowprops=dict(arrowstyle='->', color='#8888cc', lw=0.6))
        elif a >= 4e-14:
            ax.annotate(name.split('(')[0].strip(), xy=(a*1e9, g),
                        xytext=(a*1e9, g*1.4), color='white', fontsize=8)

    ax.axvline(1e-18*1e9, color='#ffaa00', linestyle=':', alpha=0.7)
    ax.text(1e-18*1e9*1.3, 3e2, 'a ≈ 1e-18 м\n(зазор → разрыв 1)',
            color='#ffaa00', fontsize=9)

    ax.set_xlabel('Зазор Казимировой полости a (нм)', fontsize=12)
    ax.set_ylabel('Разрыв = ρ_нужно / ρ_Казимир', fontsize=12)
    ax.set_title('ЗАКРЫТИЕ РАЗРЫВА: зависимость от зазора полости', fontsize=14)
    ax.grid(True, which='both', alpha=0.2)
    ax.legend(fontsize=9, facecolor='#0a0a2e', edgecolor='#333366', labelcolor='white')
    ax.set_ylim(1e-2, 1e40)
    plt.tight_layout()
    plt.savefig('/home/smboozha/portal_gun/sim_gap_closure.png', dpi=150,
                bbox_inches='tight', facecolor='#0a0a2e')
    plt.close()
    print("  ✓ График: sim_gap_closure.png")


if __name__ == "__main__":
    main()