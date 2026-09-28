#!/usr/bin/env python3
"""
PORTAL FROM JUNK — инженерный расчётчик.

Оптический релейный портал из мусора:
  две линзы Френеля (из мертвых оверхед-проекторов / заднепроекционных ТВ)
  + перископная складка из зеркал → ВЫ видите сквозь кольцо реальную сцену.

Это ПОДЛИННАЯ физика: свет реально проходит от точки B к точке A.
Проход материи — отдельная статья ([ПОСТУЛАТ], см. JUNK_PORTAL.md).

Расчёт:
  * геометрия реле (f, D, расстояние, поле зрения, увеличение)
  * кольцо из светодиодов / неона (длина, кол-во, питание)
  * датчик приближения (PIR) — мощность кольца
  * схема укладки + PNG-рендер
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch

PI = np.pi

# ============================== ВХОДНЫЕ ПАРАМЕТРЫ =================
# Оверхед-проекторская линза Френеля: D≈0.3 м, f≈0.35 м
JUNK_FRESNEL = {"name": "Френель из оверхед-проектора", "D": 0.30, "f": 0.35}
# Заднепроекционный ТВ: линза D≈0.42 м, f≈0.55 м (+ есть зеркало-складка!)
JUNK_RPTV    = {"name": "Френель из заднепроекционного ТВ", "D": 0.42, "f": 0.55}

def relay_geometry(lens: dict) -> dict:
    """Геометрия двухлинзового реле с совмещением фокальных плоскостей."""
    D, f = lens["D"], lens["f"]
    M = 1.0                     # одинаковые линзы → увеличение 1
    sep = 2 * f                 # между линзами (фокальные плоскости совпадают)
    half_fov = np.degrees(np.arctan2(D / 2.0, f))
    # Пуппиль (размер «окна», видимого наблюдателем)
    pupil = D
    # масштаб расстояния линза-зеркало-складка для перископа
    mirror_fold = sep / 2.0
    return {
        "D": D, "f": f, "M": M, "sep": sep,
        "half_fov": half_fov, "fov": 2 * half_fov,
        "pupil": pupil, "mirror_fold": mirror_fold,
    }

def led_ring(D_hoop: float, led_per_m: float = 30.0) -> dict:
    """Кольцо: окружность, число LED, мощность, питание."""
    circ = PI * D_hoop
    n_led = int(led_per_m * circ)
    mul = (led_per_m / 60.0) ** 0.0   # placeholder для неизменяемости
    # 5050 white-green: ~0.3 Вт/LED при питании 12В
    p_per_led = 0.30
    power = n_led * p_per_led
    current = power / 12.0
    return {"circ_m": circ, "n_led": n_led, "power_W": power,
            "current_A": current, "led_per_m": led_per_m}

def main():
    print("="*100)
    print("PORTAL FROM JUNK — расчёт оптического реле и кольца")
    print("="*100)

    optics = {}
    for lens in (JUNK_FRESNEL, JUNK_RPTV):
        g = relay_geometry(lens)
        optics[lens["name"]] = g
        print()
        print(f"### {lens['name']}")
        for k, v in g.items():
            print(f"    {k:<14} = {v:.3g}" if isinstance(v, float) else
                  f"    {k:<14} = {v}")
        r = led_ring(g["D"])
        print(f"    Кольцо по линзе: L={r['circ_m']:.2f} м, "
              f"LED={r['n_led']}, P={r['power_W']:.0f} Вт, "
              f"I={r['current_A']:.2f} А @12В")

    print()
    print("="*100)
    print("ИНСТРУКЦИЯ КРАТКО (подробно — JUNK_PORTAL.md):")
    print("  1) Две линзы A и B ставятся друг к другу на оси.")
    print("     Расстояние между ними ≈ 2f (у оверхеда ≈ 0.7 м).")
    print("     Между линзами — перископная складка из 2 зеркал (сохраняет ориентацию).")
    print("     Можно разносить складкой на любое расстояние (коллимированный пучок).")
    print("  2) Вокруг каждой линзы — зелёное кольцо (LED/неон).")
    print("  3) Смотреть в линзу A = видеть сцену перед B (реальный свет).")
    print("="*100)

    # ============================== PNꖅ-совместимый РЕНДЕР ================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7),
                                   gridspec_kw={'width_ratios': [1.15, 1]})
    fig.patch.set_facecolor('#0a0a2e')

    # ---- левый: схема портала из мусора (аксонометрия) ----
    ax = ax1
    ax.set_facecolor('#0a0a2e')
    ax.axis('off')
    ax.set_xlim(-1, 26)
    ax.set_ylim(-2, 20)

    # трубный "мост" между комнатами
    bridge = Rectangle((2, 1), 16, 4, color='#223', alpha=0.9,
                       edgecolor='#666', linewidth=1)
    ax.add_patch(bridge)
    ax.text(10, 3.2, 'СВЕТОВОЙ МОСТ\n(труба/короб/потолок)',
            color='#ccc', fontsize=9, ha='center', va='center')

    # стены
    for xw, name in [(0, 'КОМНАТА B'), (18, 'КОМНАТА A')]:
        ax.plot([xw, xw], [-1, 17], color='#666', linewidth=3)
        ax.text(xw + 0.4, 16, name, color='#888', fontsize=10, rotation=90)

    # портал B (лево)
    yB = 5
    ringB = Circle((2, yB), 2.2, color='#00ff41', alpha=0.15)
    ax.add_patch(ringB)
    ringB2 = Circle((2, yB), 2.0, fill=False, edgecolor='#00ff41', linewidth=4)
    ax.add_patch(ringB2)
    lensB = Circle((2, yB), 1.4, fill=False, edgecolor='#88ccff', linewidth=3,
                   linestyle='--')
    ax.add_patch(lensB)
    ax.text(2, yB + 2.8, 'ПОРТАЛ B\n(линза Френеля)', color='#00ff41',
            fontsize=9, ha='center')
    # человек у B
    ax.plot([2 - 3.4, 2 - 3.4], [yB - 1.6, yB + 0.8], color='white', linewidth=3)
    ax.plot([2 - 4.6, 2 - 3.4, 2 - 2.2], [yB + 0.8, yB + 2.2, yB + 0.8],
            color='white', linewidth=3)
    ax.plot([2 - 4.0, 2 - 3.4, 2 - 2.8, 2 - 2.2], [yB - 1.6, yB - 3.2, yB - 3.2, yB - 1.6],
            color='white', linewidth=3)

    # портал A (право)
    yA = 5
    ringA = Circle((18, yA), 2.2, color='#00ff41', alpha=0.15)
    ax.add_patch(ringA)
    ringA2 = Circle((18, yA), 2.0, fill=False, edgecolor='#00ff41', linewidth=4)
    ax.add_patch(ringA2)
    lensA = Circle((18, yA), 1.4, fill=False, edgecolor='#88ccff', linewidth=3,
                   linestyle='--')
    ax.add_patch(lensA)
    ax.text(18, yA + 2.8, 'ПОРТАЛ A\n(линза Френеля)', color='#00ff41',
            fontsize=9, ha='center')
    # глаз наблюдателя у A
    ax.text(18 + 3.0, yA + 0.4, 'ГЛАЗ', fontsize=12, ha='center',
            color='white', fontweight='bold')
    ax.text(18 + 3.0, yA - 1.6, 'смотрю в A = вижу сцену перед B',
            color='#ffaa44', fontsize=8, ha='center')

    # луч через мост
    ax.annotate('', xy=(2, yA - 0.2), xytext=(2 + 1.6, yA - 0.2),
                arrowprops=dict(arrowstyle='->', color='#00ff41', lw=2))
    ax.annotate('', xy=(18 - 1.6, yA - 0.2), xytext=(16.4, yA - 0.2),
                arrowprops=dict(arrowstyle='->', color='#00ff41', lw=2))
    ax.text(10, yA + 2.4, '→→→ реальный свет B → A ←←←', color='#00ff41',
            fontsize=10, ha='center', fontweight='bold')

    # ---- правый: кольцо и питание (электрика) ----
    ax = ax2
    ax.set_facecolor('#0a0a2e')
    ax.axis('off')
    ax.set_xlim(-3, 23)
    ax.set_ylim(0, 20)

    # Кольцо LED
    theta = np.linspace(0, 2*np.pi, 200)
    ring = plt.Circle((10, 10), 5, color='#003300', alpha=0.5)
    ax.add_patch(ring)
    ax.plot(10 + 5*np.cos(theta), 10 + 5*np.sin(theta),
            color='#00ff41', linewidth=4)
    ax.plot(10 + 5.3*np.cos(theta), 10 + 5.3*np.sin(theta),
            color='#00ff41', alpha=0.3, linewidth=2)
    for ang in np.linspace(0, 2*np.pi, 24, endpoint=False):
        x = 10 + 4.6*np.cos(ang)
        y = 10 + 4.6*np.sin(ang)
        ax.plot(x, y, color='#88ff88', marker='o', markersize=3)

    ax.text(10, 17.5, 'ЗЕЛЁНОЕ КОЛЬЦО\n(светодиодная лента 5050 или неон)',
            color='#00ff41', fontsize=11, ha='center', fontweight='bold')
    ax.text(10, 2.2, f"L={PI*0.3:.2f} м · {int(30*PI*0.3)} шт × 0.3 Вт · ~9 Вт\n"
                     "питание 12 В от БП-мусора (старый роутер/ноутбук)",
            color='#aaa', fontsize=9, ha='center')

    # цепь питания
    ax.plot([10, 10], [5, 7.4], color='#ff4444', linewidth=2)
    ax.plot([10, 10], [12.6, 15], color='#ff4444', linewidth=2)

    psu = plt.Rectangle((8.3, 6.0), 3.4, 1.6, color='#552', edgecolor='#ff8844')
    ax.add_patch(psu)
    ax.text(10, 6.8, 'БП 12 В', color='#ff8844', fontsize=9,
            ha='center', va='center', fontweight='bold')

    pir = plt.Rectangle((8.3, 14.4), 3.4, 1.6, color='#225', edgecolor='#88aaff')
    ax.add_patch(pir)
    ax.text(10, 15.2, 'PIR-датчик\n(мусор: сигналка/лампа)',
            color='#88aaff', fontsize=8, ha='center', va='center')

    okr = plt.Rectangle((13.5, 6.0), 3.4, 1.6, color='#252', edgecolor='#44ff88')
    ax.add_patch(okr)
    ax.text(15.2, 6.8, 'транзистор /\nреле (мусор)', color='#44ff88', fontsize=8,
            ha='center', va='center')
    ax.plot([11.7, 13.5], [6.8, 6.8], color='#44ff88', linewidth=2)
    ax.plot([10.6, 10.6], [6.0, 5.4], color='#ff4444', linewidth=2)
    ax.plot([10.6 + 0.15, 15.2], [5.4, 5.4], color='#ff4444', linewidth=2)
    ax.plot([15.2, 15.2], [5.4, 4.2], color='#ff4444', linewidth=2)
    ax.plot(15.2 + np.cos(theta)*0.4, 4.2 + np.sin(theta)*0.4,
            color='#ff4444', linewidth=2)
    ax.text(15.7, 4.2, '→ лента', color='#ff4444', fontsize=8)

    ax.set_title('Схема электрики кольца (всё из мусора)', color='white',
                 fontsize=12, pad=12)

    plt.tight_layout()
    plt.savefig('/home/smboozha/portal_gun/sim_junk_portal_layout.png', dpi=140,
                bbox_inches='tight', facecolor='#0a0a2e')
    plt.close()
    print("  ✓ Схема: sim_junk_portal_layout.png")

if __name__ == "__main__":
    main()