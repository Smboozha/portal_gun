#!/usr/bin/env python3
# v1.0  — PASSAGE / ПРОХОД. Итоговый честный расчёт с ПРАВИЛЬНОЙ QI-формулой SI.
#
#   ПРАВИЛЬНАЯ форма квантового неравенства (Ford-Roman, SI):
#       |ρ|_max(τ) = ħ / (8π² c³ τ⁴)      [J/м³]
#   (ранее в проекте ошибочно использовалось ħc³/(8π²τ⁴), что давало
#    завышение допуска в c⁶ ≈ 7e49 — это ошибка, см. PASSAGE.md §CORRECTION)
#
#   Пан.1 — ЗАКОН QI:  |ρ|_max(τ) vs τ, окна экспозиции и ρ_need.
#   Пан.2 — ВМОСТЯЩАЯ ДИАГРАММА горловины Морриса–Торна, b(r)=r₀.
#   Пан.3 — ВЕРДИКТ: свет/квант/человек по квантовому неравенству.

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

H = 1.054571817e-34
C = 2.99792458e8
G = 6.67430e-11
K = 8 * np.pi**2


def qi(t: float) -> float:
    return H / (K * C**3 * t**4)


def rho_need(r0: float) -> float:
    return C**4 / (8 * np.pi * G * r0**2)


fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# ---------------------------------------------------------------- Панель 1
ax = axs[0]
tau = np.logspace(-44, -8, 400)
b = qi(tau)
ax.loglog(tau, b, "b-", lw=2.5, label=r"$|\rho|_{\max}(\tau)=\hbar/(8\pi^2 c^3 \tau^4)$")

for r0, col in ((0.1, "r"), (1.0, "m")):
    ax.axhline(rho_need(r0), color=col, ls="--", lw=1.5)
    ax.text(3e-43, rho_need(r0) * 1.6, f"rho_need(r0={r0} м) = {rho_need(r0):.1e} J/m3", color=col, fontsize=10)

markers = [
    (3.2e-27, "EW-шкала\n(квант)"),
    (3.34e-10, "r₀=0.1 м / c"),
    (3.34e-9, "r₀=1 м / c"),
    (6e-9, "человек\nL/c"),
]
for t, lab in markers:
    ax.plot(t, qi(t), "ko", ms=6, zorder=5)
    ax.annotate(lab, (t, qi(t)), textcoords="offset points",
                xytext=(18, 6), fontsize=9, arrowprops=dict(arrowstyle="-"))

ax.set_xlim(1e-44, 1e-8)
ax.set_ylim(1e-34, 1e115)
ax.set_xlabel(r"время квантования τ [с]")
ax.set_ylabel(r"$\rho$ [Дж/м³]")
ax.set_title("1) ЗАКОН QI (Ford-Roman, SI-вид)")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="lower right", fontsize=9)

# ---------------------------------------------------------------- Панель 2
ax = axs[1]
r0 = 1.0
r = np.linspace(r0, 4 * r0, 300)
u = np.sqrt(r / r0 - 1)
z = (2.0 / 3.0) * r0 * u**3
rpos = np.concatenate([r[::-1], r])
zpos = np.concatenate([-z[::-1], z])
ax.fill(rpos / r0, zpos / r0, color="lime", alpha=0.35, label="отрицательная\nэнергия (горловина)")
ax.plot(r / r0, z / r0, "g-", lw=2)
ax.plot(r / r0, -z / r0, "g-", lw=2)
ax.plot([-4, 4], [0, 0], "k-", lw=1)
ax.fill_between([-4, -1], -0.1, 0.1, color="gray", alpha=0.8, label="‘верхний’ лист")
ax.fill_between([1, 4], -0.1, 0.1, color="gray", alpha=0.8)
ax.plot([-r0, r0], [0, 0], "k-", lw=4)
ax.plot(0, 0, "go", ms=10)
ax.annotate("горловина r₀", (0, 0), textcoords="offset points",
            xytext=(-6, -26), fontsize=11, ha="center")
ax.set_xlim(-4.5, 4.5)
ax.set_ylim(-4.5, 4.5)
ax.set_xlabel("r / r₀")
ax.set_ylabel("z / r₀")
ax.set_title("2) ВЛОЖЕНИЕ ГОРЛОВИНЫ  b(r)=r₀")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

# ---------------------------------------------------------------- Панель 3
ax = axs[2]
ax.set_xlim(0, 3)
ax.set_ylim(-0.05, 1.05)
rows = [
    ("СВЕТ (фотон)", "τ→0  ⇒  |ρ|_max→∞", "РАЗРЕШЕНО", "lime", "реально строим (линзы Френеля)"),
    ("КВАНТ на атометровой щели\na≈1e-18 м, τ≈3e-27 с", "|ρ|_max=4.7e44 ≈ ρ_need(0.1 м)=4.8e44",
     "РОЗНИЦА ~1.0x", "gold", "[ПОСТУЛАТ]: домен на EW-шкале"),
    ("ЧЕЛОВЕК (L/c≈6нс)", "|ρ|_max=3.8e-29  vs  ρ_need=4.8e42", "ЗАПРЕЩЕНО", "tomato",
     "разрыв ~1.3e71 — доказанное QI"),
]
for i, (name, num, verdict, col, note) in enumerate(rows):
    y = 0.98 - i * 0.34
    ax.add_patch(Rectangle((0.1, y - 0.13), 2.8, 0.3, facecolor=col, edgecolor="k", lw=1.2))
    ax.text(0.16, y + 0.03, name, fontsize=11, va="center", fontweight="bold")
    ax.text(0.16, y - 0.1, num, fontsize=9, va="center")
    ax.text(2.9, y, verdict, fontsize=10, va="center", ha="right", fontweight="bold")
    ax.text(1.55, y - 0.1, note, fontsize=8.5, ha="center", color="dimgray", va="center")
ax.axis("off")
ax.set_title("3) ВЕРДИКТ КВАНТОВОГО НЕРАВЕНСТВА", fontsize=13)

fig.suptitle("PASSAGE: честная экономика прохода (ПРАВИЛЬНАЯ QI-формула SI, фактор c⁶)",
             fontsize=14, y=0.99)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("sim_passage.png", dpi=120)
print("OK -> sim_passage.png")

# ключевые числа на stdout
print()
print("ПРАВИЛЬНАЯ QI: |ρ|_max(τ) = ħ/(8π²c³τ⁴) = 4.957e-62/τ⁴  [J/м³]")
print(f"r0=0.1м: ρ_need={rho_need(0.1):.3e},  |ρ|_QI(r0/c)={qi(0.1/C):.3e}, разрыв={rho_need(0.1)/qi(0.1/C):.2e}")
print(f"человек: |ρ|_QI(6e-9с)={qi(6e-9):.3e}, разрыв={rho_need(1.0)/qi(6e-9):.2e}")
for r0 in (0.1, 1.0):
    lgh = (4.957e-62 / rho_need(r0)) ** 0.25 * C
    print(f"легальная длина экспозиции при ρ_need(r0={r0}м): c·τ = {lgh:.3e} м")
print(f"квантовое окно τ=3.2e-27с: |ρ|_max={qi(3.2e-27):.3e} (~ρ_need(0.1)=4.8e44, запас {qi(3.2e-27)/rho_need(0.1):.3f})")