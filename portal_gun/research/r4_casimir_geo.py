"""
R4 [АГЕНТ-4, ГЕОМЕТРИЯ] — ЗНАК И ВЕЛИЧИНА ЭНЕРГИИ КАЗИМИРА ПО ТОПОЛОГИИ ОБОЛОЧКИ
==========================================================================
Вопрос: какая топология/геометрия проводящей границы даёт ОТРИЦАТЕЛЬНУЮ
энергию/плотность Казимира (ρ_vac<0 — поддержка NEC для макро-горловины)?

ЗАФИКСИРОВАНО ЧИСЛЕННО (эта сессия [ЖОН], scratch_*):
  · Окружность (1D-кольцо): E(τ) = R/τ² − 1/(12R) + τ²/(240R³)−…
    τ-экстракция базисом {τ⁻²,τ⁰,τ²,τ⁴,τ⁶,τ⁸}, τ/R∈[0.15,1.4], 16 т.:
    −1/(12R) с относительной ошибкой ~8e-9.  ЗНАК МИНУС.
  · Параллельные пластины (2D, скаляр на площадь): константа −π²/(1440a³)
    базисом {τ⁻⁴..τ⁺⁸}, τ/a∈[0.02,0.30], column-scaling lstsq:
    относительная ошибка ~4e-6.  ЗНАК МИНУС (сокращение больших членов
    между тремя слагаемыми суммы — константа рождается из отмены ~π²/480a³).
  · Сфера (PEC, Milton/Boyer): E = +(0.092353/2)/a = +0.0461765/a.
    Реализация R_l-интегралом в этом файле → ~0.01%.  ЗНАК ПЛЮС.
  · Бесконечный цилиндр (DeRaad–Milton 1981, цит. по Straley–Kolomeisky
    PRA 90, 012514 (2014)): E/L = −0.01356/a² на единицу длины.  ЗНАК МИНУС.
  · ТОР (модель «wrapped-waveguide»): ω=√((p/R)²+γ²), p∈ℤ, γ — нули
    J_l/J'_l (scipy jn_zeros/jnp_zeros). Лидирующий коэффициент
    E(τ)·τ⁴ → 3a²R (Вейль, проверен ТОЧНО: 3.0e-8 vs 3a²R=3.0e-8 — это
    закрывает прежнее расхождение ×~1.85, причиной был double-count в старом
    билдере нулей). ПРЯМОЕ извлечение константы E0 численно НЕСХОДИМО
    (обусловленность ~1e5: остаток после вычитания расходящихся τ⁻⁴..τ⁻¹
    порядка 1e5×сама константа, знак неустойчив при смене базиса/счёта).
  · Литературный якорь для ТОРА (Balian–Duplantier, SkK-2014): полная
    энергия Казимира C при минорном радиусе 1 ОТРИЦАТЕЛЬНА на всех R/a:
        R/a=2:−0.11  3:−0.16  4:−0.22  5:−0.28  6:−0.33  8:−0.44  10:−0.56
    с законом C ≈ −0.009·2πR (≈пропорциональна длине кольца).

ВЫВОД (топология ↔ знак):
  · выпуклая замкнутая сфера (g=0) — E>0, РЕПУЛЬСИЯ: НЕ держит горловину;
  · цилиндр / кольцо (тор) / компактифицированные направления — E<0,
    ρ_vac<0: ПОДДЕРЖИВАЕТ условие ρ<0 (по знаку — кандидат для стабилизации);
  · «пузырь» Казимира ρ = −π²ℏc/(720·R_b⁴) < 0 — самостоятельный источник
    отрицательной плотности (см. r4_bubound.py).
  НО: все модули — масштаб ℏc (~1e-26 Дж·м), т.е. макроскопически ничтожно;
    для оболочки на радиусе r0=0.10 м см. энергетические сравнения ниже.

МАСШТАБ ИЗМЕРЕНИЙ: ℏc = 3.161526472e-26 Дж·м; ħ=c=1 в естественных формулах,
перевод в SI умножением на ℏc и делением на характерную длину.
"""
import json
import numpy as np
import mpmath as mp
from scipy import integrate
from scipy.special import ive, kve, jn_zeros, jnp_zeros

hbar   = 1.054571817e-34
c      = 2.99792458e8
G      = 6.67430e-11
HBARC  = hbar * c          # 3.161526472e-26 Дж·м
pi     = np.pi

mp.mp.dps = 20
Z0 = 0.01                  # разделение: mpmath [0,Z0], scipy [Z0,∞)

# ======================================================================
# 1. СФЕРА (Milton/Boyer): E = −(1/2a)·(Σ_l R_l − 3/32)
# ======================================================================

def lambda_l_scipy(z, l):
    """λ_l(z) на [Z0,∞) через scipy ive/kve (амплитудно-стабильно при z→∞)."""
    nu = l + 0.5
    x = nu * z
    Ie, Ke = ive(nu, x), kve(nu, x)
    dIe = ive(nu + 1, x) + (nu / x - 1.0) * Ie
    dKe = (nu / x + 1.0) * Ke - kve(nu + 1, x)
    S = np.sqrt(np.pi * x / 2) * Ie
    E = np.sqrt(2 * x / np.pi) * Ke
    dS = S / (2 * x) + np.sqrt(np.pi * x / 2) * dIe
    dE = E / (2 * x) + np.sqrt(2 * x / np.pi) * dKe
    return E * dS + S * dE

def lambda_l_mp(z, l):
    """λ_l(z) на [0,Z0] через mpmath (малые z, x≤~0.8 для l=80)."""
    nu = l + mp.mpf('0.5')
    x = nu * z
    I = mp.besseli(nu, x); K = mp.besselk(nu, x)
    Ip = mp.besseli(nu + 1, x) + nu * I / x
    Kp = nu * K / x - mp.besselk(nu + 1, x)
    s = mp.sqrt(mp.pi * x / 2) * I
    e = mp.sqrt(2 * x / mp.pi) * K
    sp = s / (2 * x) + mp.sqrt(mp.pi * x / 2) * Ip
    ep = e / (2 * x) + mp.sqrt(2 * x / mp.pi) * Kp
    return sp * e + s * ep

def R_l(l):
    """R_l = −(1/2π)∫₀^∞ dz[(2l+1)²·ln(1−λ_l²) + (1+z²)⁻³]."""
    v0 = float(mp.quad(lambda z: (2*l+1)**2*mp.ln(1-lambda_l_mp(z,l)**2) + (1+z**2)**(-3),
                       [0, Z0], method='gauss-legendre'))
    f = lambda t: ((2*l+1)**2*np.log(1-lambda_l_scipy(np.tan(t),l)**2)
                   + (1+np.tan(t)**2)**(-3)) / np.cos(t)**2
    v1, _ = integrate.quad(f, np.arctan(Z0), pi/2 - 1e-9, limit=600,
                           epsabs=1e-12, epsrel=1e-12)
    return -(v0 + v1) / (2*pi)

def sphere_Ea(lmax=80):
    """E·a по Милтону (E = −(ΣR_l − 3/32)/2a), хвост R_l≈c/l²."""
    sums = 0.0
    for l in range(1, lmax + 1):
        sums += R_l(l)
    c_eff = R_l(lmax) * lmax * lmax
    tail = c_eff / lmax
    return -(sums + tail - 3/32) / 2.0, sums, tail

# ======================================================================
# 2. ТОР — модель wrapped-waveguide (структурная проверка Вейля)
# ======================================================================

def build_modes(a, Gmax):
    """(γ, mult) для PEC круглого волновода: γ=нули J_l/J'_l на (0,Gmax/a],
    вырождение 1 при l=0, 2 при l≥1 (обе поляризации учтены семействами)."""
    Lg = int(Gmax) + 2
    N  = int(Gmax) + 2
    out = []
    for l in range(0, Lg + 1):
        for fn in (jn_zeros, jnp_zeros):
            z = fn(l, N)
            z = z[z > 0]; z = z[z <= Gmax]
            m = 1 if l == 0 else 2
            out.extend((zz / a, m) for zz in z)
    out.sort(key=lambda t: t[0])
    return np.array(out, dtype=[("g", float), ("m", float)])

def E_tau(modes, a, R, tau, omax_units=42.0, chunk=2000000):
    """E(τ)=½Σ w·√((p/R)²+γ²)·e^{−τ·√(...)}, p∈[−pm,pm], ωmax=42/τ."""
    omax = omax_units / tau
    gs = modes["g"]; ws = modes["m"]
    d = omax**2 - gs**2
    sel = d > 0
    gs = gs[sel]; ws = ws[sel]; d = d[sel]
    pm = np.floor(np.sqrt(d) * R).astype(np.int64)
    sizes = 2*pm + 1
    tot = 0.0
    for i in range(0, len(gs), chunk):
        sl = slice(i, i+chunk)
        g = gs[sl]; w = ws[sl]; pmc = pm[sl]
        sc_ = 2*pmc + 1
        gr_ = np.repeat(g, sc_); wr_ = np.repeat(w, sc_)
        starts = np.empty(len(pmc), dtype=np.int64)
        np.cumsum(sc_[:-1], out=starts[1:]); starts[0] = 0
        p = np.arange(len(gr_)) - np.repeat(starts, sc_)
        om_ = np.sqrt((p / R)**2 + gr_**2)
        tot += float(np.sum(wr_ * om_ * np.exp(-tau * om_)))
    return 0.5 * tot

def torus_weyl_check(a=1e-3, ar=0.1, etas=(0.15, 0.20), Gmax=290):
    """Проверка лидирующего Вейлевского коэффициента: E(τ)·τ⁴ → 3a²R."""
    modes = build_modes(a, Gmax)
    R = a / ar
    rows = []
    for eta in etas:
        tau = eta * a
        E = E_tau(modes, a, R, tau)
        rows.append({"tau_over_a": eta, "E": E, "E_tau4": E * tau**4})
    theory = 3 * a*a * R
    return rows, theory

# ======================================================================
# 3. ЛИТЕРАТУРНЫЙ ЯКОРЬ ТОРА/ЦИЛИНДРА (DeRaad–Milton; SkK-2014 BD)
# ======================================================================
# C(R/a) — полная энергия Казимира тора с минорным радиусом 1
SK_TORUS = {2.0: -0.11, 3.0: -0.16, 4.0: -0.22, 5.0: -0.28,
            6.0: -0.33, 8.0: -0.44, 10.0: -0.56}
CYL_PER_LEN = -0.01356      # E/L = −0.01356/a² на единицу длины

def torus_C(Rovera):
    """интерполяция ln(R/a)→C в таблице SkK."""
    xs = np.log(np.array(sorted(SK_TORUS)))
    ys = np.array([SK_TORUS[x] for x in sorted(SK_TORUS)])
    return float(np.interp(np.log(Rovera), xs, ys))

def torus_rows(ar_list=(0.1, 0.3, 0.5), a_list=(1e-3, 1e-2, 1e-1)):
    rows = []
    for ar in ar_list:
        Rovera = 1.0 / ar
        C = torus_C(Rovera)
        for a in a_list:
            E_tot_J  = C * HBARC / a                 # полная, Дж
            R = a / Rovera
            E_L_J_m  = E_tot_J / (2*pi*R)            # на единицу длины
            cyl_tot  = CYL_PER_LEN * HBARC / a**2 * (2*pi*R)
            rho_tube = E_tot_J / (2*pi**2 * a*a * R)  # плотность в трубе, Дж/м³
            rows.append({
                "a_over_R": ar, "a_m": a,
                "C_lit": C,
                "E_tot_J": E_tot_J, "E_per_len_J_m": E_L_J_m,
                "E_cyl_est_J": cyl_tot,
                "rho_tube_J_m3": rho_tube,
                "rho_tube_kg_m3": rho_tube / c**2,
            })
    return rows

# ======================================================================
# 4. ПУЗЫРЬ (r4_bubound): ρ = −π²ℏc/(720·R_b⁴)
# ======================================================================

def bubble_rho(R_b):
    rho_J = -pi**2 * HBARC / (720.0 * R_b**4)
    return {"R_b_m": R_b, "rho_J_m3": rho_J, "rho_kg_m3": rho_J / c**2}

# ======================================================================
if __name__ == "__main__":
    results = {}

    print("=" * 110)
    print("  R4 [ГЕОМЕТРИЯ]  ЗАФИКСАЦИЯ ЗНАКА/ВЕЛИЧИНЫ КАЗИМИРА ПО ТОПОЛОГИИ;  ℏc=%.6e Дж·м" % HBARC)
    print("  Естеств. формулы в ħ=c=1; SI получен умножением на ℏc и делением на длину.")
    print("=" * 110)

    # --- 1. Сфера ---
    print("\n[1] СФЕРА PEC (Milton/Boyer):")
    Ea, sums, tail = sphere_Ea(lmax=80)
    ref = 0.0461765
    print("    E·a = %.8f   (ΣR_l=%.6f, хвост=%.2e)" % (Ea, sums, tail))
    print("    эталон 0.092353/2 = %.8f  → ratio = %.6f" % (ref, Ea / ref))
    results["sphere"] = {"Ea": Ea, "ref": ref, "ratio": Ea / ref,
                         "conclusion": "E>0, РЕПУЛЬСИЯ — сфера НЕ держит горловину"}

    # --- 2. Тор: структурная проверка Вейля ---
    print("\n[2] ТОР (wrapped-waveguide), a/R=0.1: Вейлевский коэффициент")
    wrows, theory = torus_weyl_check()
    for r in wrows:
        print("    t/a=%.2f  E·τ⁴ = %.6e  (Вейль 3a²R = %.6e, ratio %.4f)"
              % (r["tau_over_a"], r["E_tau4"], theory, r["E_tau4"] / theory))
    results["torus_wrapped_weyl"] = {
        "rows": wrows, "theory_3a2R": theory,
        "note": "лидирующий коэффициент ТОЧЕН (=3a²R); прямое извлечение константы "
                "численно не сходится (обусловленность ~1e5) → берём литературу [SkK-2014]"}

    # --- 3. Тор/цилиндр: литературный якорь ---
    print("\n[3] ТОР: полная энергия Казимира C (Balian–Duplantier, SkK-2014), минорный радиус 1")
    for Rv, C in sorted(SK_TORUS.items()):
        print("    R/a=%2d: C=%.2f" % (Rv, C))
    print("    Цилиндр (DeRaad–Milton): E/L = −0.01356/a² (единица длины)")
    rows = torus_rows()
    print("\n    Таблица для a∈{1мм,1см,10см} (SI, Дж):")
    print("      a_over_R    a_m        E_tot_J         E/L_J/m        rho_tube_J/m3")
    for r in rows:
        print("      %7.2f   %7.0e  %14.6e  %14.6e  %14.6e"
              % (r["a_over_R"], r["a_m"], r["E_tot_J"], r["E_per_len_J_m"], r["rho_tube_J_m3"]))
    results["torus"] = {"sk_table": SK_TORUS,
                        "cyl_per_len_a2": CYL_PER_LEN,
                        "rows": rows,
                        "conclusion": "E<0 для всех R/a; ρ_vac<0 в трубе — кандидат поддержки NEC"}

    # --- 4. Пузырь ---
    print("\n[4] ПУЗЫРЬ ρ = −π²ℏc/(720·R_b⁴)  (r4_bubound.py):")
    bubbles = []
    for Rb in (r0 := 0.10, 1.0, 10.0):
        b = bubble_rho(Rb)
        bubbles.append(b)
        print("    R_b=%.2f м: ρ = %.6e Дж/м³ = %.6e кг/м³" % (Rb, b["rho_J_m3"], b["rho_kg_m3"]))
    results["bubble"] = {"rows": bubbles,
                         "conclusion": "ρ<0 всегда — собственный источник отрицательной плотности"}

    # --- 5. Вывод ---
    print("\n[5] ВЫВОД (топология ↔ знак):")
    lines = [
        "  · сфера (замкнутая выпуклая, g=0): E>0 — РЕПУЛЬСИЯ, горловину НЕ держит;",
        "  · цилиндр / тор / компактификации: E<0, ρ_vac<0 — ПОДДЕРЖИВАЕТ условие ρ<0;",
        "  · пузырь Казимира: ρ<0 всегда;",
        "  · все модули ~ℏc/длина ≈ 1e-25..1e-23 Дж на единицу размера — макроскопически ничтожны",
        "    (ср. r4_bubound: достижимая гравитация только на наномасштабах).",
    ]
    for ln in lines:
        print(ln)
    results["conclusion"] = {
        "sign_by_topology": "sphere:+, cylinder/torus:−, bubble:−",
        "nec": "тор/цилиндр и пузырь дают ρ<0 (поддержка NEC); сфера даёт ρ>0 (нет)",
        "caveat": "величины ~ℏc-масштаба; тор взяты из BD-литературы (SkK-2014), "
                  "ибо прямое численное извлечение константы в wrapped-модели не сходится",
        "items": lines,
    }

    with open("/home/smboozha/portal_gun/research/r4_casimir_geo_values.json", "w") as f:
        json.dump(results, f, indent=1, ensure_ascii=False)
    print("\n[OK] r4_casimir_geo_values.json записан")