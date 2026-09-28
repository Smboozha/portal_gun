"""
R4 [АГЕНТ-4, ВЫВОД/ГРАНИЦА] — границы «пузыря» в лабораторных координатах
==========================================================================
Пузырь: 3D-область радиуса R_b с отрицательной плотностью энергии Казимира
ρ<0 внутри; снаружи — плоский фон (Минковский). Эффективная отрицательная
гравитационная масса m_eff = |ρ|·(4π/3)·R_b³/c² ∈ {0.2e-31 .. 2e-17 кг}.

ДЖОНКЦИЯ ДАРМУА–ИЗРАИЛЯ (тонкая сферическая оболочка r=R_b, статика):
  S^a_b = −(c⁴/8πG)·([K^a_b] − δ^a_b [K]),  [X] ≡ X⁺ − X⁻   (внеш−внутр),
  σ = −S^τ_τ (поверхностная плотность энергии, Дж/м²),  P_θ = S^θ_θ.
  β_± = √(1 − r_s,±/R),  r_s = 2G·M_ext/c²  →  K^θ_θ = K^φ_φ = β/R,
  K^τ_τ = −r_s/(2R²β) (+O(β-степени, несущественно при r_s/R~1e-49)).
  Для плоского внутреннего патча (M_int=0, β_−=1): K^θ_θ=K^φ_φ=1/R, K^τ_τ=0.

  СЛАБОЕ ПОЛЕ (r_s/R ≪ 1), устойчивые формы без катастрофического 1−β:
    σ  = (c⁴/4πG)·r_s/(R²(1+β))      ≈ M_ext·c²/(4πR²)
    P_θ = −(c⁴/8πG)·(r_s/R²)·[1/(2β)+1/(1+β)]  ≈ −σ

ВЕТВИ ЗНАКОВ (честно):
  (A) M_ext = +m_eff: σ>0, P_θ≈−σ<0 — натяжение, обычная оболочка.
  (B) M_ext = −m_eff (отрицательная масса ⇒ отталкивание): σ<0, P_θ>0,
      оболочка экзотическая, a_g наружу. ← главная ветвь задачи.
  (C) M_ext = 0, снаружи СТРОГО Минковский: внешнего поля нет (a_g=0),
      оболочка β_+=β_−=1 ⇒ σ=0 по Джонкции, всё уходит в скачок ∂_t: 
      внутренняя отрицательная энергия полностью экранируется (шлем-условие);
      в «объёмной» реализации σ=|ρ|·c²R/3 балансирует M_ADM=0.
  Для сетки главная таблица даёт ветвь (B); для сравнения показана и σ(A).

ЧИСЛА НЕ ТОЛЬКО «МАЛЫЕ» — НУЖНЫ УСТОЙЧИВЫЕ: r_s/R = 3e-49..3e-45, прямое
1−β коллапсирует в double (в прошлой версии σ печаталась как −0.0).

МАСШТАБ КРИВИЗНЫ ВНУТРЕННЕГО ПАТЧА:
  · реализация «однородный ρ<0» (уравнения Эйнштейна, AdS):  ℓ_ρ=√(3c²/(8πG|ρ|))
  · вложенная кривизна границы S²: k=2/R_b, характерная длина R_b/2.

ОГОВОРКИ (честно):
  1) Энергия Казимира занимает ЗАЗОР d между пластинами на оболочке, а не весь
     шар: в оболочковой реализации m_eff умножить на 3d/R_b ≪ 1 (ещё меньше).
  2) При M_ext=0 снаружи нечего измерять гравиметрией (a_g=0) — детекция только
     по локальным свойствам оболочки σ и внутренней геометрии.
  3) Масса отрицательна ⇒ ускорение пробного тела направлено НАРУЖУ.
"""
import json
import numpy as np

# ---- константы (SI) ----
hbar   = 1.054571817e-34
c      = 2.99792458e8
G      = 6.67430e-11
pi     = np.pi
g_earth = 9.80665

C_SI   = (pi**2 * hbar * c) / 720.0        # ε = −C_SI/d⁴  [Дж/м³]

R_b_list = [0.1, 1.0, 10.0]
d_list   = [1e-6, 10e-6, 100e-6]

a_thr_grav = 1e-6      # м/с² — гравиметрия достижима
a_thr_atom = 1e-9      # м/с² — атомный интерферометр


def casimir_rho(d):
    """массовая плотность ρ<0 [кг/м³] от энергии Казимира в зазоре d."""
    return -C_SI / d**4 / c**2


def m_eff_mag(R_b, rho):
    """|m_eff| = |ρ|·(4π/3)R_b³  [кг]."""
    return -rho * (4.0 / 3.0) * pi * R_b**3


def shell_junction(R_b, m_ext):
    """Дармуа–Израиль для тонкой сферы; m_ext — внешняя масса со знаком.
    Возвращает (r_s, β, σ, P_θ), σ и P_θ в Дж/м², устойчивые к 1−β."""
    rs = 2.0 * G * m_ext / c**2             # со знаком; <0 для отрицательной массы
    beta = np.sqrt(1.0 - rs / R_b)
    den = R_b**2 * (1.0 + beta)
    sigma = (c**4 / (4.0 * pi * G)) * rs / den
    P_theta = -(c**4 / (8.0 * pi * G)) * (rs / R_b**2) * (1.0 / (2.0 * beta) + 1.0 / (1.0 + beta))
    return rs, beta, sigma, P_theta


def a_field(R_b, m_eff_mag_):
    """величина пробного ускорения на границе пузыря; направление наружу, если масса<0."""
    return G * m_eff_mag_ / R_b**2


def adS_radius(rho):
    """длина кривизны внутреннего патча при реализации ρ=const (AdS-аналог)."""
    return np.sqrt(3.0 * c**2 / (8.0 * pi * G * abs(rho)))


def need_rho(R_b, a=1.0):
    """|ρ|_need для a_g = a:  a = G·|ρ|·4πR³/3/R² = 4πG|ρ|R/3 ⇒ ρ = 3a/(4πGR)."""
    return 3.0 * a / (4.0 * pi * G * R_b)


def d_need(R_b, a=1.0):
    """зазор d, дающий требуемую |ρ|: C_SI/d⁴/c² = ρ_need."""
    return (C_SI / (need_rho(R_b, a) * c**2))**0.25


def verdict(a):
    if a > a_thr_grav:
        return 'реально (гравиметрия > 1e-6 м/с²)'
    if a > a_thr_atom:
        return 'погранично (атомный интерферометр 1e-9..1e-6)'
    return 'фантастика (< 1e-9 м/с²)'


# ----------------------------------------------------------------------
rows = []
for R_b in R_b_list:
    for d in d_list:
        rho = casimir_rho(d)
        me = m_eff_mag(R_b, rho)
        # ветвь B (отрицательная масса, главная) и ветвь A (положительная, справка)
        rsB, betaB, sigB, PB = shell_junction(R_b, -me)
        rsA, betaA, sigA, PA = shell_junction(R_b, +me)
        ag = a_field(R_b, me)
        rows.append({
            'R_b_m': R_b, 'd_m': d,
            'rho_kg_m3': rho,
            'm_eff_kg': me,
            'r_s_m': rsB, 'rs_over_Rb': abs(rsB) / R_b,
            'sigma_B_J_m2': sigB,              # ветвь B: отрицательная масса
            'sigma_A_J_m2': sigA,              # ветвь A: положительная масса
            'P_theta_B_J_m2': PB, 'P_over_sigma_B': PB / sigB,
            'R_curv_adS_m': adS_radius(rho),   # масштаб кривизны внутреннего патча
            'R_embed_boundary_m': R_b / 2.0,
            'a_g_m_s2': ag,
            'a_g_over_g': ag / g_earth,
            'verdict': verdict(ag),
        })

print("=" * 108)
print("  R4 [ВЫВОД/ГРАНИЦА]  ПУЗЫРЬ С ρ<0: ИЗРАИЛЬ-ШОВ, СЕТКА R_b×d, ЛОКАЛЬНАЯ ГРАВИМЕТРИЯ")
print("  ε_Каз = −π²ℏc/(720d⁴);  |m_eff| = |ρ|·4πR_b³/3;  a_g = Gm_eff/R_b²")
print("  σ = (c⁴/4πG)·r_s/(R²(1+β)),  P_θ = −(c⁴/8πG)(r_s/R²)[1/(2β)+1/(1+β)],  β=√(1−r_s/R)")
print("  ветвь B: масса −m_eff (отталкивание, σ<0) ← главная;  ветвь A: +m_eff (σ>0) — справка")
print("=" * 108)

hdr = (("R_b/d".ljust(12))
       + ("|m_eff|[кг]".rjust(13))
       + ("r_s/R".rjust(11))
       + ("σ_B[Дж/м²]".rjust(12))
       + ("σ_A[Дж/м²]".rjust(12))
       + ("P_θB[Дж/м²]".rjust(12))
       + ("a_g[м/с²]".rjust(12))
       + ("a_g/g".rjust(10))
       + "  статус")
print(hdr)
print("-" * 108)
for r in rows:
    print(("%5.2fм/%3.0fмкм".ljust(12)) % (r['R_b_m'], r['d_m'] * 1e6)
          + ("%13.3e" % r['m_eff_kg'])
          + ("%11.2e" % r['rs_over_Rb'])
          + ("%12.3e" % r['sigma_B_J_m2'])
          + ("%12.3e" % r['sigma_A_J_m2'])
          + ("%12.3e" % r['P_theta_B_J_m2'])
          + ("%12.3e" % r['a_g_m_s2'])
          + ("%10.2e" % r['a_g_over_g'])
          + "  " + r['verdict'])

print("-" * 108)
print("  ρ(Казимир) и масштаб кривизны внутреннего патча (AdS-реализация ρ=const, ℓ_ρ):")
for d in d_list:
    rho = casimir_rho(d)
    print("    d=%5.0f мкм:  ρ=%10.3e кг/м³  ε=%10.3e Дж/м³  ℓ_ρ=%10.3e м (~%9.3f кпк)"
          % (d * 1e6, rho, rho * c**2, adS_radius(rho), adS_radius(rho) / 3.086e19))

print("-" * 108)
print("  ОГОВОРКА 1 (оболочковая реализация): энергия Казимира живёт в зазоре d,")
print("  тогда объём отрицательной энергии ≈ 4πR_b²d, а не 4πR_b³/3:")
for R_b in R_b_list:
    for d in d_list:
        f = 3.0 * d / R_b
        r = next(x for x in rows if x['R_b_m'] == R_b and x['d_m'] == d)
        print("    R_b=%5.2f м, d=%4.0f мкм: m_eff(шар)=%9.2e кг  →  m_eff(оболочка)=%9.2e кг (×%9.2e)"
              % (R_b, d * 1e6, r['m_eff_kg'], r['m_eff_kg'] * f, f))

print("=" * 108)
print("  ДО ПОРОГА a_g>1e-6 м/с² (|ρ|_need = 3·a_thr/(4πGR_b);  d_need = (C_SI/(ρ_need c²))^{1/4}):")
print("  %-16s %-14s %-14s %-14s %-14s" % ("R_b / d", "|ρ|_need кг/м³", "d_need м", "нехватка ×|ρ|", "d_need/нм"))
for r in rows:
    nrho = need_rho(r['R_b_m'], a_thr_grav)
    dn = d_need(r['R_b_m'], a_thr_grav)
    gap = nrho / abs(r['rho_kg_m3'])
    print("  %5.2f м / %4.0f мкм  %-14.3e %-14.3e %-14.3e %-14.2f"
          % (r['R_b_m'], r['d_m'] * 1e6, nrho, dn, gap, dn / 1e-9))
print("  (НЕХВАТКА ×1e22..×1e32 даже на самой плотной ячейке; d_need = 0.6..2 пм, нефизично)")

# ----------------------------------------------------------------------
# Вердикт и рекомендация
# ----------------------------------------------------------------------
best = min(rows, key=lambda r: abs(r['a_g_m_s2'] - a_thr_grav))
closest = max(rows, key=lambda r: r['a_g_m_s2'])
recommendation = {
    'measurable_config': 'НЕТ в сетке R_b=0.1..10 м × d=1..100 мкм',
    'max_a_g_m_s2': closest['a_g_m_s2'],
    'closest_config': {'R_b_m': closest['R_b_m'], 'd_m': closest['d_m']},
    'verdict_all': 'фантастика',
    'reason': ('a_g = 1.35e-40..1.35e-29 м/с² на сетке, т.е. в 1e23..1e33 раз ниже '
               'порога гравиметрии 1e-6 м/с²; даже |ρ|_need=358..3.58e4 кг/м³ (плотность '
               'воды и выше!) требует зазор d_need~0.6..2 пм (атомный/ядерный масштаб), '
               'где теория Казимира и сама идея пластин не работают'),
    'gap_rho': '×1e22..×1e32',
    'note_negative_mass': ('масса отрицательна ⇒ ускорение пробного тела направлено наружу '
                           '(отталкивание, антигравитация); скаляр |a_g| приведён в таблице'),
    'note_exterior_flat': ('при строго Minkowski-снаружи (M_ext=0) внешнее поле обнуляется, '
                           'a_g=0 снаружи — ничего не измерить гравиметрией'),
}

out = {
    'dataset': 'r4_bubound_values',
    'constants': {'hbar': hbar, 'c': c, 'G': G, 'g_earth': g_earth,
                  'C_SI_J_m3_m4': C_SI},
    'formulas': {
        'casimir_eps': '-pi^2*hbar*c/(720*d^4)  [J/m3]',
        'casimir_rho': 'eps/c^2  [kg/m3]',
        'm_eff': 'abs(rho)*(4*pi/3)*R_b^3',
        'r_s': '2*G*m_ext/c^2  (sign = sign of mass)',
        'sigma_israel_stable': '(c^4/(4*pi*G))*r_s/(R^2*(1+beta))  ~ m_ext*c^2/(4*pi*R^2)',
        'P_theta_israel_stable': '-(c^4/(8*pi*G))*(r_s/R^2)*[1/(2*beta)+1/(1+beta)] ~ -sigma',
        'beta': 'sqrt(1 - r_s/R)',
        'a_g': 'G*|m_eff|/R^2  (outward for negative mass)',
        'needed_rho': '3*a_thr/(4*pi*G*R_b)',
        'd_needed': '(C_SI/(rho_need*c^2))^(1/4)',
        'adS_radius_inner': 'sqrt(3*c^2/(8*pi*G*|rho|))',
    },
    'thresholds': {'gravimetry_m_s2': a_thr_grav, 'atom_interferometer_m_s2': a_thr_atom},
    'branches': {
        'A_positive_mass': {'sigma': '>0', 'P_theta': '<0 (tension)', 'probe': 'attractive'},
        'B_negative_mass': {'sigma': '<0 exotic shell', 'P_theta': '>0', 'probe': 'repulsive', 'used_as_primary': True},
        'C_exterior_flat_M0': {'a_g': '=0 outside', 'sigma_volume_estimate': 'abs(rho)*c^2*R/3 balancing', 'note': 'nothing to measure by gravimetry'},
    },
    'assumptions': {
        'ideal_parallel_plates': True,
        'D=3_casimir': True,
        'interior_patch': 'flat Minkowski (R_scalar=0); AdS-curvature scale listed separately',
        'exterior_patch': 'Schwarzschild with mass = sign(m_eff) regime, beta=sqrt(1-r_s/R)',
        'shell_israel': 'thin spherical shell, static, Darmois-Israel',
        'geometric_caveat': 'Casimir energy fills gap d (shell realization): rescale m_eff by 3d/R_b',
        'sign_branch': 'B (negative mass, repulsive probe) as primary, matching R4 task statement',
    },
    'values': rows,
    'recommendation': recommendation,
}

with open('/home/smboozha/portal_gun/research/r4_bubound_values.json', 'w',
          encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

print("=" * 108)
print("  ИТОГ: детектируемых (a_g≥1e-6 м/с²) конфигураций на сетке НЕТ.")
print("  Лучшая ячейка: R_b=%g м, d=%g мкм → a_g=%g м/с² (в %g раз ниже порога)."
      % (closest['R_b_m'], closest['d_m'] * 1e6, closest['a_g_m_s2'],
         a_thr_grav / max(closest['a_g_m_s2'], 1e-300)))
print("  Вывод: «фантастика» на всей сетке; рекомендации по расширению порога см. JSON.")
print("  JSON: research/r4_bubound_values.json сохранён")
print("=" * 108)