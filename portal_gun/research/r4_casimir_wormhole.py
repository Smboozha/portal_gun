"""
R4 [ОСН] — КОНСОЛИДИРОВАННЫЙ скан NEC/WEC «Пузыря Истинного Вакуума»
========================================================================
Свод 4 апгрейдов от Агентов 1–4:
  A1 r4_dce_source.py      — источник T_μν: статика Казимира (пластины/тор/пузырь/
                             сфера) + DCE-модуляция (ε, ω_d).
  A2 r4_qi_ft.py           — квантовые неравенства Форда–Романа в f(T): торсионный
                             вклад q=−22α меняет жёсткость источника,
                             τ_avail = τ_FR·(1−q)^{−1/4}; при α≤−1/22 QI не активны.
  A3 r4_casimir_geo.py     — топология: сфера даёт + (отталкивание), тор/цилиндр/
                             пузырь дают − (НЕК-поддерживающие).
  A4 r4_bubound.py         — граничные условия на краях пузыря (условия Израиля),
                             измеримость в лаборатории (a_g).

Вопросы:
  (1) В каком семействе f(T) и с каким источником Казимира выполняются NEC/WEC?
  (2) Жизнеспособность по QI-времени (с учётом кручения).
  (3) Знак энергии по топологии (нужен ρ<0).
  (4) Измерима ли горловина в лабе? (граница пузыря)

Единицы: 4π·r0²·X (X = {ρ,p_r,p_t}); r0 = 0.10 м; x = r/r0; α = T0·r0².
"""
import json
import mpmath as mp

import r4_dce_source   as A1
import r4_qi_ft        as A2
import r4_casimir_geo  as A3
import r4_bubound      as A4

r0 = 0.10
R0S = "0.10 m (10 см)"
THROAT_X = mp.mpf('1.0') + mp.mpf('1e-6')

def scan_quad_nec(alpha):
    """4πr0²·{ρ,p_r,p_t} на горловине для f=T+αT² (BHL (49)-(45), исправл. p_t)."""
    rho = -(1 + 22*alpha)/2
    pr  = -(1 + 2*alpha)/2
    pt  = (1 + 14*alpha)/2
    return rho, pr, pt

def main():
    print("="*88)
    print(" R4 [ОСН] — КОНСОЛИДИРОВАННЫЙ NEC/WEC «Пузырь Истинного Вакуума»")
    print("="*88)

    # ---------- 1) ГЕОМЕТРИЯ/ТОПОЛОГИЯ (А3): знак энергии ----------
    print("\n[1] ТОПОЛОГИЯ (А3) — знак энергии удерживающей горловину")
    print("    нужен ρ<0 вблизи горловины (NEC/флейр-аут); ρ>0 = отталкивание")
    print("    sphere  : E·a = +0.04617666 (Бойер/Милтон) → ρ>0  [НЕ подходит]")
    print("    torus   : C = -0.11..-0.56 (R/a=2..10, Balian-Duplantier) → ρ<0 [OK]")
    print("    cylinder: E/L = -0.01356/a^2 (DeRaad-Milton) → ρ<0 [OK]")
    print("    bubble  : ρ = -pi^2 hbar c /(720 R_b^4) → ρ<0 [OK]")
    print("    ВЫВОД: сфера отталкивает; тор/цилиндр/пузырь поддерживают ρ<0.")

    # ---------- 2) ИСТОЧНИК (А1): численные ρ для форм ----------
    print("\n[2] ИСТОЧНИК Казимира (А1) — 4πr0²·ρ при различных конфигурациях")
    rows = []
    for form, xarg, dm in [('plates', 1e-5, None), ('torus', 1e-3, None),
                           ('bubble', 1.0, 1e-6), ('sphere', 1.0, None)]:
        s = A1.source(form, xarg, dm)
        rows.append(dict(form=form, sign=s['sign']['rho'],
                         rho=s['rho'], pr=s['pr'], pt=s['pt'],
                         note=s['note'][:70], formula=s['formula']))
        print("    %-7s ρ=%.3e pr=%.3e pt=%.3e  [%s]" % (form, s['rho'], s['pr'],
                                                          s['pt'], s['sign']['rho']))
    print("    DCE: линейная модуляция ρ(t)=ρ_stat(1+ε sinω_d t) → ⟨ρ⟩=ρ_stat (0 вклада);")
    print("         адиабатика (1+ε sin)⁻⁴ даёт положит. фотонный член O(ε²): ⟨ρ⟩/ρ_stat=1.0513 (ε=0.1)")
    print("         реальный DCE требует 2f₁=c/d, т.е. d∈[3 см;300 м] — у нас зазоры мкм: режим квазистатический")

    # ---------- 3) NEC/WEC в f(T) (BHL) + QI-жизнеспособность (А2) ----------
    print("\n[3] NEC/WEC при горловине + QI-время жизни (А2, f=T+αT²)")
    print("    %6s | %7s%7s%7s%7s%7s | %8s | %s" %
          ("α", "ρ", "p_r", "p_t", "NECr", "NECt", "τ_life/τ_FR", "QI-верд."))
    qi_table = []
    for av in ['-5','-3','-2','-1.5','-1','-0.75','-0.5','-0.25','-0.1','0','0.1','0.25','0.5','1']:
        a = mp.mpf(av)
        rho, pr, pt = scan_quad_nec(a)
        nr, nt = rho+pr, rho+pt
        # QI по А2: τ_avail = τ_FR·(1-q)^{-1/4}, q=-22α на горловине (BHL)
        q_tor = -22*a
        tau_ratio = mp.power(mp.mpc(1 - q_tor), mp.mpf('-0.25'))
        nec_r = nr >= 0; nec_t = nt >= 0
        st = ("NECr%s NECt%s" % ("+" if nec_r else "-", "+" if nec_t else "-"))
        if q_tor >= 1:
            qi_v = "QI off (экзот. из геом.)"
        elif q_tor > 0:
            qi_v = "QI on, τ<τ_FR"
        else:
            qi_v = "QI on (α>0: хуже GR)"
        print("    %6s | %7.3f%7.3f%7.3f%7.3f%7.3f | %8.4f   | %s" %
              (av, rho, pr, pt, nr, nt, tau_ratio.real, qi_v))
        qi_table.append(dict(alpha=av, rho=float(rho), pr=float(pr), pt=float(pt),
                             nec_r=bool(nec_r), nec_t=bool(nec_t), tau_ratio=float(tau_ratio.real),
                             q_torsion=float(q_tor), qi_state=qi_v))

    print("    А2-вывод: QI Форда-Романа НЕ отменены; в GR τ~10⁻²⁷ с — GR-горловина")
    print("    нежизнеспособна. Торсионный вклад q=−22α сдвигает жёсткость: α≤−1/22")
    print("    (q≥1) → экзотика не требуется, QI неактивны; канон BHL α=−1 с запасом.")

    # ---------- 4) ГРАНИЦА ПУЗЫРЯ (А4): измеримость в лабе ----------
    print("\n[4] ГРАНИЦА ПУЗЫРЯ (А4, Израиль) — лабораторная измеримость")
    for Rb, d in [(0.1, 1e-6), (1.0, 1e-6), (10.0, 1e-6)]:
        rho_si = A4.casimir_rho(d)
        me = A4.m_eff_mag(Rb, rho_si)
        ag = A4.a_field(Rb, me)
        print("    R_b=%4.1f м d=%4.0f мкм: m_eff=%.2e кг  a_g=%.2e м/с²  (%s)" %
              (Rb, d*1e6, me, ag, "измеримо" if ag >= 1e-6 else "фантастика"))
    print("    ВЫВОД: во всех конфигурациях a_g << порога гравиметрии 1e-6 м/с²;")
    print("    для a_g=1e-6 нужен |ρ| ≈ плотность воды и зазор ~пм (нефизично).")

    # ---------- 5) КОНСОЛИДАЦИЯ ----------
    print("\n[5] КОНСОЛИДИРОВАННЫЙ ВЕРДИКТ")
    print("    ✓ NEC/WEC локально у горловины: f=T+αT² при α<0 (R2 + A1-поддержка)")
    print("    ✓ QI: торсионный сектор переносит часть жёсткости; α≤−1/22 → неактивны")
    print("    ✓ Топология: тор/пузырь/пластины дают ρ<0; сфера ρ>0 (не годится)")
    print("    ✗ Измеримость в лабе: граничные условия дают a_g ~ 10⁻³⁰..10⁻³⁹ м/с² —")
    print("      «Пузырь» на масштабах r0=10 см физически НЕДЕТЕКТИРУЕМ современными")
    print("      средствами. Реальный продукт — математически консистентная модель")
    print("      (f(T)-неэкзотичная горловина), НЕ лабораторный порт.")

    data = {
        "title": "R4 consolidated NEC/WEC scan, Casimir bubble in f(T)",
        "r0": R0S,
        "units": "4pi*r0^2 * X",
        "topology_sign": {
            "sphere": "+ (repulsive, NOT supporting)",
            "cylinder": "- (DeRaad-Milton 1981)",
            "torus": "- (Balian-Duplantier via A3; C=-0.11..-0.56)",
            "bubble": "- (plates-scale rho<0)"
        },
        "source_rho_4pir02": rows,
        "DCE": "linear mod averages to static; adiabatic adds +O(eps^2); real DCE needs cm-scale gaps",
        "nec_quad": qi_table,
        "QI": "not nullified; q=-22*alpha; alive for alpha<-1/22; GR unviable",
        "boundary": "all fantastic; a_g<<1e-6; need pm gaps",
        "verdict": ("NEC/WEC locally OK (alpha<0); QI dormant (torsion); "
                    "topology: torus/bubble support; BUT lab-detectable signal: NO")
    }
    out = '/home/smboozha/portal_gun/research/r4_consolidated.json'
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print("\n-> ", out)

if __name__ == '__main__':
    main()