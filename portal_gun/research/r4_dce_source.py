"""
R4 [АГЕНТ-1 КОД] — источники T_μν: статический Казимир + DCE (динамический Казимир)
===================================================================================
Формат как в ft2_necscan_BHL.py: функции (rho, pr, pt) в безразмерных 4π·r0²
единицах, r0 = 10 см. Конверсия физической энергии/давления [Дж/м³] в «4πr0²·X»:
  X = 4π r0² · κ · ε_физ ,  κ = 8πG/c⁴  (это источник в полевых уравнениях
  в тех же единицах, что кривизна ~ 1/r0²; пр-но: 4πr0²·G_00 ≈ 1 при r=r0).

ФИЗИКА (честно):
  1) СТАТИЧЕСКИЙ КАЗИМИР, идеальные пластины, зазор d:
       ε(d) = −π²ħc/(720 d⁴)  [Дж/м³]   (ρ однородно между пластинами)
       p_r = 3ρ  (натяжение, нормаль к пластинам, притяжение)
       p_t = −ρ  (тангенциальное давление)
     конформность EM: след T = ρ − p_r − 2p_t = 0.
     Сила на пластину: |p_r| = π²ħc/(240 d⁴) (классика ~1.3 мПа при 1 мкм).
     NEC-скан: ρ<0 (WEC нарушен); ρ+p_r=4ρ<0 нарушен радиально;
               ρ+p_t=0 — тангенц. NEC насыщен. (Казимировская анизотропная EoS.)
  2) DCE: пластина осциллирует d(t)=d(1+ε sin ω_d t), ε∈[0,1], f_d∈[1 МГц..10 ГГц].
     Мгновенно-адиабатическая модуляция (линейная по задаче): ρ(t)=ρ_stat(1+ε sinω_dt);
     честная адиабатика: ρ(t)=ρ_stat(1+ε sinω_dt)⁻⁴.
     Средняя по периоду: РОВНО 0 вклад на O(ε) (⟨sin⟩=0) → ⟨ρ⟩=ρ_stat<0.
     Настоящие фотоны DCE — положительная энергия, вклад O(ε²): квазистатически
     ⟨ρ⟩=ρ_stat·⟨(1+εsin)⁻⁴⟩ = ρ_stat(1+5ε²+105/8·ε⁴+…)>ρ_stat (менее отрицательно).
     Реальный парный DCE запускается только при резонансе f_d≈2f_n = nc/d;
     для наших зазоров (10 нм..10 мкм) 2f₁=c/d ≫ 10 ГГц → в полосе нет резонанса,
     режим адиабатический. Честные пределы скана — внизу отчёта.
ФОРМЫ (x = безразмерный параметр = размер/r0):
   'plates'  x=d/r0      ε=−π²ħc/720d⁴    знак −   (анизотропно p_r=3ρ, p_t=−ρ)
   'sphere'  x=R/r0      ε=+0.0923531ħc/(2R·4πR³/3)  знак +  (Boyer/Milton, отталкивание)
   'torus'   x=L/r0      ε=−π²ħc/45L⁴     знак −   (топологический: PBC, 1 компактное
                          направление L; EM 2 поляризации; = 16× пластин)
   'bubble'  x=R_b/r0    внутри ε=−π²ħc/720d⁴  знак − (лабораторный пузырь r4_bubound,
                          зазор d; E=ε·4π/3·R_b³)
ПРОВЕРКА коэффициентов ДВУМЯ НЕЗАВИСИМЫМИ ПУТЯМИ (прямая формула + scipy.integrate):
   прямой:  C = π²/720 (пластины), π²/45 (torus) [mpmath];
   интеграл: ζ(−3)=1/120 по Абелю–Плане 2∫₀∞ t³/(e^{2πt}−1)dt [scipy.integrate.quad],
             далее C = (π²/12)·(множитель)·2·ζ(−3). Расхождение ≤ 1e-6.
   Сфера: коэффициент 0.0923531 — ЛИТЕРАТУРНЫЙ (Boyer 1968; Balian–Duplantier 1978;
          Milton 1980; Leseduarte–Romeo 1996), самостоятельный пересчёт не проводился —
          помечено честно в JSON.
"""
import json
import numpy as np
import mpmath as mp
from scipy.integrate import quad

mp.mp.dps = 30

# ----------------------------------------------------------------- constants
r0    = 0.10                      # м   (масштаб червоточины, как в ft2_necscan_BHL)
hbar  = 1.054571817e-34           # Дж·с
c     = 2.99792458e8              # м/с
G     = 6.67430e-11               # м³/(кг·с²)
C_SI  = float(mp.pi**2 * hbar * c / 720)      # ε = −C_SI/d⁴,  C_SI [Дж·м]

kappa  = float(8*mp.pi*G/c**4)                  # 1/(кг·м) → κ·ε [1/м²]
DIMFAC = float(4*mp.pi*r0**2*kappa)             # [Дж/м³] → «4πr0²·X»

D_GAPS = [10e-9, 1e-6, 10e-6]                   # 10 нм, 1 мкм, 10 мкм
F_BAND = (1e6, 1e10)                            # [1 МГц .. 10 ГГц]
R_SHELL_COEF = 0.0923531                        # E_sph = +coef·ħc/(2R)  (Boyer/Milton)

# ---------------------------------------------------------------- static plates
def casimir_plates(d):
    """Идеальные пластины, зазор d [м]: ρ[Dж/м³], p_r, p_t (физ., анизотропн.)."""
    eps = -C_SI/d**4
    return dict(rho=eps, pr=3.0*eps, pt=-eps, F_A=-3.0*eps)

def to_dimless(si):
    """физическое значение [Дж/м³] → безразмерное 4πr0²·X."""
    return DIMFAC*si

# ---------------------- два независимых способа для коэффициентов --------------
def zeta_minus3_abelplana():
    """ζ(−3) = 1/120 через Абеля–Плану:  2∫₀∞ t³/(e^{2πt}−1) dt  (scipy.integrate)."""
    val, err = quad(lambda t: t**3/(np.exp(2*np.pi*t) - 1.0), 0.0, 60.0, limit=400)
    return 2.0*val

def slab_C(form):
    """C в ε = −C·ħc/размер⁴ для 'plates' (размер=d) и 'torus' (размер=L), EM 2 поляриз."""
    z = zeta_minus3_abelplana()                 # ≈ ζ(−3) = 1/120
    mult = 16.0 if form == 'torus' else 1.0     # torus: k→2πn/L ⇒ ×2³, ±n ⇒ ×2  => (2⁴)
    return float(mp.pi**2/12.0*2.0*mult*z)      # ×2 — две поляризации EM

C_PL_DIRECT = float(mp.pi**2/720.0)             # ε_plates = −C ħc/d⁴
C_TO_DIRECT = float(mp.pi**2/45.0)              # ε_torus  = −C ħc/L⁴
C_PL_NUM    = slab_C('plates')
C_TO_NUM    = slab_C('torus')

def _check(name, direct, num):
    rel = abs(direct-num)/abs(direct)
    ok  = rel <= 1e-6
    print("    %-22s direct=%.12e  num(quad)=%.12e  |Δ|/dir=%.2e  %s"
          % (name, direct, num, rel, "OK(≤1e-6)" if ok else "FAIL"))
    return rel, ok

# ------------------------------------------------------------ четыре формы
def source(form, x, d_m=None):
    """(rho, pr, pt) в 4πr0²-безразмерности в формате ft2_necscan_BHL.

    form : 'plates'|'sphere'|'torus'|'bubble'
    x    : d/r0 (plates), R/r0 (sphere), L/r0 (torus), R_b/r0 (bubble)
    d_m  : для bubble — зазор Казимира в метрах.
    Возвращает словарь с SI и безразмерными значениями + знак и комментарий.
    """
    if form == 'plates':
        d = x*r0
        s = casimir_plates(d)
        sign = {'rho': '-', 'pr': '-', 'pt': '+'}
        note = ('идеальные плоские пластины; ρ<0, натяжение по нормали p_r=3ρ<0; '
                'тангенциально p_t=−ρ>0')
        formula = 'ε = −π²ħc/(720 d⁴)  (d = x·r0)'
    elif form == 'sphere':
        R = x*r0
        eps = R_SHELL_COEF*hbar*c/(2.0*R*(4.0/3.0)*mp.pi*R**3)   # E/V, сфера-оболочка
        eps = float(eps)
        s = dict(rho=eps, pr=eps/3.0, pt=eps/3.0, F_A=0.0)
        sign = {'rho': '+', 'pr': '+', 'pt': '+'}
        note = ('проводящая сферическая оболочка (Boyer 1968, Milton 1980): '
                'САМОЭНЕРГИЯ ОТТАЛКИВАТЕЛЬНА, E=+0.0923531ħc/2R; усреднённая по объёму '
                'плотность >0 — НЕ даёт отрицательной энергии (в отличие от пластин)')
        formula = 'ε = +0.0923531·ħc/(2R·(4π/3)R³)  (R = x·r0)'
    elif form == 'torus':
        L = x*r0
        eps = -float(C_TO_DIRECT)*hbar*c/L**4
        s = dict(rho=eps, pr=eps/3.0, pt=eps/3.0, F_A=0.0)
        sign = {'rho': '-', 'pr': '-', 'pt': '-'}
        note = ('топологический тор: ОДНО компактное (периодическое) направление длины L, '
                'EM 2 поляризации. Результат ζ(−3): 16× пластин (знак −). '
                'Физическая кривизна тора не учитывается (плоский тор/калибровочное поле); '
                'если нужен СКРУЧЕННЫЙ тор — знак зависит от модуляции, честно предупр.')
        formula = 'ε = −π²ħc/(45 L⁴)  (L = x·r0, PBC, = 16·пластин)'
    elif form == 'bubble':
        R_b = x*r0
        if d_m is None:
            raise ValueError("bubble: задай d_m (зазор Казимира в метрах)")
        eps = -C_SI/d_m**4
        E_tot = eps*(4.0/3.0)*mp.pi*R_b**3
        s = dict(rho=eps, pr=eps/3.0, pt=eps/3.0,
                 E_tot=float(E_tot), m_eff_kg=float(-E_tot/c**2))
        sign = {'rho': '-', 'pr': '-', 'pt': '-'}
        note = ('лабораторный «пузырь» r4_bubound: объём 4π/3·R_b³ заполнен '
                'Казимировской энергией пластин (зазор d_m) → ρ<0, изотропное '
                'натяжение p=ε/3<0; E=ε·V ')
        formula = 'ε = −π²ħc/(720 d_m⁴), E = ε·(4π/3)R_b³  (R_b = x·r0)'
    else:
        raise ValueError(form)
    return dict(form=form, x=x, d_m=d_m,
                rho_SI=s['rho'], pr_SI=s['pr'], pt_SI=s['pt'],
                rho=to_dimless(s['rho']), pr=to_dimless(s['pr']), pt=to_dimless(s['pt']),
                sign=sign, formula=formula, note=note,
                E_tot_J=s.get('E_tot'), m_eff_kg=s.get('m_eff_kg'))

# ------------------------------------------------------------ DCE
def dce_ratio(t, eps, linear=True):
    """ρ(t)/ρ_stat. linear: 1+ε sin ωt (задача); exact: (1+ε sin ωt)⁻⁴ (честная адиаб.)."""
    s = np.sin(t)
    return 1.0 + eps*s if linear else (1.0 + eps*s)**(-4)

def dce_average(eps, linear=True):
    """⟨ρ⟩/ρ_stat за период (quad по ωt ∈ [0,2π])."""
    val, err = quad(lambda t: dce_ratio(t, eps, linear), 0.0, 2.0*np.pi, limit=400)
    return val/(2.0*np.pi)

def dce_series(eps):
    """mpmath-ряд ⟨(1+εsin)⁻⁴⟩ = 1 + 5ε² + 105/8·ε⁴ + (бином. coeff)·⟨sin^{2k}⟩."""
    t = mp.mpf(0)
    for k in range(1, 8):
        coef = sum(mp.binomial(-4, j)*eps**j for j in [2*k])   # член j=2k
        mean_sin = mp.binomial(2*k, k)/mp.power(4, k)
        t += coef*mean_sin
    return float(1.0 + t)

def dce_scan():
    """Честный скан: 2f₁=c/d (резонанс) против полосы [1 МГц..10 ГГц]."""
    rows = []
    for d in D_GAPS:
        f1  = c/(2*d)          # фундаментальная мода полости
        f_res = c/d            # вырожденный резонанс f_d = 2f₁
        in_band = F_BAND[0] <= f_res <= F_BAND[1]
        detune = F_BAND[1]/f_res - 1.0          # при отстройке на верхней границе полосы
        rows.append(dict(d_m=d,
                         f1_Hz=f1, f_res_Hz=f_res,
                         in_band=in_band,
                         detune_at_top=detune,
                         kappa_ad_s1=None,
                         verdict=('резонанс вне полосы: адиабатический режим, '
                                  'парный DCE подавлен') if not in_band
                                 else 'В ПОЛОСЕ: возможен параметрический резонанс'))
    d_demo = float(c/F_BAND[1])                  # 3 см: резонанс на верхней границе полосы
    for eps in (0.01, 0.1, 0.3):
        w1 = np.pi*c/d_demo
        kap = eps*w1/4.0                          # темп роста амплитуды, Law/Dodonov
        T = 1.0/F_BAND[1]
        for nper in (1, 100):
            # κ·t ≪ 1: N ≈ (κt)²; общий закон N = sinh²(κt) — только до декогерентности
            N = np.sinh(kap*nper*T)**2
            rows.append(dict(d_m=d_demo, demo=True, eps=eps,
                             kappa_s1=kap, nper=nper, N_photons=N))
    return rows

def main():
    print("="*94)
    print("  R4 [АГЕНТ-1 КОД] — ИСТОЧНИКИ T_μν: СТАТИЧЕСКИЙ КАЗИМИР + DCE")
    print("  Формат 4πr0²·{rho,pr,pt} как в ft2_necscan_BHL;  r0=0.10 м;  κ=8πG/c⁴")
    print("="*94)

    print("\n 1) ПРОВЕРКА коэффициентов двумя независимыми способами (≤1e-6):")
    r1, o1 = _check("C_plates  ε=−π²ħc/720d⁴", C_PL_DIRECT, C_PL_NUM)
    r2, o2 = _check("C_torus   ε=−π²ħc/45L⁴ ", C_TO_DIRECT, C_TO_NUM)
    z = zeta_minus3_abelplana()
    print("    ζ(−3) numer(Abel-Plana)=%.12e  exact=1/120=%.12e  |Δ|=%.2e"
          % (z, 1/120.0, abs(z-1/120.0)))

    print("\n 2) СТАТИЧЕСКИЙ КАЗИМИР (пластины): ρ, p_r=3ρ, p_t=−ρ при трёх зазорах")
    print("    %8s|%12s%12s%12s|%11s%11s%11s|%10s| %s"
          % ("d", "ε[Дж/м³]", "ρm=ε/c²[кг/м³]", "|p_r|[Па]",
             "ρ·4πr0²", "p_r·4πr0²", "p_t·4πr0²", "|F|Па", "NEC"))
    tab = []
    for d in D_GAPS:
        s = casimir_plates(d)
        dm = to_dimless(s['rho'])
        nec_r = s['rho']+s['pr']; nec_t = s['rho']+s['pt']
        tab.append(dict(d_m=d, d_mm=d*1e6,
                        eps_Jm3=s['rho'], rho_kg_m3=s['rho']/c**2,
                        F_A_Pa=s['F_A'], rho_dm=dm,
                        nec_r_Pa=nec_r, nec_t_Pa=nec_t))
        print("    %4.0f nm %9.3e %12.3e %12.3e |%12.3e%12.3e%12.3e|%10.3e| %s"
              % (d*1e9, s['rho'], s['rho']/c**2, abs(s['pr']),
                 dm, dm*3.0, -dm, s['F_A'],
                 "ρ<0 WEC−; ρ+pr<0 NEC−; ρ+pt=0" ))

    print("\n 3) ЧЕТЫРЕ ФОРМЫ (знак, формула, 4πr0²-значения):")
    forms_out = []
    for form, x, dm in [('plates', 1e-6/r0, None),
                        ('sphere', 0.1/r0, None),
                        ('torus',  1e-6/r0, None),
                        ('bubble', 1.0/r0, 1e-6)]:
        r = source(form, x, d_m=dm)
        forms_out.append(r)
        print("    %-8s x=%.3e  ρsign=%s  ρ(4πr0²)=%+.3e  pr=%+.3e  pt=%+.3e"
              % (form, r['x'], r['sign']['rho'], r['rho'], r['pr'], r['pt']))
        print("      %s" % r['formula'])
        print("      %s" % r['note'])

    print("\n 4) DCE — модуляция и честный отчёт по периоду/резонансу:")
    for eps in (0.01, 0.1, 0.3):
        a_lin = dce_average(eps, linear=True)
        a_ex  = dce_average(eps, linear=False)
        val_series = dce_series(eps)
        rel = abs(a_ex-val_series)/a_ex
        print("    ε=%.2f: ⟨ρ⟩_lin/ρ_stat=%.6e (O(ε)-вклад=0); ⟨ρ⟩_exact/ρ_stat=%.6f"
              " (Δ=%+.4e, фотоны O(ε²)); ряд=%.6f |Δряд|=%.2e %s"
              % (eps, a_lin-1.0, a_ex, a_ex-1.0, val_series, rel,
                 "OK(≤1e-6)" if rel <= 1e-6 else "FAIL"))
    dce = dce_scan()
    for r in dce:
        if r.get('demo'):
            print("    [демо d=%.3gm ε=%.2f] κ=%+.3e 1/с;  N(%d пер.)=%.3e"
                  % (r['d_m'], r['eps'], r['kappa_s1'], r['nper'], r['N_photons']))
        else:
            print("    d=%-9.2e  f₁=%+.3e ГГц  f_res=c/d=%+.3e ГГц  в полосе=%s"
                  "  отстройка@10ГГц=%+.3e  → %s"
                  % (r['d_m'], r['f1_Hz']/1e9, r['f_res_Hz']/1e9,
                     r['in_band'], r['detune_at_top'], r['verdict']))
    print("    ПРЕДЕЛЫ: резонанс 2f₁=c/d попадает в [1 МГц;10 ГГц] ТОЛЬКО при "
          "d∈[%.3g;%.3g] м — наши зазоры (10 нм..10 мкм) далеко внизу ⇒"
          % (c/F_BAND[1], c/F_BAND[0]))
    print("           режим квазистатический; реальные фотоны DCE невозможны в полосе.")
    print("           Честно: ⟨ρ(t)⟩ за период = ρ_stat + 0 (O(ε)) ; фотонная добавка "
          "положительна, O(ε²);")
    print("           зеркало теряет механическую энергию НА ЭТУ ВЕЛИЧИНУ (закон сохранения).")

    zek = (abs(C_PL_DIRECT-C_PL_NUM)/C_PL_DIRECT, abs(C_TO_DIRECT-C_TO_NUM)/C_TO_DIRECT)
    out = {
        "dataset": "r4_dce_values",
        "author": "АГЕНТ-1 КОД",
        "units": "4pi*r0^2*X как в ft2_necscan_BHL; r0=0.10 м; "
                 "физ. ε[Дж/м³]→безразмер.=4πr0²·(8πG/c⁴)·ε",
        "constants": {"hbar": hbar, "c": c, "G": G, "r0": r0,
                      "C_SI_Jm": C_SI, "kappa": kappa, "dimfac": DIMFAC},
        "static_casimir_plates": {
            "formula": "ε=−π²ħc/(720d⁴); p_r=3ε; p_t=−ε; |F_A|=π²ħc/(240d⁴)",
            "values": [dict(d_m=t['d_m'], eps_Jm3=t['eps_Jm3'],
                            rho_kg_m3=t['rho_kg_m3'], F_A_Pa=t['F_A_Pa'],
                            rho_4pi_r0_2=t['rho_dm'],
                            pr_4pi_r0_2=t['rho_dm']*3.0,
                            pt_4pi_r0_2=-t['rho_dm'])
                       for t in tab],
            "energy_conditions": "WEC ρ<0 нарушен; радиальный NEC ρ+p_r=4ρ<0 нарушен; "
                                 "тангенц. NEC ρ+p_t=0 (насыщен); след T=0 (конформно)."},
        "checks": {
            "method_A_direct": "mpmath: C_plates=π²/720, C_torus=π²/45",
            "method_B_scipy_integrate": "Abel–Plana ζ(−3)=2∫t³/(e^{2πt}−1)dt=1/120; "
                                        "C=(π²/12)·2·(16 для тор) ·ζ(−3)",
            "rel_devs": {"C_plates": zek[0], "C_torus": zek[1]},
            "all_le_1e6": bool(zek[0] <= 1e-6 and zek[1] <= 1e-6)},
        "sphere_coefficient": "E_sph=+0.0923531·ħc/(2R) — ЛИТЕРАТУРНЫЙ "
                              "(Boyer 1968; Balian–Duplantier 1978; Milton 1980; "
                              "Leseduarte–Romeo 1996). Самостоятельно НЕ пересчитан "
                              "(честно): в JSON помечен как литературный.",
        "forms": [dict(form=f['form'], x=f['x'], d_m=f['d_m'],
                       formula=f['formula'], sign=f['sign'], note=f['note'],
                       rho_4pi_r0_2=f['rho'], pr_4pi_r0_2=f['pr'], pt_4pi_r0_2=f['pt'],
                       rho_SI_Jm3=f['rho_SI'], E_tot_J=f['E_tot_J'], m_eff_kg=f['m_eff_kg'])
                  for f in forms_out],
        "dce": {
            "modulation": "задача: ρ(t)=ρ_stat(1+ε sin ω_dt); честная адиаб.: "
                          "ρ(t)=ρ_stat(1+ε sin ω_dt)⁻⁴",
            "avg_over_period": "O(ε)=0 (⟨sin⟩=0); O(ε²): ⟨ρ⟩=ρ_stat·⟨(1+εsin)⁻⁴⟩"
                               " = ρ_stat(1+5ε²+105/8ε⁴+…) — положительная фотонная добавка",
            "resonance": "f_d=2f_n=nc/d; в полосе [1 МГц,10 ГГц] ⇒ d∈[3 см,300 м]",
            "tab": dce,
            "limit_honest": "для d=10 нм..10 мкм резонанс вне полосы (2f₁=c/d≫10 ГГц), "
                            "режим квазистатический, парный DCE подавлен; "
                            "энергия фотонов берётся из работы зеркала (сохраняется)."}
    }
    with open("/home/smboozha/portal_gun/research/r4_dce_values.json", "w",
              encoding="utf-8") as fp:
        json.dump(out, fp, ensure_ascii=False, indent=2, default=str)
    print("\n  JSON: research/r4_dce_values.json сохранён")
    print("="*94)

if __name__ == "__main__":
    main()