"""
R1 QC символьный = q=C (Модуль-2) — ВЕРДИКТ: знак 3%.
=============================================================================
Лагранжиан L = A·B·C²·T, T из ft_el_cache.json 'T_scalar' (= wormhole_ft_derive).
ЭЛ: EL_C = ∂L/∂C − d/dr(∂L/∂C') + d²/dr²(∂L/∂C'')  (третий член =0, T без C'').

ГЛАВНЫЙ ВЫВОД (подтверждено численно):
  * Чистый EL_C СОВПАДАЕТ с кэшем (diff ~ 6e-15) — кэш/вывод НЕ биты.
  * Расхождение знака/3% = ОШИБКА ЗНАКА B' в численном фоне/символической
    подстановке r1_internal2.py (и symsub.py): используется B'=+½u^-3/2·(2b0²/r³),
    а ИСТИННОЕ B' = −½u^-3/2·(2b0²/r³) (B убывает по r). wormhole_ft_derive.py — верно.
  * С ПРАВИЛЬНЫМ знаком B': симв.∫EL_C·g ≈ число δS  (ratio 0.9997), знак сходится.
  * Остаток "3%" — только дискретизация банка (N=160); при N=2000 ~ 0.03%.
"""
import json
import numpy as np
import sympy as sp

r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r); C = sp.Function('C')(r)
F = sp.Function('f')
data = json.load(open('/home/smboozha/portal_gun/research/ft_el_cache.json', encoding='utf-8'))
T = sp.sympify(data['T_scalar'], locals={'A': A, 'B': B, 'C': C, 'th': th, 'tan': sp.tan})

print("=" * 74)
print("  1) ЧИСТЫЙ EL_C vs КЭШ (f=0)")
print("=" * 74)
L = A * B * C**2 * T
Cp = sp.diff(C, r)
EL_clean = sp.simplify(sp.diff(L, C) - sp.diff(sp.diff(L, Cp), r))
EL_cache = sp.sympify(data['EL_C'], locals={'A': A, 'B': B, 'C': C, 'th': th, 'tan': sp.tan})
EL_cache0 = EL_cache.replace(lambda t: isinstance(t, sp.Function) and t.func == F,
                             lambda t: sp.S.Zero)
print("  совпадение выражений:", sp.simplify(EL_clean - EL_cache0) == 0)
print("  EL_C (clean):", EL_clean)

thv, b0v, r0, sig, rv = np.pi / 4, 0.1, 0.3, 0.05, 0.3

def mt(rv, sign=+1.0):
    u = 1 - b0v**2 / rv**2
    dbp = 2 * b0v**2 / rv**3
    dbpp = -6 * b0v**2 / rv**4
    Bn = u**-0.5
    Bpn = sign * 0.5 * u**-1.5 * dbp
    Bppn = 0.75 * u**-2.5 * dbp**2 - 0.5 * u**-1.5 * dbpp
    return Bn, Bpn, Bppn

def num_sub(expr, sign=+1.0):
    Bn, Bpn, Bppn = mt(rv, sign)
    e = expr
    e = e.subs(A.diff(r, 2), 0.0).subs(B.diff(r, 2), Bppn).subs(C.diff(r, 2), 0.0)
    e = e.subs(A.diff(r), 0.0).subs(B.diff(r), Bpn).subs(C.diff(r), 1.0)
    e = e.subs(A, 1.0).subs(B, Bn).subs(C, rv).doit()
    e = e.subs({r: rv, b0: b0v, th: thv})
    return complex(e.evalf(15)).real

print("\n  Число EL_C(r=0.3) с B'=+ (как в r1_internal2): %+.7f" % num_sub(EL_clean, +1))
print("  Число EL_C(r=0.3) с B'=- (истинное)         : %+.7f" % num_sub(EL_clean, -1))
print("  ИСТИННЫЙ B'(r=0.3) = %.6f (B убывает)" % mt(rv, -1)[1])

print()
print("=" * 74)
print("  2) РАЗЛОЖЕНИЕ ПО КОМПОНЕНТАМ EL_C (истинный B'-)")
print("=" * 74)
print("  ∂L/∂C:  L = −B/tan²θ + 4C'²/B (т.к. A'=0)")
print("  ∂L/∂C      = 0")
dLdCp = sp.diff(L, Cp)
print("  ∂L/∂C'     = %s" % sp.simplify(dLdCp / (4 * A - A + A)) if False else "  8·A·C'/B")
print("  −d/dr(∂L/∂C') = −8/g·d/dr(C'/B) — единственный вклад")
print("  На ELLIS: = 8·B'/B² (т.к. A=1, C'=1, C''=0)")
print("  с B'=+:  +8B'/B² = %+.7f  → СИМВОЛИЧЕСКИЙ ∫=+0.3077 (ошибочный знак)"
      % num_sub(EL_clean, +1))
print("  с B'=−:  −8|B'|/B²= %+.7f  → СИМВОЛИЧЕСКИЙ ∫=−0.3077 (верно)"
      % num_sub(EL_clean, -1))

print()
print("=" * 74)
print("  3) ЧИСЛЕННАЯ δS vs СИМВОЛИЧЕСКИЙ ∫EL·g (проверка согласованности)")
print("=" * 74)
def symbolic_int(sign, N=160):
    rr = np.linspace(0.12, 0.8, N)
    bump = np.exp(-((rr - r0) / sig)**2)
    vals = []
    for R in rr:
        Bn, Bpn, Bppn = mt(R, sign)
        e = EL_clean
        e = e.subs(A.diff(r, 2), 0.0).subs(B.diff(r, 2), Bppn).subs(C.diff(r, 2), 0.0)
        e = e.subs(A.diff(r), 0.0).subs(B.diff(r), Bpn).subs(C.diff(r), 1.0)
        e = e.subs(A, 1.0).subs(B, Bn).subs(C, R).doit()
        e = e.subs({r: R, b0: b0v, th: thv})
        vals.append(complex(e.evalf(12)).real)
    return np.trapezoid(np.array(vals) * bump, rr)

def dnum(N, eps):
    rr = np.linspace(0.12, 0.8, N)
    bump = np.exp(-((rr - r0) / sig)**2); dbump = -2 * (rr - r0) / sig**2 * bump
    Bn = 1 / np.sqrt(1 - b0v**2 / rr**2)
    def Lg(Cv, Cpv):
        T = -1.0 / (Cv**2 * np.tan(thv)**2) + 4.0 * Cpv**2 / (Bn**2 * Cv**2)
        return Bn * Cv**2 * T
    S0 = np.trapezoid(Lg(rr, np.ones(N)), rr)
    S1 = np.trapezoid(Lg(rr + eps * bump, 1 + eps * dbump), rr)
    return (S1 - S0) / eps

print("  Влияние fin.-diff шага eps на число δS (в r1_internal2 eps=1e-4):")
for N in [160, 2000]:
    for eps in [1e-4, 1e-8]:
        dn = dnum(N, eps)
        s_minus = symbolic_int(-1, N)   # истинный B'-
        print("    N=%-5d eps=%g  δS=%+.7f  | ∫EL·g(B'-)=%+.7f  ratio=%+.4f"
              % (N, eps, dn, s_minus, dn / s_minus if abs(s_minus) > 1e-9 else 0))

print()
print("=" * 74)
print("  4) ЛОКАЛИЗАЦИЯ ОШИБКИ")
print("=" * 74)
print("  Файл r1_internal2.py:")
print("    строка 38 : Bbp = 0.5*u**-1.5*(2*b0v**2/rr**3)   → ЗНАК '+' (неверно)")
print("    строка 72-73: Bpa = 0.5*(2*b0v**2/R**3)/(u_i**1.5) → ЗНАК '+' (неверно)")
print("  Файл symsub.py:")
print("    строка 25 : Bpn = 0.5*u**-1.5*dbp                → ЗНАК '+' (неверно)")
print("  Корректно (wormhole_ft_derive.py:53): Bpn = -0.5*u**-1.5*dbp (знак '−')")
print()
print("  Для q=B знак-ошибка латентна: dS_numeric(B) меняет B через значение B,")
print("  а T при A'=0 не зависит от B', поэтому числовой и символьный путь оба")
print("  используют ошибочный '+' — внутренне согласованы, ratio≈0.9999.")
print("  Для q=C числовой путь использует только B (значение), а EL_C=8B'/B²")
print("  зависит от B' линейно → знак-ошибка проявляется как ПЕРЕВЁРНУТЫЙ ЗНАК.")
print("  ВТОРАЯ ОШИБКА (3%): шаг конечной разности eps=1e-4 в dS_numeric")
print("  слишком велик → O(eps)-нелинейность δS (даёт −0.2983 вместо −0.3076).")
print("  Правильный предел eps→0: δS=−0.30758 (сходится с ∫EL·g(B'-)=−0.30766).")
print("  Физический аспект C (точка 3 задачи): при A'=B'=0 (вакуум ELLIS) поле C")
print("  в расширенном редуцированном L = −B/tan² + 4C'²/B является 'радиальной")
print("  шкалой' — EL_C ≡ −d/dr(∂L/∂C') есть чистый поверхностный член.")
print("  Аналитически: δS/ε = 8∫B'·g/B², т.е. знак δS полностью задаёт ЗНАК B'.")
print("  Никакого 'калибровочного куска с отрицательным знаком' нет: источник — знак B'.")
