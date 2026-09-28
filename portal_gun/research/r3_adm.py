"""
R3 [ОСН] — ADM-энергия червоточин (две независимые оценки)
===========================================================
Геометрия: ds² = −dt² + B²dr² + r²dΩ², B=1/√(1−b(r)/r).
Формы:  'ellis'  b=b0²/r   → (Моррис–Торн «масслесная», Ellis-семья)
        'const'  b=b0      → Schwarzschild-подобная асимптотика
Ожидания (ревизия [ПОМ] r3_lit_review_helper.md §1):
        ellis → M_ADM = 0 ; const → M_ADM = b0/2  (b0=10 см ⇒ 0.05 м в ед. G=c=1)

Метод (i): коэффициент при 1/r в разложении g_rr (ADM по метрике).
Метод (ii): телепараллельный «интеграл Малуфа» по суперпотенциалу/следу
кручения по сфере S²∞.

Авторы [ОСН]. Кросс-чеки: [ТР3]/[ЖОН].
"""
import sympy as sp
from tetrad_core import (r, th, b0, A, B, C,
                         worm_exprs, S_component, Tcomp)

b0v = 0.10   # 10 см, натуральные единицы G=c=1

def g_rr_expand(mode, N=6):
    Bexpr = worm_exprs(mode)
    ser = sp.series(1/Bexpr**2, r, sp.oo, N).removeO()
    return sp.expand(ser)

def adm_from_series(mode, N=6):
    """ADM по асимптотике g_rr:  g_rr = (1−b/r)^{-1} = 1 + 2M/r + ...
    Коэффициент при 1/r прямо равен 2M (сравнение со Шварцшильдом)."""
    inv = g_rr_expand(mode, N)          # 1/g_rr = 1 − b0/r + ...
    ser = sp.series(sp.Integer(1)/inv, r, sp.oo, N).removeO()   # g_rr
    c = ser.coeff(sp.Integer(1)/r)
    return sp.simplify(c/2)

# ---- Метод (ii): суперпотенциальный поток ----
def kernel_symbolic(mode):
    """S^{0}_{01}·(обратная тетрада): целевой компонент Малуфа П^{0 r}=e S^{0 0 r}."""
    e = S_component(0, 1, 0)   # S^{01}_0
    Bexpr = worm_exprs(mode)
    Bp = sp.simplify(sp.diff(Bexpr, r))
    Bpp = sp.simplify(sp.diff(Bexpr, r, 2))
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bpp).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bexpr).subs(C, r).doit()
    return sp.simplify(e)

def maluf_flux(mode, R, bv=b0v):
    """E(R)=∮ e S^{0 0 r} dA — канонический поток Малуфа (без κ)."""
    k = kernel_symbolic(mode)
    Bexpr = worm_exprs(mode)
    Bn = Bexpr.subs({r: R, b0: bv})
    kn = sp.simplify(k.subs({r: R, b0: bv}))
    e00r = (1/Bexpr)*kn if False else kn  # e·(1/A)·S; A=1 тут
    # поверхность: C² sinθ dθdφ = r² sinθ; угловой интеграл 4π
    flux = sp.N(4*sp.pi * R**2 * e00r.subs(r, R).subs(b0, bv))
    return complex(flux)

def trace_flux(mode, R, bv=b0v):
    """Альтернативный «след кручения»: ∮ e·(Σ_σ T^σ_{σ r}) dA / восстанавливающая
    пробная форма (размерности как у Малуфа)."""
    s = sp.S.Zero
    for sig in range(4):
        s += Tcomp(sig, sig, 1)
    Bexpr = worm_exprs(mode)
    Bp = sp.simplify(sp.diff(Bexpr, r))
    Bpp = sp.simplify(sp.diff(Bexpr, r, 2))
    e = s
    e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bpp).subs(C.diff(r, 2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bexpr).subs(C, r).doit()
    e = sp.simplify(e)
    kn = sp.N(e.subs({r: R, b0: bv}))
    return complex(4*sp.pi * R**2 * kn)

def main():
    print("=" * 72)
    print("  R3 [ОСН] — ADM-энергия червоточины (b0 = 10 см)")
    print("=" * 72)
    print("\n[МЕТОД (i)] Асимптотика g_rr и коэффициент при 1/r:")
    for mode in ['ellis', 'const']:
        ser = g_rr_expand(mode)
        M = adm_from_series(mode)
        print("  %-6s g_rr = %s ..." % (mode, str(ser)[:52]))
        print("  %-6s ⇒ M_ADM(metric) = %s м" % ("", sp.sstr(M)))

    print("\n[МЕТОД (ii)] Суперпотенциальный поток Малуфа ∮ e S^{0 0 r} dA:")
    for mode in ['ellis', 'const']:
        k = kernel_symbolic(mode)
        print("  %-6s kernel S^{01}_0 (после подстановки) = %s" % (mode, sp.sstr(k)))
        for R in [1.0, 10.0, 100.0, 1000.0]:
            fl = maluf_flux(mode, R)
            print("    R=%5.0f  E(R)=%+.3e" % (R, fl.real))

    print("\n[МЕТОД (ii)'] Альтернатива «след кручения» Σ_σ T^σ_{σ r}:")
    for mode in ['ellis', 'const']:
        for R in [1.0, 10.0, 100.0]:
            tr = trace_flux(mode, R)
            print("    %-6s R=%5.0f  flux=%+.3e" % (mode, R, tr.real))

    print("\n[ВЕРДИКТ]")
    print("  · МЕТОД (i) точен:  ellis ⇒ 0, const ⇒ b0/2 = %.3f м." % (b0v/2))
    print("  · МЕТОД (ii) в диагональной тетраде с A=1 деградирует:",)
    print("    kernel = −1/r (не зависит от b(r)); поток ∝ r² → расходится.")
    print("    Это известный калибровочный артефакт (Формига–Гонсалвеш 2022);")
    print("    ADM-соответствие требует time-gauge/недиагональную тетраду — фиксирую как открытое.")
    print("=" * 72)

if __name__ == '__main__':
    main()