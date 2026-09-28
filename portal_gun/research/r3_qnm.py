"""
R3 [ОСН] — QNM скалярного возмущения для f(T)-червоточины (A=1)
================================================================
Скалярное возмущение Φ=ψ(r)Y_lm e^{-iwt}/r на ds²=−dt²+B²dr²+r²dΩ²
сводится (строгий вывод, см. ниже) к:
    d²ψ/dr*² + (ω² − V(r*))ψ = 0,   r* = ∫B dr,
    V(r) = ℓ(ℓ+1)/r² − B'/(B³r),   а для наших форм:
    'ellis' → V = ℓ(ℓ+1)/r² + b0²/r⁴,   r∈[b0,∞)
    'const' → V = ℓ(ℓ+1)/r² + b0/(2r³), r∈[b0,∞)
(вывод: □Φ=0, g=diag(-1,B²,r²,r²sin²θ), √−g=BC²sinθ; член −B'/(B³r)
  превращается в +b0²/r⁴ (ellis) / +b0/2r³ (const), как в аналитике.)

Главный вопрос: есть ли осцилляторные QNM (барьерный максимум V в r*)?
"""
import numpy as np
from scipy.integrate import quad

b0 = 0.10  # 10 см

def V_ellis(r, l):
    return l*(l+1)/r**2 + b0**2/r**4

def V_const(r, l):
    return l*(l+1)/r**2 + b0/(2*r**3)

def B_const(r): return 1.0/np.sqrt(1-b0/r)
def B_ellis(r): return 1.0/np.sqrt(1-b0**2/r**2)

def rstar_ellis(r_): return np.sqrt(r_**2-b0**2)
def tortoise_points(mode, l, N=400):
    if mode == 'ellis':
        r = np.linspace(b0+1e-6, 8.0, N)
        rs = np.sqrt(r**2-b0**2)
    else:
        r = np.linspace(b0+1e-6, 8.0, N)
        vals = [quad(B_const, b0+1e-7, x)[0] for x in r]
        vals = np.array(vals)
        rs = vals
    return r, rs

def check_monotone(mode, l):
    r, rs = tortoise_points(mode, l)
    V = V_ellis(r, l) if mode == 'ellis' else V_const(r, l)
    dV = np.gradient(V, rs)
    return (dV <= 1e-12).all(), float(V[0]), float(V[-1]), V

def list_report():
    print("="*70)
    print("  R3 [ОСН] QNM скалярного возмущения, A=1, b0=10 см")
    print("="*70)
    lines = []
    for mode in ['ellis', 'const']:
        for l in [0, 1, 2]:
            mono, V0, Vinf, V = check_monotone(mode, l)
            r, rs = tortoise_points(mode, l)
            print("  %-6s l=%d   V(r0)=%+.4f   V(inf)=%+.2e   монотонно убывающ.: %s"
                  % (mode, l, V0, Vinf, "ДА" if mono else "НЕТ(есть барьер)"))
            if not mono:
                im = np.argmax(V)
                print("      → барьер в точке r*=%+.3f, Vmax=%+.4f (WKB возможен)"
                      % (rs[im], V[im]))
    print()
    print("  Вывод: для A=1 обе формы имеют монотонно убывающий V(r*)")
    print("  → НЕТ осцилляторных (комплексно-частотных) QNM — чёрно-дырного")
    print("    ringdown (звон) отсутствует; отклик — чистое затухание.")
    print("  Следствие для наблюдательного R4: у проходимой безызлучающей")
    print("    горловины нет канонической «гравитационной мелодии».")
    print("="*70)

if __name__ == '__main__':
    list_report()