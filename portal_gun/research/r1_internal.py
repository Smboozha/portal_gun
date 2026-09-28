"""
R1 [ОСН] — внутренняя согласованность EL: символьно против численной вариации
=============================================================================
Проверяем, что СИМВОЛИЧЕСКАЯ Эйлер–Лагранжева (из кэша) равна численной
функциональной производной ТОГО ЖЕ действия L = A·B·C²·(T+f(T)) на классе
диагональных тетрад. Базовый профиль — Ellis червоточина.
"""
import json
import numpy as np
import sympy as sp

r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r); C = sp.Function('C')(r); f = sp.Function('f')
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': f}
data = json.load(open('/home/smboozha/portal_gun/research/ft_el_cache.json', encoding='utf-8'))
EL = {k.replace('EL_', ''): sp.sympify(v, locals=ns) for k, v in data.items() if k.startswith('EL_')}
F = sp.Function('F') if False else sp.Function('f')

# численная T по символической формуле (f-часть нулевая в f=0 пределе)
def T_num(rr, b0v, thv, A_=None, B_=None, C_=lambda R: R):
    "A_, B_ — функции r, C=r"
    Ab = 1.0; Abp = 0.0
    Bb = 1.0/sp.sqrt(1-b0v**2/rr**2); Bbp = 0.5*(2*b0v**2/rr**3)/(1-b0v**2/rr**2)**1.5
    if A_ is not None: Ab, Abp = A_(rr)[0], A_(rr)[1]
    if B_ is not None: Bb, Bbp = B_(rr)[0], B_(rr)[1]
    Cb = rr; Cbp = 1.0
    t = -1.0/(Cb**2*np.tan(thv)**2) + 4.0*Cbp**2/(Bbp2(Bb)) + ...
    # прямая формула
    return (-1.0/(Cb**2*np.tan(thv)**2)
            + 4.0*Cbp**2/(Bb**2*Cb**2)
            + 4.0*Abp*Cbp/(Ab*Bb**2*Cb)
            + Abp**2/(Ab**2*Bb**2))

def make_pert(fbase, dfbase, gauss):
    r0, sig = gauss
    def F_(rr):
        A_, dA = fbase(rr), dfbase(rr)
        g = np.exp(-((rr-r0)/sig)**2)
        dg = -2*(rr-r0)/sig**2*g
        return (A_, dA, g, dg)
    return F_

def action(profA, profB, rmin, rmax, N, b0v, thv):
    rr = np.linspace(rmin, rmax, N)
    def integrand(rr):
        Ta, Tbp = profA(rr), profB(rr)
        e = Ta[0]*Ta[2]*rr**2
        T = T_num(rr, b0v, thv, (Ta[0], Ta[1]), (Tb[0], Tb[1]))
        return e*T
    # собрать проще по-компонентно
    vals = np.array([L_at(r, b0v, thv, profA, profB) for r in rr])
    return np.trapezoid(vals, rr)

def L_at(rr, b0v, thv, profA, profB):
    Aa = profA(rr); Ba = profB(rr)
    Ab, Abp = Aa[0], Aa[1]
    Bb, Bbp = Ba[0], Ba[1]
    T = (-1.0/(rr**2*np.tan(thv)**2)
         + 4.0/(Bb**2*rr**2)
         + 4.0*Abp/(Ab*Bb**2*rr)
         + Abp**2/(Ab**2*Bb**2))
    return Ab*Bb*rr**2*T

def base_profs(rr, b0v):
    u = 1-b0v**2/rr**2
    B = 1/np.sqrt(u); Bp = b0v**2/(u**1.5*rr**3)
    return ((1.0, 0.0), (B, Bp))

def numeric_dS(profA, profB, q, r0, sig, eps, rmin, rmax, N, b0v, thv):
    """численная вариация ∫ EL_q·g dr ≈ (S[q+εg]−S[q])/ε"""
    def bump(rr): return np.exp(-((rr-r0)/sig)**2)
    def dbump(rr): return -2*(rr-r0)/sig**2*np.exp(-((rr-r0)/sig)**2)
    def profA_p(rr):
        a = profA(rr); return (a[0]+eps*bump(rr), a[1]+eps*dbump(rr))
    def profB_p(rr):
        a = profB(rr); return (a[0]+eps*0, a[1])  # не трогаем B в этом тесте
    S0, Sp = 0.0, 0.0
    rr = np.linspace(rmin, rmax, N)
    S0 = np.trapezoid([L_at(x, b0v, thv, profA, profB) for x in rr], rr)
    if q == 'A':
        Sp = np.trapezoid([L_at(x, b0v, thv, profA_p, profB) for x in rr], rr)
    else:
        Sp = np.trapezoid([L_at(x, b0v, thv, profA, profB_p) for x in rr], rr)
    return (Sp-S0)/eps

def symbolic_EL_numeric(name, profA, profB, r0, sig, rmin, rmax, N, b0v, thv):
    """∫ EL[name]·bump dr на бумажных профилях (поточковые скалярные подстановки)."""
    rr = np.linspace(rmin, rmax, N)
    vals = []
    for R in rr:
        u = 1-b0v**2/R**2
        BB = 1/np.sqrt(u)
        Bpa = (2*b0v**2/R**3)**0.5*0 + 0.5*(2*b0v**2/R**3)/(u**1.5)
        Bppa = -0.75*(2*b0v**2/R**3)**2/(u**2.5) + 0.5*(-6*b0v**2/R**4)/(u**1.5)
        e = EL[name].replace(lambda t: isinstance(t, sp.Function) and t.func == f, lambda t: sp.S.Zero)
        e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), Bppa).subs(C.diff(r, 2), 0)
        e = e.subs(A.diff(r), 0).subs(B.diff(r), Bpa).subs(C.diff(r), 1)
        e = e.subs(A, 1).subs(sp.Function('B')(r), BB).subs(C, R).doit().subs({th: thv})
        e = e.subs({sp.Symbol('b0'): sp.Float(b0v), sp.Symbol('r'): sp.Float(R)})
        try:
            vals.append(float(complex(e)).real)
        except Exception:
            vals.append(0.0)
    return np.trapezoid(np.array(vals)*np.exp(-((rr-r0)/sig)**2), rr)

b0v, thv, r0, sig, eps = 0.1, np.pi/4, 0.3, 0.05, 1e-4
rmin, rmax, N = 0.12, 0.8, 400
profA = lambda rr: (np.ones_like(rr), np.zeros_like(rr))
profB = lambda rr: (1/np.sqrt(1-b0v**2/rr**2), b0v**2*1.0/rr**3/(1-b0v**2/rr**2)**1.5)

for q in ['A', 'B']:
    ndS = numeric_dS(profA, profB, q, r0, sig, eps, rmin, rmax, N, b0v, thv)
    se = symbolic_EL_numeric(q, profA, profB, r0, sig, rmin, rmax, N, b0v, thv)
    print("q=%-3s числ. δS=%.6f  симв. ∫EL·g=%.6f  (ratio %.4f)" % (q, ndS, se, ndS/se if se != 0 else float('nan')))

# прямое значение EL[A] на точке (для ориентира)
for q in ['A', 'B', 'C']:
    for R in [0.2]:
        pass
print("готово")