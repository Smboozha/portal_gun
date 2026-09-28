"""
R1 [ОСН] — внутренняя согласованность EL (численная вариация = символьная ЭЛ)
=============================================================================
Богатство: L = A·B·C²·T (f=0 предел). Численная функциональная производная
(см. классическая граница: δS/δq·дискретизация) должна равняться ∫EL_q·g dr.
"""
import json
import numpy as np
import sympy as sp
import symsub
from symsub import r, th, b0, A, B, C, F, worm_subs, kill_f

data = json.load(open('/home/smboozha/portal_gun/research/ft_el_cache.json', encoding='utf-8'))
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': F}
EL = {k.replace('EL_', ''): sp.sympify(v, locals=ns) for k, v in data.items() if k.startswith('EL_')}

b0v = 0.1; thv = np.pi/4
rmin, rmax, N = 0.12, 0.8, 160
rr = np.linspace(rmin, rmax, N)
r0, sig, eps = 0.3, 0.05, 1e-5
bump = np.exp(-((rr-r0)/sig)**2)
dbump = -2*(rr-r0)/sig**2*bump

def T_at(R, Ab, Abp, Bb, Bbp):
    return (-1.0/(R**2*np.tan(thv)**2)
            + 4.0/(Bb**2*R**2)
            + 4.0*Abp/(Ab*Bb**2*R)
            + Abp**2/(Ab**2*Bb**2))

def L_grid(Ab, Abp, Bb, Bbp):
    e = Ab*Bb*rr**2
    T = np.array([T_at(R, Ab[i], Abp[i], Bb[i], Bbp[i]) for i, R in enumerate(rr)])
    return e*T

Ab = np.ones(N); Abp = np.zeros(N)
u = 1-b0v**2/rr**2
Bb = 1/np.sqrt(u)
Bbp = -0.5*u**-1.5*(2*b0v**2/rr**3)

def dS_numeric(q, style=None):
    """Численная функциональная производная (центральная разность, eps→точность).
    O(eps^2) вместо forward O(eps): убирает 3%-погрешность (модуль-1)."""
    he = 0.5*eps
    if q == 'A':
        bp  = Ab + eps*bump;  bm  = Ab - eps*bump
        bppp = Abp + eps*dbump; bppm = Abp - eps*dbump
        return (np.trapezoid(L_grid(bp, bppp, Bb, Bbp), rr)
                - np.trapezoid(L_grid(bm, bppm, Bb, Bbp), rr))/(2*eps)
    if q == 'B':
        Bp2p = Bb + eps*bump;  Bp2m = Bb - eps*bump
        Bpp2p = Bbp + eps*dbump; Bpp2m = Bbp - eps*dbump
        return (np.trapezoid(L_grid(Ab, Abp, Bp2p, Bpp2p), rr)
                - np.trapezoid(L_grid(Ab, Abp, Bp2m, Bpp2m), rr))/(2*eps)
    if q == 'C':
        # варьируем C→r+εg (метрика: C=C(r))
        def L_C(Cv, Cpv):
            T = np.array([
                -1.0/(Cv[i]**2*np.tan(thv)**2)
                + 4.0*Cpv[i]**2/(Bb[i]**2*Cv[i]**2)
                + 4.0*Abp[i]*Cpv[i]/(Ab[i]*Bb[i]**2*Cv[i])
                + Abp[i]**2/(Ab[i]**2*Bb[i]**2) for i, R in enumerate(rr)])
            return Ab*Bb*Cv**2*T
        C2p = rr + eps*bump;  Cp2p = 1 + eps*dbump
        C2m = rr - eps*bump;  Cp2m = 1 - eps*dbump
        return (np.trapezoid(L_C(C2p, Cp2p), rr)
                - np.trapezoid(L_C(C2m, Cp2m), rr))/(2*eps)
    raise ValueError(q)

def symbolic_int(name):
    """∫EL[name]·(функция-банка по всем узлам ELLIS)·bump dr.
    CHECK-путь (фикс 2026-09-07, модули-1/2/3):
      B^″,B^′,B ДЕРЖАТ символьными функциями r (не float!) до самого конца:
      e→ subs(B,diff(r,2),Bpp(r)) → subs(B.diff(r),Bp(r)) → subs(B,B(r))
      → координаты {r:R} ПОСЛЕДНИМИ. Старый float-путь (B.diff(r)→число)
      терял r-структуру у префакторов → ratio 2.7443 (q=A).
    Знак производных B′ — ОТРИЦАТЕЛЬНЫЙ (B убывает с r), см. symsub.worm_vals.
    """
    Bsym = B._subs_dummy if False else sp.Function('B')(r)
    # Аналитические B(r), B'(r), B''(r) как функции r (в символах b0):
    ue = 1 - b0**2/r**2
    bfunc = 1/sp.sqrt(ue)
    bpfun = sp.diff(bfunc, r)
    bppfun = sp.diff(bpfun, r)
    vals = []
    for i, R in enumerate(rr):
        if i % 40 == 0:
            print('  символ.', i, flush=True)
        e = kill_f(EL[name])
        # узлы: замена производных СИМВОЛЬНЫМИ выражениями от r
        e = e.subs(A.diff(r, 2), 0).subs(B.diff(r, 2), bppfun).subs(C.diff(r, 2), 0)
        e = e.subs(A.diff(r), 0).subs(B.diff(r), bpfun).subs(C.diff(r), 1)
        e = e.subs(A, 1).subs(B, bfunc).subs(C, r)
        # координаты последней операцией
        e = e.subs({r: R, b0: b0v, th: np.pi/4}).doit()
        vals.append(complex(e.evalf(14)).real)
    vals = np.array(vals)
    return np.trapezoid(vals*bump, rr)

for q in ['A', 'B', 'C']:
    try:
        nd = dS_numeric(q)
        se = symbolic_int(q)
        print("q=%-3s  числ.δS=%.6f  симв.∫EL·g=%.6f  ratio=%s" %
              (q, nd, se, ("%.4f" % (nd/se)) if abs(se) > 1e-9 else "undefined"))
    except Exception as ex:
        import traceback; traceback.print_exc()
print("done")