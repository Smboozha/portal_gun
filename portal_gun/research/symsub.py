"""
symsub.py [ОСН] — ЕДИНЫЙ канонический путь подстановки червоточины.
Порядок жёсткий: производные-узлы → сами функции (float) → координаты (float).
Стиль 'flt': B→1/sqrt(1-b^2/r^2) float сразу; 'sym': сначала символьный B.
"""
import sympy as sp
import numpy as np

r, th, b0 = sp.symbols('r th b0', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r); C = sp.Function('C')(r)
F = sp.Function('f')

def B_profile(mode, rv, bv):
    if mode == 'const':
        return 1.0/np.sqrt(1.0-bv/rv)
    return 1.0/np.sqrt(1.0-bv**2/rv**2)

def worm_vals(mode, rv, bv):
    u = 1.0-bv**2/rv**2 if mode != 'const' else 1.0-bv/rv
    if mode == 'const':
        dbp, dbpp = bv/rv**2, -2.0*bv/rv**3
    else:
        dbp, dbpp = 2.0*bv**2/rv**3, -6.0*bv**2/rv**4
    Bn = 1.0/np.sqrt(u)
    Bpn = -0.5*u**-1.5*dbp
    Bppn = 0.75*u**-2.5*dbp**2 - 0.5*u**-1.5*dbpp
    return dict(A=1.0, B=Bn, Bp=Bpn, Bpp=Bppn, C=rv)

def worm_subs(expr, mode, rv, bv, thv, style='flt'):
    """Каноническая подстановка. Возвращает float-complex value."""
    v = worm_vals(mode, rv, bv)
    e = expr
    e = e.subs(A.diff(r, 2), 0.0).subs(B.diff(r, 2), v['Bpp']).subs(C.diff(r, 2), 0.0)
    e = e.subs(A.diff(r), 0.0).subs(B.diff(r), v['Bp']).subs(C.diff(r), 1.0)
    e = e.subs(A, 1.0).subs(B, v['B']).subs(C, rv)
    e = e.doit()
    e = e.subs({r: rv, b0: bv, th: thv})
    return complex(e.evalf(14)).real

def kill_f(expr):
    return expr.replace(lambda t: isinstance(t, sp.Function) and t.func == F,
                        lambda t: sp.S.Zero)