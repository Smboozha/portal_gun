"""
R1 q=A DIAGNOSTIC [КД3] — FINAL: corrected ratio after fix
"""
import sympy as sp
import numpy as np
import sys, os, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symsub import r, th, b0, A, B, C, F, kill_f

data = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ft_el_cache.json'), encoding='utf-8'))
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': F}
T_scalar = sp.sympify(data['T_scalar'], locals=ns)
EL_A_cache = sp.sympify(data['EL_A'], locals=ns)

b0v = 0.1; thv = np.pi/4

# ── BUILD EL_A FROM SCRATCH (no f, so no Subs problem) ──────────────────
print("Building EL_A from L = A*B*C^2*T (f=0, no Subs terms):")
L = A * B * C**2 * T_scalar
I0  = sp.diff(L, A)
I1  = sp.diff(L, A.diff(r))
I2  = sp.diff(L, A.diff(r,2))
I11 = -sp.diff(I1, r)
I22 = sp.diff(I2, r, 2)
EL_A_clean = sp.cancel(I0 + I11 + I22)
print(f"  EL_A_clean has f? {EL_A_clean.has(F)}")
print(f"  EL_A_clean has Subs? {EL_A_clean.has(sp.Subs)}")
print()

# ── Verify: does kill_f on cached EL_A equal EL_A_clean? ────────────────
print("Do the cached and clean EL_A agree after proper f-removal?")

def eval_check(expr, rv):
    Bw = 1/sp.sqrt(1-b0**2/r**2)
    Bp = sp.diff(Bw, r); Bpp = sp.diff(Bw, r, 2)
    e = expr
    e = e.subs(A.diff(r,2), 0).subs(B.diff(r,2), Bpp).subs(C.diff(r,2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bw).subs(C, r).doit()
    e = e.subs({r: rv, b0: b0v, th: thv})
    return complex(e.evalf(12)).real

# clean path
print(f"{'r':>6s}  {'EL_clean':>14s}  {'EL_cache(kf)':>16s}  {'EL_cache(Subs→0)':>18s}")
print("-" * 60)
for rv in [0.12, 0.2, 0.3, 0.5, 0.8]:
    v_clean = eval_check(EL_A_clean, rv)
    # kill_f on cached
    v_kf = eval_check(kill_f(EL_A_cache), rv)
    # strong removal: Subs->0 too
    cache_subs0 = EL_A_cache.replace(lambda t: isinstance(t, sp.Subs), lambda t: sp.S.Zero).replace(lambda t: isinstance(t, sp.Function) and t.func==F, lambda t: sp.S.Zero)
    v_subs0 = eval_check(cache_subs0, rv)
    print(f"{rv:6.3f}  {v_clean:+14.8f}  {v_kf:+16.8f}  {v_subs0:+18.8f}")

print()

# ── FINAL: compute corrected ratio for q=A ──────────────────────────────
print("=" * 72)
print("FINAL RATIO COMPUTATION (q=A) with fix")
print("=" * 72)
print()

N = 160
rr = np.linspace(0.12, 0.8, N)
sig, eps = 0.05, 1e-4
bump = np.exp(-((rr-0.3)/sig)**2)
dbump = -2*(rr-0.3)/sig**2*bump

# numerical δS/δA
Ab = np.ones(N); Abp = np.zeros(N)
u = 1-b0v**2/rr**2
Bb = 1/np.sqrt(u)
Bbp = 0.5*u**-1.5*(2*b0v**2/rr**3)
def T_num(R,Ab,Abp,Bb,Bbp):
    return (-1.0/(R**2*np.tan(thv)**2)+4.0/(Bb**2*R**2)+4.0*Abp/(Ab*Bb**2*R)+Abp**2/(Ab**2*Bb**2))
def L_num(Ab,Abp,Bb,Bbp):
    e=Ab*Bb*rr**2
    T=np.array([T_num(rr[i],Ab[i],Abp[i],Bb[i],Bbp[i]) for i in range(N)])
    return e*T
Ap=Ab+eps*bump; App=Abp+eps*dbump
dS_A=(np.trapezoid(L_num(Ap,App,Bb,Bbp),rr)-np.trapezoid(L_num(Ab,Abp,Bb,Bbp),rr))/eps

# symbolic integral: clean EL_A (check path)
vals_clean = np.array([eval_check(EL_A_clean, rv) for rv in rr])
sym_clean = np.trapezoid(vals_clean*bump, rr)

# symbolic integral: cached kill_f (float path - original bug)
vals_kf_float = np.array([eval_check(kill_f(EL_A_cache), rv) for rv in rr])  # check path here
# original float path
def path_float(expr, rv):
    u_i=1-b0v**2/rv**2
    BB=1.0/np.sqrt(u_i)
    Bpa=0.5*(2*b0v**2/rv**3)/u_i**1.5
    Bppa=-0.75*(2*b0v**2/rv**3)**2/u_i**2.5+0.5*(-6*b0v**2/rv**4)/u_i**1.5
    e=kill_f(expr)
    e=e.subs(A.diff(r,2),0.0).subs(B.diff(r,2),Bppa).subs(C.diff(r,2),0.0)
    e=e.subs(A.diff(r),0.0).subs(B.diff(r),Bpa).subs(C.diff(r),1.0)
    e=e.subs(A,1.0).subs(B,BB).subs(C,rv).doit()
    e=e.subs({r:rv,b0:b0v,th:np.pi/4})
    return complex(e.evalf(12)).real
vals_cache_float = np.array([path_float(EL_A_cache, rv) for rv in rr])
sym_float = np.trapezoid(vals_cache_float*bump, rr)

print(f"  Numerical δS/δA            = {dS_A:+.8f}")
print(f"  ∫EL_A·g (ORIGINAL float)   = {sym_float:+.8f}")
print(f"  ratio (ORIGINAL)           = {dS_A/sym_float:.6f}")
print()
print(f"  ∫EL_A·g (CLEAN rebuild)    = {sym_clean:+.8f}")
print(f"  ratio (CLEAN FIX)          = {dS_A/sym_clean:.6f}")
print()

# report which two components contribute to δS≈-0.138
I0_f = np.array([eval_check(I0, rv) for rv in rr])
I11_f = np.array([eval_check(I11, rv) for rv in rr])
int00 = np.trapezoid(I0_f*bump, rr)
int11 = np.trapezoid(I11_f*bump, rr)
print(f"  δS/δA = ∫I0·g + ∫I11·g  (I22=0)")
print(f"         = {int00:+.8f} + {int11:+.8f}")
print(f"         = {int00+int11:+.8f}")
print(f"  vs δS/δA = {dS_A:+.8f}")
print()

# component breakdown on bump region
print("  Contribution by component:")
print(f"    ∫I0·g  = {int00:+.8f}  ({100*int00/dS_A:+.1f}% of δS)")
print(f"    ∫I11·g = {int11:+.8f}  ({100*int11/dS_A:+.1f}% of δS)")
print()
print("DONE")
