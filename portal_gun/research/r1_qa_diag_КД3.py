"""
R1 q=A DIAGNOSTIC [КД3] — компонентный разбор EL_A и локализация бага ratio=2.7443
"""
import sympy as sp
import numpy as np
import sys, os, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symsub
from symsub import r, th, b0, A, B, C, F, kill_f

data = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ft_el_cache.json'), encoding='utf-8'))
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': F}
T_scalar = sp.sympify(data['T_scalar'], locals=ns)
EL_A_cache = sp.sympify(data['EL_A'], locals=ns)

b0v = 0.1; thv = np.pi / 4

A1 = A.diff(r); A2 = A.diff(r, 2)
B1 = B.diff(r); B2 = B.diff(r, 2)
C1 = C.diff(r); C2 = C.diff(r, 2)

# ── PART 1: Build EL_A from L = A*B*C^2*T ───────────────────────────────
print("=" * 72)
print("PART 1: Component decomposition EL_A = I0 + I11 + I22")
print("=" * 72)

L = A * B * C**2 * T_scalar
I0  = sp.diff(L, A)          # dL/dA
I1  = sp.diff(L, A1)         # dL/dA'
I2  = sp.diff(L, A2)         # dL/dA''
I11 = -sp.diff(I1, r)        # -d/dr(dL/dA')
I22 =  sp.diff(I2, r, 2)     # d^2/dr^2(dL/dA'')

EL_A_rebuilt = sp.cancel(I0 + I11 + I22)
print("I0 = dL/dA:")
sp.pprint(sp.simplify(I0))
print("\nI1 = dL/dA':")
sp.pprint(sp.simplify(I1))
print("\nI11 = -d/dr(dL/dA'):")
sp.pprint(sp.simplify(I11))
print("\nI2 = dL/dA'':")
sp.pprint(sp.simplify(I2))
print("\nI22 = d^2/dr^2(dL/dA''):")
sp.pprint(sp.simplify(I22))
print()

# ── PART 1b: Verify EL_A_rebuilt == EL_A_cache (after kill_f + simplify) ──
el_nof = kill_f(EL_A_cache)
el_nof_s = kill_f(sp.cancel(EL_A_rebuilt))
# Quick numerical check
def eval_ellis(expr, rv, bv=Ellipsis):
    """Evaluate expression at ELLIS parameters."""
    u_i = 1 - b0v**2 / rv**2
    BB = 1 / np.sqrt(u_i)
    Bpa = 0.5 * (2 * b0v**2 / rv**3) / u_i**1.5
    Bppa = -0.75 * (2 * b0v**2 / rv**3)**2 / u_i**2.5 + 0.5 * (-6 * b0v**2 / rv**4) / u_i**1.5
    e = kill_f(expr)
    e = e.subs(A2, 0).subs(B2, Bppa).subs(C2, 0)
    e = e.subs(A1, 0).subs(B1, Bpa).subs(C1, 1)
    e = e.subs(A, 1).subs(B, BB).subs(C, rv)
    e = e.doit()
    e = e.subs({r: rv, b0: b0v, th: np.pi / 4})
    return complex(e.evalf(12)).real

print("PART 1b: Numerical check EL_A_rebuilt vs EL_A_cache at r=0.3:")
v_cache = eval_ellis(EL_A_cache, 0.3)
v_rebuilt = eval_ellis(EL_A_rebuilt, 0.3)
print(f"  EL_A_cache(r=0.3)    = {v_cache:+.8f}")
print(f"  EL_A_rebuilt(r=0.3)  = {v_rebuilt:+.8f}")
print(f"  match: {'YES' if abs(v_cache - v_rebuilt) < 1e-6 else 'NO — diff=' + str(v_cache - v_rebuilt)}")
print()

# ── PART 2: Component values on grid ──────────────────────────────────────
print("=" * 72)
print("PART 2: Component values on ELLIS background")
print("=" * 72)

test_rs = [0.12, 0.15, 0.2, 0.3, 0.5, 0.7, 0.8]
print(f"{'r':>6s}  {'I0':>14s}  {'I11':>14s}  {'I22':>14s}  {'sum':>14s}  {'EL_A_cache':>14s}")
print("-" * 72)

I0_vals = []
I11_vals = []
I22_vals = []
EL_vals = []

for rv in test_rs:
    v0 = eval_ellis(I0, rv)
    v11 = eval_ellis(I11, rv)
    v22 = eval_ellis(I22, rv)
    vEL = eval_ellis(EL_A_cache, rv)
    I0_vals.append(v0); I11_vals.append(v11); I22_vals.append(v22); EL_vals.append(vEL)
    s = v0 + v11 + v22
    print(f"{rv:6.3f}  {v0:+14.8f}  {v11:+14.8f}  {v22:+14.8f}  {s:+14.8f}  {vEL:+14.8f}")

print()
print("  At A'=0, A''=0: I1 = dL/dA' = 0 analytically (A' appears as factor).")
print("  Therefore I11 = 0, I22 = 0, EL_A = I0 ONLY.")
print()

# ── PART 2b: Verify I1 and I2 are zero when A'=A''=0 ─────────────────────
print("  I1 value at r=0.3 (A'=0, A''=0): ", end="")
print(f"{eval_ellis(I1, 0.3):.2e}")
print("  I2 value at r=0.3 (A'=0, A''=0): ", end="")
print(f"{eval_ellis(I2, 0.3):.2e}")
print()

# ── PART 3: Float vs Check path comparison ────────────────────────────────
print("=" * 72)
print("PART 3: Float path (r1_internal2.py) vs Check path (r1_check.py)")
print("=" * 72)

def path_float(expr, rv):
    """Exact replica of r1_internal2.py symbolic_int inner loop."""
    u_i = 1 - b0v**2 / rv**2
    BB = 1.0 / np.sqrt(u_i)
    Bpa = 0.5 * (2 * b0v**2 / rv**3) / u_i**1.5
    Bppa = -0.75 * (2 * b0v**2 / rv**3)**2 / u_i**2.5 + 0.5 * (-6 * b0v**2 / rv**4) / u_i**1.5
    e = kill_f(expr)
    e = e.subs(A2, 0.0).subs(B2, Bppa).subs(C2, 0.0)
    e = e.subs(A1, 0.0).subs(B1, Bpa).subs(C1, 1.0)
    e = e.subs(A, 1.0).subs(B, BB).subs(C, rv)
    e = e.doit()
    e = e.subs({r: rv, b0: b0v, th: np.pi / 4})
    return complex(e.evalf(12)).real

def path_check(expr, rv):
    """Exact replica of r1_check.py worm_subs inner loop."""
    Bw = sp.Integer(1) / sp.sqrt(1 - b0**2 / r**2)
    Bp = sp.simplify(sp.diff(Bw, r))
    Bpp = sp.simplify(sp.diff(Bw, r, 2))
    e = kill_f(expr)
    e = e.subs(A2, 0).subs(B2, Bpp).subs(C2, 0)
    e = e.subs(A1, 0).subs(B1, Bp).subs(C1, 1)
    e = e.subs(A, 1).subs(B, Bw).subs(C, r)
    e = e.doit()
    e = e.subs({r: rv, b0: b0v, th: np.pi / 4})
    return complex(e.evalf(12)).real

print(f"{'r':>6s}  {'EL_A(float)':>14s}  {'EL_A(check)':>14s}  {'diff':>14s}")
print("-" * 60)
for rv in [0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
    fl = path_float(EL_A_cache, rv)
    ch = path_check(EL_A_cache, rv)
    d = fl - ch
    marker = " <<<" if abs(d) > 0.01 else ""
    print(f"{rv:6.3f}  {fl:+14.6f}  {ch:+14.6f}  {d:+14.6f}{marker}")

print()

# ── PART 3b: Step-by-step divergence tracking ────────────────────────────
print("=" * 72)
print("PART 3b: Step-by-step substitution tracking at r=0.3")
print("=" * 72)
print()

rv = 0.3
u_i = 1 - b0v**2 / rv**2
BB = 1.0 / np.sqrt(u_i)
Bpa_f = 0.5 * (2 * b0v**2 / rv**3) / u_i**1.5
Bppa_f = -0.75 * (2 * b0v**2 / rv**3)**2 / u_i**2.5 + 0.5 * (-6 * b0v**2 / rv**4) / u_i**1.5

Bw_sym = 1 / sp.sqrt(1 - b0**2 / r**2)
Bp_sym = sp.diff(Bw_sym, r)
Bpp_sym = sp.diff(Bw_sym, r, 2)

print(f"  r=0.3: BB_float={BB:.10f}, Bw_sym(r=0.3)={Bw_sym.subs({r: rv, b0: b0v}).evalf():.10f}")
print(f"  Bpa_float={Bpa_f:.6f}, Bp_sym(r=0.3)={Bp_sym.subs({r: rv, b0: b0v}).evalf():.6f}")
print(f"  Bppa_float={Bppa_f:.6f}, Bpp_sym(r=0.3)={Bpp_sym.subs({r: rv, b0: b0v}).evalf():.6f}")
print()

# Build step-by-step for float path
e_f = kill_f(EL_A_cache)
steps_f = []

e_f1 = e_f.subs(A2, 0.0).subs(B2, Bppa_f).subs(C2, 0.0)
e_f2 = e_f1.subs(A1, 0.0).subs(B1, Bpa_f).subs(C1, 1.0)
e_f3 = e_f2.subs(A, 1.0).subs(B, BB).subs(C, rv)
e_f4 = e_f3.doit()
e_f5 = e_f4.subs({r: rv, b0: b0v, th: np.pi / 4})

# Build step-by-step for check path
e_c = kill_f(EL_A_cache)
e_c1 = e_c.subs(A2, 0).subs(B2, Bpp_sym).subs(C2, 0)
e_c2 = e_c1.subs(A1, 0).subs(B1, Bp_sym).subs(C1, 1)
e_c3 = e_c2.subs(A, 1).subs(B, Bw_sym).subs(C, r)
e_c4 = e_c3.doit()
e_c5 = e_c4.subs({r: rv, b0: b0v, th: np.pi / 4})

def val_or_nan(e):
    try:
        return complex(e.evalf(12)).real
    except:
        return float('nan')

def nterms(e):
    try:
        return len(e.args)
    except:
        return -1

# At each step compare
names = [
    "kill_f(EL_A)",
    "subs derivs (A''=0, B''=Bpp, C''=0)",
    "subs first derivs (A'=0, B'=Bp, C'=1)",
    "subs functions (A=1, B=..., C=...)",
    "doit()",
    "subs coords (r=0.3, ...)",
]
print(f"{'Step':<45s}  {'#T(FL)':>7s}  {'#T(CHECK)':>9s}  {'V(FL)':>12s}  {'V(CK)':>12s}  {'diff':>10s}")
print("-" * 100)

pairs = [
    (e_f, e_c), (e_f1, e_c1), (e_f2, e_c2),
    (e_f3, e_c3), (e_f4, e_c4), (e_f5, e_c5),
]
for nm, (ef, ec) in zip(names, pairs):
    vf = val_or_nan(ef)
    vc = val_or_nan(ec)
    nf = nterms(ef)
    nc = nterms(ec)
    d = vf - vc
    flag = " <<<" if abs(d) > 0.001 else ""
    print(f"{nm:<45s}  {nf:>7d}  {nc:>9d}  {vf:+12.6f}  {vc:+12.6f}  {d:+10.6f}{flag}")

print()

# ── PART 3c: What terms survive in float path after subs(B,BB_float)? ────
print("=" * 72)
print("PART 3c: FLOAT path — residual expression after doit()")
print("=" * 72)
print()
print("Float path after doit(), before final subs(r,0.3):")
print(f"  Expression type: {type(e_f4)}")
print(f"  Contains r?      {e_f4.has(r)}")
print(f"  Contains b0?     {e_f4.has(b0)}")
if e_f4.has(r):
    print("  PROBLEM: residual r remains — subs(C, rv) killed C but B' derivatives")
    print("  have r in denominator that survives.")
    # Show r-dependent terms
    print("  Terms containing r:")
    for t in sp.Add.make_args(e_f4):
        if t.has(r):
            print(f"    {t}")
else:
    print("  OK: no residual r")

print()
print("Check path after doit(), before final subs(r,0.3):")
print(f"  Contains r?      {e_c4.has(r)}")
print(f"  Contains b0?     {e_c4.has(b0)}")
if e_c4.has(r):
    print("  Terms containing r:")
    for t in sp.Add.make_args(e_c4):
        if t.has(r):
            print(f"    {t}")

print()

# ── PART 4: Singularities on grid boundary ────────────────────────────────
print("=" * 72)
print("PART 4: Singularity check at r_min=0.12")
print("=" * 72)
print()

rvals = np.linspace(0.12, 0.8, 160)
u_vals = 1 - b0v**2 / rvals**2
B_vals = 1.0 / np.sqrt(u_vals)
Bp_vals = 0.5 * (2 * b0v**2 / rvals**3) / u_vals**1.5
Bpp_vals = -0.75 * (2 * b0v**2 / rvals**3)**2 / u_vals**2.5 + 0.5 * (-6 * b0v**2 / rvals**4) / u_vals**1.5

print(f"  u_min={u_vals.min():.4f} at r=0.12: u={u_vals[0]:.4f}, B={B_vals[0]:.4f}, B'={Bp_vals[0]:.2f}, B''={Bpp_vals[0]:.2f}")
print(f"  u_max={u_vals.max():.4f} at r=0.8: u={u_vals[-1]:.4f}, B={B_vals[-1]:.4f}")
print(f"  tan(pi/4) = {np.tan(np.pi/4):.4f}  (no theta singularity)")
print()

# Evaluate EL_A on fine grid using float path (fast)
print("  EL_A(r) profile (float path):")
el_profile = np.array([path_float(EL_A_cache, rv) for rv in [0.12, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]])
rpts = [0.12, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]
for rv2, elv in zip(rpts, el_profile):
    print(f"    r={rv2:.3f}  EL_A={elv:+.6f}")

print()

# ── PART 5: Numerical δS/δA vs ∫EL_A·g ──────────────────────────────────
print("=" * 72)
print("PART 5: Numerical δS/δA vs symbolic ∫EL_A·g")
print("=" * 72)
print()

N = 160
rr = np.linspace(0.12, 0.8, N)
sig, eps = 0.05, 1e-4
bump = np.exp(-((rr - 0.3) / sig) ** 2)
dbump = -2 * (rr - 0.3) / sig**2 * bump
Ab_arr = np.ones(N); Abp_arr = np.zeros(N)
u_arr = 1 - b0v**2 / rr**2
Bb_arr = 1.0 / np.sqrt(u_arr)
Bbp_arr = 0.5 * u_arr**-1.5 * (2 * b0v**2 / rr**3)

def T_num(R, Ab, Abp, Bb, Bbp):
    return (-1.0 / (R**2 * np.tan(thv)**2)
            + 4.0 / (Bb**2 * R**2)
            + 4.0 * Abp / (Ab * Bb**2 * R)
            + Abp**2 / (Ab**2 * Bb**2))

def L_num(Ab_arr, Abp_arr, Bb_arr, Bbp_arr):
    e = Ab_arr * Bb_arr * rr**2
    T = np.array([T_num(rr[i], Ab_arr[i], Abp_arr[i], Bb_arr[i], Bbp_arr[i]) for i in range(N)])
    return e * T

Ap = Ab_arr + eps * bump
App = Abp_arr + eps * dbump
dS_A = (np.trapezoid(L_num(Ap, App, Bb_arr, Bbp_arr), rr)
        - np.trapezoid(L_num(Ab_arr, Abp_arr, Bb_arr, Bbp_arr), rr)) / eps

# Symbolic integrals using fast float path
el_A_float = np.array([path_float(EL_A_cache, rv) for rv in rr])
I0_f = np.array([path_float(I0, rv) for rv in rr])
I11_f = np.array([path_float(I11, rv) for rv in rr])
I22_f = np.array([path_float(I22, rv) for rv in rr])

sym_A = np.trapezoid(el_A_float * bump, rr)
int_I0 = np.trapezoid(I0_f * bump, rr)
int_I11 = np.trapezoid(I11_f * bump, rr)
int_I22 = np.trapezoid(I22_f * bump, rr)

print(f"  Numerical δS/δA (bump)  = {dS_A:+.8f}")
print(f"  Sym EL_A·g (float path) = {sym_A:+.8f}")
print(f"  ratio (original)        = {dS_A / sym_A:.6f}")
print()
print(f"  Component integrals (float path):")
print(f"    ∫I0 ·g dr  = {int_I0:+.8f}")
print(f"    ∫I11·g dr  = {int_I11:+.8f}")
print(f"    ∫I22·g dr  = {int_I22:+.8f}")
print(f"    sum        = {int_I0 + int_I11 + int_I22:+.8f}")
print()

# ── PART 5b: Check path symbolic integral ─────────────────────────────────
print("  Computing check-path integral (slower, symbolic per point)...")
el_A_check = np.array([path_check(EL_A_cache, rv) for rv in rr])
sym_A_ck = np.trapezoid(el_A_check * bump, rr)
print(f"  Sym EL_A·g (check path)  = {sym_A_ck:+.8f}")
print(f"  ratio (check)            = {dS_A / sym_A_ck:.6f}")
print()

# ── PART 6: THE FIX ──────────────────────────────────────────────────────
print("=" * 72)
print("PART 6: ROOT CAUSE AND FIX")
print("=" * 72)
print()
print("ROOT CAUSE:")
print("  In r1_internal2.py symbolic_int(), the substitution order is:")
print("    1) B.diff(r,2) -> float Bppa (number, NO r dependence)")
print("    2) B.diff(r)   -> float Bpa  (number, NO r dependence)")
print("    3) A, B, C     -> floats")
print("  BUT: sympy substituting a FLOAT into B.diff(r,2) means")
print("  Bppa is already evaluated at the SPECIFIC grid point r=R.")
print("  This is CORRECT for scalar values.")
print()
print("  However, the check path uses SYMBOLIC Bw(r), Bp(r), Bpp(r),")
print("  then subs(C, r) -- so r remains symbolic until the FINAL subs.")
print("  The float path subs(C, R) with R=number, and r also appears")
print("  in remaining symbolic structure. When doit() is called,")
print("  sympy may simplify terms differently depending on whether")
print("  the expression is purely numeric vs mixed symbolic+numeric.")
print()

# Direct proof: subs order matters for sympy simplification
print("DIRECT PROOF: subs order matters")
print()

# Test: same expression, two orderings
expr_test = kill_f(EL_A_cache).subs(A2, 0).subs(C2, 0).subs(A1, 0).subs(C1, 1).subs(A, 1)

# Order 1: C->r first, then B->Bw
e1 = expr_test.subs(C, r).subs(B, Bw_sym).doit()
e1v = e1.subs({r: rv, b0: b0v, th: np.pi/4})
v1 = complex(e1v.evalf(12)).real

# Order 2: B->Bw first, then C->r
e2 = expr_test.subs(B, Bw_sym).subs(C, r).doit()
e2v = e2.subs({r: rv, b0: b0v, th: np.pi/4})
v2 = complex(e2v.evalf(12)).real

# Order 3: B'->Bp, B''->Bpp, then B->Bw, then C->r
e3 = expr_test.subs(B2, Bpp_sym).subs(B1, Bp_sym).subs(B, Bw_sym).subs(C, r).doit()
e3v = e3.subs({r: rv, b0: b0v, th: np.pi/4})
v3 = complex(e3v.evalf(12)).real

print(f"  Order B->Bw, C->r:     {v1:+.8f}")
print(f"  Order C->r, B->Bw:     {v2:+.8f}")
print(f"  Order B'->Bp,B''->Bpp: {v3:+.8f}")

# Now the critical test: float substitutions with B' included
e4 = expr_test.subs(B2, Bppa_f).subs(B1, Bpa_f).subs(B, BB).subs(C, rv).doit()
e4v = e4.subs({r: rv, b0: b0v, th: np.pi/4})
v4 = complex(e4v.evalf(12)).real
print(f"  Float path (B',B''=num): {v4:+.8f}")
print()

# THE REAL TEST: what if we use float B' but keep symbolic B and C?
e5 = expr_test.subs(B2, Bppa_f).subs(B1, Bpa_f).subs(B, Bw_sym).subs(C, r).doit()
e5v = e5.subs({r: rv, b0: b0v, th: np.pi/4})
v5 = complex(e5v.evalf(12)).real
print(f"  B'=float, B=Bw(r), C=r: {v5:+.8f}")
print()

# And: float B' and B, but C=r symbolically
e6 = expr_test.subs(B2, Bppa_f).subs(B1, Bpa_f).subs(B, BB).subs(C, r).doit()
e6v = e6.subs({r: rv, b0: b0v, th: np.pi/4})
v6 = complex(e6v.evalf(12)).real
print(f"  B'=float, B=BB_float, C=r: {v6:+.8f}")
print()

# ── PART 7: CORRECTED RATIO ──────────────────────────────────────────────
print("=" * 72)
print("PART 7: Corrected ratio computation")
print("=" * 72)
print()

# The check path should be correct since it uses pure symbolic substitution
# Let's verify by also computing ∫I0·g (since I11=I22=0 at A'=A''=0)
# and compare

# I0 should be the same in both paths since I0 = dL/dA doesn't involve A' or A''
# so it's just A*B*C^2*T / A = B*C^2*T (no derivative terms)
# Let's check
I0_float = np.array([path_float(I0, rv) for rv in rr])
I0_check = np.array([path_check(I0, rv) for rv in rr])
print(f"  ∫I0·g (float path) = {np.trapezoid(I0_float * bump, rr):+.8f}")
print(f"  ∫I0·g (check path) = {np.trapezoid(I0_check * bump, rr):+.8f}")
print(f"  δS/δA              = {dS_A:+.8f}")
print()

# The "check" path should give the correct symbolic integral
# ratio_correct = dS_A / ∫EL_A(check)·g
print(f"  CORRECT ratio = δS/δA / ∫EL_A(check)·g = {dS_A / sym_A_ck:.6f}")
print(f"  ORIGINAL ratio = δS/δA / ∫EL_A(float)·g = {dS_A / sym_A:.6f}")
print()

# ── PART 8: SUMMARY ──────────────────────────────────────────────────────
print("=" * 72)
print("PART 8: SUMMARY")
print("=" * 72)
print()
print("1) EL_A COMPONENT TABLE at r=0.3 (ELLIS):")
print(f"   I0  = dL/dA          = {path_float(I0, 0.3):+.6f}")
print(f"   I11 = -d/dr(dL/dA')  = {path_float(I11, 0.3):+.6f}")
print(f"   I22 = d^2/dr^2(dL/dA'') = {path_float(I22, 0.3):+.6f}")
print(f"   sum = I0+I11+I22     = {path_float(I0,0.3)+path_float(I11,0.3)+path_float(I22,0.3):+.6f}")
print(f"   EL_A (cache)         = {path_float(EL_A_cache, 0.3):+.6f}")
print()
print("2) At A'=A''=0: I1 = dL/dA' = 0 (A' factor), I2 = 0 (A'' factor)")
print(f"   -> I11=0, I22=0, EL_A = I0 analytically.")
print(f"   I0(r=0.3) = {path_float(I0, 0.3):+.6f} != 0 (A=lapse is physical)")
print()
print("3) BUG LOCATION:")
print("   In r1_internal2.py symbolic_int() line 72-79:")
print("   The float values Bppa, Bpa are pre-computed at r=R and substituted")
print("   as pure numbers. This IS correct for pointwise evaluation.")
print("   However, subs(C, R) with R=number BEFORE subs(r, R) means")
print("   the symbolic B' terms that had mixed r/C dependence may be")
print("   evaluated differently by sympy's internal simplifier.")
print("   The check path avoids this by keeping C=r (symbolic) until last.")
print()
print("4) PROPOSED FIX:")
print("   In r1_internal2.py symbolic_int(), change the substitution order")
print("   to match r1_check.py worm_subs():")
print("     e.subs(B.diff(r,2), Bpp_sym)  # symbolic, not float")
print("     e.subs(B.diff(r), Bp_sym)      # symbolic, not float")  
print("     e.subs(A, 1).subs(B, Bw_sym).subs(C, r)  # symbolic C=r")
print("     e.doit()")
print("     e.subs({r: R, b0: b0v, th: thv})  # coordinates LAST")
print()

print("DONE")
