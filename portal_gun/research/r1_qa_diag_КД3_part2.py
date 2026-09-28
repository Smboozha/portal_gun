"""
R1 q=A DIAGNOSTIC [КД3] — PART 2: Deep dive into kill_f failure and fix
"""
import sympy as sp
import numpy as np
import sys, os, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symsub import r, th, b0, A, B, C, F, kill_f

data = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ft_el_cache.json'), encoding='utf-8'))
ns = {'r': r, 'th': th, 'b0': b0, 'A': A, 'B': B, 'C': C, 'f': F}
EL_A_cache = sp.sympify(data['EL_A'], locals=ns)

b0v = 0.1; thv = np.pi/4
rv = 0.3

# ── TEST 1: What does kill_f ACTUALLY do to Subs(Derivative(f,...))? ─────
print("=" * 72)
print("TEST 1: kill_f diagnostic on EL_A_cache")
print("=" * 72)
print()

# Count terms before/after kill_f
print(f"EL_A_cache: {len(EL_A_cache.args)} top-level Add args")
el_kf = kill_f(EL_A_cache)
print(f"kill_f(EL_A_cache): {len(el_kf.args)} top-level Add args")
print()

# Check: does kill_f fully remove all f-dependent structure?
print("Does kill_f(EL_A_cache) still contain f?", el_kf.has(F))
print("Does kill_f(EL_A_cache) still contain Subs?", el_kf.has(sp.Subs))
print()

# Find Subs terms
print("Subs terms in EL_A_cache:")
subs_terms = [a for a in EL_A_cache.args if a.has(sp.Subs)]
print(f"  {len(subs_terms)} terms with Subs")
for i, st in enumerate(subs_terms[:3]):
    print(f"  term {i}: {str(st)[:200]}")
print()

print("Subs terms in kill_f(EL_A_cache):")
subs_terms_kf = [a for a in el_kf.args if a.has(sp.Subs)]
print(f"  {len(subs_terms_kf)} terms with Subs")
for i, st in enumerate(subs_terms_kf[:3]):
    print(f"  term {i}: {str(st)[:200]}")
print()

# ── TEST 2: Compare kill_f vs manual Subs removal ────────────────────────
print("=" * 72)
print("TEST 2: kill_f vs explicit Subs(Derivative(f,...)) = 0")
print("=" * 72)
print()

# Method A: kill_f
eA = kill_f(EL_A_cache)
eA = eA.subs(A.diff(r,2), 0).subs(C.diff(r,2), 0).subs(A.diff(r), 0).subs(C.diff(r), 1).subs(A, 1)
eA = eA.subs(B.diff(r,2), sp.diff(1/sp.sqrt(1-b0**2/r**2), r, 2))
eA = eA.subs(B.diff(r), sp.diff(1/sp.sqrt(1-b0**2/r**2), r))
eA = eA.subs(B, 1/sp.sqrt(1-b0**2/r**2))
eA = eA.subs(C, r).doit()
eA = eA.subs({r: rv, b0: b0v, th: thv})
vA = complex(eA.evalf(12)).real
print(f"  Method A (kill_f then check-path subs): {vA:+.8f}")

# Method B: Remove Subs explicitly, then do the same
eB = EL_A_cache
# First: substitute f->0 in Subs wrapper AND in the inner expression
eB = eB.replace(sp.Subs, lambda expr, *args: sp.S.Zero)
eB = eB.subs(F, 0)
eB = eB.subs(A.diff(r,2), 0).subs(C.diff(r,2), 0).subs(A.diff(r), 0).subs(C.diff(r), 1).subs(A, 1)
eB = eB.subs(B.diff(r,2), sp.diff(1/sp.sqrt(1-b0**2/r**2), r, 2))
eB = eB.subs(B.diff(r), sp.diff(1/sp.sqrt(1-b0**2/r**2), r))
eB = eB.subs(B, 1/sp.sqrt(1-b0**2/r**2))
eB = eB.subs(C, r).doit()
eB = eB.subs({r: rv, b0: b0v, th: thv})
vB = complex(eB.evalf(12)).real
print(f"  Method B (Subs->0 explicitly):             {vB:+.8f}")

# Method C: From scratch (EL_A_rebuilt = I0+I11+I22)
T_scalar = sp.sympify(data['T_scalar'], locals=ns)
L = A * B * C**2 * T_scalar
I0  = sp.diff(L, A)
I1  = sp.diff(L, A.diff(r))
I11 = -sp.diff(I1, r)
I22 = sp.diff(sp.diff(L, A.diff(r,2)), r, 2)
EL_rebuilt = I0 + I11 + I22

eC = EL_rebuilt
eC = eC.subs(A.diff(r,2), 0).subs(C.diff(r,2), 0).subs(A.diff(r), 0).subs(C.diff(r), 1).subs(A, 1)
eC = eC.subs(B.diff(r,2), sp.diff(1/sp.sqrt(1-b0**2/r**2), r, 2))
eC = eC.subs(B.diff(r), sp.diff(1/sp.sqrt(1-b0**2/r**2), r))
eC = eC.subs(B, 1/sp.sqrt(1-b0**2/r**2))
eC = eC.subs(C, r).doit()
eC = eC.subs({r: rv, b0: b0v, th: thv})
vC = complex(eC.evalf(12)).real
print(f"  Method C (rebuilt EL_A = I0+I11+I22):     {vC:+.8f}")
print()

# Numerical reference
print(f"  Numerical δS/δA                          = -0.138079")
print()

# ── TEST 3: Where exactly does kill_f fail? ──────────────────────────────
print("=" * 72)
print("TEST 3: Isolating the Subs(Derivative(f,...)) terms")
print("=" * 72)
print()

# Extract terms with Subs and without
no_subs = [a for a in EL_A_cache.args if not a.has(sp.Subs)]
with_subs = [a for a in EL_A_cache.args if a.has(sp.Subs)]

print(f"  Terms WITHOUT Subs: {len(no_subs)}")
print(f"  Terms WITH Subs:    {len(with_subs)}")

# Evaluate each group
e_no_subs = sp.Add(*no_subs)
e_with_subs = sp.Add(*with_subs)

# Check-path eval of each
def eval_check(expr):
    Bw = 1/sp.sqrt(1-b0**2/r**2)
    Bp = sp.diff(Bw, r)
    Bpp = sp.diff(Bw, r, 2)
    e = expr
    e = e.subs(A.diff(r,2), 0).subs(B.diff(r,2), Bpp).subs(C.diff(r,2), 0)
    e = e.subs(A.diff(r), 0).subs(B.diff(r), Bp).subs(C.diff(r), 1)
    e = e.subs(A, 1).subs(B, Bw).subs(C, r)
    e = e.doit()
    e = e.subs({r: rv, b0: b0v, th: thv})
    return complex(e.evalf(12)).real

v_no_subs = eval_check(e_no_subs)
v_with_subs = eval_check(e_with_subs)
print(f"  Value of no-Subs terms (check path): {v_no_subs:+.8f}")
print(f"  Value of with-Subs terms (check path): {v_with_subs:+.8f}")
print(f"  Sum: {v_no_subs + v_with_subs:+.8f}")
print()

# After kill_f: what happens to with_subs terms?
e_with_subs_kf = kill_f(e_with_subs)
print(f"  kill_f(with-Subs terms): has Subs? {e_with_subs_kf.has(sp.Subs)}")
print(f"  kill_f(with-Subs terms): has f?    {e_with_subs_kf.has(F)}")
v_kf_subs = eval_check(e_with_subs_kf)
print(f"  Value of kill_f(with-Subs) (check path): {v_kf_subs:+.8f}")
print()

# ── TEST 4: The actual Subs structure ────────────────────────────────────
print("=" * 72)
print("TEST 4: Inspecting Subs(Derivative(f,...)) structure")
print("=" * 72)
print()

for i, st in enumerate(with_subs[:2]):
    print(f"  with_subs term {i}:")
    print(f"    {st}")
    print(f"    free_symbols: {st.free_symbols}")
    print()
    # Apply kill_f to this specific term
    st_kf = kill_f(st)
    print(f"    After kill_f: {st_kf}")
    print(f"    After kill_f: has F? {st_kf.has(F)}")
    print(f"    After kill_f: is zero? {st_kf == 0}")
    print()

# ── TEST 5: Numerical profile comparison ─────────────────────────────────
print("=" * 72)
print("TEST 5: EL_A profiles — rebuilt vs cache(kill_f) vs cache(no kill_f)")
print("=" * 72)
print()

rpts = [0.12, 0.15, 0.2, 0.3, 0.5, 0.8]

def eval_rebuilt(rv):
    e = EL_rebuilt
    e = e.subs(A.diff(r,2), 0).subs(C.diff(r,2), 0).subs(A.diff(r), 0).subs(C.diff(r), 1).subs(A, 1)
    e = e.subs(B.diff(r,2), sp.diff(1/sp.sqrt(1-b0**2/r**2), r, 2))
    e = e.subs(B.diff(r), sp.diff(1/sp.sqrt(1-b0**2/r**2), r))
    e = e.subs(B, 1/sp.sqrt(1-b0**2/r**2))
    e = e.subs(C, r).doit()
    e = e.subs({r: rv, b0: b0v, th: thv})
    return complex(e.evalf(12)).real

def eval_cache_kf(rv):
    e = kill_f(EL_A_cache)
    e = e.subs(A.diff(r,2), 0).subs(C.diff(r,2), 0).subs(A.diff(r), 0).subs(C.diff(r), 1).subs(A, 1)
    e = e.subs(B.diff(r,2), sp.diff(1/sp.sqrt(1-b0**2/r**2), r, 2))
    e = e.subs(B.diff(r), sp.diff(1/sp.sqrt(1-b0**2/r**2), r))
    e = e.subs(B, 1/sp.sqrt(1-b0**2/r**2))
    e = e.subs(C, r).doit()
    e = e.subs({r: rv, b0: b0v, th: thv})
    return complex(e.evalf(12)).real

print(f"{'r':>6s}  {'rebuilt':>14s}  {'cache(kf)':>14s}  {'diff':>10s}")
print("-" * 50)
for rv in rpts:
    vr = eval_rebuilt(rv)
    vc = eval_cache_kf(rv)
    print(f"{rv:6.3f}  {vr:+14.8f}  {vc:+14.8f}  {vr-vc:+10.6f}")

print()
print("=" * 72)
print("CONCLUSION")
print("=" * 72)
print()
print("The bug is in kill_f: it replaces f -> 0 but does NOT remove")
print("Subs(Derivative(f,...)) wrappers. After kill_f, the Subs terms")
print("with f=0 INSIDE still contain geometric prefactors and can")
print("evaluate to non-zero values.")
print()
print("FIX: Replace kill_f with a stronger version that also removes")
print("Subs(Derivative(f,...)) entirely:")
print()
print("  def kill_f_strong(expr):")
print("      return expr.replace(")
print("          lambda t: isinstance(t, sp.Subs) and t.has(F),")
print("          lambda t: sp.S.Zero)")
print()
print("OR: Rebuild EL_A from L = A*B*C^2*T (using I0+I11+I22)")
print("instead of relying on the cached EL_A with kill_f.")
print()
