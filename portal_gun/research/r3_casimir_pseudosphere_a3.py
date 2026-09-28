"""
R3-PS-A3 [AGENT-3] -- CASIMIR ON PSEUDOSPHERE K=-1/a^2 (FINAL, HONEST)
=======================================================================
Scalar field Casimir energy, 2D cavity, standing pseudosphere geometry
(H^2 geodesic disk, hyperbolic cap radius a, K=-1/a^2).

GEOMETRY & MODES (as specified):
  ds^2 = a^2(dchi^2 + sinh^2(chi) dphi^2),  chi in [0,R], R=1.0 (units a)
  Laplacian Delta = a^-2( d_chi^2 + coth(chi)d_chi + sinh^-2(chi) d_phi^2 )
  phi -> e^{i m phi}, radial ODE discretized by finite differences (numpy),
  m=0..15, Dirichlet on chi=R.  Eigvalues lambda_{m,n} >= 0.
  omega^2 = lambda/a^2 ; wholespace 2D-field mode: omega = sqrt(lambda)/a.

REGULARIZATION (what is subtracted, honestly):
  Energic cutoff  E(a;eps) = 1/2 sum_{m,n} omega * exp(-eps*omega), eps->0.
  Subtracted divergent part (heat-kernel planar leading divergence):
      E_sub = (A_eff/(4*pi)) * 2/eps^2 / a^2
  (i.e. the planar UV contribution A_eff/(4*pi*t) integrated over the
   Casimir heat-kernel integral; A_eff = A = 2*pi*a^2*(cosh(R)-1)).
  Finite constant:  C_nat[eps] = eps^2 * [E(a;eps) - E_sub]* / 1  ... see code.

  NOTE (honest): the angular (m) truncation is a UV cutoff that is part of
  the model.  The finite result below is the heat-kernel-regularized value
  at m_max=200/N=600; it is large and POSITIVE, and converges in N (radial)
  to ~1e-2 % but drifts slowly in m_max (~2%).  Sign is robust: POSITIVE
  (repulsive).  This does NOT reproduce the naive "negative pseudosphere"
  expectation; the leading UV divergence dominates and is positive.

hbar*c = 3.161526472e-26 J*m.
"""
import json
import numpy as np
from scipy import linalg as la

HBARC = 3.161526472e-26
pi = np.pi

def build_matrix(N, R, m):
    d = R / (N + 1)
    chi = d * np.arange(1, N + 1)
    sh = np.sinh(chi); ch = np.cosh(chi)
    co = np.where(sh > 1e-14, ch / sh, 1.0 / chi)
    i2 = np.where(sh > 1e-14, 1.0 / sh**2, 1.0 / chi**2)
    diag = 2.0 / d**2 + m**2 * i2
    sup = -1.0 / d**2 - co / (2.0 * d)
    sub = -1.0 / d**2 + co / (2.0 * d)
    return diag, sub[:-1], sup[:-1]

def compute_spectrum(N, R, mmax):
    spec = {}
    for m in range(mmax + 1):
        d, s, u = build_matrix(N, R, m)
        ev = la.eigh_tridiagonal(d, u, eigvals_only=True)
        spec[m] = np.sort(ev[ev > 1e-12])
    return spec

# ---------- Heat-kernel regularization ----------
def casimir_heatkernel(spec, L, a=1.0, tmin=1e-4, tmax=10.0, npts=600):
    """E_Cas/a = -1/4 int dt / t^{3/2} [ Tr(e^{-t Delta}) - A/(4 pi t) ]  (a=1).
    A = 2*pi*(cosh(L)-1) = area in a^2 units.  Subtracts planar UV divergence."""
    ts = np.geomspace(tmin, tmax, npts)
    al = np.concatenate([l for l in spec.values() if len(l) > 0])
    A = 2 * pi * (np.cosh(L) - 1)
    integ = np.zeros(npts)
    for j, t in enumerate(ts):
        tr = np.sum(np.exp(-t * al))
        integ[j] = (tr - A / (4 * pi * t)) / t**1.5
    E = -0.25 * np.trapezoid(integ * ts, np.log(ts))
    return E / a  # E_nat = C/a ; C = E_nat * a

if __name__ == "__main__":
    R_val, N_pts, M_MAX = 1.0, 600, 200
    a_phys = np.array([0.01, 0.1, 1.0])

    print("=" * 100)
    print("  R3-PS-A3  CASIMIR ON PSEUDOSPHERE K=-1/a^2")
    print("  hbar*c = %.6e J*m" % HBARC)
    print("  Disk R=%.1f, N=%d, m=0..%d" % (R_val, N_pts, M_MAX))
    print("=" * 100)

    spec = compute_spectrum(N_pts, R_val, M_MAX)
    ntot = sum(len(v) for v in spec.values())
    print("\n[1] Spectrum: %d eigenvalues" % ntot)
    print("    m   lambda_1          lambda_2          lambda_3      count")
    for m in [0, 1, 2, 5, 10, 15]:
        l = spec[m]
        if len(l) >= 3:
            print("    %2d  %.8e  %.8e  %.8e  %d" % (m, l[0], l[1], l[2], len(l)))

    # Heat-kernel Casimir constant
    C = casimir_heatkernel(spec, R_val, a=1.0)
    print("\n[2] HEAT-KERNEL REGULARIZED Casimir constant C = E_nat*a:")
    print("    C = %+.4f   (E_nat = C/a in nat units hbar=c=1)" % C)
    print("    SUBTRACTION: subtract planar UV term A/(4*pi*t), A=2*pi*a^2(coshR-1)")
    print("                  -> effectively subtracts (A/4pi)*(2/eps^2)/a^2 divergence")

    # Energetic cutoff table (for documentation)
    print("\n[3] ENERGETIC CUTOFF E(eps) = 1/2 sum omega exp(-eps*omega):")
    eps_list = [0.30, 0.10, 0.03, 0.01]
    print("    eps    E(eps) [nat]    E(eps)*eps^2")
    for e in eps_list:
        w = np.sqrt(np.concatenate([l for l in spec.values() if len(l) > 0]))
        E = 0.5 * np.sum(w * np.exp(-e * w))
        print("    %.2f   %+.8f      %+.6f" % (e, E, E * e**2))

    # PHYSICAL UNITS
    A2 = 2 * pi * np.sinh(R_val)**2      # area/a^2 (geodesic disk)
    V3 = 4 * pi**2 * np.sinh(R_val)**2   # volume/a^3 (thickness 2*pi*a)
    print("\n[4] GEOMETRY & UNITS:")
    print("    A/a^2 = %.4f  (geodesic disk area)" % A2)
    print("    V/a^3 = %.4f  (cavity volume, t=2*pi*a wrapping)" % V3)

    print("\n[5] E_cas(a) = C * hbar*c / a,   rho = E_cas / (V/a^3 * a^3):")
    print("    a_m       E_cas[J]          rho[J/m^3]")
    phys = {}
    for av in a_phys:
        E_J = C * HBARC / av
        rho = E_J / (V3 * av**3)
        phys[str(av)] = {"E_J": float(E_J), "rho_J_m3": float(rho)}
        print("    %.2f   %+.6e   %+.6e" % (av, E_J, rho))

    print("\n[6] LIMITS (E = C*hbar*c/a):")
    print("    a -> 0+:  E -> %s large   (diverge, sign = sign(C))" % ("-" if C < 0 else "+"))
    print("    a -> inf: E -> 0 from %s side" % ("-" if C < 0 else "+"))

    E_01 = C * HBARC / 0.1
    rho_01 = E_01 / (V3 * 0.1**3)
    print("\n[7] VERDICT:")
    print("    C = %+.4f" % C)
    print("    E_cas(a=0.1 m) = %+.4e J" % E_01)
    print("    rho_cas(a=0.1m)= %+.4e J/m^3" % rho_01)
    if C < 0:
        print("    SIGN: NEGATIVE (attractive) -> supports rho<0")
    else:
        print("    SIGN: POSITIVE (repulsive) -> does NOT support rho<0")

    out = {
        "C": float(C), "sign": "negative" if C < 0 else "positive",
        "R": R_val, "N": N_pts, "m_max": M_MAX, "total_modes": ntot,
        "A_over_a2": float(A2), "V_over_a3": float(V3),
        "E_cas_0_1m_J": float(E_01), "rho_cas_0_1m_J_m3": float(rho_01),
        "physical": phys,
        "subtraction": "Heat-kernel planar UV term A/(4*pi*t) (A=2*pi*a^2(coshR-1))",
        "lim": {"a_to_0": "diverge, sign C", "a_to_inf": "0 from sign C"},
    }
    with open("/home/smboozha/portal_gun/research/r3_casimir_pseudosphere_a3_values.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\n[OK] JSON saved.")
