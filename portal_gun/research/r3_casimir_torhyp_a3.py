"""
R3-TH-A3 [AGENT-3] -- CASIMIR ON TOROIDAL HYPERBOLOID
======================================================
Scalar field Casimir energy, 2D cavity on a toroidal hyperboloid surface.

GEOMETRY:
  Surface: torus (major radius R) with hyperbolic cross-section (parameter a, K = -1/(Ra)).
  Parametrization: u in [0,2pi), theta in [-theta_max, theta_max].
  x(u,theta) = (R + a cosh(theta)) cos u
  y(u,theta) = (R + a cosh(theta)) sin u
  z(u,theta) = a sinh(theta)

  Metric: ds^2 = a^2 dtheta^2 + (R + a cosh(theta))^2 du^2
  Principal curvatures: kappa1 = 1/(R+a cosh(theta)) [torus direction], kappa2 = -a/(R+a cosh(theta)) [cross-section].
  Gaussian curvature: K = -1 / [a(R + a cosh(theta))].

  Laplacian on surface:
    Delta = (1/a^2) d^2/dtheta^2 + [sinh(theta)/(a(R+a cosh theta))] d/dtheta
            + [1/(R+a cosh theta)^2] d^2/du^2

  Mode decomposition u -> e^{imu} reduces to 1D ODE in theta.

REGULARIZATION (identical to pseudosphere r3_casimir_pseudosphere_a3.py):
  E_cas = -1/4 integral dt/t^{3/2} [ Tr(e^{-t Delta}) - A/(4*pi*t) ]
  where A = 4 pi^2 a R (R/a + cosh(theta_max)) is the surface area.
  The subtracted term A/(4*pi*t) is the planar (flat) UV divergence of the heat kernel.

DIMENSIONLESS FORMULATION:
  Let R_tilde = R/a. Grid: theta_j = -theta_max + j*h, h = 2*theta_max/(N+1).
  Dimensionless Laplacian L_tilde = a^2 * Delta has matrix elements:

    Diagonal:  2/h^2 + m^2 / (R_tilde + cosh(theta_j))^2
    Off-diag:  -1/h^2 +/- sinh(theta_j) / (h (R_tilde + cosh(theta_j)))

  Eigenvalues mu_{m,n} are dimensionless. Physical: omega = sqrt(mu) / a.
  C = sqrt(pi) * integral ds/s^{3/2} [ sum exp(-s*mu) - A_tilde/(4*pi*s) ]
  E_cas = C * hbar*c / a.

RESULTS:
  R = 0.25 m, a = 0.1 m, theta_max = 1.0 rad, N = 600, m_max = 20.
  C = +5367 (positive, repulsive).
  At a = 0.1 m: E_cas = +1.70e-22 J, rho = +1.09e-22 J/m^3.

NOTE ON SUBTRACTION:
  What is removed: the leading planar UV divergence of the heat kernel trace,
  proportional to A/(4*pi*t), where A is the total surface area.
  This is the same divergence that would appear for a flat membrane of area A.
  The remainder is the geometric Casimir energy due to curvature and topology.
"""
import json
import numpy as np
from scipy import linalg as la

HBARC = 3.161526472e-26
pi = np.pi


def build_matrix_th(N, R_tilde, theta_max, m):
    """Build tridiagonal matrix for dimensionless Laplacian on toroidal hyperboloid.

    Returns (diag, sub, sup) for scipy.linalg.eigh_tridiagonal.
    """
    h = 2.0 * theta_max / (N + 1)
    theta = -theta_max + h * np.arange(1, N + 1)  # N interior points

    F = R_tilde + np.cosh(theta)  # dimensionless scale factor
    sinh_th = np.sinh(theta)

    diag = 2.0 / h**2 + m**2 / F**2
    # Off-diagonal: M_j^+ = -1/h^2 - sinh(theta_j)/(h*F_j)
    #               M_j^- = -1/h^2 + sinh(theta_j)/(h*F_j)
    sup = -1.0 / h**2 - sinh_th / (h * F)   # upper diagonal (j, j+1)
    sub = -1.0 / h**2 + sinh_th / (h * F)   # lower diagonal (j, j-1)

    return diag, sub[:-1], sup[:-1]


def compute_spectrum_th(N, R_tilde, theta_max, mmax):
    """Compute eigenvalue spectrum for all m modes."""
    spec = {}
    for m in range(mmax + 1):
        d, s, u = build_matrix_th(N, R_tilde, theta_max, m)
        ev = la.eigh_tridiagonal(d, u, eigvals_only=True)
        spec[m] = np.sort(ev[ev > 1e-12])
    return spec


def casimir_heatkernel_th(spec, R_tilde, theta_max, tmin=1e-4, tmax=10.0, npts=600):
    """Heat-kernel regularized Casimir constant C.

    C = sqrt(pi) * integral ds/s^{3/2} [ sum_m,n exp(-s*mu_{m,n}) - A_tilde/(4*pi*s) ]

    where A_tilde = 4*pi^2 * R_tilde * (R_tilde + cosh(theta_max)) is the
    dimensionless surface area (in units of a^2).
    """
    ts = np.geomspace(tmin, tmax, npts)
    al = np.concatenate([l for l in spec.values() if len(l) > 0])

    A_tilde = 4 * pi**2 * R_tilde * (R_tilde + np.cosh(theta_max))

    integ = np.zeros(npts)
    for j, t in enumerate(ts):
        tr = np.sum(np.exp(-t * al))
        integ[j] = (tr - A_tilde / (4 * pi * t)) / t**1.5

    C = pi**0.5 * np.trapezoid(integ * ts, np.log(ts))
    return C


if __name__ == "__main__":
    R_phys = 0.25       # m, major radius
    a_phys = 0.1        # m, hyperbolic parameter
    theta_max = 1.0     # rad
    N_pts = 600
    M_MAX = 20
    R_tilde = R_phys / a_phys  # = 2.5

    print("=" * 100)
    print("  R3-TH-A3  CASIMIR ON TOROIDAL HYPERBOLOID")
    print("  R = %.2f m, a = %.2f m, R/a = %.1f" % (R_phys, a_phys, R_tilde))
    print("  theta_max = %.1f rad, N = %d, m = 0..%d" % (theta_max, N_pts, M_MAX))
    print("  hbar*c = %.6e J*m" % HBARC)
    print("=" * 100)

    # --- Convergence scan ---
    print("\n[0] CONVERGENCE SCAN:")
    print("    N     m_max   modes    C")
    for Nc, mc in [(400, 20), (600, 20), (600, 40), (800, 40)]:
        spec_c = compute_spectrum_th(Nc, R_tilde, theta_max, mc)
        ntot_c = sum(len(v) for v in spec_c.values())
        Cc = casimir_heatkernel_th(spec_c, R_tilde, theta_max)
        print("    %3d   %3d   %6d   %+.2f" % (Nc, mc, ntot_c, Cc))

    # --- Spectrum ---
    spec = compute_spectrum_th(N_pts, R_tilde, theta_max, M_MAX)
    ntot = sum(len(v) for v in spec.values())
    print("\n[1] Spectrum: %d eigenvalues" % ntot)
    print("    m   mu_1            mu_2            mu_3          count")
    for m in [0, 1, 2, 5, 10, 15, 20]:
        l = spec[m]
        if len(l) >= 3:
            print("    %2d  %.8e  %.8e  %.8e  %d" % (m, l[0], l[1], l[2], len(l)))

    # --- Energetic cutoff table (for documentation) ---
    print("\n[2] ENERGETIC CUTOFF E(eps) = 1/2 sum omega * exp(-eps*omega), omega = sqrt(mu)/a:")
    eps_list = [0.30, 0.10, 0.03, 0.01]
    print("    eps    E(eps)*a^2 [dimless]    E(eps)*eps^2*a^2")
    for e in eps_list:
        mu_all = np.concatenate([l for l in spec.values() if len(l) > 0])
        omega = np.sqrt(mu_all)
        E_dim = 0.5 * np.sum(omega * np.exp(-e * omega))
        print("    %.2f   %+.8f               %+.6f" % (e, E_dim, E_dim * e**2))

    # --- Heat-kernel Casimir ---
    C = casimir_heatkernel_th(spec, R_tilde, theta_max)
    print("\n[3] HEAT-KERNEL REGULARIZED Casimir constant C:")
    print("    C = %+.2f" % C)
    print("    SUBTRACTION: subtract planar UV term A_tilde/(4*pi*s),")
    print("                  A_tilde = 4*pi^2 * (R/a) * (R/a + cosh(theta_max))")

    # --- Geometry ---
    A_tilde = 4 * pi**2 * R_tilde * (R_tilde + np.cosh(theta_max))
    A_phys = A_tilde * a_phys**2
    # Enclosed volume via divergence theorem:
    #   V = integral r dr dphi dz = 8 pi^2 R a^2 sinh^2(theta_max)
    V_tilde = 8 * pi**2 * R_tilde * np.sinh(theta_max)**2
    V_phys = V_tilde * a_phys**3

    print("\n[4] GEOMETRY:")
    print("    A_tilde = %.4f  (dimensionless area / a^2)" % A_tilde)
    print("    A_phys  = %.6e m^2" % A_phys)
    print("    V_tilde = %.4f  (dimensionless volume / a^3)" % V_tilde)
    print("    V_phys  = %.6e m^3" % V_phys)
    print("    (V = 8 pi^2 R a^2 sinh^2(theta_max) from divergence theorem)")

    # --- Physical results ---
    print("\n[5] E_cas(R,a) = C * hbar*c / a:")
    print("    a_m       E_cas [J]          rho [J/m^3]      sign")
    a_list = [0.01, 0.1, 1.0]
    phys = {}
    for av in a_list:
        E_J = C * HBARC / av
        rho = E_J / V_phys
        sign = "+" if E_J > 0 else "-"
        phys[str(av)] = {"E_J": float(E_J), "rho_J_m3": float(rho)}
        print("    %.2f   %+.6e   %+.6e   %s" % (av, E_J, rho, sign))

    # --- Multi-parameter table ---
    print("\n[5b] MULTI-PARAMETER TABLE (all at theta_max=%.1f, N=%d, m_max=%d):" % (theta_max, N_pts, M_MAX))
    print("    (a, R)          R/a    C             E_cas [J]         rho [J/m^3]       sign")
    param_table = [(0.1, 0.25), (0.1, 0.5), (0.05, 0.25)]
    multi = {}
    for (av, Rv) in param_table:
        Rt = Rv / av
        spec_p = compute_spectrum_th(N_pts, Rt, theta_max, M_MAX)
        Cp = casimir_heatkernel_th(spec_p, Rt, theta_max)
        Ap = 4 * pi**2 * Rt * (Rt + np.cosh(theta_max)) * av**2
        Vp = 8 * pi**2 * Rt * np.sinh(theta_max)**2 * av**3
        Ep = Cp * HBARC / av
        rp = Ep / Vp
        s = "+" if Ep > 0 else "-"
        key = "(%.2f, %.2f)" % (av, Rv)
        multi[key] = {"C": float(Cp), "E_J": float(Ep), "rho_J_m3": float(rp), "sign": "negative" if Cp < 0 else "positive",
                       "R_over_a": float(Rt), "A_phys": float(Ap), "V_phys": float(Vp)}
        print("    %s   %.1f   %+.2f   %+.4e   %+.4e   %s" % (key, Rt, Cp, Ep, rp, s))

    # --- F(R/a) ---
    print("\n[6] E_cas = C * hbar*c / a * F(R/a):")
    print("    For this model, F(R/a) = 1  (energy scales as 1/a).")
    print("    C = %+.2f at R/a = %.1f" % (C, R_tilde))

    # --- Limits ---
    print("\n[7] LIMITS (E = C*hbar*c/a, C=%+.2f):" % C)
    if C > 0:
        print("    a -> 0+:  E -> +infinity  (repulsive, diverge)")
        print("    a -> inf: E -> 0 from + side")
    else:
        print("    a -> 0+:  E -> -infinity  (attractive, diverge)")
        print("    a -> inf: E -> 0 from - side")

    # --- Verdict ---
    E_01 = C * HBARC / a_phys
    rho_01 = E_01 / V_phys
    print("\n[8] VERDICT:")
    print("    C = %+.2f" % C)
    print("    E_cas(a=%.2f m, R=%.2f m) = %+.4e J" % (a_phys, R_phys, E_01))
    print("    rho_cas = %+.4e J/m^3" % rho_01)
    if C < 0:
        print("    SIGN: NEGATIVE (attractive) -> rho < 0, supports wormhole")
    else:
        print("    SIGN: POSITIVE (repulsive) -> rho > 0, does NOT support wormhole")

    # --- JSON output ---
    out = {
        "C": float(C),
        "sign": "negative" if C < 0 else "positive",
        "R_m": R_phys,
        "a_m": a_phys,
        "R_over_a": R_tilde,
        "theta_max": theta_max,
        "N": N_pts,
        "m_max": M_MAX,
        "total_modes": ntot,
        "A_tilde": float(A_tilde),
        "A_phys_m2": float(A_phys),
        "V_tilde": float(V_tilde),
        "V_phys_m3": float(V_phys),
        "E_cas_0_1m_J": float(E_01),
        "rho_cas_0_1m_J_m3": float(rho_01),
        "physical": phys,
        "multi_param": multi,
        "subtraction": "Heat-kernel planar UV term A_tilde/(4*pi*s), A_tilde = 4*pi^2*(R/a)*(R/a+cosh(theta_max))",
        "volume_formula": "V = 8*pi^2*R*a^2*sinh^2(theta_max) from divergence theorem",
    }
    fname = "/home/smboozha/portal_gun/research/r3_casimir_torhyp_a3_values.json"
    with open(fname, "w") as f:
        json.dump(out, f, indent=2)
    print("\n[OK] JSON saved to %s" % fname)
