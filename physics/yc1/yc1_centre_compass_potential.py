"""YC1: the TVSP response compass of the SU(2) centre record.

The carrier is one Haar-distributed SU(2) face with

    c = Tr(U)/2 in [-1,1],    <c^(2n)> = Catalan(n)/4^n.

The two-source potential is

    Psi(kappa,h) = log <exp(kappa*c + h*c^2)>.

Its Hessian is the covariance matrix of (c,c^2).  Centre reflection sends
c to -c and fixes c^2.  The exact theorem and its Yang--Mills claim boundary
are in YC1_CENTRE_COMPASS_POTENTIAL.md.  This executable certificate uses
only rational formal series for the local identities and mpmath for
independent numerical replays of the Bessel formulas.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

SOURCES = (
    ROOT / "uncut/ca2/CA2_SYMMETRY_COMPASS_ATLAS.md",
    ROOT / "uncut/up1/UP1_DEGREE_OF_THE_POTENTIAL.md",
    ROOT / "physics/ct1/CT1_CENTRE_RECORD_OF_THE_TURN_CHAIN.md",
    ROOT / "tools/rw1/RW1_RETURNED_WINDING_COUNT.md",
    ROOT / "physics/tc1/TC1_CORE_OF_TURNS.md",
    ROOT / "physics/cm2/CM2_COMPACT_CORE_MATCHING.md",
)


def source_pins():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in SOURCES}


def catalan(n):
    return F(comb(2*n, n), n + 1)


def haar_moment(n):
    """Exact moment of c=Tr(U)/2 under normalized SU(2) Haar measure."""
    if n < 0:
        raise ValueError("moment order must be nonnegative")
    if n % 2:
        return F(0)
    j = n // 2
    return catalan(j) / 4**j


def sadd(a, b, sign=1, order=None):
    if order is None:
        order = max(len(a), len(b)) - 1
    return [(a[n] if n < len(a) else F(0))
            + sign*(b[n] if n < len(b) else F(0))
            for n in range(order + 1)]


def smul(a, b, order):
    out = [F(0)]*(order + 1)
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if i + j > order:
                break
            out[i + j] += x*y
    return out


def sdiv(a, b, order):
    """Formal series a/b through the requested order, with b[0] nonzero."""
    if not b or b[0] == 0:
        raise ValueError("series denominator needs a nonzero constant term")
    out = [F(0)]*(order + 1)
    for n in range(order + 1):
        value = a[n] if n < len(a) else F(0)
        value -= sum(b[j]*out[n-j] for j in range(1, min(n, len(b)-1) + 1))
        out[n] = value/b[0]
    return out


def derivative(a):
    return [(n + 1)*a[n + 1] for n in range(len(a)-1)]


def integral(a, order):
    out = [F(0)]*(order + 1)
    for n in range(min(len(a), order)):
        out[n + 1] = a[n]/(n + 1)
    return out


def tilted_raw_moment_series(power, order):
    """Series of <c^power exp(kappa c)>/<exp(kappa c)> at h=0."""
    z = [haar_moment(n)/factorial(n) for n in range(order + 1)]
    numerator = [haar_moment(n + power)/factorial(n) for n in range(order + 1)]
    return sdiv(numerator, z, order)


def response_series(order=12):
    """Exact kappa series for Psi, ST, centre cut, Hessian and compass."""
    if order < 4:
        raise ValueError("order at least four is required")
    moments = [tilted_raw_moment_series(j, order) for j in range(5)]
    r = moments[1]
    psi = integral(r, order)
    st = [F(n)*psi[n] for n in range(order + 1)]
    centre = [psi[n]-st[n]/2 for n in range(order + 1)]

    var_c = sadd(moments[2], smul(moments[1], moments[1], order), -1, order)
    var_c2 = sadd(moments[4], smul(moments[2], moments[2], order), -1, order)
    cov = sadd(moments[3], smul(moments[1], moments[2], order), -1, order)
    lost = sdiv(smul(cov, cov, order), smul(var_c, var_c2, order), order)
    chi = [-x for x in lost]
    chi[0] += 1
    return {
        "psi": psi,
        "r": r,
        "st": st,
        "centre": centre,
        "var_c": var_c,
        "var_c2": var_c2,
        "cov": cov,
        "lost": lost,
        "chi": chi,
    }


def centre_hessian():
    """Hessian at the centre-symmetric Haar point (kappa,h)=(0,0)."""
    mean_c, mean_c2 = haar_moment(1), haar_moment(2)
    var_c = haar_moment(2)-mean_c**2
    cov = haar_moment(3)-mean_c*mean_c2
    var_c2 = haar_moment(4)-mean_c2**2
    det = var_c*var_c2-cov**2
    return ((var_c, cov), (cov, var_c2)), det/(var_c*var_c2)


def face_potential(kappa):
    """Psi(kappa,0)=log(2 I_1(kappa)/kappa), continuously at zero."""
    kappa = mp.mpf(kappa)
    if not kappa:
        return mp.mpf(0)
    return mp.log(2*mp.besseli(1, kappa)/kappa)


def face_response(kappa):
    """Psi_kappa=I_2/I_1=tanh(k_c), continuously at zero."""
    kappa = mp.mpf(kappa)
    if not kappa:
        return mp.mpf(0)
    return mp.besseli(2, kappa)/mp.besseli(1, kappa)


def centre_cut(kappa):
    kappa = mp.mpf(kappa)
    return face_potential(kappa)-kappa*face_response(kappa)/2


def as_sparse(series):
    return {str(n): str(x) for n, x in enumerate(series) if x}


def run():
    checks = {}
    series = response_series(12)

    checks["SU(2) Haar moments reproduce the centre-record Catalan law"] = (
        [haar_moment(n) for n in range(9)]
        == [F(1), F(0), F(1, 4), F(0), F(1, 8), F(0), F(5, 64), F(0), F(7, 128)])

    expected_psi = {2: F(1, 8), 4: F(-1, 384), 6: F(1, 9216),
                    8: F(-1, 184320)}
    expected_r = {1: F(1, 4), 3: F(-1, 96), 5: F(1, 1536),
                  7: F(-1, 23040)}
    expected_centre = {4: F(1, 384), 6: F(-1, 4608),
                       8: F(1, 61440)}
    checks["one-face log potential and response have the exact weak series"] = (
        all(series["psi"][n] == value for n, value in expected_psi.items())
        and all(series["r"][n] == value for n, value in expected_r.items()))
    checks["the TVSP diagram centre kills the quadratic term and first retains degree four"] = (
        series["centre"][0:4] == [F(0)]*4
        and all(series["centre"][n] == value for n, value in expected_centre.items()))
    checks["the ST scalar is kappa times the centre response"] = all(
        series["st"][n] == (series["r"][n-1] if n else 0)
        for n in range(len(series["st"])))

    hessian, chi0 = centre_hessian()
    checks["centre flip diagonalizes the c and c-squared response Hessian"] = (
        hessian == ((F(1, 4), F(0)), (F(0), F(1, 16))) and chi0 == 1)
    checks["the first compass information loss away from the centre is kappa squared over four"] = (
        series["lost"][0:4] == [F(0), F(0), F(1, 4), F(0)]
        and series["chi"][0:5] == [F(1), F(0), F(-1, 4), F(0), F(1, 16)])

    # Exact formal replay of r'=1-r^2-3r/kappa.  Multiplication by kappa
    # removes the removable singularity at zero.
    r, dr = series["r"], derivative(series["r"])
    order = len(dr)-1
    lhs = [F(0)] + dr
    rhs = sadd([F(0), F(1)], [F(0)] + smul(r, r, order), -1, order + 1)
    rhs = sadd(rhs, [3*x for x in r[:order + 2]], -1, order + 1)
    checks["the Bessel response obeys its exact Riccati identity as a formal series"] = (
        lhs[:order + 2] == rhs[:order + 2])

    mp.mp.dps = 60
    density = lambda c: 2*mp.sqrt(1-c*c)/mp.pi
    integral_ok = True
    derivative_ok = True
    rw1_ok = True
    for kappa in (mp.mpf("0.125"), mp.mpf("0.75"), mp.mpf("2"), mp.mpf("9")):
        z = mp.quad(lambda c: mp.exp(kappa*c)*density(c), [-1, 0, 1])
        mean = mp.quad(lambda c: c*mp.exp(kappa*c)*density(c), [-1, 0, 1])/z
        integral_ok &= abs(mp.log(z)-face_potential(kappa)) < mp.mpf("1e-45")
        integral_ok &= abs(mean-face_response(kappa)) < mp.mpf("1e-45")
        r0 = face_response(kappa)
        rp = 1-r0*r0-3*r0/kappa
        derivative_ok &= abs(mp.diff(face_response, kappa)-rp) < mp.mpf("1e-45")
        sinh2kc = 2*r0/(1-r0*r0)
        rw1_ok &= kappa/2 < sinh2kc < 2*kappa/3
        mid_prime = (4*r0-kappa*(1-r0*r0))/2
        rw1_ok &= 0 < mid_prime < r0/2
    checks["Haar integral equals the Bessel potential and response"] = bool(integral_ok)
    checks["Bessel differentiation gives r-prime equals one minus r-squared minus 3r over kappa"] = bool(derivative_ok)
    checks["RW1 bounds make the centre-cut scalar strictly increasing"] = bool(rw1_ok)

    scalar_sandwich_ok = True
    for kappa in (mp.mpf("0.01"), mp.mpf("0.5"), mp.mpf("3"), mp.mpf("40")):
        psi, st, mid = face_potential(kappa), kappa*face_response(kappa), centre_cut(kappa)
        scalar_sandwich_ok &= 0 < mid < psi/2 and psi < st < 2*psi
    checks["the centre scalar lies below half the potential and sandwiches ST"] = bool(
        scalar_sandwich_ok)

    stored_rw1 = json.loads((ROOT / "tools/rw1/RW1_RESULT.json").read_text())
    checks["the pinned RW1 global inequality certificate is present and passing"] = bool(
        stored_rw1["checks"]["R6_centre_record_between_kappa_over_2_and_2kappa_over_3"])

    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])

    result = {
        "stage": "YC1",
        "checks": checks,
        "all_pass": True,
        "carrier": {
            "group": "one SU(2) face with normalized Haar measure",
            "centre_coordinate": "c = Tr(U)/2",
            "centre_action": "c -> -c",
            "observables": ["c", "c^2"],
            "sources": ["kappa", "h"],
            "potential": "Psi(kappa,h) = log <exp(kappa c + h c^2)>",
            "metric": "Hessian Psi = Cov(c,c^2)",
        },
        "tvsp": {
            "pre_observation": "(T,V,S,P) = (Psi_kappa,h,kappa,Psi_h)",
            "post_observation": "(V,T,S,P) = (h,Psi_kappa,kappa,Psi_h)",
            "st_scalar_at_h_zero": "kappa Psi_kappa = kappa I2(kappa)/I1(kappa) = kappa tanh(k_c)",
            "diagram_centre": "Psi_mid = Psi-(kappa Psi_kappa+h Psi_h)/2",
        },
        "centre_seam": {
            "hessian_at_origin": [[str(x) for x in row] for row in hessian],
            "chi_for_every_h_at_kappa_zero": "1",
            "reason": "centre parity makes Cov(c,c^2)=0",
        },
        "weak_series_at_h_zero": {
            name: as_sparse(series[name]) for name in ("psi", "r", "st", "centre", "lost", "chi")
        },
        "global_centre_cut_bound": {
            "input": "RW1: kappa/2 < sinh(2 k_c) < 2 kappa/3 for kappa>0",
            "derivative": "Psi_mid' = (4r-kappa(1-r^2))/2",
            "conclusion": "0 < Psi_mid' < r/2, hence 0<Psi_mid<Psi/2 and Psi<ST<2Psi for kappa>0",
        },
        "claim_boundary": {
            "gauge_invariant_two_source_response_potential_constructed": True,
            "centre_symmetry_forces_complete_compass_at_the_seam": True,
            "diagram_centre_removes_the_quadratic_face_cumulant": True,
            "first_retained_face_term_and_TC1_core_share_quartic_degree": True,
            "one_face_kappa_identified_with_CM2_theta_or_continuum_coupling": False,
            "one_face_quartic_coefficient_identified_with_TC1_e2_core": False,
            "mass_gap_or_volume_uniformity_proved": False,
        },
        "source_sha256": source_pins(),
    }
    return result


if __name__ == "__main__":
    result = run()
    (HERE / "YC1_RESULT.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
