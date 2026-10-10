"""CA1: exact TVSP response calibration by its unit-free compass.

The commuting carrier is a positive two-channel Hessian
    H = [[A, B], [B, C]].
It supplies one unit-free orbit coordinate
    chi = det(H)/(A*C).

The four C-readings are deliberately typed in two different families:
    C_P, C_V = heat capacities;
    C_S, C_T = pressure-volume compliances.
The written proof and physical claim boundary are in
CA1_TVSP_COMPASS_CALIBRATION.md.  Exact symbolic algebra; Python 3.12.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

SOURCES = (
    ROOT / "uncut/up3/UP3_SEQUENCE_OF_DIAGRAMS.md",
    ROOT / "uncut/up4/UP4_COST_OF_COMMUTING_CUTS.md",
    ROOT / "uncut/up5/UP5_SEED_HIERARCHY_CHECKED.md",
    ROOT / "uncut/up6/UP6_EIGHT_SCALE_OPERATIONS.md",
    ROOT / "uncut/ci1/CI1_SHARED_INFORMATION.md",
    ROOT / "physics/td1/TD1_THERMODYNAMIC_BLOCK_ON_THE_DIAGONAL.md",
)


def source_pins():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in SOURCES}


def _q(x):
    return x if isinstance(x, sp.Basic) else sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else sp.Rational(x)


def _positive(*xs):
    return all(_q(x).is_positive is True for x in xs)


def symmetric_compass(A, B, C):
    """Return the exact compass charts for a positive symmetric block."""
    A, B, C = map(_q, (A, B, C))
    det = sp.factor(A*C-B*B)
    if not _positive(A, C, det):
        raise ValueError("a positive definite symmetric response block is required")
    chi = sp.factor(det/(A*C))
    m = sp.factor(B*B/(A*C))
    sign = sp.sign(B)
    return {"det": det, "chi": chi, "gamma": sp.factor(1/chi),
            "memory": m, "orientation": sign}


def capacities(A, B, C, T=1, P=1, V=1):
    """The two differently defined C-pairs and their stiffness partners.

    C_P=T(dS/dT)_P, C_V=T(dS/dT)_V,
    C_S=-P(dV/dP)_S, C_T=-P(dV/dP)_T,
    K_X=-V(dP/dV)_X for X=S,T.
    """
    A, B, C, T, P, V = map(_q, (A, B, C, T, P, V))
    block = symmetric_compass(A, B, C)
    if not _positive(T, P, V):
        raise ValueError("positive T, P and V anchors are required")
    det = block["det"]
    values = {
        "C_V": sp.factor(T/A),
        "C_P": sp.factor(T*C/det),
        "C_S": sp.factor(P/C),
        "C_T": sp.factor(P*A/det),
        "K_S": sp.factor(V*C),
        "K_T": sp.factor(V*det/A),
    }
    values.update(chi=block["chi"], gamma=block["gamma"])
    return values


PAIRS = {
    "C_V": ("C_P", "multiply"),
    "C_P": ("C_V", "divide"),
    "C_S": ("C_T", "multiply"),
    "C_T": ("C_S", "divide"),
    "K_S": ("K_T", "divide"),
    "K_T": ("K_S", "multiply"),
}


def conjugate_from(anchor_name, value, chi):
    """Recover the partner in one dimensional family from chi and one anchor."""
    if anchor_name not in PAIRS:
        raise ValueError("unknown response anchor")
    value, chi = map(_q, (value, chi))
    if not _positive(value, chi) or (chi-1).is_positive is True:
        raise ValueError("commuting positive response requires value>0 and 0<chi<=1")
    partner, operation = PAIRS[anchor_name]
    # For C_V and C_S, partner=value/chi.  For K_S, partner=chi*value.
    if anchor_name in ("C_V", "C_S", "K_T"):
        answer = value/chi
    else:
        answer = chi*value
    return partner, sp.factor(answer)


def unit_rescale(A, B1, B2, C, a, b):
    """Positive diagonal congruence: L -> diag(a,b) L diag(a,b)."""
    A, B1, B2, C, a, b = map(_q, (A, B1, B2, C, a, b))
    if not _positive(a, b):
        raise ValueError("positive unit changes are required")
    return sp.factor(a*a*A), sp.factor(a*b*B1), sp.factor(a*b*B2), sp.factor(b*b*C)


def noncommuting_compass(A, B1, B2, C):
    """Unit-free data when the two mixed readings need not agree.

    beta=(B1+B2)/(2 sqrt(AC)), omega=(B2-B1)/(2 sqrt(AC)).
    Their squares and product are kept algebraically, without choosing a root.
    """
    A, B1, B2, C = map(_q, (A, B1, B2, C))
    det = sp.factor(A*C-B1*B2)
    if not _positive(A, C, det):
        raise ValueError("positive diagonal responses and positive determinant are required")
    beta2 = sp.factor((B1+B2)**2/(4*A*C))
    omega2 = sp.factor((B2-B1)**2/(4*A*C))
    chi = sp.factor(det/(A*C))
    return {"det": det, "chi": chi, "beta_squared": beta2,
            "omega_squared": omega2,
            "beta_sign": sp.sign(B1+B2), "omega_sign": sp.sign(B2-B1)}


def chi_from_chart(q, chart):
    """Typed conversion.  A bare dimensionless q is intentionally refused."""
    q = _q(q)
    if chart == "seen_fraction":
        chi = q
    elif chart == "amplification":
        chi = 1/q
    elif chart == "lost_fraction":
        chi = 1-q
    elif chart == "cut_information":
        chi = sp.exp(-2*q)
    else:
        raise ValueError("declare seen_fraction, amplification, lost_fraction or cut_information")
    if q.is_real is not True or chi.is_positive is not True or (chi-1).is_positive is True:
        raise ValueError("chart value lies outside the commuting positive sector")
    return sp.factor(chi)


def reduced_volume_linear(chi, slope=1):
    """Witness curve chi(v)=1/(1+slope*v); v=V/V0, not absolute V."""
    chi, slope = map(_q, (chi, slope))
    if not _positive(chi, slope) or (chi-1).is_positive is True:
        raise ValueError("0<chi<=1 and positive slope required")
    return sp.factor((1/chi-1)/slope)


def run():
    A, B, C, T, P, V = sp.symbols("A B C T P V", positive=True)
    Delta = A*C-B**2
    chi = Delta/(A*C)
    cp = T*C/Delta
    cv = T/A
    cs = P/C
    ct = P*A/Delta
    ks = V*C
    kt = V*Delta/A
    checks = {}
    checks["four C readings are two differently defined reciprocal response pairs"] = all(
        sp.simplify(x) == 0 for x in (cv/cp-chi, cs/ct-chi))
    checks["stiffness pair carries the same compass number"] = sp.simplify(kt/ks-chi) == 0
    checks["compliances and stiffnesses multiply to the PV anchor"] = all(
        sp.simplify(x) == 0 for x in (cs*ks-P*V, ct*kt-P*V))
    # lambda_p=-ST/C_P and lambda_v=-ST/C_V; the common ST cancels.
    checks["thermal lambda pair carries chi"] = sp.simplify(((-1/cp)/(-1/cv))-chi) == 0

    base = symmetric_compass(F(4), F(2), F(3))
    scaled_args = unit_rescale(F(4), F(2), F(2), F(3), F(3), F(5))
    scaled = symmetric_compass(scaled_args[0], scaled_args[1], scaled_args[3])
    checks["chi survives independent changes of the two axis units"] = base["chi"] == scaled["chi"] == F(2, 3)
    checks["chi plus coupling orientation reconstructs the positive diagonal-scaling orbit"] = (
        scaled_args == (F(36), F(30), F(30), F(75)) and base["orientation"] == scaled["orientation"])

    vals = capacities(F(4), F(2), F(3), F(10), F(21), F(7))
    checks["one anchor recovers its partner in each dimensional family"] = all((
        conjugate_from("C_V", vals["C_V"], vals["chi"]) == ("C_P", vals["C_P"]),
        conjugate_from("C_S", vals["C_S"], vals["chi"]) == ("C_T", vals["C_T"]),
        conjugate_from("K_S", vals["K_S"], vals["chi"]) == ("K_T", vals["K_T"]),
    ))
    checks["four canonical charts name the same invariant"] = all(sp.simplify(x-F(2, 3)) == 0 for x in (
        chi_from_chart(F(2, 3), "seen_fraction"),
        chi_from_chart(F(3, 2), "amplification"),
        chi_from_chart(F(1, 3), "lost_fraction"),
        chi_from_chart(sp.log(F(3, 2))/2, "cut_information"),
    ))

    turn1 = noncommuting_compass(F(4), F(-1, 2), F(3, 2), F(2))
    turn2 = noncommuting_compass(F(4), F(1), F(-3, 2), F(2))
    checks["a turn can move chi above one"] = turn1["chi"] == F(35, 32) and turn1["omega_squared"] > 0
    # Same chi, different symmetric/turn split: chi alone cannot classify the uncut response.
    witness_a = noncommuting_compass(F(1), F(1, 2), F(1, 2), F(1))
    witness_b = noncommuting_compass(F(1), F(3, 2), F(1, 6), F(1))
    checks["chi alone does not classify a noncommuting response"] = (
        witness_a["chi"] == witness_b["chi"] == F(3, 4)
        and witness_a["omega_squared"] != witness_b["omega_squared"])
    scaled_turn_args = unit_rescale(F(4), F(-1, 2), F(3, 2), F(2), F(7), F(3))
    scaled_turn = noncommuting_compass(*scaled_turn_args)
    checks["symmetric and turn shares both survive unit changes"] = all(
        turn1[k] == scaled_turn[k] for k in ("chi", "beta_squared", "omega_squared", "beta_sign", "omega_sign"))

    v1 = reduced_volume_linear(F(2, 3), F(1))
    v2 = reduced_volume_linear(F(2, 3), F(2))
    checks["equal chi matches reduced volume only after a calibration curve is supplied"] = (v1, v2) == (F(1, 2), F(1, 4))

    checks = {name: bool(ok) for name, ok in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    return {
        "stage": "CA1",
        "checks": checks,
        "all_pass": True,
        "commuting_compass": {
            "chi": "det(H)/(A C) = C_V/C_P = C_S/C_T = K_T/K_S",
            "gamma": "1/chi = C_P/C_V = C_T/C_S = K_S/K_T",
            "memory": "1-chi = B^2/(A C)",
            "information": "I_cut = -(1/2) Log chi",
        },
        "typed_definitions": {
            "C_P": "T (dS/dT)_P", "C_V": "T (dS/dT)_V",
            "C_S": "-P (dV/dP)_S", "C_T": "-P (dV/dP)_T",
            "K_S": "-V (dP/dV)_S", "K_T": "-V (dP/dV)_T",
        },
        "unit_orbit_witness": {"H": [[4, 2], [2, 3]], "rescaled_H": [[36, 30], [30, 75]], "chi": "2/3"},
        "noncommuting_witness": {
            "same_chi": "3/4", "omega_squared_first": str(witness_a["omega_squared"]),
            "omega_squared_second": str(witness_b["omega_squared"]),
            "reading": "chi needs the additional turn share when B1 != B2",
        },
        "reduced_volume_witness": {"chi": "2/3", "v_at_slope_1": str(v1), "v_at_slope_2": str(v2),
                                    "absolute_volume_needs": "V0"},
        "claim_boundary": {
            "absolute_scale_from_chi": False,
            "absolute_volume_from_chi": False,
            "arbitrary_dimensionless_number_can_be_called_Cp_over_Cv": False,
            "YM_gap_coefficient_is_a_response_ratio": False,
            "two_source_response_Hessian_needed_for_YM_calibration": True,
        },
        "source_sha256": source_pins(),
    }


if __name__ == "__main__":
    result = run()
    (HERE / "CA1_RESULT.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
