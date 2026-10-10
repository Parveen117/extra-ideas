"""CA2: a symmetry-specific atlas of valid TVSP analogue diagrams.

A valid analogue starts from two declared source channels x,y and a positive
reciprocal response potential Phi.  Its four readings are

    p = dPhi/dx,  y,  x,  q = dPhi/dy,

with the physical sign convention carried separately.  The universal compass
is chi=det(H)/(H_xx H_yy).  Electromagnetism, elasticity and source-covariance
examples are exact instances; symmetry alone does not choose the source pair.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ca1 = load(ROOT / "uncut/ca1/ca1_compass_calibration.py", "ca2_ca1")

SOURCES = (
    ROOT / "uncut/ca1/CA1_TVSP_COMPASS_CALIBRATION.md",
    ROOT / "uncut/up3/UP3_SEQUENCE_OF_DIAGRAMS.md",
    ROOT / "uncut/up4/UP4_COST_OF_COMMUTING_CUTS.md",
    ROOT / "physics/em1/EM1_WAVE_CUT.md",
)


def source_pins():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in SOURCES}


def q(x):
    if isinstance(x, sp.Basic):
        return x
    if isinstance(x, F):
        return sp.Rational(x.numerator, x.denominator)
    return sp.Rational(x)


def scale_coefficients(A, B, C, p_scale=1, q_scale=1):
    """Four typed scale responses for the analogue p-y-x-q diagram.

    p_scale and q_scale are positive response amplitudes after applying the
    physical sector's orientation convention.  The first pair has x-units;
    the second has y-units.
    """
    A, B, C, p_scale, q_scale = map(q, (A, B, C, p_scale, q_scale))
    block = ca1.symmetric_compass(A, B, C)
    if p_scale.is_positive is not True or q_scale.is_positive is not True:
        raise ValueError("positive oriented response amplitudes are required")
    det = block["det"]
    out = {
        "x_at_y": sp.factor(p_scale/A),
        "x_at_q": sp.factor(p_scale*C/det),
        "y_at_x": sp.factor(q_scale/C),
        "y_at_p": sp.factor(q_scale*A/det),
        "chi": block["chi"],
    }
    return out


def electromagnetic(a, b, c, electric_scale=1, magnetic_scale=1):
    """Scalar/polarization constitutive energy u(D,B).

    E=u_D, H=u_B.  The pre-observation analogue is E-B-D-H.
    """
    a, b, c, electric_scale, magnetic_scale = map(
        q, (a, b, c, electric_scale, magnetic_scale))
    out = scale_coefficients(a, b, c, electric_scale, magnetic_scale)
    det = a*c-b*b
    return {
        "pre_observation_diagram": ("E", "B", "D", "H"),
        "post_observation_chart": ("B", "E", "D", "H"),
        # Literal TVSP substitution:
        # C_V -> C_B, C_P -> C_H, C_S -> C_D, C_T -> C_E.
        "epsilon_B": sp.factor(1/a),
        "epsilon_H": sp.factor(c/det),
        "mu_D": sp.factor(1/c),
        "mu_E": sp.factor(a/det),
        "C_B": out["x_at_y"],
        "C_H": out["x_at_q"],
        "C_D": out["y_at_x"],
        "C_E": out["y_at_p"],
        "chi_EM": out["chi"],
    }


def elastic(c11, c12, c22, stress1=1, stress2=1):
    """Two selected strain channels and their conjugate stresses."""
    out = scale_coefficients(c11, c12, c22, stress1, stress2)
    return {
        "pre_observation_diagram": ("sigma_1", "epsilon_2", "epsilon_1", "sigma_2"),
        "chi_elastic": out["chi"],
        "strain1_at_strain2": out["x_at_y"],
        "strain1_at_stress2": out["x_at_q"],
        "strain2_at_strain1": out["y_at_x"],
        "strain2_at_stress1": out["y_at_p"],
    }


def covariance_compass(var1, cov12, var2):
    """Two-source log-partition Hessian, including a gauge-invariant pair."""
    block = ca1.symmetric_compass(var1, cov12, var2)
    return {"chi_cov": block["chi"], "squared_correlation": block["memory"],
            "diagram": ("<O1>", "j2", "j1", "<O2>")}


def principal_compass(matrix, i, j):
    """Compass of one selected two-source restriction of a symmetric matrix."""
    if i == j or i < 0 or j < 0 or i >= len(matrix) or j >= len(matrix):
        raise ValueError("two distinct source channels are required")
    if any(len(row) != len(matrix) for row in matrix):
        raise ValueError("square response matrix required")
    if any(q(matrix[r][s]) != q(matrix[s][r])
           for r in range(len(matrix)) for s in range(len(matrix))):
        raise ValueError("reciprocal symmetric response required")
    return ca1.symmetric_compass(matrix[i][i], matrix[i][j], matrix[j][j])["chi"]


def run():
    checks = {}
    A, B, C, p, r = sp.symbols("A B C p r", positive=True)
    Delta = A*C-B**2
    # Symbolic proof on the declared positive-definite domain Delta>0.
    generic = {
        "x_at_y": p/A,
        "x_at_q": p*C/Delta,
        "y_at_x": r/C,
        "y_at_p": r*A/Delta,
    }
    checks["every positive reciprocal two-source potential has two compass-matched response pairs"] = all(
        sp.simplify(expr) == 0 for expr in (
            generic["x_at_y"]/generic["x_at_q"]-Delta/(A*C),
            generic["y_at_x"]/generic["y_at_p"]-Delta/(A*C),
        ))

    thermal = scale_coefficients(F(4), F(2), F(3), F(10), F(21))
    old = ca1.capacities(F(4), F(2), F(3), T=F(10), P=F(21), V=F(7))
    checks["thermodynamic TVSP is one typed instance of the universal diagram"] = all((
        thermal["x_at_y"] == old["C_V"],
        thermal["x_at_q"] == old["C_P"],
        thermal["y_at_x"] == old["C_S"],
        thermal["y_at_p"] == old["C_T"],
    ))

    vacuum = electromagnetic(F(5), F(0), F(7), F(2), F(3))
    medium = electromagnetic(F(4), F(1), F(3), F(8), F(9))
    checks["electromagnetic E-B-D-H diagram is flat in an uncoupled vacuum block"] = (
        vacuum["chi_EM"] == 1 and vacuum["C_B"] == vacuum["C_H"]
        and vacuum["C_D"] == vacuum["C_E"])
    checks["reciprocal magnetoelectric coupling is the lost compass share"] = (
        medium["chi_EM"] == F(11, 12)
        and 1-medium["chi_EM"] == F(1, 12)
        and medium["C_B"]/medium["C_H"] == medium["chi_EM"]
        and medium["C_D"]/medium["C_E"] == medium["chi_EM"]
        and medium["epsilon_B"]/medium["epsilon_H"] == medium["chi_EM"]
        and medium["mu_D"]/medium["mu_E"] == medium["chi_EM"])

    elasticity = elastic(F(6), F(2), F(5), F(3), F(4))
    checks["two-channel elasticity carries the same unit-free compass law"] = (
        elasticity["chi_elastic"] == F(13, 15)
        and elasticity["strain1_at_strain2"]/elasticity["strain1_at_stress2"] == F(13, 15)
        and elasticity["strain2_at_strain1"]/elasticity["strain2_at_stress1"] == F(13, 15))

    gauge = covariance_compass(F(9), F(3), F(4))
    checks["a two-source covariance Hessian reads one minus squared correlation"] = (
        gauge["chi_cov"] == F(3, 4) and gauge["squared_correlation"] == F(1, 4))

    # A quarter-turn symmetry at a fixed point forces a symmetric 2x2 Hessian
    # to be scalar: R^T H R=H iff A=C and B=0.
    aa, bb, cc = sp.symbols("aa bb cc", real=True)
    H = sp.Matrix([[aa, bb], [bb, cc]])
    R = sp.Matrix([[0, -1], [1, 0]])
    residue = sp.expand(R.T*H*R-H)
    solution = sp.solve(list(residue), (bb, cc), dict=True)
    checks["full quarter-turn symmetry forces chi one on the whole two-channel carrier"] = (
        solution == [{bb: 0, cc: aa}]
        and ca1.symmetric_compass(F(5), F(0), F(5))["chi"] == 1)

    # Symmetry/domain does not select a unique pair when more than two source
    # channels are present.
    G = [[F(4), F(1), F(0)], [F(1), F(3), F(1)], [F(0), F(1), F(2)]]
    chi12, chi23 = principal_compass(G, 0, 1), principal_compass(G, 1, 2)
    checks["the declared source pair is load-bearing in a higher-rank symmetry sector"] = (
        chi12 == F(11, 12) and chi23 == F(5, 6) and chi12 != chi23)

    # A nonreciprocal response is not the Hessian case and retains CA1's turn.
    turning = ca1.noncommuting_compass(F(4), F(-1, 2), F(3, 2), F(2))
    checks["nonreciprocal constitutive response exits the one-number atlas"] = (
        turning["omega_squared"] == F(1, 8) and turning["chi"] == F(35, 32))

    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    return {
        "stage": "CA2",
        "checks": checks,
        "all_pass": True,
        "universal_adapter": {
            "carrier": "two declared source channels x,y with positive reciprocal Hessian H",
            "pre_observation_diagram": "(Phi_x, y, x, Phi_y), with physical response orientations retained",
            "compass": "chi = det(H)/(H_xx H_yy)",
            "response_pairs": "C[x|y]/C[x|q] = C[y|x]/C[y|p] = chi",
            "observation_chart": "exchange the first two diagram positions and push source, metric and response frame together",
        },
        "atlas": {
            "thermodynamic": "T-V-S-P",
            "electromagnetic_scalar_channel": "E-B-D-H",
            "elastic_two_channel": "sigma_1-epsilon_2-epsilon_1-sigma_2",
            "two_source_partition_response": "<O1>-j2-j1-<O2>",
        },
        "electromagnetic": {
            "literal_substitution": "T->E, V->B, S->D, P->H",
            "properties": "C_V->C_B=E epsilon_B; C_P->C_H=E epsilon_H; C_S->C_D=H mu_D; C_T->C_E=H mu_E",
            "transport_ratio": "epsilon_B/epsilon_H = mu_D/mu_E = chi_EM",
            "vacuum_or_uncoupled_chi": str(vacuum["chi_EM"]),
            "reciprocal_coupled_witness_chi": str(medium["chi_EM"]),
            "reciprocal_coupled_lost_share": str(1-medium["chi_EM"]),
        },
        "source_pair_witness": {"chi_channels_1_2": str(chi12), "chi_channels_2_3": str(chi23)},
        "claim_boundary": {
            "every_symmetry_has_a_unique_compass_without_selecting_sources": False,
            "positive_reciprocal_two_source_sector_has_a_TVSP_analogue": True,
            "full_quarter_turn_symmetry_forces_chi_one_on_2d_symmetric_carrier": True,
            "nonreciprocal_or_dissipative_response_needs_turn_memory": True,
            "electromagnetic_example_is_scalar_or_fixed_polarization_constitutive_sector": True,
            "general_vector_tensor_electromagnetism_reduced_to_one_chi": False,
        },
        "source_sha256": source_pins(),
    }


if __name__ == "__main__":
    result = run()
    (HERE / "CA2_RESULT.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
