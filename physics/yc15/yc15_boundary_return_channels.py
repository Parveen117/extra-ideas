"""YC15: full interface-source variance and return into cube-only channels.

Exact identities and rational bounds support the accompanying analytic proofs.
The cube grounds and the returning inverses are not harmonically truncated.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


y14 = load('yc15_yc14', 'physics/yc14/yc14_observer_reset_and_interface.py')
y12 = y14.y12
GAP = y12.cube_bounds(1)['actual_full_gap_lower']
SOURCE_FLOOR = Q(13, 2)  # Conservative charged-source floor retained from YC14.
INVERSE_NORM_CAP = Q(27, 640)
ETA_MAX = Q(1, 4320)
RADIUS = Q(1, 128)
EXP_UPPER = Q(16, 15)


def quaternion_product(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


def conjugate(a):
    return (a[0], -a[1], -a[2], -a[3])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def wilson(a, b, u, v):
    # Two internal links a,b and two bridges u,v, with their traversal signs.
    return quaternion_product(quaternion_product(quaternion_product(a, u),
                                                conjugate(b)), conjugate(v))[0]


def bridge_gram(a, b, aa, bb):
    # Exact second Haar moments on each of the two independent bridges.
    axes = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    return sum(wilson(a, b, u, v)*wilson(aa, bb, u, v)
               for u in axes for v in axes)/sp.Integer(16)


def edge_kinetic_bound(theta):
    theta = Q(theta)
    if not 0 <= theta <= 1:
        raise ValueError('actual correlated cube window is 0<=theta<=1')
    return theta*theta/(96*(1-theta/2)**2)


def source_bounds(theta_a, theta_b):
    variance = edge_kinetic_bound(theta_a)+edge_kinetic_bound(theta_b)
    d = SOURCE_FLOOR
    norm_squared = Q(1, 576)+variance*(1/(864*d)+1/(144*d*d))
    return dict(theta_a=Q(theta_a), theta_b=Q(theta_b),
                source_norm_squared=Q(1, 4), source_mean_energy=Q(12),
                variance_upper=variance, susceptibility_lower=Q(1, 48),
                susceptibility_upper=Q(1, 48)+variance/(144*d),
                inverse_vector_norm_squared_upper=norm_squared,
                inverse_vector_norm_rational_cap=INVERSE_NORM_CAP,
                cube_return_difference_norm_squared_upper=variance/(576*d*d))


def join_bounds(eta):
    eta = Q(eta)
    if not 0 <= eta <= ETA_MAX:
        raise ValueError('certified interface window is 0<=eta<=1/4320')
    beta = 24*eta
    cube_seed = 24*4*INVERSE_NORM_CAP*eta
    bridge_seed = 2*(4*INVERSE_NORM_CAP+Q(1, 6))*eta
    seed = max(cube_seed, bridge_seed)
    contraction = 4*beta/GAP*EXP_UPPER*(2+8*(1+2*RADIUS))
    relative = 8*beta/GAP*EXP_UPPER
    mapping = seed+contraction*RADIUS
    return dict(eta=eta, factor_gap=GAP, nonlinear_minimum_support=1,
                radius=RADIUS, beta=beta, cube_seed=cube_seed,
                bridge_seed=bridge_seed, seed=seed,
                contraction=contraction, mapping=mapping,
                mapping_reserve=RADIUS-mapping,
                relative_return=relative,
                full_gap_lower=GAP*(1-relative),
                simple_gap_floor=Q(9, 40),
                fixed_point_norm_upper=seed/(1-contraction))


def first_channel_coefficient(incident=True):
    # Independent stationary perturbation derivation. The product Wp Wf
    # resolves at energies 18,26 if the internal face contains the edge;
    # at energy 24 otherwise. Subtract the change of the reference Omega/48.
    if incident:
        return Q(1, 6)*(Q(1, 16*18)+Q(3, 16*26))-Q(1, 576)
    return Q(1, 6)*Q(1, 4*24)-Q(1, 576)


@lru_cache(None)
def symbolic_checks():
    a, b, aa, bb = (sp.symbols(prefix+'0:4', real=True) for prefix in ('a', 'b', 'c', 'd'))
    gram = bridge_gram(a, b, aa, bb)
    checks = {
        'two bridge integrations retain both internal-link dot products':
            sp.expand(gram-dot(a, aa)*dot(b, bb)/4) == 0,
    }
    x, mu = sp.symbols('x mu', positive=True)
    checks['full inverse residual identity retains the quadratic remainder'] = sp.cancel(
        1/x-(1/mu-(x-mu)/mu**2+(x-mu)**2/(mu**2*x))) == 0
    checks['inverse squared residual identity retains both positive terms'] = sp.cancel(
        1/x**2-(1/mu**2-2*(x-mu)/mu**3+
                   (x-mu)**2*(2*x+mu)/(mu**3*x*x))) == 0
    checks['inverse squared remainder decreases above a positive floor'] = sp.simplify(
        sp.diff((2*x+mu)/(mu**3*x*x), x)+2*(x+mu)/(mu**3*x**3)) == 0
    theta = sp.symbols('theta', nonnegative=True)
    kinetic = theta**2/(96*(1-theta/2)**2)
    checks['edge kinetic bound increases throughout the cube window'] = sp.factor(
        sp.diff(kinetic, theta)-theta/(48*(1-theta/2)**3)) == 0
    # One internal edge times an incident face: degree 0 / 2 on their
    # common link, three other fundamental links, hence energies 9 / 17.
    weights = {9: Q(1, 4), 17: Q(3, 4)}
    time = sp.symbols('t', nonnegative=True)
    first = sum(sp.Rational(weight)*(sp.exp(-E*time)/12+
                (sp.exp(-3*time)-sp.exp(-E*time))/(E-3))
                for E, weight in weights.items())-sp.exp(-3*time)/12
    expected = sp.exp(-3*time)/84-sp.exp(-9*time)/48+sp.exp(-17*time)/112
    checks['complete first heat response has all three electric rates'] = sp.simplify(first-expected) == 0
    checks['heat response begins only after its first time derivative'] = (
        expected.subs(time, 0) == 0 and sp.diff(expected, time).subs(time, 0) == 0
        and sp.diff(expected, time, 2).subs(time, 0) == 1)
    z = sp.symbols('z')
    positive = 3*z**5+6*z**4+9*z**3+12*z**2+8*z+4
    checks['first incident-face heat response has a positive polynomial certificate'] = sp.expand(
        4-7*z**3+3*z**7-(1-z)**2*positive) == 0
    heat_coefficient = Q(1, 4)*(Q(1, 84*12)-Q(1, 48*18)+Q(1, 112*26))
    checks['heat kernel and stationary inverse give the same low-support coefficient'] = (
        heat_coefficient == first_channel_coefficient() == Q(1, 22464))
    checks['nonincident cube faces cancel against the reference reset'] = first_channel_coefficient(False) == 0
    checks['positive response weights resolve the correct mean rate fifteen'] = sum(
        E*weight for E, weight in weights.items()) == 15
    edges, faces, _, _ = y12.cube()
    checks['each internal edge sends its first return to exactly two cube faces'] = all(
        sum(e in face for face in faces) == 2 for e in range(len(edges)))
    checks['the cube source has Haar norm squared three halves and energy twelve'] = (
        len(faces)*Q(1, 4) == Q(3, 2) and all(len(face)*3 == 12 for face in faces))
    return checks


def verify_predecessor(relative):
    record = json.loads((ROOT/relative).read_text())
    for path, expected in record['source_pins'].items():
        actual = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f'changed frozen source: {path}')


def source_pins():
    paths = ['physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md',
             'physics/yc12/yc12_correlated_cube_gap.py', 'physics/yc12/YC12_RESULT.json',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc14/yc14_observer_reset_and_interface.py', 'physics/yc14/YC14_RESULT.json',
             'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md',
             'physics/yc15/yc15_boundary_return_channels.py', 'physics/yc15/test_yc15.py']
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    for path in ('physics/yc12/YC12_RESULT.json', 'physics/yc14/YC14_RESULT.json'):
        verify_predecessor(path)
    checks = dict(symbolic_checks())
    worst = source_bounds(1, 1)
    checks['actual cube floor is retained without rounding down to one quarter'] = GAP == Q(1279511, 4672512) > Q(1, 4)
    checks['full source residual variance is bounded by one twelfth'] = worst['variance_upper'] == Q(1, 12)
    checks['susceptibility width is at most one over 11232'] = (
        worst['susceptibility_upper'] == Q(235, 11232))
    checks['inverse vector norm squared uses the complete residual'] = (
        worst['inverse_vector_norm_squared_upper'] == Q(773, 438048))
    checks['rational inverse-vector cap is outward'] = (
        worst['inverse_vector_norm_squared_upper'] < INVERSE_NORM_CAP**2)
    checks['nonconstant cube return has a full-window norm bound'] = (
        worst['cube_return_difference_norm_squared_upper'] == Q(1, 292032))
    checks['zero internal couplings return exactly the scalar forty-eighth'] = (
        source_bounds(0, 0)['susceptibility_upper'] == Q(1, 48)
        and source_bounds(0, 0)['cube_return_difference_norm_squared_upper'] == 0)
    row = join_bounds(ETA_MAX)
    checks['resolved seed costs eighty-one twentieths times eta'] = row['seed'] == Q(81, 20)*ETA_MAX
    checks['new endpoint maps the complete cluster ball strictly inward'] = row['mapping_reserve'] == Q(11417, 409443520) > 0
    checks['full nonlinear contraction remains strict'] = row['contraction'] == Q(28035072, 31987775) < 1
    checks['all excitation channels retain a relative reserve'] = row['relative_return'] == Q(5537792, 31987775) < 1
    checks['actual full-lattice gap is larger than nine fortieths'] = row['full_gap_lower'] == Q(979629, 4326400) > Q(9, 40)
    checks['interface window grows by thirty-two twenty-sevenths'] = ETA_MAX/y14.ETA_MAX == Q(32, 27)
    checks['old endpoint floor also improves'] = join_bounds(y14.ETA_MAX)['full_gap_lower'] > Q(21, 100)
    # At the new endpoint both the old coarse source and the old rounded
    # factor floor fail separately, so both identified improvements matter.
    old_seed = Q(96, 13)*ETA_MAX
    checks['old first-source bound cannot certify the new endpoint'] = old_seed+row['contraction']*RADIUS > RADIUS
    old_q = Q(5184, 5)*ETA_MAX/Q(1, 4)
    checks['discarding the certified charged-cube reserve loses this endpoint'] = row['seed']+old_q*RADIUS > RADIUS
    geometry = [y14.geometry_certificate(shape) for shape in ((4, 4, 4), (4, 6, 8))]
    for item in geometry:
        label = 'x'.join(map(str, item['shape']))
        checks[f'{label}: only a self-paired face can remove its bridge signature'] = item['bridge_signatures_unique']
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([key for key, value in checks.items() if not value])
    return y14.encode(dict(stage='YC15', base_commit='a7ae437', checks=checks,
                          full_window_source=worst,
                          first_single_cube_return_coefficient=Q(1, 22464),
                          heat_response='exp(-3t)/84-exp(-9t)/48+exp(-17t)/112',
                          first_nonincident_face_coefficient=Q(0),
                          endpoint=row, old_endpoint=join_bounds(y14.ETA_MAX),
                          full_gap_formula='gap>=1279511/4672512-(1024/5)*eta; eta<=1/4320',
                          geometry=geometry,
                          claim_boundary=dict(actual_low_support_return_derived=True,
                              complete_source_variance_and_inverse=True,
                              first_internal_coupling_derivative_exact=True,
                              finite_theta_linear_truncation_certified=False,
                              all_cluster_orders_retained_in_gap_proof=True,
                              stronger_correlated_interface_window=True,
                              effective_return_equals_scalar_potential=False,
                              isotropic_window_improved=False, spatial_RG=False,
                              continuum_mass_gap=False, formal_machine_verification=False),
                          source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC15_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC15: {len(result['checks'])} exact checks pass.")
    print('Actual full lattice gap >9/40 through eta=1/4320; first cube return coefficient 1/22464.')
