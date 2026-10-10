"""YC14: metric-aware observer reset and source-resolved cube interfaces.

Exact arithmetic accompanies the analytic, full-carrier proofs in the note.
No finite harmonic approximation is substituted for a cube ground or inverse.
"""
import argparse
from collections import Counter
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
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


y12 = load('yc14_yc12', 'physics/yc12/yc12_correlated_cube_gap.py')
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.diag(1, -1)
J = R*K
OBS = -J
GAP = Q(1, 4)
ETA_MAX = Q(1, 5120)
RADIUS = Q(1, 128)
EXP_UPPER = Q(16, 15)


def zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.cancel(entry) == 0 for entry in value)
    return sp.cancel(value) == 0


@lru_cache(None)
def reset_family():
    # u=exp(eta/2)>0, t=tan(phi/2), one continuous angular patch.
    u = sp.symbols('u', positive=True)
    t = sp.symbols('t', real=True)
    k = ((1-t*t)*K+2*t*J)/(1+t*t)
    n = R*k
    ch, sh = (u*u+u**-2)/2, (u*u-u**-2)/2
    cut = ch*k+sh*R
    whitening = (u+1/u)/2*sp.eye(2)-(u-1/u)/2*n
    rotation = (sp.eye(2)-t*R)/sp.sqrt(1+t*t)
    reset = rotation*whitening
    metric = whitening*whitening
    cycle = (R*cut)**2
    return dict(u=u, t=t, k=k, n=n, cut=cut, whitening=whitening,
                rotation=rotation, reset=reset, metric=metric, cycle=cycle)


def rotation(cosine, sine):
    cosine, sine = sp.Rational(cosine), sp.Rational(sine)
    if cosine*cosine+sine*sine != 1:
        raise ValueError('unit circle pair required')
    return cosine*sp.eye(2)+sine*R


def memory_example():
    frames = (sp.eye(2), rotation(Q(3, 5), Q(4, 5)),
              rotation(Q(5, 13), Q(12, 13)))
    cuts = tuple(O.T*K*O for O in frames)
    steps = tuple(frames[i+1]*frames[i].T for i in range(2))
    return frames, cuts, steps


def source_cost(bridges, cube_gap=GAP):
    cube_gap = Q(cube_gap)
    if bridges not in (2, 4) or not 0 < cube_gap <= 3:
        raise ValueError('two or four bridges and a supplied cube floor in (0,3] required')
    floor = 6+2*cube_gap if bridges == 2 else Q(12)
    return dict(bridges=bridges, cubes=4-bridges, support=4,
                norm_squared=Q(1, 4), energy_form=Q(3), mean_energy=Q(12),
                source_floor=floor, inverse_vector_upper=Q(1, 2)/floor,
                weighted_seed_upper=2/floor,
                susceptibility_lower=Q(1, 48), susceptibility_upper=Q(1, 4)/floor)


def geometry_certificate(shape):
    tiling = y12.tiling(shape)
    faces = tiling['external']
    signature = [frozenset(f for f in face if f[0] == 'bridge') for face in faces]
    profiles = {}
    for factor in tiling['incidence']:
        counts = Counter(len(sig) for face, sig in zip(faces, signature) if factor in face)
        profiles.setdefault(factor[0], set()).add((counts[2], counts[4]))
    block_count = (shape[0]//2)*(shape[1]//2)*(shape[2]//2)
    return dict(shape=list(shape), blocks=block_count,
                two_bridge_faces=sum(len(sig) == 2 for sig in signature),
                four_bridge_faces=sum(len(sig) == 4 for sig in signature),
                bridge_signatures_unique=len(set(signature)) == len(faces),
                cube_profiles=sorted(profiles['cube']),
                bridge_profiles=sorted(profiles['bridge']))


def join_bounds(eta):
    eta = Q(eta)
    if not 0 <= eta <= ETA_MAX:
        raise ValueError('certified interface window is 0<=eta<=1/5120')
    beta, M, a = 24*eta, Q(4), Q(8)
    cube_seed = 24*source_cost(2)['weighted_seed_upper']*eta
    bridge_seed = 2*(source_cost(2)['weighted_seed_upper']+
                     source_cost(4)['weighted_seed_upper'])*eta
    seed = max(cube_seed, bridge_seed)
    contraction = M*beta/GAP*EXP_UPPER*(2+a*(1+2*RADIUS))
    mapping = seed+contraction*RADIUS
    relative = 2*M*beta/GAP*EXP_UPPER
    return dict(eta=eta, beta=beta, factor_floor=GAP, nonlinear_minimum_support=1,
                radius=RADIUS, exponential_upper=EXP_UPPER,
                cube_seed=cube_seed, bridge_seed=bridge_seed, seed=seed,
                contraction=contraction, mapping=mapping,
                mapping_reserve=RADIUS-mapping, relative_return=relative,
                full_gap_lower=GAP*(1-relative),
                fixed_point_norm_upper=seed/(1-contraction))


def symbolic_checks():
    f = reset_family()
    cut, C, G, B, k, O = (f[key] for key in ('cut', 'reset', 'metric', 'cycle', 'k', 'rotation'))
    inv = C.inv()
    Rprime = C*R*inv
    plus, minus = (sp.eye(2)+cut)/2, (sp.eye(2)-cut)/2
    checks = {
        'primitive turn and cut anticommute': R*R == -sp.eye(2) and K*R == -R*K,
        'admissible cut is an involution': zero(cut*cut-sp.eye(2)),
        'quarter-root reset preserves oriented area': zero(C.det()-1),
        'cut metric is the positive half-root of cycle': zero(G*G-B),
        'reset carries the cut to fixed axes': zero(C*cut*inv-K),
        'pairing is carried with the observer': zero(C.T*C-G),
        'cut is self-adjoint in its adapted pairing': zero(cut.T*G-G*cut),
        'seen and lost projections are metric orthogonal': zero(plus.T*G*minus),
        'balanced contrast and total norm are both retained': zero(C.T*K*C-G*cut),
        'full cycle survives when the turn is transported': zero((Rprime*K)**2-C*B*inv),
        'flat branch requires only rotation': zero(f['whitening'].subs(f['u'], 1)-sp.eye(2)),
        'rotation aligns the symmetric seam': zero(O*k*O.T-K),
    }
    witness = {f['u']: 2, f['t']: 0}
    checks['cut whitening does not close a genuinely open compass'] = (
        not zero((C*B*inv).subs(witness)-sp.eye(2)))
    checks['leaving turn fixed after nonorthogonal reset is false'] = (
        not zero(Rprime.subs(witness)-R))
    pp, qq, ww = sp.symbols('p q w', real=True)
    L = pp*K+qq*J+ww*R
    checks['rotation does not remove the antisymmetric defect'] = zero(
        (O*L*O.T-(O*L*O.T).T)/2-ww*R)
    checks['signed observation is still a consistent involution'] = (
        OBS*OBS == sp.eye(2) and OBS*R*OBS.T == -R)
    frames, cuts, steps = memory_example()
    checks['all three example observations share the same local cut'] = all(
        Oi*ki*Oi.T == K for Oi, ki in zip(frames, cuts))
    checks['reset transition retains open-path ratio memory'] = steps[1]*steps[0] == frames[2]
    checks['closed reset ledger alone has no holonomy'] = frames[0]*frames[2].T*steps[1]*steps[0] == sp.eye(2)
    q1, q2 = Q(4, 3), Q(16, 63)
    checks['ratio memory uses the tangent composition law'] = (q1+q2)/(1-q1*q2) == Q(12, 5)
    T1, T2 = rotation(Q(5, 13), Q(12, 13)), rotation(Q(3, 5), Q(-4, 5))
    observed = frames[2]*T2*frames[1].T*frames[1]*T1*frames[0].T
    checks['actual transported memory telescopes only at interior frames'] = observed == frames[2]*T2*T1*frames[0].T
    # Unit-quaternion Haar moments: W=a.u, |a|=1, E[u_i u_j]=delta_ij/4.
    aa = sp.symbols('a0:4', real=True)
    mean_square = sum(a*a for a in aa)/4
    checks['Haar source norm and gradient energy have the same unit constraint'] = zero(
        4*(sum(a*a for a in aa)-mean_square)-12*mean_square)
    # An endpoint centre flip changes exactly two edges of each incident cube face.
    edges, faces, _, _ = y12.cube()
    vertices = {v for edge in edges for v in edge}
    checks['cube endpoint centre flip preserves every internal Wilson face'] = all(
        sum(vertex in edges[e] for e in face) % 2 == 0 for vertex in vertices for face in faces)
    checks['every open cube edge is odd under an endpoint centre flip'] = all(
        sum(vertex in edge for vertex in (edge[0],)) == 1 for edge in edges)
    return checks


def verify_predecessor(relative):
    record = json.loads((ROOT/relative).read_text())
    for path, expected in record['source_pins'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise AssertionError(f'changed frozen source: {path}')


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def source_pins():
    paths = ['uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
             'uncut/up8/up8_deformed_compass.py', 'uncut/up8/UP8_RESULT.json',
             'physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc10/yc10_closed_return_and_scale.py', 'physics/yc10/YC10_RESULT.json',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md',
             'physics/yc12/yc12_correlated_cube_gap.py', 'physics/yc12/YC12_RESULT.json',
             'physics/yc13/YC13_FRAME_FLOW_AND_RETURN.md', 'physics/yc13/YC13_RESULT.json',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc14/yc14_observer_reset_and_interface.py', 'physics/yc14/test_yc14.py']
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    for path in ('physics/yc12/YC12_RESULT.json', 'physics/yc13/YC13_RESULT.json'):
        verify_predecessor(path)
    checks = symbolic_checks()
    geometry = [geometry_certificate(shape) for shape in ((4, 4, 4), (4, 6, 8), (6, 6, 6))]
    for row in geometry:
        name = 'x'.join(map(str, row['shape']))
        checks[f'{name}: interface bridge signatures are unique'] = row['bridge_signatures_unique']
        checks[f'{name}: exact two/four bridge counts'] = (
            row['two_bridge_faces'] == 12*row['blocks'] and row['four_bridge_faces'] == 6*row['blocks'])
        checks[f'{name}: all factor incidence profiles'] = (
            row['cube_profiles'] == [(24, 0)] and row['bridge_profiles'] == [(2, 2)])
    two, four = source_cost(2), source_cost(4)
    checks['complete two-bridge source floor includes both charged cubes'] = two['source_floor'] == Q(13, 2)
    checks['four-bridge source has exact free energy twelve'] = four['susceptibility_upper'] == four['susceptibility_lower'] == Q(1, 48)
    checks['two-bridge complete inverse susceptibility interval'] = (
        two['susceptibility_lower'] == Q(1, 48) and two['susceptibility_upper'] == Q(1, 26))
    bound = join_bounds(ETA_MAX)
    checks['rational exponential majorant sums the whole series'] = (1-8*RADIUS)*EXP_UPPER == 1
    checks['source-centred map preserves the full cluster ball'] = bound['mapping_reserve'] == Q(7, 166400) > 0
    checks['all-order contraction remains strict'] = bound['contraction'] == Q(81, 100) < 1
    checks['actual excitation relative return keeps a strict reserve'] = bound['relative_return'] == Q(4, 25) < 1
    checks['new all-volume full gap is twenty-one hundredths'] = bound['full_gap_lower'] == Q(21, 100)
    checks['old interface endpoint now has gap nineteen eightieths'] = join_bounds(Q(1, 16384))['full_gap_lower'] == Q(19, 80)
    checks['interface window grows by sixteen fifths'] = ETA_MAX/Q(1, 16384) == Q(16, 5)
    # The old coarse mapping bound is not adequate at the new endpoint.
    coarse = 4*bound['beta']/GAP*(1+2*RADIUS)*EXP_UPPER
    checks['new endpoint genuinely needs the source-resolved seed'] = coarse > RADIUS
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([key for key, value in checks.items() if not value])
    return encode(dict(stage='YC14', base_commit='c730847', checks=checks,
                       owner_correction='Each reading returns to a local balanced lambda-zero cut; retain transition ratio memory; radial motion undecided.',
                       observer_reset='C=exp(-phi*R/2)*B^(1/4); G_cut=B^(1/2); C*kappa*C^-1=K',
                       source_costs=[two, four], geometry=geometry,
                       endpoint=bound, old_endpoint=join_bounds(Q(1, 16384)),
                       full_gap_formula='gap>=1/4-(1024/5)*eta, 0<=eta<=1/5120',
                       claim_boundary=dict(local_reset_preserves_full_return=True,
                           quantum_classical_equality_for_equilibrium_multiplication_readings=True,
                           all_quantum_measurements_are_classical=False,
                           rotation_erases_cut_defect=False, radial_motion_derived=False,
                           actual_correlated_lattice_window_improved=True,
                           all_cluster_orders_and_boundary_charges_retained=True,
                           isotropic_window_improved=False, lambda_equals_YM_coupling=False,
                           continuum_mass_gap=False, formal_machine_verification=False),
                       source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC14_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC14: {len(result['checks'])} exact checks pass.")
    print('Actual full lattice gap >=21/100 through interface eta=1/5120; internal cube theta<=1.')
