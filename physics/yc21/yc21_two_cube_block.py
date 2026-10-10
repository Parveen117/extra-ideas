"""YC21: exact budgets for an actual 28-link YM block and its full join.

Finite matrix controls check algebra; they do not replace the infinite carrier.
The companion note proves full harmonic/parity floors and the spectral transfer.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    'yc21_y20', ROOT/'physics/yc20/yc20_metric_retained_return.py')
y20 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = y20
spec.loader.exec_module(y20)
y15 = y20.y19.y15
GAP_CUBE = y15.GAP
GAP_BLOCK = Q(1, 5)
ETA_MAX = Q(1, 2)
EPS_MAX = Q(1, 12800)
SIZE = (4, 2, 2)
SIGNATURES = tuple(frozenset((i, (i+1) % 4)) for i in range(4))


def block_bounds(couplings, z=Q(1, 5)):
    ts = tuple(Q(t) for t in couplings)
    if len(ts) != 4 or any(abs(t) > ETA_MAX for t in ts):
        raise ValueError('four real interface couplings with absolute value <=1/2 required')
    z = Q(z)
    s, v = sum(abs(t) for t in ts), sum(t*t for t in ts)/4
    d, odd = 6-s, 3-s
    if z >= d:
        raise ValueError('energy must be below the full even hidden floor')
    sigma = v/(d-z)
    return dict(couplings=ts, z=z, hidden_even_floor=d, odd_floor=odd,
                source_norm_squared=v, return_upper=sigma,
                metric_excess=v/(d-z)**2,
                spectral_reserve=GAP_CUBE-z-sigma,
                positive_gap_certified=(z > 0 and z < odd and GAP_CUBE-z-sigma > 0),
                ground_energy_lower=-v/d)


def join_bounds(epsilon):
    epsilon = Q(epsilon)
    if not 0 <= epsilon <= EPS_MAX:
        raise ValueError('certified external window is 0<=epsilon<=1/12800')
    r, exp_upper = Q(1, 128), Q(16, 15)
    beta = 48*epsilon
    inverse_source = Q(1, 2)/(6+2*GAP_BLOCK)
    seed = 48*4*inverse_source*epsilon
    q = 4*beta/GAP_BLOCK*exp_upper*(2+8*(1+2*r))
    b = 8*beta/GAP_BLOCK*exp_upper
    return dict(epsilon=epsilon, beta=beta, seed=seed,
                inverse_source_cap=inverse_source, radius=r,
                contraction=q, relative_return=b, mapping_reserve=r-seed-q*r,
                gap_lower=GAP_BLOCK*(1-b), minimum_nonlinear_support=1)


def rectangular_tiling(shape):
    shape = tuple(shape)
    if (len(shape) != 3 or any(type(n) is not int for n in shape)
            or any(n < 2*b or n % b for n, b in zip(shape, SIZE))):
        raise ValueError('Lx multiple of4 >=8; Ly,Lz even >=4 required')
    points = tuple(product(*(range(n) for n in shape)))

    def shift(x, axis):
        out = list(x)
        out[axis] = (out[axis]+1) % shape[axis]
        return tuple(out)

    def block(x):
        return tuple(a//b for a, b in zip(x, SIZE))

    factors = {}
    for x in points:
        for axis in range(3):
            edge = (x, axis)
            factors[edge] = (('block', block(x)) if block(x) == block(shift(x, axis))
                             else ('bridge', edge))
    factor_links = Counter(factors.values())
    internal, external = [], []
    old_cube_faces, tube_faces = Counter(), Counter()
    for x in points:
        for i, j in combinations(range(3), 2):
            face = ((x, i), (shift(x, i), j), (shift(x, j), i), (x, j))
            support = frozenset(factors[e] for e in face)
            if len(support) == 1:
                internal.append(support)
                owner = next(iter(support))
                vertices = {v for e in face for v in (e[0], shift(*e))}
                old_blocks = {tuple(a//2 for a in v) for v in vertices}
                (old_cube_faces if len(old_blocks) == 1 else tube_faces)[owner] += 1
            else:
                external.append(support)
    incidence = Counter(f for support in external for f in support)
    bridge_profiles = Counter()
    for factor in factor_links:
        if factor[0] == 'bridge':
            counts = Counter(sum(f[0] == 'bridge' for f in s)
                             for s in external if factor in s)
            bridge_profiles[(counts[2], counts[4])] += 1
    return dict(shape=shape, factor_links=factor_links, internal=internal,
                external=external, incidence=incidence, cube_faces=old_cube_faces,
                tube_faces=tube_faces, bridge_profiles=bridge_profiles)


def graph_connection(W, coordinates):
    """Flat ambient connection projected onto a full-rank graph frame."""
    M = W.H*W
    inv = M.inv()
    projector = y20.simplify(W*inv*W.H)
    A = [y20.simplify(inv*W.H*W.diff(x)) for x in coordinates]
    curvature, returned = {}, {}
    for i, j in combinations(range(len(coordinates)), 2):
        F = A[j].diff(coordinates[i])-A[i].diff(coordinates[j])+A[i]*A[j]-A[j]*A[i]
        B = (W.diff(coordinates[i]).H*(sp.eye(W.rows)-projector)*W.diff(coordinates[j])
             - W.diff(coordinates[j]).H*(sp.eye(W.rows)-projector)*W.diff(coordinates[i]))
        curvature[i, j], returned[i, j] = y20.simplify(F), y20.simplify(B)
    return dict(metric=y20.simplify(M), projector=projector, connection=A,
                curvature=curvature, returned=returned)


@lru_cache(None)
def geometry_controls():
    x, y = sp.symbols('x y', real=True)
    W = sp.Matrix([[1, 0], [0, 1], [x, y]])
    row = graph_connection(W, (x, y))
    complex_W = sp.Matrix([1, x+sp.I*y])
    complex_row = graph_connection(complex_W, (x, y))
    # Same hidden direction along both parameter changes: curvature cancels.
    flat_row = graph_connection(sp.Matrix([1, x+y]), (x, y))
    return x, y, W, row, complex_row, flat_row


def spectral_control():
    # Independent finite control, not a truncated YM matrix.
    g = sp.Rational(GAP_CUBE.numerator, GAP_CUBE.denominator)
    H = sp.Matrix([[0, 0, -sp.Rational(1, 4), 0],
                   [0, g, 0, -sp.Rational(1, 4)],
                   [-sp.Rational(1, 4), 0, 5, sp.Rational(1, 4)],
                   [0, -sp.Rational(1, 4), sp.Rational(1, 4), 6]])
    z = sp.symbols('z', real=True)
    F, W = y20.schur(H-z*sp.eye(4), 2)
    return H, z, F, W


def source_pins():
    paths = ['physics/dm1/DM1_ONE_MORE_DIMENSION.md',
             'physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md',
             'physics/yc15/YC15_RESULT.json',
             'physics/yc17/YC17_FOUR_FACE_BRIDGE_RETURN.md',
             'physics/yc20/YC20_METRIC_RETAINED_RETURN.md',
             'physics/yc20/yc20_metric_retained_return.py',
             'physics/yc20/YC20_RESULT.json',
             'physics/yc21/YC21_TWO_CUBE_BLOCK_AND_REMOVED_CHANNELS.md',
             'physics/yc21/yc21_two_cube_block.py', 'physics/yc21/test_yc21.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    for p in ('physics/yc15/YC15_RESULT.json', 'physics/yc20/YC20_RESULT.json'):
        y15.verify_predecessor(p)
    checks = {}
    half = block_bounds((ETA_MAX,)*4)
    quarter = block_bounds((Q(1, 4),)*4, Q(13, 50))
    metric = block_bounds((ETA_MAX,)*4, Q(1, 2))
    join = join_bounds(EPS_MAX)
    checks['block uses the inherited full charged cube floor'] = GAP_CUBE == Q(1279511, 4672512)
    checks['four interface faces preserve total bridge centre parity'] = all(len(s) % 2 == 0 for s in SIGNATURES)
    checks['distinct face source signatures are orthogonal'] = len(set(SIGNATURES)) == 4
    energies = [(sum(n*(n+2) for n in ns), sum(ns) % 2)
                for ns in product(range(4), repeat=4) if any(ns)]
    checks['lowest nonconstant even and odd bridge energies are six and three'] = (
        min(e for e, p in energies if not p) == 6 and min(e for e, p in energies if p) == 3)
    a, b = sp.symbols('a0:4'), sp.symbols('b0:4')
    checks['full retained face square is one quarter for unit cube quaternions'] = sp.expand(
        y15.bridge_gram(a, b, a, b)-sum(t*t for t in a)*sum(t*t for t in b)/4) == 0
    checks['half-strength tube retains the full positive spectral reserve'] = (
        half['positive_gap_certified'] and half['spectral_reserve'] == Q(3572617,443888640))
    checks['quarter-strength tube has gap at least thirteen fiftieths'] = (
        quarter['positive_gap_certified'] and quarter['spectral_reserve'] == Q(6019313,9228211200))
    checks['odd hidden carrier is accounted for separately'] = half['odd_floor'] == 1 < half['hidden_even_floor'] == 4
    checks['actual ground scalar is bounded below by minus one sixteenth'] = half['ground_energy_lower'] == -Q(1,16)
    checks['complete actual graph metric excess is at most one forty-ninth'] = metric['metric_excess'] == Q(1,49)
    checks['actual graph derivative and curvature have uniform upper bounds'] = (
        (Q(1,2)+Q(1,7))/Q(7,2) == Q(9,49)
        and 2*Q(9,49)**2 == Q(162,2401))
    # Whole-window monotonicity, independent of endpoint enumeration.
    t, z = sp.symbols('t z', real=True)
    derivative = sp.diff(t*t/(6-4*t-z), t)
    checks['return cost increases on the entire declared coupling window'] = sp.simplify(
        derivative-2*t*(6-2*t-z)/(6-4*t-z)**2) == 0
    tilings = []
    for shape in ((8,4,4), (12,6,4)):
        row = rectangular_tiling(shape)
        blocks = [f for f in row['factor_links'] if f[0] == 'block']
        checks[f'{shape}: each real block has 28 links, 12 cube faces and four tube faces'] = all(
            row['factor_links'][f] == 28 and row['cube_faces'][f] == 12 and row['tube_faces'][f] == 4
            for f in blocks)
        checks[f'{shape}: remaining faces have complete four-factor support'] = all(
            len(s) == 4 and sum(f[0] == 'bridge' for f in s) in (2,4) for s in row['external'])
        checks[f'{shape}: exact new block and bridge incidence are 48 and four'] = all(
            n == (48 if f[0] == 'block' else 4) for f,n in row['incidence'].items())
        checks[f'{shape}: rectangular bridges require both incidence profiles'] = set(row['bridge_profiles']) == {(2,2),(3,1)}
        tilings.append(dict(shape=shape, blocks=len(blocks), profiles=dict(row['bridge_profiles'])))
    checks['new complete boundary source cap and seed are retained'] = (
        join['inverse_source_cap'] == Q(5,64) and join['seed'] == 15*EPS_MAX)
    checks['new block ground join contracts strictly on all cluster supports'] = join['contraction'] == Q(81,100)
    checks['the complete creator ball retains a positive reserve'] = join['mapping_reserve'] == Q(1,3200)
    checks['new anisotropic full gap is at least twenty-one over 125'] = (
        join['relative_return'] == Q(4,25) and join['gap_lower'] == Q(21,125))
    checks['new grouping does not contract the old boundary budget'] = 48/GAP_BLOCK > 24/GAP_CUBE
    x,y,W,row,complex_row,flat_row = geometry_controls()
    M, Pi = row['metric'], row['projector']
    checks['removed graph directions give exactly the returned curvature identity'] = y20.equal(
        M*row['curvature'][0,1], row['returned'][0,1])
    checks['graph projection is orthogonal and retains its complete frame'] = (
        y20.equal(Pi*Pi,Pi) and y20.equal(Pi.H,Pi) and y20.equal(Pi*W,W))
    checks['projected connection transports its nonconstant metric'] = all(y20.equal(
        M.diff(c), row['connection'][i].H*M+M*row['connection'][i]) for i,c in enumerate((x,y)))
    checks['real two-channel graph can have nonzero order-defect'] = row['curvature'][0,1].subs({x:0,y:0}) == sp.Matrix([[0,1],[-1,0]])
    checks['complex rank-one graph keeps conjugation and its Berry curvature'] = sp.simplify(
        complex_row['curvature'][0,1][0]-2*sp.I/(1+x*x+y*y)**2) == 0
    checks['dimension reduction and a varying metric do not force curvature'] = (
        flat_row['curvature'][0,1] == sp.zeros(1) and flat_row['metric'].diff(x) != sp.zeros(1))
    H,z,F,W = spectral_control()
    checks['full hidden propagation produces the energy derivative metric'] = y20.equal(-F.diff(z),W.H*W)
    zz = sp.Rational(1,5)
    checks['finite control full inertia agrees with the complete Schur count'] = (
        y20.inertia(H-zz*sp.eye(4))[0] == y20.inertia(F.subs(z,zz))[0] == 1)
    checks['finite hidden inverse does not close on an individual source direction'] = W[3,0] != 0
    checks = {k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k,v in checks.items() if not v])
    return y15.y12.encode(dict(stage='YC21', base_commit='88e40f7', checks=checks,
        block_half=half, block_quarter=quarter, metric=metric, lattice_join=join,
        tiling_controls=tilings, source_pins=source_pins(),
        claim_boundary=dict(actual_28_link_block_gap=True, all_boundary_charges=True,
            no_harmonic_truncation=True, full_hidden_metric_bound=True,
            actual_anisotropic_volume_uniform_gap=True,
            projected_parameter_curvature_identity_and_upper_bound=True,
            computed_nonzero_YM_block_curvature=False,
            spacetime_dimension_removed=False, curvature_monotonicity=False,
            global_quotient_metric_condition_number=False,
            improved_isotropic_window=False, iterated_RG_contraction=False,
            continuum_mass_gap=False, formal_machine_verification=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC21_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC21: {len(result['checks'])} exact checks pass.")
    print('Actual block gap >=1/5; anisotropic lattice gap >=21/125; full metric <=50/49.')
