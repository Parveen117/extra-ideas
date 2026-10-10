"""YC24: actual eight-sector bridge return, full metric tails and gap transfer.

The companion note proves the all-harmonic operator statements. Finite
matrices below are independent algebra controls, never a YM truncation.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    'yc24_y22', ROOT/'physics/yc22/yc22_gauge_protected_channels.py')
y22 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = y22
spec.loader.exec_module(y22)
y20 = y22.y21.y20

FACE_MASKS = (3, 6, 12, 9)
E_MASKS = FACE_MASKS
B_MASKS = (0, 15, 5, 10)
B_FLOORS = (8, 12, 6, 6)
SOURCES = (
    'uncut/ls2/LS2_RESPONSE_TOWER_LIFT.md',
    'uncut/ls2/LS2_RESULT.json',
    'physics/yc13/YC13_FRAME_FLOW_AND_RETURN.md',
    'physics/yc13/YC13_RESULT.json',
    'physics/yc17/YC17_FOUR_FACE_BRIDGE_RETURN.md',
    'physics/yc20/YC20_METRIC_RETAINED_RETURN.md',
    'physics/yc21/YC21_TWO_CUBE_BLOCK_AND_REMOVED_CHANNELS.md',
    'physics/yc21/YC21_RESULT.json',
    'physics/yc22/YC22_GAUGE_PROTECTED_CHANNELS.md',
    'physics/yc22/yc22_gauge_protected_channels.py',
    'physics/yc22/YC22_RESULT.json',
    'physics/yc23/YC23_PRIMITIVE_CUT_AND_AGITATION.md',
    'physics/yc24/YC24_PARITY_RESOLVED_RESPONSE_TOWER.md',
    'physics/yc24/yc24_parity_response_tower.py',
    'physics/yc24/test_yc24.py',
)


def parity(mask):
    return -1 if (mask & 5).bit_count() % 2 else 1


def comparison_constants(z):
    z = Q(z)
    if z >= 6:
        raise ValueError('z must be below the reference E-sector floor6')
    a = 6-z
    c = 4*sum(1/(Q(d)-z) for d in B_FLOORS)
    b = 4*sum(1/(Q(d)-z)**2 for d in B_FLOORS)
    return dict(z=z, a=a, c=c, b=b, tau=c/a)


def response_bounds(t=Q(1), z=Q(0), xi=(1, 1, 1, 1), order=2):
    t, z = Q(t), Q(z)
    xi = tuple(Q(x) for x in xi)
    if not 0 <= t <= 1 or len(xi) != 4 or any(abs(x) > 1 for x in xi):
        raise ValueError('0<=t<=1 and four fixed |xi_i|<=1 required')
    if type(order) is not int or order < 0:
        raise ValueError('nonnegative integer truncation order required')
    row = comparison_constants(z)
    q, v0 = t*t, sum(x*x for x in xi)/4
    alpha = row['a']-q*row['c']
    if alpha <= 0:
        raise ValueError('the complete hidden comparison is not positive here')
    rho = q*row['tau']
    metric_excess = q*v0*(1+q*row['b'])/alpha**2
    relative_tail = rho**(order+1)
    return dict(**row, t=t, q=q, xi=xi, v0=v0, alpha=alpha, rho=rho,
                return_cap=q*v0/alpha, metric_excess=metric_excess,
                metric_cap=1+metric_excess, order=order,
                energy_tail=q*v0/row['a']*relative_tail/(1-rho),
                metric_error=metric_excess*(2*relative_tail-relative_tail**2))


def block_bounds(theta_cap=Q(1), t=Q(1)):
    theta_cap = Q(theta_cap)
    if not 0 <= theta_cap <= 2:
        raise ValueError('inherited cube window is 0<=theta_c<=2')
    delta, z = ((Q(11), Q(12, 5)) if theta_cap <= 1 else (Q(27, 5), Q(7, 3)))
    ret = response_bounds(t, z)
    return dict(theta_cap=theta_cap, interface_cap=Q(t), cube_physical_floor=delta,
                physical_gap=z, hidden_floor=Q(5, 2), full_gap=Q(3, 4),
                alpha=ret['alpha'], return_cap=ret['return_cap'],
                spectral_reserve=delta-z-ret['return_cap'])


def join_bounds(kind):
    if kind not in ('tube', 'force'):
        raise ValueError('choose the inherited tube or force profile')
    row = y22.join_bounds(kind)
    row['predecessor_physical_gap'] = row['physical_gap']
    row['physical_reference_floor'] = Q(7, 3) if kind == 'tube' else Q(12, 5)
    row['physical_gap'] = row['physical_reference_floor']*(1-row['relative_return'])
    return row


def scalar_hidden_comparison():
    """Eight-sector norm comparison, not an eight-harmonic YM Hamiltonian."""
    c = sp.ones(4)
    return (6*sp.eye(4)).row_join(-c).col_join((-c.T).row_join(sp.diag(*B_FLOORS)))


def finite_control(t=sp.Rational(1, 2), z=sp.S.Zero, order=0):
    """A noncommuting finite graph control with two retained columns."""
    t, z = sp.sympify(t), sp.sympify(z)
    A0, B0 = sp.diag(6, 9), sp.Matrix([[8]])
    A, B = A0-z*sp.eye(2), B0-z*sp.eye(1)
    C, R = sp.Matrix([1, sp.Rational(1, 2)]), sp.eye(2)/2
    D0 = A0.row_join(-t*C).col_join((-t*C.H).row_join(B0))
    D = D0-z*sp.eye(3)
    source = R.col_join(sp.zeros(1, 2))
    Z = (D.inv()*(t*source)).applyfunc(sp.cancel)
    sigma = (t*source.H*Z).applyfunc(sp.cancel)
    U, step = A.inv()*R, A.inv()*C*B.inv()*C.H
    power, total = U, U
    for n in range(1, order+1):
        power = step*power
        total += t**(2*n)*power
    ZE = t*total
    ZN = ZE.col_join(t*B.inv()*C.H*ZE)
    sigmaN = t*R.H*ZE
    k = sp.diag(0, 11)
    H = k.row_join(-t*source.H).col_join((-t*source).row_join(D0))
    return dict(A=A, B=B, C=C, R=R, D=D, H=H, Z=Z, ZN=ZN,
                sigma=sigma, sigmaN=sigmaN, metric=sp.eye(2)+Z.H*Z,
                metricN=sp.eye(2)+ZN.H*ZN, step=step,
                F=k-z*sp.eye(2)-sigma)


def encode(obj):
    if isinstance(obj, (Q, sp.Rational)):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(x) for x in obj]
    return obj


def certificate():
    # Replay hashes, not predecessor writers or modified theorem packets.
    for path in ('physics/yc13/YC13_RESULT.json', 'physics/yc21/YC21_RESULT.json',
                 'physics/yc22/YC22_RESULT.json'):
        y22.y15.verify_predecessor(path)
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not checks[name]:
            raise AssertionError(name)

    even = {n for n in range(16) if n.bit_count() % 2 == 0}
    check('eight complete even signatures split into four plus four',
          set(E_MASKS).isdisjoint(B_MASKS) and set(E_MASKS) | set(B_MASKS) == even)
    check('alternating centre cut has the stated signs',
          all(parity(m) == -1 for m in E_MASKS) and all(parity(m) == 1 for m in B_MASKS))
    check('every actual face crosses the alternating cut exactly once',
          all(parity(m ^ f) == -parity(m) for m in even for f in FACE_MASKS))
    check('hidden signature incidence is the complete bipartite four by four graph',
          all(sum((e ^ f) == b for f in FACE_MASKS) == 1 for e in E_MASKS for b in B_MASKS))

    # Verify the cycle on the actual28-link geometry, not only supplied masks.
    graph = y22.open_block()
    bridge_order = [((1, yy, zz), 0) for yy, zz in ((0,0),(1,0),(1,1),(0,1))]
    actual_masks = []
    for face in graph['faces']:
        mask = sum(1 << j for j, bridge in enumerate(bridge_order) if bridge in face)
        if mask:
            actual_masks.append(mask)
    check('actual28-link block has precisely the four declared bridge face masks',
          len(graph['edges']) == 28 and len(actual_masks) == 4 and set(actual_masks) == set(FACE_MASKS))

    floors = {}
    for ns in product(range(4), repeat=4):
        if not any(ns) or sum(ns) % 2:
            continue
        mask = sum((n % 2) << j for j, n in enumerate(ns))
        energy = sum(n*(n+2) for n in ns)
        floors[mask] = min(floors.get(mask, energy), energy)
    check('bridge harmonic floors retain the empty-signature nonconstant sector',
          all(floors[e] == 6 for e in E_MASKS)
          and tuple(floors[b] for b in B_MASKS) == B_FLOORS)
    check('odd-length hidden returns are forbidden at every word length',
          all(parity(m) == -parity(m ^ f) for m in even for f in FACE_MASKS)
          and parity(0) == 1)

    c0, ch = comparison_constants(0), comparison_constants(Q(1,2))
    check('zero-energy positive tower norm is at most thirteen over36',
          c0['c'] == Q(13,6) and c0['b'] == Q(5,16) and c0['tau'] == Q(13,36))
    hidden = response_bounds(1, Q(5,2))
    check('whole hidden block retains a positive reserve at five halves', hidden['alpha'] == Q(193,2926))
    comparison = scalar_hidden_comparison()-sp.Rational(5,2)*sp.eye(8)
    check('independent eight-sector comparison is strictly positive at five halves',
          all(comparison[:i,:i].det() > 0 for i in range(1,9)))
    zero_row, half_row = response_bounds(), response_bounds(1,Q(1,2))
    check('complete zero-energy return and metric caps',
          zero_row['return_cap'] == Q(6,23) and zero_row['metric_excess'] == Q(189,2116))
    check('full-strength returned metric is at most eight sevenths through z one half',
          half_row['metric_excess'] == Q(78682276,576816289)
          and half_row['metric_cap'] < Q(8,7))

    z, q = sp.symbols('z q', real=True)
    constants = [8,12,6,6]
    cz = 4*sum(1/(d-z) for d in constants)
    bz = 4*sum(1/(d-z)**2 for d in constants)
    az = 6-z-q*cz
    check('whole-window monotonicity follows from exact denominator derivatives',
          sp.cancel(sp.diff(cz,z)-bz) == 0
          and sp.cancel(sp.diff(az,z)+1+q*bz) == 0
          and sp.diff(az,q) == -cz)

    weak, strong = block_bounds(2), block_bounds(1)
    check('physical block floor seven thirds has an exact positive reserve',
          weak['spectral_reserve'] == Q(29251,89115) and weak['physical_gap'] < Q(5,2))
    check('physical block floor twelve fifths has an exact positive reserve',
          strong['spectral_reserve'] == Q(7073,1555) and strong['physical_gap'] < Q(5,2))
    check('full charged floor is kept separate from the improved physical floor',
          strong['full_gap'] == weak['full_gap'] == Q(3,4))
    tube, force = join_bounds('tube'), join_bounds('force')
    check('strong-cube lattice join inherits physical gap461 over225',
          tube['physical_gap'] == Q(461,225) and tube['full_gap'] == Q(461,700))
    check('resolved-force lattice join inherits physical gap6476 over3125',
          force['physical_gap'] == Q(6476,3125) and force['full_gap'] == Q(1619,2500))
    check('full nonlinear return and creator ball reserves are unchanged and positive',
          tube['relative_return'] == Q(64,525) and force['relative_return'] == Q(256,1875)
          and tube['ball_reserve'] == Q(3,22400) and force['ball_reserve'] == Q(53,400000))
    check('boundary incidence over full gap has not become contractive',
          Q(48)/strong['full_gap'] == 64 > 24)

    control = finite_control()
    A, B, C, R = [control[k] for k in ('A','B','C','R')]
    t = sp.Rational(1,2)
    reduced = t*t*R.H*(A-t*t*C*B.inv()*C.H).inv()*R
    check('full inverse and nested bipartite elimination agree without commutation',
          y20.equal(control['sigma'], reduced)
          and A*C*B.inv()*C.H != C*B.inv()*C.H*A)
    minus = finite_control(-t)
    check('complete hidden return is even in the interface scale',
          y20.equal(control['sigma'], minus['sigma']))
    check('neutral cross-source mixing survives the complete inverse',
          control['sigma'][0,1] == sp.Rational(1,54960) and R[:,0].dot(R[:,1]) == 0)

    sigma_q = q*R.H*(A-q*C*B.inv()*C.H).inv()*R
    step = A.inv()*C*B.inv()*C.H
    for n in (1,2,3):
        derivative = (sp.factorial(n)*R.H*step**(n-1)
                      *(sp.eye(2)-q*step).inv()**(n+1)*A.inv()*R)
        check(f'actual algebra control for response derivative order{n}',
              y20.equal(sigma_q.diff(q,n), derivative))

    tail_rows = []
    for n in (0,1,2,3):
        finite = finite_control(order=n)
        bound = response_bounds(Q(1,2), 0, xi=(Q(1,2),)*4, order=n)
        tail = finite['sigma']-finite['sigmaN']
        error = finite['metric']-finite['metricN']
        ecap = sp.Rational(bound['energy_tail'])
        mcap = sp.Rational(bound['metric_error'])
        check(f'complete energy and two-sided norm tails enclose finite control order{n}',
              y20.positive_semidefinite(tail)
              and y20.positive_semidefinite(ecap*sp.eye(2)-tail)
              and y20.positive_semidefinite(mcap*sp.eye(2)-error)
              and y20.positive_semidefinite(mcap*sp.eye(2)+error))
        tail_rows.append(response_bounds(order=n))
    check('positive energy tail does not imply a one-sided metric tail',
          sp.factor((control['metric']-control['metricN']).det()) == -sp.Rational(1,978674918400))
    variable = finite_control(z=z)
    check('complete spectral derivative equals the original graph norm',
          y20.equal(-variable['F'].diff(z),variable['metric']))
    counted = finite_control(z=sp.Rational(12,5))
    check('full finite spectral inertia agrees with the retained pencil',
          y20.inertia(counted['H']-sp.Rational(12,5)*sp.eye(5))[0]
          == y20.inertia(counted['F'])[0] == 1)
    check('quoted full-strength sixth-order energy and metric tails are exact',
          tail_rows[2]['energy_tail'] == Q(2197,178848)
          and tail_rows[2]['metric_error'] == Q(1401257585,170595237888))
    check('quoted full-strength eighth-order energy and metric tails are exact',
          tail_rows[3]['energy_tail'] == Q(28561,6438528)
          and tail_rows[3]['metric_error'] == Q(665891061017,221091428302848))

    return encode(dict(
        stage='YC24', base_commit='7648401', exact_check_count=len(checks), checks=checks,
        masks=dict(E=E_MASKS,B=B_MASKS,faces=FACE_MASKS), B_floors=B_FLOORS,
        zero_response=zero_row, half_response=half_row, hidden_floor=hidden,
        cube_cap_two=weak, cube_cap_one=strong, lattice_tube=tube,lattice_force=force,
        full_strength_tail_rows=tail_rows,
        source_pins={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCES},
        classical_reference='https://arxiv.org/abs/2105.02058',
        claims=dict(actual28_link_parity_adapter=True,
                    complete_positive_coupling_squared_tower=True,
                    full_energy_and_graph_metric_remainder_bounds=True,
                    stronger_physical_block_and_lattice_floors=True,
                    all_hidden_harmonics_retained=True,
                    neutral_cross_source_mixing_retained=True,
                    wider_external_window=False,
                    positive_metric_tail_as_operator=False,
                    finite_source_span_closed=False,
                    local_tower_ratio_is_RG_contraction=False,
                    iterated_RG_contraction=False, continuum_mass_gap=False,
                    formal_machine_verification=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = HERE/'YC24_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result,indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise AssertionError('certificate or source pins changed')
    print(f"YC24: {result['exact_check_count']} exact checks pass.")
    print('Physical block gaps >=7/3 and12/5; physical lattice floors461/225 and6476/3125.')
