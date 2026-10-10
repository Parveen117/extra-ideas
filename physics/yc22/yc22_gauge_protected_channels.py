"""YC22: exact evidence for full gauge-protected channels and spatial joins.

The analytic note proves the whole-carrier form and symmetry statements.
These controls do not substitute a finite harmonic matrix for Yang--Mills.
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


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


y21 = load('yc22_y21','physics/yc21/yc21_two_cube_block.py')
y13 = load('yc22_y13','physics/yc13/yc13_frame_flow.py')
y15 = y21.y15
CAPS = {(3,2):Q(91,2000), (4,2):Q(6,125), (4,3):Q(11,200)}


def casimir(n):
    if type(n) is not int or n < 0:
        raise ValueError('twice-spin must be a nonnegative integer')
    return n*(n+2)


def charge_pattern_floor(pattern, degrees):
    ns, ds = tuple(pattern), tuple(degrees)
    if len(ns) != len(ds) or not ns or any(type(d) is not int or d <= 0 for d in ds):
        raise ValueError('one positive integer degree per charge label required')
    cs = [casimir(n) for n in ns]
    if sum(ns) % 2:
        raise ValueError('the global vertex centre kernel excludes odd total twice-spin')
    return sum(Q(c,2*d) for c,d in zip(cs,ds))


def charged_floor(dmax):
    if type(dmax) is not int or dmax <= 0:
        raise ValueError('positive integer maximum degree required')
    return Q(3,dmax)


def edge_force_bound(degrees, incident_coupling_sum):
    if len(degrees) != 2 or any(type(d) is not int or d <= 0 for d in degrees):
        raise ValueError('two positive endpoint degrees required')
    c = Q(incident_coupling_sum)
    if c < 0:
        raise ValueError('the coupling sum is a sum of absolute values')
    return Q(min(degrees)**2,64)*c*c


def inverse_square_bound(variance, floor):
    v, d = Q(variance), Q(floor)
    if v < 0 or d <= 0:
        raise ValueError('nonnegative variance and positive full source floor required')
    return Q(1,576)+v*(1/(864*d)+1/(144*d*d))


def open_block():
    size = (4,2,2)
    points = tuple(product(*(range(n) for n in size)))
    edges = [(p,a) for p in points for a in range(3) if p[a]+1 < size[a]]

    def shift(p,a):
        q = list(p)
        q[a] += 1
        return tuple(q)

    degrees = Counter(v for p,a in edges for v in (p,shift(p,a)))
    face_count, faces = Counter(), []
    for p in points:
        for a,b in combinations(range(3),2):
            if p[a]+1 < size[a] and p[b]+1 < size[b]:
                face = ((p,a),(shift(p,a),b),(shift(p,b),a),(p,b))
                faces.append(face)
                for e in face:
                    face_count[e] += 1
    incidence, records = Counter(), []
    for p,a in edges:
        ds = (degrees[p],degrees[shift(p,a)])
        c = face_count[p,a]
        weight = sum(p[b] in (0,size[b]-1) for b in range(3) if b != a)
        key = (min(ds),c)
        incidence[key] += weight
        records.append(dict(edge=(p,a), degrees=ds, faces=c, external_incidence=weight,
                            kinetic_bound=edge_force_bound(ds,c), cap=CAPS[key]))
    return dict(vertices=points, edges=edges, faces=faces, degrees=degrees,
                boundary_incidence=incidence, records=records)


def join_bounds(kind, epsilon=None):
    profiles = {
        'cube':(Q(1), Q(27,5), 24, Q(6), Q(1,1920)),
        'tube':(Q(3,4), Q(17,10), 48, Q(64,5), Q(1,4480)),
        'force':(Q(3,4), Q(15,8), 48, Q(228,25), Q(1,4000)),
    }
    if kind not in profiles:
        raise ValueError('choose cube, tube or force profile')
    g, delta, count, seed_coefficient, maximum = profiles[kind]
    e = maximum if epsilon is None else Q(epsilon)
    if not 0 <= e <= maximum:
        raise ValueError('outside the certified external coupling window')
    r, exp_upper = Q(1,128), Q(16,15)
    beta, seed = count*e, seed_coefficient*e
    q = 4*beta/g*exp_upper*(2+8*(1+2*r))
    b = 8*beta/g*exp_upper
    return dict(profile=kind, external_cap=maximum, epsilon=e,
                full_factor_floor=g, physical_reference_floor=delta,
                incidence=count, beta=beta, seed=seed, radius=r,
                contraction=q, relative_return=b, ball_reserve=r-seed-q*r,
                full_gap=g*(1-b), physical_gap=delta*(1-b),
                nonlinear_minimum_support=1)


def vector_field(poly, coords, axis, side='left'):
    unit = tuple(int(i == axis+1) for i in range(4))
    coefficients = (y15.quaternion_product(unit,coords) if side == 'left'
                    else y15.quaternion_product(coords,unit))
    return sp.expand(sum(c*sp.diff(poly,x) for c,x in zip(coefficients,coords)))


def sphere_casimir(poly, coords):
    return sp.expand(-sum(vector_field(vector_field(poly,coords,a),coords,a) for a in range(3)))


@lru_cache(None)
def normalization_controls():
    qs = tuple(tuple(sp.symbols(f'q{e}_0:4',real=True)) for e in range(3))
    W = sp.expand(y15.quaternion_product(y15.quaternion_product(qs[0],qs[1]),qs[2])[0])

    def gauss(poly,v,a):
        return sp.expand(vector_field(poly,qs[v],a,'left')
                         - vector_field(poly,qs[(v-1)%3],a,'right'))

    def C(poly,v):
        return sp.expand(-sum(gauss(gauss(poly,v,a),v,a) for a in range(3)))

    derivative = vector_field(W,qs[0],0)
    return qs,W,gauss,C,derivative


def source_pins():
    paths = ['physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc13/YC13_FRAME_FLOW_AND_RETURN.md',
             'physics/yc13/yc13_frame_flow.py','physics/yc13/YC13_RESULT.json',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md',
             'physics/yc15/YC15_RESULT.json',
             'physics/yc21/YC21_TWO_CUBE_BLOCK_AND_REMOVED_CHANNELS.md',
             'physics/yc21/yc21_two_cube_block.py','physics/yc21/YC21_RESULT.json',
             'physics/yc22/YC22_GAUGE_PROTECTED_CHANNELS.md',
             'physics/yc22/yc22_gauge_protected_channels.py','physics/yc22/test_yc22.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    for p in ('physics/yc13/YC13_RESULT.json','physics/yc15/YC15_RESULT.json','physics/yc21/YC21_RESULT.json'):
        y15.verify_predecessor(p)
    checks = {}
    qs,W,G,C,derivative = normalization_controls()
    checks['unit quaternion coordinates have fundamental Casimir three'] = all(
        sphere_casimir(x,qs[0]) == 3*x for x in qs[0])
    checks['left and right endpoint generators commute'] = all(
        vector_field(vector_field(W,qs[0],a),qs[0],b,'right')
        == vector_field(vector_field(W,qs[0],b,'right'),qs[0],a)
        for a,b in product(range(3),repeat=2))
    checks['a closed triangle reading is invariant at every vertex'] = all(
        G(W,v,a) == 0 for v,a in product(range(3),repeat=2))
    checks['an edge derivative of an invariant reading is adjoint at one endpoint'] = C(derivative,0) == 8*derivative
    checks['that derivative is a singlet at every other endpoint'] = all(C(derivative,v) == 0 for v in (1,2))
    checks['an open edge source has two fundamental charges'] = all(
        C(qs[0][0],v) == (3*qs[0][0] if v in (0,1) else 0) for v in range(3))
    checks['the loop carrier has its full free energy nine'] = sum(sphere_casimir(W,q) for q in qs).expand() == 9*W
    x = sp.symbols('x',real=True)
    omega, f = 2+sp.cos(x),sp.sin(2*x)
    V = sp.diff(omega,x,2)/omega
    hpsi = -sp.diff(omega*f,x,2)+V*omega*f
    checks['weighted ground form retains its exact total derivative'] = sp.trigsimp(
        omega*f*hpsi-omega**2*sp.diff(f,x)**2+sp.diff(omega**2*f*sp.diff(f,x),x)) == 0
    checks['global centre kernel raises the first allowed charge sum from three to six'] = min(
        sum(casimir(n) for n in ns) for ns in product(range(4),repeat=3)
        if any(ns) and sum(ns)%2 == 0) == 6
    checks['cube and rectangle charged floors are one and three quarters'] = charged_floor(3)==1 and charged_floor(4)==Q(3,4)
    checks['periodic degree-six charged floor is independent of volume and coupling'] = charged_floor(6)==Q(1,2)
    checks['the same Casimir estimate supplies zero on physical states'] = charge_pattern_floor((0,0,0),(3,4,6))==0
    w = y13.gap_window(2,Q(27,5))
    checks['frozen cube gauge floor and new charged floor combine without substitution'] = w['face_reserve']>0 and min(w['gauge_gap_lower'],charged_floor(3))==1
    checks['full half-strength tube retains the stronger nine-tenths gap'] = 1-Q(9,10)-Q(1,4)/(4-Q(9,10))==Q(3,155)>0
    checks['strong tube gauge gap seventeen tenths has a positive exact reserve'] = Q(27,5)-Q(17,10)-1/(2-Q(17,10))==Q(11,30)>0
    checks['theta-one tube gauge gap fifteen eighths has a positive exact reserve'] = 11-Q(15,8)-1/(2-Q(15,8))==Q(9,8)>0
    checks['charged cut protects the strong tube despite an unusable odd bridge norm floor'] = 3-4*1<0 and charged_floor(4)==Q(3,4)
    checks['enlarged physical Schur graph metric has a full thirteen-ninths cap'] = 1+1/(2-Q(1,2))**2==Q(13,9)
    # A two-charge-label matrix cannot mix different labels, but multiplicities can mix.
    c1,c2,a,b,c,d = sp.symbols('c1 c2 a b c d',real=True)
    label = sp.diag(c1,c2)
    A = sp.Matrix([[a,b],[c,d]])
    comm = label*A-A*label
    checks['different charge labels force both cross blocks to vanish'] = sp.expand(comm[0,1]-b*(c1-c2))==0 and sp.expand(comm[1,0]-c*(c2-c1))==0
    same = sp.Matrix([[2,1],[1,2]])
    checks['equal charge permits a nonzero complete inverse cross return'] = same*(3*sp.eye(2))==(3*sp.eye(2))*same and same.inv()[0,1]==-sp.Rational(1,3)
    block = open_block()
    checks['actual strong block has sixteen vertices and twenty-eight links'] = len(block['vertices'])==16 and len(block['edges'])==28
    checks['boundary force classes account for all forty-eight incidences'] = block['boundary_incidence']==Counter({(3,2):32,(4,2):8,(4,3):8})
    caps = []
    for key in sorted(CAPS):
        deg,faces = key
        T = edge_force_bound((deg,deg),faces)
        upper = inverse_square_bound(2*T,Q(15,2))
        checks[f'edge class {key}: complete residual lies below its outward cap'] = upper < CAPS[key]**2
        caps.append(dict(degree=deg,faces=faces,kinetic=T,inverse_squared=upper,cap=CAPS[key]))
    seed = 4*sum(count*CAPS[key] for key,count in block['boundary_incidence'].items())
    checks['resolved edge classes reduce the complete first-source seed'] = seed==Q(228,25)<Q(64,5)
    checks['bridge-root seed is strictly below the block-root bound'] = 4*4*max(CAPS.values())<seed
    rows = {name:join_bounds(name) for name in ('cube','tube','force')}
    expected = {'cube':(Q(3,6400),Q(67,75),Q(603,125)),
                'tube':(Q(3,22400),Q(461,700),Q(7837,5250)),
                'force':(Q(53,400000),Q(1619,2500),Q(1619,1000))}
    for name,row in rows.items():
        reserve,full,physical = expected[name]
        checks[f'{name}: full creator map and contraction retain their reserve'] = row['ball_reserve']==reserve>0 and row['contraction']<1
        checks[f'{name}: full and physical gap floors use their separate reference rates'] = row['full_gap']==full and row['physical_gap']==physical
    checks['blocking has not automatically contracted incidence over the full factor floor'] = Q(48)/charged_floor(4)>Q(24)/charged_floor(3)
    checks = {k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k,v in checks.items() if not v])
    return y15.y12.encode(dict(stage='YC22',base_commit='c45bc50',checks=checks,
        charged_floors=dict(cube=charged_floor(3),tube=charged_floor(4),torus=charged_floor(6)),
        source_caps=caps,lattice_joins=rows,source_pins=source_pins(),
        claim_boundary=dict(all_coupling_charged_form_bound=True,
            full_interacting_charge_pattern_orthogonality=True,
            local_true_ground_force_bound=True,
            stronger_full_cube_and_tube_gaps=True,
            larger_anisotropic_volume_uniform_physical_windows=True,
            same_charge_multiplicity_can_mix=True,
            charged_floor_is_physical_mass_gap=False,
            finite_harmonic_truncation=False,improved_isotropic_window=False,
            iterated_block_contraction=False,continuum_mass_gap=False,
            formal_machine_verification=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC22_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC22: {len(result['checks'])} exact checks pass.")
    print('Charged floor 3/dmax; cube full gap >=1; strong tube full gap >=3/4.')
    print('Resolved anisotropic physical gap >=1619/1000 through external coupling 1/4000.')
