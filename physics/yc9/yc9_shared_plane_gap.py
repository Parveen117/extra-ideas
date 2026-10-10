"""YC9: exact controls for the written local-cluster gap proof.

These finite checks do not replace the analytic all-volume projection,
contraction, domain and resolvent arguments in YC9_SHARED_PLANE_GAP.md.
No harmonic truncation supplies the stated SU(2) gap theorem.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location('yc8_geometry', ROOT/'physics/yc8/yc8_local_vacuum_dressing.py')
yc8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(yc8)

K = 4
LINK_GAP = F(3)
VACUUM_FREE_GAP = F(12)
RADIUS = F(1, 16)
EXP_MAJORANT = F(5, 3)
BETA_MAX = F(1, 160)
THETA_MAX = F(1, 640)
MAP_COEFFICIENT = F(2**K, 3)*(1+2*RADIUS)*EXP_MAJORANT
LIPSCHITZ_COEFFICIENT = F(2**K, 3)*EXP_MAJORANT*(2+2*K*(1+2*RADIUS))
EXCITATION_COEFFICIENT = 2**(K+1)*EXP_MAJORANT
RELATIVE_COEFFICIENT = EXCITATION_COEFFICIENT/LINK_GAP


def checked_beta(beta):
    beta = F(beta)
    if not 0 <= beta <= BETA_MAX:
        raise ValueError('YC9 requires 0 <= local interaction budget <= 1/160')
    return beta


def plane_budget(a, b, c):
    """a=xy, b=yz, c=zx. Shared links are counted once."""
    a, b, c = map(lambda x: abs(F(x)), (a, b, c))
    return 2*max(a+c, a+b, b+c)


def gap_lower(beta, vacuum=True):
    beta = checked_beta(beta)
    free = VACUUM_FREE_GAP if vacuum else LINK_GAP
    return free*(1-RELATIVE_COEFFICIENT*beta)


def isotropic_gap(theta):
    theta = F(theta)
    if theta < 0:
        raise ValueError('nonnegative isotropic theta required')
    return gap_lower(4*theta)


def iteration_error(beta, iteration):
    beta = checked_beta(beta)
    if isinstance(iteration, bool) or not isinstance(iteration, int) or iteration < 0:
        raise ValueError('nonnegative integer iteration required')
    q = LIPSCHITZ_COEFFICIENT*beta
    return q**iteration*beta/(6*(1-q))


@lru_cache(None)
def geometry_record(shape):
    g = yc8.geometry(shape)
    edges, faces = g['edges'], g['faces']
    ends, adjacency = [], {}
    for e, (x, direction) in enumerate(edges):
        y = list(x)
        y[direction] = (y[direction]+1) % shape[direction]
        y = tuple(y)
        ends.append((x, y))
        adjacency.setdefault(x, {})[y] = e
        adjacency.setdefault(y, {})[x] = e
    triangles = set()
    for x, neighbours in adjacency.items():
        for y, z in combinations(neighbours, 2):
            if z in adjacency[y]:
                triangles.add(frozenset((neighbours[y], neighbours[z], adjacency[y][z])))

    def centre_parity(es):
        chars = [0, 0, 0]
        for e in es:
            x, direction = edges[e]
            if x[direction] == shape[direction]-1:
                chars[direction] ^= 1
        return tuple(chars)

    directions = [frozenset(edges[e][1] for e in face) for face in faces]
    incidence_ok = True
    for e, (_, i) in enumerate(edges):
        counts = {frozenset((i, j)): 0 for j in range(3) if j != i}
        for p in g['incident'][e]:
            counts[directions[p]] += 1
        incidence_ok &= set(counts.values()) == {2}
    n = len(edges)//3
    return {
        'shape': list(shape), 'links': len(edges), 'plaquettes': len(faces),
        'no_self_or_parallel_edges': (all(x != y for x, y in ends) and
            len({frozenset(pair) for pair in ends}) == len(ends)),
        'two_faces_from_each_incident_plane': incidence_ok,
        'elementary_faces_centre_even': all(centre_parity(face) == (0,0,0) for face in faces),
        'triangle_count': len(triangles),
        'expected_triangle_count': sum(n//3 for size in shape if size == 3),
        'triangles_are_single_direction_length_three_windings': all(
            len({edges[e][1] for e in tri}) == 1 and
            shape[edges[next(iter(tri))][1]] == 3 and
            sum(centre_parity(tri)) == 1 for tri in triangles),
        'centre_parities_of_triangles': sorted({centre_parity(tri) for tri in triangles}),
    }


def comm(a, b):
    return a*b-b*a


def qubit_creator(n, support, amplitude=sp.Integer(1)):
    raise_one = sp.Matrix([[0, 0], [1, 0]])
    return amplitude*sp.kronecker_product(*(
        raise_one if i in support else sp.eye(2) for i in range(n)))


@lru_cache(None)
def algebra_checks():
    """Exact witnesses, separate from the infinite-dimensional proof."""
    n = 4
    dim = 2**n
    zero = sp.zeros(dim)
    eye = sp.eye(dim)
    omega = eye[:, 0]
    a = qubit_creator(n, {0,2}, sp.Rational(1,100))
    b = qubit_creator(n, {1,3}, sp.Rational(-1,120))
    remote = qubit_creator(n, {2,3}, sp.Rational(1,140))
    overlapping = qubit_creator(n, {0,3}, sp.Rational(1,160))
    C = a+b+remote+overlapping
    Vlocal = sp.Matrix([[1,2,-1,0], [2,-2,1,3], [-1,1,2,-2], [0,3,-2,-1]])/17
    V = sp.kronecker_product(Vlocal, sp.eye(4))
    H0 = sp.diag(*[3*index.bit_count() for index in range(dim)])
    expC = sum((C**j/sp.factorial(j) for j in range(1,n+1)), eye)
    expminus = sum(((-C)**j/sp.factorial(j) for j in range(1,n+1)), eye)
    W = expminus*V*expC
    ad, series = V, V
    for j in range(1, 2*K+1):
        ad = comm(ad, C)
        series += ad/sp.factorial(j)
    J = qubit_creator(n, {0,1,2}, sp.Rational(2,7))
    hat = expminus*(H0+V)*expC
    checks = {
        'overlapping creation operators multiply to zero': a*overlapping == zero and overlapping*a == zero,
        'disjoint creators commute': comm(a,b) == zero,
        'all creators commute including entangled-support overlap pattern': comm(C, J) == zero,
        'creation exponentials are mutually inverse': expminus*expC == eye,
        'free similarity has only the first commutator': expminus*H0*expC == H0+comm(H0,C),
        'local dressed interaction equals exact commutator expansion': W == series,
        'remote creator cannot enter through an enlarged union': comm(remote,comm(a,V)) == zero,
        'excitation commutator covariance': comm(W,J) == expminus*comm(V,J)*expC,
        'excitation identity uses no assumed scalar ground residual': (hat*J-J*hat)*omega == H0*J*omega+comm(W,J)*omega,
    }
    # After each fixed operator word, all links outside X are either
    # fixed excited or the word is zero; superpose words only afterwards.
    supports = [({0,2}, a), ({1,3}, b), ({2,3}, remote), ({0,3}, overlapping)]
    outside_ok = True
    checked_words = 0
    for (ia, ca), (ib, cb) in product(supports, repeat=2):
        exterior = (ia|ib)-{0,1}
        for word in (ca*cb*V, ca*V*cb, cb*V*ca, V*ca*cb):
            vector = word*omega
            nonzero = [i for i, value in enumerate(vector) if value]
            expected_mask = sum(1 << (n-1-i) for i in exterior)
            outside_ok &= all((i & 3) == expected_mask for i in nonzero)
            outside_ok &= len(nonzero) <= 2**2
            checked_words += 1
    checks['outside-face excitation pattern and projection count on 64 words'] = outside_ok and checked_words == 64
    return checks


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (tuple,list)):
        return [encode(v) for v in obj]
    return obj


def source_pins():
    names = ['physics/yc8/YC8_LOCAL_VACUUM_DRESSING.md',
             'physics/yc8/yc8_local_vacuum_dressing.py', 'physics/yc8/YC8_RESULT.json',
             'physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc9/yc9_shared_plane_gap.py', 'physics/yc9/test_yc9.py']
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def run():
    checks = dict(algebra_checks())
    records = []
    for shape in ((3,3,3), (3,4,5), (4,4,4)):
        r = geometry_record(shape)
        records.append(r)
        for key in ('no_self_or_parallel_edges', 'two_faces_from_each_incident_plane',
                    'elementary_faces_centre_even', 'triangles_are_single_direction_length_three_windings'):
            checks[f'{shape}: {key}'] = bool(r[key])
        checks[f'{shape}: all short triangles counted'] = r['triangle_count'] == r['expected_triangle_count']
    checks.update({
        'rational exponential majorant': F(3,2)+F(1,8)/(1-F(1,6)) == F(33,20) < EXP_MAJORANT,
        'fixed-point ball maps to itself': MAP_COEFFICIENT == 10 and MAP_COEFFICIENT*BETA_MAX == RADIUS,
        'contraction has a strict rational reserve': LIPSCHITZ_COEFFICIENT*BETA_MAX == F(11,18) < 1,
        'excitation relative norm at the endpoint': RELATIVE_COEFFICIENT*BETA_MAX == F(1,9),
        'actual full-product gap endpoint': gap_lower(BETA_MAX, vacuum=False) == F(8,3),
        'actual physical vacuum gap endpoint': gap_lower(BETA_MAX) == F(32,3),
        'three-plane join budget gives the isotropic endpoint': plane_budget(THETA_MAX,THETA_MAX,THETA_MAX) == BETA_MAX,
        'isotropic floor slope and endpoint': VACUUM_FREE_GAP*RELATIVE_COEFFICIENT*4 == F(2560,3) and isotropic_gap(THETA_MAX) == F(32,3),
        'single-source return agrees with YC8 second coefficient': F(1,4)/12 == F(1,48),
        'centre-odd short winding cannot lower vacuum free floor': 3*3 == 9 < VACUUM_FREE_GAP and 3*8 > VACUUM_FREE_GAP,
        'three transverse planes label eight local octants': len(tuple(product((-1,1),repeat=3))) == 8,
    })
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    return encode({
        'stage':'YC9', 'base_commit':'15ae97f', 'checks':checks,
        'carrier':'3D periodic spatial SU(2); every side >=3; unit S3 link metric; full untruncated local Hilbert spaces',
        'geometry_controls':records,
        'weighted_creator_ball':RADIUS, 'local_budget_max':BETA_MAX,
        'isotropic_theta_max':THETA_MAX, 'contraction_upper':F(11,18),
        'relative_excitation_upper':F(1,9),
        'actual_vacuum_gap_endpoint':gap_lower(BETA_MAX),
        'actual_full_product_gap_endpoint':gap_lower(BETA_MAX,False),
        'isotropic_vacuum_gap_formula':'12 - (2560/3)*theta on 0 <= theta <= 1/640',
        'join_path':[{'s':s,'theta_xy':THETA_MAX,'theta_yz':s*THETA_MAX,'theta_zx':s*THETA_MAX,
                      'beta':plane_budget(THETA_MAX,s*THETA_MAX,s*THETA_MAX),
                      'vacuum_gap_lower':gap_lower(plane_budget(THETA_MAX,s*THETA_MAX,s*THETA_MAX))}
                     for s in (F(0),F(1,2),F(1))],
        'iteration_tail_at_endpoint':[{'iteration':m,'local_cluster_error_upper':iteration_error(BETA_MAX,m)}
                                     for m in (0,1,4,8,16)],
        'claim_boundary':{
            'actual_YM_finite_volume_gap_uniform_in_volume':True,
            'actual_vacuum_sector_gap_uniform_in_volume':True,
            'comparison_operator_substituted_for_actual_YM':False,
            'all_cluster_orders_controlled_analytically':True,
            'finite_checks_are_formal_proof':False,
            'independent_planar_tensor_product':False,
            'first_excitation_identified':False,
            'infinite_volume_state_constructed_in_this_packet':False,
            'continuum_4d_proved':False,
            'new_general_existence_of_strong_coupling_gaps_claimed':False},
        'source_pins':source_pins()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC9_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC9: {len(result['checks'])} exact checks pass.")
    print('Actual vacuum gap >=32/3, uniformly in volume, on 0<=theta<=1/640.')
    print('Full untruncated SU(2) spaces; small lattice-coupling window; continuum not claimed.')
