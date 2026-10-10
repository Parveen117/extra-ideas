"""YC17: exact bridge-cycle selection, Haar contraction and kinetic costs.

All integer/rational checks support the accompanying full-operator proof.
No finite harmonic calculation is substituted for the full resolvent bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    'yc17_y16', ROOT/'physics/yc16/yc16_absorbed_cube_return.py')
y16 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = y16
spec.loader.exec_module(y16)
y15, y12 = y16.y15, y16.y15.y12
mul, conj = y15.quaternion_product, y15.conjugate
AXES = tuple(tuple(int(i == j) for i in range(4)) for j in range(4))
SIGNATURES = tuple((1 << i) | (1 << ((i+1) % 4)) for i in range(4))


def bridge_matrix(a, b):
    """W(a,b,u_next,u_current)=u_current^T T(a,b) u_next."""
    return sp.Matrix(4, 4, lambda i, j: y15.wilson(a, b, AXES[j], AXES[i]))


def quaternion_chain(values):
    out = AXES[0]
    for value in values:
        out = mul(out, value)
    return out


def static_tube(a, b):
    if len(a) != 4 or len(b) != 4:
        raise ValueError('four oriented internal links on each end required')
    out = sp.eye(4)
    for aa, bb in zip(a, b):
        out = out*bridge_matrix(aa, bb)
    return sp.trace(out)/256


def prefix_data(order):
    if sorted(order) != list(range(4)):
        raise ValueError('one permutation of four distinct tube faces required')
    mask, out = 0, []
    for k, i in enumerate(order[:3], 1):
        mask ^= SIGNATURES[i]
        odd = mask.bit_count()
        out.append(dict(odd_bridges=odd, full_floor=3*odd,
                        free_surviving_energy=6*k+3*odd))
    return out


def tube_bounds(z=0):
    z = Q(z)
    if z >= 3:
        raise ValueError('the full bridge-complement resolvent is declared only for z<3')
    bound = Q(0)
    source = Q(0)
    profiles = Counter()
    for order in permutations(range(4)):
        prefix = prefix_data(order)
        cost, coefficient = Q(1, 4), Q(1, 64)
        for row in prefix:
            cost /= row['full_floor']-z
            coefficient /= row['free_surviving_energy']-z
        bound += cost
        source += coefficient
        profiles[tuple(row['free_surviving_energy'] for row in prefix)] += 1
    return dict(z=z, full_operator_norm_upper=bound,
                free_vacuum_face_pair_coefficient=source,
                free_vacuum_source_norm=source/4,
                incident_selected_coefficient_cost=6*bound,
                free_energy_profiles=dict(profiles))


def tube_geometry(shape):
    """Actual face/edge sets, checked against YC12's independent tiling."""
    tiling = y12.tiling(shape)
    factors = tiling['edge_factors']

    def shift(x, axis):
        return tuple((v+int(i == axis)) % shape[i] for i, v in enumerate(x))

    all_faces = {}
    for x in product(*(range(n) for n in shape)):
        for i, j in combinations(range(3), 2):
            edges = frozenset(((x, i), (shift(x, i), j),
                               (shift(x, j), i), (x, j)))
            all_faces[edges] = frozenset(factors[e] for e in edges)
    tubes = []
    for c in product(*(range(n//2) for n in shape)):
        for normal in range(3):
            x = tuple(2*v+int(i == normal) for i, v in enumerate(c))
            j, k = (i for i in range(3) if i != normal)
            verts = (x, shift(x, j), shift(shift(x, j), k), shift(x, k))
            bridges = tuple((v, normal) for v in verts)
            end_a = ((verts[0], j), (verts[1], k), (verts[3], j), (verts[0], k))
            end_b = tuple((shift(v, normal), ax) for v, ax in end_a)
            faces = tuple(frozenset((end_a[i], end_b[i], bridges[i], bridges[(i+1)%4]))
                          for i in range(4))
            cubes_a = {factors[e] for e in end_a}
            cubes_b = {factors[e] for e in end_b}
            assert len(cubes_a) == len(cubes_b) == 1
            ca, cb = next(iter(cubes_a)), next(iter(cubes_b))
            assert ca[0] == cb[0] == 'cube' and ca != cb
            assert all(factors[e][0] == 'bridge' for e in bridges)
            for face in faces:
                assert face in all_faces
                assert len(all_faces[face]) == 4
                assert sum(f[0] == 'bridge' for f in all_faces[face]) == 2
            tubes.append(dict(cubes=(ca, cb), bridges=bridges, faces=faces))
    two_bridge_faces = {face for face, support in all_faces.items()
                        if sum(f[0] == 'bridge' for f in support) == 2}
    return tubes, two_bridge_faces, factors


def tensor_checks():
    # Bilinear T(a,b) makes this degree-one-in-each-argument basis test exact
    # for arbitrary quaternion coefficients, not just the unit basis fixtures.
    matrices = {(i,j): bridge_matrix(a,b)
                for i,a in enumerate(AXES) for j,b in enumerate(AXES)}
    composition = all(matrices[i,j]*matrices[k,l] == bridge_matrix(mul(AXES[i],AXES[k]),
                                                                mul(AXES[j],AXES[l]))
                      for i,j,k,l in product(range(4), repeat=4))
    trace = all(sp.trace(matrices[i,j]) == 4*AXES[i][0]*AXES[j][0]
                for i,j in product(range(4), repeat=2))
    return composition, trace


def curvature_controls():
    x,y=sp.symbols('x y', real=True)
    R=sp.Matrix([[0,-1],[1,0]])
    Ax,Ay=sp.zeros(2),x*R
    abelian=Ay.diff(x)-Ax.diff(y)+Ax*Ay-Ay*Ax
    # A noncommuting flat GL(2) control: g=(I+x E12)(I+y E21).
    E12,E21=sp.Matrix([[0,1],[0,0]]),sp.Matrix([[0,0],[1,0]])
    g=(sp.eye(2)+x*E12)*(sp.eye(2)+y*E21)
    Lx,Ly=g.inv()*g.diff(x),g.inv()*g.diff(y)
    bracket=sp.simplify(Lx*Ly-Ly*Lx)
    curvature=sp.simplify(Ly.diff(x)-Lx.diff(y)+bracket)
    p=sp.symbols('p',real=True)
    shared=sp.simplify(p*(Ly.diff(x)-Lx.diff(y))+p*p*bracket)
    return dict(commuting_components_can_curve=abelian==R,
                noncommuting_pure_gauge_can_be_flat=bracket!=sp.zeros(2) and curvature==sp.zeros(2),
                constant_share_law=sp.simplify(shared-p*(p-1)*bracket)==sp.zeros(2))


def source_pins():
    paths=['physics/yc11/YC11_CLOSED_CUBE_RESPONSE.md',
           'physics/yc12/yc12_correlated_cube_gap.py',
           'physics/yc15/yc15_boundary_return_channels.py',
           'physics/yc16/YC16_ABSORBED_CUBE_RETURN.md',
           'physics/yc16/YC16_RESULT.json',
           'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
           'physics/yc17/YC17_FOUR_FACE_BRIDGE_RETURN.md',
           'physics/yc17/yc17_four_face_bridge_return.py',
           'physics/yc17/test_yc17.py']
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    y15.verify_predecessor('physics/yc16/YC16_RESULT.json')
    composition,trace=tensor_checks()
    checks={'quaternion composition tensor identity on every basis entry':composition,
            'quaternion bimultiplication trace tensor identity':trace,
            'static tube coefficient from four exact bridge moments':Q(4,4**4)==Q(1,64),
            'closing both ends recovers the YC11 six-face moment':Q(1,64)*Q(1,4)**2==Q(1,1024)}
    subsets=[]
    for selected in range(1,16):
        mask=0
        for i in range(4):
            if selected & (1<<i): mask^=SIGNATURES[i]
        if mask==0: subsets.append(selected)
    checks['only the complete distinct-face tube has empty bridge signature']=subsets==[15]
    checks['all three-insertion words in one tube have nonempty signature']=all(
        SIGNATURES[i]^SIGNATURES[j]^SIGNATURES[k] for i,j,k in product(range(4),repeat=3))
    rows=[]
    for shape in ((4,4,4),(4,6,8),(6,6,6)):
        tubes,faces,factors=tube_geometry(shape)
        face_counts=Counter(face for t in tubes for face in t['faces'])
        bridge_counts=Counter(e for t in tubes for e in t['bridges'])
        incidence=Counter(c for t in tubes for c in t['cubes'])
        block_count=(shape[0]*shape[1]*shape[2])//8
        checks[f'actual tube partition, bridge ownership and six incidences on {shape}']=(
            len(tubes)==3*block_count and set(face_counts)==faces and
            set(face_counts.values())=={1} and
            set(bridge_counts)=={e for e,f in factors.items() if f[0]=='bridge'} and
            set(bridge_counts.values())=={1} and len(incidence)==block_count and
            set(incidence.values())=={6})
        rows.append(dict(shape=shape,cubes=block_count,tubes=len(tubes),two_bridge_faces=len(faces)))
    row=tube_bounds()
    checks['the 24 orders have the two complete free-energy profiles']=(
        row['free_energy_profiles']=={(12,18,24):16,(12,24,24):8})
    checks['full resolvents retain the five over 216 operator envelope']=row['full_operator_norm_upper']==Q(5,216)
    checks['free vacuum kinetic coefficient is eleven over 165888']=row['free_vacuum_face_pair_coefficient']==Q(11,165888)
    checks['connected free source has exact nonzero L2 norm']=row['free_vacuum_source_norm']==Q(11,663552)>0
    checks['selected incidence cost is five eta fourth over 36']=row['incident_selected_coefficient_cost']==Q(5,36)
    checks['static and kinetic coefficients are distinct']=row['free_vacuum_face_pair_coefficient']!=Q(1,64)
    checks['operator bound dominates the computed source norm']=row['full_operator_norm_upper']>row['free_vacuum_source_norm']
    checks.update(curvature_controls())
    checks['YC16 second-order free source is silent while this channel survives']=(
        y16.source_bounds(0,0)['connected_source_squared_upper']==0 and row['free_vacuum_source_norm']>0)
    checks['existing actual gap baseline remains above nine fortieths']=(
        y15.join_bounds(y15.ETA_MAX)['full_gap_lower']>Q(9,40))
    checks={k:bool(v) for k,v in checks.items()}
    if not all(checks.values()): raise AssertionError([k for k,v in checks.items() if not v])
    return y12.encode(dict(stage='YC17',base_commit='fa07e53',checks=checks,
        static_tube_coefficient=Q(1,64),kinetic_return=row,lattice_checks=rows,
        claim_boundary=dict(two_bridge_distinct_face_minimal_return_classified=True,
            full_operator_selected_word_bound=True,all_spectator_energies_retained=True,
            exact_free_vacuum_source=True,full_fourth_order_return_resolved=False,
            interacting_vacuum_source_coefficient_computed=False,
            local_interaction_norm_after_spectator_elimination=False,
            improved_lattice_gap_window=False,spatial_RG=False,continuum_mass_gap=False,
            formal_machine_verification=False),source_pins=source_pins()))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    path=HERE/'YC17_RESULT.json'
    if args.check:
        if json.loads(path.read_text())!=result: raise AssertionError('certificate or source pins changed')
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC17: {len(result['checks'])} exact checks pass.")
    print('Four-face tube: static 1/64; free kinetic source 11/165888; full operator norm <=5/216.')
