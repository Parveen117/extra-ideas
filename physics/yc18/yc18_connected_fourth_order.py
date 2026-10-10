"""YC18: complete low-order bridge selection and connected creator recursion.

Exact finite checks support the full-space and uniform analytic proofs in the
note. Finite matrix controls do not replace the infinite harmonic inverses.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('yc18_y17',ROOT/'physics/yc17/yc17_four_face_bridge_return.py')
y17=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=y17
spec.loader.exec_module(y17)
y15,y12=y17.y15,y17.y12
RHO=Q(1,4320)
BALL=Q(1,128)


def lattice_audit(shape):
    faces=y12.tiling(shape)['external']
    bridges=sorted({f for support in faces for f in support if f[0]=='bridge'})
    ids={b:i for i,b in enumerate(bridges)}
    masks=[sum(1<<ids[f] for f in support if f[0]=='bridge') for support in faces]
    if len(set(masks))!=len(masks):raise AssertionError('distinct faces share a bridge signature')
    lookup={mask:i for i,mask in enumerate(masks)}
    pairs=defaultdict(list)
    triples=set()
    max_shared=0
    for i,j in combinations(range(len(faces)),2):
        overlap=(masks[i]&masks[j]).bit_count()
        max_shared=max(max_shared,overlap)
        value=masks[i]^masks[j]
        pairs[value].append((i,j))
        if value in lookup:triples.add(tuple(sorted((i,j,lookup[value]))))
    quads=set()
    for bucket in pairs.values():
        for a,b in combinations(bucket,2):
            indices=tuple(sorted(a+b))
            assert len(set(indices))==4
            quads.add(indices)
    actual_masks={frozenset(masks[i] for i in quad) for quad in quads}
    tubes,_,_=y17.tube_geometry(shape)
    expected_masks=set()
    for tube in tubes:
        signature_set=[]
        for face in tube['faces']:
            signature_set.append(sum(1<<ids[('bridge',edge)] for edge in face
                                     if ('bridge',edge) in ids))
        expected_masks.add(frozenset(signature_set))
    return dict(shape=shape,external_faces=len(faces),max_shared_bridges=max_shared,
                third_order_distinct_returns=len(triples),fourth_order_distinct_returns=len(quads),
                all_fourth_distinct_are_exactly_tubes=actual_masks==expected_masks,
                all_fourth_distinct_have_only_two_bridge_faces=all(
                    all(masks[i].bit_count()==2 for i in quad) for quad in quads))


def self_four_bound(m):
    if m not in (2,4):raise ValueError('external faces have two or four bridges')
    return Q(1,288*m*m)


def paired_bound(m,n,shared):
    if m not in (2,4) or n not in (2,4) or shared not in (0,1):
        raise ValueError('ordinary distinct external face types and overlap required')
    s=m+n-2*shared
    return (Q(1,144*m*n) if shared else Q(0))+Q((m+n)**2,108*s*m*m*n*n)


def paired_bound_by_words(m,n,shared):
    paired_bound(m,n,shared)  # Validate the physical parameter range.
    masks=[(1<<m)-1, ((1<<n)-1)<<(m-shared)]
    bound=Q(0)
    kept=0
    for order in sorted(set(permutations((0,0,1,1)))):
        if not shared and order[0]==order[1]:continue
        mask=0
        word=Q(1,4)
        for i in order[:3]:
            mask^=masks[i]
            floor=3*mask.bit_count() if mask else 8
            word/=floor
        bound+=word
        kept+=1
    return bound,kept


def tail_bound(magnitude,order=4):
    magnitude=Q(magnitude)
    if not 0<=magnitude<RHO or type(order) is not int or order<0:
        raise ValueError('nonnegative magnitude below the analytic radius and integer order required')
    x=magnitude/RHO
    return BALL*x**(order+1)/(1-x)


def qubit_reference(energies):
    if any(e<=0 for e in energies):raise ValueError('positive factor gaps required')
    dim=1<<len(energies)
    diagonal=[sum(sp.Rational(e) for i,e in enumerate(energies) if mask&(1<<i))
              for mask in range(dim)]
    omega=sp.zeros(dim,1)
    omega[0]=1
    inverse=sp.diag(0,*(1/e for e in diagonal[1:]))
    return sp.diag(*diagonal),inverse,omega


def creator(vector):
    dim=len(vector)
    if dim<2 or dim&(dim-1) or vector[0]!=0:
        raise ValueError('excited vector in an ordered binary factor carrier required')
    out=sp.zeros(dim)
    for mask in range(1,dim):
        if vector[mask]==0:continue
        for start in range(dim):
            if not start&mask:out[start|mask,start]+=vector[mask]
    return out


def factor_matrix(matrix,site,count):
    if not 0<=site<count:raise ValueError('invalid binary factor site')
    factors=[matrix if i==site else sp.eye(2) for i in reversed(range(count))]
    return sp.kronecker_product(*factors)


def commutator_coefficients(energies,V):
    _,inverse,omega=qubit_reference(energies)
    if V.shape!=(len(omega),len(omega)):raise ValueError('wrong interaction carrier')
    def B(*vectors):
        operator=V
        for vector in vectors:
            C=creator(vector)
            operator=operator*C-C*operator
        return sp.simplify(inverse*operator*omega)
    c1=B()
    c2=B(c1)
    c3=B(c2)+B(c1,c1)/2
    c4=B(c3)+B(c1,c2)+B(c1,c1,c1)/6
    return [sp.simplify(c) for c in (c1,c2,c3,c4)]


def rayleigh_coefficients(energies,V,order=4):
    """Independent intermediate-normalized recursion for H0-zeta V."""
    _,inverse,omega=qubit_reference(energies)
    psi=[omega]
    energy=[sp.Integer(0)]
    for n in range(1,order+1):
        energy.append(sp.expand(-(omega.T*V*psi[n-1])[0]))
        rhs=V*psi[n-1]
        for k in range(1,n+1):rhs+=energy[k]*psi[n-k]
        psi.append(sp.simplify(inverse*rhs))
    return psi,energy


def compositions(total,length):
    if length==1:
        if total>=1:yield (total,)
        return
    for first in range(1,total-length+2):
        for rest in compositions(total-first,length-1):yield (first,)+rest


def creation_log(psi):
    """Independent log(1+Y) in the commuting nilpotent creation algebra."""
    omega=psi[0]
    Y=[None]+[creator(v) for v in psi[1:]]
    result=[]
    for n in range(1,len(psi)):
        C=sp.zeros(len(omega))
        for length in range(1,n+1):
            for indices in compositions(n,length):
                term=sp.eye(len(omega))
                for i in indices:term=term*Y[i]
                C+=sp.Rational((-1)**(length+1),length)*term
        result.append(sp.simplify(C*omega))
    return result


def chain_fixture():
    X=sp.Matrix([[0,1],[1,0]])
    Z=sp.diag(1,-1)
    x=[factor_matrix(X,i,3) for i in range(3)]
    z1=factor_matrix(Z,1,3)
    V=x[0]*x[1]+z1*x[2]/2+x[1]*x[2]/3
    return (2,3,5),V


def disconnected_control():
    a,b=sp.symbols('a b',real=True)
    X=sp.Matrix([[0,1],[1,0]])
    V=(a*factor_matrix(X,0,2)+b*factor_matrix(X,1,2))/2
    _,D,omega=qubit_reference((12,12))
    second=(omega.T*V*D*V*omega)[0]
    derivative=(omega.T*V*D*D*V*omega)[0]
    fourth=(omega.T*V*D*V*D*V*D*V*omega)[0]
    cross=lambda value:sp.expand(value).coeff(a,2).coeff(b,2)
    c=commutator_coefficients((12,12),V)
    return dict(raw_fourth_cross=cross(fourth),energy_argument_cross=cross(second*derivative),
                net_energy_cross=cross(second*derivative-fourth),
                all_mixed_creator_components_zero=all(vector[3]==0 for vector in c))


def source_pins():
    paths=['physics/yc9/YC9_SHARED_PLANE_GAP.md',
           'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
           'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md',
           'physics/yc15/yc15_boundary_return_channels.py',
           'physics/yc15/YC15_RESULT.json',
           'physics/yc17/YC17_FOUR_FACE_BRIDGE_RETURN.md',
           'physics/yc17/yc17_four_face_bridge_return.py','physics/yc17/YC17_RESULT.json',
           'physics/yc18/YC18_CONNECTED_FOURTH_ORDER.md',
           'physics/yc18/yc18_connected_fourth_order.py','physics/yc18/test_yc18.py']
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    y15.verify_predecessor('physics/yc17/YC17_RESULT.json')
    y15.verify_predecessor('physics/yc15/YC15_RESULT.json')
    checks={}
    audits=[]
    for shape in ((4,4,4),(4,6,8),(6,6,6)):
        audit=lattice_audit(shape)
        audits.append(audit)
        checks[f'complete distinct low-order bridge classification on {shape}']=(
            audit['max_shared_bridges']==1 and audit['third_order_distinct_returns']==0 and
            audit['all_fourth_distinct_are_exactly_tubes'] and
            audit['all_fourth_distinct_have_only_two_bridge_faces'])
    pair_rows=[]
    expected={(2,2,0):Q(1,432),(2,2,1):Q(11,1728),
              (2,4,0):Q(1,1152),(2,4,1):Q(5,2304),
              (4,4,0):Q(1,3456),(4,4,1):Q(17,20736)}
    for (m,n,r),value in expected.items():
        by_words,count=paired_bound_by_words(m,n,r)
        checks[f'complete repeated pair bound for {m},{n},{r}']=(
            paired_bound(m,n,r)==by_words==value and count==(6 if r else 4))
        pair_rows.append(dict(m=m,n=n,shared_bridges=r,norm_upper=value,nonzero_orders_upper=count))
    checks['even bridge floor gives both fourth self-return bounds']=(
        self_four_bound(2)==Q(1,1152) and self_four_bound(4)==Q(1,4608))
    z=sp.symbols('z')
    envelope=4/(6-z)**3+2/((6-z)**2*(12-z))
    checks['spectator derivative retains all three resolvent insertions']=(
        sp.diff(envelope,z).subs(z,0)==Q(29,2592))
    dc=disconnected_control()
    checks['disconnected raw return has a nonzero mixed fourth coefficient']=dc['raw_fourth_cross']==Q(1,13824)
    checks['energy dependence exactly cancels disconnected scalar return']=(
        dc['energy_argument_cross']==Q(1,13824) and dc['net_energy_cross']==0)
    checks['creation logarithm removes disconnected mixed components']=dc['all_mixed_creator_components_zero']
    energies,V=chain_fixture()
    direct=commutator_coefficients(energies,V)
    psi,energy=rayleigh_coefficients(energies,V)
    logged=creation_log(psi)
    checks['four complete commutator coefficients equal independent Rayleigh creation logarithms']=direct==logged
    checks['connected calculation is not the unlogged wavefunction']=any(direct[i]!=psi[i+1] for i in range(4))
    endpoint=y15.join_bounds(RHO)
    checks['uniform complex disk has the pinned contraction and mapping reserve']=(
        endpoint['contraction']==Q(28035072,31987775)<1 and
        endpoint['mapping_reserve']==Q(11417,409443520)>0)
    checks['complete fourth-order tail at half radius is one over 2048']=tail_bound(RHO/2)==Q(1,2048)
    checks['complete fourth-order tail at quarter radius is one over 98304']=tail_bound(RHO/4)==Q(1,98304)
    tube=Q(y17.tube_bounds()['free_vacuum_face_pair_coefficient'])/24
    checks['tube source inverse gives the local connected creation coefficient']=tube==Q(11,3981312)
    checks['tube creator Hilbert and weighted-root norms have distinct normalization']=(
        tube/4==Q(11,15925248) and tube/2==Q(11,7962624))
    checks['actual full gap baseline is preserved']=endpoint['full_gap_lower']>Q(9,40)
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise AssertionError([k for k,v in checks.items() if not v])
    return y12.encode(dict(stage='YC18',base_commit='16bb4a0',checks=checks,
        lattice_audits=audits,paired_word_bounds=pair_rows,
        spectator_tube_relative_constant=Q(29,2592),
        disconnected_control={k:str(v) if isinstance(v,sp.Basic) else v for k,v in dc.items()},
        analytic_radius=RHO,creation_ball=BALL,half_radius_tail=tail_bound(RHO/2),
        quarter_radius_tail=tail_bound(RHO/4),tube_creation_coefficient=tube,
        claim_boundary=dict(all_bridge_word_types_through_four_classified=True,
            each_allowed_fourth_multiset_has_full_operator_bound=True,
            ground_creator_recursion_through_four_complete=True,
            connected_local_ground_coefficients=True,uniform_complete_ground_tail=True,
            all_interacting_cube_coefficients_numerically_evaluated=False,
            full_excitation_effective_operator_locally_closed=False,
            improved_gap_window=False,iterated_spatial_RG=False,continuum_mass_gap=False,
            formal_machine_verification=False),source_pins=source_pins()))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    path=HERE/'YC18_RESULT.json'
    if args.check:
        if json.loads(path.read_text())!=result:raise AssertionError('certificate or source pins changed')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC18: {len(result['checks'])} exact checks pass.")
    print('Complete fourth-order bridge selection; connected local ground recursion; uniform full tail.')
