"""YC19: local vacuum-annihilating interactions and full excitation errors.

The full-carrier analytic proofs are in the note. Finite exact controls check
creation subtraction, operator coefficients, transported metrics and budgets.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import product
import json
from math import prod
from pathlib import Path
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('yc19_y18',ROOT/'physics/yc18/yc18_connected_fourth_order.py')
y18=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=y18
spec.loader.exec_module(y18)
y15,y12=y18.y15,y18.y12
RHO=y18.RHO
BALL=y18.BALL
GAP=y15.GAP
LOCAL=Q(130)*RHO
RELATIVE=Q(5537792,31987775)


def vacuum_lift(vector,dims):
    """Full creation lift, including the scalar, on finite controls only."""
    if not dims or any(type(d) is not int or d<2 for d in dims) or len(vector)!=prod(dims):
        raise ValueError('vector and complete factor dimensions must agree')
    out=sp.zeros(prod(dims))
    for row,digits in enumerate(product(*(range(d) for d in dims))):
        value=vector[row]
        if value==0:continue
        factors=[]
        for digit,dim in zip(digits,dims):
            if digit==0:factor=sp.eye(dim)
            else:
                factor=sp.zeros(dim)
                factor[digit,0]=1
            factors.append(factor)
        out+=value*sp.kronecker_product(*factors)
    return out


def normal_order(A,dims):
    if A.shape!=(prod(dims),prod(dims)):raise ValueError('wrong local operator carrier')
    return A-vacuum_lift(A[:,0],dims)


def nilpotent_exp(C):
    if C.rows!=C.cols:raise ValueError('square creator required')
    term=sp.eye(C.rows)
    out=term
    for n in range(1,C.rows+1):
        term=term*C/n
        out+=term
        if term==sp.zeros(C.rows):return out
    raise ValueError('the supplied operator is not nilpotent')


def series_product(a,b,order):
    return [sum((a[k]*b[n-k] for k in range(n+1)),sp.zeros(a[0].rows)) for n in range(order+1)]


def series_exp(C,order):
    dim=C[0].rows
    power=[sp.eye(dim)]+[sp.zeros(dim) for _ in range(order)]
    out=[v.copy() for v in power]
    for k in range(1,order+1):
        power=[v/k for v in series_product(power,C,order)]
        out=[a+b for a,b in zip(out,power)]
    return out


def excitation_coefficients(energies,V):
    c=y18.commutator_coefficients(energies,V)
    c1,c2,c3=(y18.creator(v) for v in c[:3])
    comm=lambda a,b:a*b-b*a
    L=[V,comm(V,c1),comm(V,c2)+comm(comm(V,c1),c1)/2,
       comm(V,c3)+comm(comm(V,c1),c2)+comm(comm(comm(V,c1),c1),c1)/6]
    return [-normal_order(value,(2,)*len(energies)) for value in L]


def similarity_coefficients(energies,V):
    """Independent RS ground, creation logarithm and full operator similarity."""
    H0,_,_=y18.qubit_reference(energies)
    psi,energy=y18.rayleigh_coefficients(energies,V)
    c=y18.creation_log(psi)
    C=[sp.zeros(H0.rows)]+[y18.creator(v) for v in c]
    plus=series_exp(C,4)
    minus=series_exp([-v for v in C],4)
    h=[H0,-V-energy[1]*sp.eye(H0.rows)]+[-energy[n]*sp.eye(H0.rows) for n in range(2,5)]
    full=series_product(series_product(minus,h,4),plus,4)
    return full[1:]


def exact_metric_control():
    c=sp.Matrix([0,sp.Rational(1,5),sp.Rational(1,7),sp.Rational(1,11)])
    C=y18.creator(c)
    S=nilpotent_exp(C)
    omega=sp.eye(4)[:,0]
    psi=S*omega
    P=psi*psi.T/(psi.T*psi)[0]
    D=sp.diag(1,3,5,8)
    H=(sp.eye(4)-P)*D*(sp.eye(4)-P)
    H0=sp.diag(0,2,3,5)
    W=S.inv()*(H-H0)*S
    dressed=S.inv()*H*S
    metric=S.T*S
    A=dressed[1:,1:]
    quotient=metric[1:,1:]-metric[1:,0:1]*metric[0:1,1:]/metric[0,0]
    return dict(S=S,H=H,H0=H0,W=W,dressed=dressed,metric=metric,A=A,quotient=quotient)


def full_local_bound(magnitude):
    magnitude=Q(magnitude)
    if not 0<=magnitude<=RHO:raise ValueError('full reference disk is |zeta|<=1/4320')
    return 130*magnitude


def truncation_bounds(magnitude,order=4,block_size=1,range_base=1):
    magnitude=Q(magnitude)
    range_base=Q(range_base)
    if not 0<=magnitude<RHO or type(order) is not int or order<1:
        raise ValueError('Taylor estimate needs |zeta|<rho and positive integer order')
    if type(block_size) is not int or block_size<1:raise ValueError('positive integer block size required')
    x=magnitude/RHO
    if range_base<1 or range_base*x>=1:raise ValueError('exponential range estimate needs b>=1 and bx<1')
    tail=x**(order+1)/(1-x)
    polynomial=sum((x**n for n in range(1,order+1)),Q(0))
    return dict(magnitude=magnitude,x=x,order=order,
                exact_full_local_upper=full_local_bound(magnitude),
                local_tail_upper=LOCAL*tail,relative_tail_upper=RELATIVE*tail,
                truncated_relative_upper=RELATIVE*polynomial,
                size_weighted_upper=4*LOCAL*x/(1-x)**2,
                exponential_range_upper=LOCAL*range_base*x/(1-range_base*x),
                range_base=range_base,
                regrouped_local_upper=block_size*full_local_bound(magnitude),
                regrouped_tail_upper=block_size*LOCAL*tail)


def resolvent_transfer(magnitude,energy):
    row=truncation_bounds(magnitude)
    energy=Q(energy)
    denominator=1-energy/GAP-row['truncated_relative_upper']
    if energy<0 or denominator<=0:raise ValueError('the declared approximate-resolvent bound is unavailable')
    approximate=1/denominator
    reserve=1-row['relative_tail_upper']*approximate
    return dict(energy=energy,approximate_relative_inverse_upper=approximate,
                error_product=row['relative_tail_upper']*approximate,
                certified=reserve>0,
                exact_relative_inverse_upper=approximate/reserve if reserve>0 else None)


def incidence_budget(terms,partition=None):
    """Nonnegative labelled norm budgets; partition maps factors to blocks."""
    sums={}
    for support,norm in terms:
        norm=Q(norm)
        if norm<0:raise ValueError('operator norm budgets cannot be negative')
        labels=set(support) if partition is None else {partition[i] for i in support}
        for i in labels:sums[i]=sums.get(i,Q(0))+norm
    return max(sums.values(),default=Q(0))


def source_pins():
    paths=['physics/yc9/YC9_SHARED_PLANE_GAP.md',
           'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
           'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md','physics/yc15/YC15_RESULT.json',
           'physics/yc18/YC18_CONNECTED_FOURTH_ORDER.md',
           'physics/yc18/yc18_connected_fourth_order.py','physics/yc18/YC18_RESULT.json',
           'physics/yc19/YC19_LOCAL_EXCITATION_OPERATOR.md',
           'physics/yc19/yc19_local_excitation_operator.py','physics/yc19/test_yc19.py']
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    y15.verify_predecessor('physics/yc18/YC18_RESULT.json')
    y15.verify_predecessor('physics/yc15/YC15_RESULT.json')
    checks={}
    A=sp.Matrix([[1,2,3,4],[2,5,1,0],[3,1,6,2],[4,0,2,7]])/100
    N=normal_order(A,(2,2))
    omega=sp.eye(4)[:,0]
    checks['complete local creation subtraction annihilates the reference']=N*omega==sp.zeros(4,1)
    checks['subtracting only the scalar would leave an excited source']=(A-A[0,0]*sp.eye(4))*omega!=sp.zeros(4,1)
    checks['creation lift is not a whole-support rank-one source replacement']=vacuum_lift(A[:,0],(2,2))!=A[:,0]*omega.T
    checks['normal order has the exact action on every excited column']=all(
        N*sp.eye(4)[:,j]==(A*y18.creator(sp.eye(4)[:,j])-y18.creator(sp.eye(4)[:,j])*A)*omega
        for j in range(1,4))
    energies,V=y18.chain_fixture()
    coefficients=excitation_coefficients(energies,V)
    checks['all four excitation coefficients equal the independent complete similarity']=coefficients==similarity_coefficients(energies,V)
    checks['every excitation coefficient annihilates the reference']=all(K[:,0]==sp.zeros(8,1) for K in coefficients)
    checks['a transformed coefficient need not be Hermitian in the old pairing']=any(K!=K.T for K in coefficients)
    control=exact_metric_control()
    dressed,metric=control['dressed'],control['metric']
    checks['exact similarity equals the locally vacuum-subtracted interaction']=(
        dressed==control['H0']+normal_order(control['W'],(2,2)))
    checks['the full transported positive metric restores adjointness']=dressed.T*metric==metric*dressed
    G=control['quotient']
    checks['the ground quotient keeps its Schur metric']=control['A'].T*G==G*control['A']
    checks['the quotient metric is positive and not the old Euclidean pairing']=(
        all(G[:n,:n].det()>0 for n in (1,2,3)) and G!=sp.eye(3) and control['A']!=control['A'].T)
    checks['local root-count envelope is exactly 130 times the interface amplitude']=(
        5*24*(1+2*BALL)*Q(16,15)==130 and LOCAL==Q(13,432))
    checks['full relative excitation envelope agrees with the pinned lattice proof']=(
        8*24*RHO/GAP*Q(16,15)==RELATIVE==y15.join_bounds(RHO)['relative_return'])
    half=truncation_bounds(RHO/2)
    quarter=truncation_bounds(RHO/4,range_base=2)
    checks['full local operator tail at half radius is thirteen over 6912']=half['local_tail_upper']==Q(13,6912)
    checks['full relative excitation tail at half radius is below 0.01083']=(
        half['relative_tail_upper']==Q(346112,31987775)<Q(1083,100000))
    checks['quarter-radius local tail is thirteen over 331776']=quarter['local_tail_upper']==Q(13,331776)
    checks['quarter-radius interaction has a summable exponential range weight']=quarter['exponential_range_upper']==LOCAL
    checks['quarter-radius interaction has a summable support-size weight']=quarter['size_weighted_upper']==Q(13,243)
    transfer=resolvent_transfer(RHO/2,GAP/2)
    checks['complete error transfers the approximate excitation resolvent']=(
        transfer['error_product']==Q(692224,21604415)<1 and transfer['certified'])
    terms=[({0},Q(1)),({1},Q(1))]
    checks['regrouping can cost the entire block cardinality instead of contracting']=(
        incidence_budget(terms)==1 and incidence_budget(terms,{0:0,1:0})==2)
    checks['existing actual full gap is retained without claiming a larger window']=y15.join_bounds(RHO)['full_gap_lower']>Q(9,40)
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise AssertionError([key for key,value in checks.items() if not value])
    return y12.encode(dict(stage='YC19',base_commit='9a63c80',checks=checks,
        full_local_coefficient=130,local_cauchy_envelope=LOCAL,
        relative_cauchy_envelope=RELATIVE,half_radius=half,quarter_radius=quarter,
        resolvent_control=transfer,
        claim_boundary=dict(full_excitation_operator_local_family=True,
            every_local_term_annihilates_reference=True,uniform_size_and_range_bounds=True,
            complete_operator_truncation_tail=True,transported_full_and_quotient_metrics=True,
            controlled_finite_block_regrouping=True,automatic_budget_contraction=False,
            factorized_block_metric=False,uniform_Hilbert_similarity_condition_number=False,
            improved_gap_window=False,iterated_spatial_RG=False,continuum_mass_gap=False,
            formal_machine_verification=False),source_pins=source_pins()))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    path=HERE/'YC19_RESULT.json'
    if args.check:
        if json.loads(path.read_text())!=result:raise AssertionError('certificate or source pins changed')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC19: {len(result['checks'])} exact checks pass.")
    print('Full local excitation dressing, inherited metric, range decay and complete relative tail certified.')
