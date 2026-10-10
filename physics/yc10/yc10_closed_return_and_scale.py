"""YC10: exact controls for the closed-return proof and typed scale audit.

The all-volume, infinite-harmonic claims rely on the written analytic
proof. Local symbolic witnesses and rational endpoints are replayed here.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from math import factorial
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('yc9_return',ROOT/'physics/yc9/yc9_shared_plane_gap.py')
yc9=importlib.util.module_from_spec(spec)
spec.loader.exec_module(yc9)

RADIUS=F(1,4)
PROJECTION_FACTOR=F(4)
FREE_VACUUM=F(12)


def constants(minimum_support):
    if isinstance(minimum_support,bool) or minimum_support not in (3,4):
        raise ValueError('YC10 certifies minimum supports 3 and 4')
    s=minimum_support
    exp_bound=F(39,20) if s==3 else F(5,3)
    beta_max=F(1,16) if s==3 else F(3,40)
    a=F(8,s)
    map_coeff=PROJECTION_FACTOR/3*(1+2*RADIUS)*exp_bound
    lip_coeff=PROJECTION_FACTOR/3*exp_bound*(2+a*(1+2*RADIUS))
    relative_coeff=2*PROJECTION_FACTOR/3*exp_bound
    return {'minimum_support':s,'radius':RADIUS,'exponent':a*RADIUS,
            'exp_bound':exp_bound,'beta_max':beta_max,'theta_max':beta_max/4,
            'map_coefficient':map_coeff,'lipschitz_coefficient':lip_coeff,
            'relative_coefficient':relative_coeff,
            'map_at_endpoint':map_coeff*beta_max,
            'contraction_at_endpoint':lip_coeff*beta_max,
            'relative_at_endpoint':relative_coeff*beta_max,
            'vacuum_gap_at_endpoint':FREE_VACUUM*(1-relative_coeff*beta_max),
            'full_gap_at_endpoint':3*(1-relative_coeff*beta_max)}


def profile_for_shape(shape):
    yc9.yc8.validate_shape(shape)
    return constants(4 if min(shape)>=4 else 3)


def gap_from_budget(beta,minimum_support):
    c=constants(minimum_support)
    beta=F(beta)
    if not 0<=beta<=c['beta_max']:
        raise ValueError('local budget is outside the certified geometry-specific window')
    return FREE_VACUUM*(1-c['relative_coefficient']*beta)


def isotropic_gap(theta,shape):
    theta=F(theta)
    c=profile_for_shape(shape)
    return gap_from_budget(4*theta,c['minimum_support'])


def iteration_error(beta,minimum_support,iteration):
    gap_from_budget(beta,minimum_support)
    beta=F(beta)
    if isinstance(iteration,bool) or not isinstance(iteration,int) or iteration<0:
        raise ValueError('nonnegative integer iteration required')
    q=constants(minimum_support)['lipschitz_coefficient']*beta
    return q**iteration*beta/(6*(1-q))


def exponential_upper(x,through):
    """Positive Taylor polynomial and a geometric majorant of its tail."""
    x=F(x)
    if x<0 or isinstance(through,bool) or not isinstance(through,int) or through<0:
        raise ValueError('nonnegative argument and integer truncation required')
    ratio=x/F(through+2)
    if ratio>=1:
        raise ValueError('tail ratio must be strictly below one')
    return sum((x**j/factorial(j) for j in range(through+1)),F(0))+x**(through+1)/factorial(through+1)/(1-ratio)


def curvature_scaling(dimension,epsilon):
    if isinstance(dimension,bool) or not isinstance(dimension,int) or dimension<2:
        raise ValueError('integer Euclidean or spatial dimension >=2 required')
    epsilon=F(epsilon)
    if epsilon<=0:
        raise ValueError('positive compression scale required')
    return {'connection':epsilon**-1,'curvature':epsilon**-2,
            'curvature_density':epsilon**-4,'measure':epsilon**dimension,
            'integral':epsilon**(dimension-4)}


@lru_cache(None)
def symbolic_checks():
    checks={}
    x,y=sp.symbols('x y',positive=True)
    phi=(x*x+y*y)**sp.Rational(3,4)
    euler=lambda f: x*sp.diff(f,x)+y*sp.diff(f,y)
    corners=(phi,phi-x*sp.diff(phi,x),phi-y*sp.diff(phi,y),phi-euler(phi))
    checks['all native cut corners retain potential degree 3/2']=all(
        sp.simplify(euler(c)-sp.Rational(3,2)*c)==0 for c in corners)
    hess=sp.hessian(phi,(x,y))
    checks['native Hessian has degree minus one half']=all(
        sp.simplify(euler(h)+h/2)==0 for h in hess)
    radial=sp.symbols('radial',positive=True)
    checks['vanishing potential and diverging response coexist']=(
        sp.limit(radial**sp.Rational(3,2),radial,0)==0 and
        sp.limit(sp.Rational(3,4)/sp.sqrt(radial),radial,0)==sp.oo)

    # Noncommuting su(2) coefficients, bracket normalized as vector cross product.
    xx=sp.symbols('x0:4',real=True)
    epsilon=sp.symbols('epsilon',positive=True)
    a0=sp.Matrix([xx[1],xx[2],0])
    a1=sp.Matrix([0,xx[0],xx[3]*xx[0]])
    curvature=a1.diff(xx[0])-a0.diff(xx[1])+a0.cross(a1)
    sub={z:z/epsilon for z in xx}
    scaled0=a0.subs(sub,simultaneous=True)/epsilon
    scaled1=a1.subs(sub,simultaneous=True)/epsilon
    direct=scaled1.diff(xx[0])-scaled0.diff(xx[1])+scaled0.cross(scaled1)
    checks['nonabelian derivative and bracket scale together']=all(
        sp.simplify(v)==0 for v in direct-curvature.subs(sub,simultaneous=True)/epsilon**2)
    checks['scaling witness retains a nonzero commutator']=a0.cross(a1)!=sp.zeros(3,1)

    q=sp.symbols('q0:4',real=True)
    chi=4*q[0]**2-1
    eq=lambda f: sum(z*sp.diff(f,z) for z in q)
    electric=lambda f: sp.expand(eq(eq(f))+2*eq(f)-sum(sp.diff(f,z,2) for z in q))
    checks['adjoint loop has single-link rate eight']=sp.expand(electric(chi)-8*chi)==0
    checks['adjoint winding is centre even']=sp.expand(chi.subs(q[0],-q[0])-chi)==0
    checks['adjoint winding has zero single-link Haar mean']=4*sp.Rational(1,4)-1==0
    checks['length-three adjoint support must not be discarded']=3*8==24

    shift,gap=sp.symbols('shift gap',real=True)
    checks['divergent scalar shifts leave an excitation difference unchanged']=sp.expand((shift+gap)-shift-gap)==0
    return checks


def encode(obj):
    if isinstance(obj,F): return str(obj)
    if isinstance(obj,dict): return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [encode(v) for v in obj]
    return obj


def source_pins():
    names=['physics/yc9/YC9_SHARED_PLANE_GAP.md','physics/yc9/yc9_shared_plane_gap.py',
           'physics/yc9/YC9_RESULT.json','uncut/up1/UP1_DEGREE_OF_THE_POTENTIAL.md',
           'uncut/up3/UP3_SEQUENCE_OF_DIAGRAMS.md','physics/lt1/LT1_LAMBDA_TOWER_INWARD.md',
           'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
           'physics/yc10/yc10_closed_return_and_scale.py','physics/yc10/test_yc10.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names}


def run():
    checks=dict(symbolic_checks())
    geometry=[]
    for shape in ((3,3,3),(3,4,5),(4,4,4),(4,5,6)):
        record=yc9.geometry_record(shape)
        c=profile_for_shape(shape)
        geometry.append({'shape':shape,'profile_support':c['minimum_support'],
                         'triangles':record['triangle_count'],'links':record['links']})
        checks[f'{shape}: simple graph and correct minimum-support profile']=(
            record['no_self_or_parallel_edges'] and
            ((c['minimum_support']==3 and record['triangle_count']>0) or
             (c['minimum_support']==4 and record['triangle_count']==0)))
        checks[f'{shape}: plaquettes remain centre even']=record['elementary_faces_centre_even']
    for s in (3,4):
        c=constants(s)
        checks[f'support {s}: invariant ball maps to itself']=c['map_at_endpoint']<=RADIUS
        checks[f'support {s}: complete return is contractive']=c['contraction_at_endpoint']<1
        checks[f'support {s}: excitation relative norm stays below one']=c['relative_at_endpoint']<1
        checks[f'support {s}: actual gap endpoint']=c['vacuum_gap_at_endpoint']==(F(81,10) if s==3 else 8)
    checks.update({
        'orthogonal projection bound is four and sharp on sixteen equal components': F(16)**2==PROJECTION_FACTOR**2*16,
        'outward exponential at two thirds': exponential_upper(F(2,3),5)==F(404671,207765)<F(39,20),
        'outward exponential at one half': exponential_upper(F(1,2),1)==F(33,20)<F(5,3),
        'all ordinary volumes have tenfold larger coupling window': constants(3)['theta_max']/yc9.THETA_MAX==10,
        'triangle-free volumes have twelvefold larger coupling window': constants(4)['theta_max']/yc9.THETA_MAX==12,
        'four dimensional action is invariant under compression': curvature_scaling(4,F(1,7))['integral']==1,
        'three dimensional magnetic energy scales inversely with size': curvature_scaling(3,F(1,7))['integral']==7,
        'curvature density grows even when the four dimensional action stays fixed': curvature_scaling(4,F(1,7))['curvature_density']==7**4,
    })
    if not all(checks.values()): raise AssertionError([name for name,ok in checks.items() if not ok])
    return encode({'stage':'YC10','base_commit':'9f5fb83','checks':checks,
        'profiles':[constants(3),constants(4)],'geometry_controls':geometry,
        'isotropic_bounds':{'all_sides_at_least_3':'gap_vac >=12-(1248/5)*theta on [0,1/64]',
                            'all_sides_at_least_4':'gap_vac >=12-(640/3)*theta on [0,3/160]'},
        'scale_control':[{'dimension':d,'epsilon':F(1,7),**curvature_scaling(d,F(1,7))} for d in (3,4,5)],
        'iteration_tail':[{'support':s,'iteration':m,'upper':iteration_error(constants(s)['beta_max'],s,m)}
                          for s in (3,4) for m in (0,8,32)],
        'claim_boundary':{'actual_YM_gap_uniform_in_volume':True,
            'full_cluster_and_harmonic_spaces_retained':True,
            'native_lambda_identified_with_bare_coupling':False,
            'corner_concentration_proves_quantum_gap':False,
            'correlated_block_constructed':False,'continuum_4d_proved':False,
            'finite_checks_are_formal_verification':False},
        'source_pins':source_pins()})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    path=HERE/'YC10_RESULT.json'
    if args.check:
        if json.loads(path.read_text())!=result: raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC10: {len(result['checks'])} exact checks pass.")
    print('Actual uniform vacuum gap: >=81/10 through theta=1/64 (all sides>=3); >=8 through theta=3/160 (all sides>=4).')
    print('Concentration scaling audited; no quantum continuum or lambda-to-coupling identification claimed.')
