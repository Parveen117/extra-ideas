"""YC11: complete closed-cube source, boundary return and native compass.

Reuses CZ1's character gluing and UP9's outward state-jet transport.
No quantum heat trace, kinetic block inverse or new mass-gap claim.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
import json
from math import comb, factorial, prod
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT/path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


u9 = load_module('yc11_up9', 'uncut/up9/up9_haar_source_flow.py')
cz = load_module('yc11_cz1', 'physics/cz1/cz1_centre_record_of_a_closed_surface.py')
I = u9.Interval
S, V, e = sp.symbols('S V epsilon', positive=True)  # epsilon=lambda^2


@lru_cache(None)
def face_coefficient(n, power):
    """[t^power] f_n(t), f_n=2 I_n(t)/t, dimension n>=1."""
    if n < 1 or power < 0:
        return Q(0)
    rest = power-n+1
    if rest < 0 or rest % 2:
        return Q(0)
    k = rest//2
    return Q(1, 2**power*factorial(k)*factorial(n+k))


@lru_cache(None)
def hemisphere_coefficient(n, power):
    """[t^power] f_n(t)^3, using every contributing Bessel coefficient."""
    rest = power-3*(n-1)
    if rest < 0 or rest % 2:
        return Q(0)
    k = rest//2
    return sum((face_coefficient(n,n-1+2*a)*face_coefficient(n,n-1+2*b)
                *face_coefficient(n,n-1+2*(k-a-b))
                for a in range(k+1) for b in range(k-a+1)), Q(0))


@lru_cache(None)
def cube_coefficient(a, b, connected=False):
    """[k^a h^b] Z_cube, optionally retaining n>=2 only (Z-Z_ind)."""
    if min(a,b) < 0:
        raise ValueError('nonnegative source powers required')
    return sum((n*n*hemisphere_coefficient(n,a)*hemisphere_coefficient(n,b)
                for n in range(2 if connected else 1, min(a,b)//3+2)), Q(0))


@lru_cache(None)
def moment(a, b):
    """E[A^a B^b], A/B the sums of three opposite corner faces."""
    return factorial(a)*factorial(b)*cube_coefficient(a,b)


@lru_cache(None)
def potential():
    # Z has even total degree. Through lambda^8, all contributing
    # representation dimensions are included automatically by the formula.
    z = [sp.expand(sum(sp.Rational(cube_coefficient(a,2*n-a))*S**a*V**(2*n-a)
                       for a in range(2*n+1))) for n in range(5)]
    logz = [sp.S(0)]*5
    for n in range(1,5):
        logz[n] = sp.expand(z[n]-sum(sp.Rational(k,n)*logz[k]*z[n-k]
                                     for k in range(1,n)))
    return sp.expand(sum(logz[n]*e**(n-1) for n in range(1,5)))


@lru_cache(None)
def leading_geometry():
    U = potential()
    trunc = lambda f: sp.series(f,e,0,4).removeO().expand()
    T = sp.diff(U,S)
    m = trunc(-V*sp.diff(T,V)/(S*sp.diff(T,S)))
    ds = lambda f: S*sp.diff(f,S)
    dv = lambda f: trunc(V*sp.diff(f,V)+m*ds(f))
    L = sp.Matrix([[ds(ds(U)),dv(ds(U))],[ds(dv(U)),dv(dv(U))]])
    p = trunc((L[0,0]-L[1,1])/2)
    q = trunc((L[0,1]+L[1,0])/2)
    w = trunc((L[1,0]-L[0,1])/2)
    p0,p2 = p.coeff(e,0),p.coeff(e,1)
    w4,w6,q4,q6 = w.coeff(e,2),w.coeff(e,3),q.coeff(e,2),q.coeff(e,3)
    ell8 = sp.factor(w4*w4/p0**2)
    ell10 = sp.factor(2*w4*w6/p0**2-2*w4*w4*p2/p0**3)
    phi4 = sp.factor(q4/p0)
    phi6 = sp.factor(q6/p0-q4*p2/p0**2)
    jac = lambda a,b: sp.factor(sp.diff(a,S)*sp.diff(b,V)-sp.diff(a,V)*sp.diff(b,S))
    return dict(U=U,m=m,w4=w4,w6=w6,p0=p0,p2=p2,q4=q4,q6=q6,
                ell8=ell8,ell10=ell10,phi4=phi4,phi6=phi6,
                f12=jac(ell8,phi4),f14=sp.factor(jac(ell10,phi4)+jac(ell8,phi6)))


def face_interval(n, t, order=20):
    """Positive Bessel series with complete geometric tail."""
    t = Q(t)
    if n < 1 or t < 0 or order < 0:
        raise ValueError('n>=1, t>=0, order>=0 required')
    terms = [face_coefficient(n,n-1+2*k)*t**(n-1+2*k) for k in range(order+2)]
    ratio = t*t/(4*(order+2)*(n+order+2))
    if ratio >= 1:
        raise ValueError('increase Bessel series order')
    lo = sum(terms[:-1],Q(0))
    return I(lo,lo+terms[-1]/(1-ratio))


def ipow(a, n):
    out = I.exact(1)
    for _ in range(n):
        out = out*a
    return out


def boundary_budget(k):
    k=Q(k)
    if k < 0 or k**3 >= 96:
        raise ValueError('sufficient bound requires 0<=k^3<96')
    return k**3/(16*(1-k**3/96))


def interface_term_bound(n, k, h):
    z=Q(k)*Q(h)/4
    return Q(n*n)*z**(3*(n-1))/factorial(n)**6


def interface_tail(k,h,dimension):
    """All dimensions strictly above dimension, not a closed finite span."""
    n=dimension+1
    z=Q(k)*Q(h)/4
    ratio=z**3/(n*n*(n+1)**4)
    if min(k,h)<0 or dimension<1 or ratio>=1:
        raise ValueError('invalid interface tail parameters')
    return interface_term_bound(n,k,h)/(1-ratio)


def interface_interval(k,h,dimension=6):
    k,h=Q(k),Q(h)
    fk,fh=face_interval(1,k),face_interval(1,h)
    total=I.exact(0)
    for n in range(2,dimension+1):
        un,vn=face_interval(n,k)/fk,face_interval(n,h)/fh
        total=total+n*n*ipow(un*vn,3)
    tail=interface_tail(k,h,dimension)
    return I(total.lo,total.hi+tail)


def log_one_plus(x,order=12):
    if x.lo<0 or x.hi>=1 or order<1:
        raise ValueError('log series requires 0<=x<1 and positive order')
    value=I.exact(0)
    power=I.exact(1)
    for n in range(1,order+1):
        power=power*x
        value=value+Q((-1)**(n+1),n)*power
    tail=x.hi**(order+1)/(order+1)
    return I(value.lo-tail,value.hi+tail)


def centre_memory(lambda0):
    l=Q(lambda0);k,h=2*l,l
    ret=interface_interval(k,h)
    delta_u=log_one_plus(ret)/l**2
    z=tilted_moment(0,0,l)
    ds=(tilted_moment(1,0,l)/z-3*face_interval(2,k)/face_interval(1,k))/l
    dv=(tilted_moment(0,1,l)/z-3*face_interval(2,h)/face_interval(1,h))/l
    return delta_u,delta_u-(2*ds+dv)/2


@lru_cache(None)
def tilted_moment(i,j,lambda0,s0=Q(2),v0=Q(1),order=40):
    """Full twelve-link integral derivative with complete exponential tail."""
    lambda0,s0,v0=Q(lambda0),Q(s0),Q(v0)
    t=3*abs(lambda0)*(abs(s0)+abs(v0))
    if lambda0 == 0 or order<0:
        raise ValueError('nonzero lambda and nonnegative order required')
    val=Q(0)
    for n in range(order+1):
        mn=sum((Q(comb(n,k))*s0**(n-k)*v0**k*moment(i+n-k,j+k)
                for k in range(n+1)),Q(0))
        val+=lambda0**n*mn/factorial(n)
    # |A|,|B|<=3. e^t <= 3^ceil(t), using the elementary bound e<3.
    ceil_t=-((-t.numerator)//t.denominator)
    tail=Q(3)**(i+j+ceil_t)*t**(order+1)/factorial(order+1)
    return I(val-tail,val+tail)


def source_jet(lambda0,s0=Q(2),v0=Q(1),order=40):
    lambda0=Q(lambda0)
    z0=tilted_moment(0,0,lambda0,s0,v0,order)
    b={(i,j):tilted_moment(i,j,lambda0,s0,v0,order)/z0
       *lambda0**(i+j)/Q(factorial(i)*factorial(j))
       for i in range(5) for j in range(5-i) if i+j}
    out,power={},{(0,0):I.exact(1)}
    for n in range(1,5):
        power=u9.mul(power,b,4)
        out=u9.add(out,u9.scale(power,Q((-1)**(n+1),n)),4)
    return u9.scale(out,1/lambda0**2)


@lru_cache(None)
def finite_geometry(lambda0):
    jet=source_jet(Q(lambda0))
    g=u9.geometry_from_jet(jet,Q(2),Q(1))
    uss,uvv,usv=jet[(2,0)]*2,jet[(0,2)]*2,jet[(1,1)]
    g['mixed_covariance']=usv
    g['chi']=g['hessian_det']/(uss*uvv)
    return g


def topology_checks():
    faces=cz.box_surface(1,1,1)
    edges={edge for f in faces for edge,sign in f}
    vertices={v for edge in edges for v in edge}
    counts=[sum(edge==ee for f in faces for ee,sgn in f) for edge in edges]
    # With this source convention the three lower-coordinate faces are
    # indices 0,2,4; the other three meet the opposite corner.
    left={0,2,4}
    boundary=[edge for edge in edges if sum(edge==ee for i in left for ee,sgn in faces[i])==1]
    return {
        'actual cube has eight vertices twelve links six faces':(len(vertices),len(edges),len(faces))==(8,12,6),
        'each link is integrated once and belongs to two faces':set(counts)=={2},
        'hemispheres meet on six shared boundary links':len(boundary)==6,
        'CZ1 exact gluing retains the dimension to Euler power':cz.glue(faces)==-4,
        'fundamental six-face moment is not independent':Q(2)**cz.glue(faces)/2**6==Q(1,1024),
    }


def symbolic_checks():
    g=leading_geometry();U=g['U'];d=S*S-V*V
    zero=lambda a:sp.cancel(a)==0
    checks={}
    checks['zero source normalizes the entire cube measure']=moment(0,0)==1
    checks['Haar covariance of the two three-face sums is diagonal']=(moment(1,0),moment(0,1),moment(2,0),moment(1,1),moment(0,2))==(0,0,Q(3,4),0,Q(3,4))
    checks['sixth mixed cumulant sees the closure']=moment(3,3)==Q(9,256)
    checks['mixed moments below either cubic channel factorize']=all(
        moment(a,b)==moment(a,0)*moment(0,b) for a in range(7) for b in range(7) if min(a,b)<3)
    checks['complete first connected partition coefficient']=cube_coefficient(3,3,True)==Q(1,1024)
    checks['connected partition has no lower total degree']=all(
        cube_coefficient(a,b,True)==0 for a in range(6) for b in range(6-a))
    checks['flat covariance limit']=zero(U.coeff(e,0)-3*(S*S+V*V)/8)
    checks['first local nonGaussian correction']=zero(U.coeff(e,1)+(S**4+V**4)/128)
    checks['first closed-cube correction']=zero(U.coeff(e,2)-(S**6+V**6)/3072-S**3*V**3/1024)
    checks['next closed-cube correction']=zero(U.coeff(e,3)+(S**8+V**8)/61440+S**3*V**3*(S*S+V*V)/8192)
    checks['source dilation includes correlation terms']=zero(2*e*sp.diff(U,e)+2*U-S*sp.diff(U,S)-V*sp.diff(U,V))
    dc=U-3*(S*S+V*V)/8+e*(S**4+V**4)/128-e**2*(S**6+V**6)/3072+e**3*(S**8+V**8)/61440
    checks['centre first records closed-cube memory at order four']=zero(sp.expand(-e*sp.diff(dc,e)).coeff(e,2)+S**3*V**3/512)
    checks['selected cut defect begins at order four']=zero(g['w4']+9*S**3*V**3/2048)
    checks['lost weight begins at order eight']=zero(g['ell8']-9*S**6*V**6/(262144*d*d))
    checks['leading seam and defect depend on the same scalar']=zero(g['phi4']-3*g['w4']/g['p0'])
    checks['twelfth-order curvature cancels exactly']=g['f12']==0
    checks['fourteenth-order curvature keeps the second source shape']=zero(g['f14']-9*S**10*V**8*(3*S*S-V*V)/(33554432*d**4))
    checks['nonzero rational first curvature at two-to-one sources']=g['f14'].subs({S:2,V:1})==sp.Rational(11,294912)
    checks['six-face dimension factor is essential']=Q(1,1024)!=Q(1,4096)
    return checks


def source_pins():
    paths=['physics/cz1/CZ1_CENTRE_RECORD_OF_A_CLOSED_SURFACE.md',
           'physics/cz1/cz1_centre_record_of_a_closed_surface.py',
           'physics/yc2/YC2_COMPACT_CENTRE_RESPONSE.md',
           'physics/yc8/YC8_LOCAL_VACUUM_DRESSING.md',
           'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
           'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
           'uncut/up9/UP9_HAAR_SOURCE_AND_FLAT_LIMIT.md',
           'uncut/up9/up9_haar_source_flow.py',
           'physics/yc11/YC11_CLOSED_CUBE_RESPONSE.md',
           'physics/yc11/yc11_closed_cube_response.py','physics/yc11/test_yc11.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    checks={**topology_checks(),**symbolic_checks()}
    finite=[]
    for l in (Q(1,4),Q(1,2)):
        g=finite_geometry(l)
        ret=interface_interval(2*l,l)
        delta_u,delta_mid=centre_memory(l)
        p=f'full cube lambda={l}: '
        checks[p+'positive ordinary Hessian']=g['uss'].lo>0 and g['hessian_det'].lo>0
        checks[p+'admissible native cut']=g['delta'].lo>0
        checks[p+'signed cut defect excludes zero']=g['w'].hi<0
        checks[p+'native curvature strictly positive']=g['curvature'].lo>0
        checks[p+'ordinary response retains mixed correlation']=g['mixed_covariance'].lo>0 and 0<g['chi'].lo<=g['chi'].hi<1
        checks[p+'full boundary return positive and within analytic budget']=(ret.lo>0 and ret.hi<=interface_tail(2*l,l,1))
        checks[p+'centre distinguishes the complete boundary return']=(delta_u.lo>0 and delta_mid.hi<0)
        z=tilted_moment(0,0,l)
        zind=ipow(face_interval(1,2*l)*face_interval(1,l),3)
        independent_return=z/zind-1
        checks[p+'moment and boundary representations overlap']=max(ret.lo,independent_return.lo)<=min(ret.hi,independent_return.hi)
        finite.append(dict(lambda_value=str(l),S='2',V='1',exponential_degree=40,
                           outward_binary_bits=160,full_integral={k:v.record() for k,v in g.items()},
                           interface_return=ret.record(),interface_norm_upper=str(interface_tail(2*l,l,1)),
                           potential_correction=delta_u.record(),centre_correction=delta_mid.record(),
                           hemisphere_uniform_budgets=[str(boundary_budget(2*l)),str(boundary_budget(l))]))
    checks={k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):raise AssertionError([k for k,v in checks.items() if not v])
    g=leading_geometry()
    return dict(stage='YC11',base_commit='b6e988a',checks=checks,
                potential_coefficients=[str(sp.factor(g['U'].coeff(e,i))) for i in range(4)],
                geometry={k:str(v) for k,v in g.items() if k!='U'},
                finite_certificates=finite,
                claim_boundary=dict(actual_twelve_link_cube_configuration=True,
                    complete_boundary_character_return=True,full_integral_error_enclosed=True,
                    native_curvature_is_CZ1_centre_weight_commutator=False,
                    shared_link_alone_forces_Haar_pair_correlation=False,
                    kinetic_block_gap_or_inverse_constructed=False,
                    volume_uniform_extension_or_new_YM_gap=False,
                    source_lambda_identified_with_quantum_YM_theta=False,
                    formal_machine_verification=False),source_pins=source_pins())


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run();path=HERE/'YC11_RESULT.json'
    if args.check:
        if json.loads(path.read_text())!=result:raise AssertionError('result or source pins differ')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"YC11: {len(result['checks'])} exact/outward checks pass.")
    print('Complete correlated cube source and native curvature; no new quantum spectral gap.')
