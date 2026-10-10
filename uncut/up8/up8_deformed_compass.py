"""UP8: centre-sourced cut deformation and native seam transport.

Exact algebra and complete polynomial sign certificates; no physical mass/gap.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.diag(1, -1)
J = R*K
I = sp.eye(2)
x, y = sp.symbols('x y', nonnegative=True)


def zero(expr):
    if isinstance(expr, sp.MatrixBase):
        return all(sp.cancel(v) == 0 for v in expr)
    return sp.cancel(expr) == 0


def parts(L):
    return tuple(sp.factor(e) for e in (
        (L[0, 0]+L[1, 1])/2, (L[0, 0]-L[1, 1])/2,
        (L[0, 1]+L[1, 0])/2, (L[1, 0]-L[0, 1])/2))


def cycle(L):
    _, p, q, w = parts(L)
    delta = sp.factor(p*p+q*q-w*w)
    X = p*K+q*J+w*R
    return (R*X*R*X/delta).applyfunc(sp.factor)


def polar_tangent(M):
    return sp.factor((M[1, 0]-M[0, 1])/(M[0, 0]+M[1, 1]))


@lru_cache(None)
def family():
    h0 = 1+2*x+y
    m = (3+4*x+6*y)/(6*h0)
    # These act on f after factoring out z0=S^3/V.
    ds = lambda f: 3*f+x*sp.diff(f, x)
    dv = lambda f: (-1+3*m)*f+m*x*sp.diff(f, x)-y*sp.diff(f, y)
    u = 1+x+y
    L = sp.Matrix([[ds(ds(u)), dv(ds(u))], [ds(dv(u)), dv(dv(u))]]).applyfunc(sp.factor)
    _, p, q, w = parts(L)
    rho2 = sp.factor(p*p+q*q)
    delta = sp.factor(rho2-w*w)
    ell = sp.factor(w*w/delta)
    phix = sp.factor((p*sp.diff(q, x)-q*sp.diff(p, x))/rho2)
    phiy = sp.factor((p*sp.diff(q, y)-q*sp.diff(p, y))/rho2)
    fxy = sp.factor(sp.diff(ell, x)*phiy-sp.diff(ell, y)*phix)
    H = sp.cancel(5184*h0**6*delta)
    G = sp.cancel(fxy*H**2/(10368*x*(4*y+1)*h0**3*(4*x+3*y+3)**2))
    return dict(m=m, L=L, p=p, q=q, w=w, delta=delta, ell=ell,
                phix=phix, phiy=phiy, fxy=fxy, H=H, G=G)


def at(expr, xx, yy):
    out = expr.subs({x: sp.Rational(xx), y: sp.Rational(yy)})
    return out.applyfunc(sp.factor) if isinstance(out, sp.MatrixBase) else sp.factor(out)


def polynomial_record(expr):
    poly = sp.Poly(expr, x, y)
    return {'degree': int(poly.total_degree()), 'terms': len(poly.terms()),
            'constant': str(poly.eval({x: 0, y: 0})),
            'coefficients': [{'powers': list(powers), 'coefficient': str(c)}
                             for powers, c in poly.terms()]}


@lru_cache(None)
def connection_checks():
    checks = {}
    t = sp.symbols('t', positive=True)  # t=e^eta, not a physical clock
    sh, ch = (t-1/t)/2, (t+1/t)/2
    # Evaluate at phi=0. All other orientations follow by common SO(2) conjugation.
    kap = ch*K+sh*R
    B = (R*kap)**2
    h = ch*I-sh*J
    dphi = lambda M: (R*M-M*R)/2
    checks['cut squares to identity'] = zero(kap*kap-I)
    checks['positive root gives the complete compass return'] = zero(h*h-B)
    checks['return has determinant one'] = zero(B.det()-1)
    checks['flat cut closes exactly'] = B.subs(t, 1) == I
    plus, minus = sp.Matrix([t+1,t-1]), sp.Matrix([t-1,t+1])
    checks['oblique cut retains both signed eigenlines'] = zero(kap*plus-plus) and zero(kap*minus+minus)
    Xi = B.inv()*sp.diff(B, t)
    Xj = B.inv()*dphi(B)
    Ai, Aj = Xi/2, Xj/2
    ai = h*Ai*h.inv()-sp.diff(h, t)*h.inv()
    aj = h*Aj*h.inv()-dphi(h)*h.inv()
    checks['orthonormal connection has no radial component'] = zero(ai)
    checks['seam connection equals normalized defect squared'] = zero(aj-sh**2*R)
    checks['native transport preserves its declared metric'] = (
        zero(sp.diff(B, t)-Ai.T*B-B*Ai) and zero(dphi(B)-Aj.T*B-B*Aj))
    checks['full Maurer Cartan reading is flat'] = zero(
        sp.diff(Xj, t)-dphi(Xi)+Xi*Xj-Xj*Xi)
    Fij = sp.diff(Aj, t)-dphi(Ai)+Ai*Aj-Aj*Ai
    checks['half share curvature retains commutator factor'] = zero(Fij+(Xi*Xj-Xj*Xi)/4)
    checks['curvature in orthonormal frame is d lost wedge d angle'] = zero(
        h*Fij*h.inv()-sp.diff(sh**2, t)*R)
    a, p, q, w, tau = sp.symbols('a p q w tau', real=True)
    L = a*I+p*K+q*J+w*R
    checks['polynomial tower rescales one fixed traceless direction'] = zero(
        L+tau*L*L-(a+tau*(a*a+p*p+q*q-w*w))*I
        -(1+2*tau*a)*(p*K+q*J+w*R))
    a1,b1,c1,a2,b2,c2 = sp.symbols('a1 b1 c1 a2 b2 c2', real=True)
    M = (a2*I+b2*K+c2*J)*(a1*I+b1*K+c1*J)
    mm, _, _, ww = parts(M)
    checks['noncollinear cycles generate the skew coefficient'] = zero(ww-(c2*b1-b2*c1))
    unnormalized_positive = (mm*I-ww*R)*M
    checks['polar rotation removes the whole skew part'] = zero(unnormalized_positive-unnormalized_positive.T)
    return checks


@lru_cache(None)
def source_checks():
    checks = {}
    f = family()
    S, V, lam = sp.symbols('S V lambda', positive=True)
    U = S**3/V+lam*S**4/V+lam*S**3/V**2
    T = sp.diff(U, S)
    m = sp.factor(-V*sp.diff(T, V)/(S*sp.diff(T, S)))
    ds = lambda z: S*sp.diff(z, S)
    dv = lambda z: V*sp.diff(z, V)+m*ds(z)
    direct = sp.Matrix([[ds(ds(U)), dv(ds(U))], [ds(dv(U)), dv(dv(U))]])
    sub = {x: lam*S, y: lam/V}
    checks['response is obtained from the centre derivatives'] = zero(
        direct-S**3/V*f['L'].subs(sub))
    checks['frame defect source is retained'] = zero(
        direct[1, 0]-direct[0, 1]-ds(m)*S*T)
    checks['exact deformation defect'] = zero(
        x*sp.diff(f['m'], x)+x*(1+4*y)/(3*(1+2*x+y)**2))
    for aa, bb in ((3, -1), (4, -1), (3, -2)):
        u = S**aa*V**bb
        hes = sp.hessian(u, (S, V))
        checks[f'convex monomial determinant {aa},{bb}'] = zero(
            hes.det()-aa*bb*(1-aa-bb)*S**(2*aa-2)*V**(2*bb-2)) and aa*(aa-1)>0 and aa*bb*(1-aa-bb)>=0
    origin = {x: 0, y: 0}
    checks['flat source cross response'] = f['L'].subs(origin) == sp.Matrix([[9, sp.Rational(3,2)], [sp.Rational(3,2), sp.Rational(1,4)]])
    checks['flat source discriminant'] = f['delta'].subs(origin) == sp.Rational(1369, 64)
    checks['flat source return and connection weight'] = cycle(f['L'].subs(origin)) == I and f['ell'].subs(origin) == 0
    for name in ('H', 'G'):
        poly = sp.Poly(f[name], x, y)
        checks[f'complete positive coefficient certificate {name}'] = (
            all(c>0 for c in poly.coeffs()) and poly.eval(origin)>0)
    checks['cut discriminant factorization'] = zero(f['delta']-f['H']/(5184*(1+2*x+y)**6))
    checks['native curvature factorization'] = zero(f['fxy']-
        10368*x*(4*y+1)*(1+2*x+y)**3*(4*x+3*y+3)**2*f['G']/f['H']**2)
    checks['seam angle is smooth in the positive quadrant'] = all(
        all(c>0 for c in sp.Poly(sp.fraction(sp.factor(f[k]))[0], x, y).coeffs()) for k in ('p', 'q'))
    checks['first defect jet'] = sp.diff(f['w'], x).subs(origin) == -sp.Rational(1,2)
    checks['second lost weight jet'] = sp.diff(f['ell'], x, 2).subs(origin)/2 == sp.Rational(16,1369)
    checks['two independent first seam jets'] = (
        f['phix'].subs(origin) == -sp.Rational(652,1369) and f['phiy'].subs(origin) == sp.Rational(420,1369))
    checks['first nonzero native curvature is cubic in deformation'] = (
        sp.cancel(f['fxy']/x).subs(origin) == sp.Rational(13440,1874161))
    checks['nonzero exact curvature control one'] = at(f['fxy'], 1, 1) == sp.Rational(76489856000,2084237405595601)
    checks['nonzero exact curvature control two'] = at(f['fxy'], sp.Rational(1,2), sp.Rational(1,3)) == sp.Rational(483327,4355431112)
    b1, b2 = cycle(at(f['L'],1,1)), cycle(at(f['L'],2,2))
    angle = polar_tangent(b2*b1)
    checks['two sourced cycles have a nonzero polar turn'] = angle == -sp.Rational(273483623347200,8726504489924666083)
    checks['reversing ordered cycles reverses their turn'] = polar_tangent(b1*b2) == -angle
    checks['large deformation does not force monotone curvature growth'] = (
        sp.limit(lam**2*f['fxy'].subs({x:lam,y:lam}),lam,sp.oo) == 0)
    return checks


def source_pins():
    paths = [
        'uncut/up1/UP1_DEGREE_OF_THE_POTENTIAL.md',
        'uncut/up4/UP4_COST_OF_COMMUTING_CUTS.md',
        'uncut/up5/UP5_SEED_HIERARCHY_CHECKED.md',
        'uncut/up6/UP6_EIGHT_SCALE_OPERATIONS.md',
        'uncut/up7/UP7_TURN_AND_CUT.md', 'uncut/up7/up7_turn_and_cut.py',
        'uncut/tt1/TT1_TIME_IN_THE_THERMO_DIAGRAM.md',
        'uncut/dw1/DW1_DEFECT_IS_THE_TURN_PART.md',
        'uncut/cp1/CP1_SELF_SIMILAR_CENTRE.md',
        'uncut/on1/ON1_ORIGIN_OF_NONCOMMUTATION.md',
        'uncut/ss1/SS1_SCALE_IS_THE_SOURCE.md',
        'physics/ph3/PH3_TWO_READINGS_TRANSPORT.md',
        'physics/nc1/NC1_NATIVE_CURVATURE_OF_THE_FIELD.md',
        'physics/lc1/LC1_LOST_IS_COUNTED.md',
        'physics/qc4/QC4_RETURN_IS_MEMORY.md',
        'physics/qc5/QC5_MOVING_SHARE.md',
        'physics/hc1/HC1_SIGNED_HORIZON_AND_NATIVE_CURVATURE.md',
        'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
        'uncut/up8/up8_deformed_compass.py', 'uncut/up8/test_up8.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    checks = {name:bool(ok) for name,ok in {**connection_checks(), **source_checks()}.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    f = family()
    return {'stage':'UP8', 'base_commit':'d1ac633', 'checks':checks,
        'source':'U_lambda = S^3/V + lambda*S^4/V + lambda*S^3/V^2',
        'seam_connection':'a = ell R dphi, ell = w^2/(p^2+q^2-w^2)',
        'curvature':'F = R d ell wedge d phi',
        'positive_polynomials':{k:polynomial_record(f[k]) for k in ('H','G')},
        'point_controls':[{'x':str(xx),'y':str(yy),
            'response':[[str(z) for z in row] for row in at(f['L'],xx,yy).tolist()],
            'ell':str(at(f['ell'],xx,yy)), 'curvature_xy':str(at(f['fxy'],xx,yy))}
            for xx,yy in ((0,0),(1,1),(sp.Rational(1,2),sp.Rational(1,3)))],
        'claim_boundary':{'explicit_flat_zero_deformation':True,
            'strict_native_curvature_for_positive_deformation_in_this_family':True,
            'transport_protocol_declared':True,
            'curvature_monotone_in_deformation':False,
            'deformation_is_the_polynomial_tower_parameter':False,
            'isotropic_contraction_or_arrow_proved':False,
            'all_prior_physics_identified_as_lambda_zero':False,
            'full_SU2_YM_connection_or_new_mass_gap':False,
            'finite_checks_are_formal_verification':False},
        'source_pins':source_pins()}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    path=HERE/'UP8_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result: raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(f"UP8: {len(result['checks'])} exact checks pass.")
    print('Flat lambda=0; native seam curvature nonzero for every positive lambda in the declared family.')
