"""YC13: lambda-frame transport, observed covariance and the full cube return.

Exact identities and outward surface-integral bounds accompany the written
proof. Native lambda is not identified with physical time or YM theta.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
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


u8 = load('yc13_up8', 'uncut/up8/up8_deformed_compass.py')
u9 = load('yc13_up9', 'uncut/up9/up9_haar_source_flow.py')
y12 = load('yc13_yc12', 'physics/yc12/yc12_correlated_cube_gap.py')
I = u9.Interval
R = u8.R
OBS = sp.Matrix([[0, -1], [-1, 0]])
x, y = sp.symbols('x y', nonnegative=True)
S, V, lam, ratio = sp.symbols('S V lambda r', positive=True)


def zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.cancel(a) == 0 for a in value)
    return sp.cancel(value) == 0


@lru_cache(None)
def family():
    # U=S^3/V + lambda(S^2+V^2)=z0*(1+x+y),
    # x=lambda V/S, y=lambda(V/S)^3, z0=S^3/V.
    m = sp.Rational(3, 2)/(3+x)
    ds = lambda f: 3*f-x*sp.diff(f, x)-3*y*sp.diff(f, y)
    dv = lambda f: (-1+3*m)*f+(1-m)*(x*sp.diff(f, x)+3*y*sp.diff(f, y))
    U = 1+x+y
    L = sp.Matrix([[ds(ds(U)), dv(ds(U))], [ds(dv(U)), dv(dv(U))]]).applyfunc(sp.factor)
    _, p, q, w = u8.parts(L)
    delta = sp.factor(p*p+q*q-w*w)
    ell = sp.cancel(w*w/delta)
    phi = [sp.cancel((p*sp.diff(q, t)-q*sp.diff(p, t))/(p*p+q*q)) for t in (x, y)]
    fxy = sp.factor(sp.diff(ell, x)*phi[1]-sp.diff(ell, y)*phi[0])
    return dict(m=m, L=L, p=p, q=q, w=w, delta=delta, ell=ell,
                phix=phi[0], phiy=phi[1], fxy=fxy)


def positive_power(value, n):
    if not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer power required')
    result = I.exact(1)
    for _ in range(n):
        result = result*value
    return result


def curvature_interval(lambdas, ratios):
    """Outward F_lambda,r/R for lambda>=0, r>0; exact rational operations."""
    if lambdas.lo < 0 or ratios.lo <= 0:
        raise ValueError('nonnegative lambda and positive ratio required')
    xx = lambdas*ratios
    yy = lambdas*positive_power(ratios, 3)
    a = xx+3
    p0 = ((((16*xx+188)*xx+756)*xx+1350)*xx+945)/(8*positive_power(a, 3))
    pp = p0-2*yy
    pp2 = pp*pp
    if pp.lo <= 0 <= pp.hi:
        pp2 = I(Q(0), max(pp.lo**2, pp.hi**2))
    # Manifestly positive discriminant avoids an artificial interval zero.
    dd = pp2+9*positive_power(2*xx+3, 3)/(4*positive_power(a, 3))
    fxy = 81*xx*positive_power(2*xx+3, 3)/(8*positive_power(a, 6)*dd*dd)
    return 2*lambdas*positive_power(ratios, 3)*fxy


@lru_cache(None)
def loop_certificate(partitions=16):
    """Whole rectangle [0,1/4] x [1/2,1], including between all sample nodes."""
    if type(partitions) is not int or partitions <= 0:
        raise ValueError('positive integer partition count required')
    dl, dr = Q(1, 4*partitions), Q(1, 2*partitions)
    integral = I.exact(0)
    for i in range(partitions):
        for j in range(partitions):
            cell = curvature_interval(I(i*dl, (i+1)*dl), I(Q(1, 2)+j*dr, Q(1, 2)+(j+1)*dr))
            integral = integral+dl*dr*cell
    return -integral  # Counterclockwise loop: Theta=-integral F/R.


def uniform_mode_return(zeta):
    zeta = Q(zeta)
    if zeta >= 18:
        raise ValueError('zeta must be below 18')
    return Q(1, 2)/(18-zeta)+Q(1, 2)/(24-zeta)+Q(3, 2)/(26-zeta)+Q(1, 4)/(32-zeta)


def return_step(theta1, theta2, zeta):
    theta1, theta2, zeta = Q(theta1), Q(theta2), Q(zeta)
    if not 0 <= theta1 <= theta2 or zeta >= 18:
        raise ValueError('ordered nonnegative couplings and zeta<18 required')
    u1, u2 = 6*theta1/(18-zeta), 6*theta2/(18-zeta)
    if u2 >= 1:
        raise ValueError('full hidden reserve is exhausted')
    return (1-u1*u1)/(1-u2*u2)


def gap_window(theta_max, zeta):
    theta_max, zeta = Q(theta_max), Q(zeta)
    if theta_max < 0 or not 0 < zeta < 12:
        raise ValueError('positive candidate gap below 12 required')
    factor = return_step(0, theta_max, zeta)
    cost = theta_max**2*factor*uniform_mode_return(zeta)
    reserve = 12-zeta-cost
    if reserve <= 0:
        raise ValueError('face block is not certified positive')
    return dict(theta_max=theta_max, centered_threshold=zeta,
                full_hidden_ratio=6*theta_max/(18-zeta),
                return_upper=cost, face_reserve=reserve, gauge_gap_lower=zeta)


def verified_certificate(relative):
    old = json.loads((ROOT/relative).read_text())
    for name, digest in old['source_pins'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise AssertionError(f'frozen source pin changed: {name}')
    return old


def cube_mixed_certificates():
    # Reuse only a verified frozen certificate, not rounded displayed decimals.
    old = verified_certificate('physics/yc11/YC11_RESULT.json')
    result = []
    for record in old['finite_certificates']:
        l, ss, vv = Q(record['lambda_value']), Q(record['S']), Q(record['V'])
        fs = I(*(Q(a) for a in record['full_integral']['curvature']))
        result.append(dict(lambda_value=l, S=ss, V=vv, state=fs.record(),
                           mixed_lambda_S=(-vv/l*fs).record(), mixed_lambda_V=(ss/l*fs).record()))
    return result


def symbolic_checks():
    f = family()
    checks = {f'UP8 connection: {k}': bool(v) for k, v in u8.connection_checks().items()}
    U = S**3/V+lam*(S*S+V*V)
    T = sp.diff(U, S)
    m = sp.factor(-V*sp.diff(T, V)/(S*sp.diff(T, S)))
    ds = lambda g: S*sp.diff(g, S)
    dv = lambda g: V*sp.diff(g, V)+m*ds(g)
    direct = sp.Matrix([[ds(ds(U)), dv(ds(U))], [ds(dv(U)), dv(dv(U))]])
    subs = {x: lam*V/S, y: lam*(V/S)**3}
    checks['response comes from the declared scalar potential'] = zero(direct-S**3/V*f['L'].subs(subs))
    checks['potential is homogeneous of degree two at every lambda'] = zero(S*sp.diff(U, S)+V*sp.diff(U, V)-2*U)
    checks['ordinary Hessian is strictly positive by base determinant and positive addition'] = (
        sp.hessian(S**3/V, (S, V)).det() == 3*S**4/V**4 and
        sp.hessian(U, (S, V))-sp.hessian(S**3/V, (S, V)) == 2*lam*sp.eye(2))
    checks['zero member has the previously certified closed flat response'] = (
        f['w'].subs({x: 0, y: 0}) == 0 and f['delta'].subs({x: 0, y: 0}) == sp.Rational(1369, 64))
    checks['strict cut positivity holds globally without cancellation'] = zero(
        f['q']**2-f['w']**2-9*(2*x+3)**3/(4*(x+3)**3))
    checks['mixed curvature has a strictly positive factored numerator'] = zero(
        f['fxy']*f['delta']**2-81*x*(2*x+3)**3/(8*(x+3)**6))
    jac = sp.Matrix([[sp.diff(subs[a], b) for b in (S, V)] for a in (x, y)]).det()
    checks['every fixed-lambda state surface has zero curvature'] = zero(jac)
    sx = {x: lam*ratio, y: lam*ratio**3}
    jac_lr = sp.Matrix([[sp.diff(sx[a], b) for b in (lam, ratio)] for a in (x, y)]).det()
    checks['lambda and shape give a genuine mixed area'] = zero(jac_lr-2*lam*ratio**3)
    leading_xy = sp.factor(sp.limit(f['fxy']/x, x, 0).subs(y, 0))
    checks['first mixed curvature coefficient is exact'] = leading_xy == sp.Rational(1536, 1874161)
    actual = sp.factor((2*lam*ratio**3*f['fxy'].subs(sx)).subs({lam: sp.Rational(1, 4), ratio: sp.Rational(1, 2)}))
    checks['finite mixed witness is nonzero'] = actual == sp.Rational(5184000000000000, 959239754252475521293)
    a, b, c = sp.symbols('a b c', real=True, nonzero=True)
    quadratic = (a*S*S+2*b*S*V+c*V*V)/2
    mq = -b*V/(a*S)
    dq = lambda g: V*sp.diff(g, V)+mq*ds(g)
    qL = sp.Matrix([[ds(ds(quadratic)), dq(ds(quadratic))], [ds(dq(quadratic)), dq(dq(quadratic))]])
    _, pp, qq, ww = u8.parts(qL)
    checks['arbitrary quadratic coefficients retain q equals minus w'] = zero(qq+ww)
    checks['quadratic families reduce lost weight to squared angle ratio'] = zero(ww*ww/(pp*pp+qq*qq-ww*ww)-(qq/pp)**2)
    aa, pp, qq, ww = sp.symbols('aa pp qq ww', real=True)
    LL = aa*sp.eye(2)+pp*u8.K+qq*u8.J+ww*R
    changed = OBS*LL*OBS.T
    checks['observation is the signed involution with transported orientation'] = OBS*OBS == sp.eye(2) and OBS*R*OBS.T == -R
    checks['observed response keeps magnitude and reverses signed turn'] = u8.parts(changed) == (aa, -pp, qq, -ww)
    checks['full compass return transforms with the observed frame'] = zero(u8.cycle(changed)-OBS*u8.cycle(LL)*OBS.T)
    # Pullback from k=lambda S, h=lambda V: arbitrary first jets suffice.
    elx, ely, phx, phy = sp.symbols('ell_k ell_h phi_k phi_h')
    ell_l, phi_l = S*elx+V*ely, S*phx+V*phy
    state = lam**2*(elx*phy-ely*phx)
    mixed_s, mixed_v = ell_l*lam*phx-lam*elx*phi_l, ell_l*lam*phy-lam*ely*phi_l
    checks['source amplitude mixed S curvature is fixed by state curvature'] = zero(mixed_s+V*state/lam)
    checks['source amplitude mixed V curvature is fixed by state curvature'] = zero(mixed_v-S*state/lam)
    checks['amplitude family has an exact null scale direction'] = zero(lam*ell_l-S*lam*elx-V*lam*ely) and zero(lam*phi_l-S*lam*phx-V*lam*phy)
    v, t1, t2 = sp.symbols('v t1 t2', nonnegative=True)
    step_difference = 1/(1-t2*t2*v)-1/(1-t1*t1*v)
    checks['full-hidden even return step has positive scalar spectral difference'] = zero(
        step_difference-(t2*t2-t1*t1)*v/((1-t1*t1*v)*(1-t2*t2*v)))
    return checks


def encode(obj):
    if isinstance(obj, Q):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [encode(v) for v in obj]
    return obj


def source_pins():
    paths = ['physics/SIGNED_COMPASS_YM_HANDOFF.md',
             'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
             'uncut/up8/up8_deformed_compass.py',
             'uncut/up9/UP9_HAAR_SOURCE_AND_FLAT_LIMIT.md',
             'uncut/up9/up9_haar_source_flow.py',
             'physics/yc11/YC11_CLOSED_CUBE_RESPONSE.md', 'physics/yc11/YC11_RESULT.json',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md',
             'physics/yc12/yc12_correlated_cube_gap.py', 'physics/yc12/YC12_RESULT.json',
             'physics/yc13/YC13_FRAME_FLOW_AND_RETURN.md',
             'physics/yc13/yc13_frame_flow.py', 'physics/yc13/test_yc13.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    verified_certificate('physics/yc12/YC12_RESULT.json')
    checks = symbolic_checks()
    loop = loop_certificate()
    checks['entire frame-state loop has a strict returned angle'] = loop.lo < loop.hi < 0
    checks['certified loop cannot hide a full turn'] = -Q(1, 100) < loop.lo
    windows = [gap_window(1, 11), gap_window(Q(3, 2), Q(42, 5)), gap_window(2, Q(27, 5))]
    for row in windows:
        t, z = row['theta_max'], row['centered_threshold']
        pencil = y12.pencil(t, z, True)
        face = [a[1:] for a in pencil[1:]]
        checks[f'theta<={t}: exact face inertia gives a full-window gap {z}'] = (
            y12.y4.inertia(face) == (0, 0, 6) and y12.y4.inertia(pencil) == (1, 0, 6))
    checks['doubling step at centered energy zero retains its full inverse'] = return_step(Q(1, 2), 1, 0) == Q(35, 32)
    checks['resolved cube source norm is eleven quarters'] = all(
        sum(sum(m[i][j] for m in y12.source_grams().values()) for j in range(1, 7)) == Q(11, 4)
        for i in range(1, 7))
    mixed = cube_mixed_certificates()
    for row in mixed:
        checks[f"actual cube lambda={row['lambda_value']}: mixed source curvatures have certified signs"] = (
            Q(row['mixed_lambda_S'][1]) < 0 < Q(row['mixed_lambda_V'][0]))
    checks = {k: bool(v) for k, v in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    f = family()
    return encode(dict(stage='YC13', base_commit='89769b8', checks=checks,
                       owner_convention='TVSP before observation; observed VTSP film at lambda=0',
                       family='U_lambda=S^3/V+lambda*(S^2+V^2)',
                       normalized_geometry={k: str(f[k]) for k in ('m', 'p', 'q', 'w')},
                       factored_fxy_times_delta_squared='81*x*(2*x+3)^3/(8*(x+3)^6)',
                       first_mixed_curvature='F_lambda,r/R=(3072/1874161)*lambda^2*r^4+O(lambda^3)',
                       loop=dict(lambda_interval=['0', '1/4'], ratio_interval=['1/2', '1'],
                                 subdivisions_per_axis=16, outward_binary_bits=160,
                                 returned_angle=loop.record()),
                       actual_cube_mixed_certificates=mixed, full_window_gauge_gaps=windows,
                       claim_boundary=dict(frame_state_mixed_curvature=True,
                           signed_observation_is_a_covariant_chart=True,
                           observation_swap_alone_generates_curvature=False,
                           physical_measurement_or_time_law=False,
                           lambda_is_theta_or_spatial_RG=False,
                           new_full_window_actual_cube_gauge_gap=True,
                           improved_charged_block_or_uniform_lattice_window=False,
                           energy_or_continuum_mass_derived=False, formal_machine_verification=False),
                       source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC13_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC13: {len(result['checks'])} exact/outward checks pass.")
    print('Native frame transport and complete cube return; gauge gap >=27/5 through theta=2.')
