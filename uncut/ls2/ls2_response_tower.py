"""LS2: exact typed jet lift and the actual UP8 rational response tower.

Analytic and all-order conclusions have written proofs in the companion note.
This program checks their algebra and the supplied family's source adapter.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'uncut/up8'))
import up8_deformed_compass as up8

I, R, K, J = up8.I, up8.R, up8.K, up8.J
x, y = up8.x, up8.y
t = sp.symbols('t', real=True)
A, B = sp.symbols('A B', positive=True)
S, V = sp.symbols('S V', positive=True)
SOURCES = (
    'uncut/up3/UP3_SEQUENCE_OF_DIAGRAMS.md',
    'uncut/up6/UP6_EIGHT_SCALE_OPERATIONS.md',
    'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
    'uncut/up8/up8_deformed_compass.py',
    'uncut/ls1/LS1_COMPASS_SPACE.md',
    'uncut/ls1/ls1_compass_space.py',
    'physics/yc23/YC23_PRIMITIVE_CUT_AND_AGITATION.md',
    '02-relational-response/R14_SCOPE_CORRECTION.md',
    'uncut/ls2/LS2_RESPONSE_TOWER_LIFT.md',
    'uncut/ls2/ls2_response_tower.py',
    'uncut/ls2/test_ls2.py',
)
T42 = dict(
    commit='3cc5a33b05c16d59c90994ddda69dedc0d392424',
    path='theorum/42_cut_graded_lambda_jacobian_tower_theorem.md',
    repository='https://github.com/Parveen117/Recognition-Kernel-Framework',
    sha256='e5304cc970bfa0e18140aa1b5373f16733cd5613e25ae4f29642d10f5b6431ed',
)


def clean(value):
    return value.applyfunc(sp.cancel) if isinstance(value, sp.MatrixBase) else sp.cancel(value)


def zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(zero(z) for z in value)
    return sp.simplify(value) == 0


def matrix_sum(values):
    return sum(values, sp.zeros(2))


def normalize_jets(coefficients):
    """Taylor coefficients -> retained mean, positive size and cut coefficients.

    Inputs are 2x2 matrices in one fixed trivialization, not higher tensors
    without a declared contraction. Positivity of Delta at the base is required.
    """
    if not coefficients or any(z.shape != (2, 2) for z in coefficients):
        raise ValueError('a nonempty list of 2x2 response coefficients is required')
    means = [sp.trace(z)/2 for z in coefficients]
    traceless = [z-a*I for z, a in zip(coefficients, means)]
    delta = [sp.simplify(sum(sp.trace(traceless[j]*traceless[n-j])/2
                            for j in range(n+1))) for n in range(len(coefficients))]
    if delta[0].is_positive is not True:
        raise ValueError('the base discriminant must be certified positive')
    size = [sp.sqrt(delta[0])]
    cuts = [clean(traceless[0]/size[0])]
    for n in range(1, len(coefficients)):
        size.append(sp.simplify((delta[n]-sum(size[j]*size[n-j]
                                             for j in range(1, n)))/(2*size[0])))
        cuts.append(clean((traceless[n]-matrix_sum(size[j]*cuts[n-j]
                                                   for j in range(1, n+1)))/size[0]))
    return means, size, cuts


@lru_cache(None)
def polynomial_numerator():
    h = 1+2*x+y
    return clean(h**3*up8.family()['L'])


def response_coefficients(alpha, beta, order):
    """Exact source ray x=alpha*t,y=beta*t; coefficients include 1/n!."""
    alpha, beta = sp.sympify(alpha), sp.sympify(beta)
    c = 2*alpha+beta
    numerator = polynomial_numerator().subs({x: alpha*t, y: beta*t})
    polynomials = [sp.Poly(z, t) for z in numerator]
    out = []
    for n in range(order+1):
        nn = sp.Matrix(2, 2, [z.nth(n) for z in polynomials])
        correction = matrix_sum(sp.binomial(3, j)*c**j*out[n-j]
                                for j in range(1, min(3, n)+1))
        out.append(clean(nn-correction))
    return out


def reconstruct_from_five(coefficients, c):
    """Reconstruct given the independently known denominator (1+c*t)^3."""
    if len(coefficients) != 5:
        raise ValueError('five Taylor coefficients and the known c are required')
    numerator = []
    for n in range(5):
        numerator.append(matrix_sum(sp.binomial(3, j)*c**j*coefficients[n-j]
                                    for j in range(min(3, n)+1)))
    return clean(matrix_sum(z*t**n for n, z in enumerate(numerator))/(1+c*t)**3)


@lru_cache(None)
def first_compass_jets():
    f = up8.family()
    at_zero = {x: 0, y: 0}
    sigma0 = sp.sqrt(f['delta'].subs(at_zero))
    r1 = (f['w'].diff(x)*S+f['w'].diff(y)/V).subs(at_zero)/sigma0
    phi1 = (f['phix']*S+f['phiy']/V).subs(at_zero)
    return sp.factor(r1), sp.factor(phi1)


def certificate():
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not checks[name]:
            raise AssertionError(name)

    f = up8.family()
    origin = {x: 0, y: 0}
    a, p, q, w = sp.symbols('a p q w', real=True)
    raw = a*I+p*K+q*J+w*R
    check('quadratic response factorization on the declared carrier',
          zero((raw-a*I)**2-(p*p+q*q-w*w)*I))
    ps, qs, ws = sp.symbols('ps qs ws', real=True)
    D = p*p+q*q-w*w
    derivative_D = 2*(p*ps+q*qs-w*ws)
    # Tangency after normalized differentiation, with positive sqrt(D) suppressed.
    tangent_numerator = (ps*K+qs*J+ws*R)*D-(p*K+q*J+w*R)*derivative_D/2
    check('normalization differential is tangent to the cut',
          zero((raw-a*I)*tangent_numerator+tangent_numerator*(raw-a*I)))

    cut_path = sp.sqrt(1+t*t)*K+t*R
    check('valid analytic cut has a derivative outside the real-cut sector',
          zero(cut_path**2-I) and cut_path.diff(t).subs(t, 0) == R and R**2 == -I)
    input_jets = response_coefficients(1, 1, 4)
    means, sizes, cuts = normalize_jets(input_jets)
    check('retained jet lift reconstructs every response coefficient through four',
          all(zero(input_jets[n]-means[n]*I-matrix_sum(sizes[j]*cuts[n-j]
                                                     for j in range(n+1)))
              for n in range(5)))
    check('normalized jet products obey the involution convolution through four',
          all(zero(matrix_sum(cuts[j]*cuts[n-j] for j in range(n+1))
                   -(I if n == 0 else sp.zeros(2))) for n in range(5)))

    ep, eq, ew = sp.symbols('ep eq ew', real=True)
    change = (p+ep)**2+(q+eq)**2-(w+ew)**2-D
    check('sector-error identity has the stated linear and quadratic terms',
          zero(change-2*(p*ep+q*eq-w*ew)-(ep*ep+eq*eq-ew*ew)))

    N = polynomial_numerator()
    h = 1+2*x+y
    check('source response equals a quartic numerator over a cubic denominator',
          zero(N-h**3*f['L']) and
          all(sp.Poly(z, x, y).total_degree() == 4 for z in N))
    generic_coeffs = response_coefficients(A, B, 4)
    reconstructed = reconstruct_from_five(generic_coeffs, 2*A+B)
    check('five known-denominator jets reconstruct the entire actual response',
          zero(reconstructed-f['L'].subs({x: A*t, y: B*t})))
    coeffs = response_coefficients(1, 1, 9)
    check('third-order coefficient recurrence after degree four',
          all(zero(coeffs[n]+9*coeffs[n-1]+27*coeffs[n-2]+27*coeffs[n-3])
              for n in range(5, 10)))
    check('derivative recurrence includes factorial weights',
          all(zero(sp.factorial(n)*coeffs[n]
                   +9*n*sp.factorial(n-1)*coeffs[n-1]
                   +27*n*(n-1)*sp.factorial(n-2)*coeffs[n-2]
                   +27*n*(n-1)*(n-2)*sp.factorial(n-3)*coeffs[n-3])
              for n in range(5, 10)))
    line = f['L'].subs({x: t, y: t})
    expected22 = -(2*t+1)*(4*t**3-81*t**2-54*t-9)/(36*(3*t+1)**3)
    check('actual diagonal-ray lower-right response', zero(line[1, 1]-expected22))
    check('negative pole has nonzero third-order numerator',
          sp.cancel((1+3*t)**3*line[1, 1]).subs(t, -sp.Rational(1, 3)) == sp.Rational(1, 729))
    check('rational reconstruction reaches a regular point beyond the Taylor disk',
          reconstruct_from_five(input_jets, 3).subs(t, 1) == up8.at(f['L'], 1, 1))
    check('Taylor truncation is not substituted for exact continuation',
          sum((z for z in input_jets), sp.zeros(2)) != up8.at(f['L'], 1, 1))
    for name in ('H', 'G'):
        poly = sp.Poly(f[name], x, y)
        check(f'complete inherited positive polynomial {name}',
              all(z > 0 for z in poly.coeffs()) and poly.eval(origin) > 0)
    check('positive-axis discriminant bound uses the exact constant coefficient',
          zero(f['delta']-f['H']/(5184*h**6)) and
          sp.Rational(f['H'].subs(origin), 5184) == sp.Rational(1369, 64))

    r1, phi1 = first_compass_jets()
    check('first signed-depth response jet', zero(r1+4*S/37))
    check('first seam-angle response jet', zero(phi1-(-652*S+420/V)/1369))
    jet_jacobian = sp.Matrix([r1, phi1]).jacobian((S, V))
    check('first jet has full state rank everywhere on the positive patch',
          zero(jet_jacobian.det()-1680/(37**3*V**2)))
    recovered_S = -37*r1/4
    recovered_V = 420/(1369*phi1+652*recovered_S)
    check('first jet decoder recovers both original coordinates exactly',
          zero(recovered_S-S) and zero(recovered_V-V))
    check('zero-value normalized compass is constant despite varying raw size',
          f['w'].subs(origin) == 0 and f['q'].subs(origin)/f['p'].subs(origin) == sp.Rational(12, 35)
          and (S**3/V).diff(S) != 0)
    G2 = clean(jet_jacobian.T*jet_jacobian)
    expected_G2 = sp.Matrix([[447008, 273840/V**2], [273840/V**2, 176400/V**4]])/37**4
    check('second-order pullback metric equals the first-jet Gram matrix', zero(G2-expected_G2))
    check('renormalized metric has a strictly positive determinant',
          zero(G2.det()-(1680/(37**3*V**2))**2))
    F3 = sp.factor(2*r1*jet_jacobian.det())
    check('cubic native curvature is determined by the first response jet',
          zero(F3+13440*S/(37**4*V**2)))
    inherited_leading = sp.cancel(f['fxy']/x).subs(origin)
    check('independent exact source-curvature formula has the same cubic term',
          zero(-S*inherited_leading/V**2-F3))
    xy_jacobian = sp.Matrix([t*S, t/V]).jacobian((S, V)).det()
    check('orientation of the source pullback is retained', zero(xy_jacobian+t*t/V**2))
    check('all-positive curvature factorization is inherited without truncation',
          zero(f['fxy']-10368*x*(4*y+1)*h**3*(4*x+3*y+3)**2*f['G']/f['H']**2))
    check('nonzero finite source curvature after exact continuation',
          -up8.at(f['fxy'], 1, 1) == -sp.Rational(76489856000, 2084237405595601))

    # A two-parameter higher-order onset control: m=2, n=3.
    ss, vv = sp.symbols('ss vv', real=True)
    depth, angle = t**2*ss, t**3*vv
    control = 2*depth*(depth.diff(ss)*angle.diff(vv)-depth.diff(vv)*angle.diff(ss))
    check('general spatial onset exponent is twice depth order plus angle order',
          control == 2*t**7*ss)

    return dict(
        stage='LS2', base_commit='625434c', exact_check_count=len(checks), checks=checks,
        source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCES},
        external_source=T42,
        response_recurrence='A_n=-3c A_(n-1)-3c^2 A_(n-2)-c^3 A_(n-3), n>=5, c=2S+1/V',
        first_jet=dict(depth=str(r1), angle=str(phi1), state_jacobian=str(sp.factor(jet_jacobian.det()))),
        metric_leading=[[str(z) for z in row] for row in G2.tolist()],
        curvature_leading=str(F3),
        claims=dict(lossless_admissible_response_jet_lift=True,
                    actual_family_rational_tower=True,
                    all_positive_source_axis_admissible=True,
                    first_jet_repairs_two_state_directions=True,
                    metric_second_curvature_third_order=True,
                    derivative_layers_are_independent_cuts=False,
                    capacity_recursion_identified_with_derivative_tower=False,
                    lambda_is_physical_RG_scale=False,
                    full_primitive_carrier_selected=False,
                    new_yang_mills_gap=False, continuum_mass_gap=False),
        proof_scope='written all-order and analytic proofs, exact algebra controls; not formal verification',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = HERE/'LS2_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise AssertionError('result or source pins differ')
    print(f"LS2: {result['exact_check_count']} exact checks pass.")
    print('Actual response tower reconstructed; first jet repairs both state directions.')
