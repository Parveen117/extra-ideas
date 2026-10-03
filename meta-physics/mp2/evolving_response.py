"""MP-2 rational controls; reuse the pinned native engine, Python 3.12 only.

The time-continuous assertions are proved in THEOREM.md. No floating point,
ODE integrator, physical calibration or replacement operator engine here.
"""
from fractions import Fraction as F
import importlib
from math import isqrt
from pathlib import Path
import sys

u = None


def configure(publications_root):
    global u
    directory = Path(publications_root).resolve() / 'papers/yang-mills-certified-benchmark/certificates'
    if not (directory / 'ym54_response_protocol.py').is_file():
        raise ValueError('Pinned Publications checkout required')
    sys.path.insert(0, str(directory))
    module = importlib.import_module('ym54_response_protocol')
    if Path(module.__file__).resolve().parent != directory:
        raise ValueError('An engine from a different checkout is already loaded')
    u = module


def inner(A, B):
    A, B = u.matrix2(A), u.matrix2(B)
    return sum((x*y for row, other in zip(A, B) for x, y in zip(row, other)), F(0))


def cofactor(V):
    V = u.matrix2(V)
    return [[V[1][1], -V[1][0]], [-V[0][1], V[0][0]]]


def state(V):
    V = u.matrix2(V)
    B, d = inner(V, V), u.det2(V)
    if not B:
        raise ValueError('Nonzero response required')
    plus = u.a.scale(u.a.add(V, cofactor(V)), F(1, 2))
    minus = u.a.scale(u.a.add(V, cofactor(V), -1), F(1, 2))
    return dict(V=V, B=B, d=d, plus=plus, minus=minus,
                a_plus=inner(plus, plus), a_minus=inner(minus, minus))


def velocity(V, kappa=1, omega=0, error=None):
    s = state(V)
    kappa, omega = u.exact(kappa), u.exact(omega)
    if kappa < 0:
        raise ValueError('Relaxation requires kappa>=0')
    radial = u.a.add(s['V'], u.a.scale(cofactor(V), -2*s['d']/s['B']))
    out = u.a.add(u.a.scale(u.mm(u.R, s['V']), omega), u.a.scale(radial, -kappa))
    return out if error is None else u.a.add(out, u.matrix2(error))


def balances(V, kappa=1, omega=0, error=None):
    s = state(V)
    v = velocity(V, kappa, omega, error)
    loss = u.exact(kappa)*(s['B']-4*s['d']**2/s['B'])
    return dict(d_dot=inner(cofactor(V), v), B_dot=2*inner(V, v),
                m_dot=loss, open_energy=inner(V, v)+loss)


def hessian_jets(V):
    V = u.matrix2(V)
    (p1, p2), (q1, q2) = V
    return [[[2*p1+q2, q1], [q1, q2]], [[q1, q2], [q2, q1-2*p2]]]


def hessian(V, x=0, y=0):
    A, B = hessian_jets(V)
    return u.a.add(u.I, u.a.add(u.a.scale(A, u.exact(x)), u.a.scale(B, u.exact(y))))


def response(V):
    return u.response_data(u.I, hessian_jets(V), u.I)


def tensor(V):
    V = u.matrix2(V)
    S = u.a.scale(u.mm(V, u.a.transpose(V)), F(1, 2))
    return [row+[F(0)] for row in S]+[[F(0)]*3]


def tensor_dot(V, v):
    V, v = u.matrix2(V), u.matrix2(v)
    S = u.a.scale(u.a.add(u.mm(v, u.a.transpose(V)),
                         u.mm(V, u.a.transpose(v))), F(1, 2))
    return [row+[F(0)] for row in S]+[[F(0)]*3]


def floor(d_floor, B_ceiling):
    d, B = u.exact(d_floor), u.exact(B_ceiling)
    if d <= 0 or B < 2*d:
        raise ValueError('Require D>0 and B>=2D')
    beta = d*d/(2*B)
    return dict(even=beta, full=min(d/4, beta))


def robust_budget(V, r0, epsilon):
    s = state(V)
    r0, epsilon = u.exact(r0), u.exact(epsilon)
    if r0 <= 0 or r0*r0 < s['B'] or epsilon < 0:
        raise ValueError('Invalid initial norm or integrated error budget')
    D = abs(s['d'])-r0*epsilon-epsilon*epsilon/2
    B = (r0+epsilon)**2
    return dict(D=D, B=B, **floor(D, B))


def rotation(O):
    O = u.matrix2(O)
    if u.mm(u.a.transpose(O), O) != u.I or u.det2(O) != 1:
        raise ValueError('Oriented orthogonal frame required')
    return O


def pair_step(V, scale_plus, scale_minus, O=None):
    """Check a rational point on a finite relaxation orbit, not an ODE stepper."""
    s = state(V)
    if not s['d']:
        raise ValueError('Nonzero oriented area required')
    ap, am = u.exact(scale_plus), u.exact(scale_minus)
    if not (0 < ap <= 1 and 0 < am <= 1):
        raise ValueError('Positive finite-time shrink factors required')
    O = rotation(u.I if O is None else O)
    W = u.mm(O, u.a.add(u.a.scale(s['plus'], ap), u.a.scale(s['minus'], am)))
    t = state(W)
    if t['d'] != s['d']:
        raise ValueError('Pair scales do not preserve oriented area')
    m = (s['B']-t['B'])/2
    return dict(V=W, O=O, m=m, B=t['B'], d=t['d'])


def rational_root(x):
    x = u.exact(x)
    if x < 0:
        raise ValueError('Nonnegative root required')
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    if a*a != x.numerator or b*b != x.denominator:
        raise ValueError('This finite checker requires a rational square root')
    return F(a, b)


def recover_rational(V, O, m):
    """Exact finite witness for T3; general roots use the written completion."""
    m = u.exact(m)
    if m < 0:
        raise ValueError('Nonnegative retained record required')
    W = u.mm(u.a.transpose(rotation(O)), u.matrix2(V))
    s = state(W)
    result = u.a.scale(u.I, 0)
    for branch, weight in [('plus', 'a_plus'), ('minus', 'a_minus')]:
        a = s[weight]
        if not a:
            if m:
                raise ValueError('A lost direction cannot be reconstructed from a scalar weight')
        else:
            result = u.a.add(result, u.a.scale(s[branch], rational_root((a+m)/a)))
    return result


def gamma_linear(f, g, C):
    """Linear product form for symmetric C, including indefinite Cdot."""
    if len(C) != 3 or any(len(row) != 3 for row in C):
        raise ValueError('Three native generator coefficients required')
    C = [[u.exact(x) for x in row] for row in C]
    if C != u.a.transpose(C):
        raise ValueError('Symmetric coefficients required')
    result = {}
    for i in range(3):
        for j in range(3):
            result = u.p.add(result, u.p.scale(u.p.mul(u.y.native_generator(f, i+1),
                               u.y.native_generator(g, j+1)), C[i][j]))
    return result


def work_balance(f, C, Cdot):
    Lf = u.q.lap(f, C)
    energy = u.z.energy(f, C)
    work = u.y.phi(gamma_linear(f, f, Cdot))
    dissipation = 2*u.y.phi(u.p.mul(Lf, Lf))
    return dict(energy=energy, work=work, dissipation=dissipation,
                energy_dot=work-dissipation, squared_norm_dot=-2*energy)


def work_witness():
    V = u.a.diag((4, 1))
    q0, q1 = u.p.COORD[:2]
    f = u.p.add(u.p.add(u.p.mul(q0, q0), u.p.mul(q1, q1)),
                u.p.scale(u.y.RADIUS, F(-1, 2)))
    C = tensor(V)
    Cdot = tensor_dot(V, velocity(V, 10))
    return f, C, Cdot, work_balance(f, C, Cdot)


def resolvent(C, degree, h):
    """Exact rational contraction proxy; this is not the heat exponential."""
    h = u.exact(h)
    if h <= 0:
        raise ValueError('Positive step required')
    L = u.q.operator(degree, C)
    return u.a.inverse(u.a.add(u.a.identity(len(L)), u.a.scale(L, h)))


def exact_examples():
    V = u.a.diag((4, 1))
    O = [[F(3, 5), F(-4, 5)], [F(4, 5), F(3, 5)]]
    step = pair_step(V, F(41, 50), F(3, 10), O)
    recovered = recover_rational(step['V'], O, step['m'])
    return dict(initial=dict(d=F(4), B=F(17)), pair_step=step,
                initial_recovered=recovered == V, full_floor=floor(4, 17)['full'],
                moving_work=work_witness()[3], robust=robust_budget(V, 5, F(1, 10)),
                balanced_full=F(1), nearby_full=F(881, 800),
                optimal_full_squared=F(4, 3))
