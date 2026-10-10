"""UP9: SU(2) Haar source selects the deformed compass.

Exact cumulants, written analytic identities, and outward rational enclosures
of the full integral. The lambda flow is source dilation, not spatial RG.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S, V, lam = sp.symbols('S V lambda', positive=True)


@lru_cache(None)
def haar(n):
    if not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer moment order required')
    return Q(0) if n % 2 else Q(comb(n, n//2), (n//2+1)*4**(n//2))


@lru_cache(None)
def mixed(i, j):
    """E[c^i (c^2-1/4)^j]."""
    if min(i, j) < 0:
        raise ValueError('nonnegative mixed powers required')
    return sum((Q(comb(j, k))*(-Q(1, 4))**(j-k)*haar(i+2*k)
                for k in range(j+1)), Q(0))


@lru_cache(None)
def cumulant_potential(order=6):
    """Taylor coefficients of U_lambda through lambda^(order-2)."""
    mu = [sp.expand(sum(sp.binomial(n, j)*sp.Rational(mixed(n-j, j))
                       *S**(n-j)*V**j for j in range(n+1)))
          for n in range(order+1)]
    cumulants = [sp.S(0)]*(order+1)
    for n in range(1, order+1):
        cumulants[n] = sp.expand(mu[n]-sum(
            sp.binomial(n-1, j-1)*cumulants[j]*mu[n-j] for j in range(1, n)))
    return sp.expand(sum(cumulants[n]*lam**(n-2)/sp.factorial(n)
                         for n in range(2, order+1)))


@lru_cache(None)
def leading_geometry():
    U = cumulant_potential()
    T = sp.diff(U, S)
    m = sp.series(-V*sp.diff(T, V)/(S*sp.diff(T, S)), lam, 0, 5).removeO().expand()
    w = sp.series(S*sp.diff(m, S)*S*T/2, lam, 0, 5).removeO().expand()
    ds = lambda f: S*sp.diff(f, S)
    dv = lambda f: V*sp.diff(f, V)+m*ds(f)
    # Only the constant diagonal and first off-diagonal coefficients are needed.
    p0 = sp.factor((ds(ds(U))-dv(dv(U))).subs(lam, 0)/2)
    q1 = sp.factor(sp.expand(dv(ds(U))+ds(dv(U))).coeff(lam, 1)/2)
    ell8 = sp.factor(w.coeff(lam, 4)**2/p0**2)
    phi1 = sp.factor(q1/p0)
    f9 = sp.factor(sp.diff(ell8, S)*sp.diff(phi1, V)
                   -sp.diff(ell8, V)*sp.diff(phi1, S))
    return dict(U=U, m=m, w=w, p0=p0, phi1=phi1, ell8=ell8, f9=f9)


# Every arithmetic operation rounds outwards to this fixed rational grid.
# No floating point is used to obtain a sign or an enclosure.
GRID = 2**160


def floor_grid(x):
    x = Q(x)*GRID
    return Q(x.numerator//x.denominator, GRID)


def ceil_grid(x):
    return -floor_grid(-Q(x))


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('reversed interval')
        object.__setattr__(self, 'lo', floor_grid(self.lo))
        object.__setattr__(self, 'hi', ceil_grid(self.hi))

    @classmethod
    def exact(cls, x):
        return cls(Q(x), Q(x))

    def __add__(self, other):
        b = as_interval(other)
        return Interval(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-as_interval(other)

    def __rsub__(self, other):
        return as_interval(other)+-self

    def __mul__(self, other):
        b = as_interval(other)
        p = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return Interval(min(p), max(p))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval crosses zero')
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*as_interval(other).reciprocal()

    def __rtruediv__(self, other):
        return as_interval(other)*self.reciprocal()

    def contains(self, value):
        return self.lo <= value <= self.hi

    def record(self):
        return [str(self.lo), str(self.hi)]


def as_interval(x):
    return x if isinstance(x, Interval) else Interval.exact(x)


ZERO = Interval.exact(0)
ONE = Interval.exact(1)


def trunc(a, degree):
    return {ij: v for ij, v in a.items() if sum(ij) <= degree}


def add(a, b, degree):
    return {ij: a.get(ij, ZERO)+b.get(ij, ZERO)
            for ij in a.keys() | b.keys() if sum(ij) <= degree}


def scale(a, k):
    return {ij: v*k for ij, v in a.items()}


def mul(a, b, degree):
    out = {}
    for (i, j), v in a.items():
        for (k, l), w in b.items():
            ij = (i+k, j+l)
            if sum(ij) <= degree:
                out[ij] = out.get(ij, ZERO)+v*w
    return out


def diff(a, axis):
    out = {}
    for ij, v in a.items():
        if ij[axis]:
            key = list(ij)
            key[axis] -= 1
            out[tuple(key)] = v*ij[axis]
    return out


def inverse(a, degree):
    inv0 = 1/a[(0, 0)]
    # The constant of the nonconstant part is exactly zero, not an interval
    # subtraction of a0/a0; this is the formal Taylor identity.
    b = scale({ij: v for ij, v in a.items() if sum(ij)}, -inv0)
    term = {(0, 0): ONE}
    out = term
    for _ in range(degree):
        term = mul(term, b, degree)
        out = add(out, term, degree)
    return scale(out, inv0)


def divide(a, b, degree):
    return mul(a, inverse(b, degree), degree)


@lru_cache(None)
def tilted_moment(i, j, lambda0, s0=Q(1), v0=Q(1), order=24):
    """Full Haar integral enclosed by polynomial moments + exponential tail."""
    lambda0, s0, v0 = Q(lambda0), Q(s0), Q(v0)
    t = abs(lambda0)*(abs(s0)+3*abs(v0)/4)
    if not 0 < t < 1 or order < 0:
        raise ValueError('this tail implementation requires 0 < source norm < 1')
    val = Q(0)
    for n in range(order+1):
        mn = sum((Q(comb(n, k))*s0**(n-k)*v0**k*mixed(i+n-k, j+k)
                  for k in range(n+1)), Q(0))
        val += lambda0**n*mn/factorial(n)
    # e^t <= 1/(1-t), |c|<=1 and |c^2-1/4|<=3/4.
    tail = Q(3, 4)**j*t**(order+1)/(factorial(order+1)*(1-t))
    return Interval(val-tail, val+tail)


def source_jet(lambda0, s0=Q(1), v0=Q(1), order=24):
    lambda0 = Q(lambda0)
    z0 = tilted_moment(0, 0, lambda0, s0, v0, order)
    b = {(i, j): tilted_moment(i, j, lambda0, s0, v0, order)/z0
         *lambda0**(i+j)/Q(factorial(i)*factorial(j))
         for i in range(5) for j in range(5-i) if i+j}
    # log(Z/Z0), to fourth order in state increments. Constant U is not
    # needed by any response below; it is not asserted to vanish.
    out, power = {}, {(0, 0): ONE}
    for k in range(1, 5):
        power = mul(power, b, 4)
        out = add(out, scale(power, Q((-1)**(k+1), k)), 4)
    return scale(out, 1/lambda0**2)


def geometry_from_jet(U, s0=Q(1), v0=Q(1)):
    sj, vj = {(0, 0): as_interval(s0), (1, 0): ONE}, {(0, 0): as_interval(v0), (0, 1): ONE}
    ds = lambda f, n: mul(sj, diff(f, 0), n)
    T = diff(U, 0)
    m = scale(divide(mul(vj, diff(T, 1), 2), mul(sj, diff(T, 0), 2), 2), -1)
    f = ds(U, 2)
    g = add(mul(vj, diff(U, 1), 2), mul(m, f, 2), 2)
    dv = lambda f: add(mul(vj, diff(f, 1), 1), mul(m, ds(f, 1), 1), 1)
    L = [[ds(f, 1), dv(f)], [ds(g, 1), dv(g)]]
    p = scale(add(L[0][0], scale(L[1][1], -1), 1), Q(1, 2))
    q = scale(add(L[0][1], L[1][0], 1), Q(1, 2))
    w = scale(add(L[1][0], scale(L[0][1], -1), 1), Q(1, 2))
    rho2 = add(mul(p, p, 1), mul(q, q, 1), 1)
    delta = add(rho2, scale(mul(w, w, 1), -1), 1)
    ell = divide(mul(w, w, 1), delta, 1)
    phi = [(p[(0, 0)]*q.get(ij, ZERO)-q[(0, 0)]*p.get(ij, ZERO))/rho2[(0, 0)]
           for ij in ((1, 0), (0, 1))]
    curvature = ell.get((1, 0), ZERO)*phi[1]-ell.get((0, 1), ZERO)*phi[0]
    uss, usv, uvv = U[(2, 0)]*2, U[(1, 1)], U[(0, 2)]*2
    # An independent response formula for the antisymmetric source.
    direct_w = mul(ds(m, 1), f, 1).get((0, 0), ZERO)/2
    return dict(delta=delta[(0, 0)], ell=ell[(0, 0)], w=w[(0, 0)],
                p=p[(0, 0)], q=q[(0, 0)], curvature=curvature,
                hessian_det=uss*uvv-usv*usv, uss=uss, direct_w=direct_w)


@lru_cache(None)
def finite_geometry(lambda0, order=24):
    return geometry_from_jet(source_jet(Q(lambda0), order=order))


def symbolic_checks():
    g = leading_geometry()
    U, m, w = g['U'], g['m'], g['w']
    d = 4*S**2-V**2
    z = lambda a: sp.cancel(a) == 0
    checks = {}
    checks['normalized Haar and centered mixed covariance'] = (
        [haar(n) for n in range(9)] == [Q(1), Q(0), Q(1,4), Q(0), Q(1,8), Q(0), Q(5,64), Q(0), Q(7,128)]
        and mixed(0,1) == mixed(1,1) == 0 and mixed(0,2) == Q(1,16))
    checks['Gaussian response is the source-selected flat limit'] = z(U.coeff(lam,0)-(4*S*S+V*V)/32)
    checks['cubic source coefficient is fixed by Haar moments'] = z(U.coeff(lam,1)-V*(12*S*S+V*V)/384)
    checks['fourth cumulant coefficient'] = z(U.coeff(lam,2)+S*S*(2*S*S-3*V*V)/768)
    checks['fifth cumulant coefficient'] = z(U.coeff(lam,3)+V*(60*S**4+V**4)/30720)
    checks['sixth cumulant coefficient'] = z(U.coeff(lam,4)-(32*S**6-192*S**4*V**2-24*S**2*V**4-V**6)/294912)
    checks['exact source dilation at every retained degree'] = z(lam*sp.diff(U,lam)+2*U-S*sp.diff(U,S)-V*sp.diff(U,V))
    regular_heat = sp.expand(lam*sp.diff(U,V)-sp.diff(U,S,2)+sp.Rational(1,4)-lam**2*sp.diff(U,S)**2)
    checks['centered heat source identity through all available orders'] = all(z(regular_heat.coeff(lam,j)) for j in range(5))
    checks['Haar integration by parts determines moment recursion'] = all(
        (n+3)*haar(n+1) == n*haar(n-1) for n in range(1,18))
    checks['frame response cancels to fourth order in its S dependence'] = z(
        m+lam*V/4-lam**3*V**3/128+lam**4*V**2*d/1536)
    checks['first cut defect has exact fourth-order source'] = z(w+lam**4*S**4*V**2/1536)
    checks['flat cut discriminant avoids V equals 2S'] = z(g['p0']-d/16)
    checks['first seam angle coefficient'] = z(g['phi1']+S*S*V/d)
    checks['first lost weight is eighth order'] = z(g['ell8']-S**8*V**4/(9216*d**2))
    checks['first native curvature is ninth order'] = z(g['f9']+S**9*V**4*(2*S*S+V*V)/(1152*d**4))
    checks['nonzero rational curvature coefficient at unit source'] = g['f9'].subs({S:1,V:1}) == -sp.Rational(1,31104)
    b = sp.symbols('b', positive=True)
    checks['independent aggregation step holds coefficient by coefficient'] = z(
        U.subs(lam,lam/sp.sqrt(b))-b*U.subs({S:S/sp.sqrt(b),V:V/sp.sqrt(b)}, simultaneous=True))
    centre = sp.expand(U-(S*sp.diff(U,S)+V*sp.diff(U,V))/2)
    checks['diagram centre is exactly the source dilation response'] = z(centre+lam*sp.diff(U,lam)/2)
    checks['centre retains the cubic mixed record before curvature'] = z(centre.coeff(lam,1)+V*(12*S*S+V*V)/768)
    return checks


def source_pins():
    paths = ['physics/yc1/YC1_CENTRE_COMPASS_POTENTIAL.md',
             'physics/yc1/yc1_centre_compass_potential.py',
             'physics/ct1/CT1_CENTRE_RECORD_OF_THE_TURN_CHAIN.md',
             'uncut/up1/UP1_DEGREE_OF_THE_POTENTIAL.md',
             'uncut/up6/UP6_EIGHT_SCALE_OPERATIONS.md',
             'uncut/ss1/SS1_SCALE_IS_THE_SOURCE.md',
             'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
             'uncut/up8/up8_deformed_compass.py',
             'physics/lc1/LC1_LOST_IS_COUNTED.md',
             'uncut/up9/UP9_HAAR_SOURCE_AND_FLAT_LIMIT.md',
             'uncut/up9/up9_haar_source_flow.py', 'uncut/up9/test_up9.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    checks = symbolic_checks()
    finite = []
    for l0 in (Q(1,4), Q(1,2)):
        f = finite_geometry(l0)
        prefix = f'full integral lambda={l0}: '
        checks[prefix+'positive ordinary Hessian enclosure'] = f['uss'].lo > 0 and f['hessian_det'].lo > 0
        checks[prefix+'positive cut discriminant enclosure'] = f['delta'].lo > 0
        checks[prefix+'signed cut defect excludes zero'] = f['w'].hi < 0
        checks[prefix+'lost weight strictly positive'] = f['ell'].lo > 0
        checks[prefix+'native curvature excludes zero'] = f['curvature'].hi < 0
        checks[prefix+'source commutator and matrix readings overlap'] = (
            max(f['w'].lo,f['direct_w'].lo) <= min(f['w'].hi,f['direct_w'].hi))
        finite.append({'lambda':str(l0), 'S':'1', 'V':'1', 'exponential_order':24,
                       'outward_binary_bits':160, 'enclosures':{k:v.record() for k,v in f.items()}})
    checks = {k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k,v in checks.items() if not v])
    g = leading_geometry()
    return dict(stage='UP9', base_commit='78157ee', checks=checks,
        potential='U_lambda = log E_Haar exp(lambda*(S*c+V*(c^2-1/4)))/lambda^2',
        coefficients=[str(sp.factor(g['U'].coeff(lam,j))) for j in range(5)],
        geometry={k:str(g[k]) for k in ('m','w','p0','phi1','ell8','f9')},
        finite_integral_certificates=finite,
        claim_boundary=dict(source_is_actual_SU2_Haar=True,
            source_dilation_and_independent_aggregation=True,
            flat_zero_limit_on_nondegenerate_cut_patch=True,
            nonzero_curvature_for_sufficiently_small_positive_lambda=True,
            all_positive_lambda_curvature_proved=False,
            finite_integrals_include_full_exponential_tails=True,
            lambda_is_YM_coupling_or_spatial_RG=False,
            new_YM_gap_or_full_SU2_connection=False,
            formal_machine_verified_proof=False), source_pins=source_pins())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'UP9_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"UP9: {len(result['checks'])} exact/outward checks pass.")
    print('SU(2) source-selected flat limit; full-integral curvature witnesses; no new spectral gap.')
