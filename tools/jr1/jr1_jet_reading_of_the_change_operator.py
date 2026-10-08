"""JR1: the jet reading of the change operator.

Source object: the "lambda-residue"  Res_lam(f, z0) = e^(-i(lam-1) z0) Res(f e^(i(lam-1) z), z0)  of the
catalog draft, stated there for "lambda-analytic" f  (C f = lam f,  C = I + i d/dz).
Sources read before building: RKF F00-E (cut-complex field, iota^2 = -1, factorial polynomials E_N),
Publications GE1-T2 (the change operator C = I + D and its kernel), RKF theorum/46 (first visible jet).

Carrier: Laurent tails at a point over the exact cut-complex rationals (rad + iota turn),
    f = sum_(j=1..m) a_(-j) t^(-j) + (regular part) ,   t = z - z0 ,   D t^k = k t^(k-1) .
No contour and no limit is used: every sum below is finite.

    reading      Rd_c(f) = sum_(j>=1) a_(-j) c^(j-1)/(j-1)!          (c = iota (lam - 1))

J1  Rd_c(f) is the t^(-1) coefficient of f * E_N(c t), any N >= m-1: it is the draft's expression.
J2  Rd_c(D f) = -c Rd_c(f).  Hence  Rd_c(C f) = lam Rd_c(f)  and  Rd_c(p(C) f) = p(lam) Rd_c(f):
    the reading is an eigen-reading of the change operator on every tail, with eigenvalue lam.
J3  At a simple pole Rd_c(f) = a_(-1) for every lam: one jet carries no lam.
J4  A pole of order m is recovered from m readings at distinct lam, and not from fewer.
    The single reading at lam = 1 is the first jet a_(-1) alone.
J5  The eigenfunctions of C (C f = lam f) are multiples of Exp(-c z); a tail with a pole is never one
    ((C - lam) f has the term -iota m a_(-m) t^(-m-1)).  So on the draft's stated domain every reading is zero.
    The reading and the eigenfunction are partners: Rd_c(f) is the first jet of f Exp(c t).
J6  Shift:  Rd_c(f E(b t)) = Rd_(c+b)(f).   Resolvent (a corollary of J2):  (C - mu) g = f  gives
    Rd_c(g) = Rd_c(f)/(lam - mu).

Python 3.12, standard library only.
"""
from fractions import Fraction as F
from math import factorial
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


class G:
    """exact cut-complex rational  rad + iota*turn ,  iota^2 = -1"""
    __slots__ = ('r', 't')

    def __init__(self, r=0, t=0):
        self.r, self.t = F(r), F(t)

    def __add__(self, o):
        o = g(o)
        return G(self.r + o.r, self.t + o.t)
    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.t)

    def __sub__(self, o):
        return self + (-g(o))

    def __rsub__(self, o):
        return g(o) - self

    def __mul__(self, o):
        o = g(o)
        return G(self.r*o.r - self.t*o.t, self.r*o.t + self.t*o.r)
    __rmul__ = __mul__

    def inv(self):
        n = self.r*self.r + self.t*self.t
        return G(self.r/n, -self.t/n)

    def __truediv__(self, o):
        return self*g(o).inv()

    def __pow__(self, n):
        out = G(1)
        for _ in range(n):
            out = out*self
        return out

    def __eq__(self, o):
        o = g(o)
        return self.r == o.r and self.t == o.t

    def __hash__(self):
        return hash((self.r, self.t))

    def __repr__(self):
        return f'({self.r} + {self.t} iota)'


def g(x):
    return x if isinstance(x, G) else G(x)


IOTA = G(0, 1)


# ---------------------------------------------------------------- Laurent tails: dict exponent -> coefficient
def clean(f):
    return {k: v for k, v in f.items() if v != 0}


def lmul(f, h):
    out = {}
    for i, x in f.items():
        for j, y in h.items():
            out[i + j] = out.get(i + j, G(0)) + x*y
    return clean(out)


def ladd(f, h, s=1):
    out = dict(f)
    for k, v in h.items():
        out[k] = out.get(k, G(0)) + s*v
    return clean(out)


def lscale(f, c):
    return clean({k: c*v for k, v in f.items()})


def deriv(f):
    return clean({k - 1: k*v for k, v in f.items() if k != 0})


def exp_poly(c, n):
    """E_N(c t) = sum_(k<=N) c^k t^k / k!   (RKF F00-E factorial polynomial)"""
    return {k: (g(c)**k)/factorial(k) for k in range(n + 1)}


def order_of_pole(f):
    neg = [-k for k in f if k < 0]
    return max(neg) if neg else 0


def reading(f, c):
    c = g(c)
    return sum((f[-j]*(c**(j - 1))/factorial(j - 1) for j in range(1, order_of_pole(f) + 1) if -j in f), G(0))


def c_of(lam):
    return IOTA*(g(lam) - 1)


def change(f):
    """C f = f + iota D f"""
    return ladd(f, lscale(deriv(f), IOTA))


def random_tail(rng, m, reg=4):
    f = {}
    for k in range(-m, reg + 1):
        f[k] = G(rng.randint(-5, 5), rng.randint(-5, 5))
    if f[-m] == 0:
        f[-m] = G(1)
    return clean(f)


def solve(mat, rhs):
    """exact Gauss elimination over the cut-complex rationals; returns None if singular"""
    n = len(mat)
    a = [row[:] + [rhs[i]] for i, row in enumerate(mat)]
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            return None
        a[col], a[piv] = a[piv], a[col]
        inv = a[col][col].inv()
        a[col] = [x*inv for x in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                fac = a[r][col]
                a[r] = [x - fac*y for x, y in zip(a[r], a[col])]
    return [a[i][n] for i in range(n)]


def run():
    out = {}
    rng = random.Random(20261009)
    tails = [random_tail(rng, m) for m in (1, 2, 3, 4, 5) for _ in range(6)]
    lams = [G(2), G(0), G(F(1, 2), F(3, 4)), G(1, F(1, 3)), G(-2, 5)]     # includes the draft's lam = 1 + iota/m

    # J1: two routes
    out['J1_reading_is_the_first_jet_of_f_times_factorial_polynomial'] = all(
        reading(f, c_of(l)) == lmul(f, exp_poly(c_of(l), order_of_pole(f) + extra)).get(-1, G(0))
        for f in tails for l in lams for extra in (-1, 0, 3))

    # J2: eigen-reading
    out['J2_reading_of_the_derivative'] = all(
        reading(deriv(f), c_of(l)) == -c_of(l)*reading(f, c_of(l)) for f in tails for l in lams)
    out['J2_change_operator_acts_as_lambda'] = all(
        reading(change(f), c_of(l)) == l*reading(f, c_of(l)) for f in tails for l in lams)
    ok = True
    for f in tails[:10]:
        for l in lams:
            cf = change(f)
            ccf = change(cf)
            pf = ladd(ladd(lscale(ccf, G(3)), lscale(cf, G(-2, 1))), lscale(f, G(5)))    # p(C) = 3C^2 + (-2+iota)C + 5
            ok &= reading(pf, c_of(l)) == (3*l*l + G(-2, 1)*l + 5)*reading(f, c_of(l))
    out['J2_polynomials_of_the_change_operator'] = ok

    # two routes for D on an actual rational function:  f = 1/((z-z0)^2 (z-z1)),  d = z1 - z0
    d = G(F(3, 2), -2)
    nreg = 8
    fa = clean({n - 2: -(d**(n + 1)).inv() for n in range(nreg)})
    # quotient rule:  f' = -2/((z-z0)^3 (z-z1)) - 1/((z-z0)^2 (z-z1)^2)
    fb = ladd(clean({n - 3: 2*(d**(n + 1)).inv() for n in range(nreg)}),
              clean({n - 2: -(n + 1)*(d**(n + 2)).inv() for n in range(nreg)}))
    polar = lambda f: {k: v for k, v in f.items() if k < 0}
    out['J2_derivative_rule_agrees_with_the_quotient_rule'] = polar(deriv(fa)) == polar(fb)
    out['J2_rational_function_example'] = all(
        reading(change(fa), c_of(l)) == l*reading(fa, c_of(l)) for l in lams)

    # J3: one jet carries no lambda
    simple = [f for f in tails if order_of_pole(f) == 1]
    out['J3_simple_pole_reading_is_the_first_jet_for_every_lambda'] = all(
        reading(f, c_of(l)) == f[-1] for f in simple for l in lams)
    double = [f for f in tails if order_of_pole(f) == 2]
    out['J3_higher_poles_do_depend_on_lambda'] = all(
        reading(f, c_of(lams[0])) != reading(f, c_of(lams[1])) for f in double)

    # J4: recovery from m readings, not from fewer
    ok = True
    for f in tails:
        m = order_of_pole(f)
        cs = [c_of(l) for l in lams[:m]]
        mat = [[(c**(j - 1))/factorial(j - 1) for j in range(1, m + 1)] for c in cs]
        sol = solve(mat, [reading(f, c) for c in cs])
        ok &= sol == [f.get(-j, G(0)) for j in range(1, m + 1)]
    out['J4_m_readings_recover_a_pole_of_order_m'] = ok
    c1, c2 = c_of(lams[0]), c_of(lams[2])
    # a non-zero tail of order 3 whose readings at c1 and c2 both vanish:  reading = (c - c1)(c - c2)
    hidden = {-3: G(2), -2: -(c1 + c2), -1: c1*c2}
    out['J4_fewer_readings_do_not'] = (reading(hidden, c1) == 0 and reading(hidden, c2) == 0
                                        and clean(hidden) != {})
    out['J4_reading_at_lambda_one_is_the_first_jet_only'] = all(
        reading(f, c_of(G(1))) == f.get(-1, G(0)) for f in tails)

    # J5: eigenfunctions of C have no pole
    ok = True
    for l in lams:
        c = c_of(l)
        n = 12
        e = exp_poly(-c, n)                                            # truncated Exp(-c t)
        lhs, rhs = change(e), lscale(e, l)
        ok &= all(lhs.get(k, G(0)) == rhs.get(k, G(0)) for k in range(n))    # exact below the cut-off
    out['J5_exp_of_minus_c_is_an_eigenfunction'] = ok
    # no tail with a pole is an eigenfunction: (C - lam) f has the term  -iota m a_(-m) t^(-m-1)
    ok = True
    for f in tails:
        m = order_of_pole(f)
        for l in lams:
            defect = ladd(change(f), lscale(f, l), -1)
            ok &= defect.get(-m - 1, G(0)) == -IOTA*m*f[-m] and defect.get(-m - 1, G(0)) != 0
    out['J5_a_tail_with_a_pole_is_never_an_eigenfunction'] = ok
    # a general regular series with C f = lam f is fixed by its value at t = 0
    ok = True
    for l in lams:
        c = c_of(l)
        a = [G(7, -3)]
        for k in range(1, 10):                                         # (k) a_k * iota = (lam - 1) a_(k-1)
            a.append((g(l) - 1)*a[k - 1]/(IOTA*k))
        ok &= all(a[k] == a[0]*((-c)**k)/factorial(k) for k in range(10))
    out['J5_eigenfunction_is_unique_up_to_a_factor'] = ok

    # J6: shift and resolvent
    b = G(F(2, 3), -1)
    out['J6_shift_law'] = all(
        reading(lmul(f, exp_poly(b, order_of_pole(f) + 2)), c_of(l)) == reading(f, c_of(l) + b)
        for f in tails for l in lams)
    mu = G(F(1, 5), 2)
    out['J6_resolvent'] = all(
        reading(gg, c_of(l)) == reading(ladd(change(gg), lscale(gg, mu), -1), c_of(l))/(l - mu)
        for gg in tails for l in lams)

    # controls ----------------------------------------------------------------------------------------------
    f = double[0]
    l = lams[2]
    no_factorial = lambda f, c: sum((f[-j]*(g(c)**(j - 1)) for j in range(1, order_of_pole(f) + 1) if -j in f), G(0))
    f5 = [t for t in tails if order_of_pole(t) == 5][0]
    out['control_without_factorial_weights_the_law_fails'] = (
        no_factorial(deriv(f5), c_of(l)) != -c_of(l)*no_factorial(f5, c_of(l)))
    out['control_opposite_sign_of_c_fails'] = reading(change(f), -c_of(l)) != l*reading(f, -c_of(l))

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out


if __name__ == '__main__':
    out = run()
    with open(os.path.join(HERE, 'JR1_RESULT.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
