"""Direct exact controls for YC34's new source/support algebra only."""

from itertools import combinations
from math import comb
import sympy as s


def creator(n, support, amplitudes):
    """Vacuum-to-excited map, with identity on untouched two-level factors."""
    ans = s.zeros(2**n)
    mask = sum(1 << j for j in support)
    for col in range(2**n):
        if col & mask == 0:
            for pattern, value in amplitudes.items():
                assert pattern & ~mask == 0
                ans[col | pattern, col] += value
    return ans


def gamma(n, source):
    ans = s.zeros(2**n)
    for pattern, value in enumerate(source):
        support = [j for j in range(n) if pattern & (1 << j)]
        ans += creator(n, support, {pattern: value})
    return ans


def normal(n, a):
    assert a == a.conjugate().T and a[0, 0] == 0
    lift = gamma(n, a[:, 0])
    return a - lift - lift.conjugate().T


def restricted_psd(a):
    assert a == a.conjugate().T
    for size in range(1, a.rows + 1):
        for indices in combinations(range(a.rows), size):
            assert a.extract(indices, indices).det() >= 0


# A local-vacuum sandwich fails to preserve a cancellation across supports.
n = 2
c0 = creator(n, [0], {1: 1})
x0 = c0 + c0.T
q0 = s.diag(0, 1, 0, 1)
q01 = s.diag(0, 1, 1, 1)
wrong = q0*x0*q0 + q01*(-x0)*q01
assert wrong != s.zeros(4)
assert wrong[2, 3] == -1 and wrong[3, 2] == -1
assert normal(n, x0) + normal(n, -x0) == s.zeros(4)
# Enlarging a source support leaves its creator's spectator identity intact.
assert gamma(n, x0[:, 0]) == c0

# Nontrivial source cancellation across three overlapping support patches.
n = 3
c = [creator(n, [j], {1 << j: 1}) for j in range(n)]
x = [a + a.T for a in c]
q = [s.diag(*[int(bool(i & (1 << j))) for i in range(8)])
     for j in range(n)]
patches = [x[0] + q[0]*x[1],
           -x[0] + x[2] + q[1]*x[2],
           -x[2] + q[2]*x[0]]
total = sum(patches, s.zeros(8))
assert total[:, 0] == s.zeros(8, 1)
ordered = [normal(n, a) for a in patches]
assert all(a == a.T and a[:, 0] == s.zeros(8, 1) for a in ordered)
assert sum(ordered, s.zeros(8)) == total

# Creator products: disjoint supports compose, overlapping ones vanish.
assert c[0]*c[2] == creator(n, [0, 2], {5: 1})
assert c[0]*c[2] == c[2]*c[0]
assert c[0]*creator(n, [0, 2], {5: 1}) == s.zeros(8)

# The finite-model Hilbert bound, with complex coefficients and all supports.
v = s.Matrix([1, 2, -1, s.I, 3, -2*s.I, 1, -1]) / 7
g = gamma(n, v)
vnorm2 = (v.conjugate().T*v)[0]
for k in (0, 1, 2):
    indices = [j for j in range(8) if j.bit_count() <= k]
    gp = g[:, indices]
    d = sum(comb(n, j) for j in range(k + 1))
    restricted_psd(d*vnorm2*s.eye(len(indices)) - gp.conjugate().T*gp)

# Even vacuum-annihilating errors need a factor-count budget.
for volume in (2, 5):
    for k in range(volume + 1):
        restricted_norm = max(j.bit_count() for j in range(2**volume)
                              if j.bit_count() <= k)
        assert restricted_norm == min(k, volume)

# Zero mixed coupling: energies add before normalized functional calculus.
tau = s.Rational(4)
pair_return = s.Rational(2 + 3)/(2 + 3 + tau)
product_return = s.Rational(2)/(2 + tau)*s.Rational(3)/(3 + tau)
assert pair_return == s.Rational(5, 9)
assert product_return == s.Rational(1, 7)
assert pair_return != product_return
assert s.Rational(1, 2 + 3) != s.Rational(1, 2)*s.Rational(1, 3)

# Exponents from the incidence, local lift and output-number factors.
nu, k = 32, 2
summand_power = 3 + (7 - nu) + s.Rational(3*k, 2) + s.Rational(3, 2)
assert summand_power + 1 == -s.Rational(33, 2)

print("Spectator cancellation, overlap products and symmetric source lift: pass")
print("Finite-model Hilbert bound and excitation-count dependence: pass")
print("Additive zero return and k=2 spatial-tail exponent: pass")
print("No numerical Yang--Mills radius, gap or harmonic truncation was tested.")
