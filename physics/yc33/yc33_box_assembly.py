"""Exact YC33 convention and coefficient-assembly controls only."""

from fractions import Fraction as Fr
from itertools import combinations
import sympy as s


def equal(a, b):
    assert all(s.simplify(x) == 0 for x in a - b)


def psd(a):
    equal(a, a.conjugate().T)
    for size in range(1, a.rows + 1):
        for indices in combinations(range(a.rows), size):
            assert s.simplify(a.extract(indices, indices).det()) >= 0


# Moving vacuum split, with nonzero diagonal energy change.
x = s.symbols("x", real=True)
rot = s.Matrix([[1 - x*x, -2*x], [2*x, 1 - x*x]]) / (1 + x*x)
diag = s.diag(x, 5 + 2*x)
h = rot * diag * rot.T
pi = rot * s.diag(1, 0) * rot.T
d = 2 * s.Matrix([[0, s.I], [-s.I, 0]]) / (1 + x*x)
observed_source = s.diff(h, x) + s.I * (h*d - d*h)
equal(observed_source, rot * s.diag(1, 2) * rot.T)
equal(observed_source*pi, pi*observed_source)
equal(s.diff(pi, x), s.I * (d*pi - pi*d))
equal(s.diff(h, x)*pi - pi*s.diff(h, x),
      -(h*s.diff(pi, x) - s.diff(pi, x)*h))
equal(s.diff(rot.T*h*rot, x), s.diag(1, 2))

# Independent algebraic control of unequal source-box errors.
t = s.Matrix([[s.Rational(3, 10), s.Rational(1, 40), 0],
              [s.Rational(1, 40), s.Rational(1, 4), s.Rational(1, 50)],
              [0, s.Rational(1, 50), s.Rational(1, 5)]])
error = s.Matrix([[0, s.I, 1], [-2*s.I, 0, 1], [0, -1, 0]]) / 100
raw = t + error
assert raw != raw.conjugate().T
epsilon = s.Rational(1, 20)
assert sum(s.simplify(s.conjugate(z)*z) for z in error) <= epsilon**2
sym = (raw + raw.conjugate().T) / 2
positive = sym + epsilon * s.eye(3)
psd(t)
psd(s.eye(3)/2 - t)
psd(positive - t)
psd(t + 2*epsilon*s.eye(3) - positive)
tau = 4
exact = (s.eye(3) + tau*t).inv()
boxed = (s.eye(3) + tau*positive).inv()
psd(exact - boxed)
psd(2*tau*epsilon*s.eye(3) - (exact - boxed))

for label, m in (("28", Fr(592, 175)), ("64", Fr(11188, 3325))):
    # The geometric/boundary proof must deliver these two independent parts.
    tail = Fr(1, 1600)
    box = Fr(1, 1600)
    total = tail + box
    q = 1 - 1 / (1 + 4 * (1/m + 2*total))
    assert total == Fr(1, 800)
    assert q**9 < Fr(1, 200)
    assert 8*total + q**9 < Fr(3, 200)
    print(f"{label}-link: combined budget 1/800, degree-8 error < 3/200")

print("Moving-frame cancellation, asymmetric assembly and positive correction: pass")
print("No open-box spectrum, practical radii or harmonic tail were numerically certified.")
