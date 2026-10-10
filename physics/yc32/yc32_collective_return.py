"""Exact direct controls for YC32; not a Yang--Mills harmonic truncation."""

from fractions import Fraction as Fr
import sympy as s


def equal(a, b):
    assert all(s.simplify(x) == 0 for x in a - b)


def psd2(a):
    equal(a, a.T)
    assert a[0, 0] >= 0 and a[1, 1] >= 0 and a.det() >= 0


# The three coordinates are an independent finite algebraic control model.
a = s.Matrix([[5, -1, 0], [-1, 4, -1], [0, -1, 6]])
v = s.Matrix([[1, 0], [0, 0], [0, 1]])
hidden = s.Matrix([0, 1, 0])
p = v * v.T
qh = s.eye(3) - p
tau = s.Integer(4)
y = a.inv() * v
t = v.T * y
f = a * (a + tau * s.eye(3)).inv()
schur = v.T * f * v - (v.T * f * hidden) * (hidden.T * f * hidden).inv() * (hidden.T * f * v)
reduced = (s.eye(2) + tau * t).inv()
equal(schur, reduced)
assert schur != v.T * f * v

# Radius-zero column localization in this model; all discarded entries kept
# in the exact error calculation. Frobenius gives a rational norm certificate.
yr = s.Matrix([[y[0, 0], 0], [0, 0], [0, y[2, 1]]])
tr = v.T * yr
err = y - yr
epsilon = s.Rational(1, 10)
assert sum(x*x for x in err) <= epsilon**2
tr_plus = tr + epsilon * s.eye(2)
psd2(tr_plus - t)
psd2(t + 2 * epsilon * s.eye(2) - tr_plus)
sr = (s.eye(2) + tau * tr_plus).inv()
psd2(reduced - sr)
psd2(2 * tau * epsilon * s.eye(2) - (reduced - sr))

wf = v + tau * qh * y * reduced
wfr = v + tau * qh * yr * sr
equal(v.T * wf, s.eye(2))
equal(v.T * wfr, s.eye(2))
equal(wf.T * f * wf, reduced)
equal(wf, f.inv() * v * reduced)
source = s.Matrix([2, -3, 1])
solution = wf * reduced.inv() * wf.T * source + hidden * (hidden.T * f * hidden).inv() * hidden.T * source
equal(solution, f.inv() * source)

wa = y * t.inv()
war = v + qh * yr * tr_plus.inv()
equal(v.T * wa, s.eye(2))
equal(v.T * war, s.eye(2))
equal(wa.T * a * wa, t.inv())

# A positive matrix can lose positivity after spatial band truncation.
positive = s.Matrix([[1, s.Rational(3, 4), s.Rational(3, 4)],
                     [s.Rational(3, 4), 1, s.Rational(3, 4)],
                     [s.Rational(3, 4), s.Rational(3, 4), 1]])
banded = positive.copy()
banded[0, 2] = banded[2, 0] = 0
assert all(positive[:i, :i].det() > 0 for i in (1, 2, 3))
assert banded.det() < 0

for label, gap, expected_q, expected_delta in (
    ("28", Fr(592, 175), Fr(1103, 2028), Fr(925, 2028)),
    ("64", Fr(11188, 3325), Fr(335297, 614997), Fr(279700, 614997)),
):
    correction = Fr(1, 800)
    b = 1 + 4 * (1 / gap + 2 * correction)
    q = 1 - 1 / b
    assert q == expected_q and 1 / b == expected_delta
    assert q**9 < Fr(1, 200)
    assert 8 * correction + q**9 < Fr(3, 200)
    print(f"{label}-link: q <= {q}, reserve >= {1/b}, degree-8 error < 3/200")

print("Exact full Schur return, positive correction, lifts and complete source: pass")
print("The collective all-volume estimate is a written proof; no radius was numerically certified.")
