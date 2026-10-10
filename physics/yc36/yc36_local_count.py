"""Direct exact controls of YC36's changed geometry and bounded algebra.

No Yang--Mills harmonic truncation, numerical box solve, or practical
filter/radius/count certificate is supplied by these finite controls.
"""

from itertools import combinations
from math import comb
import sympy as s


def psd(a):
    assert a == a.T
    for size in range(1, a.rows + 1):
        for indices in combinations(range(a.rows), size):
            assert a.extract(indices, indices).det() >= 0


def cycle_distance(x, y, length):
    return min(abs(x-y), length-abs(x-y))


def hausdorff(a, b, length):
    def directed(x, y):
        return max(min(cycle_distance(i, j, length) for j in y) for i in x)
    return max(directed(a, b), directed(b, a))


def repair(c, error, q):
    return q*(c + error*s.eye(c.rows))/(q + 2*error)


def reduction(f, retained):
    a, b = f[:retained, :retained], f[:retained, retained:]
    hidden = f[retained:, retained:]
    lift = s.eye(retained).col_join(-hidden.inv()*b.T)
    return a-b*hidden.inv()*b.T, lift


# A nearest-pair cut misses unmatched spectators and acquires volume degree.
for volume in (9, 17):
    supports = [subset for count in (1, 2)
                for subset in combinations(range(volume), count)]
    near_pair = [a for a in supports if 0 in a]
    near_hausdorff = [a for a in supports if hausdorff((0,), a, volume) <= 1]
    assert len(near_pair) == volume
    assert len(near_hausdorff) == 6
    assert len(near_hausdorff) <= sum(comb(3, count) for count in (1, 2))
    assert hausdorff((0,), (0, 4), volume) == 4
    # A far-separated input still has a bounded two-sided support ball.
    source = (0, 4)
    neighbours = [a for a in supports if hausdorff(source, a, volume) <= 1]
    union_size = len({x for x in range(volume)
                      if min(cycle_distance(x, y, volume) for y in source) <= 1})
    assert len(neighbours) <= sum(comb(union_size, count) for count in (1, 2))

# Exact local normalization, including a deliberately asymmetric assembly.
q, delta, error = s.Rational(1, 2), s.Rational(1, 2), s.Rational(1, 80)
k = s.Matrix([[s.Rational(1, 5), s.Rational(1, 20), s.Rational(1, 30)],
              [s.Rational(1, 20), s.Rational(1, 4), s.Rational(1, 40)],
              [s.Rational(1, 30), s.Rational(1, 40), s.Rational(1, 6)]])
psd(k)
psd(q*s.eye(3)-k)
raw = k + s.Matrix([[1, 3, 0], [-1, -2, 1], [0, -1, 0]])/200
assert raw != raw.T
c = (raw+raw.T)/2
psd(error*s.eye(3)-(c-k))
psd(error*s.eye(3)+(c-k))
kt = repair(c, error, q)
epsilon_local = 2*q*error/(q+2*error)
psd(kt)
psd(q*s.eye(3)-kt)
psd(epsilon_local*s.eye(3)-(kt-k))
psd(epsilon_local*s.eye(3)+(kt-k))
f, ft = s.eye(3)-k, s.eye(3)-kt
psd(ft-delta*s.eye(3))
psd(s.eye(3)-ft)

# The affine error bound is sharp, even on two commuting eigenchannels.
k_edge = s.diag(0, q)
c_edge = s.diag(error, q-error)
assert repair(c_edge, error, q)-k_edge == s.diag(epsilon_local, -epsilon_local)

# Actual reduced error and lift error obey the reserve-based comparison.
schur, lift = reduction(f, 1)
schur_t, lift_t = reduction(ft, 1)
psd(epsilon_local*s.eye(1)/delta-(schur_t-schur))
psd(epsilon_local*s.eye(1)/delta+(schur_t-schur))
dlift = lift_t-lift
psd(epsilon_local**2*s.eye(1)/delta**3-dlift.T*dlift)

# The finite hidden inverse has an explicit residual; it is not an exact lift.
p = 2
a, b, hidden = ft[:1, :1], ft[:1, 1:], ft[1:, 1:]
ch = s.eye(2)-hidden
rp = sum((ch**n for n in range(p+1)), s.zeros(2))
sp = a-b*rp*b.T
wp = s.eye(1).col_join(-rp*b.T)
assert hidden.inv()-rp == ch**(p+1)*hidden.inv()
psd(sp-schur_t)
psd(q**(p+3)*s.eye(1)/delta-(sp-schur_t))
psd(sp-delta*s.eye(1))
psd(s.eye(1)-sp)
residual = (ft*wp)[1:, :]
assert residual == ch**(p+1)*b.T and residual != s.zeros(2, 1)
assert (ft*wp)[:1, :] == sp
assert wp.T*ft*wp != sp
dw = wp-lift_t
psd(q**(2*p+4)*s.eye(1)/delta**2-dw.T*dw)
psd(q**(2*p+4)*s.eye(1)-residual.T*residual)

# Mixed-zero specialization: additive support energies, exact local blocks,
# identity affine repair, and a polynomial with zero hidden residual.
tau = s.Rational(3)
energies = s.diag(4, 5, 4+5)
k_zero = tau*(energies+tau*s.eye(3)).inv()
assert repair(k_zero, s.Rational(0), q) == k_zero
f_zero = s.eye(3)-k_zero
assert f_zero[2, 2] == s.Rational(3, 4)
assert f_zero[2, 2] != f_zero[0, 0]*f_zero[1, 1]
s_zero, w_zero = reduction(f_zero, 2)
b_zero = f_zero[:2, 2:]
assert b_zero == s.zeros(2, 1)
assert s_zero == f_zero[:2, :2]
assert w_zero == s.eye(3)[:, :2]

print("Hausdorff support degree versus unmatched nearest-pair spectators: pass")
print("Affine reserve repair, sharp error and reduced source/lift comparison: pass")
print("Finite hidden elimination, retained residual and typed zero recovery: pass")
print("No numerical Yang--Mills box, harmonic tail or practical cutoff certified.")
