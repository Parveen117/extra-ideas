"""Direct exact controls of YC35's new count and Schur identities.

These finite matrices do not certify a Yang--Mills harmonic truncation,
the infinite filter integral, or a practical excitation-count cutoff.
"""

from itertools import combinations
import sympy as s


def psd(a):
    """All principal minors, in exact rational arithmetic."""
    assert a == a.T
    for size in range(1, a.rows + 1):
        for indices in combinations(range(a.rows), size):
            assert a.extract(indices, indices).det() >= 0


def vacuum_expectation(a, factor):
    """Two qubits; keep one factor, put its spectator in vacuum."""
    bit = 1 << factor
    local = a.extract([0, bit], [0, bit])
    embedded = s.zeros(4)
    for row in range(4):
        for col in range(4):
            if row & ~bit == col & ~bit:
                embedded[row, col] = local[int(bool(row & bit)),
                                          int(bool(col & bit))]
    return local, embedded


def reduced(f, retained):
    inclusion = s.eye(f.rows)[:, :retained]
    schur = (inclusion.T*f.inv()*inclusion).inv()
    lift = f.inv()*inclusion*schur
    return schur, lift


# Positive vacuum-fixed observables stay locally vacuum-annihilating.
u = s.eye(4)
u[1, 1] = u[3, 3] = s.Rational(3, 5)
u[1, 3], u[3, 1] = -s.Rational(4, 5), s.Rational(4, 5)
assert u.T*u == s.eye(4) and u[:, 0] == s.eye(4)[:, 0]
q = [s.diag(0, 1, 0, 1), s.diag(0, 0, 1, 1)]
number = q[0] + q[1]
for projection in q:
    observed = u*projection*u.T
    assert observed*observed == observed
    for factor in (0, 1):
        local, embedded = vacuum_expectation(observed, factor)
        psd(local)
        psd(s.eye(2) - local)
        assert local[:, 0] == s.zeros(2, 1)
        assert embedded[:, 0] == s.zeros(4, 1)
        psd(q[factor] - embedded)
        shell = observed - embedded
        assert shell[:, 0] == s.zeros(4, 1)
psd(2*number - u.T*number*u)

# A fixed energy gap alone does not control reference excitation count.
large_count = s.symbols("L", positive=True)
gap_example = s.Matrix([[s.Rational(3, 2), s.Rational(1, 2)],
                        [s.Rational(1, 2), s.Rational(3, 2)]])
assert set(gap_example.eigenvals()) == {1, 2}
response = gap_example.inv()[:, 0]
assert (response.T*s.diag(1, large_count)*response)[0] == (9 + large_count)/16

# Four excited states with count 1,1,2,3. The k=1 carrier has two states.
number = s.diag(1, 1, 2, 3)
a = s.diag(2, 3, 4, 5) + s.ones(4)
tau, delta = s.Rational(3), s.Rational(2, 5)
psd(a - 2*s.eye(4))
f = a*(a + tau*s.eye(4)).inv()
psd(f - delta*s.eye(4))
psd(s.eye(4) - f)
source = s.eye(4)[:, :2]
response = a.inv()*source
schur, lift = reduced(f, 2)
assert schur == (s.eye(2) + tau*source.T*response).inv()
assert f*lift == source*schur
assert source.T*lift == s.eye(2)
psd(lift.T*lift - s.eye(2))
psd(s.eye(2)/delta - lift.T*lift)

# Count j=2 retains the first three states; it cuts the full bounded f.
inclusion_j = s.eye(4)[:, :3]
p_j = inclusion_j*inclusion_j.T
f_j = inclusion_j.T*f*inclusion_j
schur_j, local_lift = reduced(f_j, 2)
lift_j = inclusion_j*local_lift
difference = lift_j - lift
tail = (s.eye(4) - p_j)*lift
error = schur_j - schur
assert error == difference.T*f*difference
psd(error)
psd(tail.T*f*tail - error)
psd(error - delta*difference.T*difference)
psd(schur - delta*s.eye(2))
psd(s.eye(2) - schur_j)
psd(f[:2, :2] - schur_j)  # The smaller j=1 carrier has greater Schur return.
psd(s.eye(2)/delta - lift_j.T*lift_j)
assert (p_j*lift).T*f*(p_j*lift) == schur + tail.T*f*tail

# A rigorous toy B^2 from the weighted map's squared Frobenius norm.
# This is not YC35's Yang--Mills filter constant.
b_squared = s.trace(number.inv()*a.inv()*number*a.inv())
psd(b_squared*number - a.inv()*number*a.inv())
eta_squared = b_squared/3  # k/(j+1), with k=1 and j=2.
response_tail = (s.eye(4) - p_j)*response
psd(eta_squared*s.eye(2) - response_tail.T*response_tail)
psd(tau**2*eta_squared*s.eye(2) - error)
psd(tau**2*eta_squared*s.eye(2)/delta - difference.T*difference)

# Compressing first and normalizing afterwards is a different operator.
a_j = inclusion_j.T*a*inclusion_j
normalize_after_cut = a_j*(a_j + tau*s.eye(3)).inv()
assert normalize_after_cut != f_j
assert reduced(normalize_after_cut, 2)[0] != schur_j

# Exact zero/mixing-free specialization preserves each number support.
a_zero = s.diag(2, 3, 4, 5)
f_zero = a_zero*(a_zero + tau*s.eye(4)).inv()
s_zero, w_zero = reduced(f_zero, 2)
s_zero_j, _ = reduced(f_zero[:3, :3], 2)
assert s_zero_j == s_zero
assert w_zero == source
assert (s.eye(4) - p_j)*a_zero.inv()*source == s.zeros(4, 2)

print("Positive localization and the gap-alone count counterexample: pass")
print("Stationary Schur lift, monotonicity, reserve and square-tail bound: pass")
print("Compression/normalization distinction and exact zero compatibility: pass")
print("No Yang--Mills harmonic truncation or numerical cutoff was certified.")
