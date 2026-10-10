"""Direct YC31 convention and arithmetic checks; no lattice truncation."""

from fractions import Fraction as F

import sympy as sp


def zero(matrix):
    assert all(sp.simplify(entry) == 0 for entry in matrix)


t = sp.symbols("t", real=True)
eye = sp.eye(2)
vac = sp.Matrix([1, 0])
rot = sp.Matrix([[1 - t*t, -2*t], [2*t, 1 - t*t]]) / (1 + t*t)
ground = rot * vac
proj = ground * ground.T
ham = rot * sp.diag(0, 5) * rot.T
deriv_in_frame = sp.simplify(rot.T * sp.diff(ham, t) * rot)
energies = [0, 5]
# F_W(omega)=i/omega, and U'=iDU. Zero-frequency entries vanish.
gen_in_frame = sp.Matrix(2, 2, lambda i, j:
    0 if i == j else sp.I * deriv_in_frame[i, j] / (energies[i] - energies[j]))
gen = sp.simplify(rot * gen_in_frame * rot.T)
zero(gen - gen.conjugate().T)
zero(sp.diff(rot, t) - sp.I * gen * rot)
zero(sp.diff(proj, t) - sp.I * (gen * proj - proj * gen))
zero(rot.T * rot - eye)

# f(0)=0 and f(5)=4/9: exact centered return in its physical pairing.
tau = sp.Integer(4)
kernel = rot * sp.diag(0, tau / (5 + tau)) * rot.T
source_observable = sp.Matrix([[2, 1], [1, -1]])
excited_source = (eye - proj) * source_observable * ground
zero(kernel * excited_source - kernel * source_observable * ground)
zero(rot.T * kernel * rot - sp.diag(0, sp.Rational(4, 9)))
zero(kernel * ground)

# At H(ell)=diag(0,12)-ell V the ground derivative is +V Omega/12.
v = sp.Matrix([[0, 1], [1, 0]])
jet_generator = sp.Matrix([[0, sp.I / 12], [-sp.I / 12, 0]])
zero(sp.I * jet_generator * vac - v * vac / 12)

for label, incidence, cap, gap, expected_velocity, expected_q in (
    ("28", 48, F(1, 3200), F(592, 175), F(3, 25), F(175, 323)),
    ("64", 88, F(1, 5700), F(11188, 3325), F(176, 1425), F(3325, 6122)),
):
    assert gap > F(1, 2) > F(1, 4)
    assert 8 * incidence * cap == expected_velocity
    assert F(4) / (gap + 4) == expected_q
    print(f"{label}-link: v/e <= {expected_velocity}; q <= {expected_q}")

print("Exact flow sign, projection derivative, pairing, source return and zero jet: pass")
print("All-volume bounds are written proofs; filter constants were not numerically certified.")
