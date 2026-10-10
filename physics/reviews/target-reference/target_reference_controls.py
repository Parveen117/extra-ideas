"""Exact controls for the target-reference conversation audit, not YM spectra."""

import sympy as s


def eq(a, b):
    assert all(s.simplify(x) == 0 for x in a - b)


def inner(a, b):
    return s.trace(a.T*b)


x, y, lam = s.symbols("x y lam", real=True)
eye = s.eye(2)

# The scalar centre vanishes even with nonzero mixed response and Cp/Cv > 1.
u = (x*x + 2*x*y + 4*y*y)/2
h = s.hessian(u, (x, y))
assert s.expand(u - (x*s.diff(u, x) + y*s.diff(u, y))/2) == 0
assert h.det() == 3
assert h[0, 0]*h[1, 1]/h.det() == s.Rational(4, 3)

# Eigenspace agreement does not solve the one-parameter step equations.
h = s.diag(1, 2)
target = s.diag(2, 3)
eq(h*target, target*h)
assert s.solve([1 + lam - 2, 2 + 4*lam - 3], [lam]) == []

# A separable, flat response can still change under an integrable square.
u = s.exp(x) + s.exp(y)
hfield = s.hessian(u, (x, y))
ut = u + lam*(s.exp(2*x) + s.exp(2*y))/4
eq(s.hessian(ut, (x, y)), hfield + lam*hfield**2)
assert hfield**2 != s.zeros(2)

# Target decomposition, independently in a rotated rational eigenframe.
rot = s.Matrix([[3, -4], [4, 3]])/5
eq(rot.T*rot, eye)
h = rot*s.diag(1, 4)*rot.T
target = s.Matrix([[3, 1], [1, 2]])
m, mt = s.trace(h)/2, s.trace(target)/2
xx, xt = h - m*eye, target - mt*eye
beta = inner(xx, xt)/inner(xx, xx)
alpha = mt - beta*m
perp = target - alpha*eye - beta*h
assert inner(perp, eye) == 0
assert inner(perp, h) == 0
bracket = h*target - target*h
assert inner(bracket, bracket) == (4-1)**2*inner(perp, perp)
spectral_part = rot*s.diag(*(rot.T*target*rot).diagonal())*rot.T
eq(alpha*eye + beta*h, spectral_part)

# Force information supplies a third symmetric-matrix direction when it mixes.
h = s.diag(1, 2)
g = s.Matrix([1, 1])
force = g*g.T
assert h*force - force*h != s.zeros(2)
assert s.Matrix.hstack(eye.reshape(4, 1), h.reshape(4, 1),
                      force.reshape(4, 1)).rank() == 3
aligned = s.Matrix([1, 0])
eq(h*(aligned*aligned.T), (aligned*aligned.T)*h)

# The whole chain-rule family, including its force calibration.
u = (x*x + 2*y*y)/2 + x + y
grad = s.Matrix([s.diff(u, t) for t in (x, y)])
hfield = s.hessian(u, (x, y))
eq(s.hessian(u + lam*u*u/2, (x, y)),
   (1 + lam*u)*hfield + lam*grad*grad.T)

# Least-squares direction on a declared affine response family.
directions = [h, force]
target = s.Matrix([[4, 1], [1, 3]])
residual = target - h
gram = s.Matrix([[inner(a, b) for b in directions] for a in directions])
rhs = s.Matrix([inner(a, residual) for a in directions])
coeff = gram.inv()*rhs
move = sum((c*a for c, a in zip(coeff, directions)), s.zeros(2))
assert -inner(residual, move) == -inner(move, move) < 0
assert all(inner(a, residual - move) == 0 for a in directions)

# Exact Schur elimination brings data not present in powers of A_PP.
a = s.Matrix([[2, 0, 1], [0, 3, 1], [1, 1, 4]])
assert [a[:n, :n].det() for n in (1, 2, 3)] == [2, 6, 19]
ap = a[:2, :2]
bridge = a[:2, 2:3]
effective = ap - bridge*bridge.T/4
eq(effective, a.inv()[:2, :2].inv())
eq(ap*effective - effective*ap,
   s.Matrix([[0, s.Rational(1, 4)], [-s.Rational(1, 4), 0]]))

print("Silent centre, constrained spectral reach and active flat update: pass")
print("Target projection, force channel and admitted local direction: pass")
print("Source-selected Schur return beyond retained spectral powers: pass")
print("No YM continuum conclusion or predecessor-suite replay.")
