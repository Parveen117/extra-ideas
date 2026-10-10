"""Exact algebra controls for SIGNED_COMPASS_YM_HANDOFF.md; no gap certificate."""
import sympy as sp


def check():
    a, b, c, b1, b2 = sp.symbols("a b c b1 b2", real=True)
    j = sp.Matrix([[0, -1], [-1, 0]])
    quadrants = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
    images = [tuple(j * sp.Matrix(q)) for q in quadrants]
    assert images == [quadrants[2], quadrants[1], quadrants[0], quadrants[3]]
    assert j * j == sp.eye(2) and j.det() == -1

    # Generic tangent data dT=a dS+b dV; dP=-b dS-c dV.
    # These rational identities apply wherever their derivative charts exist.
    t, v, s, p = [sp.Matrix(row) for row in
                  [(0, -1), (-a, -b), (1, 0), (-b, -c)]]

    def bracket(f, g):
        return sp.det(sp.Matrix.hstack(f, g))

    def derivative(f, g, held):
        return sp.cancel(bracket(f, held) / bracket(g, held))

    residuals = [
        derivative(v, t, s) + derivative(p, s, t),
        derivative(s, t, v) - derivative(p, v, t),
        derivative(v, p, s) - derivative(t, s, p),
        derivative(s, p, v) + derivative(t, v, p),
    ]
    assert all(sp.cancel(r) == 0 for r in residuals)

    T, V, S, P = sp.symbols("T V S P")
    t0, v0, s0, p0 = -V, -T, S, P
    assert sp.expand(S*T+V*P+s0*v0+t0*p0) == 0

    response = sp.Matrix([[a, b1], [b2, c]])
    moved = j*response*j.T
    assert moved == sp.Matrix([[c, b2], [b1, a]])
    assert moved.det() == response.det()
    assert moved[0, 0]*moved[1, 1] == a*c
    assert (moved[1, 0]-moved[0, 1])/2 == -(b2-b1)/2

    # A correlated hidden block: changing a component sign alone must fail.
    d = sp.Matrix([[5, 1], [1, 4]])
    r = sp.Matrix([[1, 2], [3, -1]])
    z = sp.Rational(1, 2)
    uq = sp.diag(-1, 1)
    returned = r.T*(d-z*sp.eye(2)).inv()*r
    moved_r = uq*r*j.T
    moved_d = uq*d*uq.T
    moved_return = moved_r.T*(moved_d-z*sp.eye(2)).inv()*moved_r
    assert moved_return == j*returned*j.T
    assert (-r).T*(d-z*sp.eye(2)).inv()*(-r) == returned
    wrong_return = (uq*r).T*(d-z*sp.eye(2)).inv()*(uq*r)
    assert wrong_return != returned
    print("PASS: quadrant map; four generic Maxwell identities; centre pairing;")
    print("      response orientation; full-return covariance; incorrect-sign control.")


if __name__ == "__main__":
    check()
