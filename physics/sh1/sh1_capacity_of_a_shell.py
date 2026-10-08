"""SH1: how many readings a shell holds: 2 n^2, from least-cost readings and independence.

MC1: least cost = no Laplacian.  QC1-T3: a reading of given degree splits into contents.  RC1-R7: the first-order bound
histories of the 1/r attraction are great circles of the 3-sphere in four dimensions.  IN1: a reading has two components.
MM1: the memory of a record is an elementary symmetric function of its tensor.
Exact linear algebra (sympy, rationals).  Python 3.12."""
import itertools, json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def monomials(d, deg):
    return [e for e in itertools.product(range(deg + 1), repeat=d) if sum(e) == deg]


def least_cost_count(d, deg):
    """dimension of the homogeneous polynomial readings of degree deg in d variables with zero Laplacian."""
    if deg < 2:
        return len(monomials(d, deg))
    src, dst = monomials(d, deg), monomials(d, deg - 2)
    idx = {m: i for i, m in enumerate(dst)}
    M = sp.zeros(len(dst), len(src))
    for j, m in enumerate(src):
        for a in range(d):
            if m[a] >= 2:
                t = list(m); t[a] -= 2
                M[idx[tuple(t)], j] += m[a]*(m[a] - 1)
    return len(src) - M.rank()


def e_k(M, k):
    n = M.shape[0]
    return sum(M.extract(list(c), list(c)).det() for c in itertools.combinations(range(n), k))


def run():
    out, rec = {}, {}
    # S1: least-cost readings of degree l: two cuts: 2 (the contents +-l of QC1-T3); three cuts: 2l + 1; four variables: (l+1)^2
    two = [least_cost_count(2, l) for l in range(1, 7)]
    three = [least_cost_count(3, l) for l in range(0, 7)]
    four = [least_cost_count(4, l) for l in range(0, 6)]
    rec['two_cuts'] = two; rec['three_cuts'] = three; rec['four_variables'] = four
    out['S1_two_cuts'] = two == [2]*6
    out['S1_three_cuts'] = three == [2*l + 1 for l in range(0, 7)]
    out['S1_four_variables'] = four == [(l + 1)**2 for l in range(0, 6)]
    # S2: shell n = degrees 0 .. n-1 in three cuts = degree n-1 on the 3-sphere of RC1-R7
    out['S2_shell'] = all(sum(three[:n]) == n*n == four[n - 1] for n in range(1, 6))
    # S3: two components per reading (IN1): 2 n^2
    caps = [2*n*n for n in range(1, 5)]
    rec['capacities'] = caps
    out['S3_capacities'] = caps == [2, 8, 18, 32]
    # S4: independence.  IN1-T5's memory of two readings is the squared wedge: n1 n2 - r1.r2 = 2 |psi1 ^ psi2|^2
    a, b, c, d = sp.symbols('a b c d')
    I = sp.I
    C = [sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]])]
    p1, p2 = sp.Matrix([a, b]), sp.Matrix([c, d])
    nn = lambda p: (p.H*p)[0]
    rr = lambda p: [(p.H*Ci*p)[0] for Ci in C]
    lhs = sp.expand(nn(p1)*nn(p2) - sum(x*y for x, y in zip(rr(p1), rr(p2))))
    wedge = a*d - b*c
    out['S4_two_readings_wedge'] = sp.simplify(lhs - 2*sp.expand(wedge*sp.conjugate(wedge))) == 0
    # N readings in a space of D: the record tensor rho = sum p_s psi_s psi_s^+ has e_N = (prod p_s) x Gram determinant: zero iff dependent
    vs = [sp.Matrix([1, 2, 0, 1]), sp.Matrix([0, 1, 1, 3]), sp.Matrix([2, 0, 1, 1]), sp.Matrix([1, 1, 1, 0])]
    ps = [sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 8), sp.Rational(1, 8)]
    ok = True
    for N in (2, 3, 4):
        rho = sum((ps[s]*vs[s]*vs[s].T for s in range(N)), sp.zeros(4))
        gram = sp.Matrix(N, N, lambda i, j: vs[i].dot(vs[j])).det()
        prod = sp.prod(ps[:N])
        ok &= (e_k(rho, N) == prod*gram) and gram != 0
    out['S4_memory_is_gram_determinant'] = ok
    dep = [vs[0], vs[1], vs[0] + 2*vs[1]]
    rho = sum((sp.Rational(1, 3)*v*v.T for v in dep), sp.zeros(4))
    out['S4_dependent_readings_leave_none'] = e_k(rho, 3) == 0
    five = vs + [sp.Matrix([3, 1, 4, 1])]
    rho5 = sum((sp.Rational(1, 5)*v*v.T for v in five), sp.zeros(4))
    out['S4_no_more_than_the_dimension'] = sp.Matrix(5, 5, lambda i, j: five[i].dot(five[j])).det() == 0 and rho5.rank() == 4
    # S5: the rows of the table of elements have lengths 2, 8, 8, 18, 18, 32 (known): every length is a capacity 2 n^2
    rows = [2, 8, 8, 18, 18, 32]
    noble = [2, 10, 18, 36, 54, 86]
    out['S5_rows_are_capacities'] = all(r in caps for r in rows) and [sum(rows[:i + 1]) for i in range(6)] == noble
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, record=rec), open(os.path.join(HERE, 'SH1_RESULT.json'), 'w'), indent=1)
    return out, rec


if __name__ == '__main__':
    print(run())
