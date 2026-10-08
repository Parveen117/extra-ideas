"""WD2: the carrier holds exactly three cuts, and their readings close the boundary in whole turns. Exact (Gaussian rationals)."""
import json, os, itertools, sympy as sp
I = sp.I; I2 = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); S = R*K
C = [K, S, I*R]                                           # the three cuts (FR1/IN1)
out = {}
anti = lambda A, B: A*B + B*A
out['T1_three_cuts'] = all(c*c == I2 for c in C) and all(anti(C[i], C[j]) == sp.zeros(2) for i in range(3) for j in range(i))
# T2: no fourth cut: anything anticommuting with all three is zero
a, b, c_, d_ = sp.symbols('a b c d')
X = sp.Matrix([[a, b], [c_, d_]])
eqs = [e for Cc in C for e in anti(X, Cc)]
out['T2_no_fourth'] = sp.solve(eqs, [a, b, c_, d_], dict=True) == [{a: 0, b: 0, c_: 0, d_: 0}]
out['T2_product_is_iota'] = sp.simplify(C[0]*C[1]*C[2] - I*I2) == sp.zeros(2)
# readings (FR1): rho(n) = (1 + n.C)/2 along +-cut directions
P = lambda i, s: (I2 + s*C[i])/2
# T3: triple product of readings round one eighth of the boundary
t = (P(0, 1)*P(1, 1)*P(2, 1)).trace()
out['T3_octant'] = sp.simplify(t - (1 + I)/4) == 0        # one-eighth turn; modulus (1/2)(1/sqrt 2)
out['T3_modulus'] = sp.simplify(sp.Abs(t) - sp.Rational(1, 2)/sp.sqrt(2)) == 0
# T4: all eight octants, each taken with outward sense: the turn-parts multiply to one full turn
prod = sp.Integer(1); phases = []
for s in itertools.product((1, -1), repeat=3):
    order = (0, 1, 2) if s[0]*s[1]*s[2] == 1 else (0, 2, 1)      # outward sense
    tr = (P(order[0], s[order[0]])*P(order[1], s[order[1]])*P(order[2], s[order[2]])).trace()
    ph = sp.simplify(tr/sp.Abs(tr)); phases.append(ph); prod *= ph
out['T4_each_one_eighth'] = all(sp.simplify(p**8 - 1) == 0 and sp.simplify(p - (1 + I)/sp.sqrt(2)) == 0 for p in phases)
out['T4_whole_boundary'] = sp.simplify(prod - 1) == 0 and sp.simplify(sp.prod(phases[:4]) + 1) == 0   # half the boundary: -1
# T5: with two cuts only (K, S) every such product is real: no turn is carried
t2 = (P(0, 1)*P(1, 1)*P(0, -1)*P(1, -1)).trace()
out['T5_two_cuts_real'] = sp.im(t2) == 0 and all(sp.im((P(0, s1)*P(1, s2)).trace()) == 0 for s1 in (1, -1) for s2 in (1, -1))
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'WD2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
