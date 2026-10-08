"""DO1: the observer is on the diagonal.  For the observer whose cuts all read alike, the law is  seen = lost  (F = R - D = 0),
and for a centre that alone gives rho = -1/2.

T24: S = R + D, F = R - D with Z the block and P a cut.  The cut is not fixed: the diagonal observer is the cut frame
in which every cut reads the same.  Exact rational arithmetic and sympy.  Python 3.12."""
import itertools, json, os
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, rec = {}, {}
    # P1: one plane.  Block with rates l1, l2 read by cuts turned by chi from its own axes
    l1, l2, chi = sp.symbols('lambda1 lambda2 chi', real=True)
    Rot = sp.Matrix([[sp.cos(chi), -sp.sin(chi)], [sp.sin(chi), sp.cos(chi)]])
    X = (Rot.T*sp.diag(l1, l2)*Rot).applyfunc(sp.simplify)
    a, dl = (l1 + l2)/2, (l1 - l2)/2
    out['P1_entries'] = z(X[0, 0] - (a + dl*sp.cos(2*chi))) and z(X[1, 1] - (a - dl*sp.cos(2*chi))) and z(X[0, 1]**2 - dl**2*sp.sin(2*chi)**2)
    Xd = X.subs(chi, sp.pi/4).applyfunc(sp.simplify)
    out['P1_diagonal_reads_alike'] = z(Xd[0, 0] - a) and z(Xd[1, 1] - a) and z(Xd[0, 1]**2 - dl**2)
    # T24 on the diagonal: for each cut R = a^2, D = b^2, and F = R - D is the determinant of the block (EMK-1)
    R_, D_ = Xd[0, 0]**2, Xd[0, 1]**2
    out['P1_F_is_the_determinant'] = z(R_ - D_ - l1*l2)
    # for any other observer the single-cut F is not the determinant
    Fgen = sp.simplify(X[0, 0]**2 - X[0, 1]**2 - l1*l2)
    out['P1_only_on_the_diagonal'] = z(Fgen.subs(chi, sp.pi/4)) and not z(Fgen.subs(chi, sp.pi/6)) and not z(Fgen.subs(chi, 0))
    # maximum: the lost part D(chi) = dl^2 sin^2(2 chi) is largest on the diagonal, zero on the block's own axes
    Dchi = dl**2*sp.sin(2*chi)**2
    out['P1_lost_is_largest_on_the_diagonal'] = z(sp.diff(Dchi, chi).subs(chi, sp.pi/4)) and z(Dchi.subs(chi, sp.pi/4) - dl**2) and z(Dchi.subs(chi, 0))
    # complete: on the diagonal one cut's (seen, lost) = (mean, half-difference): both rates follow
    out['P1_complete'] = z((a + dl) - l1) and z((a - dl) - l2)
    # P2: three cuts (and d cuts) equally inclined to the fall: K = c [ rho n n^T + (1 - n n^T) ], n = (1,...,1)/sqrt(d)
    rho, c = sp.symbols('rho c', real=True)
    ok = True
    for d in (2, 3, 4, 5, 6):
        J = sp.ones(d, d)
        K = c*(sp.eye(d) + (rho - 1)*J/d)
        seen, lost = K[0, 0], K[0, 1]
        ok &= z(seen - c*(rho + d - 1)/d) and z(lost - c*(rho - 1)/d)
        sol = sp.solve(sp.Eq(seen**2, lost**2), rho)
        ok &= (sol == [-sp.Rational(d - 2, 2)])
        # then every plane determinant vanishes, the lost part of a cut is (d-1) x its seen part, and e2 = 0
        Kg = K.subs(rho, -sp.Rational(d - 2, 2))
        ok &= all(z(Kg[i, i]*Kg[j, j] - Kg[i, j]**2) for i in range(d) for j in range(i + 1, d))
        ok &= z(sum(Kg[0, j]**2 for j in range(1, d)) - (d - 1)*Kg[0, 0]**2)
    out['P2_seen_equals_lost_gives_rho'] = ok
    K3 = (sp.eye(3) + (sp.Rational(-1, 2) - 1)*sp.ones(3, 3)/3)
    out['P2_three_cuts_block'] = K3 == sp.Rational(1, 2)*sp.Matrix([[1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) and sorted(K3.eigenvals().items()) == sorted({sp.Rational(-1, 2): 1, sp.Integer(1): 2}.items())
    # P3: the converse.  All blocks with every cut reading a and every pair losing a (signs free):
    def spectra(d):
        found = {}
        pairs = list(itertools.combinations(range(d), 2))
        for signs in itertools.product((1, -1), repeat=len(pairs)):
            M = sp.eye(d)
            for (i, j), s_ in zip(pairs, signs):
                M[i, j] = M[j, i] = s_
            key = tuple(sorted(sp.nsimplify(v) for v, mult in M.eigenvals().items() for _ in range(mult)))
            prod = 1
            if d == 3:
                prod = signs[0]*signs[1]*signs[2]
            found.setdefault(str(key), set()).add(prod)
        return found
    s2, s3, s4 = spectra(2), spectra(3), spectra(4)
    rec['spectra_two_cuts'] = sorted(s2)
    rec['spectra_three_cuts'] = {k: sorted(v) for k, v in s3.items()}
    rec['spectra_four_cuts'] = sorted(s4)
    out['P3_two_cuts_single_reading_only'] = sorted(s2) == ['(0, 2)']
    out['P3_three_cuts_nothing_or_the_centre'] = s3 == {'(0, 0, 3)': {1}, '(-1, 2, 2)': {-1}}
    out['P3_four_cuts_has_more'] = len(s4) > 2
    # P4: three cuts, any block with every cut reading the same a: e2 = sum over planes of (a^2 - K_ij^2) = (seen total) - (lost total)/2
    a_, b12, b13, b23 = sp.symbols('a b12 b13 b23', real=True)
    Kd = sp.Matrix([[a_, b12, b13], [b12, a_, b23], [b13, b23, a_]])
    e2 = sum(Kd[i, i]*Kd[j, j] - Kd[i, j]**2 for i in range(3) for j in range(i + 1, 3))
    seenT = 3*a_**2; lostT = 2*(b12**2 + b13**2 + b23**2)
    out['P4_law_for_the_diagonal_observer'] = z(e2 - (seenT - lostT/2)) and z(e2 - sum(a_**2 - b**2 for b in (b12, b13, b23)))
    # such a frame exists for every block: example with rates 2, 5, -1 (mean 2): two steps, the second a diagonal in one plane
    D0 = sp.diag(5, -1, 2)
    t1 = sp.atan(sp.sqrt(sp.Rational(1, 1)))                        # first plane (5, -1): turn so that one cut reads the mean 2: cos^2 = 1/2
    G1 = sp.Matrix([[sp.cos(t1), -sp.sin(t1), 0], [sp.sin(t1), sp.cos(t1), 0], [0, 0, 1]])
    M1 = (G1.T*D0*G1).applyfunc(sp.simplify)
    t2 = sp.pi/4
    G2 = sp.Matrix([[1, 0, 0], [0, sp.cos(t2), -sp.sin(t2)], [0, sp.sin(t2), sp.cos(t2)]])
    M2 = (G2.T*M1*G2).applyfunc(sp.simplify)
    out['P4_frame_exists_example'] = all(z(M2[i, i] - 2) for i in range(3))
    # P5: a record of turns: the diagonal is the equal share of 1 and R (a quarter turn apart): seen = lost = 1/2, mean of size 1/sqrt 2
    mean = (Fr(1, 2), Fr(1, 2))                                     # (1 + R)/2 as a + bR
    Rturn = mean[0]**2 + mean[1]**2
    out['P5_turn_record_on_the_diagonal'] = (Rturn == Fr(1, 2)) and (1 - Rturn == Rturn)
    th = sp.Symbol('Theta', positive=True)
    out['P5_quarter_turn'] = sp.solveset(sp.Eq(sp.cos(th/2)**2, sp.Rational(1, 2)), th, sp.Interval(0, sp.pi)) == sp.FiniteSet(sp.pi/2)
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, record=rec), open(os.path.join(HERE, 'DO1_RESULT.json'), 'w'), indent=1)
    return out, rec


if __name__ == '__main__':
    print(run())
