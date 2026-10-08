"""LC1: the part lost to a cut is the part that is counted.

Sources: LN1-N1 (lost = (H^2)_ii - (H_ii)^2), QC1-T3 (content q silent iff q Theta in 2 pi Z), QC1-T4 (on a circle
of hyperbolic radius d: Theta = -pi (cosh d - 1); silent circles cosh d_n = 1 + 2n/q; heights k + n, k = q/2),
RMG1-T2 (f = (1/4) sinh(l/2) dl ^ dphi, d = l/2), WQ1-W3 (levels kappa omega (n + 1/2)), GR2 (once-around turn
2 pi (1/N - 1)), LT1-T3 (level below: tanh^2(eta/2) = x = r_s / 4s).
Normalised element H = Exp(d n), its root M = Exp((d/2) n); both read in a cut across the axis n.
    seen(H) = cosh d ,      lost(M) = sinh^2(d/2) =: L ,      seen(H) = 1 + 2 L .
C1  Return around a circle:  Theta = -pi (cosh d - 1) = -2 pi L.   The return, in turns, is the lost part.
C2  Connection:  A = L dphi ,  dA = dL ^ dphi = (1/4) sinh(l/2) dl ^ dphi  (RMG1's curvature).
C3  Content q is silent  <=>  q L is a whole number:  L_n = n / q.   The count is the lost part.
C4  Height of the ladder:  k cosh d_n = k (1 + 2 L_n) = k + n :  k times the seen part.  The lowest value k = q/2 is
    what stays seen when nothing is lost; for q = 1 the levels (n + 1/2) of WQ1 are (1/2) seen, n = lost.
C5  Gravity: the once-around turn of GR2 is  2 pi (1/N - 1) = 4 pi x / (1 - x) = 4 pi sinh^2(eta/2) :
    4 pi times the lost part at the lower level of the tower.
Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../qc1')
sys.path.insert(0, '../ln1')
import qc1_reading_content_closure as qc
import ln1_lost_and_returned as ln


def element(t):
    """Exp(u L) read in the K cut, with e^u = t."""
    ch, sh = (t + 1/t)/2, (t - 1/t)/2
    return [[ch, sh], [sh, ch]]


def laurent_derivative(p):
    return {e-1: c*e for e, c in p.items() if e != 0}


def laurent_eval(p, t):
    return sum(c*t**e for e, c in p.items())


def run():
    rows = []
    for t in (F(2), F(3, 2), F(5, 4), F(11, 10)):          # t = e^{d/2}
        M = element(t)
        H = ln.mul(M, M)
        if H[0][0]*H[1][1] - H[0][1]**2 != 1:
            raise ValueError('normalised element must have unit determinant')
        L = ln.lost(M)[0]
        seen = H[0][0]                                      # cosh d
        if seen != 1 + 2*L:
            raise ValueError('seen(H) = 1 + 2 lost(M) failed')
        turns = (seen - 1)/2                                # -Theta / 2 pi  (QC1-T4)
        if turns != L:
            raise ValueError('the return in turns must be the lost part')
        rows.append((str(t), str(seen), str(L)))
    # C2: with l = 4 ln t,  L = (t^2 - 2 + t^-2)/4 ;  dL/dl = (t/4) dL/dt  must equal (1/4) sinh(l/2) = (t^2 - t^-2)/8
    Lp = {2: F(1, 4), 0: F(-1, 2), -2: F(1, 4)}
    dL = laurent_derivative(Lp)
    for t in (F(2), F(7, 5), F(1, 3)):
        if (t/4)*laurent_eval(dL, t) != (t**2 - t**-2)/8:
            raise ValueError('dA is not the curvature of RMG1')
    # C3, C4
    ladder_rows = []
    for q in (1, 2, 3, 4):
        k = F(q, 2)
        theirs = qc.ladder(q, 3)
        for n in (1, 2, 3):
            L = F(n, q)
            cosh_d = 1 + 2*L
            if (q*L).denominator != 1 or k*cosh_d != k + n:
                raise ValueError('ladder from the lost part failed')
            circle = theirs['circles'][n-1]
            if circle['cosh_d'] != str(cosh_d) or circle['height'] != str(k + n):
                raise ValueError('disagreement with the certified ladder of QC1')
        for bad in (F(1, 2*q + 1), F(2*q + 1, 2*q)):
            if (q*bad).denominator == 1:
                raise ValueError('a fractional lost part must not be silent')
        ladder_rows.append((q, str(k), [str(k + n) for n in (1, 2, 3)], [str(F(n, q)) for n in (1, 2, 3)]))
    # WQ1: q = 1 levels n + 1/2 = seen/2 at lost = n
    for n in range(0, 5):
        if F(1, 2)*(1 + 2*n) != n + F(1, 2):
            raise ValueError('half the seen part failed')
    # C5
    grav = []
    for x in (F(1, 4), F(1, 9), F(1, 100), F(1, 3)):
        N = (1 - x)/(1 + x)
        turn = 1/N - 1                                      # GR2, in units of 2 pi
        lost_low = x/(1 - x)                                # sinh^2(eta/2) with tanh^2(eta/2) = x
        if turn != 2*lost_low:
            raise ValueError('once-around turn is not twice the lost part of the lower level')
        grav.append((str(x), str(N), str(turn)))
    silent_radii = []
    for n in (1, 2, 3):
        N = F(1, n + 1)                                     # 1/N - 1 = n
        r_over_rs = 1/(1 - N*N)
        silent_radii.append((n, str(r_over_rs)))
    return dict(elements=rows, ladder=ladder_rows, gravity=grav, whole_turn_radii_in_r_s=silent_radii)


if __name__ == '__main__':
    res = run()
    json.dump(res, open('LC1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
