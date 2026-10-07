"""PR2: mass and antimass as the two sheets; what an accelerated observer can and cannot read.

EMK elements X = x H + y K + t R (H, K split, R circular; R28 roles) have X^2 = x^2 + y^2 - t^2.
Frame changes act by similarity with boosts Exp(eta A / 2) and turns Exp(phi R / 2).
  T1  the form X^2 and, for timelike X, the sign of t are invariant: two sheets.
  T2  the boosted mass generator g R is (energy, momentum): g (cosh eta R + sinh eta A').
  T3  antimass = the reversed turn Exp(-theta X) = the other sheet.
  T4  the pair record {Exp(theta X), Exp(-theta X)} with weights (p, 1-p):
        mean = cos(theta) + (2p-1) sin(theta) X,   memory = 4 p (1-p) sin^2(theta);
      the even part and the memory are scalars (same in every frame, along every
      acceleration history); the odd part is a frame-dependent vector, zero iff p = 1/2.
  T5  the generator alpha A + g R has square alpha^2 - g^2 (three sectors).
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../pr1')
import pr1_sector_speed_and_mass as pr1

ONE, H, K, R, mm, add, sc = pr1.ONE, pr1.H, pr1.K, pr1.R, pr1.mm, pr1.add, pr1.sc
ZERO = sc(F(0), ONE)


def element(x, y, t):
    return add(add(sc(x, H), sc(y, K)), sc(t, R))


def coords(X):
    """(scalar, x, y, t) of a 2x2 element in the basis 1, H, K, R."""
    s = (X[0][0]+X[1][1])/2
    x = (X[0][0]-X[1][1])/2
    y = (X[0][1]+X[1][0])/2
    t = (X[1][0]-X[0][1])/2
    if add(add(sc(s, ONE), element(x, y, t)), X, -1) != ZERO:
        raise ValueError('decomposition failed')
    return s, x, y, t


def boost(axis, ch, sh):
    """Exp(eta A / 2) and its inverse, A = a H + b K a unit split element."""
    a, b = axis
    if a*a+b*b != 1 or ch*ch-sh*sh != 1 or ch <= 0:
        raise ValueError('not a boost')
    A = add(sc(a, H), sc(b, K))
    return add(sc(ch, ONE), sc(sh, A)), add(sc(ch, ONE), sc(-sh, A))


def turn(c, s):
    if c*c+s*s != 1:
        raise ValueError('not a turn')
    return add(sc(c, ONE), sc(s, R)), add(sc(c, ONE), sc(-s, R))


def act(g, X):
    return mm(mm(g[0], X), g[1])


FRAMES = [boost((F(1), F(0)), F(5, 4), F(3, 4)), boost((F(3, 5), F(4, 5)), F(13, 12), F(5, 12)),
          boost((F(0), F(1)), F(17, 8), F(15, 8)), turn(F(4, 5), F(3, 5)),
          boost((F(-5, 13), F(12, 13)), F(5, 3), F(-4, 3))]


def form_control():
    rows = []
    for x, y, t in ((F(1, 2), F(1, 3), F(1)), (F(2), F(-1), F(-3)), (F(0), F(0), F(1)), (F(3), F(1), F(1))):
        X = element(x, y, t)
        q = x*x+y*y-t*t
        if mm(X, X) != sc(q, ONE):
            raise ValueError('square law failed')
        cur = X
        for g in FRAMES:
            cur = act(g, cur)
            s, xx, yy, tt = coords(cur)
            if s != 0 or xx*xx+yy*yy-tt*tt != q:
                raise ValueError('form is not invariant')
            if q < 0 and (tt > 0) != (t > 0):
                raise ValueError('a frame change moved a timelike element to the other sheet')
        rows.append(dict(start=[str(x), str(y), str(t)], form=str(q),
                         kind='timelike' if q < 0 else 'spacelike' if q > 0 else 'null',
                         end=[str(v) for v in coords(cur)[1:]]))
    return rows


def energy_momentum_control(g):
    rows = []
    for ch, sh in ((F(5, 4), F(3, 4)), (F(13, 12), F(5, 12)), (F(17, 8), F(15, 8))):
        B = boost((F(1), F(0)), ch, sh)
        C, S = ch*ch+sh*sh, 2*ch*sh                         # cosh eta, sinh eta
        s, x, y, t = coords(act(B, sc(g, R)))
        energy, momentum2 = t, x*x+y*y
        if energy != g*C or momentum2 != g*g*S*S or s != 0:
            raise ValueError('boosted mass generator is not (g cosh, g sinh)')
        if energy*energy-momentum2 != g*g:
            raise ValueError('energy^2 - momentum^2 = mass^2 failed')
        # the reversed turn sits on the other sheet with the same form
        s2, x2, y2, t2 = coords(act(B, sc(-g, R)))
        if t2 != -energy or x2*x2+y2*y2 != momentum2:
            raise ValueError('antimass is not the mirror sheet')
        rows.append(dict(cosh_eta=str(C), energy=str(energy), momentum_squared=str(momentum2),
                         velocity_squared=str(momentum2/(energy*energy))))
    return rows


def pair_control(c, s, p):
    """Record {Exp(theta R), Exp(-theta R)} with weights (p, 1-p), carried through all frames."""
    if c*c+s*s != 1 or not 0 <= p <= 1:
        raise ValueError('bad record')
    U, V = add(sc(c, ONE), sc(s, R)), add(sc(c, ONE), sc(-s, R))
    hist = []
    odd_seen = set()
    for n in range(len(FRAMES)+1):
        Un, Vn = U, V
        for g in FRAMES[:n]:
            Un, Vn = act(g, Un), act(g, Vn)
        mean = add(sc(p, Un), sc(1-p, Vn))
        sm, x, y, t = coords(mean)
        if sm != c:
            raise ValueError('even part is not the frame-independent scalar cos(theta)')
        odd = element(x, y, t)
        if mm(odd, odd) != sc(-((2*p-1)*s)**2, ONE):
            raise ValueError('odd part does not have the invariant square')
        # memory = variance of the two records about their mean: p |U - S|^2 + (1-p) |V - S|^2
        dU, dV = add(Un, mean, -1), add(Vn, mean, -1)
        var = add(sc(p, mm(dU, dU)), sc(1-p, mm(dV, dV)))
        if var != sc(-4*p*(1-p)*s*s, ONE):
            raise ValueError('memory is not the scalar 4 p (1-p) sin^2')
        if mm(Un, Vn) != ONE:
            raise ValueError('the two records are not reverse to each other in this frame')
        odd_seen.add((x, y, t))
        hist.append([str(x), str(y), str(t)])
    balanced = (p == F(1, 2))
    if balanced != (odd_seen == {(F(0), F(0), F(0))}) and s != 0:
        raise ValueError('odd part must vanish in every frame exactly for the balanced pair')
    if not balanced and s != 0 and len(odd_seen) == 1:
        raise ValueError('an unbalanced record should be frame dependent')
    return dict(weight=str(p), even_part=str(c), memory=str(4*p*(1-p)*s*s),
                odd_coefficient=str((2*p-1)*s), odd_part_by_frame=hist)


def sector_control(alpha, g):
    G = add(sc(alpha, H), sc(g, R))
    q = alpha*alpha-g*g
    if mm(G, G) != sc(q, ONE):
        raise ValueError('square law failed')
    return dict(split_rate=str(alpha), turn_rate=str(g), square=str(q),
                sector='circular (turns)' if q < 0 else 'dual (stops)' if q == 0 else 'split (runs away)')


def run():
    return dict(form=form_control(), energy_momentum=energy_momentum_control(F(3, 2)),
                pair=[pair_control(F(3, 5), F(4, 5), p) for p in (F(1), F(3, 4), F(1, 2), F(0))],
                sectors=[sector_control(a, F(1)) for a in (F(1, 2), F(1), F(2))])


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PR2_RESULT.json', 'w'), indent=1)
    for r in res['form']:
        print(r)
    for r in res['energy_momentum']:
        print(r)
    for r in res['pair']:
        print({k: v for k, v in r.items() if k != 'odd_part_by_frame'}, r['odd_part_by_frame'][:3])
    print(res['sectors'])
