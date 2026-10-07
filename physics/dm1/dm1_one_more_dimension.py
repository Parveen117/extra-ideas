"""DM1: one more dimension reads the curvature - the hypotenuse, the unseen leg, and mass as the third momentum.

Carrier and cuts as in IN1 / FR1: C1 = K, C2 = RK, C3 = iota R; reading tensor rho, readings (n; r1, r2, r3).
T1  Single reading: n^2 = r1^2 + r2^2 + r3^2.  An observer with the two real cuts sees the speed
    cos(a) = sqrt(r1^2 + r2^2)/n and misses sin(a) = r3/n;  for a real reading sin(a) = 0.
T2  r3 is the invariant of the two-cut frame: every real frame change (det 1) changes n, r1, r2 and keeps r3.
    The two-cut observer's rest mass is the reading of the cut they do not have.
T3  Mass is the third momentum:  c (xi1 C1 + xi2 C2) + c (iota k3) C3 = c (xi1 K + xi2 RK) - c k3 R,
    so the three-cut massless law on a mode of third wave number k3 is the two-cut law with flip rate
    g = -c k3; the sign of k3 is the sheet.
T4  The ladder: with one cut the unseen part is r2^2 + r3^2, with two it is r3^2, with three it is 0.
T5  One level up: a record rho = A A^dagger has n^2 - r.r = 4 |det A|^2; it is null (a single reading)
    on the doubled carrier, and its unrecoverable memory is the pairing with the second sheet.
Exact arithmetic over the Gaussian rationals, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../in1')
import in1_the_invariant as in1

c, Z, O, IOTA = in1.c, in1.Z, in1.O, in1.IOTA
mm, madd, scale, dagger, det, tensor, readings, form = (in1.mm, in1.madd, in1.scale, in1.dagger, in1.det,
                                                        in1.tensor, in1.readings, in1.form)
C1, C2, C3, ONE = in1.C1, in1.C2, in1.C3, in1.ONE
RMAT = [[Z, c(-1)], [O, Z]]                                   # the turn R on the cut-complex carrier


def cscale(z, a):
    return [[in1.cmul(z, v) for v in row] for row in a]


STATES = [(c(3), c(1)), (c(2), c(0, 1)), (c(2, 1), c(1, -3)), (c(1), c(F(3, 4), F(1))), (c(4), c(3))]


def pythagoras_control():
    rows = []
    for psi in STATES:
        n, r = readings(tensor(psi))
        if n*n != r[0]**2+r[1]**2+r[2]**2:
            raise ValueError('hypotenuse law failed')
        cos2 = (r[0]**2+r[1]**2)/(n*n)
        sin2 = (r[2]/n)**2
        if cos2+sin2 != 1:
            raise ValueError('seen and unseen legs do not complete the hypotenuse')
        real = all(x[1] == 0 for x in psi)
        if real != (sin2 == 0):
            raise ValueError('the unseen leg must vanish exactly for real readings')
        rows.append(dict(uncut=str(n), seen_speed_squared=str(cos2), unseen_squared=str(sin2), real=real,
                         ladder=dict(one_cut=str((r[1]**2+r[2]**2)/(n*n)), two_cuts=str(sin2), three_cuts='0')))
    return rows


REAL_FRAMES = [[[c(2), Z], [Z, c(F(1, 2))]], [[c(F(5, 4)), c(F(3, 4))], [c(F(3, 4)), c(F(5, 4))]],
               [[c(F(3, 5)), c(F(-4, 5))], [c(F(4, 5)), c(F(3, 5))]], [[O, c(3)], [Z, O]]]
COMPLEX_FRAME = [[c(F(3, 5), F(4, 5)), Z], [Z, c(F(3, 5), F(-4, 5))]]


def third_reading_control():
    rows = []
    for psi in STATES[1:4]:
        rho = tensor(psi)
        n0, r0 = readings(rho)
        moved = False
        for g in REAL_FRAMES:
            if det(g) != O or any(v[1] != 0 for row in g for v in row):
                raise ValueError('not a real frame change of determinant one')
            cur = mm(mm(g, rho), dagger(g))
            n, r = readings(cur)
            if r[2] != r0[2]:
                raise ValueError('the third reading changed under a two-cut frame change')
            if n*n-r[0]**2-r[1]**2 != r0[2]**2:
                raise ValueError('two-cut rest mass is not the third reading')
            moved = moved or (n, r[0], r[1]) != (n0, r0[0], r0[1])
        if not moved:
            raise ValueError('the seen readings should change')
        n, r = readings(mm(mm(COMPLEX_FRAME, rho), dagger(COMPLEX_FRAME)))
        rows.append(dict(third_reading=str(r0[2]), after_a_cut_complex_frame_change=str(r[2])))
    if all(r['third_reading'] == r['after_a_cut_complex_frame_change'] for r in rows):
        raise ValueError('a three-cut frame change should move the third reading')
    return rows


def mass_control():
    rows = []
    for speed, xi1, xi2, k3 in ((F(1), F(2), F(-1), F(3)), (F(1, 2), F(1, 3), F(5), F(-2)), (F(3), F(0), F(0), F(7, 4))):
        three = madd(madd(scale(speed*xi1, C1), scale(speed*xi2, C2)), cscale(c(0, speed*k3), C3))
        two = madd(madd(scale(speed*xi1, C1), scale(speed*xi2, C2)), scale(-speed*k3, RMAT))
        if three != two:
            raise ValueError('the third propagation term is not the mass turn')
        want = speed*speed*(xi1*xi1+xi2*xi2-k3*k3)
        if mm(two, two) != scale(want, ONE):
            raise ValueError('square law failed')
        rows.append(dict(speed=str(speed), third_wave_number=str(k3), flip_rate=str(-speed*k3),
                         sheet='upper' if k3 < 0 else 'lower'))
    # reversing the third momentum reverses the turn and nothing else
    a = madd(scale(F(1), C1), cscale(c(0, 5), C3))
    b = madd(scale(F(1), C1), cscale(c(0, -5), C3))
    if madd(a, scale(F(-1), b)) != scale(F(-10), RMAT):
        raise ValueError('sheet is not the direction along the third cut')
    return rows


def doubled_control():
    rows = []
    for A in ([[c(3), c(1)], [c(1), c(0, 2)]], [[c(1, 1), c(2)], [c(0, 1), c(1, -1)]], [[c(2), c(4)], [c(1), c(2)]],
              [[c(F(3, 5)), Z], [Z, c(F(4, 5))]]):
        rho = mm(A, dagger(A))
        n, r = readings(rho)
        dA = det(A)
        mod2 = dA[0]**2+dA[1]**2
        if form(n, r) != 4*mod2:
            raise ValueError('unrecoverable memory is not 4 |det A|^2')
        # on the doubled carrier the reading is single: its tensor has rank one, i.e. all 2x2 minors of the
        # coefficient vector vanish; here that is the statement that A itself is the reading
        rank_one = (dA == Z)
        if rank_one != (form(n, r) == 0):
            raise ValueError('a record is null exactly when it is a single reading')
        rows.append(dict(unrecoverable=str(form(n, r)), pairing_with_second_sheet_squared=str(mod2), single=rank_one))
    return rows


def run():
    return dict(pythagoras=pythagoras_control(), third_reading=third_reading_control(), mass=mass_control(),
                doubled=doubled_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('DM1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        for row in v:
            print(k, row)
