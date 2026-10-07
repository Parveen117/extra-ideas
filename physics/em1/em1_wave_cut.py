"""EM1: the electromagnetic wave under a cut - observed, lost, and what stays.

Reading: F = E + iota B (cut-complex three-vector, units with c = 1).  The cut is conjugation:
E = (F + conj F)/2 is one sheet, iota B = (F - conj F)/2 the other.
T1  u = (E.E + B.B)/2 = conj(F).F / 2,  S = E x B,  and  u^2 - S.S = |F.F|^2 / 4,
    with F.F = (E.E - B.B) + 2 iota E.B.
T2  A single wave has F.F = 0, hence u = |S| (nothing lost).  Two waves that are not parallel
    have F.F != 0: a part that no frame can turn into flow.
T3  Frame change along x with rapidity eta acts on F as a turn through the cut-complex angle iota eta:
    F_y' = cosh F_y + iota sinh F_z,  F_z' = cosh F_z - iota sinh F_y.  It is the usual change of
    E and B, it keeps F.F, and it changes u and S.
T4  What the cut reads: E.E and B.B change with the frame; E.E - B.B and E.B do not.
T5  A wave along the frame change: amplitude x D, u x D^2, frequency x D, length of a train of a
    fixed number of crests x 1/D  =>  (energy of the train)/(frequency) is unchanged.
Exact arithmetic over the Gaussian rationals, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json

Z = (F(0), F(0))


def c(re, im=0):
    return (F(re), F(im))


def cadd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def conj(a):
    return (a[0], -a[1])


def field(E, B):
    return [(F(e), F(b)) for e, b in zip(E, B)]


def split(Fv):
    return [x[0] for x in Fv], [x[1] for x in Fv]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def cdot(a, b):
    s = Z
    for x, y in zip(a, b):
        s = cadd(s, cmul(x, y))
    return s


def energy_flow(Fv):
    E, B = split(Fv)
    u = (dot(E, E)+dot(B, B))/2
    S = cross(E, B)
    return u, S


def invariant(Fv):
    return cdot(Fv, Fv)


def identity_control():
    rows = []
    for E, B in (((1, 2, 0), (0, 1, 3)), ((F(1, 2), -1, 4), (2, F(1, 3), 1)), ((0, 1, 0), (0, 0, 1)),
                 ((3, 0, 0), (0, 0, 0)), ((1, 1, 1), (1, 1, 1))):
        Fv = field(E, B)
        u, S = energy_flow(Fv)
        ff = invariant(Fv)
        Ev, Bv = split(Fv)
        if ff != (dot(Ev, Ev)-dot(Bv, Bv), 2*dot(Ev, Bv)):
            raise ValueError('F.F is not (E.E - B.B) + 2 iota E.B')
        if u*u-dot(S, S) != (ff[0]**2+ff[1]**2)/4:
            raise ValueError('u^2 - S.S = |F.F|^2/4 failed')
        half = sum((x[0]**2+x[1]**2) for x in Fv)/2
        if half != u:
            raise ValueError('u is not conj(F).F / 2')
        rows.append(dict(u=str(u), flow_squared=str(dot(S, S)), unrecoverable_squared=str((ff[0]**2+ff[1]**2)/4),
                         single=(ff == Z)))
    if not rows[2]['single'] or rows[0]['single']:
        raise ValueError('a single wave must be null; a general field must not be')
    return rows


def boost(Fv, ch, sh):
    """Frame change along x: a turn about x through iota eta."""
    if ch*ch-sh*sh != 1:
        raise ValueError('not a boost')
    fy = cadd(cmul(c(ch), Fv[1]), cmul(c(0, sh), Fv[2]))
    fz = cadd(cmul(c(ch), Fv[2]), cmul(c(0, -sh), Fv[1]))
    return [Fv[0], fy, fz]


def frame_control():
    rows = []
    for E, B in (((1, 2, 0), (0, 1, 3)), ((0, 1, 0), (0, 0, 1)), ((F(1, 2), -1, 4), (2, F(1, 3), 1)), ((0, 1, 0), (0, 0, 0))):
        Fv = field(E, B)
        for ch, sh in ((F(5, 4), F(3, 4)), (F(13, 12), F(-5, 12)), (F(17, 8), F(15, 8))):
            G = boost(Fv, ch, sh)
            Ev, Bv = split(Fv)
            E2, B2 = split(G)
            # the usual change of the fields, v = sh/ch
            want_E = [Ev[0], ch*Ev[1]-sh*Bv[2], ch*Ev[2]+sh*Bv[1]]
            want_B = [Bv[0], ch*Bv[1]+sh*Ev[2], ch*Bv[2]-sh*Ev[1]]
            if (E2, B2) != (want_E, want_B):
                raise ValueError('the turn through iota eta is not the frame change of E and B')
            if invariant(G) != invariant(Fv):
                raise ValueError('F.F changed with the frame')
            u0, S0 = energy_flow(Fv)
            u1, S1 = energy_flow(G)
            if u1*u1-dot(S1, S1) != u0*u0-dot(S0, S0):
                raise ValueError('u^2 - S.S changed with the frame')
        Gl = boost(Fv, F(5, 4), F(3, 4))
        Ev, Bv = split(Fv)
        E2, B2 = split(Gl)
        rows.append(dict(E_squared=[str(dot(Ev, Ev)), str(dot(E2, E2))], B_squared=[str(dot(Bv, Bv)), str(dot(B2, B2))],
                         E2_minus_B2=str(dot(Ev, Ev)-dot(Bv, Bv)), E_dot_B=str(dot(Ev, Bv))))
    if rows[0]['E_squared'][0] == rows[0]['E_squared'][1]:
        raise ValueError('what the cut reads should depend on the frame')
    # a purely electric reading acquires a magnetic sheet in another frame
    if rows[3]['B_squared'] != ['0', '9/16']:
        raise ValueError('pure electric reading should show a magnetic part after the frame change')
    return rows


def two_wave_control():
    def wave(direction, amp):
        """Circular wave along +x or -x: F = amp (0, 1, -+iota)."""
        return [Z, c(amp), c(0, -amp if direction > 0 else amp)]
    a, b = wave(+1, 2), wave(+1, 3)
    same = [cadd(x, y) for x, y in zip(a, b)]
    opposite = [cadd(x, y) for x, y in zip(wave(+1, 2), wave(-1, 3))]
    if invariant(a) != Z or invariant(same) != Z:
        raise ValueError('single and parallel waves must be null')
    if invariant(opposite) == Z:
        raise ValueError('counter-running waves must carry an unrecoverable part')
    u, S = energy_flow(opposite)
    u_a, _ = energy_flow(wave(+1, 2))
    u_b, _ = energy_flow(wave(-1, 3))
    if u*u-dot(S, S) != 2*u_a*u_b*(1-(-1)):                  # 2 u1 u2 (1 - cos), opposite directions
        raise ValueError('two-wave invariant is not 2 u1 u2 (1 - cos)')
    return dict(counter_running=dict(u=str(u), flow=[str(x) for x in S],
                                     unrecoverable_squared=str(u*u-dot(S, S))),
                parallel_unrecoverable='0')


def doppler_control():
    rows = []
    for ch, sh in ((F(5, 4), F(3, 4)), (F(17, 8), F(15, 8)), (F(13, 12), F(-5, 12))):
        D = ch+sh
        Fv = [Z, c(1), c(0, -1)]                              # circular wave along x
        G = boost(Fv, ch, sh)
        if G != [Z, c(D), c(0, -D)]:
            raise ValueError('amplitude does not scale by the Doppler factor')
        u0, _ = energy_flow(Fv)
        u1, _ = energy_flow(G)
        if u1 != D*D*u0:
            raise ValueError('energy density does not scale by D^2')
        freq = D                                              # phase rate of the same wave, PR1-T4: omega' = omega (C + S)
        length = 1/D                                          # a train of a fixed number of crests
        if (u1*length)/freq != u0:
            raise ValueError('energy of the train over frequency is not invariant')
        rows.append(dict(doppler=str(D), energy_density_factor=str(D*D), train_energy_over_frequency='unchanged'))
    return rows


def run():
    return dict(identity=identity_control(), frames=frame_control(), two_waves=two_wave_control(), doppler=doppler_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('EM1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
