"""WQ1: a wave is a pair of readings; its enclosed area comes in whole units; energy = count x unit x frequency.

The owner's statement: every curvature is made of waves; in discrete form it gives h nu.
Sources: QC1-T1/T2 (two readings are a canonical pair; the return is their bracket), QC1-T3 (content q is silent
iff q Theta in 2 pi Z), QC1-T4 (ladder k + n, k = q/2), QC2 (p acts as kappa d/dw), EM1 (energy of a wave train /
frequency is the same in every frame; the frame supplies the unit), DC1-K1 (one unit for frames that exchange).
W1  One wave mode = one pair (w, p) turning at rate omega, energy E = (omega/2)(w^2 + p^2).
    Its orbit encloses area A with  A / 2 pi = E / omega.
W2  The return around the orbit is the bracket: Theta = A / kappa.  A recordable mode (content 1 silent) needs
        A = 2 pi kappa n ,   i.e.   E = n kappa omega = n h nu   with  h = 2 pi kappa ,  nu = omega / 2 pi .
W3  Operator form (exact, on polynomials): raising z, lowering kappa d/dz, bracket kappa; H = kappa omega (z d/dz + 1/2)
    has H z^n = kappa omega (n + 1/2) z^n.  The 1/2 is QC1's k = q/2.
W4  Under a change of frame energy and frequency scale by the same factor (EM1): the count n does not change.
W5  A superposition of modes: the area (curvature flux) of each mode is counted separately; energies add.
Exact rational arithmetic; floats only for the illustration.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json


def action(E, omega):
    """A / 2 pi for the orbit of energy E: radius^2 = 2E/omega, area = pi radius^2."""
    radius2 = 2*E/omega
    area_over_pi = radius2
    return area_over_pi/2


def silent(q, J, kappa):
    """content q returns unchanged iff q * (A/kappa) is a whole number of turns: q J / kappa integer."""
    return (q*J/kappa).denominator == 1


def orbit_control():
    rows = []
    for omega, kappa in ((F(3), F(1, 2)), (F(7, 2), F(2, 5))):
        for n in (1, 2, 5):
            E = n*kappa*omega
            J = action(E, omega)
            if J != E/omega or J != n*kappa or not silent(1, J, kappa):
                raise ValueError('whole-count orbit failed')
        for E in (kappa*omega/2, F(7, 3)*kappa*omega):
            if silent(1, action(E, omega), kappa):
                raise ValueError('a fractional area must not be silent for content 1')
        # content 2 (the covariance) is silent already at half counts (QC1-T3)
        if not silent(2, action(kappa*omega/2, omega), kappa):
            raise ValueError('content 2 at half count failed')
        rows.append((str(omega), str(kappa)))
    return rows


# polynomials in z as coefficient lists
def raise_(p):
    return [F(0)] + p


def lower(p, kappa):
    return [kappa*k*p[k] for k in range(1, len(p))] or [F(0)]


def sub(a, b):
    n = max(len(a), len(b))
    a = a + [F(0)]*(n-len(a)); b = b + [F(0)]*(n-len(b))
    return [x-y for x, y in zip(a, b)]


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def hamiltonian(p, kappa, omega):
    """kappa omega (z d/dz + 1/2) = omega (z . kappa d/dz + kappa/2)."""
    zd = raise_(lower(p, kappa))
    n = max(len(zd), len(p))
    zd = zd + [F(0)]*(n-len(zd)); pp = p + [F(0)]*(n-len(p))
    return [omega*(a + kappa*b/2) for a, b in zip(zd, pp)]


def ladder_control():
    kappa, omega = F(2, 5), F(7, 2)
    rows = []
    for p in ([F(3), F(-1), F(2, 3), F(5)], [F(0), F(0), F(1)]):
        comm = trim(sub(lower(raise_(p), kappa), raise_(lower(p, kappa))))
        if comm != trim([kappa*x for x in p]):
            raise ValueError('bracket of lowering and raising is not kappa')
    for n in range(0, 6):
        mono = [F(0)]*n + [F(1)]
        got = trim(hamiltonian(mono, kappa, omega))
        want = trim([kappa*omega*(n + F(1, 2))*x for x in mono])
        if got != want:
            raise ValueError('ladder failed')
        rows.append(str(kappa*omega*(n + F(1, 2))))
    steps = {F(rows[i+1]) - F(rows[i]) for i in range(len(rows)-1)}
    if steps != {kappa*omega}:
        raise ValueError('steps must all equal kappa omega')
    return dict(levels=rows, step=str(kappa*omega), lowest=str(kappa*omega/2))


def frame_control():
    rows = []
    for D in (F(2), F(3, 2), F(1, 5)):           # Doppler factors of rational boosts (beta = 3/5, 5/13, -12/13)
        kappa, omega, n = F(1, 3), F(9, 2), 4
        E = n*kappa*omega
        E2, w2 = E*D, omega*D                      # EM1: train energy and frequency scale alike
        if E2/w2 != E/omega or E2/(kappa*w2) != n:
            raise ValueError('the count must not depend on the frame')
        rows.append((str(D), str(E2), str(w2)))
    return rows


def superposition_control():
    kappa = F(1, 2)
    modes = [(F(3), 2), (F(5), 0), (F(7, 2), 3)]              # (omega, count)
    total = sum(n*kappa*w for w, n in modes)
    areas = [action(n*kappa*w, w) for w, n in modes]
    if areas != [n*kappa for _, n in modes] or total != kappa*(2*3 + 0 + 3*F(7, 2)):
        raise ValueError('modes must be counted separately and energies add')
    return dict(total_energy=str(total), counts=[n for _, n in modes])


def illustration():
    h, c = 6.62607015e-34, 299792458.0
    return dict(green_light_500nm_J=h*c/500e-9, sound_1kHz_J=h*1000.0, wave_100Hz_J=h*100.0)


def run():
    return dict(orbit=orbit_control(), ladder=ladder_control(), frame=frame_control(),
                superposition=superposition_control(), illustration=illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('WQ1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
