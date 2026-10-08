"""ZP1: an indeterminate sum resolved -- the seen floor of every mode, and the force it leaves.

The owner's programme: what remains is either solved or shown indeterminate; indeterminacy is handled by the
Bindu-Lopa operators.  First case taken: the sum of the floors of all modes between two walls.
Sources: LC1-C4 / WQ1-W3 (every mode has the floor (1/2) kappa omega: what stays seen when nothing is lost),
PS1-T4 (sum n w^n in the two Exp sectors: scale sector 1/(4 sinh^2(x/2)) = 1/x^2 - 1/12 + ..., phase sector
-1/(4 sin^2(phi/2)) = -1/phi^2 - 1/12 - ...; the constant is sector-blind, the pole flips), EM1 (the wave),
FR1 (three cuts).
Z1  Exact: sum_{n=1}^{N} n w^n = w (1 - (N+1) w^N + N w^(N+1)) / (1 - w)^2 ; limit w/(1-w)^2 for |w| < 1.
Z2  Series (exact rational):  1/(e^y - 1) = 1/y - 1/2 + y/12 - y^3/720 + ... ;
    sum n e^(-n x) = 1/x^2 - 1/12 + x^2/240 - ... ;  phase sector: -1/phi^2 - 1/12 - phi^2/240 - ...
    The constant -1/12 is the same in both sectors; the pole and the quadratic term change sign.
Z3  One cut (a line of length L, modes omega_n = n pi c / L, floor (1/2) kappa omega_n, weight e^(-n eps/L)):
        E = pi kappa c L / (2 eps^2)  -  pi kappa c / (24 L)  + O(eps^2) .
    The indeterminate part is proportional to L.  For a wall at a inside a fixed length it does not depend on a:
    it is silent in the force.  What remains:  E = -pi kappa c / (24 L).
Z4  Three cuts (two plates a distance d apart, two transverse cuts free):
        E / area = 3 d kappa c / (pi^2 eps^4) - kappa c / (2 pi eps^3)  -  pi^2 kappa c / (720 d^3)  + O(eps^2) .
    First term proportional to d (silent in the force by the same argument), second independent of d.
    Force per area:  - pi^2 kappa c / (240 d^4) .   With kappa = hbar: 13.0 Pa at 100 nm.
Exact rational arithmetic for all coefficients; floats only for the illustration.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json
import math

ORDER = 10


def s_mul(a, b):
    out = [F(0)]*ORDER
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < ORDER:
                    out[i+j] += x*y
    return out


def s_inv(a):
    """reciprocal of a power series with a[0] != 0"""
    out = [F(0)]*ORDER
    out[0] = 1/a[0]
    for n in range(1, ORDER):
        out[n] = -sum(a[k]*out[n-k] for k in range(1, n+1))/a[0]
    return out


def exp_series(sign=1):
    return [F(sign**n, math.factorial(n)) for n in range(ORDER + 1)]


def bernoulli_generating():
    """y / (e^y - 1) as a power series"""
    e = exp_series()
    denom = e[1:ORDER+1]                      # (e^y - 1)/y
    return s_inv(denom)


def scale_sector():
    """x^2 * sum n e^(-n x) = x^2 e^-x / (1 - e^-x)^2"""
    e = exp_series(-1)
    one_minus = [-c for c in e[1:ORDER+1]]    # (1 - e^-x)/x
    return s_mul(e[:ORDER], s_inv(s_mul(one_minus, one_minus)))


def phase_sector():
    """phi^2 * ( -1 / (4 sin^2(phi/2)) )"""
    # sin(phi/2)/(phi/2) as a series in phi
    sinc = [F(0)]*ORDER
    for k in range(0, ORDER, 2):
        sinc[k] = F((-1)**(k//2), math.factorial(k+1))/F(2)**k
    sq = s_mul(sinc, sinc)                    # (2 sin(phi/2)/phi)^2 = 4 sin^2(phi/2) / phi^2
    return [-c for c in s_inv(sq)]


def run():
    # Z1
    for w in (F(1, 2), F(2, 3), F(-1, 4)):
        for N in (1, 5, 12):
            direct = sum(n*w**n for n in range(1, N+1))
            closed = w*(1 - (N+1)*w**N + N*w**(N+1))/(1-w)**2
            if direct != closed:
                raise ValueError('finite sum failed')
    # Z2
    bern = bernoulli_generating()
    if bern[:5] != [F(1), F(-1, 2), F(1, 12), F(0), F(-1, 720)]:
        raise ValueError('series of y/(e^y - 1) failed')
    sc = scale_sector()
    ph = phase_sector()
    if sc[:5] != [F(1), F(0), F(-1, 12), F(0), F(1, 240)]:
        raise ValueError('scale sector failed')
    if ph[:5] != [F(-1), F(0), F(-1, 12), F(0), F(-1, 240)]:
        raise ValueError('phase sector failed')
    if sc[2] != ph[2] or sc[0] != -ph[0] or sc[4] != -ph[4]:
        raise ValueError('constant must be sector-blind; pole and quadratic term must flip')
    # Z3: E/(pi kappa c) = (1/(2L)) [ (L/eps)^2 - 1/12 ] = L/(2 eps^2) - 1/(24 L)
    one = []
    for L, eps in ((F(3), F(1, 100)), (F(7, 2), F(1, 50))):
        pole = (1/(2*L))*sc[0]*(L/eps)**2
        const = (1/(2*L))*sc[2]
        if pole != L/(2*eps**2) or const != -1/(24*L):
            raise ValueError('one-cut energy failed')
        one.append((str(L), str(const)))
    Ltot, eps = F(10), F(1, 100)
    poles = {a: a/(2*eps**2) + (Ltot-a)/(2*eps**2) for a in (F(1), F(3), F(7, 2))}
    if len(set(poles.values())) != 1:
        raise ValueError('the indeterminate part must not depend on where the wall is')
    fin = lambda a: -1/(24*a) - 1/(24*(Ltot-a))
    a, h = F(2), F(1, 1000)
    slope = (fin(a+h) - fin(a-h))/(2*h)
    exact_slope = lambda a_: 1/(24*a_*a_) - 1/(24*(Ltot-a_)**2)
    if abs(slope - exact_slope(a)) > F(1, 10**6) or exact_slope(a) <= 0:
        raise ValueError('the wall must be pulled toward the nearer end')
    # Z4: (1/2 pi) d^2/d eps^2 [ (1/eps) / (e^(p eps) - 1) ],  p = pi/d.  Laurent terms as {power of eps: (coefficient, power of p)}
    inner = {-2: (bern[0], -1), -1: (bern[1], 0), 0: (bern[2], 1), 1: (bern[3], 2), 2: (bern[4], 3)}
    second = {e-2: (cf*e*(e-1), pw) for e, (cf, pw) in inner.items() if cf*e*(e-1) != 0}
    if second != {-4: (F(6), -1), -3: (F(-1), 0), 0: (F(-1, 360), 3)}:
        raise ValueError('three-cut expansion failed')
    # finite part: (1/(2 pi)) * (-1/360) p^3 = -pi^2 / (720 d^3)
    finite_coeff = F(1, 2)*second[0][0]                 # multiplies pi^2 / d^3
    bulk_coeff = F(1, 2)*second[-4][0]                  # multiplies d / (pi^2 eps^4)
    if finite_coeff != F(-1, 720) or bulk_coeff != F(3):
        raise ValueError('plate energy coefficients failed')
    pressure_coeff = 3*finite_coeff                     # d/dd of d^-3 gives -3 d^-4 ; force = -dE/dd
    if pressure_coeff != F(-1, 240):
        raise ValueError('plate force coefficient failed')
    hbar, cl = 1.054571817e-34, 299792458.0
    pres = lambda d: math.pi**2*hbar*cl/(240*d**4)
    return dict(series=dict(y_over_exp_minus_1=[str(x) for x in bern[:5]], scale=[str(x) for x in sc[:5]], phase=[str(x) for x in ph[:5]]),
                one_cut=one, plate_energy_coefficient=str(finite_coeff), plate_force_coefficient=str(pressure_coeff),
                pressure_Pa={'100 nm': pres(100e-9), '1 micrometre': pres(1e-6)},
                one_cut_energy_J_at_1_micrometre=-math.pi*hbar*cl/(24*1e-6))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('ZP1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
