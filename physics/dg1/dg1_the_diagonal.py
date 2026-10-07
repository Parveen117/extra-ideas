"""DG1: the diagonal -- equal share of observed and lost -- in the gravity field.

The owner's remark: 1/2 is the diagonal observer, sin 45 = cos 45.
Sources: GR1 (N^2 + m = 1; horizon = all memory), MO1 (circular histories: energy and speed), NC1 (ratio rho),
NT-3 / QC5 (share law q(1-q)), R38 / HB1 (eight-mark clock, 1/sqrt 2), MC1.
With x = r_s / r:   observed N^2 = 1 - x,   lost m = x.
D1  Equal share m = N^2 = 1/2 is the radius r = 2 r_s; there N = 1/sqrt(2), the value of the eight-mark clock.
D2  At that radius, and only there:  (a) the speed of a circular history equals the speed of the fall frame;
    (b) a circular history has exactly the energy of rest far away (no binding).
        circular: E^2 = (1-x)^2 / (1 - 3x/2),  v^2 = (x/2)/(1-x);   fall: v^2 = x.
D3  The share law q(1-q) is largest at q = 1/2, where it is 1/4: the factor of the native curvature F = -(1/4)[X,X].
D4  rho = (1/2) d ln m / d ln r: the 1/2 in rho = -1/2 is the step from the amplitude tanh(eta) to its square m;
    the -1 is the fall-off of m in a frame of three cuts (MC1).
Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json


def energy2(x):
    return (1-x)**2/(1-F(3, 2)*x)


def orbit_speed2(x):
    return (x/2)/(1-x)


def run():
    half = F(1, 2)
    # D1
    if 1-half != half:
        raise ValueError('equal share failed')
    # D2: uniqueness on a grid of rationals below the last circular radius x < 2/3
    hits_speed, hits_energy = [], []
    for num in range(1, 200):
        x = F(num, 300)                      # 0 < x < 2/3
        if orbit_speed2(x) == x:
            hits_speed.append(x)
        if energy2(x) == 1:
            hits_energy.append(x)
    if hits_speed != [half] or hits_energy != [half]:
        raise ValueError('the diagonal radius must be the only one')
    # algebraic form: (x/2)/(1-x) = x  <=>  x (1 - 2x) = 0 ;  (1-x)^2 = 1 - 3x/2  <=>  x (2x - 1) = 0
    for x in (F(1, 7), F(2, 5), half, F(3, 5)):
        if (orbit_speed2(x)-x)*(1-x) != -x*(1-2*x)/2:
            raise ValueError('speed identity failed')
        if (energy2(x)-1)*(1-F(3, 2)*x) != x*(2*x-1)/2:
            raise ValueError('energy identity failed')
    side = 'bound' if energy2(F(1, 4)) < 1 else 'not bound'
    inner = 'unbound' if energy2(F(3, 5)) > 1 else 'bound'
    # D3
    best = max((F(k, 100)*(1-F(k, 100)), F(k, 100)) for k in range(0, 101))
    if best != (F(1, 4), half):
        raise ValueError('share law maximum failed')
    # D4: m = A r^p  =>  rho = p/2 ; p = -1 gives -1/2
    rho = {p: F(p, 2) for p in (-1, -2, 0)}
    return dict(diagonal_radius_in_r_s='2', clock_factor_squared=str(half), fall_speed_squared=str(half),
                orbit_speed_squared=str(orbit_speed2(half)), orbit_energy_squared=str(energy2(half)),
                outside=side, inside=inner, share_law_maximum=str(best[0]), rho_for_power={str(k): str(v) for k, v in rho.items()})


if __name__ == '__main__':
    res = run()
    json.dump(res, open('DG1_RESULT.json', 'w'), indent=1)
    print(res)
