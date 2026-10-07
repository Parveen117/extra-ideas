"""GB1: where gravity sits in the block, and where it does not.

Question (OB1): is gravity the law of the scalar / volume part of the block, as light is of the
traceless part?
B1  Every invertible block is (scale) x (phase turn) x (unit block), det = scale^2 * phase^2.
    On a reading rho -> M rho M^dagger:  the phase does nothing;  the scale multiplies n and r and
    changes n^2 - r.r;  the unit block changes n, r and keeps n^2 - r.r.
B2  A scalar-only gravity bends light too little.  Clock factor only (flat space):
    u'' + u = r_s E^2 / (2 L^2 (1 - r_s u)^2), first-order deflection r_s / b.
    A common factor on time and space bends it not at all.  The form of MO1 gives 2 r_s / b.
B3  The field of MO1 as a block: H = (1 + beta rhat.C)/N is a self-dagger unit block; it carries the
    resting unit reading to (1/N ; beta rhat / N) - the density and proper density of GR1.
B4  Light acts on a reading through the algebra (F rho + rho F^dagger, times the content q);
    the field of B3 acts through the group (g rho g^dagger), the same on every reading.
Exact arithmetic over the Gaussian rationals, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../in1')
import in1_the_invariant as in1

c, Z, O = in1.c, in1.Z, in1.O
mm, madd, scale, dagger, det, tensor, readings, form = (in1.mm, in1.madd, in1.scale, in1.dagger, in1.det,
                                                        in1.tensor, in1.readings, in1.form)
C1, C2, C3, ONE = in1.C1, in1.C2, in1.C3, in1.ONE
CUTS = [C1, C2, C3]


def cscale(z, a):
    return [[in1.cmul(z, v) for v in row] for row in a]


def vec(v):
    out = [[Z, Z], [Z, Z]]
    for x, Cm in zip(v, CUTS):
        out = madd(out, scale(x, Cm))
    return out


def factor_control():
    rho = madd(scale(F(1, 3), tensor((c(3), c(1, 2)))), scale(F(2, 3), tensor((c(1, -1), c(2)))))
    n0, r0 = readings(rho)
    q0 = form(n0, r0)
    rows = []
    unit = [[c(F(5, 4)), c(0, F(3, 4))], [c(0, F(-3, 4)), c(F(5, 4))]]       # a boost along C3
    if det(unit) != O:
        raise ValueError('not a unit block')
    for name, M in (('phase turn', cscale(c(F(3, 5), F(4, 5)), ONE)), ('scale', scale(F(3), ONE)), ('unit block', unit)):
        moved = mm(mm(M, rho), dagger(M))
        n, r = readings(moved)
        rows.append(dict(part=name, uncut=[str(n0), str(n)], invariant=[str(q0), str(form(n, r))],
                         reading_changed=(n, r) != (n0, r0)))
    if rows[0]['reading_changed']:
        raise ValueError('a reading must be blind to the phase turn')
    if rows[1]['invariant'][1] != str(q0*81) or rows[1]['uncut'][1] != str(n0*9):
        raise ValueError('the scale must multiply n by 9 and the invariant by 81')
    if rows[2]['invariant'][0] != rows[2]['invariant'][1] or not rows[2]['reading_changed']:
        raise ValueError('a unit block must move the reading and keep the invariant')
    # a general block: det = scale^2 * phase^2
    M = cscale(in1.cmul(c(3), c(F(3, 5), F(4, 5))), unit)
    d = det(M)
    want = in1.cmul(c(9), in1.cmul(c(F(3, 5), F(4, 5)), c(F(3, 5), F(4, 5))))
    if d != want:
        raise ValueError('determinant is not scale^2 phase^2')
    return rows


def light_control():
    rows = []
    for r_s, EoverL2, u in ((F(1), F(1, 100), F(1, 20)), (F(2), F(1, 9), F(1, 50)), (F(1, 2), F(4), F(1, 3))):
        # clock factor only: (du/dphi)^2 = (E^2/L^2)/(1 - r_s u) - u^2
        dF = EoverL2*r_s/(1-r_s*u)**2-2*u
        if dF/2 != r_s*EoverL2/(2*(1-r_s*u)**2)-u:
            raise ValueError('clock-only orbit equation failed')
        # full form (MO1): (du/dphi)^2 = E^2/L^2 - (1 - r_s u) u^2 for light
        dG = -(-r_s*u*u+2*u*(1-r_s*u))
        if dG/2 != -u+F(3, 2)*r_s*u*u:
            raise ValueError('full-form light equation failed')
        rows.append(dict(r_s=str(r_s)))
    # first order: clock only u1 = r_s/(2 b^2) (constant) -> total r_s/b; full form -> 2 r_s/b (MO1-M6)
    r_s, b = F(1), F(1000)
    clock_only = 2*(r_s/(2*b*b))*b                           # delta each side = b * u1, two sides
    full = 2*r_s/b
    if clock_only != r_s/b or full != 2*clock_only:
        raise ValueError('deflection values failed')
    # a common factor on time and space: which displacements are light-like does not depend on it
    for dt, dx, dy in ((F(5), F(3), F(4)), (F(13), F(5), F(12)), (F(2), F(1), F(1)), (F(1), F(1), F(1))):
        base = dt*dt-dx*dx-dy*dy
        for phi2 in (F(1, 4), F(9)):
            if (phi2*base == 0) != (base == 0) or (phi2*base > 0) != (base > 0):
                raise ValueError('a common factor changed which displacements are light-like')
    arcsec = 1.7512432813682448
    return dict(points=len(rows), clock_factor_only=str(clock_only), common_factor='0', form_of_MO1=str(full),
                at_the_sun_arcsec=dict(clock_only=arcsec/2, common_factor=0.0, form_of_MO1=arcsec))


def field_block_control():
    rows = []
    for beta, N, rhat in ((F(3, 5), F(4, 5), (F(1), F(0), F(0))), (F(4, 5), F(3, 5), (F(2, 3), F(-1, 3), F(2, 3))),
                          (F(5, 13), F(12, 13), (F(0), F(3, 5), F(4, 5)))):
        if N*N != 1-beta*beta or sum(x*x for x in rhat) != 1:
            raise ValueError('bad field point')
        H = scale(1/N, madd(ONE, vec([beta*x for x in rhat])))
        if dagger(H) != H or det(H) != O:
            raise ValueError('field block is not a self-dagger unit block')
        n, r = readings(scale(F(1, 2), H))                    # image of the resting unit reading rho = 1/2
        if n != 1/N or r != [beta*x/N for x in rhat] or form(n, r) != 1:
            raise ValueError('resting reading is not carried to (1/N ; beta rhat / N)')
        rows.append(dict(memory=str(beta*beta), density=str(n), proper_density=str(beta/N)))
    return rows


def action_control():
    """Algebra action (light) scales with the content q and flips on the other sheet; group action does not."""
    rho = scale(F(1, 2), madd(scale(F(5), ONE), vec([F(1), F(-2), F(2)])))
    Fb = madd(vec([F(1), F(2), F(0)]), cscale(c(0, 1), vec([F(0), F(1), F(3)])))
    base = madd(mm(Fb, rho), mm(rho, dagger(Fb)))
    out = {}
    for q in (1, -1, 2):
        change = scale(F(q), base)
        out[str(q)] = [str(x) for x in ([readings_raw(change)[0]]+readings_raw(change)[1])]
    if out['1'] == out['-1'] or [F(x) for x in out['2']] != [2*F(x) for x in out['1']]:
        raise ValueError('light must act in proportion to the content')
    g = [[c(F(5, 4)), c(F(3, 4))], [c(F(3, 4)), c(F(5, 4))]]
    moved = mm(mm(g, rho), dagger(g))
    n, r = readings(moved)
    n0, r0 = readings(rho)
    if form(n, r) != form(n0, r0):
        raise ValueError('group action must keep the invariant')
    return dict(light_change_by_content=out, group_action='the same for every reading; invariant kept')


def readings_raw(a):
    n = in1.trace(a)
    r = [in1.trace(mm(a, Cm)) for Cm in CUTS]
    return n[0], [x[0] for x in r]


def run():
    return dict(factors=factor_control(), light=light_control(), field_block=field_block_control(),
                action=action_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('GB1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
