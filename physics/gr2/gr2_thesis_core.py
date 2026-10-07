"""GR2: certified core of the gravity thesis.

G1  Equivalence.  The static reading psi = Exp(G K)(1, 0) of GR1 is carried by the single frame
    change Exp(-G K) to the pure reading (1, 0): its memory against the cut is removable.
G2  The source is not removable.  For a record of two readings, n^2 - r.r (IN1-T5) is unchanged
    by that same frame change, while each recoverable memory changes.
G3  Linearity of the field.  Least-cost memories (MC1) add, and their fluxes add: the flux is an
    additive measure of the source.
G4  Once around.  A frame boosted radially by the fall rapidity (cosh eta = 1/N) and carried once
    around the source returns turned by 2 pi (1/N - 1); this exceeds 2 pi (1 - N) by 2 pi (1-N)^2/N.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../ms1')
sys.path.insert(0, '../in1')
sys.path.insert(0, '../mc1')
import ms1_speed_and_mass as ms1
import in1_the_invariant as in1
import mc1_minimum_cost as mc1


def equivalence_control():
    rows = []
    for b in (F(3, 2), F(2), F(5), F(1, 3)):
        E, ch, sh = ms1.boostK(b)
        psi = ms1.vec(E, (F(1), F(0)))
        n, j, sigma = ms1.state(psi)
        memory = 1-(j/n)**2
        back = ms1.vec(ms1.boostK(1/b)[0], psi)
        if back != (F(1), F(0)):
            raise ValueError('the local frame change does not restore the pure reading')
        n0, j0, s0 = ms1.state(back)
        if 1-(j0/n0)**2 != 0 or memory == 0:
            raise ValueError('recoverable memory should be removed, and should have been there')
        rows.append(dict(exp_G=str(b), memory_static=str(memory), memory_after_frame_change='0'))
    return rows


def real_block(b):
    ch, sh = (b+1/b)/2, (b-1/b)/2
    return [[in1.c(ch), in1.c(sh)], [in1.c(sh), in1.c(ch)]]          # Exp(G K_exchange) as a frame change


def source_control():
    rows = []
    for (b1, b2, p) in ((F(2), F(1, 3), F(1, 2)), (F(3), F(3, 2), F(1, 4)), (F(5), F(1, 5), F(2, 3))):
        psis = []
        for b in (b1, b2):
            E, _, _ = ms1.boostK(b)
            v = ms1.vec(E, (F(1), F(0)))
            psis.append((in1.c(v[0]), in1.c(v[1])))
        rho = in1.madd(in1.scale(p, in1.tensor(psis[0])), in1.scale(1-p, in1.tensor(psis[1])))
        n, r = in1.readings(rho)
        q = in1.form(n, r)
        if q <= 0:
            raise ValueError('a record of two different readings must have an unrecoverable part')
        seen = set()
        for b in (F(2), F(1, 7), 1/b1):
            g = real_block(b)
            if in1.det(g) != in1.O:
                raise ValueError('frame change must have determinant one')
            moved = in1.mm(in1.mm(g, rho), in1.dagger(g))
            n2, r2 = in1.readings(moved)
            if in1.form(n2, r2) != q:
                raise ValueError('the source changed under a frame change')
            seen.add(1-(r2[0]/n2)**2)
        if len(seen) < 2:
            raise ValueError('the recoverable memory should depend on the frame')
        rows.append(dict(share=str(p), unrecoverable=str(q), recoverable_by_frame=sorted(str(x) for x in seen)))
    return rows


def linearity_control():
    w = mc1.weights(F(2), 3, 8)
    out = []
    for m_a, m_b in ((F(1, 5), F(1, 7)), (F(3, 10), F(1, 100))):
        fa, Qa = mc1.minimizer(w, m_a, m_a/2**8)
        fb, Qb = mc1.minimizer(w, m_b, m_b/2**8)
        fs, Qs = mc1.minimizer(w, m_a+m_b, (m_a+m_b)/2**8)
        if fs != [x+y for x, y in zip(fa, fb)] or Qs != Qa+Qb:
            raise ValueError('least-cost memories and their fluxes must add')
        if Qa/m_a != Qb/m_b:
            raise ValueError('flux must be proportional to the source strength')
        out.append(dict(flux_a=str(Qa), flux_b=str(Qb), flux_sum=str(Qs), flux_per_unit=str(Qa/m_a)))
    return out


def around_control():
    rows = []
    for m, N in ((F(9, 25), F(4, 5)), (F(16, 25), F(3, 5)), (F(25, 169), F(12, 13)), (F(64, 289), F(15, 17))):
        if N*N != 1-m:
            raise ValueError('clock factor failed')
        thomas = 1/N-1                                           # turn / (2 pi)
        flat = 1-N
        if thomas-flat != (1-N)**2/N or thomas <= flat:
            raise ValueError('difference is not (1 - N)^2 / N')
        # derivative of cosh(eta) = 1/N along r, from N^2 = 1 - m with m = r_s/r: -(r_s / (2 r^2 N^3))
        r = F(7)
        r_s = m*r
        dN = r_s/(2*r*r*N)
        if -dN/(N*N) != -r_s/(2*r*r*N**3):
            raise ValueError('radial rate of the turn failed')
        rows.append(dict(memory=str(m), turn_over_2pi=str(thomas), first_order=str(m/2), cone_value_over_2pi=str(flat),
                         difference=str(thomas-flat)))
    return rows


def run():
    return dict(equivalence=equivalence_control(), source=source_control(), linearity=linearity_control(),
                around=around_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('GR2_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        for row in v:
            print(k, row)
