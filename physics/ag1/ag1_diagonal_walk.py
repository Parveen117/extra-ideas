"""AG1: the Riemann line's seam on the Yang-Mills record, then the diagonal walk.

The owner's instruction: carry the Riemann line's results to Yang-Mills -- the same along the seam on the
real and the imaginary axis, and after that the diagonal walk.

What is done.  The record is a pair of turns with the heat weight.  Its contents are the whole cut-complex
numbers alpha = a + b iota (PM1's carrier); the weight of a content is q^N(alpha), q = Exp(-pi t), t the heat
time of the cell.  Three sums over the lattice (S = source, F = flip, L = the half-shifted, lost, sum):
        S(t) = sum q^(a^2+b^2) ,   F(t) = sum (-1)^(a+b) q^(a^2+b^2) ,   L(t) = sum q^((a+1/2)^2+(b+1/2)^2) .
The flip is theorum/24's F = R - D for the cut "content on the diagonal sublattice (1 + iota) Z[iota]" (even a+b)
against the rest: the centre mark of the record.

Sources read before building: tools DS1 (F = Z* j Z; the diagonal vector 1 + iota is the prime 2; mirror form;
no crossing), PM1 (the lattice, its seam law at k = 0), RW1, Publications LAM-1 T3 / LAM-2 / LAM-3 (theta seam by
enclosures; the character mod 4; the wall of cancellation crossed by exact arithmetic), YM75 (the Riemann line's
Schur and Birman-Schwinger tools already carried to the actual Yang-Mills gap -- not repeated here), physics
DU1 (doubling the cell; the record that is its own dual), CT1, CZ1 (the weight of a twist; weak-end exponent
pi^2/2F), MG1 (the weak-coupling wall), RKF theorum/24, theorum/28.

A1  THE SEAM, THE SAME ON BOTH AXES.  The mirror t -> 1/t (contents <-> windings; the Riemann line's theta seam)
        S(1/t) = t S(t) ,     F(1/t) = t L(t) ,     L(1/t) = t F(t) .
    At t = 1 the lattice of contents and the lattice of windings are the same square lattice.
A2  A RIGHT TRIANGLE.   S^2 = F^2 + L^2  (exact, as series).  The mirror exchanges the two legs; on the seam they
    are equal:  F = L,  (F/S)^2 = 1/2 -- the eight-mark value, an angle of 45 degrees.
A3  THE DIAGONAL STEP.  Multiplying the lattice by the diagonal vector 1 + iota doubles every norm: the seen part
    of the record is the record of the doubled cell,   R(q) = (S + F)/2 = S(q^2) ,  and the count of contents
    of norm 2n equals that of norm n.  The mirror J (exchange of a and b) cancels the odd coset in the mixed
    flip, which leaves  F(q^2)^2 = S(q) F(q).
A4  THE WALK IS THE ARITHMETIC-GEOMETRIC MEAN.   S' = (S + F)/2 ,   F' = sqrt(S F)   for t -> 2t.
    So  F < F' < S' < S :  the flip only grows, the source only falls, they never cross, and the
    arithmetic-geometric mean of (S, F) is the same at every cell size: it is 1, the count of the origin.
    S'^2 - F'^2 = ((S - F)/2)^2 :  the doubled cell's source is the hypotenuse over its flip and the lost part.
A5  THE WEAK END IS THE DIAGONAL, AND THE MIRROR CROSSES THE WALL.  As t -> 0,  F/S -> 0 (seen = lost).
        F(t) = (4/t) Exp(-pi/(2t)) (1 + Exp(-2 pi/t) + ...)^2 :   smaller than every power of t.
    In the description by contents this number is what is left of a sum of terms of size 1 (39 digits cancel at
    t = 1/60); in the mirror description it is one positive term.  Both are enclosed and agree.
A6  THE TURN BLOCK (SU(2) contents m = 2j + 1, weight m^2, heat weight; u = t/4pi).  Closed surface with heat
    weight: flip = K_-(u)/K_+(u), K_+ = sum m^2 Exp(-pi u m^2), K_- = sum (-1)^(m-1) m^2 Exp(-pi u m^2).
    Mirror description:  K_+ = (1/(4 pi u^(3/2))) sum (1 - 2 pi n^2/u) Exp(-pi n^2/u) ,
                         K_- = (1/(4 pi u^(3/2))) sum (2 pi (n+1/2)^2/u - 1) Exp(-pi (n+1/2)^2/u) .
    Weak end: flip = 2 (pi/(2u) - 1) Exp(-pi/(4u)) (1 + ...): exponent pi^2/T in the heat time T of the whole
    surface; with T = 2F/kappa this is CZ1's weak-end exponent pi^2 kappa/(2F).
A8  THE TURN BLOCK WALKS AT LEAST AS FAST.  With rho = flip of the closed surface at heat time T and rho' at 2T,
        rho' (1 + rho) >= 2 sqrt(rho)          (the pair record has equality: A4)
    at every coupling tried: enclosures at nine couplings from u = 1/100 to u = 2, strict.  At the strong end
    rho' - 2 sqrt(rho)/(1 + rho) = 64 q^12 + ...  (q = Exp(-pi u));  at the weak end the ratio grows without bound.
    Not proved for all couplings.
A7  THE REFLECTION FORM IS THE MIRROR FORM.  For a transfer record c(m) = sum w lambda^m, the form
    sum f_s f_t-dagger c(s+t) has inertia (real values + pairs, pairs) (DS1-D4 with the mirror dagger): it is
    non-negative exactly when every transfer value is on the rad axis.  The gap is the same count with the line
    moved:  n_-(r c(s+t) - c(s+t+1)) = number of values above r.

Exact integer series; directed rational enclosures (tools/pm1).  Python 3.12, standard library only.
"""
from fractions import Fraction as F
from math import isqrt
import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(folder, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, '..', '..', 'tools', folder, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pm1 = load('pm1', 'pm1_prime_turn_series')
ds1 = load('ds1', 'ds1_one_diagonal')
down, up, BITS, SCALE = pm1.down, pm1.up, pm1.BITS, pm1.SCALE
TINY = F(1, 1 << (BITS - 20))


# ---------------------------------------------------------------- intervals (pairs of rationals, directed rounding)
def iv(x, y=None):
    return (F(x), F(x if y is None else y))


def iadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def isub(x, y):
    return (x[0] - y[1], x[1] - y[0])


def imul(x, y):
    ps = [x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1]]
    return (down(min(ps)), up(max(ps)))


def iscale(x, c):
    c = F(c)
    return (x[0]*c, x[1]*c) if c >= 0 else (x[1]*c, x[0]*c)


def idiv(x, y):
    assert y[0] > 0
    return imul(x, (down(1/y[1]), up(1/y[0])))


def isqrt_iv(x):
    assert x[0] >= 0
    lo = F(isqrt((x[0].numerator*SCALE*SCALE)//x[0].denominator), SCALE)
    hi = F(isqrt(-((-x[1].numerator*SCALE*SCALE)//x[1].denominator)) + 1, SCALE)
    return (lo, hi)


def overlap(x, y):
    assert x[0] <= x[1] and y[0] <= y[1], 'not an enclosure'
    return max(x[0], y[0]) <= min(x[1], y[1])


def width(x):
    return x[1] - x[0]


def below(x, y):
    return x[1] < y[0]


def powers_of_squares(e, shift=0):
    """yield enclosures of e^((n + shift/2)^2 * 4/4): for shift = 0 the numbers e^(n^2), n = 1, 2, ...;
    for shift = 1, with e an enclosure of the fourth root, the numbers e^((2n+1)^2), n = 0, 1, ..."""
    if shift == 0:
        cur, step, e2 = e, imul(e, imul(e, e)), imul(e, e)        # e^1, then multiply by e^3, e^5, ...
        while True:
            yield cur
            cur = imul(cur, step)
            step = imul(step, e2)
    else:
        e8 = imul(imul(imul(e, e), imul(e, e)), imul(imul(e, e), imul(e, e)))
        cur, step = e, e8                                          # e^1, then multiply by e^8, e^16, e^24, ...
        while True:
            yield cur
            cur = imul(cur, step)
            step = imul(step, e8)


def theta3(e):
    """1 + 2 sum_(n>=1) e^(n^2)"""
    total = iv(1)
    for term in powers_of_squares(e):
        total = iadd(total, iscale(term, 2))
        if term[1] < TINY:
            tail = 2*term[1]/(1 - e[1])
            return (total[0], total[1] + tail)


def theta4(e):
    """1 + 2 sum_(n>=1) (-1)^n e^(n^2)"""
    total, sign = iv(1), -1
    for term in powers_of_squares(e):
        total = iadd(total, iscale(term, 2*sign))
        sign = -sign
        if term[1] < TINY:
            tail = 2*term[1]/(1 - e[1])
            return (total[0] - tail, total[1] + tail)


def theta2(e4):
    """2 sum_(n>=0) e^((n+1/2)^2), given an enclosure e4 of the fourth root of e"""
    total = iv(0)
    for term in powers_of_squares(e4, shift=1):
        total = iadd(total, iscale(term, 2))
        if term[1] < TINY:
            tail = 2*term[1]/(1 - e4[1])
            return (total[0], total[1] + tail)


def sfl(t, pi_iv):
    """enclosures of S(t), F(t), L(t)"""
    e = pm1.exp_neg_pi(F(t), pi_iv)
    e4 = pm1.exp_neg_pi(F(t)/4, pi_iv)
    t3, t4, t2 = theta3(e), theta4(e), theta2(e4)
    return imul(t3, t3), imul(t4, t4), imul(t2, t2)


def k_contents(u, pi_iv):
    """K_+ and K_- by contents: sum m^2 e^(-pi u m^2), sum (-1)^(m-1) m^2 e^(-pi u m^2)"""
    e = pm1.exp_neg_pi(F(u), pi_iv)
    kp, km, m = iv(0), iv(0), 0
    for term in powers_of_squares(e):
        m += 1
        kp = iadd(kp, iscale(term, m*m))
        km = iadd(km, iscale(term, m*m if m % 2 else -m*m))
        if (m + 1)**2*term[1] < TINY:
            tail = 2*(m + 1)**2*term[1]/(1 - e[1])
            return (kp[0], kp[1] + tail), (km[0] - tail, km[1] + tail)


def k_windings(u, root_u_cubed, pi_iv):
    """K_+ and K_- by windings; root_u_cubed = u^(3/2), rational for the couplings used"""
    u = F(u)
    e = pm1.exp_neg_pi(1/u, pi_iv)
    e4 = pm1.exp_neg_pi(1/(4*u), pi_iv)
    pref = idiv(iv(1), iscale(pi_iv, 4*root_u_cubed))
    sp, n = iv(1), 0
    for term in powers_of_squares(e):
        n += 1
        coef = isub(iv(1), iscale(pi_iv, F(2*n*n)/u))                 # 1 - 2 pi n^2/u
        sp = iadd(sp, iscale(imul(coef, term), 2))
        if (n + 1)**2*term[1]/u < TINY:
            tail = 8*(n + 1)**2*term[1]/(u*(1 - e[1]))
            sp = (sp[0] - tail, sp[1] + tail)
            break
    sm, n = iv(0), 0
    for term in powers_of_squares(e4, shift=1):
        odd = 2*n + 1                                                  # (n + 1/2)^2 = odd^2/4
        coef = isub(iscale(pi_iv, F(odd*odd, 2)/u), iv(1))            # 2 pi (n+1/2)^2/u - 1
        sm = iadd(sm, iscale(imul(coef, term), 2))
        n += 1
        if (2*n + 1)**2*term[1]/u < TINY:
            tail = 8*(2*n + 1)**2*term[1]/(u*(1 - e4[1]))
            sm = (sm[0] - tail, sm[1] + tail)
            break
    return imul(pref, sp), imul(pref, sm)


# ---------------------------------------------------------------- exact integer series in q
def series(order):
    s, f, g = [0]*(order + 1), [0]*(order + 1), [0]*(order + 1)
    r = isqrt(order)
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a*a + b*b
            if n <= order:
                s[n] += 1
                f[n] += -1 if (a + b) % 2 else 1
                g[n] += -1 if b % 2 else 1
    return s, f, g


def smul(a, b, order):
    out = [0]*(order + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j > order:
                    break
                out[i + j] += x*y
    return out


def dilate(a, order):
    """a(q) -> a(q^2)"""
    out = [0]*(order + 1)
    for i, x in enumerate(a):
        if 2*i <= order:
            out[2*i] = x
    return out


def sadd1(a):
    return [a[0] + 1] + a[1:]


def run():
    out, num = {}, {}
    pi_iv = pm1.pi_interval()
    order = 400
    s, f, g = series(order)
    d = [(x - y)//2 for x, y in zip(s, f)]

    # A3 (exact series) -------------------------------------------------------------------------------------
    out['A3_seen_part_is_the_record_of_the_doubled_cell'] = (
        [(x + y)//2 for x, y in zip(s, f)] == dilate(s, order) and all((x + y) % 2 == 0 for x, y in zip(s, f)))
    out['A3_count_of_norm_2n_equals_count_of_norm_n'] = all(s[2*n] == s[n] for n in range(1, order//2 + 1))
    ok = True
    for a in range(-6, 7):
        for b in range(-6, 7):
            z = ds1.G(a, b)*ds1.G(1, 1)                                   # times the diagonal vector
            ok &= z.norm() == 2*(a*a + b*b) and (z.r + z.t) % 2 == 0
    out['A3_diagonal_vector_doubles_every_norm_and_lands_on_the_even_coset'] = ok
    out['A3_mixed_flip_is_the_flip_of_the_doubled_cell'] = g == dilate(f, order)
    # the odd coset cancels under the exchange of a and b
    odd_part = [0]*(order + 1)
    for a in range(-20, 21):
        for b in range(-20, 21):
            if (a + b) % 2 and a*a + b*b <= order:
                odd_part[a*a + b*b] += -1 if b % 2 else 1
    out['A3_mirror_cancels_the_odd_coset'] = all(x == 0 for x in odd_part)
    out['A3_square_of_the_mixed_flip'] = smul(g, g, order) == smul(s, f, order)
    out['A3_counts_are_the_character_mod_four'] = all(
        s[n] == 4*sum(pm1.chi4(k) for k in range(1, n + 1) if n % k == 0) for n in range(1, order + 1))

    # A2 (exact series):  S(q^2)^2 = F(q^2)^2 + D(q)^2 ,  D(q) = L at the doubled cell
    s2, f2 = dilate(s, order), dilate(f, order)
    out['A2_right_triangle_as_series'] = smul(s2, s2, order) == [x + y for x, y in zip(smul(f2, f2, order), smul(d, d, order))]
    out['A4_walk_as_series'] = (
        [2*x for x in s2] == [x + y for x, y in zip(s, f)] and smul(f2, f2, order) == smul(s, f, order)
        )
    # S'^2 - F'^2 = ((S - F)/2)^2 is the same statement
    out['A4_hypotenuse_form'] = [x - y for x, y in zip(smul(s2, s2, order), smul(f2, f2, order))] == smul(d, d, order)

    # A1, A2 (enclosures) -----------------------------------------------------------------------------------
    ok_mirror = ok_tri = True
    widths = []
    for t in (F(2), F(3, 2), F(5)):
        s_a, f_a, l_a = sfl(t, pi_iv)
        s_b, f_b, l_b = sfl(1/t, pi_iv)
        ok_mirror &= overlap(s_b, iscale(s_a, t)) and overlap(f_b, iscale(l_a, t)) and overlap(l_b, iscale(f_a, t))
        ok_tri &= overlap(imul(s_a, s_a), iadd(imul(f_a, f_a), imul(l_a, l_a)))
        ok_tri &= overlap(imul(s_b, s_b), iadd(imul(f_b, f_b), imul(l_b, l_b)))
        widths += [width(s_b), width(f_b), width(l_b)]
    out['A1_mirror_exchanges_flip_and_lost_sum'] = ok_mirror
    out['A2_right_triangle_at_the_points'] = ok_tri
    out['A1_enclosures_narrower_than_1e_minus_100'] = max(widths) < F(1, 10**100)
    wrong = sfl(F(2), pi_iv)
    half = sfl(F(1, 2), pi_iv)
    out['A1_wrong_partner_is_separated'] = not overlap(half[1], iscale(wrong[1], 2))      # F(1/t) against t F(t)
    s1, f1, l1 = sfl(F(1), pi_iv)
    out['A2_on_the_seam_the_legs_are_equal'] = overlap(f1, l1) and overlap(iscale(imul(f1, f1), 2), imul(s1, s1))
    off = sfl(F(11, 10), pi_iv)
    out['A2_off_the_seam_they_are_not'] = not overlap(off[1], off[2]) and below(off[2], off[1])
    num['A2_F_over_S_on_the_seam'] = str(round(float(f1[0]/s1[0]), 15))

    # A4 (enclosures): the walk, the order, the invariant -----------------------------------------------------
    ok_walk = ok_order = True
    for t in (F(1, 4), F(1, 2), F(1), F(2)):
        s_a, f_a, _ = sfl(t, pi_iv)
        s_b, f_b, _ = sfl(2*t, pi_iv)
        ok_walk &= overlap(s_b, iscale(iadd(s_a, f_a), F(1, 2))) and overlap(imul(f_b, f_b), imul(s_a, f_a))
        ok_order &= below(f_a, f_b) and below(f_b, s_b) and below(s_b, s_a)
    out['A4_walk_at_the_points'] = ok_walk
    out['A4_never_across'] = ok_order
    a, b = sfl(F(1, 8), pi_iv)[:2]                                     # a weak-coupling start: F/S about 1e-4
    num['A4_start_F_over_S'] = str(round(float(b[0]/a[0]), 8))
    steps = 0
    while a[0] - b[1] > F(1, 10**50) and steps < 30:
        a, b = iscale(iadd(a, b), F(1, 2)), isqrt_iv(imul(a, b))
        steps += 1
    for _ in range(2):
        a, b = iscale(iadd(a, b), F(1, 2)), isqrt_iv(imul(a, b))
    out['A4_arithmetic_geometric_mean_is_one'] = a[0] <= 1 <= a[1] and width(a) < F(1, 10**50) and overlap(a, b)
    num['A4_steps_to_fifty_digits'] = steps
    other = (iv(F(3, 2)), iv(F(1, 2)))
    for _ in range(8):
        other = (iscale(iadd(*other), F(1, 2)), isqrt_iv(imul(*other)))
    out['A4_control_another_pair_has_another_mean'] = not (other[0][0] <= 1 <= other[0][1])

    # A5: the weak end --------------------------------------------------------------------------------------
    ok_two = ok_lower = True
    weak = {}
    for t in (F(1, 20), F(1, 60)):
        e = pm1.exp_neg_pi(t, pi_iv)
        t4 = theta4(e)
        by_contents = imul(t4, t4)
        q4 = pm1.exp_neg_pi(1/(4*t), pi_iv)                           # Exp(-pi/(4t))
        by_windings = iscale(imul(theta2(q4), theta2(q4)), 1/t)
        ok_two &= overlap(by_contents, by_windings) and by_windings[0] > 0
        first = iscale(imul(q4, q4), 4/t)                             # (4/t) Exp(-pi/(2t))
        ok_lower &= first[0] <= by_windings[1] and by_windings[1] < first[1]*F(1001, 1000)
        digits = len(str(int(1/by_windings[1]))) - 1
        weak[str(t)] = {'flip': f'{float(by_windings[0]):.6e}', 'digits_cancelled_in_the_content_sum': digits,
                        'width_by_contents': f'{float(width(by_contents)):.1e}',
                        'width_by_windings': f'{float(width(by_windings)):.1e}'}
    out['A5_two_descriptions_agree_at_the_weak_end'] = ok_two
    out['A5_first_winding_term_is_the_flip_to_a_thousandth'] = ok_lower
    num['A5_weak_end'] = weak
    # smaller than every power: t^-k F(t) still below 1 for k up to 20 at t = 1/60
    fl = iscale(imul(theta2(pm1.exp_neg_pi(F(15), pi_iv)), theta2(pm1.exp_neg_pi(F(15), pi_iv))), 60)
    out['A5_smaller_than_twenty_powers'] = fl[1]*F(60)**20 < 1
    s_w, f_w, _ = sfl(F(1, 20), pi_iv)
    out['A5_weak_end_is_the_diagonal'] = f_w[1]/s_w[0] < F(1, 10**11)

    # A6: the turn block ------------------------------------------------------------------------------------
    ok_two = ok_pos = ok_weak = True
    flips = {}
    for u, root in ((F(1), F(1)), (F(1, 4), F(1, 8)), (F(1, 16), F(1, 64)), (F(1, 36), F(1, 216))):
        kp_c, km_c = k_contents(u, pi_iv)
        kp_w, km_w = k_windings(u, root, pi_iv)
        ok_two &= overlap(kp_c, kp_w) and overlap(km_c, km_w)
        ok_pos &= km_w[0] > 0 and below(km_w, kp_w)
        flip = idiv(km_w, kp_w)
        lead = iscale(imul(isub(iscale(pi_iv, 1/(2*u)), iv(1)), pm1.exp_neg_pi(1/(4*u), pi_iv)), 2)
        if u <= F(1, 16):
            # K_+ = (1/(4 pi u^(3/2))) (1 + ...), so flip = lead * (1 + small): within a thousandth of the first term
            ok_weak &= flip[0] > lead[0]*F(999, 1000) and flip[1] < lead[1]*F(1001, 1000)
        flips[str(u)] = f'{float(flip[0]):.6e}'
    out['A6_two_descriptions_of_the_turn_block_agree'] = ok_two
    out['A6_flip_of_a_closed_surface_is_positive_and_below_one'] = ok_pos
    out['A6_weak_end_is_the_first_winding_term'] = ok_weak
    num['A6_flip_at_u'] = flips
    # (the exponent pi/(4u), u = T/(4 pi), is pi^2/T; with t = 2/kappa per face and F faces it is CZ1's
    #  pi^2 kappa/(2F): algebra, stated in the text, not a check)

    # A8: the turn block against the walk of the pair -------------------------------------------------------
    ok, ratios = True, {}
    for u in (F(1, 100), F(1, 72), F(1, 36), F(1, 16), F(1, 8), F(1, 4), F(1, 2), F(1), F(2)):
        kp1, km1 = k_contents(u, pi_iv)
        kp2, km2 = k_contents(2*u, pi_iv)
        rho, rho2 = idiv(km1, kp1), idiv(km2, kp2)
        lhs = imul(imul(rho2, rho2), imul(iadd(iv(1), rho), iadd(iv(1), rho)))     # rho'^2 (1 + rho)^2
        rhs = iscale(rho, 4)                                                       # 4 rho
        ok &= rho[0] > 0 and below(rhs, lhs)
        ratios[str(u)] = f'{float(lhs[0]/rhs[1])**0.5:.6g}'
    out['A8_turn_block_flip_walks_at_least_as_fast_as_the_mean'] = ok
    num['A8_ratio_to_the_mean_step'] = ratios
    # strong end: series in q of  K_-/K_+  to order 12 and of the mean step
    def ser(coeffs, order=14):
        return [F(coeffs.get(k, 0)) for k in range(order + 1)]

    def mul(a, b, order=14):
        o = [F(0)]*(order + 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if i + j > order:
                        break
                    o[i + j] += x*y
        return o

    def inv(a, order=14):
        o = [F(0)]*(order + 1)
        o[0] = 1/a[0]
        for n in range(1, order + 1):
            o[n] = -sum(a[k]*o[n - k] for k in range(1, n + 1))/a[0]
        return o
    kp_s = ser({0: 1, 3: 4, 8: 9})                       # K_+/q = 1 + 4 q^3 + 9 q^8 + 16 q^15 ...
    km_s = ser({0: 1, 3: -4, 8: 9})
    rho_s = mul(km_s, inv(kp_s))
    rho2_s = [rho_s[k//2] if k % 2 == 0 else F(0) for k in range(15)]              # q -> q^2
    # mean step m = 2 sqrt(rho)/(1 + rho):  m^2 (1 + rho)^2 = 4 rho ; solve for the series m with m(0) = 1
    target = mul([4*x for x in rho_s], inv(mul([F(1)] + [F(0)]*14, mul(sadd1(rho_s), sadd1(rho_s)))))
    msq = target
    mser = [F(1)] + [F(0)]*14
    for n in range(1, 15):
        mser[n] = (msq[n] - sum(mser[k]*mser[n - k] for k in range(1, n)))/2
    diff = [x - y for x, y in zip(rho2_s, mser)]
    out['A8_strong_end_first_difference_is_64_q12'] = all(x == 0 for x in diff[:12]) and diff[12] == 64
    G = ds1.G

    def hankel(values, weights, size, shift=0):
        # c(s + t + shift) = sum w lambda^(s+t+shift): DS1-D4's mirror form for the mirror lambda -> lambda-dagger
        return [[sum((w*(v**(i + j + shift)) for v, w in zip(values, weights)), G(0))
                 for j in range(size)] for i in range(size)]
    real = [G(F(1, 2)), G(F(1, 3)), G(F(1, 5)), G(F(1, 7))]
    out['A7_real_transfer_values_give_a_positive_reflection_form'] = ds1.inertia(hankel(real, [F(1)]*4, 6)) == (4, 0, 2)
    pair = [G(F(1, 2)), G(F(1, 3), F(1, 4)), G(F(1, 3), F(-1, 4))]
    out['A7_a_pair_off_the_axis_is_one_count'] = ds1.inertia(hankel(pair, [F(1)]*3, 5)) == (2, 1, 2)
    ok = True
    for r in (F(3, 5), F(2, 5), F(1, 4), F(1, 6), F(1, 10)):
        h0, h1 = hankel(real, [F(1)]*4, 4), hankel(real, [F(1)]*4, 4, shift=1)
        count = ds1.inertia(ds1.madd([[r*x for x in row] for row in h0], h1, -1))[1]
        ok &= count == sum(1 for v in real if v.r > r)
    out['A7_gap_is_the_same_count_with_the_line_moved'] = ok

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out, num


if __name__ == '__main__':
    out, num = run()
    with open(os.path.join(HERE, 'AG1_RESULT.json'), 'w') as fh:
        json.dump({'checks': out, 'numbers': num}, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
    print(json.dumps(num, indent=1))
