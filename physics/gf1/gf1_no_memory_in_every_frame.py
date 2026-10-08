"""GF1: the two no-memory conditions of MM1 for every frame, not only on flat slices.

A coordinate-free quadratic law in the order defect c_ab^k of a frame is Q = a1 I1 + a2 I2 + a3 I3 (TP1).
The conditions are algebraic in the order defect at a place, so they can be put to every frame:
  (P) an order defect that lies in one plane (two directions only) leaves no memory;
  (T) an order defect that is a pure turn rate of the cuts leaves no memory.
Symbolic (sympy), with the functions of CV1 and TP1.  Python 3.12."""
import itertools, json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)
z = lambda e: sp.simplify(e) == 0
a1, a2, a3 = sp.symbols('a1 a2 a3')


def general_defect():
    """all 24 components c_ab^k = -c_ba^k as free symbols."""
    c = [[[sp.Integer(0)]*4 for _ in range(4)] for _ in range(4)]
    syms = {}
    for a in range(4):
        for b in range(a + 1, 4):
            for k in range(4):
                s = sp.Symbol('c%d%d_%d' % (a, b, k))
                syms[(a, b, k)] = s
                c[a][b][k] = s; c[b][a][k] = -s
    return c, syms


def Qgen(c):
    I1, I2, I3 = tp.invariants(c)
    return sp.expand(a1*I1 + a2*I2 + a3*I3)


def flat(E):
    c = cv.structure(E); w = cv.connection(c); R = cv.curvature(E, c, w)
    return all(v == 0 for v in R.values()), c


def run():
    out = {}
    c, syms = general_defect()
    Q = Qgen(c)
    planes = list(itertools.combinations(range(4), 2))
    # P: order defect confined to one plane (a, b): only c_ab^a and c_ab^b
    ok, forms = True, {}
    for (a, b) in planes:
        keep = {syms[(a, b, a)], syms[(a, b, b)]}
        q = Q.subs({s: 0 for s in syms.values() if s not in keep})
        fac = sp.factor(q)
        ratio = sp.simplify(q/(2*a1 + a2 + a3))
        ok &= (not ratio.has(a1, a2, a3)) and ratio != 0
        forms['%d%d' % (a, b)] = str(ratio)
    out['P_one_plane'] = ok
    # T: pure turn rate of the cuts along e_0: c_0i^j = -c_0j^i
    w1, w2, w3 = sp.symbols('w1 w2 w3')
    sub = {s: 0 for s in syms.values()}
    sub.update({syms[(0, 1, 2)]: w3, syms[(0, 2, 1)]: -w3, syms[(0, 2, 3)]: w1, syms[(0, 3, 2)]: -w1, syms[(0, 3, 1)]: w2, syms[(0, 1, 3)]: -w2})
    qt = Q.subs(sub)
    ratio_t = sp.simplify(qt/(2*a1 - a2))
    out['T_turn_rate'] = (not ratio_t.has(a1, a2, a3)) and z(ratio_t - 2*(w1**2 + w2**2 + w3**2))
    # both conditions: 1 : 2 : -4, for every frame
    sol = sp.solve([2*a1 + a2 + a3, 2*a1 - a2], [a2, a3], dict=True)[0]
    out['PT_select_TP1'] = sol[a2] == 2*a1 and sol[a3] == -4*a1
    # the selected law has no square of any in-plane component: memory only between different planes
    Qs = Q.subs({a1: sp.Rational(1, 4), a2: sp.Rational(1, 2), a3: -1})
    inplane = [syms[(a, b, k)] for (a, b) in planes for k in (a, b)]
    out['S_no_squares_of_in_plane_parts'] = all(Qs.coeff(s, 2) == 0 for s in inplane)
    # two in-plane parts meet in Q only if their planes share exactly one direction
    ok = True
    for (p, q) in itertools.combinations([(a, b, k) for (a, b) in planes for k in (a, b)], 2):
        co = Qs.coeff(syms[p], 1).coeff(syms[q], 1)
        shared = len(set(p[:2]) & set(q[:2]))
        if co != 0:
            ok &= (shared == 1)
    out['S_cross_terms_only_between_planes_sharing_a_direction'] = ok
    # frames of flat space (cv.X used as Cartesian or cylindrical names)
    t, x, y, zz = cv.X
    g, om = sp.symbols('g omega', positive=True)
    sel = {a1: sp.Rational(1, 4), a2: sp.Rational(1, 2), a3: -1}
    # (a) uniformly accelerated frame: order defect in the plane (time, x)
    Ea = sp.Matrix([[1/(1 + g*x), 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    fa, ca = flat(Ea)
    qa = Qgen(ca)
    out['F_accelerated_frame'] = fa and z(qa.subs(sel)) and not z(qa.subs({a1: 1, a2: 0, a3: 0}))
    # (b) frame turning in time about z
    Eb = sp.Matrix([[1, 0, 0, 0], [0, sp.cos(om*t), sp.sin(om*t), 0], [0, -sp.sin(om*t), sp.cos(om*t), 0], [0, 0, 0, 1]])
    fb, cb = flat(Eb)
    qb = Qgen(cb)
    out['F_turning_frame'] = fb and z(qb.subs(sel)) and z(sp.simplify(qb/(2*a1 - a2)).diff(a1)) and not z(qb.subs({a1: 1, a2: 0, a3: 0}))
    # (c) polar frame in one plane of space (x -> radius, y -> angle): order defect in that plane
    Ec = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1/x, 0], [0, 0, 0, 1]])
    cc = cv.structure(Ec)
    qc = Qgen(cc)
    out['F_polar_frame_one_plane'] = z(qc.subs(sel)) and not z(qc.subs({a1: 1, a2: 0, a3: 0}))
    # (d) spherical frame (two planes): the selected law is not zero there, it is a pure boundary term: Q = -2 div(trace vector)
    r, th, ph = x, y, zz
    Ed = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1/r, 0], [0, 0, 0, 1/(r*sp.sin(th))]])
    cd = cv.structure(Ed)
    qd = sp.simplify(Qgen(cd).subs(sel))
    tr = [sum(cd[a][b][a] for a in range(4)) for b in range(4)]
    vol = r**2*sp.sin(th)
    V = [sum(cv.ETA[b, b]*tr[b]*Ed[b, mu] for b in range(4)) for mu in range(4)]
    div = sp.simplify(sum(sp.diff(vol*V[mu], cv.X[mu]) for mu in range(4))/vol)
    out['F_spherical_frame_boundary_term'] = (not z(qd)) and (z(qd + 2*div) or z(qd - 2*div))
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, one_plane_forms=forms, turn_form=str(ratio_t)), open(os.path.join(HERE, 'GF1_RESULT.json'), 'w'), indent=1)
    return out, forms, str(ratio_t)


if __name__ == '__main__':
    print(run())
