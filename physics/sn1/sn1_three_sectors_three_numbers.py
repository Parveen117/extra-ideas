"""SN1: the spread of a count through a cut, in the three sectors: 1/4, 1/2, 1/3.

FD1-F4: an exclusive record seen through a cut counts like independent coins whose biases q are the readings of
the seen kernel; mean = sum q, spread = sum q(1 - q) (QD1's F = R - D, the share law).  The pure number is
spread / mean.  Here the seen share is taken in each of the three sectors of the algebra (R^2 = -1, eps^2 = 0,
K^2 = +1) with its angle spread evenly.  sympy.  Python 3.12."""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def run():
    out, num = {}, {}
    zz = lambda e: sp.simplify(e) == 0
    ph, x, eta, s, X = sp.symbols('phi x eta s X', positive=True)

    sectors = {
        'turn':  (sp.cos(ph)**2, ph, sp.pi/2),          # seen cos^2, lost sin^2 (IN1); angle on the quarter circle
        'shear': (1/(1 + x**2), x, X),                  # lost / seen = x^2
        'boost': (1/sp.cosh(eta)**2, eta, s),           # seen 1 - tanh^2, lost tanh^2 (GR1, DG1-D4)
    }
    res = {}
    for name, (q, var, top) in sectors.items():
        if name == 'boost':                                  # antiderivatives in u = tanh(eta), checked by differentiation
            T = sp.tanh(var)
            prim = (T, T**3/3, 2*T**5/5 - T**3/3)
            integrands = (q, q*(1 - q), q*(1 - q)*(1 - 2*q))
            assert all(sp.simplify((sp.diff(P_, var) - f_).rewrite(sp.exp)) == 0 for P_, f_ in zip(prim, integrands))
            Ts = sp.tanh(top)
            res[name] = (Ts**2/3, 2*Ts**4/5 - Ts**2/3)
            continue
        mean = sp.integrate(q, (var, 0, top))
        second = sp.integrate(sp.simplify(q*(1 - q)), (var, 0, top))
        third = sp.integrate(sp.simplify(q*(1 - q)*(1 - 2*q)), (var, 0, top))
        res[name] = (sp.simplify(second/mean), sp.simplify(third/mean))

    # all three are one form: lost / seen = t^2 with t = tan(phi), x, sinh(eta)
    out['S0_one_form'] = zz(1/(1 + sp.tan(ph)**2) - sp.cos(ph)**2) and zz(1/(1 + sp.sinh(eta)**2) - 1/sp.cosh(eta)**2)

    out['S1_turn_sector_one_quarter'] = res['turn'][0] == sp.Rational(1, 4) and res['turn'][1] == 0
    f2, f3 = res['shear']
    out['S2_shear_sector_one_half'] = sp.limit(f2, X, sp.oo) == sp.Rational(1, 2) and sp.limit(f3, X, sp.oo) == sp.Rational(1, 4)
    out['S2_shear_closed_form'] = zz(f2 - (1 - X/((1 + X**2)*sp.atan(X)))/2)
    g2, g3 = res['boost']
    out['S3_boost_sector_one_third'] = g2 == sp.tanh(s)**2/3 and sp.limit(g2, s, sp.oo) == sp.Rational(1, 3) and sp.limit(g3, s, sp.oo) == sp.Rational(1, 15)
    # the law of the seen share in the turn sector: the arcsine law, mean 1/2
    qv = sp.symbols('q', positive=True)
    dens = 1/(sp.pi*sp.sqrt(qv*(1 - qv)))
    out['S1_arcsine_law'] = sp.integrate(dens, (qv, 0, 1)) == 1 and sp.integrate(qv*dens, (qv, 0, 1)) == sp.Rational(1, 2) and sp.integrate(qv*(1 - qv)*dens, (qv, 0, 1)) == sp.Rational(1, 8)
    # a single cut on the diagonal: q = 1/2, spread / mean = 1/2 ; a cut along the reading: 0 ; nothing seen: 1
    out['S4_single_cut'] = [(qq*(1 - qq)/qq) for qq in (sp.Rational(1, 2), sp.Integer(1))] == [sp.Rational(1, 2), 0] and sp.limit((qv*(1 - qv))/qv, qv, 0) == 1
    # how fast the boost sector reaches 1/3
    num['boost_sector_at_angle'] = {str(v): float((sp.tanh(v)**2/3)) for v in (1, 2, 3, 5)}
    num['numbers'] = {'turn': '1/4, 0', 'shear': '1/2, 1/4', 'boost': '1/3, 1/15'}
    num['measured'] = {'turn (open cavities, arXiv:cond-mat/0009087)': '1/4', 'boost (long wires, arXiv:cond-mat/9808042)': '1/3, reached for large reservoirs',
                       'shear (symmetric double barrier)': '1/2 (recalled)'}

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'SN1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
