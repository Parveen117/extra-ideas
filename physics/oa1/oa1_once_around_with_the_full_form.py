"""OA1: the once-around turn with the full count form of MO1 (thesis, 'Next' item 3).

A direction carried once around a circle of radius r at angular rate Omega, in the form
    dtau^2 = dt^2 - (dr + beta dt)^2 - r^2 dphi^2 ,  beta^2 = r_s / r        (MO1; MA1),
with the connection of the form that has no torsion (CV1's reading), the carried direction kept across the history.
Pure numbers: x = r_s / r (the memory m), y = (r Omega)^2.   Symbolic (sympy)."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def transport():
    t, r, ph = sp.symbols('t r phi', real=True)
    rs, Om, U = sp.symbols('r_s Omega U', positive=True)          # U = (u^t)^2
    beta = sp.sqrt(rs/r); X = [t, r, ph]
    g = sp.Matrix([[1 - beta**2, -beta, 0], [-beta, -1, 0], [0, 0, -r**2]])
    gi = g.inv()
    G = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(3))/2)
           for c in range(3)] for b in range(3)] for a in range(3)]
    ut = sp.sqrt(U)
    u = sp.Matrix([ut, 0, ut*Om])
    acc = sp.Matrix([sum(G[m][i][j]*u[i]*u[j] for i in range(3) for j in range(3)) for m in range(3)])
    al = g*acc
    M = sp.zeros(3)                                                # dS/dtau = M S : transport, direction kept across the history
    for m in range(3):
        for n in range(3):
            M[m, n] = -sum(G[m][i][n]*u[i] for i in range(3)) - u[m]*al[n]
    lam = sp.symbols('lambda')
    cp = sp.expand((M - lam*sp.eye(3)).det())
    w2 = sp.simplify(-cp.coeff(lam, 1))                            # -lambda^3 - w2 lambda
    Uv = 1/(1 - rs/r - r**2*Om**2)                                 # unit count along the history
    w2 = sp.simplify(w2.subs(U, Uv))
    per2 = sp.simplify(w2/(Om**2*Uv))                              # (angle per circuit in the co-turning basis / 2 pi)^2
    return r, rs, Om, sp.simplify(cp.coeff(lam, 0).subs(U, Uv)), sp.simplify(cp.coeff(lam, 2).subs(U, Uv)), per2


def run():
    out = {}
    z = lambda e: sp.simplify(e) == 0
    r, rs, Om, c0, c2, per2 = transport()
    x, y = sp.symbols('x y', positive=True)
    p2 = sp.simplify(per2.subs({rs: x*r, Om: sp.sqrt(y)/r}))
    out['T0_pure_turn'] = z(c0) and z(c2)
    # T1: angle per circuit in the co-turning basis = 2 pi (1 - 3x/2)/sqrt(1 - x - y)
    out['T1_general'] = z(p2 - (1 - sp.Rational(3, 2)*x)**2/(1 - x - y))
    turn = lambda xx, yy: 1 - (1 - sp.Rational(3, 2)*xx)/sp.sqrt(1 - xx - yy)        # turn against far directions, in whole turns
    # T2: in the line's quantities: N^2 = 1 - x, static acceleration g = x/(2 r N) (GR1), local speed^2 = y/N^2
    N, v = sp.symbols('N v', positive=True)
    lhs = turn(1 - N**2, v**2*N**2)
    out['T2_line_form'] = z(lhs - (1 - (N - (1 - N**2)/(2*N))/sp.sqrt(1 - v**2)))   # 1 - gamma (N - r g)
    # T3: free circle (MO1's period law y = x/2): 1 - sqrt(1 - 3x/2); first order (3/2) pi x ; second (9/16) pi x^2
    free = turn(x, x/2)
    out['T3_free_circle'] = z(free - (1 - sp.sqrt(1 - sp.Rational(3, 2)*x)))
    ser = sp.series(2*sp.pi*free, x, 0, 3).removeO()
    out['T3_orders'] = z(ser - (sp.Rational(3, 2)*sp.pi*x + sp.Rational(9, 16)*sp.pi*x**2))
    # T4: GR2-G4 (radial boost only, static frames): 2 pi (1/N - 1) = pi x + (3/4) pi x^2 : the first order of the full form is 3/2 of it
    gr2 = sp.series(2*sp.pi*(1/sp.sqrt(1 - x) - 1), x, 0, 3).removeO()
    out['T4_GR2_orders'] = z(gr2 - (sp.pi*x + sp.Rational(3, 4)*sp.pi*x**2))
    out['T4_ratio_first_order'] = sp.limit(2*sp.pi*free/(2*sp.pi*(1/sp.sqrt(1 - x) - 1)), x, 0) == sp.Rational(3, 2)
    # T5: carried slowly on a held circle: (1 - N) + m/(2N): the cone part plus radius x static acceleration
    slow = turn(x, 0)
    out['T5_slow'] = z(slow.subs(x, 1 - N**2) - ((1 - N) + (1 - N**2)/(2*N)))
    # T6: no field: 1 - gamma (the turn of a direction carried round a circle at speed v)
    out['T6_no_field'] = z(turn(0, v**2) - (1 - 1/sp.sqrt(1 - v**2)))
    # T7: the carried direction keeps step with the circle exactly at x = 2/3 (r = 3 r_s / 2), for every speed
    out['T7_light_circle'] = z(turn(sp.Rational(2, 3), y) - 1)
    # numbers: circuit of the Earth at 7020 km (the thesis's example)
    GM, c, rr = 3.986004418e14, 299792458.0, 7.020e6
    xe = 2*GM/c**2/rr
    per_circuit_mas = 1.5*3.141592653589793*xe*206264806.247
    circuits_per_year = 31557600.0/(2*3.141592653589793*(rr**3/GM)**0.5)
    out_num = dict(x=xe, per_circuit_mas=round(per_circuit_mas, 4), GR2_per_circuit_mas=round(per_circuit_mas/1.5, 4),
                   per_year_arcsec=round(per_circuit_mas*circuits_per_year/1000, 3))
    # the same formula at the radius recalled for the mission's orbit (7027.4 km: from memory, not verified at source)
    out_num['per_year_arcsec_at_7027p4km'] = round(out_num['per_year_arcsec']*(7020.0/7027.4)**2.5, 3)
    # the Earth-Moon pair as a carried direction on its circuit of the Sun (1 circuit per year)
    x_sun = 2*1.32712440018e20/c**2/1.495978707e11
    out_num['sun_circuit_mas_per_year'] = round(2*3.141592653589793*(1 - (1 - 1.5*x_sun)**0.5)*206264806.247, 2)
    out = {k: bool(val) for k, val in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=out_num), open(os.path.join(HERE, 'OA1_RESULT.json'), 'w'), indent=1)
    return out, out_num


if __name__ == '__main__':
    print(run())
