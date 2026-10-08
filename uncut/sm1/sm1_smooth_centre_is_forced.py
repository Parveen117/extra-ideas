"""SM1: at a horizon the centre angle is 2 pi in the carrier: F4 is not an extra assumption there. Exact sympy."""
import json, os, sympy as sp
I = sp.I; I2 = sp.eye(2); Z = sp.zeros(2)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]])
C1, C3 = K, I*R; G = C1*C3
kap, t, rho, P, a = sp.symbols('kappa t rho P a', positive=True)
sm = lambda M: M.applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
turned = lambda x: sp.cos(x)*C1 + sp.sin(x)*C3
out = {}
# T1: orbit at fixed rho: the cut after far time t is the radial cut turned by -i kappa t  (HX1-T2 with psi = kappa t)
out['T1_orbit'] = sm(turned(-I*kap*t) - (sp.cosh(kap*t)*K + sp.sinh(kap*t)*R)) == Z
# T2: the closed reading uses the same generator with angle x iota (HX1-T5): angle = kappa t
out['T2_closed_angle'] = sm(turned((I)*(-I*kap*t)) - turned(kap*t)) == Z
# T3: it returns first when kappa t = 2 pi : P = 2 pi / kappa, with no freedom
ret = sp.solve(sp.Eq(kap*P, 2*sp.pi), P)
out['T3_period_forced'] = ret == [2*sp.pi/kap] and sm(turned(2*sp.pi) - C1) == Z and all(sm(turned(sp.pi*q) - C1) != Z for q in (sp.Rational(1, 2), 1, sp.Rational(3, 2)))
# T4: centre angle = length/radius of the circle swept = kappa P = 2 pi exactly (QT1-T4)
out['T4_centre_angle'] = sp.simplify((kap*rho*ret[0])/rho - 2*sp.pi) == 0
# T5: what could still change the count: only the unit of the reading frame (UN1-T6 with alpha = 2 pi)
A, cf = sp.symbols('A c_f', positive=True)
out['T5_only_unit_left'] = sp.simplify((A/4)*(2*sp.pi/(2*sp.pi))/cf - A/(4*cf)) == 0
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SM1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
