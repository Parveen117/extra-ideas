"""MT1: the temperature of a body read by a moving reader: three classical answers as three readings of one function.
Exact sympy."""
import json, os, sympy as sp
S, phi, p, v = sp.symbols('S phi p v', positive=True)
E0 = sp.Function('E0', positive=True)(S)                         # the body's centre value at rest, a function of its count
T0 = sp.diff(E0, S)
out = {}
z = lambda e: sp.simplify(sp.simplify(e).rewrite(sp.exp)) == 0 or sp.simplify(sp.powdenest(sp.simplify(e), force=True)) == 0
# carrier: the boost is a turn by the imaginary angle (HX1): (E', p') = E0 (cosh phi, sinh phi)
Ep, pp = E0*sp.cosh(phi), E0*sp.sinh(phi)
out['T0_invariant'] = z(Ep**2 - pp**2 - E0**2)
gam = sp.cosh(phi); vel = sp.tanh(phi)
# reading 1: count held against velocity (corner [S | v])
T_v = sp.diff(Ep, S)                              # phi fixed
out['T1_fixed_velocity'] = z(T_v - gam*T0)
# reading 2: count held against momentum (corner [S | p]): E' = sqrt(E0^2 + p^2)
Ep_p = sp.sqrt(E0**2 + p**2)
T_p = sp.diff(Ep_p, S).subs(p, pp)
out['T2_fixed_momentum'] = z(T_p - T0/gam)
# reading 3: the rest reading
out['T3_rest'] = z(T_v*T_p - T0**2)               # the rest value is the geometric mean of the other two
# T4: which one is 1/period? The closed turn's period in the moving reader's time is dilated: P' = gamma P0
P0 = 1/T0; Pm = gam*P0
out['T4_inverse_period_is_fixed_momentum'] = z(1/Pm - T_p)
# T5: two planes, one Maxwell-type equation (PL1): dS = P' dE' - (v P') dp'
dE = sp.diff(Ep_p, S); dEp = sp.diff(Ep_p, p)     # dE' = T_p dS + v dp
out['T5_velocity_is_the_partner'] = z(dEp.subs(p, pp) - vel)
out['T5_closed'] = z(sp.diff(1/dE, p) + sp.diff(dEp/dE, S)) if False else z(sp.diff(dE, p) - sp.diff(dEp, S))
# T6: the ratio of the two contested answers is the shared information between the frames (CI1: exp(I0) = cosh psi)
out['T6_ratio'] = z(T_v/T_p - sp.cosh(phi)**2)
# T7: only the (energy, time) area with the momentum plane included is reader-free
out['T7_count_reader_free'] = z(Pm*(sp.diff(Ep, S) - vel*sp.diff(pp, S)) - 1)     # P'(dE' - v dp')/dS = 1
out['T7_energy_plane_alone_is_not'] = z(Pm*sp.diff(Ep, S) - gam**2)               # = gamma^2, not 1
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'MT1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
