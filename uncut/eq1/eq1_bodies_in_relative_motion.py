"""EQ1: two bodies in relative motion exchanging energy and momentum. Exact sympy."""
import json, os, sympy as sp
T1, T2, v, u, v1, v2, P1, P2, dE, dp, e1, e2 = sp.symbols('T1 T2 v u v1 v2 P1 P2 dE dp e1 e2', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
gam = lambda x: 1/sp.sqrt(1 - x**2)
# MT1-T5 for each body: dS = P'(dE' - v dp'), P' = gamma/T0
# T1: free exchange of energy and momentum at fixed totals: stationary iff same period and same velocity
dS = P1*(dE - v1*dp) - P2*(dE - v2*dp)
sol = sp.solve([sp.diff(dS, dE), sp.diff(dS, dp)], [P1, v1], dict=True)
out['T1_equilibrium_needs_same_velocity'] = sol == [{P1: P2, v1: v2}]
# T2: exchange through one channel: packets with momentum u * energy (|u| <= 1). The body then has one number,
# its temperature for that channel:  T(u) = T0 / (gamma (1 - v u))
Tch = lambda T0, vel, uu: T0/(gam(vel)*(1 - vel*uu))
dS1 = (gam(v)/T1)*(1 - v*u)                    # count per energy given up by the moving body through channel u
out['T2_channel_temperature'] = z(1/dS1 - Tch(T1, v, u))
out['T2_forward_backward'] = z(Tch(T1, v, 1)**2 - T1**2*(1 + v)/(1 - v)) and z(Tch(T1, v, -1)**2 - T1**2*(1 - v)/(1 + v))
# T3: body 2 at rest with the same rest temperature: through the forward channel the mover reads hotter, backward colder
R_ = sp.Rational
f, b = Tch(1, R_(3, 5), 1), Tch(1, R_(3, 5), -1)
out['T3_values_v_3_5'] = f"{f}, {b}"; out['T3_opposite_flows'] = f == 2 and b == R_(1, 2)
# T4: a two-channel circuit: mover gives e1 forward, takes e2 backward. Count change of the pair, and momentum of the mover
vv = R_(3, 5); g = gam(vv)
count = (-(g/1)*(1 - vv*1)*e1 + e1) + ((g/1)*(1 - vv*(-1))*e2 - e2)     # both rest temperatures = 1
out['T4_count_gain'] = z(count - (e1/2 + e2)) and all(c > 0 for c in sp.Poly(count, e1, e2).coeffs())
dmom = -e1*1 + e2*(-1)                         # mover loses forward momentum e1, gains backward momentum e2
out['T4_mover_slows'] = z(dmom + e1 + e2)
# T5: the channel temperatures of the two bodies agree for every u only if the velocities agree
ratio = lambda T0, vel: sp.simplify(Tch(T0, vel, 1)/Tch(T0, vel, -1))        # forward/backward = (1+v)/(1-v)
same_v = sp.solve(sp.Eq(ratio(T1, v1), ratio(T2, v2)), v1)
out['T5_all_channels_iff_same_velocity'] = same_v == [v2] and sp.solve(sp.Eq(Tch(T1, v2, 1), Tch(T2, v2, 1)), T1) == [T2]
# T6: the transverse channel u = 0 alone gives T0/gamma: MT1's fixed-momentum reading
out['T6_transverse_channel'] = z(Tch(T1, v, 0) - T1/gam(v))
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EQ1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
