"""NA1: a plane whose turns do not commute (the three cuts). The frame part of the master equation gains a
quadratic term; its flux over one eighth of the boundary is WD2's one-eighth turn. Exact sympy."""
import json, os, itertools, sympy as sp
I = sp.I; I2 = sp.eye(2); Z = sp.zeros(2)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); S = R*K
C = [K, S, I*R]
out = {}
zm = lambda M: M.applyfunc(lambda e: sp.simplify(sp.expand(sp.simplify(e).rewrite(sp.exp)))) == Z
eps = lambda a, b, c: sp.LeviCivita(a, b, c)
# T1: the three cuts close among themselves: [C_a, C_b] = 2 iota s eps_abc C_c, one sign s for all
comm = lambda A, B: A*B - B*A
s_ = sp.simplify((comm(C[0], C[1])*C[2])[0, 0]/(2*I))
out['T1_closure'] = s_**2 == 1 and all(zm(comm(C[a], C[b]) - 2*I*s_*sum((eps(a, b, c)*C[c] for c in range(3)), Z)) for a in range(3) for b in range(3))
# T2: partner one-forms a^c_i(q) on two base variables; A_i = (iota/2) sum_c a^c_i C_c.
x, y = sp.symbols('x y', real=True)
a = [[sp.Function(f'a{c}{i}')(x, y) for i in range(2)] for c in range(3)]     # a[c][i]
A = [sum(((I/2)*a[c][i]*C[c] for c in range(3)), Z) for i in range(2)]
F = sp.diff(A[1], x) - sp.diff(A[0], y) + comm(A[0], A[1])                    # frame part: dA + A^A
# component form: F = (iota/2) sum_c f^c C_c with f^c = d a^c - s * eps_cab a^a_x a^b_y
f = [sp.diff(a[c][1], x) - sp.diff(a[c][0], y) - s_*sum(eps(c, p, q)*a[p][0]*a[q][1] for p in range(3) for q in range(3)) for c in range(3)]
out['T2_quadratic_term'] = zm(F - sum(((I/2)*f[c]*C[c] for c in range(3)), Z))
# T3: commuting sub-case: only one cut active -> the quadratic term vanishes, ME1's pair equation is left
sub1 = {a[1][0]: 0, a[1][1]: 0, a[2][0]: 0, a[2][1]: 0}
out['T3_one_cut_is_ME1'] = sp.simplify(f[0].subs(sub1) - (sp.diff(a[0][1], x) - sp.diff(a[0][0], y))) == 0 and all(sp.simplify(f[c].subs(sub1)) == 0 for c in (1, 2))
# T4: a pure change of reader g(x,y) = Exp(iota x C1/2) Exp(iota y C2/2): A = g^-1 dg has no frame part, although its pieces do not commute
Ex = lambda t, Cc: sp.cos(t/2)*I2 + I*sp.sin(t/2)*Cc
g = Ex(x, C[0])*Ex(y, C[1]); gi = g.inv()
Ag = [(gi*sp.diff(g, x)).applyfunc(sp.simplify), (gi*sp.diff(g, y)).applyfunc(sp.simplify)]
Fg = sp.diff(Ag[1], x) - sp.diff(Ag[0], y) + comm(Ag[0], Ag[1])
out['T4_change_of_reader_is_flat'] = zm(Fg) and not zm(comm(Ag[0], Ag[1]))
# T5: the reading along a direction n(theta, phi) keeps one partner form a = (1 - cos theta)/2 dphi; its frame part
#     is half the area form, and its flux over one eighth of the boundary is pi/4 - WD2-T3's one-eighth turn
th, ph = sp.symbols('theta phi', real=True)
flux_density = sp.diff((1 - sp.cos(th))/2, th)
out['T5_half_area_form'] = sp.simplify(flux_density - sp.sin(th)/2) == 0
octant = sp.integrate(sp.integrate(flux_density, (th, 0, sp.pi/2)), (ph, 0, sp.pi/2))
out['T5_octant_is_one_eighth_turn'] = sp.simplify(octant - sp.pi/4) == 0 and sp.simplify(sp.expand_complex(sp.exp(I*octant)) - (1 + I)/sp.sqrt(2)) == 0
whole = sp.integrate(sp.integrate(flux_density, (th, 0, sp.pi)), (ph, 0, 2*sp.pi))
out['T5_whole_boundary_one_turn'] = sp.simplify(whole - 2*sp.pi) == 0
# T6: count of flatness conditions for m base variables and a plane of g non-commuting turns: g * m(m-1)/2
out['T6_counts'] = [3*m*(m - 1)//2 for m in (2, 3, 4)] == [3, 9, 18]
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'NA1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
