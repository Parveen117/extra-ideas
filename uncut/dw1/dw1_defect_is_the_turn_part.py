"""DW1: the ON1 defect, read on the centre potential, is the turn-part w of UP5. Exact sympy."""
import json, os, sympy as sp
S, V, a, b = sp.symbols('S V a b', positive=True)
Ug = sp.Function('U')(S, V); f = sp.Function('f')(S, V)
def stage(U):
    T, P = sp.diff(U, S), -sp.diff(U, V)
    DS = lambda g: S*sp.diff(g, S); DV = lambda g: V*sp.diff(g, V)
    eps = DS(P)/DV(P)
    ESP = lambda g: DS(g) - eps*DV(g)          # E[S|P]
    EVP = lambda g: DV(g) - DS(g)/eps          # E[V|P]
    return T, P, DS, DV, eps, ESP, EVP
z = lambda e: sp.simplify(e.doit()) == 0
out = {}
T, P, DS, DV, eps, ESP, EVP = stage(Ug)
c1 = lambda g: ESP(DV(g)) - DV(ESP(g))          # [E[S|P], E[V|S]]
c2 = lambda g: EVP(DS(g)) - DS(EVP(g))          # [E[V|P], E[S|V]]
# T1: UP5-T5 form with alpha=0, F=0, beta = running
out['T1_frame_part'] = z(c1(f) - DV(eps)*DV(f)) and z(c2(f) - DS(1/eps)*DS(f))
w1, w2 = c1(Ug)/2, c2(Ug)/2
# T2: the turn-part is (running) x (work reading)
out['T2_w'] = z(w1 + DV(eps)*P*V/2) and z(w2 - DS(1/eps)*T*S/2)
# self-similar centre
G = sp.Function('G'); Us = S**(1/a)*G(V/S**(b/a))
T, P, DS, DV, eps, ESP, EVP = stage(Us)
D = lambda g: a*DS(g) + b*DV(g)
w1 = -DV(eps)*P*V/2; w2 = DS(1/eps)*T*S/2
# T3: the two turn-parts are one: a eps^2 PV w2 + ... tied by the weights
out['T3_one_w'] = z(a*eps**2*P*V*w2 + b*T*S*w1)
# T4: the pure turn-part w/U is constant along the scaling flow
out['T4_flow_invariant'] = z(D(w1/Us)) and z(D(w2/Us))
# T5: with a scale it runs along every flow a:b  (witness)
Ub = S**2/(V-1) + S**3
T, P, DS, DV, eps, ESP, EVP = stage(Ub)
wt = -DV(eps)*P*V/(2*Ub)
A_, B_ = sp.symbols('A_ B_')
run = sp.simplify(A_*DS(wt) + B_*DV(wt))
vals = [run.subs({S: s0, V: v0}) for s0, v0 in ((1, 2), (2, 3), (1, 3))]
out['T5_runs'] = sp.solve(vals, [A_, B_], dict=True) in ([], [{A_: 0, B_: 0}])
# T6: sign. exchanging the order of the two readings flips w; value at a point
num = sp.simplify(wt.subs({S: 1, V: 2})); out['T6_value'] = str(num); out['T6_sign'] = num != 0
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'DW1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
