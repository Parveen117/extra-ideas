"""CP1: a self-similar centre. U(l^a S, l^b V) = l U. Exact sympy."""
import json, os, sympy as sp
S, V, a, b = sp.symbols('S V a b', positive=True)
G = sp.Function('G')
U = S**(1/a)*G(V/S**(b/a))                 # general solution of a S U_S + b V U_V = U
DS = lambda g: S*sp.diff(g, S); DV = lambda g: V*sp.diff(g, V)
D = lambda g: a*DS(g) + b*DV(g)            # the scaling flow
T, P = sp.diff(U, S), -sp.diff(U, V)
z = lambda e: sp.simplify(e.doit()) == 0
out = {}
out['T0_weights'] = z(D(U)-U) and z(D(T)-(1-a)*T) and z(D(P)-(1-b)*P)
eps = DS(P)/DV(P)                           # response ratio of NC1
bS, bV = DS(eps), DV(eps)                   # the two defects of NC1
out['T1_one_defect'] = z(a*bS + b*bV)       # defects tied: bS : bV = -b : a
out['T2_weights_read'] = z(a*DS(P) + b*DV(P) - (1-b)*P)   # a*eps_S + b*eps_V = 1-b on Log P
# T3: the defect does not vanish in general (so criticality is NOT commutation, it is one defect)
ex = sp.simplify(bV.subs(G, sp.Lambda(S, 1/(1+S)+S**2)).doit().subs({a: 2, b: 3, S: 1, V: 2}))
out['T3_defect_value'] = str(ex); out['T3_nonzero'] = ex != 0
# T4: without self-similarity the two defects are independent
P2 = S/(V-1) + S**2; e2 = DS(P2)/DV(P2)
r = sp.simplify(DS(e2)/DV(e2))
out['T4_independent'] = sp.simplify(sp.diff(r, S)) != 0   # ratio not a constant -> no weights exist
# T5: same flow acts on the Maxwell 2-form readings: T and P ratio T S/(P V) is flow-invariant
out['T5_invariant'] = z(D(T*S/(P*V)))
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CP1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
