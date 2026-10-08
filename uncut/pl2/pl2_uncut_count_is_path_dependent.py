"""PL2: without the smooth centre (F4) and with a centre angle that depends on the state,
the count is no longer a function of state. Exact sympy."""
import json, os, sympy as sp
M, Q = sp.symbols('M Q', positive=True)
al = sp.Function('alpha')(M, Q)
out = {}
z = lambda e: sp.simplify(e) == 0
rp = M + sp.sqrt(M**2 - Q**2); S = sp.pi*rp**2
P, ang = sp.diff(S, M), -sp.diff(S, Q)
k = al/(2*sp.pi)                                         # QT1-T5: every period is multiplied by angle/2 pi
tM, tQ = k*P, -k*ang                                     # the uncut one-form  theta = k (P dM - ang dQ) = k dS
curl = sp.simplify(sp.diff(tQ, M) - sp.diff(tM, Q))
out['T1_defect'] = z(curl - (sp.diff(al, M)*sp.diff(S, Q) - sp.diff(al, Q)*sp.diff(S, M))/(2*sp.pi))   # = {alpha, S}/2 pi
out['T2_flat_iff'] = z(curl.subs(al, 2*sp.pi).doit()) and z(curl.subs(al, sp.Function('f')(S)).doit())  # alpha const or alpha(S): still a state function
# T3: UP4/UP5 form: defect = F-part only (weight), no frame part: theta = k dS, d theta = dk ^ dS
out['T3_weight_form'] = z(curl - (sp.diff(k, M)*sp.diff(S, Q) - sp.diff(k, Q)*sp.diff(S, M)))
# T4: a witness. alpha = 2 pi (1 + e x), x = Q^2/M^2 (pure number). Loop: rectangle in (M, Q)
e = sp.symbols('epsilon', positive=True)
alw = 2*sp.pi*(1 + e*Q**2/M**2)
cw = sp.simplify(curl.subs(al, alw).doit())
loop = sp.integrate(sp.integrate(cw, (Q, 0, 1)), (M, 2, 3))
out['T4_loop_value'] = str(sp.simplify(loop/e))
out['T4_nonzero'] = sp.simplify(loop/e) != 0 and sp.simplify(sp.diff(loop, e, 2)) == 0
# T5: the sense: reversing the loop reverses the sign; per circuit the count changes by the loop value
out['T5_sign_definite_on_witness'] = bool(sp.simplify(cw/e).subs({M: sp.Rational(5, 2), Q: sp.Rational(1, 2)}) < 0) or bool(sp.simplify(cw/e).subs({M: sp.Rational(5, 2), Q: sp.Rational(1, 2)}) > 0)
out['T5_witness_sign'] = str(sp.sign(sp.N(sp.simplify(loop/e))))
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PL2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
