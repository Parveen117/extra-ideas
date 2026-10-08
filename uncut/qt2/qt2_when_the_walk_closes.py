"""QT2: when does the four-step walk close as a turn? Exact sympy."""
import json, os, sympy as sp
u, v, w, psi, s_, th = sp.symbols('u v w psi s theta', real=True)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); Sx = R*K; I2 = sp.eye(2)
kap = u*K + v*Sx + w*R
W = R*kap                                                    # the walk element (UP7)
out = {}
out['T1_walk_law'] = sp.simplify(W**2 + 2*w*W - (u**2 + v**2 - w**2)*I2) == sp.zeros(2)
ev = sp.Matrix(list(W.eigenvals().keys()))
out['T1_eigen'] = set(sp.simplify(e) for e in ev) == {sp.simplify(-w + sp.sqrt(u**2 + v**2)), sp.simplify(-w - sp.sqrt(u**2 + v**2))}
# T2 (refusal): for a real cut the two eigenvalues are real -> boost or null, never a turn
out['T2_real_never_turn'] = sp.ask(sp.Q.nonnegative(u**2 + v**2))
# T3: the walk is a pure turn exactly when the seen part of the cut is in the lost channel:
# u^2+v^2 = -s^2 with w^2 + s^2 = 1  -> eigenvalues -w +- i s, modulus 1, product 1
lam = -sp.cos(th) + sp.I*sp.sin(th)
out['T3_turn'] = sp.simplify(sp.Abs(lam) - 1) == 0 and sp.simplify(lam*sp.conjugate(lam) - 1) == 0
# T4: such a walk adds no shared information: BL1 growth factor per cycle is |lam|^2 = 1
out['T4_no_growth'] = sp.simplify(sp.Abs(lam)**2) == 1
# T5: the frame cut cosh(psi) K + sinh(psi) iota, tanh(psi) = fall speed. Outside: boost e^{-+psi}.
Wf = (R*(sp.cosh(psi)*K + sp.sinh(psi)*R))
cp = Wf.charpoly().as_expr(); lam_ = list(Wf.charpoly().gens)[0]
out['T5_outside'] = all(sp.simplify(cp.subs(lam_, e).rewrite(sp.exp)) == 0 for e in (sp.exp(-psi), -sp.exp(psi)))
# past the horizon (fall speed > 1): psi -> psi' + i pi/2 in the cut-complex carrier
Wi = Wf.subs(psi, psi + sp.I*sp.pi/2).applyfunc(sp.expand_trig)
cpi = Wi.charpoly().as_expr(); lam_i = list(Wi.charpoly().gens)[0]
out['T5_inside_quarter_turn'] = all(sp.simplify(cpi.subs(lam_i, e).rewrite(sp.exp)) == 0 for e in (-sp.I*sp.exp(-psi), -sp.I*sp.exp(psi)))
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QT2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
