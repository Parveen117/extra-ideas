"""ON1: where non-commutation comes from. Exact sympy; log coordinates s=Log S, v=Log V, p=Log P(s,v)."""
import json, sympy as sp
s, v = sp.symbols('s v', real=True)
p = sp.Function('p')(s, v)
f = sp.Function('f')(s, v)
ps, pv = sp.diff(p, s), sp.diff(p, v)
# scale readings: E[S|V]=d_s, E[V|S]=d_v, E[S|P]=d_s-(ps/pv)d_v, E[V|P]=d_v-(pv/ps)d_s
ESV = lambda g: sp.diff(g, s)
EVS = lambda g: sp.diff(g, v)
ESP = lambda g: sp.diff(g, s) - ps/pv*sp.diff(g, v)
EVP = lambda g: sp.diff(g, v) - pv/ps*sp.diff(g, s)
eps = ps/pv                       # = -(dLogV/dLogS)_P : a pure ratio of responses
comm = lambda A, B, g: sp.simplify(A(B(g)) - B(A(g)))
out = {}
# T1: the defect is the scale-derivative (running) of the pure ratio
out['T1a'] = sp.simplify(comm(ESV, ESP, f) + ESV(eps)*EVS(f)) == 0
out['T1b'] = sp.simplify(comm(EVS, EVP, f) + EVS(1/eps)*ESV(f)) == 0
out['T1c'] = comm(ESV, EVS, f) == 0 and sp.simplify(ESP(p)) == 0
# T2: all defects vanish  <=>  eps constant  <=>  p = F(c*s+v) (one scaling variable)
c = sp.symbols('c', positive=True); F = sp.Function('F')
sub = lambda e: sp.simplify(e.subs(p, F(c*s+v)).doit())
out['T2_forward'] = sub(eps) == c and sub(ESV(eps)) == 0 and sub(EVS(1/eps)) == 0
# converse: d_s eps = 0 and d_v eps = 0 give ps - c pv = 0, whose general solution is F(c s + v)
x_, y_ = sp.symbols('x_ y_'); q = sp.Function('q')
pq = q(c*s+v, s)                                   # any p, written in variables (c s+v, s)
resid = sp.simplify(sp.diff(pq, s) - c*sp.diff(pq, v) - sp.diff(q(x_, y_), y_).subs({x_: c*s+v, y_: s}))
out['T2_converse'] = sp.simplify(resid.doit()) == 0   # ps - c pv = q_y : zero iff p depends on c s+v alone
# T3: examples. power law commutes; a potential with a length scale does not
ex_flat = sp.simplify(ESV(eps).subs(p, 3*s-2*v).doit())
b = sp.symbols('b', positive=True)
pvdw = sp.log(sp.exp(s)/(sp.exp(v)-b))            # P = S/(V-b): excluded volume b is a scale
ex_b = sp.simplify(EVS(1/eps).subs(p, pvdw).doit())
out['T3_flat'] = ex_flat == 0
out['T3_scale_value'] = str(ex_b)
out['T3_scale'] = ex_b != 0 and sp.limit(ex_b, b, 0) == 0
# T4: the running is the same object as the level-rule ratio: along P the pair (S,V) has degree
# k = 1 - eps ... flow equation d eps/d s = beta(eps, s, v): beta is the commutator coefficient
out['T4_beta'] = sp.simplify(comm(ESV, ESP, v) + ESV(eps)) == 0   # read on f=v: defect = -beta itself
out['pass'] = all(val for k_, val in out.items() if isinstance(val, bool))
json.dump({k_: (bool(x) if isinstance(x, (bool, sp.logic.boolalg.BooleanAtom)) else x) for k_, x in out.items()},
          open(__file__.rsplit('/',1)[0] + '/ON1_RESULT.json', 'w'), indent=1)
if __name__ == '__main__':
    print(out)
