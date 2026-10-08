"""SW2, wide family: the fall function free in both r and x = cos(theta):  b = K/rho,  K^2 = P(r, x).

Same frame and same formulas as sw2_turning_centre_exact.py.  Writes the contracted curvature to
sw2_wide_family.json; the main script reads that file and draws the consequences.  Takes about twenty minutes,
so it is not part of the routine test.  Symbolic (sympy).  Python 3.12."""
import json, os, time
import sympy as sp
r, x, a = sp.symbols('r x a', positive=True)
rho, s, sig, K = sp.symbols('rho s sigma K', positive=True)
J = {}
names = ['', 'r', 'x', 'rr', 'rx', 'xx', 'rrr', 'rrx', 'rxx', 'xxx']
for n in names:
    J[n] = sp.Symbol('P' + n)
P = J['']
key = lambda n: ''.join(sorted(n))
SQ = {rho: r**2 + a**2*x**2, s: r**2 + a**2, sig: 1 - x**2, K: P}
ETA = sp.diag(1, -1, -1, -1)

def red(e):
    e = sp.cancel(e)
    num, den = sp.fraction(e)
    def lower(p):
        p = sp.expand(p)
        for sym, sq in SQ.items():
            pp = sp.Poly(p, sym)
            p = sp.expand(sum(co*sq**(k//2)*sym**(k % 2) for (k,), co in pp.terms()))
        return p
    return sp.cancel(lower(num)/lower(den))

def D(f, v):
    out = sp.diff(f, r if v == 'r' else x)
    if v == 'r':
        out += sp.diff(f, rho)*r/rho + sp.diff(f, s)*r/s
    else:
        out += sp.diff(f, rho)*a**2*x/rho - sp.diff(f, sig)*x/sig
    out += sp.diff(f, K)*J[v]/(2*K)
    for n in names:
        if len(n) < 3:
            out += sp.diff(f, J[n])*J[key(n + v)]
    return out

def build(b):
    E = sp.Matrix([[1, -b*s/rho, 0, 0], [0, s/rho, 0, 0], [0, 0, -sig/rho, 0], [0, a*sig*b/rho, 0, 1/(s*sig)]])
    TH = sp.Matrix([[1, 0, 0, 0], [b, rho/s, 0, -a*sig**2*b], [0, 0, -rho/sig, 0], [0, 0, 0, s*sig]])
    ap = lambda v, f: v[1]*D(f, 'r') + v[2]*D(f, 'x')
    c = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(i+1, 4):
            comm = [ap(E.row(i), E[j, m]) - ap(E.row(j), E[i, m]) for m in range(4)]
            for k in range(4):
                c[i][j][k] = red(sum(comm[m]*TH[k, m] for m in range(4)))
                c[j][i][k] = -c[i][j][k]
    low = lambda i, j, k: ETA[k, k]*c[i][j][k]
    w = [[[red(ETA[i, i]*sp.Rational(1, 2)*(low(k, j, i) - low(j, i, k) + low(i, k, j))) for k in range(4)] for j in range(4)] for i in range(4)]
    def riem(i, j, k, l):
        v = ap(E.row(k), w[i][j][l]) - ap(E.row(l), w[i][j][k])
        v += sum(w[i][e][k]*w[e][j][l] - w[i][e][l]*w[e][j][k] for e in range(4))
        v -= sum(c[k][l][e]*w[i][j][e] for e in range(4))
        return v
    Ric = sp.Matrix(4, 4, lambda j, l: red(sum(riem(i, j, i, l) for i in range(4))))
    return E, c, w, Ric

if __name__ == '__main__':
    t0 = time.time()
    E, c, w, Ric = build(K/rho)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump({'%d%d' % (j, l): sp.srepr(Ric[j, l]) for j in range(4) for l in range(j, 4)},
              open(os.path.join(here, 'sw2_wide_family.json'), 'w'), indent=0)
    print('seconds', round(time.time() - t0))
