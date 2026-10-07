"""MC1: 1/r from minimum cost - the least total variance between neighbouring shells.

Shells at radii r_k = b^k (equal steps of the scale s = log r), d space dimensions.
A reading f_k on shell k.  Cost = sum over neighbouring shells of (cells on the shell) x
(squared difference per unit distance) x (thickness):
        C[f] = sum_k w_k (f_{k+1} - f_k)^2 ,      w_k = r_k^(d-2)   (up to one common factor).
T1  the minimizer with fixed ends has a constant flux Q = w_k (f_{k+1} - f_k) and is exactly
        f_k = A + B r_k^-(d-2)        (d != 2);      f_k = A + B k   (d = 2).
T2  it is the unique minimum: cost(f + h) - cost(f) = sum_k w_k (h_{k+1} - h_k)^2 > 0.
T3  gradient x area is the same on every shell: (f_{k+1} - f_k)/(r_{k+1} - r_k) x r_k^(d-1) = constant.
T4  d = 3 with f = memory: m = r_s / r, and GR1's clock factor at every node.
T5  the choice of reading matters: least-cost N instead of least-cost memory differs at second order.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json


def weights(b, d, n):
    return [b**((d-2)*k) for k in range(n)]                  # w_k for the link k -> k+1


def cost(f, w):
    return sum(wk*(f[k+1]-f[k])**2 for k, wk in enumerate(w))


def minimizer(w, f0, fn):
    """Solve the stationarity equations w_{k-1}(f_k - f_{k-1}) = w_k (f_{k+1} - f_k) exactly."""
    n = len(w)
    # constant flux Q: f_{k+1} - f_k = Q / w_k, sum = fn - f0
    total = sum(1/wk for wk in w)
    Q = (fn-f0)/total
    f = [f0]
    for wk in w:
        f.append(f[-1]+Q/wk)
    return f, Q


def tridiagonal_check(f, w):
    for k in range(1, len(f)-1):
        if w[k-1]*(f[k]-f[k-1]) != w[k]*(f[k+1]-f[k]):
            raise ValueError('stationarity failed')


def shell_control(b, d, n, f0, fn):
    w = weights(b, d, n)
    f, Q = minimizer(w, f0, fn)
    tridiagonal_check(f, w)
    r = [b**k for k in range(n+1)]
    if d != 2:
        # closed form A + B r^-(d-2) through the two ends
        e = d-2
        B = (f0-fn)/(1-r[n]**(-e))
        A = f0-B
        if any(f[k] != A+B*r[k]**(-e) for k in range(n+1)):
            raise ValueError('minimizer is not the power law r^-(d-2)')
    else:
        if any(f[k+1]-f[k] != f[1]-f[0] for k in range(n)):
            raise ValueError('d = 2 minimizer is not linear in the scale')
    # T2: every perturbation vanishing at the ends raises the cost by exactly its own cost
    base = cost(f, w)
    for h in ([F(0)]+[F((-1)**k, k+1) for k in range(1, n)]+[F(0)],
              [F(0)]+[F(1, 7)]*(n-1)+[F(0)]):
        g = [a+c for a, c in zip(f, h)]
        if cost(g, w)-base != cost(h, w) or cost(h, w) <= 0:
            raise ValueError('second variation failed')
    # T3: gradient times area constant
    vals = {(f[k+1]-f[k])/(r[k+1]-r[k])*r[k]**(d-1) for k in range(n)}
    if len(vals) != 1:
        raise ValueError('gradient times area is not constant')
    return dict(dimension=d, ratio=str(b), shells=n, flux=str(Q), cost=str(base),
                exponent=-(d-2), profile=[str(x) for x in f[:5]])


def gravity_control():
    """d = 3, reading = memory, ends m(r_0) = m0 and m(infinity) -> 0 approached on a long chain."""
    b, m0 = F(2), F(16, 25)                                   # r_s / r_0 = 16/25, N_0 = 3/5
    rows = []
    for n in (4, 8, 16, 32):
        w = weights(b, 3, n)
        target = m0/b**n                                      # the exact 1/r value at the outer end
        f, Q = minimizer(w, m0, target)
        if any(f[k] != m0/b**k for k in range(n+1)):
            raise ValueError('least-cost memory is not r_s / r')
        rows.append(dict(shells=n, flux=str(Q)))
    if len({r['flux'] for r in rows}) != 1:
        raise ValueError('the flux must not depend on where the chain is cut')
    # GR1 at the nodes where N is rational: m = 16/25 -> N = 3/5; m = 9/25 needs another start
    nodes = []
    for m, N in ((F(16, 25), F(3, 5)), (F(9, 25), F(4, 5)), (F(25, 169), F(12, 13))):
        if N*N != 1-m:
            raise ValueError('clock factor failed')
        nodes.append(dict(memory=str(m), N=str(N), share=str((1+N)/2)))
    # T5: least-cost N (harmonic N with the same far field) against least-cost memory
    m = F(9, 25)
    N_memory = F(4, 5)                                        # sqrt(1 - m)
    N_harmonic = 1-m/2                                        # N = 1 - r_s/(2r)
    if N_harmonic == N_memory or N_harmonic**2-(1-m) != m*m/4:
        raise ValueError('second-order difference is not (r_s/2r)^2')
    return dict(chains=rows, flux_is='-(b - 1)/b times m0 per unit of the common weight factor', nodes=nodes,
                least_cost_N_vs_memory=dict(memory=str(m), N_from_memory=str(N_memory),
                                            N_from_harmonic_N=str(N_harmonic), difference_in_N_squared=str(m*m/4)))


def tail_control():
    """QC5's logistic share in four dimensions: the deficit 1 - p falls as r^-2 = r^-(d-2)."""
    lam2 = F(4, 9)
    vals = []
    for r2 in (F(1), F(100), F(10000), F(10**6)):
        p = r2/(r2+lam2)
        vals.append((1-p)*r2/lam2)
    if not all(a < b < 1 for a, b in zip(vals, vals[1:])):
        raise ValueError('deficit times r^2 must increase to lambda^2')
    return dict(dimension=4, exponent=-2, last_ratio=str(vals[-1]))


def run():
    return dict(shells=[shell_control(F(2), d, 6, F(1), F(0)) for d in (1, 2, 3, 4, 5)]
                + [shell_control(F(3, 2), 3, 7, F(5), F(2)), shell_control(F(3), 4, 5, F(0), F(1))],
                gravity=gravity_control(), tail=tail_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('MC1_RESULT.json', 'w'), indent=1)
    for r in res['shells'][:5]:
        print({k: r[k] for k in ('dimension', 'exponent', 'flux', 'profile')})
    print(res['gravity']['chains'], res['gravity']['least_cost_N_vs_memory'])
    print(res['tail'])
