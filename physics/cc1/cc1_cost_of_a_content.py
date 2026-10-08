"""CC1: the cost of carrying a content, and its ratios.

MG1: the mass count itself is open; ratios in which the cell and the coupling cancel are not.
The free rate of a content in the Yang-Mills line is its Casimir number (YM-22, YM-40, YM95, YM96).
Here the Casimir numbers are computed exactly from the completeness form on tensor powers
(integer matrices, exact ranks), and their ratios are set against a published simulation.
Python 3.12."""
import itertools, json, os
from fractions import Fraction as Fr
import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

HERE = os.path.dirname(os.path.abspath(__file__))


def scaled_casimir(N, p, q):
    """2N x Casimir on (fundamental)^p (x) (dual)^q, as an integer matrix.
    C = (1/2) [ sum_ij E_ij E_ji - (1/N) (sum_i E_ii)^2 ]: fundamental (N^2-1)/(2N) (YM96)."""
    slots = p + q
    dim = N**slots
    def unit(i, j):
        e = np.zeros((N, N), dtype=np.int64); e[i, j] = 1; return e
    def act(i, j):
        tot = np.zeros((dim, dim), dtype=np.int64)
        for s in range(slots):
            mats = [np.eye(N, dtype=np.int64)]*slots
            mats = list(mats)
            mats[s] = unit(i, j) if s < p else -unit(j, i)
            M = mats[0]
            for t in mats[1:]:
                M = np.kron(M, t)
            tot += M
        return tot
    E = {(i, j): act(i, j) for i in range(N) for j in range(N)}
    S = sum(E[i, j] @ E[j, i] for i in range(N) for j in range(N))
    T = sum(E[i, i] for i in range(N))
    return N*S - T @ T, dim


def nullity(M, val):
    A = M - val*np.eye(M.shape[0], dtype=np.int64)
    return M.shape[0] - DomainMatrix(A.tolist(), A.shape, ZZ).rank()


def c3(p, q):                       # SU(3)
    return Fr(p*p + q*q + p*q, 3) + p + q


def dim3(p, q):
    return (p + 1)*(q + 1)*(p + q + 2)//2


# continuum-extrapolated ratios V_D / V_F, arXiv:hep-lat/0006022 Tables X and XI, rows with r < 2 r0 (the range of its 5 % statement)
SIM = {
    '8':   [(2.24, .09), (2.27, .03), (2.21, .04), (2.23, .04), (2.24, .02), (2.23, .03), (2.24, .03), (2.23, .04), (2.22, .04), (2.24, .04), (2.27, .05), (2.28, .04), (2.24, .04), (2.24, .05), (2.29, .05), (2.28, .06), (2.33, .06), (2.26, .10), (2.27, .09)],
    '6':   [(2.48, .11), (2.53, .04), (2.50, .04), (2.45, .05), (2.50, .03), (2.45, .04), (2.45, .05), (2.50, .05), (2.49, .04), (2.48, .05), (2.49, .07), (2.43, .07), (2.55, .05), (2.49, .08), (2.60, .08), (2.61, .10), (2.55, .09), (2.67, .12), (2.66, .11)],
    '15a': [(3.97, .29), (3.99, .10), (3.86, .15), (3.96, .11), (3.97, .08), (3.95, .11), (4.00, .08), (3.87, .12), (3.85, .11), (3.88, .12), (4.11, .16), (4.05, .16), (3.86, .13), (4.11, .19), (4.09, .13), (3.84, .17), (4.00, .14), (4.10, .19), (4.10, .18)],
    '10':  [(4.45, .37), (4.39, .18), (4.38, .16), (4.37, .13), (4.45, .11), (4.43, .12), (4.43, .10), (4.34, .14), (4.37, .14), (4.48, .10), (4.46, .11), (4.50, .10), (4.44, .10), (4.55, .12), (4.48, .13), (4.50, .16), (4.79, .14), (4.74, .20), (4.71, .19)],
    '27':  [(6.23, .65), (6.10, .33), (5.96, .29), (6.10, .16), (6.21, .15), (6.18, .16), (6.17, .16), (6.28, .14), (6.21, .12), (6.06, .17), (5.76, .24), (6.00, .21), (6.21, .16), (6.50, .23), (6.02, .21), (6.39, .31), (6.53, .27), (6.51, .33), (6.54, .32)],
    '24':  [(5.98, .57), (5.94, .22), (5.82, .23), (5.86, .14), (5.97, .14), (5.92, .14), (5.95, .14), (5.95, .15), (5.91, .14), (5.82, .19), (5.52, .24), (5.99, .23), (5.96, .18), (6.52, .23), (6.05, .29), (6.19, .28), (6.09, .22), (6.21, .30), (6.16, .29)],
    '15s': [(6.84, .83), (6.95, .27), (6.87, .23), (6.86, .19), (6.90, .17), (6.89, .17), (6.81, .21), (7.00, .13), (6.95, .12), (6.71, .19), (6.99, .16), (7.08, .14), (7.19, .12), (7.27, .15), (7.23, .14), (7.33, .17), (7.51, .16), (7.16, .48), (7.33, .35)],
}
REP = {'8': (1, 1), '6': (2, 0), '15a': (2, 1), '10': (3, 0), '27': (2, 2), '24': (3, 1), '15s': (4, 0)}


def run():
    out, num = {}, {}

    # C1: three colours.  On each tensor power the scaled Casimir is diagonal with the values c(p', q'); the top value has the dimension of (p, q)
    ok_split, ok_top, cas = True, True, {}
    for (p, q) in [(1, 0), (2, 0), (1, 1), (3, 0), (2, 1), (4, 0), (3, 1), (2, 2)]:
        M, dim = scaled_casimir(3, p, q)
        vals = sorted({6*c3(a, b) for a in range(p + q + 1) for b in range(p + q + 1) if a + b <= p + q})
        nul = {int(v): nullity(M, int(v)) for v in vals}
        ok_split &= sum(nul.values()) == dim
        top = int(6*c3(p, q))
        ok_top &= nul[top] == dim3(p, q) and all(n_ == 0 for v, n_ in nul.items() if v > top)
        cas[(p, q)] = c3(p, q)
    out['C1_casimir_splits_every_tensor_power'] = ok_split
    out['C1_top_value_has_the_dimension_of_the_content'] = ok_top
    out['C1_fundamental_is_four_thirds_adjoint_three'] = cas[(1, 0)] == Fr(4, 3) and cas[(1, 1)] == 3

    # C2: two colours: contents 1/2, 1, 3/2, 2 have j(j+1); YM-40's pinch 5/4 is C_1 - C_1/2
    ok2 = True
    for n_ in (1, 2, 3, 4):
        M, dim = scaled_casimir(2, n_, 0)
        ok2 &= nullity(M, n_*(n_ + 2)) == n_ + 1 and sum(nullity(M, k*(k + 2)) for k in range(n_ % 2, n_ + 1, 2)) == dim
    out['C2_two_colours_j_j_plus_1'] = ok2
    out['C2_pinch_five_quarters'] = Fr(2) - Fr(3, 4) == Fr(5, 4) and Fr(2)/Fr(3, 4) == Fr(8, 3)

    # C3: the ratios
    ratio = {k: cas[v]/cas[(1, 0)] for k, v in REP.items()}
    expect = {'8': Fr(9, 4), '6': Fr(5, 2), '15a': Fr(4), '10': Fr(9, 2), '27': Fr(6), '24': Fr(25, 4), '15s': Fr(7)}
    out['C3_ratios'] = ratio == expect
    num['ratios'] = {k: str(v) for k, v in ratio.items()}

    # C4: against the simulation
    tab, worst_mean, worst_point, npts, inside2 = {}, 0.0, 0.0, 0, 0
    for k, rows in SIM.items():
        w = [1/e**2 for _, e in rows]
        mean = sum(v*wi for (v, _), wi in zip(rows, w))/sum(w)
        err = sum(w)**-0.5
        e0 = float(ratio[k])
        tab[k] = dict(ratio=e0, simulation_mean=round(mean, 3), sigma_uncorrelated=round(err, 3), relative=round(mean/e0 - 1, 4))
        worst_mean = max(worst_mean, abs(mean/e0 - 1))
        for v, e in rows:
            npts += 1; inside2 += abs(v - e0) <= 2*e
            worst_point = max(worst_point, abs(v/e0 - 1))
    num['against_simulation'] = tab
    num['points'] = npts; num['points_within_two_sigma'] = int(inside2)
    num['worst_mean_deviation'] = worst_mean; num['worst_single_point_deviation'] = worst_point
    out['C4_every_mean_within_five_percent'] = worst_mean < 0.05
    out['C4_most_points_within_two_sigma'] = inside2 >= 0.9*npts

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'CC1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k, v in o.items(): print(k, v)
    for k, v in n.items(): print(k, v)
