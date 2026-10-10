"""YC20: exact finite controls for metric-aware full hidden returns.

The full-carrier, domain and spectral hypotheses are in the companion note.
These rational/complex matrix controls are not a finite approximation of YM.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    'yc20_y19', ROOT / 'physics/yc19/yc19_local_excitation_operator.py')
y19 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = y19
spec.loader.exec_module(y19)

DRAFT_SOURCES = {
    'Transforma — Projective Upgrade (λ‑eigengeometry).pdf':
        'ed31fd0908451a8e83f38d9ed5f7bc9187e305189acdea4817e74b2fa8aa2795',
    '1Monti Operator & Winding Number — Cheat Sheet.pdf':
        'bd0343d2676da7fcc5d1a911f8499001e72ff809289947e88f14cc19b424f848',
    'ilovepdf_merged.pdf':
        '16984f77e8ed150e141425c55916ea4227804301069df4fc55855df15a1c7d74',
}


def simplify(M):
    return M.applyfunc(sp.simplify)


def equal(A, B):
    return A.shape == B.shape and simplify(A-B) == sp.zeros(*A.shape)


def require_hermitian(M):
    if M.rows != M.cols or not equal(M, M.H):
        raise ValueError('a Hermitian square matrix is required')
    if M.atoms(sp.Float):
        raise ValueError('exact finite controls do not accept floating entries')


def positive_definite(M):
    require_hermitian(M)
    return all(sp.simplify(M[:k, :k].det()).is_positive is True
               for k in range(1, M.rows+1))


def positive_semidefinite(M):
    """All principal minors: deliberately small finite certificate controls."""
    require_hermitian(M)
    return all(sp.simplify(M.extract(I, I).det()).is_nonnegative is True
               for n in range(1, M.rows+1)
               for I in combinations(range(M.rows), n))


def schur(M, retained):
    """Eliminate the entire finite complement, returning its harmonic lift."""
    if M.rows != M.cols or type(retained) is not int or not 0 < retained < M.rows:
        raise ValueError('a nontrivial retained cut of a square matrix is required')
    A, B = M[:retained, :retained], M[:retained, retained:]
    C, D = M[retained:, :retained], M[retained:, retained:]
    if sp.simplify(D.det()) == 0:
        raise ValueError('the hidden block must be invertible')
    Z = -D.inv()*C
    W = sp.eye(retained).col_join(Z)
    return simplify(A+B*Z), simplify(W)


def metric_frame(K, G, retained):
    require_hermitian(K)
    require_hermitian(G)
    if K.shape != G.shape or not positive_definite(G):
        raise ValueError('energy and positive metric must share the carrier')
    if type(retained) is not int or not 0 < retained < G.rows:
        raise ValueError('a nontrivial retained cut is required')
    p, q = retained, G.rows-retained
    Gp, N, Gh = G[:p, :p], G[:p, p:], G[p:, p:]
    Kp, J, Kh = K[:p, :p], K[:p, p:], K[p:, p:]
    X = Gh.inv()*N.H
    T = sp.eye(p).row_join(sp.zeros(p, q)).col_join(
        (-X).row_join(sp.eye(q)))
    g = simplify(Gp-N*X)
    k = simplify(Kp-J*X-X.H*J.H+X.H*Kh*X)
    b = simplify(J-X.H*Kh)
    return dict(K=K, G=G, retained=p, Gp=Gp, N=N, Gh=Gh, Kh=Kh,
                X=X, T=T, g=g, k=k, b=b)


def retained_return(frame, z):
    z = sp.sympify(z)
    if z.is_real is not True:
        raise ValueError('the energy-step identities use real energy')
    p, Gh, Kh = frame['retained'], frame['Gh'], frame['Kh']
    g, k, b, X = (frame[key] for key in ('g', 'k', 'b', 'X'))
    D = Kh-z*Gh
    if sp.simplify(D.det()) == 0:
        raise ValueError('energy hits a hidden pole')
    V = D.inv()*b.H
    W = sp.eye(p).col_join(-X-V)
    F = k-z*g-b*V
    slope = g+V.H*Gh*V
    second = -2*V.H*Gh*D.inv()*Gh*V
    return {key: simplify(value) for key, value in
            dict(D=D, F=F, W=W, slope=slope, second=second).items()}


def compatible_cut(frame):
    """Transport the diagonal two-channel cut into the original metric frame."""
    p, n = frame['retained'], frame['G'].rows
    J0 = sp.diag(sp.eye(p), -sp.eye(n-p))
    Kcut = simplify(frame['T']*J0*frame['T'].inv())
    return dict(cut=Kcut, plus=(sp.eye(n)+Kcut)/2, minus=(sp.eye(n)-Kcut)/2)


def residual_enclosure(frame, z, hidden_floor, Y):
    z, d = sp.sympify(z), sp.sympify(hidden_floor)
    if (d-z).is_positive is not True:
        raise ValueError('energy must be strictly below the proved hidden floor')
    Gh, Kh, b = frame['Gh'], frame['Kh'], frame['b']
    if not positive_semidefinite(Kh-d*Gh):
        raise ValueError('the proposed full hidden floor is false')
    if Y.shape != b.H.shape or Y.atoms(sp.Float):
        raise ValueError('exact trial columns must have the complete hidden carrier')
    row = retained_return(frame, z)
    D = row['D']
    R = b.H-D*Y
    trial = b*Y+Y.H*b.H-Y.H*D*Y
    error = R.H*Gh.inv()*R/(d-z)
    upper = frame['k']-z*frame['g']-trial
    return {key: simplify(value) for key, value in
            dict(residual=R, trial=trial, error=error, lower=upper-error,
                 upper=upper, exact=row['F'], exact_error=R.H*D.inv()*R).items()}


def inertia(M):
    """Exact congruence count (negative, zero, positive), with 2x2 pivots."""
    require_hermitian(M)
    n = M.rows
    if n == 0:
        return (0, 0, 0)
    for i in range(n):
        a = sp.simplify(M[i, i])
        if a != 0:
            if a.is_positive is not True and a.is_negative is not True:
                raise ValueError('pivot sign is not exactly decidable')
            order = [i]+[j for j in range(n) if j != i]
            N = M.extract(order, order)
            rest = simplify(N[1:, 1:]-N[1:, :1]*N[:1, 1:]/a)
            count = inertia(rest)
            extra = (1, 0, 0) if a.is_negative else (0, 0, 1)
            return tuple(x+y for x, y in zip(count, extra))
    for i, j in combinations(range(n), 2):
        if sp.simplify(M[i, j]) != 0:
            order = [i, j]+[k for k in range(n) if k not in (i, j)]
            N = M.extract(order, order)
            rest = simplify(N[2:, 2:]-N[2:, :2]*N[:2, :2].inv()*N[:2, 2:])
            count = inertia(rest)
            return (count[0]+1, count[1], count[2]+1)
    return (0, n, 0)


def iteration_reserve(squared_correlations):
    values = [Q(c) for c in squared_correlations]
    if any(c < 0 or c >= 1 for c in values):
        raise ValueError('each squared correlation must be in [0,1)')
    product = Q(1)
    for c in values:
        product *= 1-c
    return dict(product=product, sum=sum(values, Q(0)),
                additive_lower=max(Q(0), 1-sum(values, Q(0))))


def complex_fixture():
    S = sp.Matrix([[1, sp.I/3, sp.Rational(1,5)],
                   [0, 1, (1+sp.I)/4], [0, 0, 1]])
    H = sp.Matrix([[2, sp.Rational(1,3), sp.I/5],
                   [sp.Rational(1,3), 5, (1-sp.I)/7],
                   [-sp.I/5, (1+sp.I)/7, 9]])
    return simplify(S.H*H*S), simplify(S.H*S)


def nested_fixture():
    S = sp.Matrix([[1, sp.Rational(1,3), 0, sp.Rational(1,5)],
                   [0, 1, sp.Rational(1,4), 0],
                   [0, 0, 1, sp.Rational(1,7)], [0, 0, 0, 1]])
    G = S.T*S
    K = S.T*sp.diag(2, 4, 7, 11)*S
    return K, G


@lru_cache(maxsize=1)
def projector_control():
    theta, phi = sp.symbols('theta phi', real=True)
    cuts = [sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-sp.I],[sp.I,0]]), sp.diag(1,-1)]
    n = [sp.sin(theta)*sp.cos(phi), sp.sin(theta)*sp.sin(phi), sp.cos(theta)]
    N = sum((x*C for x, C in zip(n, cuts)), sp.zeros(2))
    P = (sp.eye(2)+N)/2
    curvature = sp.trigsimp(-sp.I*sp.trace(
        P*(P.diff(theta)*P.diff(phi)-P.diff(phi)*P.diff(theta))))
    return theta, phi, N, P, curvature


def source_pins():
    paths = ['uncut/ca1/CA1_TVSP_COMPASS_CALIBRATION.md',
             'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc19/YC19_LOCAL_EXCITATION_OPERATOR.md',
             'physics/yc19/yc19_local_excitation_operator.py',
             'physics/yc19/YC19_RESULT.json',
             'physics/yc20/YC20_METRIC_RETAINED_RETURN.md',
             'physics/yc20/yc20_metric_retained_return.py',
             'physics/yc20/test_yc20.py']
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    y19.y15.verify_predecessor('physics/yc19/YC19_RESULT.json')
    checks = {}
    control = y19.exact_metric_control()
    A, G = control['A'], control['quotient']
    K = G*A
    frame = metric_frame(K, G, 1)
    T, g, Gh, k, b, Kh = (frame[key] for key in ('T','g','Gh','k','b','Kh'))
    cut = compatible_cut(frame)
    Kcut, Pplus, Pminus = (cut[key] for key in ('cut','plus','minus'))
    checks['primitive represented cut is an involution self-adjoint in the carried metric'] = (
        equal(Kcut*Kcut, sp.eye(3)) and equal(Kcut.H*G,G*Kcut))
    checks['cut generates complementary metric-orthogonal channels'] = (
        equal(Pplus+Pminus,sp.eye(3)) and equal(Pplus*Pminus,sp.zeros(3))
        and equal(Pplus.H*G*Pminus,sp.zeros(3)))
    source = sp.Matrix([2,1,-3])
    plus, minus = Pplus*source, Pminus*source
    checks['cut source weights recover the total and signed reading exactly'] = (
        equal(plus.H*G*plus+minus.H*G*minus,source.H*G*source)
        and equal(plus.H*G*plus-minus.H*G*minus,source.H*G*Kcut*source))
    checks['Schur metric is the retained cut channel pairing'] = equal(Pplus.H*G*Pplus,sp.diag(g,sp.zeros(2)))
    balanced = sp.ones(2,1)
    native_cut = sp.diag(1,-1)
    checks['balanced centre reading does not annihilate the cut or force equilibrium'] = (
        (balanced.T*native_cut*balanced)[0] == 0
        and native_cut*balanced != sp.zeros(2,1)
        and (-(balanced.T*(sp.diag(1,2)*native_cut+native_cut*sp.diag(1,2))*balanced))[0] == 2)
    checks['YC19 quotient energy is Hermitian although its old-coordinate operator is not'] = (
        equal(K, K.H) and not equal(A, A.H))
    checks['metric reset retains its exact Schur complement'] = equal(T.H*G*T, sp.diag(g, Gh))
    checks['metric reset changes both kinetic block and hidden source'] = (
        equal(T.H*K*T, k.row_join(b).col_join(b.H.row_join(Kh)))
        and not equal(b, K[:1,1:])
        and equal(Pminus*A*Pplus, Pminus*(A*Kcut-Kcut*A)*Pplus/2)
        and equal((T.inv()*A*T)[1:,:1], Gh.inv()*b.H))
    z = sp.Rational(1,2)
    row = retained_return(frame, z)
    direct, direct_lift = schur(K-z*G, 1)
    checks['corrected return equals direct generalized-pencil elimination'] = equal(row['F'], direct)
    checks['complete hidden lift solves the exact generalized equation'] = (
        equal(row['W'], direct_lift)
        and equal((K-z*G)*row['W'], row['F'].col_join(sp.zeros(2,1))))
    checks['returned slope is the physical norm of the complete lifted reading'] = equal(
        row['slope'], row['W'].H*G*row['W'])
    checks['slope contains a positive hidden addition and curvature in energy is concave'] = (
        positive_semidefinite(row['slope']-g) and positive_semidefinite(-row['second']))
    x = sp.Rational(1,3)
    other = retained_return(frame, x)
    checks['two-energy step is exact without commuting metric and energy'] = equal(
        row['F']-other['F'], -(z-x)*row['W'].H*G*other['W'])
    checks['the complete finite hidden floor one is valid'] = positive_semidefinite(Kh-Gh)
    Y = Kh.inv()*b.H  # Deliberately reuse the zero-energy solve at z=1/2.
    enclosure = residual_enclosure(frame, z, 1, Y)
    checks['trial hidden solve leaves a nonzero fully retained residual'] = enclosure['residual'] != sp.zeros(2,1)
    checks['completed-square return error is exact'] = equal(
        b*row['D'].inv()*b.H, enclosure['trial']+enclosure['exact_error'])
    checks['full residual Gram bounds the hidden error in the correct metric'] = positive_semidefinite(
        enclosure['error']-enclosure['exact_error'])
    checks['retained lower and upper bounds enclose the complete energy pencil'] = (
        positive_semidefinite(enclosure['exact']-enclosure['lower'])
        and positive_semidefinite(enclosure['upper']-enclosure['exact']))
    checks['residual certificate excludes all excitations through one half in the finite control'] = positive_definite(enclosure['lower'])
    checks['independent full inertia confirms the retained excitation count'] = (
        inertia(K-z*G)[0] == inertia(row['F'])[0] == 0)
    checks['generalized determinant factors through the complete hidden block'] = (
        sp.simplify((K-z*G).det()-row['D'].det()*row['F'].det()) == 0)

    Kc, Gc = complex_fixture()
    fc = metric_frame(Kc, Gc, 1)
    rc = retained_return(fc, sp.Rational(1,4))
    checks['complex frame requires conjugate transpose and preserves the exact norm'] = (
        equal(rc['slope'], rc['W'].H*Gc*rc['W']) and not equal(Gc, Gc.T))
    Bp, Bq = sp.Matrix([[2]]), sp.Matrix([[1,sp.I/5],[0,2]])
    Z = sp.diag(Bp, Bq)
    changed = metric_frame(Z.H*Kc*Z, Z.H*Gc*Z, 1)
    changed_row = retained_return(changed, sp.Rational(1,4))
    checks['separate invertible frame changes transport the retained metric and pencil'] = (
        equal(changed['g'], Bp.H*fc['g']*Bp)
        and equal(changed_row['F'], Bp.H*rc['F']*Bp))
    checks['projective one-channel metric reserve survives changes of units and hidden basis'] = (
        sp.simplify(changed['g'][0,0]/changed['Gp'][0,0]-fc['g'][0,0]/fc['Gp'][0,0]) == 0)

    Kn, Gn = nested_fixture()
    L = Kn-x*Gn
    F1, W1 = schur(L, 2)
    F2, W2 = schur(F1, 1)
    Ftotal, Wtotal = schur(L, 1)
    G1, _ = schur(Gn, 2)
    G2, _ = schur(G1, 1)
    Gtotal, _ = schur(Gn, 1)
    checks['nested energy elimination is associative'] = equal(F2, Ftotal)
    checks['nested metric elimination is associative'] = equal(G2, Gtotal)
    checks['nested lifts retain every hidden component'] = equal(W1*W2, Wtotal)
    checks['energy-dependent derivative metric composes through the whole nested lift'] = equal(
        W2.H*(W1.H*Gn*W1)*W2, Wtotal.H*Gn*Wtotal)
    checks['energy-return metric is not generally the metric-only Schur complement'] = not equal(
        Wtotal.H*Gn*Wtotal, Gtotal)
    reserve = iteration_reserve([Q(1,16), Q(1,32), Q(1,64)])
    checks['conditional iteration product has the stated additive floor'] = (
        reserve['product'] >= reserve['additive_lower'] == Q(57,64))

    shear = sp.Matrix([[1,2],[0,1]])
    Hs = sp.diag(1,4)
    As, Gs = shear.inv()*Hs*shear, shear.T*shear
    checks['discarding the metric can manufacture a negative direction'] = (
        inertia((As+As.T)/2)[0] == 1 and As.charpoly().as_expr() == Hs.charpoly().as_expr())
    short, _ = schur(Gs, 1)
    checks['metric reserve can shrink while the physical spectrum is unchanged'] = short[0,0]/Gs[0,0] == sp.Rational(1,5)
    theta, phi, Np, P, curvature = projector_control()
    checks['three-cut projector is idempotent and has nonzero frame curvature'] = (
        equal(Np*Np, sp.eye(2)) and equal(P*P, P)
        and sp.trigsimp(curvature-sp.sin(theta)/2) == 0)
    epsilon, aa = sp.symbols('epsilon aa', positive=True)
    Hp = aa*sp.eye(2)+epsilon*Np
    checks['frame curvature does not determine spectral scale'] = (
        equal(Hp*P, (aa+epsilon)*P)
        and sp.diff(curvature, epsilon) == 0
        and sp.limit(2*epsilon, epsilon, 0, dir='+') == 0)
    t = sp.symbols('t', real=True)
    checks['integer winding is erased by its endpoint exponential'] = all(
        sp.exp(2*sp.pi*sp.I*n) == 1 for n in (-2,-1,0,1,2))
    xx, yy = sp.symbols('xx yy', positive=True)
    f = xx**2+xx*yy+yy**2
    checks['a smooth scalar logarithmic connection has zero local curvature'] = (
        sp.simplify(sp.diff(sp.diff(f,yy)/f,xx)-sp.diff(sp.diff(f,xx)/f,yy)) == 0)
    Gr = sp.Matrix([[1,t],[t,1]])
    checks['two-channel metric reserve is one minus squared correlation'] = (
        sp.simplify(Gr.det()-(1-t*t)) == 0
        and equal(Gr*sp.Matrix([1,-1]), (1-t)*sp.Matrix([1,-1])))
    checks['actual YM baseline is inherited without enlargement'] = y19.y15.join_bounds(y19.RHO)['full_gap_lower'] > Q(9,40)
    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name,value in checks.items() if not value])
    finite = dict(energy=str(z), hidden_floor='1',
                  retained_metric=str(g[0,0]),
                  corrected_source=[str(v) for v in b],
                  full_residual=[str(v) for v in enclosure['residual']],
                  lower=str(enclosure['lower'][0,0]),
                  exact=str(enclosure['exact'][0,0]),
                  upper=str(enclosure['upper'][0,0]),
                  negative_count=0, actual_YM_computation=False)
    return y19.y12.encode(dict(
        stage='YC20', base_commit='c51ac1e', checks=checks,
        finite_control=finite, conditional_iteration_control=reserve,
        inherited_actual_YM_gap_lower=y19.y15.join_bounds(y19.RHO)['full_gap_lower'],
        claim_boundary=dict(
            actual_YC19_quotient_adapter_under_stated_domain_and_floor_hypotheses=True,
            cut_first_complementary_source_split=True,
            equilibrium_from_cut_alone=False, potential_postulated_as_primitive=False,
            complete_metric_and_source_reset=True, exact_energy_step=True,
            full_residual_Gram_certificate=True, nested_elimination_consistency=True,
            conditional_projective_metric_reserve=True,
            new_YM_block_metric_bound=False, new_YM_block_spectrum=False,
            local_metric_square_root=False, uniform_Hilbert_condition_number=False,
            improved_interaction_budget=False, improved_gap_window=False,
            iterated_spatial_RG=False, continuum_mass_gap=False,
            formal_machine_verification=False),
        idea_source_hashes=DRAFT_SOURCES, source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC20_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC20: {len(result['checks'])} exact checks pass.")
    print('Metric/source reset, full residual enclosure and nested energy-step controls verified.')
