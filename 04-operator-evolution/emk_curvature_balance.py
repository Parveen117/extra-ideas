"""R9: a scoped native cut-balance protocol and its information ledger.

The K cut-swap is inherited from EMK-T1. Equal recognition of its two
branches selects the balanced endpoint; a physical rate is not selected.
Generic algebra laws need K squared = I. Quantum channel statements use
the admitted Hermitian/unitary two-mode chart and exact QTH-1 source.
"""

from fractions import Fraction as Q
from aghora_return import transpose
from emk_curvature_observation import (
    K, R, RK, add, commutator, identity, is_zero, matrix, mul, same_carrier,
    scalar, scale, sub, trace,
)


def parameter(theta):
    theta = scalar(theta)
    if not 0 <= theta <= 1:
        raise ValueError('A recognition branch weight must lie in [0,1]')
    return theta


def cut_action(a, cut=K):
    a, cut = same_carrier(a, cut)
    if mul(cut, cut) != identity(len(cut)):
        raise ValueError('The declared cut must square to identity')
    return mul(mul(cut, a), cut)


def grade(a, cut=K):
    a = matrix(a)
    changed = cut_action(a, cut)
    return scale(add(a, changed), Q(1, 2)), scale(sub(a, changed), Q(1, 2))


def recognize(a, theta=Q(1, 2), cut=K):
    """Phi_theta(A)=(1-theta) A + theta K A K."""
    theta = parameter(theta)
    a = matrix(a)
    return add(scale(a, 1-theta), scale(cut_action(a, cut), theta))


def iterate(a, branch_weights, cut=K):
    """Clock-free ordered recognition events, with an exact closed form."""
    a = matrix(a)
    even, odd = grade(a, cut)
    factor, value = Q(1), a
    weights = tuple(parameter(x) for x in branch_weights)
    for theta in weights:
        value = recognize(value, theta, cut)
        factor *= 1-2*theta
    predicted = add(even, scale(odd, factor))
    return {'value': value, 'even': even, 'initial_odd': odd,
            'odd_multiplier': factor, 'weights': weights,
            'closed_form_residual': sub(value, predicted),
            'cut_balanced': is_zero(grade(value, cut)[1])}


def curvature_selection(a, b, cut=K):
    a, b = same_carrier(a, b)
    ae, ao = grade(a, cut)
    be, bo = grade(b, cut)
    full = commutator(a, b)
    observed = recognize(full, cut=cut)
    recomputed = commutator(ae, be)
    odd_pair_memory = commutator(ao, bo)
    return {'raw': full, 'recognized_raw': observed,
            'recomputed_from_recognized_generators': recomputed,
            'odd_pair_memory': odd_pair_memory,
            'identity_residual': sub(observed, add(recomputed, odd_pair_memory))}


def validate_real_state(rho):
    rho = matrix(rho)
    if len(rho) != 2 or transpose(rho) != rho or trace(rho) != 1:
        raise ValueError('This state adapter requires a real symmetric two-mode density matrix')
    if rho[0][0] < 0 or rho[1][1] < 0 or rho[0][0]*rho[1][1]-rho[0][1]**2 < 0:
        raise ValueError('The declared state must be positive semidefinite')
    return rho


def axis_state(m):
    m = scalar(m)
    if abs(m) > 1:
        raise ValueError('Axis imbalance must lie in [-1,1]')
    return matrix([[(1+m)/2, 0], [0, (1-m)/2]])


def full_rank_axis(m):
    m = scalar(m)
    if abs(m) >= 1:
        raise ValueError('The Fisher/SLD adapter requires |m|<1')
    return m


def curvature_budget(rho, hermitian_probe):
    """Signed mean plus variance equals the second moment, with fixed probes."""
    rho = validate_real_state(rho)
    rho, h = same_carrier(rho, hermitian_probe)
    if transpose(h) != h:
        raise ValueError('The variance budget requires a real Hermitian probe')
    mean = trace(mul(rho, h))
    second = trace(mul(rho, mul(h, h)))
    variance = second-mean*mean
    if variance < 0:
        raise ValueError('A positive-state Hermitian variance cannot be negative')
    return {'mean': mean, 'mean_squared': mean*mean, 'second_moment': second,
            'variance': variance, 'budget_residual': second-mean*mean-variance}


def information_pushforward(m, theta):
    """The SAME family rho_m(x,y) is pushed through the declared cut channel.

    Values AND derivatives transform: (x,y,m)->(x,lambda*y,lambda*m).
    Recomputed QTH-1 tensor is g'=diag(1,lambda^2), B'_xy=lambda^2*m.
    A retained branch flag permits branchwise inverse K^b.
    """
    m, theta = full_rank_axis(m), parameter(theta)
    lam = 1-2*theta
    before = identity(2)
    after = matrix([[1, 0], [0, lam*lam]])
    return {'initial_imbalance': m, 'theta': theta, 'lambda': lam,
            'information_before': before, 'information_after': after,
            'discarded_information': sub(before, after),
            'flagged_branch_information': before,
            'fixed_probe_curvature_mean_after': lam*m,
            'recomputed_tensor_curvature_after': lam**2*m,
            'post_channel_SLD_pair_noncommutes': lam != 0}


def source_information_probe(qth, m, theta):
    """Replay analytic pushforward values through unchanged QTH-1 operations."""
    report = information_pushforward(m, theta)
    m, lam = report['initial_imbalance'], report['lambda']
    rho = qth.bloch_state((Q(0), Q(0), lam*m))
    dx = qth.bloch_derivative((Q(1), Q(0), Q(0)))
    dy = qth.bloch_derivative((Q(0), lam, Q(0)))
    gxx, _ = qth.geometric_tensor(rho, dx, dx)
    gyy, _ = qth.geometric_tensor(rho, dy, dy)
    gxy, curvature = qth.geometric_tensor(rho, dx, dy)
    report['source_information_after'] = matrix([[gxx, gxy], [gxy, gyy]])
    report['source_curvature_after'] = curvature
    report['source_SLD_commutator_nonzero'] = (
        qth.msub(qth.mmul(qth.sld(rho, dx), qth.sld(rho, dy)),
                  qth.mmul(qth.sld(rho, dy), qth.sld(rho, dx))) != qth.mzero())
    branches = []
    for branch, probability in ((0, 1-report['theta']), (1, report['theta'])):
        sign = 1-2*branch
        state = qth.bloch_state((Q(0), Q(0), sign*m))
        bx, by = dx, qth.bloch_derivative((Q(0), Q(sign), Q(0)))
        a, _ = qth.geometric_tensor(state, bx, bx)
        b, _ = qth.geometric_tensor(state, by, by)
        ab, _ = qth.geometric_tensor(state, bx, by)
        branches.append((probability, matrix([[a, ab], [ab, b]])))
    total = scale(identity(2), 0)
    for probability, information in branches:
        total = add(total, scale(information, probability))
    report['source_flagged_branch_information'] = total
    return report


def measurement_information(m, effects):
    """Finite one-copy qubit effects E_a=w_a(I+n_a.sigma).

    Inputs explicitly certify positivity and completeness. The two tangent
    directions are x,y at (0,0,m), with QTH information I_2.
    """
    m = full_rank_axis(m)
    effects = tuple((scalar(w), tuple(map(scalar, n))) for w, n in effects)
    if (not effects or any(w < 0 or len(n) != 3 or sum(x*x for x in n) > 1 for w, n in effects)
            or sum(w for w, n in effects) != 1
            or any(sum(w*n[i] for w, n in effects) != 0 for i in range(3))):
        raise ValueError('Positive complete rational qubit effects are required')
    c = scale(identity(2), 0)
    probabilities = []
    for w, n in effects:
        nx, ny, nz = n
        denominator = 1+m*nz
        probabilities.append(w*denominator)
        c = add(c, scale(matrix([[nx*nx, nx*ny], [nx*ny, ny*ny]]), w/denominator))
    return {'classical_information': c, 'quantum_information': identity(2),
            'discarded_information': sub(identity(2), c),
            'classical_trace': trace(c), 'joint_trace_gap': 2-trace(c),
            'probabilities': tuple(probabilities)}


def projective_effects(direction):
    direction = tuple(map(scalar, direction))
    if len(direction) != 3 or sum(x*x for x in direction) != 1:
        raise ValueError('A rational unit measurement direction is required')
    return ((Q(1, 2), direction), (Q(1, 2), tuple(-x for x in direction)))


def exact_witnesses(qth, master):
    z = scale(K, 0)
    radial = curvature_selection(K, R)
    tangential = curvature_selection(R, RK)
    # A cut-balanced state can retain an even tangential curvature mean.
    rho_even = scale(add(identity(2), scale(K, Q(1, 2))), Q(1, 2))
    h = scale(commutator(K, R), Q(1, 2))
    initial = axis_state(Q(3, 5))
    balanced = recognize(initial)
    unflagged = source_information_probe(qth, Q(3, 5), Q(1, 2))
    master_channels = [(K, z), (z, R)]
    full_master = matrix(master.master_curvature(master_channels))
    return {
        'radial_tangential': radial, 'tangential_tangential': tangential,
        'balanced_even_state': rho_even,
        'surviving_even_curvature_mean': trace(mul(rho_even, tangential['raw'])),
        'nonzero_raw_master_with_zero_individual_channels': {
            'raw_master': full_master,
            'individual_channels': tuple(matrix(x) for x in master.channel_curvatures(master_channels)),
            'mixed': matrix(master.mixed_curvature(master_channels)),
            'recognized_raw': recognize(full_master)},
        'finite_balance_sequence': iterate(initial, (Q(1, 4), Q(3, 4), Q(1, 2))),
        'initial_curvature_budget': curvature_budget(initial, h),
        'balanced_curvature_budget': curvature_budget(balanced, h),
        'balanced_information_pushforward': unflagged,
        'balanced_projective_measurement': measurement_information(Q(0), projective_effects((1, 0, 0))),
        'balanced_equal_axis_measurement': measurement_information(Q(0), (
            (Q(1, 4), (1, 0, 0)), (Q(1, 4), (-1, 0, 0)),
            (Q(1, 4), (0, 1, 0)), (Q(1, 4), (0, -1, 0)))),
    }
