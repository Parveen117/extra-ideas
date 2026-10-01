"""R12: derive a future observer and operational metric from a native event law.

No metric is primitive. A finite real native representation, calibrated Gaussian
record noise, admitted transports and an event stopping law specify an experiment.
Its tagged Fisher form solves an exact positive Stein equation. The finite
linear algebra here is an experiment adapter; the canonical native engine and
its N03 observer completion remain in RKF/operator_foundation.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
import native_curvature_descent as n
import spectral_curvature_observer as spectral


def rect(rows):
    value = tuple(tuple(n.scalar(x) for x in row) for row in rows)
    if not value or not value[0] or any(len(row) != len(value[0]) for row in value):
        raise ValueError('A nonempty rectangular response map is required')
    return value


def product(*factors):
    value = factors[0]
    for later in factors[1:]:
        value = n.rproduct(value, later)
    return value


def positive_semidefinite(value):
    """Exact symmetric Schur test, including zero pivots and singular forms."""
    a = n.matrix(value)
    if a != n.transpose(a):
        return False
    while a:
        if any(a[i][i] < 0 for i in range(len(a))):
            return False
        pivot = next((i for i in range(len(a)) if a[i][i] > 0), None)
        if pivot is None:
            return n.is_zero(a)
        rest = [i for i in range(len(a)) if i != pivot]
        a = tuple(tuple(a[i][j]-a[i][pivot]*a[pivot][j]/a[pivot][pivot]
                        for j in rest) for i in rest)
    return True


def positive_definite(value):
    a = n.matrix(value)
    return (a == n.transpose(a) and all(spectral.determinant(
        tuple(row[:k] for row in a[:k])) > 0 for k in range(1, len(a)+1)))


def independent_rows(rows):
    """Choose quotient coordinates from a solved form, without word saturation."""
    work = [list(row) for row in rect(rows)]
    position = 0
    for col in range(len(work[0])):
        pivot = next((i for i in range(position, len(work)) if work[i][col]), None)
        if pivot is None:
            continue
        work[position], work[pivot] = work[pivot], work[position]
        divisor = work[position][col]
        work[position] = [x/divisor for x in work[position]]
        for i in range(len(work)):
            if i != position:
                factor = work[i][col]
                work[i] = [x-factor*y for x, y in zip(work[i], work[position])]
        position += 1
        if position == len(work):
            break
    return tuple(tuple(row) for row in work[:position])


@dataclass(frozen=True)
class EventExperiment:
    readout: tuple
    precision: tuple
    transports: tuple
    probabilities: tuple
    continuation: Q
    stability_form: tuple
    contraction: Q

    def __post_init__(self):
        d, w = rect(self.readout), n.matrix(self.precision)
        size = len(d[0])
        ts = tuple(n.matrix(t) for t in self.transports)
        ps = tuple(n.scalar(p) for p in self.probabilities)
        r, beta = n.scalar(self.continuation), n.scalar(self.contraction)
        h = n.matrix(self.stability_form)
        if len(w) != len(d) or not positive_definite(w):
            raise ValueError('Record precision must be positive definite on every record channel')
        if not ts or len(ts) != len(ps) or sum(ps) != 1 or any(p <= 0 for p in ps):
            raise ValueError('Every admitted event has positive state-independent probability')
        if not 0 < r < 1 or not 0 <= beta < 1:
            raise ValueError('Stopping and contraction must give a convergent event law')
        if len(h) != size or not positive_definite(h) or any(len(t) != size for t in ts):
            raise ValueError('The stability witness must cover the full native carrier')
        for t in ts:
            n.inverse(t)  # R7 lawful transport; a resetting map is a different contract.
        averaged = n.scale(h, 0)
        for p, t in zip(ps, ts):
            averaged = n.add(averaged, n.scale(product(n.transpose(t), h, t), r*p))
        margin = n.sub(n.scale(h, beta), averaged)
        if not positive_semidefinite(margin):
            raise ValueError('The supplied stability inequality is false')
        for name, value in (('readout', d), ('precision', w), ('transports', ts),
                            ('probabilities', ps), ('continuation', r),
                            ('stability_form', h), ('contraction', beta)):
            object.__setattr__(self, name, value)

    @property
    def seed(self):
        return product(n.transpose(self.readout), self.precision, self.readout)


def event_operator(experiment, form):
    result = n.scale(form, 0)
    for p, t in zip(experiment.probabilities, experiment.transports):
        result = n.add(result, n.scale(product(n.transpose(t), form, t),
                                       experiment.continuation*p))
    return result


def stein_solve(experiment, rhs):
    """Solve (I-event_operator)G=rhs on symmetric forms over exact Q."""
    rhs = n.matrix(rhs)
    size = len(experiment.stability_form)
    if len(rhs) != size or rhs != n.transpose(rhs):
        raise ValueError('A symmetric same-carrier source is required')
    pairs = [(i, j) for i in range(size) for j in range(i, size)]
    basis = tuple(n.matrix([[int((a, b) == (i, j) or (a, b) == (j, i))
                              for b in range(size)] for a in range(size)]) for i, j in pairs)
    columns = tuple(n.sub(b, event_operator(experiment, b)) for b in basis)
    system = tuple(tuple(col[i][j] for col in columns) for i, j in pairs)
    vector = tuple((rhs[i][j],) for i, j in pairs)
    coefficients = product(n.inverse(system), vector)
    result = n.scale(rhs, 0)
    for coefficient, b in zip(coefficients, basis):
        result = n.add(result, n.scale(b, coefficient[0]))
    if n.sub(result, event_operator(experiment, result)) != rhs:
        raise ArithmeticError('The exact Stein solve did not reproduce its source')
    return result


def tagged_metric(experiment):
    return stein_solve(experiment, n.scale(experiment.seed, 1-experiment.continuation))


def observer_from_metric(form, experiment):
    c = independent_rows(form)
    size = len(form)
    kernel = n.nullspace(c, size)
    if not c:
        return {'observer': (), 'rank': 0, 'kernel': kernel, 'section': (),
                'quotient_metric': (), 'descended_transports': ()}
    section = product(n.transpose(c), n.inverse(product(c, n.transpose(c))))
    g = product(n.transpose(section), form, section)
    actions = tuple(product(c, t, section) for t in experiment.transports)
    if (not positive_definite(g) or product(n.transpose(c), g, c) != form
            or product(product(experiment.readout, section), c) != experiment.readout
            or any(product(a, c) != product(c, t)
                   for a, t in zip(actions, experiment.transports))):
        raise ArithmeticError('The solved response form did not descend to its observer quotient')
    return {'observer': c, 'rank': len(c), 'kernel': kernel, 'section': section,
            'quotient_metric': g, 'descended_transports': actions}


def derive_observer_metric(experiment):
    form = tagged_metric(experiment)
    if not positive_semidefinite(form):
        raise ArithmeticError('A stable tagged Fisher form cannot be indefinite')
    result = observer_from_metric(form, experiment)
    result.update({'metric': form, 'seed': experiment.seed,
                   'stein_residual': n.sub(n.sub(form, event_operator(experiment, form)),
                                          n.scale(experiment.seed, 1-experiment.continuation))})
    return result


def origin_untagged_metric(experiment):
    """Fisher form after forgetting every event tag, at the zero state only.

All conditional Gaussian means coincide there. The mixture score is the
average conditional score. At nonzero states this mean formula is not valid.
"""
    size = len(experiment.stability_form)
    mean = n.scale(n.identity(size), 0)
    for p, t in zip(experiment.probabilities, experiment.transports):
        mean = n.add(mean, n.scale(t, p))
    history_mean = n.scale(n.inverse(n.sub(n.identity(size),
                                          n.scale(mean, experiment.continuation))),
                          1-experiment.continuation)
    response = product(experiment.readout, history_mean)
    form = product(n.transpose(response), experiment.precision, response)
    discard = n.sub(tagged_metric(experiment), form)
    if not positive_semidefinite(discard):
        raise ArithmeticError('Forgetting a tag cannot create local Fisher information')
    return {'mean_transport': history_mean, 'mean_response': response,
            'metric': form, 'discard': discard, 'state_scope': 'zero state only'}


def finite_tagged_metric(experiment, maximum_length):
    if not isinstance(maximum_length, int) or maximum_length < 0:
        raise ValueError('Supply a nonnegative event aperture')
    term = n.scale(experiment.seed, 1-experiment.continuation)
    result = term
    for _ in range(maximum_length):
        term = event_operator(experiment, term)
        result = n.add(result, term)
    return result


def response_two_jet(experiment, first_transports, second_transports,
                     first_seed=None, second_seed=None):
    """Differentiate the same Stein selection equation twice.

Probabilities, continuation and record calibration are constant. Transport
derivatives include all product terms; no metric derivative is supplied.
"""
    count, catalogue = len(first_transports), len(experiment.transports)
    if (not count or len(second_transports) != count
            or any(len(row) != count for row in second_transports)
            or any(len(row) != catalogue for row in first_transports)
            or any(len(row) != catalogue for plane in second_transports for row in plane)):
        raise ValueError('Declare a complete transport two-jet')
    first = tuple(tuple(n.matrix(t) for t in row) for row in first_transports)
    second = tuple(tuple(tuple(n.matrix(t) for t in row) for row in plane)
                   for plane in second_transports)
    n.same_carrier(experiment.stability_form, *(t for row in first for t in row),
                   *(t for plane in second for row in plane for t in row))
    if any(second[i][j] != second[j][i] for i in range(count) for j in range(count)):
        raise ValueError('Coordinate transport second partials must commute')
    zero = n.scale(experiment.seed, 0)
    qi = tuple(n.matrix(t) for t in first_seed) if first_seed is not None else (zero,)*count
    qij = (tuple(tuple(n.matrix(t) for t in row) for row in second_seed)
           if second_seed is not None else tuple((zero,)*count for _ in range(count)))
    if (len(qi) != count or len(qij) != count or any(len(row) != count for row in qij)
            or any(qij[i][j] != qij[j][i] for i in range(count) for j in range(count))):
        raise ValueError('Declare symmetric seed derivatives of the same two-jet')
    g = tagged_metric(experiment)
    gi = []
    for i in range(count):
        rhs = n.scale(qi[i], 1-experiment.continuation)
        for a, (p, t) in enumerate(zip(experiment.probabilities, experiment.transports)):
            ti = first[i][a]
            rhs = n.add(rhs, n.scale(n.add(product(n.transpose(ti), g, t),
                                           product(n.transpose(t), g, ti)),
                                    experiment.continuation*p))
        gi.append(stein_solve(experiment, rhs))
    gij = []
    for i in range(count):
        row = []
        for j in range(count):
            rhs = n.scale(qij[i][j], 1-experiment.continuation)
            for a, (p, t) in enumerate(zip(experiment.probabilities, experiment.transports)):
                ti, tj, tij = first[i][a], first[j][a], second[i][j][a]
                terms = (product(n.transpose(tij), g, t), product(n.transpose(t), g, tij),
                         product(n.transpose(ti), gi[j], t), product(n.transpose(t), gi[j], ti),
                         product(n.transpose(tj), gi[i], t), product(n.transpose(t), gi[i], tj),
                         product(n.transpose(ti), g, tj), product(n.transpose(tj), g, ti))
                for term in terms:
                    rhs = n.add(rhs, n.scale(term, experiment.continuation*p))
            row.append(stein_solve(experiment, rhs))
        gij.append(tuple(row))
    return {'value': g, 'first': tuple(gi), 'second': tuple(gij)}


def tangent_metric(jet, tangent_map, directions):
    """An admitted constant tangent soldering pulls back the derived response form."""
    e = rect(tangent_map)
    if len(e) != len(jet['value']) or len(e[0]) != len(directions):
        raise ValueError('Tangent modes and response carrier do not match')
    pull = lambda g: product(n.transpose(e), g, e)
    value = pull(jet['value'])
    if not positive_definite(value):
        raise ValueError('The admitted tangent modes contain an invisible direction')
    return n.MetricJet(tuple(directions), value, tuple(pull(g) for g in jet['first']),
                       tuple(tuple(pull(g) for g in row) for row in jet['second']))


def seam_experiment(coupling, normal, continuation):
    c, v, r = n.scalar(coupling), n.scalar(normal), n.scalar(continuation)
    if not 0 < r < 1:
        raise ValueError('The native stopping law requires 0<continuation<1')
    nilpotent = spectral.nilpotent_carrier(((1, 0), (0, 0)))
    ts = (n.add(n.identity(4), n.scale(nilpotent, c*v)),
          n.sub(n.identity(4), n.scale(nilpotent, c*v)))
    kappa = c*c*r/(1-r)
    h = n.identity(4)
    h = n.add(h, n.matrix([[2*kappa*v*v, 0, 0, 0], [0, 0, 0, 0],
                           [0, 0, 0, 0], [0, 0, 0, 0]]))
    return EventExperiment(n.projection(3, 4), n.identity(3), ts, (Q(1, 2),)*2,
                           r, h, (1+r)/2)


def seam_foundation(coupling=Q(2), normal=Q(1, 3), continuation=Q(1, 2)):
    c, v, r = n.scalar(coupling), n.scalar(normal), n.scalar(continuation)
    experiment = seam_experiment(c, v, r)
    nilpotent = spectral.nilpotent_carrier(((1, 0), (0, 0)))
    zero = n.scale(nilpotent, 0)
    first = ((zero, zero), (n.scale(nilpotent, c), n.scale(nilpotent, -c)))
    second = tuple(tuple((zero, zero) for _ in range(2)) for _ in range(2))
    jet = response_two_jet(experiment, first, second)
    tangent = n.transpose(n.projection(2, 4))
    g = tangent_metric(jet, tangent, ('u', 'v'))
    native = spectral.native_connection(((c, 0), (0, 0)))
    return {'experiment': experiment, 'observer_metric': derive_observer_metric(experiment),
            'untagged_at_origin': origin_untagged_metric(experiment), 'metric_jet': g,
            'kappa': c*c*r/(1-r), 'levi_civita': g.levi_civita(),
            'raw_native_curvature': native.curvature(0, 1),
            'native_visible_decomposition': n.visible_decomposition(native, g),
            'native_connection_is_levi_civita': n.riemann_reduction(native, g, n.projection(2, 4))[
                'Riemann_sector_certified']}


def feedback_experiment():
    """R11 markers admitted as actions, with both hidden modes accessible."""
    lower = (((1, 0), (0, 0)), ((0, 0), (0, 1)))
    generators = tuple(spectral.nilpotent_carrier(data) for data in lower) + (
        spectral.marker_bank(2)[0], spectral.marker_bank(2)[3])
    transports = tuple(n.add(n.identity(4), n.scale(a, sign))
                       for a in generators for sign in (1, -1))
    return EventExperiment(n.projection(2, 4), n.identity(2), transports,
                           (Q(1, 8),)*8, Q(1, 3), n.identity(4), Q(1, 2))


def exact_witnesses():
    s = seam_foundation()
    g = s['metric_jet']
    d = s['native_visible_decomposition']
    return {'balanced_sheet': {
        'coupling': Q(2), 'normal_coordinate': Q(1, 3), 'continuation': Q(1, 2),
        'derived_kappa': s['kappa'], 'derived_observer': s['observer_metric']['observer'],
        'derived_response_metric': s['observer_metric']['metric'], 'tangent_metric': g.value,
        'tangent_metric_first': g.first, 'tangent_metric_second': g.second,
        'gaussian_curvature': -s['kappa']/(g.value[0][0]**2),
        'forgotten_tag_metric_at_origin': s['untagged_at_origin']['metric'],
        'forgotten_information': s['untagged_at_origin']['discard'],
        'raw_native_curvature': s['raw_native_curvature'],
        'riemann_curvature': d['Riemann'], 'distortion_curvature': d['distortion'],
        'excursion_curvature': d['excursion'],
        'native_connection_is_levi_civita': s['native_connection_is_levi_civita']},
        'feedback_catalogue': derive_observer_metric(feedback_experiment())}
