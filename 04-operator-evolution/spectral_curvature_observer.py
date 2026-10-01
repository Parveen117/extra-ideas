"""R11: a spectral-blind curvature family and its model-relative observer repair.

All carriers and probe insertions are admitted, finite rational adapters.
Native composition and the tangent-quotient contract reuse R7-R10.
The spectral paper's insertion convention is trace(M D), not trace(M.T D).
"""

from itertools import combinations
from math import comb
import native_curvature_descent as n
from emk_curvature_balance import recognize
from emk_curvature_observation import trace

Q = n.Q


def lower_data(lower):
    lower = tuple(tuple(n.scalar(x) for x in row) for row in lower)
    if not lower:
        raise ValueError('Declare at least one hidden mode')
    return n.rectangular(lower, len(lower), 2)


def elementary(size, row, column):
    if not 0 <= row < size or not 0 <= column < size:
        raise ValueError('Matrix-unit index lies outside its carrier')
    return n.matrix([[int(i == row and j == column) for j in range(size)] for i in range(size)])


def hidden_projector(hidden):
    if not isinstance(hidden, int) or hidden < 1:
        raise ValueError('Declare a positive finite hidden dimension')
    return n.block(n.scale(n.identity(2), 0), n.rzero(2, hidden), n.rzero(hidden, 2), n.identity(hidden))


def nilpotent_carrier(lower):
    lower = lower_data(lower)
    h = len(lower)
    return n.block(n.scale(n.identity(2), 0), n.rzero(2, h), lower, n.scale(n.identity(h), 0))


def cut(hidden):
    return n.sub(n.identity(2+hidden), n.scale(hidden_projector(hidden), 2))


def native_connection(lower):
    lower = lower_data(lower)
    return n.ConnectionJet.constant(('u', 'v'), (hidden_projector(len(lower)), nilpotent_carrier(lower)))


def classical_report(lower):
    lower = lower_data(lower)
    return n.riemann_reduction(native_connection(lower), n.flat_metric(directions=('u', 'v')),
                              n.projection(2, 2+len(lower)))


def extract_lower(operator):
    operator = n.matrix(operator)
    h = len(operator)-2
    if h < 1:
        raise ValueError('The declared family needs two visible modes and a hidden sector')
    lower = tuple(row[:2] for row in operator[2:])
    if nilpotent_carrier(lower) != operator:
        raise ValueError('Operator lies outside the declared strict-lower curvature family')
    return lower


def order_factors(lower, a=Q(1, 3), b=Q(2, 5)):
    lower, a, b = lower_data(lower), n.scalar(a), n.scalar(b)
    if a == -1:
        raise ValueError('The finite return protocol requires invertible I+aP')
    unit = n.identity(2+len(lower))
    return n.add(unit, n.scale(hidden_projector(len(lower)), a)), n.add(unit, n.scale(nilpotent_carrier(lower), b))


def order_loop(lower, a=Q(1, 3), b=Q(2, 5)):
    """Four labelled moves have return X Y X^-1 Y^-1 = I+a*b*N."""
    x, y = order_factors(lower, a, b)
    maps = (n.inverse(y), n.inverse(x), y, x)
    vertices = ('p', 'a', 'b', 'c', 'p')
    moves = [n.Move(vertices[i], vertices[i+1], (('V', maps[i]),)) for i in range(4)]
    loop = moves[0]
    for move in moves[1:]:
        loop = loop.then(move)
    return loop


def determinant(operator):
    """Independent finite Gaussian determinant evaluator over Q."""
    work = [list(row) for row in n.matrix(operator)]
    value, size = Q(1), len(work)
    for column in range(size):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        if pivot is None:
            return Q(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        diagonal = work[column][column]
        value *= diagonal
        for i in range(column+1, size):
            ratio = work[i][column]/diagonal
            for j in range(column+1, size):
                work[i][j] -= ratio*work[column][j]
            work[i][column] = Q(0)
    return value


def spectral_coefficients(operator):
    """det(I-q*operator), via exact Newton trace identities."""
    operator = n.matrix(operator)
    size, power, moments = len(operator), n.identity(len(operator)), [Q(0)]
    for _ in range(size):
        power = n.mul(power, operator)
        moments.append(trace(power))
    coefficients = [Q(1)]
    for k in range(1, size+1):
        coefficients.append(-sum((coefficients[k-j]*moments[j] for j in range(1, k+1)), Q(0))/k)
    return tuple(coefficients)


def spectral_value(operator, parameter):
    operator, parameter = n.matrix(operator), n.scalar(parameter)
    return determinant(n.sub(n.identity(len(operator)), n.scale(operator, parameter)))


def marker_bank(hidden):
    hidden_projector(hidden)
    return tuple(elementary(2+hidden, column, 2+row) for row in range(hidden) for column in range(2))


def symmetric_marker_bank(hidden):
    return tuple(n.add(marker, n.transpose(marker)) for marker in marker_bank(hidden))


def nonfeedback_catalogue(hidden):
    hidden_projector(hidden)
    size = 2+hidden
    return tuple(elementary(size, row, column) for row in range(size) for column in range(size)
                 if not (row < 2 and column >= 2))


def marked_trace(marker, operator):
    marker, operator = n.same_carrier(marker, operator)
    return trace(n.mul(marker, operator))


def marked_determinant(marker, operator, parameter):
    marker, operator = n.same_carrier(marker, operator)
    return spectral_value(n.mul(marker, operator), parameter)


def feedback_report(lower, marker):
    """If the probe is admitted as a transport, its excursion is visible."""
    lower = lower_data(lower)
    marker, operator = n.same_carrier(marker, nilpotent_carrier(lower))
    connection = n.ConnectionJet.constant(('u', 'v'), (marker, operator))
    report = n.visible_decomposition(connection, n.flat_metric(directions=('u', 'v')))
    report['trace_visible_excursion'] = trace(report['excursion'])
    report['marked_curvature_response'] = marked_trace(marker, operator)
    return report


def analysis_rows(markers, hidden):
    basis = [nilpotent_carrier([[int(r == row and c == column) for c in range(2)] for r in range(hidden)])
             for row in range(hidden) for column in range(2)]
    return tuple(tuple(marked_trace(marker, operator) for operator in basis) for marker in markers)


def rank(rows, columns):
    return columns-len(n.nullspace(rows, columns))


def observer_report(rows, target_rows, columns):
    """Apply the EXISTING spectral paper's target-relative minimum-lift law."""
    rows = tuple(tuple(n.scalar(x) for x in row) for row in rows)
    target_rows = tuple(tuple(n.scalar(x) for x in row) for row in target_rows)
    if columns < 1 or any(len(row) != columns for row in rows+target_rows):
        raise ValueError('Observation and target rows must have the declared common domain')
    kernel = n.nullspace(rows, columns)
    restriction = tuple(tuple(sum((x*y for x, y in zip(row, vector)), Q(0)) for vector in kernel)
                        for row in target_rows)
    target_blindness = rank(restriction, len(kernel))
    return {'blind_kernel': kernel, 'blind_dimension': len(kernel),
            'target_blindness': target_blindness, 'minimum_arbitrary_scalar_supplements': target_blindness,
            'target_faithful': target_blindness == 0}


def catalogue_minimum(markers, hidden, target_rows):
    """Small exact catalogue search; counts scalar trace channels, not devices."""
    rows = analysis_rows(markers, hidden)
    columns = 2*hidden
    for count in range(len(rows)+1):
        for indices in combinations(range(len(rows)), count):
            if observer_report(tuple(rows[i] for i in indices), target_rows, columns)['target_faithful']:
                return {'minimum_catalogue_channels': count, 'indices': indices}
    return {'minimum_catalogue_channels': None, 'indices': None}


def balance_response_report(hidden, theta):
    theta = n.scalar(theta)
    if not 0 <= theta <= 1:
        raise ValueError('A declared branch weight must lie in [0,1]')
    factor, columns = 1-2*theta, 2*hidden
    rows = n.scale(n.identity(columns), factor)
    return {'odd_multiplier': factor, 'response_rows': rows,
            'response_Gram': n.mul(n.transpose(rows), rows),
            'observer': observer_report(rows, n.identity(columns), columns)}


def blindness_ledger(stages, target_rows, columns):
    """Successive row spaces must shrink; no lost target is charged twice."""
    reports = [observer_report(rows, target_rows, columns) for rows in stages]
    for prior, later in zip(stages, stages[1:]):
        if rank(tuple(prior)+tuple(later), columns) != rank(prior, columns):
            raise ValueError('A restriction ledger cannot silently introduce new observation rows')
    indices = tuple(report['target_blindness'] for report in reports)
    return {'target_blindness_by_stage': indices,
            'initial_blindness': indices[0] if indices else 0,
            'increments': tuple(b-a for a, b in zip(indices, indices[1:]))}


def decode_determinants(values, hidden, parameter, factor=1):
    parameter, factor = n.scalar(parameter), n.scalar(factor)
    values = tuple(n.scalar(x) for x in values)
    if len(values) != 2*hidden or parameter*factor == 0:
        raise ValueError('Full-bank decoding needs 2h scalar values and a nonzero calibrated gain')
    decoded = tuple((1-value)/(parameter*factor) for value in values)
    return n.rectangular((decoded[2*i:2*i+2] for i in range(hidden)), hidden, 2)


def classify_determinants(values, hidden, parameter, error=0, factor=1):
    error, parameter, factor = n.scalar(error), n.scalar(parameter), n.scalar(factor)
    if error < 0:
        raise ValueError('The supplied deterministic error bound must be nonnegative')
    decoded = decode_determinants(values, hidden, parameter, factor)
    radius = error/abs(parameter*factor)
    intervals = tuple((value-radius, value+radius) for row in decoded for value in row)
    if any(lo > 0 or hi < 0 for lo, hi in intervals):
        status = 'CERTIFIED_NONZERO_CURVATURE_IN_DECLARED_FAMILY'
    elif error == 0:
        status = 'CERTIFIED_ZERO_CURVATURE_IN_DECLARED_FAMILY'
    else:
        status = 'UNRESOLVED_CURVATURE_AT_DECLARED_ERROR'
    return {'status': status, 'lower_estimate': decoded, 'entry_intervals': intervals,
            'squared_Frobenius_error_bound': 2*hidden*radius*radius}


def cubic_connected_difference(a, b, c, parameter):
    """Independent coefficient [xyz](L_CBA-L_CAB), not a Magnus approximation.

    Polynomial matrices use x^2=y^2=z^2=0, which preserves this coefficient.
    Only W, W^2 and W^3 can contribute to trace log(I-rW).
    """
    a, b, c = n.same_carrier(a, b, c)
    parameter = n.scalar(parameter)
    if parameter == 1:
        raise ValueError('The normalized connected coefficient is singular at q=1')
    unit, zero = n.identity(len(a)), n.scale(a, 0)

    def product(left, right):
        out = {}
        for mask, value in left.items():
            for other_mask, other in right.items():
                if not mask & other_mask:
                    key = mask | other_mask
                    out[key] = n.add(out.get(key, zero), n.mul(value, other))
        return out

    pa, pb, pc = {0: unit, 1: a}, {0: unit, 2: b}, {0: unit, 4: c}
    endpoints = (product(pc, product(pb, pa)), product(pc, product(pa, pb)))
    responses = []
    ratio = parameter/(1-parameter)
    for endpoint in endpoints:
        w = {mask: value for mask, value in endpoint.items() if mask}
        power, response = w, Q(0)
        for k in range(1, 4):
            response -= ratio**k*trace(power.get(7, zero))/k
            power = product(power, w)
        responses.append(response)
    return responses[0]-responses[1]


def exact_witnesses():
    lower = lower_data(((2, Q(-3, 4)), (Q(5, 7), 1)))
    hidden, operator = len(lower), nilpotent_carrier(lower)
    projector, markers = hidden_projector(hidden), marker_bank(hidden)
    a, b, parameter = Q(1, 3), Q(2, 5), Q(2)
    loop = dict(order_loop(lower, a, b).operators)['V']
    values = tuple(marked_determinant(marker, operator, parameter) for marker in markers)
    zero_rows = (tuple(Q(0) for _ in range(2*hidden)),)
    return {'lower_coupling': lower, 'full_curvature': native_connection(lower).curvature(0, 1),
            'curvature_spectral_coefficients': spectral_coefficients(operator),
            'Riemann_reduction': classical_report(lower),
            'finite_order_loop': {'a': a, 'b': b, 'carrier': loop,
                                  'spectral_coefficients': spectral_coefficients(loop),
                                  'identity_spectral_coefficients': tuple(Q((-1)**k*comb(2+hidden, k)) for k in range(3+hidden))},
            'marked_determinants': {'parameter': parameter, 'values': values,
                                    'decoded_lower': decode_determinants(values, hidden, parameter)},
            'minimal_scalar_repair': observer_report(zero_rows, n.identity(2*hidden), 2*hidden),
            'no_feedback_catalogue': observer_report(analysis_rows(nonfeedback_catalogue(hidden), hidden),
                                                      n.identity(2*hidden), 2*hidden),
            'cubic_source_channel': {'computed': cubic_connected_difference(projector, operator, markers[0], parameter),
                                     'source_formula': parameter/(1-parameter)**2*marked_trace(markers[0], operator)},
            'probe_excursion_connection': feedback_report(lower, markers[0]),
            'balanced_response': balance_response_report(hidden, Q(1, 2)),
            'loss_ledger': blindness_ledger((n.identity(2*hidden), zero_rows, zero_rows, zero_rows),
                                            n.identity(2*hidden), 2*hidden)}
