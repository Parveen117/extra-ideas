"""R13: complementary native transport generates memory, noise and event data.

Exact hidden-state elimination needs no stochastic primitive. A specified
balanced unit-norm preparation supplies endogenous force covariance. Additive
native norm weights for a monitored cut supply event probabilities and a
non-Gaussian event Fisher metric. Coherent retention and monitored reset are
different protocols; the canonical operator engine remains unchanged in RKF.
"""

from fractions import Fraction as Q
import observer_metric_foundation as o
from aghora_return import cayley_flow

n = o.n


def steps(value):
    if type(value) is not int or value < 0:
        raise ValueError('Use a nonnegative integer event count')
    return value


def power(value, exponent):
    result = n.identity(len(value))
    for _ in range(steps(exponent)):
        result = n.mul(value, result)
    return result


def rotation(cosine, sine):
    c, s = n.scalar(cosine), n.scalar(sine)
    if c*c+s*s != 1:
        raise ValueError('The native rotation chart requires cosine^2+sine^2=1')
    return n.add(n.scale(n.identity(2), c), n.scale(n.R, s))


def cayley_rotation(parameter):
    a = n.scalar(parameter)
    c, s = (1-a*a)/(1+a*a), 2*a/(1+a*a)
    value = rotation(c, s)
    if value != cayley_flow(n.R, a):
        raise ArithmeticError('The exact native Cayley chart did not match')
    return {'cosine': c, 'sine': s, 'transport': value}


def split(transport, visible):
    t = n.matrix(transport)
    n.inverse(t)
    if type(visible) is not int or not 0 < visible < len(t):
        raise ValueError('Declare nonempty visible and complementary modes')
    return {'transport': t, 'A': tuple(row[:visible] for row in t[:visible]),
            'B': tuple(row[visible:] for row in t[:visible]),
            'C': tuple(row[:visible] for row in t[visible:]),
            'D': tuple(row[visible:] for row in t[visible:])}


def memory_kernels(transport, visible, count):
    m = split(transport, visible)
    return tuple(o.product(m['B'], power(m['D'], lag), m['C'])
                 for lag in range(steps(count)))


def reduced_history(transport, visible, initial_visible, initial_hidden, count):
    """Compare full evolution with exact A + memory + hidden-initial-force steps."""
    m = split(transport, visible)
    x, z = o.rect(initial_visible), o.rect(initial_hidden)
    if (len(x) != visible or len(z) != len(m['D']) or len(x[0]) != 1 or len(z[0]) != 1):
        raise ValueError('Initial visible and hidden columns must match the split')
    state = x+z
    full, xs, forces, memory, markov, residuals = [state], [x], [], [], [], []
    kernels = memory_kernels(transport, visible, count)
    for k in range(steps(count)):
        a = o.product(m['A'], xs[k])
        history = n.rzero(visible, 1)
        for j in range(k):
            history = tuple(tuple(u+v for u, v in zip(ar, br))
                            for ar, br in zip(history, o.product(kernels[k-1-j], xs[j])))
        force = o.product(m['B'], power(m['D'], k), z)
        predicted = tuple(tuple(u+v+w for u, v, w in zip(ar, br, cr))
                          for ar, br, cr in zip(a, history, force))
        state = o.product(m['transport'], state)
        actual = state[:visible]
        full.append(state)
        xs.append(actual)
        forces.append(force)
        memory.append(history)
        markov.append(a)
        residuals.append(n.rsubtract(actual, predicted))
    return {'full_states': tuple(full), 'visible_states': tuple(xs),
            'markov_terms': tuple(markov), 'memory_terms': tuple(memory),
            'complement_forces': tuple(forces), 'residuals': tuple(residuals), 'kernels': kernels}


def complementary_covariance(transport, visible, initial_covariance, count):
    """Cov(xi_k,xi_l)=B D^k S (D^l)^T B^T, with no Gaussian assumption."""
    m = split(transport, visible)
    sigma = n.matrix(initial_covariance)
    if len(sigma) != len(m['D']) or not o.positive_semidefinite(sigma):
        raise ValueError('Supply a positive hidden initial covariance on the complementary carrier')
    maps = tuple(o.product(m['B'], power(m['D'], k)) for k in range(steps(count)))
    return tuple(tuple(o.product(a, sigma, n.transpose(b)) for b in maps) for a in maps)


def balanced_unit_pair(visible_amplitude, hidden_amplitude):
    """Native fixed-norm preparation and exchange balance determine all moments."""
    x, z = n.scalar(visible_amplitude), n.scalar(hidden_amplitude)
    if x*x+z*z != 1:
        raise ValueError('The selected native initial preparation has unit total norm')
    return {'visible': x, 'hidden_branches': (z, -z), 'weights': (Q(1, 2),)*2,
            'hidden_mean': Q(0), 'hidden_variance': 1-x*x,
            'hidden_fourth_moment': (1-x*x)**2,
            'gaussian_hidden_law': False}


def endogenous_noise(cosine, sine, visible_amplitude, hidden_amplitude, count):
    t = rotation(cosine, sine)
    c, s = t[0][0], t[1][0]
    preparation = balanced_unit_pair(visible_amplitude, hidden_amplitude)
    variance = preparation['hidden_variance']
    covariance = complementary_covariance(t, 1, ((variance,),), count)
    kernels = memory_kernels(t, 1, count)
    noises = tuple(tuple(-s*(c**k)*z for k in range(steps(count)))
                   for z in preparation['hidden_branches'])
    mean = tuple(sum(p*row[k] for p, row in zip(preparation['weights'], noises))
                 for k in range(count))
    covariance_by_branches = tuple(tuple(sum(p*row[k]*row[l]
                                               for p, row in zip(preparation['weights'], noises))
                                          for l in range(count)) for k in range(count))
    return {'preparation': preparation, 'noise_paths': noises, 'noise_mean': mean,
            'covariance': covariance, 'covariance_by_branches': covariance_by_branches,
            'memory_kernels': kernels,
            'memory_noise_residual': tuple(covariance[k][0][0][0]+variance*kernels[k][0][0]
                                           for k in range(count)),
            'temporal_rank': len(o.independent_rows(covariance_by_branches)) if count else 0}


def recover_complement(cosine, sine, first_visible, next_visible):
    t = rotation(cosine, sine)
    c, s = t[0][0], t[1][0]
    if s == 0:
        raise ValueError('The complementary mode never reaches this visible channel')
    x0, x1 = n.scalar(first_visible), n.scalar(next_visible)
    return (c*x0-x1)/s


def complement_error_radius(cosine, sine, first_error, next_error):
    t = rotation(cosine, sine)
    e0, e1 = n.scalar(first_error), n.scalar(next_error)
    if e0 < 0 or e1 < 0 or t[1][0] == 0:
        raise ValueError('Use nonnegative response errors and a visible complementary gain')
    return (abs(t[0][0])*e0+e1)/abs(t[1][0])


def branch_weights(cosine, sine):
    t = rotation(cosine, sine)
    c, s = t[0][0], t[1][0]
    return {'continue': c*c, 'exit': s*s,
            'binary_event_variance': c*c*s*s, 'native_split_norm': c*c+s*s}


def cut_incompatibility(cosine, sine):
    t = rotation(cosine, sine)
    j = n.scale(n.RK, -1)
    p = n.scale(n.add(n.identity(2), j), Q(1, 2))
    q = n.sub(n.identity(2), p)
    mismatch = n.commutator(p, t)
    retained = n.mul(n.mul(p, t), p)
    leaking = n.mul(n.mul(q, t), p)
    return {'cut_transport_commutator': mismatch,
            'commutator_square': n.mul(mismatch, mismatch),
            'complementary_norm': n.mul(n.transpose(leaking), leaking),
            'retained_norm_loss': n.sub(p, n.mul(n.transpose(retained), retained))}


def first_exit(cosine, sine, maximum_index):
    """Monitored P/Q first exit: norm of Q T (P T P)^k P, k starts at 0."""
    weights = branch_weights(cosine, sine)
    r, q = weights['continue'], weights['exit']
    maximum_index = steps(maximum_index)
    probabilities = tuple(q*(r**k) for k in range(maximum_index+1))
    tail = r**(maximum_index+1)
    if sum(probabilities)+tail != 1:
        raise ArithmeticError('The native first-exit ledger did not close')
    return {'continuation': r, 'exit': q, 'probabilities': probabilities,
            'unresolved_tail': tail, 'total_mass_with_tail': sum(probabilities)+tail,
            'terminates_almost_surely': r < 1,
            'mean_continuations': r/q if q else None,
            'variance_continuations': r/(q*q) if q else None}


def varying_first_exit(rotations):
    alive, exits = Q(1), []
    for c, s in rotations:
        w = branch_weights(c, s)
        exits.append(alive*w['exit'])
        alive *= w['continue']
    return {'probabilities': tuple(exits), 'unresolved_tail': alive,
            'total_mass_with_tail': sum(exits)+alive}


def coherent_vs_monitored(cosine, sine, count):
    t = rotation(cosine, sine)
    count = steps(count)
    c, s = t[0][0], t[1][0]
    p = n.matrix(((1, 0), (0, 0)))
    complement = n.sub(n.identity(2), p)
    full_two = n.mul(t, t)
    reset_two = n.mul(n.mul(n.mul(p, t), p), n.mul(n.mul(p, t), p))
    return {'coherent_visible_amplitude': power(t, count)[0][0],
            'monitored_survival_amplitude': c**count,
            'coherent_final_visible_weight': power(t, count)[0][0]**2,
            'monitored_survival_weight': c**(2*count),
            'two_step_compression_defect': n.sub(n.mul(n.mul(p, full_two), p), reset_two),
            'two_step_excursion': n.mul(n.mul(n.mul(n.mul(p, t), complement), t), p),
            'two_step_excursion_coefficient': -s*s}


def event_information(cosine, sine):
    """Native angular tangent c'=-s, s'=c; no Gaussian metric calibration."""
    weights = branch_weights(cosine, sine)
    c, s = n.scalar(cosine), n.scalar(sine)
    r, q = weights['continue'], weights['exit']
    if not 0 < r < 1:
        raise ValueError('Regular event Fisher information requires both outcomes to have positive weight')
    dr = -2*c*s
    binary = dr*dr*(1/r+1/q)
    waiting_r = 1/(r*q*q)
    return {'binary_fisher_angular': binary, 'waiting_fisher_continuation': waiting_r,
            'waiting_fisher_angular': dr*dr*waiting_r,
            'expected_trials_to_exit': 1/q,
            'information_per_expected_trial': dr*dr*waiting_r*q,
            'likelihood': 'native norm-weighted monitored event record; not Gaussian'}


def native_binary_response(amplitude):
    """Rational native meter: rotate the balanced pair by Cayley(amplitude/4)."""
    a = n.scalar(amplitude)
    data = cayley_rotation(a/4)
    c, s = data['cosine'], data['sine']
    p0, p1 = (c-s)**2/2, (c+s)**2/2
    return {'probabilities': (p0, p1), 'contrast_mean': p1-p0,
            'contrast_variance': 1-(p1-p0)**2,
            'origin_probability_derivatives': (Q(-1, 2), Q(1, 2)),
            'origin_contrast_derivative': Q(1)}


def native_binary_score_metric(analysis_rows):
    """Product balanced native meters give D^T D at the zero contrast state.

Each meter has outcomes +/-1, norm weights 1/2 and derivatives +/-D_j/2.
Their product preparation determines independence and covariance I; no
Gaussian likelihood or externally calibrated noise precision is used.
"""
    rows = o.rect(analysis_rows)
    size = len(rows[0])
    form = n.scale(n.identity(size), 0)
    for row in rows:
        for sign in (-1, 1):
            dp = tuple(sign*x/2 for x in row)
            contribution = n.matrix([[x*y/Q(1, 2) for y in dp] for x in dp])
            form = n.add(form, contribution)
    return {'metric': form, 'record_covariance_at_origin': n.identity(len(rows)),
            'outcome_probabilities_at_origin': (Q(1, 2), Q(1, 2)),
            'Gaussian_assumption': False, 'product_native_meter_preparation': True,
            'likelihood': 'native balanced binary norm readout'}


def native_binary_seam(cosine, sine, normal_coordinate):
    """Non-Gaussian realization of R12's completed seam selection equation.

The monitored rotational gate sets r=c^2. The same complementary gain s
drives the paired sheet record. At an anchored zero contrast, independent
native balanced binary meters generate the seed Fisher form. The R12 Stein
and derivative routines are reused as algebra, not as a Gaussian likelihood.
"""
    w = branch_weights(cosine, sine)
    if not 0 < w['continue'] < 1:
        raise ValueError('The completed binary seam experiment requires a regular stopping law')
    v, s = n.scalar(normal_coordinate), n.scalar(sine)
    geometry = o.seam_foundation(s, v, w['continue'])
    meters = native_binary_score_metric(geometry['experiment'].readout)
    derived = geometry['observer_metric']
    if meters['metric'] != derived['seed']:
        raise ArithmeticError('The native binary seed did not match the completed selection equation')
    return {'derived_continuation': w['continue'], 'derived_complement_gain': s,
            'derived_kappa': geometry['kappa'], 'seed_information': meters['metric'],
            'record_covariance_at_origin': meters['record_covariance_at_origin'],
            'future_observer': derived['observer'], 'response_metric': derived['metric'],
            'untagged_response_metric_at_origin': derived['seed'],
            'discarded_tag_information': n.sub(derived['metric'], derived['seed']),
            'tangent_metric': geometry['metric_jet'].value,
            'tangent_metric_first': geometry['metric_jet'].first,
            'tangent_metric_second': geometry['metric_jet'].second,
            'Gaussian_assumption': False, 'external_record_noise_precision': False,
            'likelihood': 'tagged product of native balanced binary meters',
            'gate_information': event_information(cosine, sine),
            'gate_response_cross_score_at_origin': (Q(0),)*4,
            'state_scope': 'anchored zero contrast',
            'protocol_identification': 'sheet gain is the monitored complementary sine leg'}


def r12_event_bridge(cosine, sine, normal_coordinate):
    """Derive R12 continuation and match its record gain to the complementary leg.

This is an additional explicit identification of protocols. R12's Gaussian
calibration is retained, not derived from the binary hidden preparation.
"""
    w = branch_weights(cosine, sine)
    if not 0 < w['continue'] < 1:
        raise ValueError('The R12 bridge requires a nondegenerate stopping law')
    report = o.seam_foundation(n.scalar(sine), n.scalar(normal_coordinate), w['continue'])
    return {'derived_continuation': w['continue'], 'derived_record_gain': n.scalar(sine),
            'derived_kappa': report['kappa'], 'tangent_metric': report['metric_jet'].value,
            'Gaussian_calibration_retained_from_R12': True,
            'Gaussian_calibration_derived_from_complement': False,
            'protocol_identification': 'R12 record gain equals native QTP sine amplitude'}


def exact_witnesses():
    c, s, x, z = Q(3, 5), Q(4, 5), Q(5, 13), Q(12, 13)
    noise = endogenous_noise(c, s, x, z, 4)
    history = reduced_history(rotation(c, s), 1, ((x,),), ((z,),), 4)
    return {'rotation': rotation(c, s), 'unit_preparation': noise['preparation'],
            'visible_history': history['visible_states'], 'endogenous_noise': noise,
            'memory_equation_residuals': history['residuals'],
            'first_exit': first_exit(c, s, 4), 'event_information': event_information(c, s),
            'cut_incompatibility': cut_incompatibility(c, s),
            'coherent_vs_monitored': coherent_vs_monitored(c, s, 2),
            'native_binary_seam': native_binary_seam(c, s, Q(1, 3)),
            'R12_stopping_and_gain_bridge': r12_event_bridge(c, s, Q(1, 3))}
