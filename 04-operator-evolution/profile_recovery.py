"""R5: exact depth-prefix recovery from calibrated weak-response coefficients.

Input c[k] is the coefficient of z**(2*k+1) in x_0(z), NOT its
unscaled derivative and NOT samples near the closed cut z=1. The native
paired products are multiplied by the same known z. Continued-fraction
inversion is classical; the proof and scope are in PROFILE_RECOVERY_R5.md.
"""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess

from native_bond_bridge import check_runtime_pins


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Q)):
        raise TypeError('Use integers, rational strings or Fraction, not floats')
    return Q(value)


def positive_cells(values, *, allow_empty=False):
    values = tuple(exact(value) for value in values)
    if (not values and not allow_empty) or any(value <= 0 for value in values):
        raise ValueError('Positive paired cell products are required')
    return values


def reciprocal_jet(coefficients):
    """Reciprocal modulo t**n; retain exactly the n supplied coefficients."""
    values = tuple(exact(value) for value in coefficients)
    if not values or not values[0]:
        raise ValueError('A nonzero constant coefficient is required')
    answer = [1 / values[0]]
    for k in range(1, len(values)):
        answer.append(-sum(values[i] * answer[k-i]
                           for i in range(1, k+1)) / values[0])
    return tuple(answer)


def recover_prefix(coefficients):
    """Recover m positive products from m odd coefficients; tail stays unknown.

    Each reciprocal-and-shift consumes one coefficient. A nonpositive
    stripped constant rejects a positive profile of the requested depth.
    In particular, an apparent zero tail is not silently called an infinite
    strictly positive profile or evidence for physical termination.
    """
    current = tuple(exact(value) for value in coefficients)
    if not current:
        raise ValueError('At least one odd response coefficient is required')
    recovered = []
    while current:
        if current[0] <= 0:
            raise ValueError(f'Nonpositive recovered product at depth {len(recovered)}')
        recovered.append(current[0])
        current = reciprocal_jet(current)[1:] if len(current) > 1 else ()
    return tuple(recovered)


def finite_fraction(cells):
    """N(t)/D(t)=x(z)/z for a FINITE zero-tail aperture, t=z**2.

    This polynomial companion computes formal coefficients. Real responses
    and tail enclosures are checked through the unchanged native solver.
    Coefficients are in ascending order; no cancellation is needed.
    """
    cells = positive_cells(cells)
    numerator, denominator = (Q(0),), (Q(1),)
    for cell in reversed(cells):
        new_numerator = tuple(cell * value for value in denominator)
        new_denominator = list(denominator) + [Q(0)] * max(
            0, len(numerator) + 1 - len(denominator))
        for i, value in enumerate(numerator):
            new_denominator[i+1] += cell * value
        while len(new_denominator) > 1 and new_denominator[-1] == 0:
            new_denominator.pop()
        numerator, denominator = new_numerator, tuple(new_denominator)
    return numerator, denominator


def finite_jet(cells, count):
    """Finite-aperture odd coefficients; arbitrary deeper tails need count cells.

    If len(cells) >= count, these count coefficients are shared by every
    positive continuation. Otherwise they describe ONLY the zero-tail model.
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError('The coefficient count must be a positive integer')
    numerator, denominator = finite_fraction(cells)
    coefficients = []
    for k in range(count):
        value = numerator[k] if k < len(numerator) else Q(0)
        value -= sum(denominator[i] * coefficients[k-i]
                     for i in range(1, min(k+1, len(denominator))))
        coefficients.append(value / denominator[0])
    return tuple(coefficients)


def first_difference(prefix, left_next, right_next):
    """First differing odd coefficient after a common prefix (left - right)."""
    prefix = positive_cells(prefix, allow_empty=True)
    left_next, right_next = positive_cells((left_next, right_next))
    if left_next == right_next:
        raise ValueError('The first unequal cell products must differ')
    coefficient = left_next - right_next
    for cell in prefix:
        coefficient *= -cell**2
    return {'power': 2*len(prefix)+1, 'coefficient': coefficient}


def weak_tail_bound(prefix, probe, tail_cell_bound):
    """Upper bound for variation across tails with their first cell <= Q.

    The deeper tails must be nonnegative. This is a forward uncertainty
    bound, not a stability bound for recovery from noisy coefficients.
    """
    prefix = positive_cells(prefix)
    probe, tail_cell_bound = positive_cells((probe, tail_cell_bound))
    bound = probe * tail_cell_bound
    for cell in prefix:
        bound *= (cell * probe)**2
    return bound


def native_queries(rkf_root, requests):
    """Hash-check the pinned runtime before executing its original solver."""
    checked = check_runtime_pins(rkf_root)
    completed = subprocess.run(
        ['node', str(Path(__file__).with_name('native_profile_probe.cjs')),
         str(Path(rkf_root).resolve())],
        input=json.dumps(requests), check=True, capture_output=True,
        text=True, timeout=30)
    packet = json.loads(completed.stdout)
    packet['checked_runtime_sha256'] = checked
    return packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coefficients_json', help='Example: \'[3, -18, 216]\'')
    args = parser.parse_args()
    values = json.loads(args.coefficients_json)
    if not isinstance(values, list):
        parser.error('Supply a JSON array of exact coefficients')
    try:
        recovered = recover_prefix(values)
    except (ValueError, TypeError, ZeroDivisionError) as error:
        parser.error(str(error))
    print(json.dumps({
        'paired_cell_prefix': list(map(str, recovered)),
        'last_input_power_of_z': 2*len(recovered)-1,
        'probe_calibration_required': True,
        'unresolved_tail': True,
        'infinite_periodicity_certified': False,
        'universal_constant_selected': False,
    }, indent=2))


if __name__ == '__main__':
    main()
