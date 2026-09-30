"""R6: certified depth recovery from supplied coefficient error intervals.

Intervals enclose calibrated coefficients of z, z**3, ... at zero coupling.
Their experimental validity is a supplied contract. All endpoints are exact
rationals. Correlations may widen the answers; they never justify shrinking
an enclosure or replacing an unresolved sign by its midpoint.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
import json

from profile_recovery import exact


@dataclass(frozen=True)
class Interval:
    lower: Q
    upper: Q

    def __post_init__(self):
        object.__setattr__(self, 'lower', exact(self.lower))
        object.__setattr__(self, 'upper', exact(self.upper))
        if self.lower > self.upper:
            raise ValueError('Interval endpoints must be ordered')

    @staticmethod
    def of(value):
        return value if isinstance(value, Interval) else Interval(value, value)

    def __add__(self, other):
        other = Interval.of(other)
        return Interval(self.lower+other.lower, self.upper+other.upper)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + (-Interval.of(other))

    def __mul__(self, other):
        other = Interval.of(other)
        values = [a*b for a in (self.lower, self.upper)
                  for b in (other.lower, other.upper)]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lower <= 0 <= self.upper:
            raise ZeroDivisionError('An interval containing zero has no finite reciprocal')
        return Interval(1/self.upper, 1/self.lower)

    def contains(self, value):
        return self.lower <= exact(value) <= self.upper

    @property
    def width(self):
        return self.upper-self.lower


def boxes(values):
    return tuple(value if isinstance(value, Interval) else Interval(*value)
                 for value in values)


def reciprocal_intervals(coefficients):
    current = boxes(coefficients)
    if not current:
        raise ValueError('At least one coefficient interval is required')
    inverse = current[0].reciprocal()
    result = [inverse]
    for k in range(1, len(current)):
        result.append(-sum(current[i]*result[k-i] for i in range(1, k+1))*inverse)
    return tuple(result)


def recover_intervals(coefficient_boxes):
    """Enclose every admissible prefix; stop before a nonpositive/unknown pivot.

    CERTIFIED_PREFIX: all supplied coefficient realizations have positive
    products through the requested depth. UNRESOLVED_POSITIVITY: enclosure
    crosses zero; a compatible profile may exist. INCOMPATIBLE_POSITIVE_PROFILE:
    the attempted product is nonpositive for every realization in the box.
    Neither stopping status asserts that a physical chain terminates.
    """
    current = boxes(coefficient_boxes)
    if not current:
        raise ValueError('At least one coefficient interval is required')
    count, recovered = len(current), []
    while current:
        candidate = current[0]
        if candidate.lower <= 0:
            return {
                'status': ('INCOMPATIBLE_POSITIVE_PROFILE' if candidate.upper <= 0
                           else 'UNRESOLVED_POSITIVITY'),
                'cell_intervals': tuple(recovered),
                'requested_depth': count,
                'stopped_depth': len(recovered),
                'attempted_cell_interval': candidate,
                'unresolved_tail': True,
            }
        recovered.append(candidate)
        current = reciprocal_intervals(current)[1:] if len(current) > 1 else ()
    return {'status': 'CERTIFIED_PREFIX', 'cell_intervals': tuple(recovered),
            'requested_depth': count, 'stopped_depth': None,
            'attempted_cell_interval': None, 'unresolved_tail': True}


def predict_interval(cell_boxes, probe, tail_cell_bound):
    """Boundary response enclosure, including prefix uncertainty and hidden tail.

    The tail's first cell is at most Q=tail_cell_bound; deeper responses are
    nonnegative, so its response is in [0,zQ]. The result is exact for the
    independent product box and that tail interval, and conservative for
    correlated products returned by recovery. Real forward corners are
    independently checked with the original native solver.
    """
    cells = boxes(cell_boxes)
    z, q = exact(probe), exact(tail_cell_bound)
    if not cells or z <= 0 or q <= 0 or any(cell.lower <= 0 for cell in cells):
        raise ValueError('Require a positive prefix, probe and supplied tail bound')
    response = Interval(0, z*q)
    for cell in reversed(cells):
        response = Interval(
            z*cell.lower/(1+z*cell.lower*response.upper),
            z*cell.upper/(1+z*cell.upper*response.lower))
    return response


def prediction_corners(cell_boxes, probe, tail_cell_bound):
    """Two alternating corners attaining the box range, for native replay."""
    cells = boxes(cell_boxes)
    predict_interval(cells, probe, tail_cell_bound)  # validate the same contract
    bound = exact(probe)*exact(tail_cell_bound)
    return {
        'lower': {'cells': tuple(c.lower if j % 2 == 0 else c.upper
                                 for j, c in enumerate(cells)),
                  'tail': Q(0) if len(cells) % 2 == 0 else bound},
        'upper': {'cells': tuple(c.upper if j % 2 == 0 else c.lower
                                 for j, c in enumerate(cells)),
                  'tail': bound if len(cells) % 2 == 0 else Q(0)},
    }


def serializable(value):
    if isinstance(value, Interval):
        return [str(value.lower), str(value.upper)]
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coefficient_boxes_json', help='JSON array of [lower,upper] pairs')
    parser.add_argument('--probe', help='Optional positive common product multiplier')
    parser.add_argument('--tail-bound', help='Supplied upper bound on the first unseen product')
    args = parser.parse_args()
    try:
        result = recover_intervals(json.loads(args.coefficient_boxes_json))
        if bool(args.probe) != bool(args.tail_bound):
            raise ValueError('Supply --probe and --tail-bound together')
        if args.probe:
            if result['status'] == 'INCOMPATIBLE_POSITIVE_PROFILE':
                raise ValueError('Incompatible data have no admitted positive profile to predict')
            result['predicted_response_interval'] = predict_interval(
                result['cell_intervals'], args.probe, args.tail_bound)
            result['prediction_prefix_depth'] = len(result['cell_intervals'])
        result['coefficient_error_bounds_supplied'] = True
        result['physical_experiment_performed'] = False
        print(json.dumps(serializable(result), indent=2))
    except (ValueError, TypeError, ZeroDivisionError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
