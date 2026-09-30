"""Run the focused exact checks and write a reproducible verification record."""

import hashlib
import json
from pathlib import Path
import sys
import unittest

from response_transport import triangle


def main():
    here = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(str(here), 'test_response_transport.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    paths = [here / 'response_transport.py', here / 'test_response_transport.py',
             here / 'verify_r1.py',
             here.parent / '03-lambda-reference' / 'RETURN_INVARIANTS_R1.md']
    net = triangle()
    changed = net.rescale({'A': 5, 'B': 11, 'C': 13})
    record = {
        'status': 'PASS_EXACT_FINITE_CHECKS' if result.wasSuccessful() else 'FAIL',
        'python': sys.version.split()[0],
        'arithmetic': 'fractions.Fraction; exact rational operations',
        'tests_run': result.testsRun,
        'failures': len(result.failures),
        'errors': len(result.errors),
        'scope': 'Finite reproducible examples and counterexamples. General claims are proved in RETURN_INVARIANTS_R1.md; this run is not a formal proof assistant or physical validation.',
        'example_inputs_are_chosen': True,
        'triangle_return': str(net.holonomy(('A', 'B', 'C', 'A'))),
        'rescaled_triangle_return': str(changed.holonomy(('A', 'B', 'C', 'A'))),
        'triangle_cycle_rank': net.cycle_rank,
        'universal_lambda_selected': False,
        'sha256': {str(p.relative_to(here.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in paths},
    }
    (here / 'R1_VERIFICATION.json').write_text(json.dumps(record, indent=2) + '\n')
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
