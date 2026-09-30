"""Run the focused exact R2 checks and record scope and file hashes."""

import hashlib
import json
from pathlib import Path
import sys
import unittest

from aghora_return import compressed_return, orthogonal_example


def main():
    here = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(str(here), 'test_aghora_return.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    _, _, _, returned, bond = orthogonal_example()
    compressed, defect, _ = compressed_return(returned, bond)
    paths = [here / 'aghora_return.py', here / 'test_aghora_return.py', here / 'verify_r2.py',
             here.parent / '03-lambda-reference' / 'AGHORA_RETURN_R2.md']
    record = {
        'status': 'PASS_EXACT_FINITE_CHECKS' if result.wasSuccessful() else 'FAIL',
        'python': sys.version.split()[0],
        'arithmetic': 'fractions.Fraction; exact rational matrices',
        'tests_run': result.testsRun,
        'failures': len(result.failures),
        'errors': len(result.errors),
        'scope': 'Exact finite examples and hypothesis counterexamples. General exponential-flow and projection identities are proved in AGHORA_RETURN_R2.md. No numerical exponential approximation or physical experiment.',
        'return_eigenvalue_constraint': 'r^2=1',
        'unique_return_sign_selected': False,
        'universal_lambda_selected': False,
        'example_inputs_are_chosen': True,
        'orthogonal_example_compressed_coefficient': str(compressed[0][0]),
        'orthogonal_example_squared_norm_defect': str(defect[0][0]),
        'sha256': {str(p.relative_to(here.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in paths},
    }
    (here / 'R2_VERIFICATION.json').write_text(json.dumps(record, indent=2) + '\n')
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
