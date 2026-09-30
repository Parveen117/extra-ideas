"""Record the exact R3 checks and verify the unchanged R1/R2 evidence."""

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest

from seam_bond import rational_seam_return


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    here = Path(__file__).resolve().parent
    root = here.parent
    suite = unittest.defaultTestLoader.discover(str(here), 'test_seam_bond.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    preserved = {}
    for revision in ('R1', 'R2'):
        earlier = json.loads((here / f'{revision}_VERIFICATION.json').read_text())
        preserved[revision] = {
            'recorded_status': earlier['status'],
            'files_checked': len(earlier['sha256']),
            'all_recorded_hashes_match': all(
                (root / name).is_file() and digest(root / name) == expected
                for name, expected in earlier['sha256'].items()),
        }
    example = rational_seam_return(2, 3, Q(1, 9))
    paths = [
        '4ways.tex',
        '01-intrinsic-numbers/SEAM_COMPLEX_COORDINATES_R3.md',
        '02-relational-response/TYPED_SEAM_R3.md',
        '03-lambda-reference/SEAM_BOND_COMPLEX_STRUCTURE_R3.md',
        '04-operator-evolution/aghora_return.py',
        '04-operator-evolution/seam_bond.py',
        '04-operator-evolution/test_seam_bond.py',
        '04-operator-evolution/verify_r3.py',
    ]
    passed = result.wasSuccessful() and all(
        entry['all_recorded_hashes_match'] for entry in preserved.values())
    record = {
        'status': 'PASS_EXACT_FINITE_CHECKS' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'arithmetic': 'fractions.Fraction; exact rational matrices',
        'tests_run': result.testsRun,
        'failures': len(result.failures),
        'errors': len(result.errors),
        'scope': 'Exact two-mode examples, coordinate covariance, and hypothesis counterexamples. General identities are proved in SEAM_BOND_COMPLEX_STRUCTURE_R3.md. No approximate exponential or physical experiment.',
        'added_bond_rule': 'range(B)=L and kernel(B)=A L, with V=L direct-sum A L',
        'bond_unique_under_added_rule': True,
        'complex_structure_identity': '[A,B]^2=-I',
        'relative_metric_unique_up_to_scale_under_added_conditions': True,
        'compatible_generator': 'G=omega[A,B]; omega remains free',
        'unique_return_sign_selected': False,
        'universal_lambda_selected': False,
        'source_full_closure_identified': False,
        'physical_validation': False,
        'example': {
            'inputs_are_chosen': True,
            'gamma': '2', 'omega': '3', 'cayley_parameter': '1/9',
            'compressed_coefficient': str(example['mu']),
            'retained_squared_norm_fraction': str(example['mu']**2),
            'removed_squared_norm_fraction': str(example['defect_coefficient']),
        },
        'earlier_verification_integrity': preserved,
        'sha256': {name: digest(root / name) for name in paths},
    }
    (here / 'R3_VERIFICATION.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f"R3 record: {record['status']}; earlier recorded file hashes preserved: "
          f"{all(entry['all_recorded_hashes_match'] for entry in preserved.values())}")
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
