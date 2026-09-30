"""Check uncertain profile recovery and preserve R1-R5 verification evidence."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import unittest

from interval_profile_recovery import predict_interval, recover_intervals, serializable


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', type=Path, required=True)
    args = parser.parse_args()
    os.environ['RKF_SOURCE_ROOT'] = str(args.rkf_root.resolve())
    import test_interval_profile_recovery as tests
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.IntervalProfileRecoveryTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    here = Path(__file__).resolve().parent
    root = here.parent
    previous = {}
    for revision in range(1, 6):
        record = json.loads((here/f'R{revision}_VERIFICATION.json').read_text())
        previous[f'R{revision}'] = {
            'recorded_status': record['status'],
            'files_checked': len(record['sha256']),
            'all_recorded_hashes_match': all(
                (root/path).is_file() and digest(root/path) == wanted
                for path, wanted in record['sha256'].items()),
        }
    packet = getattr(tests.IntervalProfileRecoveryTests, 'packet', {})
    native = getattr(tests.IntervalProfileRecoveryTests, 'native', {})
    passed = result.wasSuccessful() and result.testsRun == 10 and bool(packet) and all(
        item['all_recorded_hashes_match'] and item['recorded_status'].startswith('PASS')
        for item in previous.values())
    examples = {}
    for name, last in (('base', 216), ('changed', 180)):
        inputs = tests.observations(last)
        recovered = recover_intervals(inputs)
        examples[name] = {
            'coefficient_intervals': inputs,
            'recovery': recovered,
            'probe': '1/10', 'first_unrecovered_product_upper_bound': '3',
            'prediction_interval': predict_interval(recovered['cell_intervals'], '1/10', 3),
            'completed_native_response_interval': native.get(name, {}).get('interval'),
        }
    paths = ['4ways.tex', '03-lambda-reference/INTERVAL_PROFILE_RECOVERY_R6.md',
             '04-operator-evolution/interval_profile_recovery.py',
             '04-operator-evolution/test_interval_profile_recovery.py',
             '04-operator-evolution/verify_r6.py',
             '04-operator-evolution/profile_recovery.py',
             '04-operator-evolution/native_profile_probe.cjs',
             '04-operator-evolution/native_bond_bridge.py',
             '04-operator-evolution/R4_SOURCE_PINS.json',
             '04-operator-evolution/R5_SOURCE_PINS.json'] + [
                 f'04-operator-evolution/R{revision}_VERIFICATION.json' for revision in range(1, 6)]
    record = {
        'status': 'PASS_CERTIFIED_INTERVAL_PROFILE_CHECKS' if passed else 'FAIL',
        'python': sys.version.split()[0], 'node': packet.get('runtime'),
        'arithmetic': 'exact rational interval endpoints and unchanged native forward solver',
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'native_prediction_corner_comparison_count': 8,
        'native_query_count': len(packet.get('results', [])),
        'upstream_runtime_sha256': packet.get('checked_runtime_sha256', {}),
        'coefficient_error_examples': examples,
        'undecided_depth_example': recover_intervals([(3, 3), (-18, -18), (107, 109)]),
        'incompatible_depth_example': recover_intervals([(3, 3), (-18, -18), (100, 107)]),
        'native_prediction_corners': {key: value['value'] for key, value in native.items()
                                      if key.startswith('corner_')},
        'supplied_coefficient_error_bounds_propagated': True,
        'experimental_coefficient_error_bounds_established': False,
        'unknown_calibration_gains_removed': False,
        'whole_source_selected': False,
        'universal_lambda_selected': False,
        'physical_experiment_performed': False,
        'earlier_verification_integrity': previous,
        'sha256': {path: digest(root/path) for path in paths},
    }
    (here/'R6_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R6: {record['status']}; {result.testsRun} tests; 8 native corner comparisons")
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
