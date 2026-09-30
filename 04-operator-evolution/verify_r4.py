"""Verify R4 against the unchanged source checkout and preserve earlier hashes."""

import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import sys
import unittest

from native_bond_bridge import closure_budget


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strings(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: strings(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [strings(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    args = parser.parse_args()
    os.environ['RKF_SOURCE_ROOT'] = str(args.rkf_root.resolve())
    here = Path(__file__).resolve().parent
    root = here.parent
    import test_native_bond_bridge as tests
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.NativeBondBridgeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    packet = getattr(tests.NativeBondBridgeTests, 'packet', {})
    prior = {}
    for revision in ('R1', 'R2', 'R3'):
        earlier = json.loads((here/f'{revision}_VERIFICATION.json').read_text())
        prior[revision] = {
            'recorded_status': earlier['status'],
            'files_checked': len(earlier['sha256']),
            'all_recorded_hashes_match': all(
                (root/name).is_file() and sha256(root/name) == expected
                for name, expected in earlier['sha256'].items()),
        }
    passed = (result.wasSuccessful() and result.testsRun == 12 and bool(packet)
              and all(entry['all_recorded_hashes_match'] for entry in prior.values()))
    budgets = {}
    for name, prefix in packet.get('prefixes', {}).items():
        interval = prefix['interval']
        budgets[name] = {
            'cells': prefix['count'],
            'source_interval': [interval['lower'], interval['upper']],
            'actual_zero_tail': prefix['actual_zero_tail'],
            **strings(closure_budget(Q(interval['lower']), Q(interval['upper']))),
        }
    paths = [
        '4ways.tex', 'CROSS_REPO_LINEAGE.md',
        '03-lambda-reference/NATIVE_BOND_BRIDGE_R4.md',
        '04-operator-evolution/aghora_return.py',
        '04-operator-evolution/seam_bond.py',
        '04-operator-evolution/native_bond_bridge.py',
        '04-operator-evolution/native_bond_probe.cjs',
        '04-operator-evolution/test_native_bond_bridge.py',
        '04-operator-evolution/verify_r4.py',
        '04-operator-evolution/R4_SOURCE_PINS.json',
    ]
    record = {
        'status': 'PASS_EXACT_SOURCE_BRIDGE_CHECKS' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': packet.get('runtime'),
        'tests_run': result.testsRun,
        'failures': len(result.failures), 'errors': len(result.errors),
        'arithmetic': 'fractions.Fraction and unchanged native exact cut arithmetic',
        'upstream_runtime_sha256': packet.get('checked_runtime_sha256', {}),
        'native_presentation_sha256': packet.get('presentation_sha256'),
        'native_polynomial_replay_count': sum(
            len(case['proofs']) for case in packet.get('algebra_checks', [])),
        'native_polynomial_replays': packet.get('algebra_checks', []),
        'source_schur_proofs_replayed': packet.get('source_schur_proof_names', []),
        'finite_budgets': budgets,
        'bridge_scope': 'The supplied positive paired-depth source, A=K, and its normalized boundary response; exact completed cut under the existing source criterion a=b+1.',
        'mirror_complement_follows_in_this_source_realization': True,
        'native_iota_newly_derived': False,
        'normalization_certifies_raw_response_closure': False,
        'finite_exact_cut_certifies_completed_cut': False,
        'universal_source_family_selected': False,
        'universal_lambda_selected': False,
        'dynamical_return_sign_selected': False,
        'physical_validation': False,
        'earlier_verification_integrity': prior,
        'sha256': {name: sha256(root/name) for name in paths},
    }
    (here/'R4_VERIFICATION.json').write_text(json.dumps(record, indent=2)+'\n')
    print(f"R4 record: {record['status']}; native polynomial replays: "
          f"{record['native_polynomial_replay_count']}")
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
