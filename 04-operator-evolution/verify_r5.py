"""Verify R5 profile recovery against the unchanged native paired-depth source."""

import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import sys
import unittest

from profile_recovery import finite_jet, first_difference, recover_prefix, weak_tail_bound


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


def outward_decimal(interval, digits=9):
    """Round outward using only integer division, retaining trailing zeros."""
    scale = 10**digits
    lower, upper = Q(interval['lower']), Q(interval['upper'])
    lo = (lower.numerator*scale)//lower.denominator
    hi = -((-upper.numerator*scale)//upper.denominator)
    if lo < 0:
        raise ValueError('This formatter is for the positive witness intervals')
    return [f'{value//scale}.{value%scale:0{digits}d}' for value in (lo, hi)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    args = parser.parse_args()
    os.environ['RKF_SOURCE_ROOT'] = str(args.rkf_root.resolve())
    here, root = Path(__file__).resolve().parent, Path(__file__).resolve().parent.parent
    import test_profile_recovery as tests
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.ProfileRecoveryTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    packet = getattr(tests.ProfileRecoveryTests, 'packet', {})
    native = getattr(tests.ProfileRecoveryTests, 'native', {})
    prior = {}
    for revision in ('R1', 'R2', 'R3', 'R4'):
        earlier = json.loads((here/f'{revision}_VERIFICATION.json').read_text())
        prior[revision] = {
            'recorded_status': earlier['status'],
            'files_checked': len(earlier['sha256']),
            'all_recorded_hashes_match': all(
                (root/name).is_file() and sha256(root/name) == expected
                for name, expected in earlier['sha256'].items()),
        }
    passed = (result.wasSuccessful() and result.testsRun == 14 and bool(packet)
              and all(item['all_recorded_hashes_match'] and
                      item['recorded_status'].startswith('PASS') for item in prior.values()))
    heldout, budgets = {}, []
    for name in ('base', 'changed'):
        if f'{name}_heldout' in native:
            interval = native[f'{name}_heldout']['interval']
            heldout[name] = {'probe': '1/2', 'interval': interval,
                             'outward_decimal_interval': outward_decimal(interval)}
    for i, (prefix, probe) in enumerate(tests.BOUND_CASES):
        if f'bound_{i}' in native:
            budgets.append({'prefix': prefix, 'probe': probe,
                            'tail_first_cell_upper_bound': '3',
                            'native_interval': native[f'bound_{i}'],
                            'weak_variation_bound': weak_tail_bound(prefix, probe, 3)})
    paths = [
        '4ways.tex', 'CROSS_REPO_LINEAGE.md',
        '03-lambda-reference/PROFILE_RECOVERY_R5.md',
        '04-operator-evolution/aghora_return.py',
        '04-operator-evolution/native_bond_bridge.py',
        '04-operator-evolution/R4_SOURCE_PINS.json',
        '04-operator-evolution/profile_recovery.py',
        '04-operator-evolution/native_profile_probe.cjs',
        '04-operator-evolution/test_profile_recovery.py',
        '04-operator-evolution/verify_r5.py',
        '04-operator-evolution/R5_SOURCE_PINS.json',
    ] + [f'04-operator-evolution/{revision}_VERIFICATION.json'
         for revision in ('R1', 'R2', 'R3', 'R4')]
    record = {
        'status': 'PASS_EXACT_PROFILE_RECOVERY_CHECKS' if passed else 'FAIL',
        'python': sys.version.split()[0], 'node': packet.get('runtime'),
        'tests_run': result.testsRun, 'failures': len(result.failures),
        'errors': len(result.errors),
        'arithmetic': 'fractions.Fraction and unchanged native exact cut arithmetic',
        'upstream_runtime_sha256': packet.get('checked_runtime_sha256', {}),
        'native_query_count': len(packet.get('results', [])),
        'native_finite_polynomial_comparison_count': len(tests.PROFILES)*len(tests.PROBES),
        'native_forward_values': {key: value['value'] for key, value in native.items()
                                  if key.startswith('poly_')},
        'native_schur': native.get('native_schur'),
        'native_schur_proof_replay_count': len(native.get('native_schur', {}).get('proofs', {})),
        'input_coefficient_convention': 'c_k is the coefficient of z^(2k+1) at z=0, with calibrated common product multiplier z',
        'recovery_example': {'input_coefficients': [3, -18, 216],
                             'recovered_prefix': recover_prefix([3, -18, 216]),
                             'unresolved_tail': True},
        'equal_cut_witness': {
            'base_profile': '(3,2) repeated',
            'changed_profile': 'prefix (3,2), then (2,1) repeated',
            'base_coefficients': finite_jet([3, 2, 3], 3),
            'changed_coefficients': finite_jet([3, 2, 2], 3),
            'first_difference': first_difference([3, 2], 3, 2),
            'closure_basis': 'Original RKF SY1 exact tail synthesis, followed by the shared prefix; not inferred from a finite interval containing one',
            'native_tail_designs': {parity: native.get(f'design_{parity}')
                                   for parity in ('even', 'odd')},
            'native_prefix_closure_checks': {str(m): native.get(f'closed_prefix_{m}', {}).get('value')
                                             for m in range(1, 7)},
            'heldout_probe': heldout,
        },
        'unresolved_tail_budgets': budgets,
        'classical_inversion_newly_discovered': False,
        'whole_profile_selected_by_finite_jet_and_exact_cut': False,
        'noisy_observation_recovery_certified': False,
        'individual_directional_amplitudes_recovered': False,
        'universal_lambda_selected': False,
        'physical_validation': False,
        'earlier_verification_integrity': prior,
        'sha256': {name: sha256(root/name) for name in paths},
    }
    (here/'R5_VERIFICATION.json').write_text(json.dumps(strings(record), indent=2)+'\n')
    print(f"R5 record: {record['status']}; tests: {result.testsRun}; "
          f"native finite comparisons: {record['native_finite_polynomial_comparison_count']}")
    if heldout:
        print('Held-out outward intervals: ' + json.dumps({
            key: value['outward_decimal_interval'] for key, value in heldout.items()}))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
