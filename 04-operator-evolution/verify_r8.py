"""Verify R8 exact claims and unchanged canonical-source certificates.

Run with Python 3.11/3.12 and Node, against the pinned RKF checkout.
No archive source or canonical engine is installed/copied into this repo.
"""

import argparse
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_blob(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


def check_sources(root, files):
    checked = {}
    for item in files:
        source = root/item['path']
        if not source.is_file() or sha256(source) != item['sha256']:
            raise ValueError('Pinned RKF SHA256 mismatch: '+item['path'])
        if git_blob(source.read_bytes()) != item['git_blob_sha']:
            raise ValueError('Pinned RKF Git blob mismatch: '+item['path'])
        checked[item['path']] = item['sha256']
    return checked


def serializable(value):
    from fractions import Fraction
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {'/'.join(k) if isinstance(k, tuple) else k: serializable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(v) for v in value]
    return value


def preserved_evidence(root, here):
    result = {}
    paths = [here/f'R{i}_VERIFICATION.json' for i in range(1, 8)] + [here/'EMK_MASTER_REVIEW_VERIFICATION.json']
    for path in paths:
        record = json.loads(path.read_text())
        matches = {p: (root/p).is_file() and sha256(root/p) == expected for p, expected in record['sha256'].items()}
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha256(path),
                             'recorded_paths_checked': len(matches), 'all_hashes_match': all(matches.values())}
    return result


def mutation_controls(tests, c):
    """Run independent witness tests against six deliberately wrong variants."""
    original_report, original_recover = c.compression_report, c.recover_odd

    def wrong_excursion_sign(*args):
        report = original_report(*args)
        report['excursions'] = c.scale(report['excursions'], -1)
        report['identity_residual'] = c.sub(report['reduced'], c.sub(report['compressed_full'], report['excursions']))
        return report

    def omit_excursion(*args):
        report = original_report(*args)
        report['reduced'] = report['compressed_full']
        report['reduced_flat'] = c.is_zero(report['reduced'])
        return report

    def wrong_recovery(q_plus, q_minus):
        return c.add(c.scale(c.R, (q_plus-q_minus)/2), c.scale(c.RK, (q_plus+q_minus)/2))

    import emk_tensor_calculus as tc
    original_return = tc.return_report

    def alias_sheet(*args, **kwargs):
        report = original_return(*args, **kwargs)
        report['sheet_residue'] %= 2
        report['full_unledgered_return'] = report['full_carrier_return'] and report['sheet_residue'] == 0
        return report

    variants = [
        ('wrong_excursion_sign', c, 'compression_report', wrong_excursion_sign,
         'test_full_curved_reduced_flat_with_excursion_cancellation'),
        ('omit_excursion_memory', c, 'compression_report', omit_excursion,
         'test_full_curved_reduced_flat_with_excursion_cancellation'),
        ('wrong_odd_reconstruction_sign', c, 'recover_odd', wrong_recovery,
         'test_two_readouts_recover_every_declared_odd_target'),
        ('omit_connection_derivatives', c.ConnectionJet, 'curvature',
         lambda self, i, j: c.commutator(self.operators[i], self.operators[j]),
         'test_noncommuting_coefficients_can_have_zero_connection_curvature'),
        ('trace_zero_means_flat', c, 'is_zero', lambda a: c.trace(a) == 0,
         'test_native_nonzero_curvature_and_trace_blindness'),
        ('sheet_alias_mod_two', tc, 'return_report', alias_sheet,
         'test_existing_sheet_memory_survives_local_readout_flatness'),
    ]
    outcomes = {}
    for label, target, name, replacement, test_name in variants:
        original = getattr(target, name)
        try:
            setattr(target, name, replacement)
            stream = io.StringIO()
            suite = unittest.TestSuite([tests.CurvatureObservationTests(test_name)])
            with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                result = unittest.TextTestRunner(stream=stream).run(suite)
            outcomes[label] = {'rejected': len(result.failures) > 0 and len(result.errors) == 0,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(target, name, original)
    return outcomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--node', default='node')
    args = parser.parse_args()
    if sys.version_info[:2] not in ((3, 11), (3, 12)):
        parser.error('Use Python 3.11 or 3.12')
    here = Path(__file__).resolve().parent
    root = here.parent
    pins = json.loads((here/'R8_SOURCE_PINS.json').read_text())
    source_root = args.rkf_root.resolve()
    try:
        checked = check_sources(source_root, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    # All upstream byte checks precede imports and Node execution.
    os.environ['RKF_SOURCE_ROOT'] = str(source_root)
    import emk_curvature_observation as c
    import test_emk_curvature_observation as tests
    import emk_tensor_calculus as tc
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.CurvatureObservationTests)
    expected_test_count = suite.countTestCases()
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    mutations = mutation_controls(tests, c)
    command = [args.node, str(here/'r8_native_curvature_probe.cjs'), '--rkf-root', str(source_root),
               '--expected-input-sha256', pins['native_job_input_sha256'],
               '--expected-compression-sha256', pins['compression_contract_input_sha256']]
    completed = subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
    native = json.loads(completed.stdout)
    native_path = here/'R8_NATIVE_CERTIFICATE.json'
    native_path.write_text(json.dumps(native, indent=2)+'\n')
    native_observer = native['packet']['results'][2]['result']
    native_ok = (native['replay']['status'] == 'REPLAY_MATCH'
                 and native['replay']['verified_job']
                 and native['packet']['input_sha256'] == pins['native_job_input_sha256']
                 and native_observer['initial_rank'] == 2 and native_observer['completed_rank'] == 4
                 and native['compression_certificate']['input_sha256'] == pins['compression_contract_input_sha256']
                 and native['compression_certificate']['result']['status'] == 'EQUAL_IN_DECLARED_QUOTIENT'
                 and native['compression_certificate']['replay'] == 'REPLAY_MATCH'
                 and native['compression_certificate']['finite_carrier_assumed'] is False
                 and all(native['negative_controls'].values()))
    previous = preserved_evidence(root, here)
    witness = c.exact_witnesses()
    sheet = tc.Move('p', 'p', (('V', c.identity(2)),), 2)
    sheet_audit = tc.return_report(sheet, ((tc.Slot('V', 1),),))
    passed = (result.wasSuccessful() and result.testsRun == expected_test_count == 23
              and all(x['rejected'] for x in mutations.values()) and native_ok
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/CURVATURE_OBSERVATION_R8.md',
        '04-operator-evolution/emk_curvature_observation.py',
        '04-operator-evolution/test_emk_curvature_observation.py',
        '04-operator-evolution/r8_native_curvature_probe.cjs',
        '04-operator-evolution/verify_r8.py', '04-operator-evolution/R8_SOURCE_PINS.json',
        '04-operator-evolution/R8_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/aghora_return.py', '04-operator-evolution/emk_tensor_calculus.py',
    ] + [str(path.relative_to(root)) for path in sorted(here.glob('R[1-7]_VERIFICATION.json'))] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R8', 'date': '2026-10-01',
        'status': 'PASS_EXACT_CURVATURE_OBSERVATION_CERTIFICATES' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'arithmetic': 'exact rational matrices; unchanged canonical cut-field symbolic engine',
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'compression_cases': 100, 'odd_target_reconstructions': 25,
                               'source_curvature_comparisons': 27, 'constant_response_cases': 81},
        'mutation_controls': mutations, 'native_negative_controls': native['negative_controls'],
        'native_job': {'input_sha256': pins['native_job_input_sha256'], 'replay': native['replay'],
                       'equality_proofs_replayed': 2, 'observer': native_observer,
                       'engine_home': pins['canonical_engine_home'], 'source_reused_unchanged': True},
        'general_compression_proof': {
            'input_sha256': pins['compression_contract_input_sha256'],
            'replay': native['compression_certificate']['replay'],
            'presentation_audit': native['compression_certificate']['audit']['status'],
            'finite_carrier_assumed': False,
            'rewrites_replayed': len(native['compression_certificate']['result']['certificate']['steps']),
        },
        'upstream': {'repository': pins['rkf']['repository'], 'commit': pins['rkf']['commit'],
                     'source_count': len(checked), 'sha256': checked},
        'exact_witnesses': witness, 'independent_sheet_audit': sheet_audit,
        'preserved_evidence': previous,
        'claim_levels': {'general_results': 'written algebraic proofs with explicit assumptions',
                         'finite_results': 'exact rational tests and source proof replay',
                         'Lean_formalization_added': False, 'external_peer_review': False,
                         'physical_experiment': False, 'global_PDE_existence': False},
        'sha256': {path: sha256(root/path) for path in paths},
    }
    (here/'R8_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R8: {record['status']}; {result.testsRun} tests; 6 math mutations rejected; "
          '1 general identity and 2 KIR proof replays; 3 native tampering controls; R1-R7 and master review preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
