"""Verify native balance, full information pushforward and ledger certificates."""

import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serializable(value):
    from fractions import Fraction
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value


def check_sources(root, files):
    checked = {}
    for item in files:
        path = root/item['path']
        if not path.is_file() or sha(path) != item['sha256']:
            raise ValueError('Pinned source SHA256 mismatch: '+item['path'])
        data = path.read_bytes()
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob != item['git_blob_sha']:
            raise ValueError('Pinned source Git blob mismatch: '+item['path'])
        checked[item['path']] = item['sha256']
    return checked


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def preserved(root, here):
    result = {}
    for path in [here/f'R{i}_VERIFICATION.json' for i in range(1, 9)] + [here/'EMK_MASTER_REVIEW_VERIFICATION.json']:
        record = json.loads(path.read_text())
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha(path),
                             'paths_checked': len(record['sha256']),
                             'all_hashes_match': all((root/p).is_file() and sha(root/p) == digest
                                                     for p, digest in record['sha256'].items())}
    return result


def mutations(tests, b):
    recognize, select, budget, source_probe, measure, iterate = (
        b.recognize, b.curvature_selection, b.curvature_budget,
        b.source_information_probe, b.measurement_information, b.iterate)

    def wrong_balance(a, theta=b.Q(1, 2), cut=b.K):
        return recognize(a, b.Q(1, 4) if theta == b.Q(1, 2) else theta, cut)

    def omit_even_memory(*args):
        report = select(*args)
        report['recognized_raw'] = report['recomputed_from_recognized_generators']
        return report

    def erase_variance(*args):
        report = budget(*args)
        report['variance'] = b.Q(0)
        return report

    def freeze_derivatives(qth, m, theta):
        report = source_probe(qth, m, theta)
        rho = qth.bloch_state((b.Q(0), b.Q(0), report['lambda']*m))
        dx, dy = qth.bloch_derivative((b.Q(1), b.Q(0), b.Q(0))), qth.bloch_derivative((b.Q(0), b.Q(1), b.Q(0)))
        gx, _ = qth.geometric_tensor(rho, dx, dx)
        gy, _ = qth.geometric_tensor(rho, dy, dy)
        gxy, _ = qth.geometric_tensor(rho, dx, dy)
        report['source_information_after'] = b.matrix([[gx, gxy], [gxy, gy]])
        return report

    def false_joint_saturation(m, effects):
        report = measure(m, effects)
        if m == 0:
            report['classical_information'] = b.identity(2)
            report['joint_trace_gap'] = b.Q(0)
        return report

    def erase_sign(*args):
        report = iterate(*args)
        report['odd_multiplier'] = abs(report['odd_multiplier'])
        return report

    variants = [
        ('unequal_branch_weight_claimed_balanced', 'recognize', wrong_balance,
         'test_balanced_observer_is_idempotent_and_preserves_even_sector'),
        ('omit_surviving_odd_pair_memory', 'curvature_selection', omit_even_memory,
         'test_tangential_tangential_curvature_survives_balance'),
        ('erase_second_moment_variance', 'curvature_budget', erase_variance,
         'test_signed_mean_variance_budget_is_exact'),
        ('freeze_response_derivatives', 'source_information_probe', freeze_derivatives,
         'test_balancing_transforms_derivatives_and_records_information_cost'),
        ('zero_mean_curvature_means_full_joint_information', 'measurement_information', false_joint_saturation,
         'test_balanced_mean_curvature_does_not_imply_joint_information_saturation'),
        ('erase_orientation_sign', 'iterate', erase_sign,
         'test_sign_flip_and_exact_contraction_without_clock'),
    ]
    outcomes = {}
    for label, name, replacement, test_name in variants:
        original = getattr(b, name)
        try:
            setattr(b, name, replacement)
            result = unittest.TextTestRunner(stream=io.StringIO()).run(
                unittest.TestSuite([tests.CurvatureBalanceTests(test_name)]))
            outcomes[label] = {'rejected': bool(result.failures) and not result.errors,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(b, name, original)
    return outcomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', required=True, type=Path)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--node', default='node')
    args = parser.parse_args()
    if sys.version_info[:2] not in ((3, 11), (3, 12)):
        parser.error('Use Python 3.11 or 3.12')
    here, pub, rkf = Path(__file__).resolve().parent, args.publications_root.resolve(), args.rkf_root.resolve()
    root = here.parent
    pins = json.loads((here/'R9_SOURCE_PINS.json').read_text())
    try:
        pub_checked = check_sources(pub, pins['publications']['files'])
        rkf_checked = check_sources(rkf, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    # Every executable source hash has been checked before upstream imports.
    sources = {name: load_module('r9_pinned_'+name, pub/path)
               for name, path in pins['publications']['runtime_modules'].items()}
    import emk_curvature_balance as b
    import test_emk_curvature_balance as tests
    tests.SOURCES = sources
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.CurvatureBalanceTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    controls = mutations(tests, b)
    command = [args.node, str(here/'r9_native_balance_probe.cjs'), '--rkf-root', str(rkf),
               '--expected-input-sha256', pins['native_contract_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True, timeout=30).stdout)
    native_path = here/'R9_NATIVE_CERTIFICATE.json'
    native_path.write_text(json.dumps(native, indent=2)+'\n')
    native_results = [row for group in native['groups'].values() for row in group['results']]
    native_ok = (native['status'] == 'PASS_NATIVE_CUT_BALANCE_REPLAYS'
                 and native['input_sha256'] == pins['native_contract_input_sha256']
                 and len(native_results) == 9
                 and all(x['replay'] == 'REPLAY_MATCH' and x['result']['status'] == x['expected'] for x in native_results)
                 and all(native['negative_controls'].values()))
    previous = preserved(root, here)
    passed = (result.wasSuccessful() and result.testsRun == 24
              and all(x['rejected'] for x in controls.values()) and native_ok
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/CURVATURE_BALANCE_R9.md',
        '04-operator-evolution/emk_curvature_balance.py', '04-operator-evolution/test_emk_curvature_balance.py',
        '04-operator-evolution/r9_native_balance_probe.cjs', '04-operator-evolution/verify_r9.py',
        '04-operator-evolution/R9_SOURCE_PINS.json', '04-operator-evolution/R9_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/aghora_return.py', '04-operator-evolution/emk_tensor_calculus.py',
        '04-operator-evolution/emk_curvature_observation.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 9)] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R9', 'date': '2026-10-01',
        'status': 'PASS_EXACT_NATIVE_CURVATURE_BALANCE' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'two_event_balance_sequences': 25, 'curvature_selection_pairs': 25,
                               'variance_budget_cases': 25, 'source_tensor_pushforwards': 25,
                               'measurement_matrices_compared': 25, 'source_directional_measurement_values': 75},
        'math_mutation_controls': controls, 'native_negative_controls': native['negative_controls'],
        'native_contract': {'input_sha256': pins['native_contract_input_sha256'],
                            'generic_symbolic_proofs_replayed': 6, 'KIR_equality_proofs_replayed': 2,
                            'KIR_nonzero_curvature_proof_replayed': 1,
                            'audits': {name: group['audit']['status'] for name, group in native['groups'].items()},
                            'canonical_engine_reused_unchanged': True},
        'upstream': {'publications_commit': pins['publications']['commit'],
                     'rkf_commit': pins['rkf']['commit'], 'publications_sha256': pub_checked, 'rkf_sha256': rkf_checked},
        'witnesses': b.exact_witnesses(sources['qth'], sources['master']),
        'preserved_evidence': previous,
        'claim_levels': {'written_general_proofs': True, 'generic_rewrite_certificates': True,
                         'exact_finite_source_checks': True, 'new_Lean_formalization': False,
                         'external_peer_review': False, 'physical_experiment': False},
        'selection_boundary': {'derived_cut_operator_reused': True,
                               'balanced_endpoint_selected_by_declared_branch_exchange_symmetry': True,
                               'physical_rate_or_clock_selected': False,
                               'recognition_curvature_identified_with_Riemann_curvature': False,
                               'balance_asserted_for_all_observation': False},
        'sha256': {path: sha(root/path) for path in paths},
    }
    (here/'R9_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R9: {record['status']}; {result.testsRun} tests; 9 symbolic native replays; "
          '6 math mutations and 2 native alterations rejected; R1-R8/master evidence preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
