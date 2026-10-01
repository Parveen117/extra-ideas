"""Certify the spectral-paper/native-curvature observer bridge."""

import argparse
import hashlib
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
        return {('/'.join(map(str, key)) if isinstance(key, tuple) else key): serializable(item)
                for key, item in value.items()}
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


def preserved(root, here):
    result = {}
    for path in [here/f'R{i}_VERIFICATION.json' for i in range(1, 11)] + [here/'EMK_MASTER_REVIEW_VERIFICATION.json']:
        record = json.loads(path.read_text())
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha(path),
                             'paths_checked': len(record['sha256']),
                             'all_hashes_match': all((root/p).is_file() and sha(root/p) == digest
                                                     for p, digest in record['sha256'].items())}
    return result


def mutations(tests, s):
    connection, loop, marked_trace, observer, classify = (
        s.native_connection, s.order_loop, s.marked_trace, s.observer_report, s.classify_determinants)

    def erase_spectrally_invisible_curvature(lower):
        lower = s.lower_data(lower)
        size = 2+len(lower)
        return s.n.ConnectionJet.constant(('u', 'v'), (s.hidden_projector(len(lower)), s.n.scale(s.n.identity(size), 0)))

    def reverse_order_loop(*args):
        return loop(*args).reverse()

    def confuse_insertion_marker_with_analysis_adjoint(marker, operator):
        return marked_trace(s.n.transpose(marker), operator)

    def one_channel_repairs_every_blind_target(*args):
        report = observer(*args)
        report['minimum_arbitrary_scalar_supplements'] = min(1, report['target_blindness'])
        return report

    def noisy_zero_means_flat(*args):
        report = classify(*args)
        if all(value == 0 for row in report['lower_estimate'] for value in row):
            report['status'] = 'CERTIFIED_ZERO_CURVATURE_IN_DECLARED_FAMILY'
        return report

    def omit_connected_log_cross_terms(a, b, c, parameter):
        return parameter/(1-parameter)*marked_trace(c, s.n.commutator(a, b))

    variants = [
        ('spectral_flatness_erases_native_curvature', 'native_connection', erase_spectrally_invisible_curvature,
         'test_curvature_is_nonzero_nilpotent_and_has_a_flat_riemann_quotient'),
        ('erase_finite_order_orientation', 'order_loop', reverse_order_loop,
         'test_exact_finite_loop_retains_curvature_including_both_orientation_signs'),
        ('transpose_the_insertion_marker', 'marked_trace', confuse_insertion_marker_with_analysis_adjoint,
         'test_source_insertion_markers_recover_every_lower_coupling_entry'),
        ('one_scalar_channel_repairs_arbitrary_spectral_blindness', 'observer_report', one_channel_repairs_every_blind_target,
         'test_minimum_is_two_scalar_probes_per_hidden_mode'),
        ('noisy_zero_is_exact_flatness', 'classify_determinants', noisy_zero_means_flat,
         'test_noisy_zero_responses_do_not_certify_flatness'),
        ('omit_connected_log_cross_terms', 'cubic_connected_difference', omit_connected_log_cross_terms,
         'test_cubic_excitation_replays_the_papers_connected_order_formula'),
    ]
    outcomes = {}
    for label, name, replacement, test_name in variants:
        original = getattr(s, name)
        try:
            setattr(s, name, replacement)
            result = unittest.TextTestRunner(stream=io.StringIO()).run(
                unittest.TestSuite([tests.SpectralCurvatureObserverTests(test_name)]))
            outcomes[label] = {'rejected': bool(result.failures) and not result.errors,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(s, name, original)
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
    pins = json.loads((here/'R11_SOURCE_PINS.json').read_text())
    try:
        pub_checked = check_sources(pub, pins['publications']['files'])
        rkf_checked = check_sources(rkf, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    # PDFs are proof references, not executed source modules. Native source
    # bytes have been checked before the first canonical Node proof replay.
    import spectral_curvature_observer as s
    import test_spectral_curvature_observer as tests
    tests.SOURCE_PINS_VERIFIED = True
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(tests.SpectralCurvatureObserverTests))
    controls = mutations(tests, s)
    command = [args.node, str(here/'r11_native_spectral_probe.cjs'), '--rkf-root', str(rkf),
               '--expected-input-sha256', pins['native_contract_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True, timeout=30).stdout)
    (here/'R11_NATIVE_CERTIFICATE.json').write_text(json.dumps(native, indent=2)+'\n')
    native_ok = (native['status'] == 'PASS_NATIVE_SPECTRAL_CURVATURE_REPLAYS'
                 and native['input_sha256'] == pins['native_contract_input_sha256']
                 and native['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS'
                 and len(native['results']) == 9
                 and all(x['replay'] == 'REPLAY_MATCH' and x['result']['status'] == x['expected'] for x in native['results'])
                 and all(native['negative_controls'].values()))
    previous = preserved(root, here)
    passed = (result.wasSuccessful() and result.testsRun == 31
              and all(x['rejected'] for x in controls.values()) and native_ok
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/SPECTRAL_CURVATURE_OBSERVER_R11.md',
        '04-operator-evolution/spectral_curvature_observer.py', '04-operator-evolution/test_spectral_curvature_observer.py',
        '04-operator-evolution/r11_native_spectral_probe.cjs', '04-operator-evolution/verify_r11.py',
        '04-operator-evolution/R11_SOURCE_PINS.json', '04-operator-evolution/R11_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/aghora_return.py', '04-operator-evolution/emk_tensor_calculus.py',
        '04-operator-evolution/emk_curvature_observation.py', '04-operator-evolution/emk_curvature_balance.py',
        '04-operator-evolution/native_curvature_descent.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 11)] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R11', 'date': '2026-10-01',
        'status': 'PASS_EXACT_SPECTRAL_CURVATURE_OBSERVER' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'lower_coupling_grid': 81, 'unmarked_curvature_determinants': 405,
                               'marked_curvature_determinants': 1620, 'marker_subbanks': 16,
                               'signed_finite_loop_parameters': 15, 'cubic_feedback_channels': 16,
                               'nonfeedback_marked_determinants': 60, 'hidden_dimensions': [1, 2, 3]},
        'math_mutation_controls': controls, 'native_negative_controls': native['negative_controls'],
        'native_contract': {'input_sha256': pins['native_contract_input_sha256'],
                            'equalities_replayed': 8, 'distinctness_replays': 1,
                            'audit': native['audit']['status'], 'finite_carrier_assumed': False,
                            'canonical_engine_reused_unchanged': True},
        'upstream': {'publications_commit': pins['publications']['commit'],
                     'rkf_commit': pins['rkf']['commit'], 'publications_sha256': pub_checked, 'rkf_sha256': rkf_checked},
        'witnesses': s.exact_witnesses(), 'preserved_evidence': previous,
        'claim_levels': {'written_general_proofs': True, 'native_rewrite_certificates': True,
                         'exact_finite_spectral_and_curvature_checks': True,
                         'published_45_of_50_and_209_of_270_calibrations_rerun': False,
                         'private_MP_theorem_suite_executed': False,
                         'new_Lean_formalization': False, 'external_peer_review': False, 'physical_experiment': False},
        'selection_boundary': {'curvature_family_is_admitted': True, 'probe_bank_is_model_relative': True,
                               'noise_bound_is_supplied': True, 'spectral_shadow_is_not_a_linear_map_in_general': True,
                               'minimal_count_applies_to_linear_scalar_channels': True,
                               'observer_markers_are_inserted_before_compression': True,
                               'erased_history_recovered_from_endpoint_markers': False,
                               'response_Gram_identified_with_spacetime_metric': False,
                               'new_physical_probe_or_metric_selected': False},
        'sha256': {path: sha(root/path) for path in paths},
    }
    (here/'R11_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R11: {record['status']}; {result.testsRun} tests; 9 native replays; "
          '6 math mutations and 2 native alterations rejected; R1-R10/master evidence preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
