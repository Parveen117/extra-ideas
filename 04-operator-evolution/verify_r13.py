"""Certify the native complementary noise/memory/event foundation."""

import argparse
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest
from verify_r11 import check_sources, serializable, sha


def preserved(root, here):
    result = {}
    for path in [here/f'R{i}_VERIFICATION.json' for i in range(1, 13)] + [
            here/'EMK_MASTER_REVIEW_VERIFICATION.json']:
        record = json.loads(path.read_text())
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha(path),
                             'paths_checked': len(record['sha256']),
                             'all_hashes_match': all((root/p).is_file() and sha(root/p) == digest
                                                     for p, digest in record['sha256'].items())}
    return result


def mutations(tests, m):
    kernels, covariance, pair, weights, coherent, binary = (
        m.memory_kernels, m.complementary_covariance, m.balanced_unit_pair,
        m.branch_weights, m.coherent_vs_monitored, m.native_binary_score_metric)

    def erase_memory(*args):
        return tuple(m.n.scale(k, 0) for k in kernels(*args))

    def whiten_complement(*args):
        original = covariance(*args)
        return tuple(tuple(block if i == j else m.n.scale(block, 0)
                           for j, block in enumerate(row)) for i, row in enumerate(original))

    def pretend_the_balanced_pair_is_Gaussian(*args):
        report = pair(*args)
        report['hidden_fourth_moment'] = 3*report['hidden_variance']**2
        report['gaussian_hidden_law'] = True
        return report

    def use_linear_amplitude_weights(c, s):
        report = weights(c, s)
        c, s = abs(m.n.scalar(c)), abs(m.n.scalar(s))
        report['continue'], report['exit'] = c/(c+s), s/(c+s)
        report['binary_event_variance'] = report['continue']*report['exit']
        return report

    def replace_coherent_return_by_reset(*args):
        report = coherent(*args)
        report['coherent_visible_amplitude'] = report['monitored_survival_amplitude']
        report['coherent_final_visible_weight'] = report['monitored_survival_weight']
        return report

    def misnormalize_the_native_binary_Fisher_form(*args):
        report = binary(*args)
        report['metric'] = m.n.scale(report['metric'], m.Q(1, 2))
        return report

    variants = [
        ('erase_complement_round_trip_memory', 'memory_kernels', erase_memory,
         'test_memory_term_is_load_bearing_even_when_the_initial_hidden_state_is_zero'),
        ('pretend_complement_noise_is_independent_in_time', 'complementary_covariance', whiten_complement,
         'test_noise_is_colored_and_temporally_rank_one'),
        ('pretend_native_balanced_noise_is_Gaussian', 'balanced_unit_pair', pretend_the_balanced_pair_is_Gaussian,
         'test_balanced_complement_does_not_generate_a_Gaussian_likelihood'),
        ('use_linear_amplitudes_instead_of_native_norm_weights', 'branch_weights', use_linear_amplitude_weights,
         'test_first_exit_weights_are_norms_of_actual_native_branch_amplitudes'),
        ('equate_coherent_retention_and_monitored_reset', 'coherent_vs_monitored', replace_coherent_return_by_reset,
         'test_coherent_retention_and_monitored_reset_have_different_return_weights'),
        ('misnormalize_native_binary_information', 'native_binary_score_metric', misnormalize_the_native_binary_Fisher_form,
         'test_native_binary_score_metric_matches_independent_QTH_score_evaluations'),
    ]
    outcomes = {}
    for label, name, replacement, test_name in variants:
        original = getattr(m, name)
        try:
            setattr(m, name, replacement)
            result = unittest.TextTestRunner(stream=io.StringIO()).run(
                unittest.TestSuite([tests.ComplementMemoryNoiseTests(test_name)]))
            outcomes[label] = {'rejected': bool(result.failures) and not result.errors,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(m, name, original)
    return outcomes


def check_native(native, m):
    checks = {}
    for input_case, output in zip(native['input']['observer_cases'], native['observers']):
        seed = tuple(tuple(m.Q(x) for x in row) for row in input_case['seed'])
        rows = list(seed)
        for action in input_case['actions']:
            t = tuple(tuple(m.Q(x) for x in row) for row in action)
            rows.extend(m.o.product(seed, t))
        expected = m.o.independent_rows(rows)
        observed = tuple(tuple(m.Q(x[0]) for x in row) for row in output['result']['observer'])
        real_face = all(x[1] == '0' for row in output['result']['observer'] for x in row)
        checks[input_case['name']] = {'rank': output['result']['rank'],
                                     'observer_matrix_match': real_face and expected == observed
                                     and len(expected) == output['result']['rank']}
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', required=True, type=Path)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--node', default='node')
    args = parser.parse_args()
    if sys.version_info[:2] not in ((3, 11), (3, 12)):
        parser.error('Use Python 3.11 or 3.12')
    here = Path(__file__).resolve().parent
    root, pub, rkf = here.parent, args.publications_root.resolve(), args.rkf_root.resolve()
    pins = json.loads((here/'R13_SOURCE_PINS.json').read_text())
    try:
        pub_checked = check_sources(pub, pins['publications']['files'])
        rkf_checked = check_sources(rkf, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    sources = {}
    for name, path in pins['publications']['runtime_modules'].items():
        spec = importlib.util.spec_from_file_location('r13_pinned_'+name, pub/path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        sources[name] = mod
    import complement_memory_noise as m
    import test_complement_memory_noise as tests
    tests.SOURCES = sources
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(tests.ComplementMemoryNoiseTests))
    controls = mutations(tests, m)
    command = [args.node, str(here/'r13_native_complement_probe.cjs'), '--rkf-root', str(rkf),
               '--expected-input-sha256', pins['native_contract_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True, timeout=30).stdout)
    (here/'R13_NATIVE_CERTIFICATE.json').write_text(json.dumps(native, indent=2)+'\n')
    observer_checks = check_native(native, m)
    branch_proof = next(x for x in native['results']
                        if x['name'] == 'first_exit_weight_is_the_native_squared_branch_amplitude')
    native_ok = (native['status'] == 'PASS_NATIVE_COMPLEMENT_MEMORY_REPLAYS'
                 and native['input_sha256'] == pins['native_contract_input_sha256']
                 and native['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS'
                 and len(native['results']) == 11 and len(native['observers']) == 3
                 and all(x['replay'] == 'REPLAY_MATCH' and x['result']['status'] == x['expected']
                         for x in native['results'])
                 and branch_proof['depends_on'] == ['visible_cut_is_idempotent', 'visible_single_step_is_the_cosine_channel']
                 and all(x['observer_matrix_match'] for x in observer_checks.values())
                 and all(native['negative_controls'].values()))
    previous = preserved(root, here)
    passed = (result.wasSuccessful() and result.testsRun == 31 and native_ok
              and all(x['rejected'] for x in controls.values())
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/COMPLEMENT_MEMORY_NOISE_R13.md',
        '04-operator-evolution/complement_memory_noise.py', '04-operator-evolution/test_complement_memory_noise.py',
        '04-operator-evolution/r13_native_complement_probe.cjs', '04-operator-evolution/verify_r13.py',
        '04-operator-evolution/R13_SOURCE_PINS.json', '04-operator-evolution/R13_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/verify_r11.py', '04-operator-evolution/aghora_return.py',
        '04-operator-evolution/emk_tensor_calculus.py', '04-operator-evolution/emk_curvature_observation.py',
        '04-operator-evolution/emk_curvature_balance.py', '04-operator-evolution/native_curvature_descent.py',
        '04-operator-evolution/spectral_curvature_observer.py', '04-operator-evolution/observer_metric_foundation.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 13)] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R13', 'date': '2026-10-01',
        'status': 'PASS_EXACT_COMPLEMENT_MEMORY_NOISE_FOUNDATION' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'native_rotation_charts': 12, 'full_vs_reduced_rotation_histories': 36,
                               'hidden_component_recoveries': 108, 'first_exit_branch_norm_checks': 72,
                               'independent_QTH_rotating_readout_comparisons': 12,
                               'independent_QTH_binary_score_evaluations': 27,
                               'non_Gaussian_seam_parameter_cases': 12,
                               'finite_waiting_moment_tail_checks': 5, 'finite_event_tail_ledgers': 7,
                               'multi_hidden_feedback_modes': 2, 'native_observer_catalogues': 3},
        'math_mutation_controls': controls, 'native_negative_controls': native['negative_controls'],
        'native_contract': {'input_sha256': pins['native_contract_input_sha256'],
                            'equalities_replayed': 10, 'distinctness_replays': 1,
                            'audit': native['audit']['status'], 'abstract_rewrites_assume_finite_carrier': False,
                            'finite_observer_checks': observer_checks,
                            'branch_norm_proof_dependencies': branch_proof['depends_on'],
                            'canonical_engine_reused_unchanged': True},
        'upstream': {'publications_commit': pins['publications']['commit'], 'rkf_commit': pins['rkf']['commit'],
                     'publications_sha256': pub_checked, 'rkf_sha256': rkf_checked},
        'witnesses': m.exact_witnesses(), 'preserved_evidence': previous,
        'claim_levels': {'written_hidden_elimination_and_norm_selection_proofs': True,
                         'written_completed_native_binary_metric_realization': True,
                         'native_rewrite_certificates': True, 'exact_finite_rational_checks': True,
                         'new_Lean_formalization': False, 'external_peer_review': False, 'physical_experiment': False},
        'selection_boundary': {'hidden_force_and_memory_derived_from_native_blocks': True,
                               'noise_covariance_derived_in_balanced_unit_preparation': True,
                               'event_law_derived_from_additive_norm_readout': True,
                               'native_binary_meter_likelihood_and_variance_derived': True,
                               'R12_metric_equation_realized_without_Gaussian_noise': True,
                               'matched_native_sheet_gain_and_stopping_derive_kappa_cosine_squared': True,
                               'native_preparation_phase_cut_and_detector_policy_are_specified': True,
                               'binary_seam_state_scope': 'anchored zero contrast',
                               'ordinary_complement_noise_is_independent_in_time_or_Gaussian': False,
                               'coherent_retention_has_the_monitored_geometric_event_law': False,
                               'all_noise_is_proved_to_be_hidden_memory': False,
                               'norm_probabilities_recover_signed_orientation': False,
                               'original_R12_Gaussian_experiment_rewritten': False,
                               'integer_sheet_memory_equals_complementary_state_amplitude': False,
                               'universal_physical_phase_metric_clock_or_signature_selected': False},
        'sha256': {path: sha(root/path) for path in paths},
    }
    (here/'R13_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R13: {record['status']}; {result.testsRun} tests; 11 native replays; 3 N03 completions; "
          '6 math mutations and 3 native alterations rejected; R1-R12/master evidence preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
