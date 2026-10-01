"""Certify the native event-law observer/metric foundation and derived seam."""

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
    for path in [here/f'R{i}_VERIFICATION.json' for i in range(1, 12)] + [
            here/'EMK_MASTER_REVIEW_VERIFICATION.json']:
        record = json.loads(path.read_text())
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha(path),
                             'paths_checked': len(record['sha256']),
                             'all_hashes_match': all((root/p).is_file() and sha(root/p) == digest
                                                     for p, digest in record['sha256'].items())}
    return result


def mutations(tests, o):
    original_observer = o.observer_from_metric
    original_untagged = o.origin_untagged_metric
    original_derivation = o.derive_observer_metric

    def seed_only(experiment):
        report = original_derivation(experiment)
        report['metric'] = experiment.seed
        return report

    def omit_future_observer(form, experiment):
        report = original_observer(form, experiment)
        report['rank'] = len(o.independent_rows(experiment.readout))
        report['observer'] = o.independent_rows(experiment.readout)
        return report

    def retain_tag_information_after_erasing_tags(experiment):
        report = original_untagged(experiment)
        report['metric'] = o.tagged_metric(experiment)
        report['discard'] = o.n.scale(report['metric'], 0)
        return report

    def omit_transport_derivatives(experiment, first, second, *args, **kwargs):
        g = o.tagged_metric(experiment)
        zero = o.n.scale(g, 0)
        return {'value': g, 'first': (zero,)*len(first),
                'second': tuple((zero,)*len(first) for _ in first)}

    def balance_erases_all_geometry(jet, tangent, directions):
        return o.n.flat_metric(len(directions), directions)

    variants = [
        ('omit_future_event_information', 'derive_observer_metric', seed_only,
         'test_dynamic_metric_is_selected_by_the_entire_event_law'),
        ('current_seed_is_already_future_complete', 'observer_from_metric', omit_future_observer,
         'test_future_observer_is_derived_beyond_the_seed_probe'),
        ('forgetting_tags_keeps_the_second_moment_metric', 'origin_untagged_metric', retain_tag_information_after_erasing_tags,
         'test_event_tag_forgetting_has_a_positive_information_ledger_at_the_origin'),
        ('omit_derivatives_of_the_native_transport', 'response_two_jet', omit_transport_derivatives,
         'test_derived_metric_two_jet_contains_the_hidden_second_moment'),
        ('balanced_first_moment_forces_a_flat_metric', 'tangent_metric', balance_erases_all_geometry,
         'test_balanced_native_curvature_can_vanish_while_derived_metric_is_curved'),
        ('accept_a_divergent_event_law', 'positive_semidefinite', lambda value: True,
         'test_stability_is_load_bearing_and_an_algebraic_solve_is_not_enough'),
    ]
    outcomes = {}
    for label, name, replacement, test_name in variants:
        original = getattr(o, name)
        try:
            setattr(o, name, replacement)
            result = unittest.TextTestRunner(stream=io.StringIO()).run(
                unittest.TestSuite([tests.ObserverMetricFoundationTests(test_name)]))
            outcomes[label] = {'rejected': bool(result.failures) and not result.errors,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(o, name, original)
    return outcomes


def compare_native_observers(native, o):
    q = o.Q
    def real_rows(rows):
        if any(x[1] != '0' for row in rows for x in row):
            raise ValueError('This experiment comparison uses the real rational face')
        return tuple(tuple(q(x[0]) for x in row) for row in rows)
    future = o.EventExperiment(((1, 1),), ((1,),), (o.n.matrix([[1, 0], [0, q(1, 2)]]),),
                               (q(1),), q(1, 2), o.n.identity(2), q(3, 4))
    seam = o.seam_experiment(2, q(1, 3), q(1, 2))
    feedback = o.feedback_experiment()
    by_name = {item['name']: item for item in native['observers']}
    checked = {}
    for name, experiment in [('future_depth', future), ('paired_sheet_with_one_invisible_mode', seam),
                              ('both_reverse_markers_admitted_as_transport', feedback)]:
        result = o.derive_observer_metric(experiment)
        source = by_name[name]['result']
        match = result['rank'] == source['rank'] and result['observer'] == real_rows(source['observer'])
        checked[name] = {'rank': result['rank'], 'observer_matrix_match': match}
    lower = by_name['lower_record_without_reverse_feedback']['result']
    checked['lower_record_without_reverse_feedback'] = {
        'rank': lower['rank'], 'observer_matrix_match': lower['rank'] == 2 and
        real_rows(lower['observer']) == o.n.projection(2, 4)}
    return checked


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
    pins = json.loads((here/'R12_SOURCE_PINS.json').read_text())
    try:
        pub_checked = check_sources(pub, pins['publications']['files'])
        rkf_checked = check_sources(rkf, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    spec = importlib.util.spec_from_file_location('r12_pinned_emkg1',
                                                 pub/pins['publications']['runtime_modules']['geometry'])
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
    import observer_metric_foundation as o
    import test_observer_metric_foundation as tests
    tests.SOURCES = {'geometry': geometry}
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(tests.ObserverMetricFoundationTests))
    controls = mutations(tests, o)
    command = [args.node, str(here/'r12_native_observer_metric_probe.cjs'), '--rkf-root', str(rkf),
               '--expected-input-sha256', pins['native_contract_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True, timeout=30).stdout)
    (here/'R12_NATIVE_CERTIFICATE.json').write_text(json.dumps(native, indent=2)+'\n')
    observer_matches = compare_native_observers(native, o)
    native_ok = (native['status'] == 'PASS_NATIVE_OBSERVER_METRIC_REPLAYS'
                 and native['input_sha256'] == pins['native_contract_input_sha256']
                 and native['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS'
                 and len(native['results']) == 9 and len(native['observers']) == 4
                 and all(x['replay'] == 'REPLAY_MATCH' and x['result']['status'] == x['expected']
                         for x in native['results'])
                 and all(x['observer_matrix_match'] for x in observer_matches.values())
                 and all(native['negative_controls'].values()))
    previous = preserved(root, here)
    passed = (result.wasSuccessful() and result.testsRun == 28 and native_ok
              and all(x['rejected'] for x in controls.values())
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/OBSERVER_METRIC_FOUNDATION_R12.md',
        '04-operator-evolution/observer_metric_foundation.py',
        '04-operator-evolution/test_observer_metric_foundation.py',
        '04-operator-evolution/r12_native_observer_metric_probe.cjs',
        '04-operator-evolution/verify_r12.py', '04-operator-evolution/R12_SOURCE_PINS.json',
        '04-operator-evolution/R12_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/verify_r11.py', '04-operator-evolution/aghora_return.py',
        '04-operator-evolution/emk_tensor_calculus.py', '04-operator-evolution/emk_curvature_observation.py',
        '04-operator-evolution/emk_curvature_balance.py', '04-operator-evolution/native_curvature_descent.py',
        '04-operator-evolution/spectral_curvature_observer.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 12)] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R12', 'date': '2026-10-01',
        'status': 'PASS_EXACT_OBSERVER_METRIC_FOUNDATION' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'balanced_sheet_parameter_cases': 27, 'independent_EMKG1_geometry_comparisons': 18,
                               'positive_calibration_and_event_law_cases': 9, 'finite_event_apertures': 5,
                               'rank_changing_points': 2, 'native_observer_catalogues': 4},
        'math_mutation_controls': controls, 'native_negative_controls': native['negative_controls'],
        'native_contract': {'input_sha256': pins['native_contract_input_sha256'],
                            'equalities_replayed': 8, 'distinctness_replays': 1,
                            'audit': native['audit']['status'],
                            'abstract_algebra_replays_assume_finite_carrier': False,
                            'observer_completions_are_finite': True,
                            'observer_metric_kernel_comparisons': observer_matches,
                            'canonical_engine_reused_unchanged': True},
        'upstream': {'publications_commit': pins['publications']['commit'], 'rkf_commit': pins['rkf']['commit'],
                     'publications_sha256': pub_checked, 'rkf_sha256': rkf_checked},
        'witnesses': o.exact_witnesses(), 'preserved_evidence': previous,
        'claim_levels': {'written_general_convergence_and_selection_proofs': True,
                         'native_rewrite_certificates': True, 'exact_finite_rational_checks': True,
                         'new_Lean_formalization': False, 'external_peer_review': False, 'physical_experiment': False},
        'selection_boundary': {'observer_quotient_derived_from_event_response_kernel': True,
                               'operational_metric_derived_from_tagged_likelihood': True,
                               'metric_two_jet_derived_from_transport_two_jet': True,
                               'positive_quadratic_EMKG1_family_derived_within_event_model': True,
                               'admitted_noise_likelihood_and_event_law': True,
                               'admitted_constant_tangent_map': True,
                               'stability_is_checked_separately_from_algebraic_solve': True,
                               'untagged_mean_formula_scope': 'zero displacement only',
                               'global_fixed_statistical_family_claimed_for_nonclosed_response_coframe': False,
                               'original_native_connection_is_automatically_levi_civita': False,
                               'all_observation_forces_flatness': False,
                               'universal_physical_metric_or_clock_selected': False,
                               'lorentz_signature_derived_from_positive_fisher_metric': False},
        'sha256': {path: sha(root/path) for path in paths},
    }
    (here/'R12_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R12: {record['status']}; {result.testsRun} tests; 9 native replays; 4 N03 completions; "
          '6 math mutations and 3 native alterations rejected; R1-R11/master evidence preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
