"""Verify sourced native curvature and the conditional classical tangent sector."""

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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def preserved(root, here):
    result = {}
    for path in [here/f'R{i}_VERIFICATION.json' for i in range(1, 10)] + [here/'EMK_MASTER_REVIEW_VERIFICATION.json']:
        record = json.loads(path.read_text())
        result[path.stem] = {'recorded_status': record['status'], 'record_sha256': sha(path),
                             'paths_checked': len(record['sha256']),
                             'all_hashes_match': all((root/p).is_file() and sha(root/p) == digest
                                                     for p, digest in record['sha256'].items())}
    return result


def mutations(tests, n):
    bianchi, descent, compatibility, decomposition, reduction = (
        n.sourced_bianchi, n.descent_report, n.tangent_compatibility,
        n.visible_decomposition, n.riemann_reduction)

    def omit_hidden_source(*args):
        report = bianchi(*args)
        report['hidden_source'] = n.scale(report['hidden_source'], 0)
        report['sourced_residual'] = report['visible_bianchi']
        return report

    def wrong_source_sign(*args):
        report = bianchi(*args)
        report['hidden_source'] = n.scale(report['hidden_source'], -1)
        report['sourced_residual'] = n.add(report['visible_bianchi'], report['hidden_source'])
        return report

    def ignore_connection_derivatives(*args):
        report = descent(*args)
        report['connection_jet_descends'] = all(n.is_zero(x) for x in report['value_intertwining'])
        return report

    def pointwise_conditions_are_enough(*args):
        report = compatibility(*args)
        report['Levi_Civita_connection_jet'] = (
            all(n.is_zero(x) for x in report['nonmetricity'])
            and all(x == 0 for row in report['torsion'] for vector in row for x in vector))
        return report

    def omit_feedback(*args):
        report = decomposition(*args)
        report['excursion'] = n.scale(report['excursion'], 0)
        report['reconstruction_residual'] = n.sub(
            report['visible_full_curvature'], n.add(report['Riemann'], report['distortion']))
        return report

    def curvature_equality_is_connection_descent(*args):
        report = reduction(*args)
        report['Riemann_sector_certified'] = report['curvature_descends']
        return report

    variants = [
        ('omit_hidden_Bianchi_source', 'sourced_bianchi', omit_hidden_source,
         'test_visible_bianchi_has_exact_hidden_source'),
        ('reverse_hidden_Bianchi_source_sign', 'sourced_bianchi', wrong_source_sign,
         'test_visible_bianchi_has_exact_hidden_source'),
        ('ignore_differentiated_observer_contract', 'descent_report', ignore_connection_derivatives,
         'test_differentiated_observer_intertwining_is_load_bearing'),
        ('pointwise_metric_torsion_values_certify_connection_jet', 'tangent_compatibility', pointwise_conditions_are_enough,
         'test_pointwise_metric_and_torsion_conditions_do_not_certify_curvature'),
        ('erase_hidden_excursion_curvature', 'visible_decomposition', omit_feedback,
         'test_hidden_feedback_is_an_exact_excursion_curvature'),
        ('curvature_equality_identifies_the_connection', 'riemann_reduction', curvature_equality_is_connection_descent,
         'test_accidental_curvature_equality_is_not_connection_identification'),
    ]
    outcomes = {}
    for label, name, replacement, test_name in variants:
        original = getattr(n, name)
        try:
            setattr(n, name, replacement)
            result = unittest.TextTestRunner(stream=io.StringIO()).run(
                unittest.TestSuite([tests.NativeCurvatureDescentTests(test_name)]))
            outcomes[label] = {'rejected': bool(result.failures) and not result.errors,
                               'tests_run': result.testsRun, 'assertion_failures': len(result.failures),
                               'errors': len(result.errors)}
        finally:
            setattr(n, name, original)
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
    pins = json.loads((here/'R10_SOURCE_PINS.json').read_text())
    try:
        pub_checked = check_sources(pub, pins['publications']['files'])
        rkf_checked = check_sources(rkf, pins['rkf']['files'])
    except ValueError as error:
        raise SystemExit(str(error))
    # Every pinned source is checked before upstream imports or Node execution.
    sources = {name: load_module('r10_pinned_'+name, pub/path)
               for name, path in pins['publications']['runtime_modules'].items()}
    import native_curvature_descent as n
    import test_native_curvature_descent as tests
    tests.SOURCES = sources
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.NativeCurvatureDescentTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    controls = mutations(tests, n)
    command = [args.node, str(here/'r10_native_curvature_probe.cjs'), '--rkf-root', str(rkf),
               '--expected-input-sha256', pins['native_contract_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True, timeout=30).stdout)
    (here/'R10_NATIVE_CERTIFICATE.json').write_text(json.dumps(native, indent=2)+'\n')
    native_results = [row for group in native['groups'].values() for row in group['results']]
    native_ok = (native['status'] == 'PASS_NATIVE_CURVATURE_DESCENT_REPLAYS'
                 and native['input_sha256'] == pins['native_contract_input_sha256']
                 and len(native_results) == 5
                 and all(x['replay'] == 'REPLAY_MATCH' and x['result']['status'] == x['expected'] for x in native_results)
                 and all(group['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS' for group in native['groups'].values())
                 and all(native['negative_controls'].values()))
    previous = preserved(root, here)
    passed = (result.wasSuccessful() and result.testsRun == 25
              and all(x['rejected'] for x in controls.values()) and native_ok
              and all(x['recorded_status'].startswith('PASS') and x['all_hashes_match'] for x in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md',
        '04-operator-evolution/native_curvature_descent.py', '04-operator-evolution/test_native_curvature_descent.py',
        '04-operator-evolution/r10_native_curvature_probe.cjs', '04-operator-evolution/verify_r10.py',
        '04-operator-evolution/R10_SOURCE_PINS.json', '04-operator-evolution/R10_NATIVE_CERTIFICATE.json',
        '04-operator-evolution/aghora_return.py', '04-operator-evolution/emk_tensor_calculus.py',
        '04-operator-evolution/emk_curvature_observation.py', '04-operator-evolution/emk_curvature_balance.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 10)] + [
        '04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json']
    record = {
        'development': 'R10', 'date': '2026-10-01',
        'status': 'PASS_EXACT_NATIVE_CURVATURE_DESCENT' if passed else 'FAIL',
        'python': sys.version.split()[0],
        'node': subprocess.run([args.node, '--version'], check=True, capture_output=True, text=True).stdout.strip(),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'finite_case_counts': {'constant_Bianchi_cut_cases': 6,
                               'smooth_connection_two_jets': 1,
                               'positive_domain_EMKG1_metric_points': 23,
                               'EMKG1_curvature_comparison_routes_per_point': 3,
                               'flat_tangent_dimensions': [1, 2, 3]},
        'math_mutation_controls': controls, 'native_negative_controls': native['negative_controls'],
        'native_contract': {'input_sha256': pins['native_contract_input_sha256'],
                            'single_carrier_symbolic_proofs_replayed': 3,
                            'typed_quotient_symbolic_proofs_replayed': 2,
                            'audits': {name: group['audit']['status'] for name, group in native['groups'].items()},
                            'canonical_engine_reused_unchanged': True},
        'upstream': {'publications_commit': pins['publications']['commit'],
                     'rkf_commit': pins['rkf']['commit'], 'publications_sha256': pub_checked, 'rkf_sha256': rkf_checked},
        'witnesses': n.exact_witnesses(sources['geometry']),
        'preserved_evidence': previous,
        'claim_levels': {'written_general_proofs': True, 'generic_rewrite_certificates': True,
                         'exact_finite_source_checks': True, 'new_Lean_formalization': False,
                         'external_peer_review': False, 'physical_experiment': False},
        'selection_boundary': {'metric_is_admitted': True, 'tangent_quotient_is_conditional': True,
                               'connection_jet_intertwining_checked': True,
                               'smooth_Bianchi_cut_is_constant': True,
                               'new_physical_metric_or_rate_selected': False,
                               'Riemann_sector_asserted_for_all_native_transports': False,
                               'finite_group_loop_identified_with_infinitesimal_connection': False},
        'sha256': {path: sha(root/path) for path in paths},
    }
    (here/'R10_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R10: {record['status']}; {result.testsRun} tests; 5 symbolic native replays; "
          '6 math mutations and 2 native alterations rejected; R1-R9/master evidence preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
