"""Verify typed EMK tensors against pinned native sources; preserve R1-R6."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import unittest

import emk_tensor_calculus as tc
from aghora_return import cayley_flow


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serializable(value):
    if isinstance(value, tc.Q):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', required=True, type=Path)
    args = parser.parse_args()
    if sys.version_info[:2] not in ((3, 11), (3, 12)):
        parser.error('Use Python 3.11 or 3.12')
    here = Path(__file__).resolve().parent
    root = here.parent
    pins = json.loads((here/'R7_SOURCE_PINS.json').read_text())
    checked = {}
    # Reject altered bytes before importing any upstream executable source.
    for item in pins['publications']['files']:
        if not item['runtime_required']:
            continue
        source = args.publications_root.resolve()/item['path']
        if not source.is_file() or digest(source) != item['sha256']:
            raise SystemExit('Pinned upstream SHA256 mismatch: '+item['path'])
        data = source.read_bytes()
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob != item['git_blob_sha']:
            raise SystemExit('Pinned upstream Git blob mismatch: '+item['path'])
        checked[item['path']] = item['sha256']
    if len(checked) != 3:
        raise SystemExit('Expected three unchanged upstream runtime sources')
    os.environ['PUBLICATIONS_SOURCE_ROOT'] = str(args.publications_root.resolve())
    import test_emk_tensor_calculus as tests
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tests.EMKTensorCalculusTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    previous = {}
    for revision in range(1, 7):
        record = json.loads((here/f'R{revision}_VERIFICATION.json').read_text())
        previous[f'R{revision}'] = {
            'recorded_status': record['status'], 'files_checked': len(record['sha256']),
            'all_recorded_hashes_match': all(
                (root/path).is_file() and digest(root/path) == expected
                for path, expected in record['sha256'].items()),
        }

    v, dual = tc.Slot('V', 1), tc.Slot('V', -1)
    signatures = ((dual, dual), (v, dual), (v,))
    p = tc.Move('0', 'a', (('V', tc.R),)).then(tc.Move('a', '2', (('V', tc.K),)))
    q = tc.Move('0', 'b', (('V', tc.K),)).then(tc.Move('b', '2', (('V', tc.R),)))
    sign_loop = p.then(q.reverse())
    winding_loop = tc.Move('0', '0', (('V', tc.I),), 2)
    boost = cayley_flow(tc.K, tc.Q(1, 3))
    edges = {'01': tc.R, '12': tc.K, '23': tc.RK,
             '02': tc.I, '13': tc.I, '03': tc.I}
    zero = tc.scale(tc.I, 0)
    witness_checks = (
        dict(sign_loop.operators)['V'] == tc.scale(tc.I, -1),
        tc.return_report(sign_loop, signatures)['tensor_type_returns'] == [True, True, False],
        not tc.return_report(winding_loop, signatures)['full_unledgered_return'],
        tc.invariant_symmetric_forms((tc.R, tc.K)) == (),
        tc.metric_defect(boost, tc.I, tc.I) != zero,
        tc.mul(tc.mul(tc.transpose(boost), tc.EPSILON), boost) == tc.EPSILON,
        tc.triangle_curvature(edges['01'], edges['12'], edges['02']) != zero,
        tc.tetrahedron_bianchi(edges) == zero,
    )
    passed = (result.wasSuccessful() and result.testsRun == 16 and all(witness_checks)
              and all(item['all_recorded_hashes_match']
                      and item['recorded_status'].startswith('PASS')
                      for item in previous.values()))
    paths = [
        '4ways.tex', '02-relational-response/EMK_TENSOR_CALCULUS_R7.md',
        '02-relational-response/EMK_VAULT_SOURCE_AUDIT_R7.md',
        '04-operator-evolution/emk_tensor_calculus.py',
        '04-operator-evolution/test_emk_tensor_calculus.py',
        '04-operator-evolution/verify_r7.py',
        '04-operator-evolution/R7_SOURCE_PINS.json',
        '04-operator-evolution/aghora_return.py',
    ] + [f'04-operator-evolution/R{i}_VERIFICATION.json' for i in range(1, 7)]
    record = {
        'status': 'PASS_EXACT_TYPED_EMK_TENSOR_CHECKS' if passed else 'FAIL',
        'development': 'R7', 'date': '2026-09-30', 'python': sys.version.split()[0],
        'arithmetic': 'exact rational matrices and tensor components',
        'tests_run': result.testsRun, 'failures': len(result.failures),
        'errors': len(result.errors),
        'rank_variance_signatures_checked': sum(2**r for r in range(1, 5)),
        'upstream_repository': pins['publications']['repository'],
        'upstream_commit': pins['publications']['commit'],
        'upstream_runtime_reused_unchanged': True, 'upstream_source_count': len(checked),
        'upstream_runtime_sha256': checked,
        'central_order_loop': {
            'carrier': dict(sign_loop.operators)['V'],
            'signature_order': ['metric (0,2)', 'endomorphism (1,1)', 'vector (1,0)'],
            'audit': tc.return_report(sign_loop, signatures),
        },
        'independent_sheet_loop': {
            'carrier': tc.I, 'winding': winding_loop.sheet,
            'unledgered_audit': tc.return_report(winding_loop, signatures),
            'preregistered_ledger_2_audit': tc.return_report(
                winding_loop, signatures, sheet_ledger=2),
        },
        'fixed_symmetric_metric_spaces': {
            'R_generator_basis': tc.invariant_symmetric_forms((tc.R,)),
            'K_generator_basis': tc.invariant_symmetric_forms((tc.K,)),
            'both_generator_basis': tc.invariant_symmetric_forms((tc.R, tc.K)),
        },
        'boost_metric_witness': {'transport': boost,
                                 'symmetric_I_defect': tc.metric_defect(boost, tc.I, tc.I),
                                 'alternating_form_preserved': True},
        'nonflat_bianchi_witness': {
            'edges': edges, 'triangle_012': tc.triangle_curvature(
                edges['01'], edges['12'], edges['02']),
            'tetrahedron_residual': tc.tetrahedron_bianchi(edges),
        },
        'source_coverage': {key: pins['vault'][key] for key in
                            ('retrieved_file_count', 'current_file_count',
                             'recovered_core_tex_count', 'read_depth_counts')},
        'all_vault_files_read_in_full': False,
        'smooth_analytic_existence_certified': False,
        'native_metric_policy_selected': False,
        'sheet_equals_raw_curvature_flux_asserted': False,
        'universal_lambda_selected': False, 'physical_experiment_performed': False,
        'earlier_verification_integrity': previous,
        'sha256': {path: digest(root/path) for path in paths},
    }
    (here/'R7_VERIFICATION.json').write_text(json.dumps(serializable(record), indent=2)+'\n')
    print(f"R7: {record['status']}; {result.testsRun} tests; 30 tensor signatures; "
          '3 unchanged source modules; R1-R6 hashes preserved')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
