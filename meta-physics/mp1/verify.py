"""MP-1 replay: Python 3.12, stdlib only, no private credentials."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import unittest
import axis_memory as am

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def json_value(value):
    if isinstance(value, am.Q):
        return str(value)
    if isinstance(value, dict):
        return {k: json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return value


def source_checks():
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for row in pins['consumed_unchanged_files']:
        p = (ROOT/row['path']).resolve()
        require(p.is_relative_to(ROOT), 'Source path escapes repository')
        b = p.read_bytes()
        require(sha(b) == row['sha256'], 'Source changed: '+row['path'])
        gitsha = hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        require(gitsha == row['git_blob_sha1'], 'Git source digest changed')
    spec = importlib.util.spec_from_file_location(
        'old_register_checker', ROOT/'04-operator-evolution/verify_repository.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    stages = old.verify_register(json.loads((ROOT/'CURRENT_MANIFEST.json').read_text()))
    require(stages == 46, 'MP-1 must not rewrite the historical R1-R46 register')
    return {'unchanged_dependencies': len(pins['consumed_unchanged_files']),
            'frozen_registered_stages': stages}


def inventory_checks():
    inventory = json.loads((HERE.parent/'SOURCE_INVENTORY.json').read_text())
    files = inventory['files']
    require([x['id'] for x in files] == [f'{i:03d}' for i in range(1, 99)], 'Incomplete PDF inventory')
    require(sum(x['pages'] for x in files) == inventory['page_count'] == 299, 'Page count mismatch')
    require(inventory['file_count'] == len(files) == 98, 'File count mismatch')
    require(len(set(x['extracted_text_sha256'] for x in files)) ==
            inventory['distinct_extracted_text_count'] == 95, 'Text duplicate count mismatch')
    require(inventory['raw_pdfs_published'] is False, 'Raw PDF publication is outside this packet')
    require(not list(HERE.parent.rglob('*.pdf')), 'Raw PDFs do not belong in this packet')
    for row in files:
        require(row['audit_gate'] and row['status'] in
                ('pending', 'scoped MP-1 development; remainder pending'), 'Unscoped source promotion')
        for field, n in [('git_blob_sha1', 40), ('pdf_sha256', 64), ('extracted_text_sha256', 64)]:
            require(len(row[field]) == n and set(row[field]) <= set('0123456789abcdef'), 'Malformed hash')
        if row['same_extracted_text_as']:
            previous = files[int(row['same_extracted_text_as']) - 1]
            require(previous['extracted_text_sha256'] == row['extracted_text_sha256'], 'False text duplicate')
            require(previous['pdf_sha256'] != row['pdf_sha256'], 'Binary/text distinction changed')
    selected = [x['id'] for x in files if x['status'] != 'pending']
    require(selected == ['017', '078'], 'Changed source certification scope')
    return {'files_read': 98, 'pages_read': 299, 'distinct_extracted_texts': 95,
            'scoped_development_sources': selected, 'complete_pdfs_certified': 0,
            'private_source_hashes_rechecked_by_public_ci': False}


def build_report():
    require(sys.version_info[:2] == (3, 12), 'This packet uses Python 3.12 only')
    integrity = source_checks()
    inventory = inventory_checks()
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(HERE), pattern='test_*.py')
    results = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    require(results.wasSuccessful(), stream.getvalue())
    require(results.testsRun == 23 and not results.skipped, 'Unexpected test coverage')
    examples = {str(n): json_value(am.polygon_report(am.diamond(n))) for n in [-2, -1, 0, 1, 2, 3]}
    proof = (HERE/'THEOREM.md').read_text()
    require(all(f'## T{i} ' in proof for i in range(1, 8)), 'Missing written statement')
    inputs = ['meta-physics/README.md', 'meta-physics/SOURCE_AUDIT.md',
              'meta-physics/SOURCE_INVENTORY.json', 'meta-physics/mp1/THEOREM.md',
              'meta-physics/mp1/axis_memory.py', 'meta-physics/mp1/test_axis_memory.py',
              'meta-physics/mp1/SOURCE_PINS.json', 'meta-physics/mp1/verify.py',
              '.github/workflows/meta-physics.yml']
    report = {
        'schema': 'MP1-VERIFICATION-1',
        'status': 'PASS_SCOPED_WRITTEN_PROOFS_AND_EXACT_CHECKS',
        'written_results': 7, 'python': '3.12', 'tests_passed': results.testsRun,
        'finite_vertex_perturbations_checked': 6561,
        'universal_protection_basis': 'T4 written proof, not extrapolation from finite tests',
        'formal_proof_assistant_verified': False,
        'physical_law_or_empirical_effect_certified': False,
        'integrity': integrity, 'inventory': inventory, 'exact_loop_examples': examples,
        'noise_example': json_value(am.polygon_report(am.diamond(1), vertex_error=am.Q(1, 8),
                                                     interpolation_error=am.Q(1, 8))),
        'boundary_control': json_value(am.polygon_report(am.diamond(1), vertex_error=am.Q(1, 2))),
        'reference_probe': {str(sign): str(am.reference_intensity((sign, 0))) for sign in [1, -1]},
        'unresolved_gates': ['physical observer/source selection', 'between-sample error bound',
                             'entropy cost', 'scalar event threshold', 'probabilistic noise law'],
        'input_sha256': {p: sha((ROOT/p).read_bytes()) for p in inputs},
    }
    report['canonical_report_sha256'] = sha(canonical(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', action='store_true', help='write a reviewed new certificate')
    modes.add_argument('--check', action='store_true', help='replay and compare without changing files')
    args = parser.parse_args()
    report = build_report()
    path = HERE/'VERIFICATION.json'
    data = json.dumps(report, indent=2, ensure_ascii=False)+'\n'
    if args.write:
        path.write_text(data)
    else:
        require(path.read_text() == data, 'Certificate differs from replay; review the change')
    print(json.dumps({'status': report['status'], 'tests': report['tests_passed'],
                      'certificate_sha256': sha(data.encode()), 'private_pdfs_needed': False}))


if __name__ == '__main__':
    main()
