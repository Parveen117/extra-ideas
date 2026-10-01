#!/usr/bin/env python3
"""Replay R23's native future-response quotient and observer completion and its immutable source chain."""
import argparse
import copy
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from verify_r18 import canonical, check_file, digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def derivation_gate(ledger, pins):
    """Audit pinned declarations; this is not a semantic proof assistant."""
    nodes = {row['id']: row for row in ledger['nodes']}
    if len(nodes) != len(ledger['nodes']):
        raise ValueError('Duplicate derivation node')
    known = {('local', row['path']): row['sha256']
             for row in pins['local_inputs'] + pins['preserved_evidence']}
    for home, key in [('rkf', 'rkf'), ('publications', 'publications_comparison_only')]:
        known.update({(home, row['path']): row['sha256'] for row in pins[key]['files']})
    bindings = {(row['home'], row['path'], row['commit'], row['sha256'])
                for row in pins['source_citation_bindings']}
    for row in nodes.values():
        if 'source' in row:
            src = row['source']
            if known.get((src['home'], src['path'])) != src['sha256']:
                raise ValueError('Unpinned source: ' + row['id'])
            if (src['home'], src['path'], src['commit'], src['sha256']) not in bindings:
                raise ValueError('Changed source citation: ' + row['id'])
    roots = ledger['certified_claims']
    if roots != [f'R23.{n}' for n in range(1, 10)]:
        raise ValueError('Changed certified claim set')
    allowed = {'NATIVE_SOURCE', 'NATIVE_CONSTRUCTION', 'NATIVE_TARGET_CONSTRUCTION', 'NATIVE_DERIVED'}
    visited, active = set(), set()

    def walk(node_id):
        if node_id not in nodes:
            raise ValueError('Missing dependency: ' + node_id)
        if node_id in active:
            raise ValueError('Cyclic derivation')
        if node_id in visited:
            return
        row = nodes[node_id]
        if row['kind'] not in allowed:
            raise ValueError('Forbidden proof premise: ' + row['kind'])
        if row.get('claim_class') == 'PHYSICAL_SELECTION':
            raise ValueError('Physical selection is not proved in R23')
        if row['kind'] == 'NATIVE_TARGET_CONSTRUCTION' and row.get('physical_selection_proved') is not False:
            raise ValueError('Target construction promoted to physical selection')
        if row['kind'] == 'NATIVE_SOURCE' and 'source' not in row:
            raise ValueError('Native source missing citation')
        active.add(node_id)
        for dep in row['depends_on']:
            walk(dep)
        active.remove(node_id)
        visited.add(node_id)

    for root in roots:
        walk(root)
    return {'status': 'PASS_DECLARED_NATIVE_DERIVATION_GRAPH',
            'certified_claims': len(roots), 'proof_path_nodes': len(visited),
            'imported_classical_premises_on_proof_paths': 0,
            'comparison_references_excluded_from_proof_paths': sorted(
                key for key, row in nodes.items() if row['kind'] == 'COMPARISON_ONLY'),
            'target_definitions_retained_as_definitions': sorted(
                key for key in visited if nodes[key]['kind'] == 'NATIVE_TARGET_CONSTRUCTION'),
            'automatic_semantic_or_formal_proof_verification': False}


def mutation_checks(ledger, pins):
    checks = {}

    def reject(name, change):
        changed = copy.deepcopy(ledger)
        change(changed)
        try:
            derivation_gate(changed, pins)
        except ValueError:
            checks[name] = True
        else:
            raise AssertionError('Invalid graph accepted: ' + name)

    def node(book, name):
        return next(row for row in book['nodes'] if row['id'] == name)

    for kind in ['IMPORTED_CLASSICAL', 'ADMITTED_INPUT', 'OPEN', 'COMPARISON_ONLY']:
        def plant(book, category=kind):
            book['nodes'].append({'id': 'planted', 'kind': category, 'depends_on': []})
            node(book, 'R23.1')['depends_on'].append('planted')
        reject('reject_' + kind.lower() + '_premise', plant)
    reject('reject_physical_claim_promotion', lambda book:
           node(book, 'R23.1').update(claim_class='PHYSICAL_SELECTION'))
    reject('reject_protocol_definition_promotion', lambda book:
           node(book, 'allocation_words').update(physical_selection_proved=True))
    reject('reject_dependency_cycle', lambda book:
           node(book, 'R23.1')['depends_on'].append('R23.1'))
    reject('reject_unpinned_source', lambda book:
           node(book, 'source_roles')['source'].update(sha256='0' * 64))
    reject('reject_wrong_citation_commit', lambda book:
           node(book, 'source_roles')['source'].update(commit='0' * 40))
    return checks


def main():
    if not __debug__:
        raise RuntimeError('Certificate assertions require Python without -O.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rkf-root', required=True, type=Path)
    ap.add_argument('--publications-root', required=True, type=Path)
    ap.add_argument('--output', type=Path, default=HERE / 'R23_VERIFICATION.json')
    ap.add_argument('--native-output', type=Path, default=HERE / 'R23_NATIVE_CERTIFICATE.json')
    args = ap.parse_args()
    pin_bytes = (HERE / 'R23_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    for row in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, row)
    for key, home in [('rkf', args.rkf_root), ('publications_comparison_only', args.publications_root)]:
        for row in pins[key]['files']:
            check_file(home, row)
    proof = (ROOT / '02-relational-response/FUTURE_RESPONSE_QUOTIENT_R23.md').read_text()
    ledger = json.loads((HERE / 'R23_DERIVATION_LEDGER.json').read_text())
    sections = list(re.finditer(r'### R23\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)', proof, re.S))
    assert len(sections) == len(ledger['claims']) == 9
    nodes = {row['id']: row for row in ledger['nodes']}
    for n, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        assert int(section[1]) == n and claim['id'] == f'R23.{n}'
        assert claim['title'] == section[2] and '**Proof.**' in section[3]
        assert claim['written_section_sha256'] == digest(section[0].encode())
        assert nodes[claim['id']]['written_section_sha256'] == claim['written_section_sha256']
        assert claim['evidence'] == 'WRITTEN_NATIVE_PROOF_WITH_SCOPED_EXACT_CHECKS'
    assert ledger['no_classical_premise_rule'] is True
    assert ledger['physical_observer_selected'] is False
    graph, graph_negatives = derivation_gate(ledger, pins), mutation_checks(ledger, pins)

    with tempfile.TemporaryDirectory(prefix='r23-source-') as tmp:
        report_path, native_path = Path(tmp) / 'r22.json', Path(tmp) / 'r22-native.json'
        subprocess.run([sys.executable, '-B', str(HERE / 'verify_r22.py'),
                        '--rkf-root', str(args.rkf_root), '--publications-root', str(args.publications_root),
                        '--output', str(report_path), '--native-output', str(native_path)],
                       check=True, capture_output=True, text=True)
        previous, frozen = json.loads(report_path.read_text()), json.loads((HERE / 'R22_VERIFICATION.json').read_text())
        for key in ['status', 'source_pins_sha256', 'native_certificate_sha256', 'previous_source_replayed',
                    'counters', 'derivation_graph', 'derivation_graph_mutations_rejected']:
            assert previous[key] == frozen[key], 'Frozen R22 replay mismatch: ' + key
        assert digest(native_path.read_bytes()) == frozen['native_certificate_sha256']

        row = next(row for row in pins['local_inputs'] if row['path'].endswith('future_response_quotient.cjs'))
        altered = Path(tmp) / row['path']
        altered.parent.mkdir(parents=True, exist_ok=True)
        altered.write_bytes((ROOT / row['path']).read_bytes() + b'\n// altered\n')
        try:
            check_file(Path(tmp), row)
        except ValueError:
            altered_source_rejected = True
        else:
            raise AssertionError('Altered source accepted')
    assert previous['status'] == 'PASS_R22_NATIVE_MEETING_GEOMETRY'

    command = ['node', str(HERE / 'future_response_quotient.cjs'),
               '--rkf-root', str(args.rkf_root), '--input', str(HERE / 'R17_NATIVE_INPUT.json'),
               '--expected-input-sha256', pins['native_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_R23_NATIVE_FUTURE_RESPONSE'
    assert native['check_count'] == 10 and all(row['passed'] for row in native['exact_checks'])
    assert native['symbolic_replay_count'] == 12
    assert all(row['replay'] == 'REPLAY_MATCH' for row in native['symbolic_replays'])
    check_names = {row['name'] for row in native['exact_checks']}
    assert all(set(claim['exact_check_groups']) <= check_names for claim in ledger['claims'])
    assert len(native['rejected_false_alternatives']) == 16 and all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']
    assert native['native_future_response_quotient_derived'] and native['native_record_probe_completion_derived']
    for flag in ['ordinary_complex_measurement_stochastic_or_classical_observer_premise',
                 'physical_observer_selected', 'physical_metric_c_alpha_derived',
                 'formal_proof_assistant_verified']:
        assert native[flag] is False
    wrong = list(command)
    wrong[-1] = '0' * 64
    result = subprocess.run(wrong, capture_output=True, text=True)
    assert result.returncode != 0 and 'input pin mismatch' in result.stderr
    for row in pins['preserved_evidence']:
        check_file(ROOT, row)

    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False) + '\n').encode()
    report = {
        'schema': 'extra-ideas.r23.verification.v1', 'date': '2026-10-01',
        'status': 'PASS_R23_NATIVE_FUTURE_RESPONSE',
        'written_results': 9, 'exact_check_groups': 10, 'native_symbolic_replays': 12,
        'counters': {row['name']: {k: v for k, v in row.items() if k not in ('name', 'passed')}
                     for row in native['exact_checks']},
        'rejected_false_alternatives': native['rejected_false_alternatives'],
        'derivation_graph': graph, 'derivation_graph_mutations_rejected': graph_negatives,
        'altered_source_rejected': altered_source_rejected,
        'altered_native_certificate_rejected': True, 'wrong_input_pin_rejected': True,
        'previous_source_replayed': {'commit': pins['source_commit'], 'status': previous['status'],
            'native_certificate_sha256': previous['native_certificate_sha256'],
            'old_certificate_or_source_pin_modified': False,
            'earlier_source_chain': previous['previous_source_replayed']},
        'upstream_canonical_commit': pins['rkf']['commit'],
        'upstream_files_verified': len(pins['rkf']['files']),
        'preserved_prior_files': len(pins['preserved_evidence']),
        'source_pins_sha256': digest(pin_bytes), 'native_certificate_sha256': digest(native_bytes),
        'comparison_references_are_proof_dependencies': False,
        'publications_comparison_certificates_rerun': False,
        'native_future_response_quotient_derived': True, 'native_record_probe_completion_derived': True,
        'physical_observer_selected': False,
        'physical_metric_c_alpha_derived': False,
        'formal_proof_assistant_verified': False,
        'runtime': {'python': sys.version.split()[0], 'node': subprocess.run(
            ['node', '--version'], check=True, capture_output=True, text=True).stdout.strip()},
        'scope': 'Native translation symmetry and future-response quotient on the real signed pair target. Complete exchange-probe rank, minimum finite-catalogue observer, initial-role return, endpoint decoding and native record-probe completion. Meeting geometry does not determine signed response. Written proofs and exact finite checks; no physical observer or metric selection.'}
    report['canonical_report_sha256'] = digest(canonical(report))
    args.native_output.parent.mkdir(parents=True, exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(report['status'])
    print('9 written results; 10 exact groups; 12 native replays; 16 false alternatives and 9 graph mutations rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
