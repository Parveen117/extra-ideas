#!/usr/bin/env python3
"""Replay R28's cut-curvature exchange, hidden memory and signed propagation and its immutable source chain."""
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
    if roots != [f'R28.{n}' for n in range(1, 7)]:
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
            raise ValueError('Physical selection is not proved in R28')
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
            node(book, 'R28.1')['depends_on'].append('planted')
        reject('reject_' + kind.lower() + '_premise', plant)
    reject('reject_physical_claim_promotion', lambda book:
           node(book, 'R28.1').update(claim_class='PHYSICAL_SELECTION'))
    reject('reject_protocol_definition_promotion', lambda book:
           node(book, 'curvature_response_target').update(physical_selection_proved=True))
    reject('reject_dependency_cycle', lambda book:
           node(book, 'R28.1')['depends_on'].append('R28.1'))
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
    ap.add_argument('--output', type=Path, default=HERE / 'R28_VERIFICATION.json')
    ap.add_argument('--native-output', type=Path, default=HERE / 'R28_NATIVE_CERTIFICATE.json')
    args = ap.parse_args()
    pin_bytes = (HERE / 'R28_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    for row in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, row)
    for key, home in [('rkf', args.rkf_root), ('publications_comparison_only', args.publications_root)]:
        for row in pins[key]['files']:
            check_file(home, row)
    proof = (ROOT / '02-relational-response/NATIVE_CUT_CURVATURE_PROPAGATION_R28.md').read_text()
    ledger = json.loads((HERE / 'R28_DERIVATION_LEDGER.json').read_text())
    sections = list(re.finditer(r'### R28\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)', proof, re.S))
    assert len(sections) == len(ledger['claims']) == 6
    nodes = {row['id']: row for row in ledger['nodes']}
    for n, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        assert int(section[1]) == n and claim['id'] == f'R28.{n}'
        assert claim['title'] == section[2] and '**Proof.**' in section[3]
        assert claim['written_section_sha256'] == digest(section[0].encode())
        assert nodes[claim['id']]['written_section_sha256'] == claim['written_section_sha256']
        assert claim['evidence'] == 'WRITTEN_NATIVE_PROOF_WITH_SCOPED_EXACT_CHECKS'
    assert ledger['no_classical_premise_rule']
    for flag in ['physical_vacuum_identified','electromagnetism_derived','physical_metric_c_alpha_derived','formal_proof_assistant_verified']:
        assert ledger[flag] is False
    graph, graph_negatives = derivation_gate(ledger, pins), mutation_checks(ledger, pins)

    with tempfile.TemporaryDirectory(prefix='r28-source-') as tmp:
        report_path, native_path = Path(tmp) / 'r27.json', Path(tmp) / 'r27-native.json'
        subprocess.run([sys.executable, '-B', str(HERE / 'verify_r27.py'),
                        '--rkf-root', str(args.rkf_root), '--publications-root', str(args.publications_root),
                        '--output', str(report_path), '--native-output', str(native_path)],
                       check=True, capture_output=True, text=True)
        previous, frozen = json.loads(report_path.read_text()), json.loads((HERE / 'R27_VERIFICATION.json').read_text())
        for key in ['status', 'source_pins_sha256', 'native_certificate_sha256', 'previous_source_replayed',
                    'counters', 'derivation_graph', 'derivation_graph_mutations_rejected']:
            assert previous[key] == frozen[key], 'Frozen R27 replay mismatch: ' + key
        assert digest(native_path.read_bytes()) == frozen['native_certificate_sha256']
        row = next(row for row in pins['local_inputs'] if row['path'].endswith('native_cut_curvature_propagation.cjs'))
        altered = Path(tmp) / row['path']
        altered.parent.mkdir(parents=True, exist_ok=True)
        altered.write_bytes((ROOT / row['path']).read_bytes() + b'\n// altered\n')
        try:
            check_file(Path(tmp), row)
        except ValueError:
            altered_source_rejected = True
        else:
            raise AssertionError('Altered source accepted')
    assert previous['status'] == 'PASS_R27_NATIVE_REPLICA_SELECTION'
    command = ['node', str(HERE / 'native_cut_curvature_propagation.cjs'),
               '--rkf-root', str(args.rkf_root), '--input', str(HERE / 'R17_NATIVE_INPUT.json'),
               '--expected-input-sha256', pins['native_input_sha256'],
               '--r27-certificate', str(HERE / 'R27_NATIVE_CERTIFICATE.json'),
               '--expected-r27-sha256', pins['r27_native_certificate_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_R28_NATIVE_CUT_CURVATURE_PROPAGATION'
    assert native['check_count'] == 8 and all(row['passed'] for row in native['exact_checks'])
    assert native['symbolic_replay_count'] == 10
    assert all(row['replay'] == 'REPLAY_MATCH' for row in native['symbolic_replays'])
    check_names = {row['name'] for row in native['exact_checks']}
    assert all(set(claim['exact_check_groups']) <= check_names for claim in ledger['claims'])
    assert len(native['rejected_false_alternatives']) == 9 and all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']
    for flag in ['native_curvature_exchange_derived','hidden_memory_kernel_derived','signed_curvature_response_wave_derived','native_chart_free_cyclic_ledger_derived']:
        assert native[flag] is True
    generic = next(row for row in native['symbolic_presentations'] if row['name']=='native_associative_arrows_and_involution')
    assert generic['specification']['tokens'] == ['A','B','C','J']
    assert [row['lhs'] for row in generic['specification']['rules']] == [['J','J']]
    assert generic['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS'
    for flag in ['ordinary_complex_or_classical_field_premise','physical_vacuum_identified','electromagnetism_derived','physical_metric_c_alpha_derived','formal_proof_assistant_verified']:
        assert native[flag] is False
    for option, message in [('--expected-input-sha256','input pin mismatch'),('--expected-r27-sha256','R27 certificate pin mismatch')]:
        wrong = list(command)
        wrong[wrong.index(option)+1] = '0'*64
        result = subprocess.run(wrong, capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr
    assert native['r27_native_certificate_sha256'] == pins['r27_native_certificate_sha256']
    for row in pins['preserved_evidence']:
        check_file(ROOT, row)
    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False)+'\n').encode()
    report = {
        'schema':'extra-ideas.r28.verification.v1','date':'2026-10-01',
        'status':'PASS_R28_NATIVE_CUT_CURVATURE_PROPAGATION',
        'written_results':6,'exact_check_groups':8,'native_symbolic_replays':10,
        'counters':{row['name']:{k:v for k,v in row.items() if k not in ('name','passed')} for row in native['exact_checks']},
        'rejected_false_alternatives':native['rejected_false_alternatives'],
        'derivation_graph':graph,'derivation_graph_mutations_rejected':graph_negatives,
        'altered_source_rejected':altered_source_rejected,'altered_native_certificate_rejected':True,
        'wrong_input_pin_rejected':True,'wrong_r27_certificate_pin_rejected':True,
        'previous_source_replayed':{'commit':pins['source_commit'],'status':previous['status'],
            'native_certificate_sha256':previous['native_certificate_sha256'],'old_certificate_or_source_pin_modified':False,
            'earlier_source_chain':previous['previous_source_replayed']},
        'upstream_canonical_commit':pins['rkf']['commit'],'upstream_files_verified':len(pins['rkf']['files']),
        'preserved_prior_files':len(pins['preserved_evidence']),
        'source_pins_sha256':digest(pin_bytes),'native_certificate_sha256':digest(native_bytes),
        'comparison_references_are_proof_dependencies':False,'comparison_source_certificates_rerun':False,
        'native_curvature_exchange_derived':True,'hidden_memory_kernel_derived':True,
        'signed_curvature_response_wave_derived':True,'native_chart_free_cyclic_ledger_derived':True,
        'physical_vacuum_identified':False,'electromagnetism_derived':False,'physical_metric_c_alpha_derived':False,
        'formal_proof_assistant_verified':False,
        'runtime':{'python':sys.version.split()[0],'node':subprocess.run(['node','--version'],check=True,capture_output=True,text=True).stdout.strip()},
        'scope':'Six written proofs bind native cut curvature, transition exchange, exact hidden-memory kernel, signed response wave, chart-free cyclic ledger and non-uniqueness of dynamics from algebra alone. Eight finite exact groups and ten native word replays support this scope; no physical vacuum/EM/c/alpha or proof-assistant claim.'}
    report['canonical_report_sha256'] = digest(canonical(report))
    args.native_output.parent.mkdir(parents=True,exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(report['status'])
    print('6 written proofs; 8 exact groups; 10 native replays; 9 false alternatives and 9 graph mutations rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
