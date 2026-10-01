#!/usr/bin/env python3
"""Verify native current retention and enforce its declared derivation graph."""
import argparse
import copy
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from verify_r18 import canonical, check_file, digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def derivation_gate(ledger, pins):
    """Audit declared dependency origins, not arbitrary mathematical prose."""
    nodes = {row['id']: row for row in ledger['nodes']}
    if len(nodes) != len(ledger['nodes']):
        raise ValueError('Duplicate derivation node')
    whitelisted = {
        ('local', row['path']): row['sha256']
        for row in pins['local_inputs'] + pins['preserved_evidence']
    }
    for home, key in [('rkf', 'rkf'), ('publications', 'publications_comparison_only')]:
        whitelisted.update({(home, row['path']): row['sha256'] for row in pins[key]['files']})
    roots = ledger['certified_claims']
    if roots != [f'R19.{i}' for i in range(1, 10)]:
        raise ValueError('Incomplete or changed certified claim set')
    for row in nodes.values():
        if 'source' in row:
            source = row['source']
            key = (source['home'], source['path'])
            if whitelisted.get(key) != source['sha256']:
                raise ValueError('Unpinned source node: ' + row['id'])
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
            raise ValueError('Forbidden proof dependency: ' + row['kind'])
        if row.get('claim_class') == 'PHYSICAL_SELECTION':
            # This packet contains no physical-selection proof. It may not
            # promote a target definition by changing a status field.
            raise ValueError('Physical selection remains open in R19')
        if row['kind'] == 'NATIVE_TARGET_CONSTRUCTION' and row.get('physical_selection_proved') is not False:
            raise ValueError('Target construction promoted without derivation')
        if row['kind'] == 'NATIVE_SOURCE' and 'source' not in row:
            raise ValueError('Native source missing pin')
        active.add(node_id)
        for dep in row['depends_on']:
            walk(dep)
        active.remove(node_id)
        visited.add(node_id)

    for root in roots:
        walk(root)
    return {
        'status': 'PASS_DECLARED_NATIVE_DERIVATION_GRAPH',
        'certified_claims': len(roots), 'proof_path_nodes': len(visited),
        'imported_classical_premises_on_proof_paths': 0,
        'comparison_references_excluded_from_proof_paths': sorted(
            key for key, row in nodes.items() if row['kind'] == 'COMPARISON_ONLY'),
        'target_definitions_retained_as_definitions': sorted(
            key for key in visited if nodes[key]['kind'] == 'NATIVE_TARGET_CONSTRUCTION'),
        'automatic_semantic_or_formal_proof_verification': False,
    }


def mutated_graph_checks(ledger, pins):
    results = {}
    for kind in ['IMPORTED_CLASSICAL', 'OPEN', 'COMPARISON_ONLY']:
        changed = copy.deepcopy(ledger)
        changed['nodes'].append({'id': 'planted', 'kind': kind, 'depends_on': []})
        next(row for row in changed['nodes'] if row['id'] == 'R19.1')['depends_on'].append('planted')
        try:
            derivation_gate(changed, pins)
        except ValueError:
            results['reject_' + kind.lower() + '_proof_dependency'] = True
        else:
            raise AssertionError('Forbidden dependency accepted: ' + kind)
    changed = copy.deepcopy(ledger)
    next(row for row in changed['nodes'] if row['id'] == 'R19.1')['claim_class'] = 'PHYSICAL_SELECTION'
    try:
        derivation_gate(changed, pins)
    except ValueError:
        results['reject_physical_selection_promotion'] = True
    else:
        raise AssertionError('Physical promotion accepted')
    changed = copy.deepcopy(ledger)
    next(row for row in changed['nodes'] if row['id'] == 'R19.1')['depends_on'].append('R19.1')
    try:
        derivation_gate(changed, pins)
    except ValueError:
        results['reject_derivation_cycle'] = True
    else:
        raise AssertionError('Derivation cycle accepted')
    changed = copy.deepcopy(ledger)
    next(row for row in changed['nodes'] if row['id'] == 'source_roles')['source']['sha256'] = '0' * 64
    try:
        derivation_gate(changed, pins)
    except ValueError:
        results['reject_unpinned_native_source'] = True
    else:
        raise AssertionError('Unpinned source accepted')
    return results


def main():
    if not __debug__:
        raise RuntimeError('Certificate assertions require Python without -O.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rkf-root', required=True, type=Path)
    ap.add_argument('--publications-root', required=True, type=Path)
    ap.add_argument('--output', type=Path, default=HERE / 'R19_VERIFICATION.json')
    ap.add_argument('--native-output', type=Path, default=HERE / 'R19_NATIVE_CERTIFICATE.json')
    args = ap.parse_args()
    pin_bytes = (HERE / 'R19_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    for row in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, row)
    for key, home in [('rkf', args.rkf_root), ('publications_comparison_only', args.publications_root)]:
        for row in pins[key]['files']:
            check_file(home, row)
    proof = (ROOT / '02-relational-response/CURRENT_MEMORY_RETENTION_R19.md').read_text()
    ledger = json.loads((HERE / 'R19_DERIVATION_LEDGER.json').read_text())
    pattern = r'### R19\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)'
    sections = list(re.finditer(pattern, proof, re.S))
    assert len(sections) == len(ledger['claims']) == 9
    nodes = {row['id']: row for row in ledger['nodes']}
    for n, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        assert claim['id'] == f'R19.{n}' and int(section[1]) == n
        assert claim['title'] == section[2] and '**Proof.**' in section[3]
        assert claim['written_section_sha256'] == digest(section[0].encode())
        assert nodes[claim['id']]['written_section_sha256'] == claim['written_section_sha256']
        assert claim['evidence'] == 'WRITTEN_NATIVE_PROOF_WITH_SCOPED_EXACT_CHECKS'
    assert ledger['physical_interaction_selection_proved'] is False
    assert ledger['no_classical_premise_rule'] is True
    graph = derivation_gate(ledger, pins)
    graph_negative = mutated_graph_checks(ledger, pins)

    historical = json.loads((HERE / 'R19_HISTORICAL_REPLAY_INPUTS.json').read_text())
    old_pins = json.loads((HERE / 'R17_SOURCE_PINS.json').read_text())
    nav_paths = {'README.md', '02-relational-response/README.md', '04-operator-evolution/README.md'}
    assert historical['commit'] == 'f561681959f5af93389db6bc77763c6794d76a87'
    assert {row['path'] for row in historical['files']} == nav_paths
    old_nav = {row['path']: row for row in old_pins['local_inputs'] if row['path'] in nav_paths}
    assert len(old_nav) == 3
    for row in historical['files']:
        assert row['sha256'] == old_nav[row['path']]['sha256'] == digest(row['content'].encode())
        assert row['git_blob_sha'] == old_nav[row['path']]['git_blob_sha']

    with tempfile.TemporaryDirectory(prefix='r19-source-') as tmp:
        # R17 pinned its three then-current navigation READMEs. R18 changed
        # those files but its old nested verifier reused the active checkout.
        # Reproduce the immutable proof inputs in isolation instead of changing
        # either the old pins or the user's current navigation.
        snapshot = Path(tmp) / 'source-snapshot'
        parent_pins = json.loads((HERE / 'R18_SOURCE_PINS.json').read_text())
        required = {row['path'] for row in parent_pins['local_inputs'] + parent_pins['preserved_evidence']}
        required.add('04-operator-evolution/R18_SOURCE_PINS.json')
        for path in required:
            dest = snapshot / path
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / path, dest)
        for row in historical['files']:
            dest = snapshot / row['path']
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(row['content'].encode())
            check_file(snapshot, old_nav[row['path']])
        report_path, native_path = Path(tmp) / 'r18.json', Path(tmp) / 'r18-native.json'
        subprocess.run([sys.executable, '-B', str(snapshot / '04-operator-evolution/verify_r18.py'),
                        '--rkf-root', str(args.rkf_root), '--publications-root', str(args.publications_root),
                        '--output', str(report_path), '--native-output', str(native_path)],
                       check=True, capture_output=True, text=True)
        previous = json.loads(report_path.read_text())
        frozen = json.loads((HERE / 'R18_VERIFICATION.json').read_text())
        for key in ['status', 'source_pins_sha256', 'native_certificate_sha256',
                    'previous_source_replayed', 'counters']:
            assert previous[key] == frozen[key], 'R18 replay mismatch: ' + key
    assert previous['status'] == 'PASS_R18_CUT_ADDRESS_TRANSPORT'
    assert not previous['physical_metric_selection_proved']

    command = ['node', str(HERE / 'current_memory_retention.cjs'),
               '--rkf-root', str(args.rkf_root), '--input', str(HERE / 'R17_NATIVE_INPUT.json'),
               '--expected-input-sha256', pins['native_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_R19_NATIVE_CURRENT_MEMORY'
    assert native['check_count'] == 11 and all(row['passed'] for row in native['exact_checks'])
    assert native['symbolic_replay_count'] == 9
    assert all(row['replay'] == 'REPLAY_MATCH' for row in native['symbolic_replays'])
    assert len(native['rejected_false_alternatives']) == 12
    assert all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']
    assert not native['ordinary_complex_hilbert_probability_or_physics_law_input']
    assert not native['physical_interaction_selection_proved']
    assert not native['formal_proof_assistant_verified']
    wrong = list(command)
    wrong[-1] = '0' * 64
    result = subprocess.run(wrong, capture_output=True, text=True)
    assert result.returncode != 0 and 'input pin mismatch' in result.stderr
    for row in pins['preserved_evidence']:
        check_file(ROOT, row)

    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False) + '\n').encode()
    report = {
        'schema': 'extra-ideas.r19.verification.v1', 'date': '2026-10-01',
        'status': 'PASS_R19_NATIVE_CURRENT_RETENTION',
        'written_results': 9, 'exact_check_groups': 11, 'native_symbolic_replays': 9,
        'counters': {row['name']: {k: v for k, v in row.items() if k not in ('name', 'passed')}
                     for row in native['exact_checks']},
        'rejected_false_alternatives': native['rejected_false_alternatives'],
        'derivation_graph': graph, 'derivation_graph_mutations_rejected': graph_negative,
        'altered_native_certificate_rejected': True, 'wrong_input_pin_rejected': True,
        'previous_source_replayed': {'commit': pins['source_commit'], 'status': previous['status'],
            'native_certificate_sha256': previous['native_certificate_sha256'],
            'historical_navigation_inputs_restored_in_temporary_snapshot': 3,
            'old_certificate_or_source_pin_modified': False,
            'earlier_source_chain': previous['previous_source_replayed']},
        'upstream_canonical_commit': pins['rkf']['commit'],
        'upstream_files_verified': len(pins['rkf']['files']),
        'preserved_prior_files': len(pins['preserved_evidence']),
        'source_pins_sha256': digest(pin_bytes), 'native_certificate_sha256': digest(native_bytes),
        'comparison_references_are_proof_dependencies': False,
        'physical_interaction_selection_proved': False,
        'physical_metric_c_alpha_or_new_force_derived': False,
        'formal_proof_assistant_verified': False,
        'inherited_target_constructions_physically_selected': False,
        'runtime': {'python': sys.version.split()[0], 'node': subprocess.run(
            ['node', '--version'], check=True, capture_output=True, text=True).stdout.strip()},
        'scope': 'Native pair-record propagation and memory derived from the inherited R18 target. Written arbitrary-horizon proofs and exact finite checks. The origin gate audits declared dependencies; it is not an automated semantic proof or physical selection.'}
    report['canonical_report_sha256'] = digest(canonical(report))
    args.native_output.parent.mkdir(parents=True, exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(report['status'])
    print('9 written results; 11 exact check groups; 9 native replays; 12 false alternatives and 6 dependency mutations rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
