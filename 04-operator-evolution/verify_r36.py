#!/usr/bin/env python3
"""Replay R36 native curvature-observer geometry and current with the frozen source chain."""
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
    if roots != [f'R36.{n}' for n in range(1, 8)]:
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
            raise ValueError('Physical selection is not proved in R36')
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
            node(book, 'R36.1')['depends_on'].append('planted')
        reject('reject_' + kind.lower() + '_premise', plant)
    reject('reject_physical_claim_promotion', lambda book:
           node(book, 'R36.1').update(claim_class='PHYSICAL_SELECTION'))
    reject('reject_protocol_definition_promotion', lambda book:
           node(book, 'native_curvature_observer_target').update(physical_selection_proved=True))
    reject('reject_dependency_cycle', lambda book:
           node(book, 'R36.1')['depends_on'].append('R36.1'))
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
    ap.add_argument('--output', type=Path, default=HERE / 'R36_VERIFICATION.json')
    ap.add_argument('--native-output', type=Path, default=HERE / 'R36_NATIVE_CERTIFICATE.json')
    args = ap.parse_args()
    pin_bytes = (HERE / 'R36_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    for row in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, row)
    for key, home in [('rkf', args.rkf_root), ('publications_comparison_only', args.publications_root)]:
        for row in pins[key]['files']:
            check_file(home, row)
    proof = (ROOT / '02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md').read_text()
    ledger = json.loads((HERE / 'R36_DERIVATION_LEDGER.json').read_text())
    sections = list(re.finditer(r'### R36\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)', proof, re.S))
    assert len(sections) == len(ledger['claims']) == 7
    nodes = {row['id']: row for row in ledger['nodes']}
    for n, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        assert int(section[1]) == n and claim['id'] == f'R36.{n}'
        assert claim['title'] == section[2] and '**Proof.**' in section[3]
        assert claim['written_section_sha256'] == digest(section[0].encode())
        assert nodes[claim['id']]['written_section_sha256'] == claim['written_section_sha256']
        assert claim['evidence'] == 'WRITTEN_NATIVE_PROOF_WITH_SCOPED_EXACT_CHECKS'
    assert ledger['no_classical_premise_rule']
    reading = ledger['inherited_combined_reading_register']
    assert [row['paper'] for row in reading] == [f'R{n}' for n in range(1, 34)]
    assert reading == json.loads((HERE / 'R35_DERIVATION_LEDGER.json').read_text())['combined_reading_register']
    preserved_hashes = {row['path']: row['sha256'] for row in pins['preserved_evidence']}
    for entry in reading:
        assert entry['source_commit'] == ledger['inherited_reading_source_commit'] and entry['files']
        assert entry['status'] == 'READ_FOR_SYNTHESIS_NOT_AUTOMATIC_PROOF_PREMISE'
        for row in entry['files']:
            assert preserved_hashes[row['path']] == row['sha256']
    for flag in ['physical_vacuum_identified','electromagnetism_derived','physical_metric_c_alpha_derived','physical_electric_charge_or_particle_mass_identified','physical_spatial_dimension_or_gauge_group_selected','primitive_physical_force_derived','formal_proof_assistant_verified']:
        assert ledger[flag] is False
    graph, graph_negatives = derivation_gate(ledger, pins), mutation_checks(ledger, pins)

    with tempfile.TemporaryDirectory(prefix='r36-source-') as tmp:
        report_path, native_path = Path(tmp) / 'r35.json', Path(tmp) / 'r35-native.json'
        subprocess.run([sys.executable, '-B', str(HERE / 'verify_r35.py'),
                        '--rkf-root', str(args.rkf_root), '--publications-root', str(args.publications_root),
                        '--output', str(report_path), '--native-output', str(native_path)],
                       check=True, capture_output=True, text=True)
        previous, frozen = json.loads(report_path.read_text()), json.loads((HERE / 'R35_VERIFICATION.json').read_text())
        for key in ['status', 'source_pins_sha256', 'native_certificate_sha256', 'previous_source_replayed',
                    'counters', 'derivation_graph', 'derivation_graph_mutations_rejected']:
            assert previous[key] == frozen[key], 'Frozen R35 replay mismatch: ' + key
        assert digest(native_path.read_bytes()) == frozen['native_certificate_sha256']
        row = next(row for row in pins['local_inputs'] if row['path'].endswith('native_curvature_observer.cjs'))
        altered = Path(tmp) / row['path']
        altered.parent.mkdir(parents=True, exist_ok=True)
        altered.write_bytes((ROOT / row['path']).read_bytes() + b'\n// altered\n')
        try:
            check_file(Path(tmp), row)
        except ValueError:
            altered_source_rejected = True
        else:
            raise AssertionError('Altered source accepted')
    assert previous['status'] == 'PASS_R35_NATIVE_SIGNED_ENVELOPE'
    command = ['node', str(HERE / 'native_curvature_observer.cjs'),
               '--rkf-root', str(args.rkf_root), '--input', str(HERE / 'R17_NATIVE_INPUT.json'),
               '--expected-input-sha256', pins['native_input_sha256'],
               '--r35-certificate', str(HERE / 'R35_NATIVE_CERTIFICATE.json'),
               '--expected-r35-sha256', pins['r35_native_certificate_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_R36_NATIVE_CURVATURE_OBSERVER'
    assert native['check_count'] == 9 and all(row['passed'] for row in native['exact_checks'])
    assert native['symbolic_replay_count'] == 14
    assert all(row['replay'] == 'REPLAY_MATCH' for row in native['symbolic_replays'])
    check_names = {row['name'] for row in native['exact_checks']}
    assert all(set(claim['exact_check_groups']) <= check_names for claim in ledger['claims'])
    assert len(native['rejected_false_alternatives']) == 15 and all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']
    for flag in ['native_response_metric_and_area_curvature_derived','native_minimal_curvature_observer_derived','native_information_current_cone_derived','native_observed_dynamics_and_limit_derived']:
        assert native[flag] is True
    canonical_source = native['symbolic_presentations'][0]
    assert canonical_source['name'] == 'canonical_native_source'
    assert set(canonical_source['specification']['tokens']) == {'R','K'}
    assert canonical_source['audit']['status'] == 'CONFLUENT_BY_CHECKED_DIAMONDS'
    for flag in ['ordinary_complex_or_classical_field_premise','physical_vacuum_identified','electromagnetism_derived','physical_metric_c_alpha_derived','physical_electric_charge_or_particle_mass_identified','physical_spatial_dimension_or_gauge_group_selected','primitive_physical_force_derived','formal_proof_assistant_verified']:
        assert native[flag] is False
    for option, message in [('--expected-input-sha256','input pin mismatch'),('--expected-r35-sha256','R35 certificate pin mismatch')]:
        wrong = list(command)
        wrong[wrong.index(option)+1] = '0'*64
        result = subprocess.run(wrong, capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr
    assert native['r35_native_certificate_sha256'] == pins['r35_native_certificate_sha256']
    for row in pins['preserved_evidence']:
        check_file(ROOT, row)
    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False)+'\n').encode()
    report = {
        'schema':'extra-ideas.r36.verification.v1','date':'2026-10-02',
        'status':'PASS_R36_NATIVE_CURVATURE_OBSERVER',
        'inherited_combined_reading_papers':len(reading),
        'written_results':7,'exact_check_groups':9,'native_symbolic_replays':14,
        'counters':{row['name']:{k:v for k,v in row.items() if k not in ('name','passed')} for row in native['exact_checks']},
        'rejected_false_alternatives':native['rejected_false_alternatives'],
        'derivation_graph':graph,'derivation_graph_mutations_rejected':graph_negatives,
        'altered_source_rejected':altered_source_rejected,'altered_native_certificate_rejected':True,
        'wrong_input_pin_rejected':True,'wrong_r35_certificate_pin_rejected':True,
        'previous_source_replayed':{'commit':pins['source_commit'],'status':previous['status'],
            'native_certificate_sha256':previous['native_certificate_sha256'],'old_certificate_or_source_pin_modified':False,
            'earlier_source_chain':previous['previous_source_replayed']},
        'upstream_canonical_commit':pins['rkf']['commit'],'upstream_files_verified':len(pins['rkf']['files']),
        'preserved_prior_files':len(pins['preserved_evidence']),
        'source_pins_sha256':digest(pin_bytes),'native_certificate_sha256':digest(native_bytes),
        'comparison_references_are_proof_dependencies':False,'comparison_source_certificates_rerun':False,
        'native_response_metric_and_area_curvature_derived':True,'native_minimal_curvature_observer_derived':True,
        'native_information_current_cone_derived':True,'native_observed_dynamics_and_limit_derived':True,
        'physical_vacuum_identified':False,'electromagnetism_derived':False,'physical_metric_c_alpha_derived':False,
        'physical_electric_charge_or_particle_mass_identified':False,'physical_spatial_dimension_or_gauge_group_selected':False,'primitive_physical_force_derived':False,'formal_proof_assistant_verified':False,
        'runtime':{'python':sys.version.split()[0],'node':subprocess.run(['node','--version'],check=True,capture_output=True,text=True).stdout.strip()},
        'scope':'Seven written proofs: native response metric and oriented curvature, minimal complete curvature observer, exact paired dynamics and memory, sharp information-current cone and coherence deficit, controlled modewise flow, and frame/Riemann/physical boundaries. The constructed fields are readouts of one source, not independently identified physical fields.'}
    report['canonical_report_sha256'] = digest(canonical(report))
    args.native_output.parent.mkdir(parents=True,exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(report['status'])
    print('7 written proofs; 9 exact groups; 14 native replays; 15 false alternatives and 9 graph mutations rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
