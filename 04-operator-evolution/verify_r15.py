#!/usr/bin/env python3
"""Reproduce R15 without mutating upstream sources or their certificate pins."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PACKET = HERE / 'certificates/r15'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def check_file(path, record):
    data = path.read_bytes()
    if digest(data) != record['sha256']:
        raise ValueError('SHA-256 mismatch: ' + str(path))
    if 'git_blob_sha' in record:
        git_sha = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if git_sha != record['git_blob_sha']:
            raise ValueError('Git blob mismatch: ' + str(path))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--publications-root', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=HERE / 'R15_VERIFICATION.json')
    args = parser.parse_args()
    pins = json.loads((PACKET / 'SOURCE_PINS.json').read_text())
    for record in pins['local_files']:
        check_file(ROOT / record['path'], record)
    for name, path in [('rkf', args.rkf_root), ('publications', args.publications_root)]:
        for record in pins[name]['files']:
            check_file(path / record['path'], record)
    for record in pins['preserved_previous_files']:
        check_file(ROOT / record['path'], record)

    source = (PACKET / 'emk_topology.original.tex').read_bytes()
    ledger = json.loads((PACKET / 'CLAIM_LEDGER.json').read_text())
    assert digest(source) == ledger['original_source_sha256']
    normalized = source.decode().replace('\r\n', '\n')
    pattern = r'\\begin\{(theorem|lemma|proposition|corollary)\}(?:\[([^\]]*)\])?(.*?)\\end\{\1\}'
    environments = list(re.finditer(pattern, normalized, re.S))
    assert len(environments) == len(ledger['claims']) == 48
    for i, (m, row) in enumerate(zip(environments, ledger['claims']), 1):
        assert row['id'] == f'ET{i:02}'
        assert row['environment'] == m[1] and row['title'] == (m[2] or '')
        assert row['start_line'] == normalized[:m.start()].count('\n') + 1
        assert row['end_line'] == normalized[:m.end()].count('\n') + 1
        assert row['normalized_environment_sha256'] == digest(m[0].encode())
        assert row['status'] and row['reason'] and row['replacement_or_evidence']
    assert ledger['whole_source_verdict'] == 'NOT_CERTIFIED_AS_WRITTEN'

    # Import only after checking the local verifier and implementation pins.
    import emk_topology_core as core
    finite = core.exhaustive_checks()
    assert finite['finite_systems_and_observations'] == 3678
    witnesses = core.boundary_witnesses()
    mutations = core.mutation_controls()
    phase = core.phase_operator_checks()

    command = ['node', str(HERE / 'r15_native_topology_probe.cjs'),
               '--rkf-root', str(args.rkf_root), '--expected-input-sha256', pins['native_input_sha256']]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    native = json.loads(result.stdout)
    assert native == json.loads((HERE / 'R15_NATIVE_CERTIFICATE.json').read_text())
    replays = [row for group in native['groups'].values() for row in group['results']]
    assert len(replays) == 10 and all(row['replay'] == 'REPLAY_MATCH' for row in replays)

    ugd_path = args.publications_root / 'papers/emk-ugd-algebra/certificates/ugd1_numerals.py'
    spec = importlib.util.spec_from_file_location('r15_pinned_ugd1', ugd_path)
    ugd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ugd)
    # Call the pure builder, never main(), which writes upstream result/pin files.
    ugd_cert = ugd.build_certificate()
    ugd_bytes = (json.dumps(ugd_cert, indent=2, sort_keys=True) + '\n').encode()
    ugd_folder = ugd_path.parent
    assert ugd_bytes == (ugd_folder / 'UGD1_RESULT.json').read_bytes()
    assert digest(ugd_bytes) == (ugd_folder / 'EXPECTED_UGD1.sha256').read_text().strip()
    a = ugd.Numeral({0: (7, 0)}, ledger=0)
    b = ugd.Numeral({0: (7, 0)}, ledger=2)
    assert a.classical_projection() == b.classical_projection() == 7
    assert a != b and a.total_seam() == 0 and b.total_seam() == 2

    report = {
        'schema': 'extra-ideas.r15.verification.v1', 'date': '2026-10-01',
        'status': 'SCOPED_CHECKS_PASS_ORIGINAL_NOT_CERTIFIED',
        'whole_original_manuscript_certified': False,
        'formal_proof_assistant': False,
        'scope': 'Written constructive proofs, exact finite checks, original-claim audit, native operator replay and UGD-1 reuse. No three-axiom universal sufficiency or physical validation.',
        'source_pins_sha256': digest((PACKET / 'SOURCE_PINS.json').read_bytes()),
        'original_source_sha256': digest(source), 'source_claim_environments': len(environments),
        'claim_status_counts': dict(sorted(Counter(r['status'] for r in ledger['claims']).items())),
        'finite_exhaustive_checks': finite, 'boundary_witnesses': witnesses,
        'mathematical_mutation_controls_rejected': mutations,
        'signed_phase_operator_derivation_checks': phase,
        'native_symbolic_replays': len(replays), 'native_input_sha256': native['input_sha256'],
        'native_packet_sha256': digest((HERE / 'R15_NATIVE_CERTIFICATE.json').read_bytes()),
        'native_negative_controls': native['negative_controls'], 'native_regular_basis': native['regular_basis'],
        'ugd1_upstream_replay': {'status': 'BYTE_AND_EXPECTED_HASH_MATCH', 'sha256': digest(ugd_bytes),
            'source_commit': pins['publications']['commit'], 'upstream_files_rewritten': False,
            'projection_blindness_witness': {'shared_projection': 7, 'seam_totals': [0, 2], 'full_states_differ': True}},
        'preserved_previous_file_count': len(pins['preserved_previous_files']),
        'prior_R1_R14_evidence_changed': False,
        'r14_global_foundational_interpretation': 'SUPERSEDED_BY_R14_SCOPE_CORRECTION',
        'r14_full_UGD_nonuniqueness_inference': 'WITHDRAWN; restricted examples remain restricted',
        'python': sys.version.split()[0],
        'node': subprocess.run(['node', '--version'], check=True, text=True, capture_output=True).stdout.strip(),
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ['status', 'source_claim_environments',
        'native_symbolic_replays', 'preserved_previous_file_count']}, indent=2))
    print('finite cases:', finite['finite_systems_and_observations'],
          '| signed phase basis checks:', phase['signed_basis_vectors'],
          '| UGD-1:', report['ugd1_upstream_replay']['status'])


if __name__ == '__main__':
    main()
