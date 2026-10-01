#!/usr/bin/env python3
"""Verify native cut-address transport without overwriting prior certificates."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def check_file(root, record):
    data = (root / record['path']).read_bytes()
    if digest(data) != record['sha256']:
        raise ValueError('Source hash mismatch: ' + record['path'])
    git_hash = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    if git_hash != record['git_blob_sha']:
        raise ValueError('Git blob mismatch: ' + record['path'])


def main():
    if not __debug__:
        raise RuntimeError('Certificate assertions require Python without -O.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rkf-root', required=True, type=Path)
    ap.add_argument('--publications-root', required=True, type=Path)
    ap.add_argument('--output', type=Path, default=HERE / 'R18_VERIFICATION.json')
    ap.add_argument('--native-output', type=Path, default=HERE / 'R18_NATIVE_CERTIFICATE.json')
    args = ap.parse_args()
    pin_bytes = (HERE / 'R18_SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    for record in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, record)
    for key, home in [('rkf', args.rkf_root), ('publications_comparison_only', args.publications_root)]:
        for record in pins[key]['files']:
            check_file(home, record)

    # The entire constructive source chain is replayed into temporary files.
    # Hash agreement with the frozen R17 report binds this adapter to that chain.
    with tempfile.TemporaryDirectory(prefix='r18-source-') as tmp:
        report_path, native_path = Path(tmp) / 'r17.json', Path(tmp) / 'r17-native.json'
        subprocess.run([sys.executable, '-B', str(HERE / 'verify_r17.py'),
                        '--rkf-root', str(args.rkf_root), '--output', str(report_path),
                        '--native-output', str(native_path)], check=True, capture_output=True, text=True)
        previous = json.loads(report_path.read_text())
        frozen = json.loads((HERE / 'R17_VERIFICATION.json').read_text())
        # Runtime versions are informative, not mathematical dependencies.
        for key in ['status', 'source_pins_sha256', 'native_certificate_sha256',
                    'foundation_reproduced', 'counters']:
            assert previous[key] == frozen[key], 'R17 replay mismatch: ' + key
    assert previous['status'] == 'PASS_R17_SOURCE_COUNT_MEMORY'
    assert previous['new_formal_proof_assistant_verified'] is False
    assert previous['core_ordinary_complex_or_gaussian_input'] is False

    proof = (ROOT / '02-relational-response/CUT_TRANSPORT_METRIC_R18.md').read_text()
    ledger = json.loads((HERE / 'R18_CLAIM_LEDGER.json').read_text())
    pattern = r'### R18\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)'
    sections = list(re.finditer(pattern, proof, re.S))
    assert len(sections) == len(ledger['claims']) == 9
    for n, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        assert claim['id'] == f'R18.{n}' and int(section[1]) == n
        assert claim['title'] == section[2]
        assert claim['written_section_sha256'] == digest(section[0].encode())
        assert '**Proof.**' in section[3]
        assert claim['evidence'] == 'WRITTEN_PROOF_WITH_SCOPED_EXACT_CERTIFICATE'
    assert not ledger['physical_metric_selection_proved']
    assert not ledger['continuum_convergence_proved']
    assert not ledger['full_ugd_replaced_by_two_roles']
    assert ledger['address_and_aggregation_are_explicit_target_definitions']

    command = ['node', str(HERE / 'cut_transport.cjs'), '--rkf-root', str(args.rkf_root),
               '--input', str(HERE / 'R17_NATIVE_INPUT.json'),
               '--expected-input-sha256', pins['native_input_sha256']]
    native = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_R18_NATIVE_CUT_TRANSPORT'
    assert native['check_count'] == 10 and all(x['passed'] for x in native['exact_checks'])
    assert native['symbolic_replay_count'] == 8
    assert all(x['replay'] == 'REPLAY_MATCH' for x in native['symbolic_replays'])
    assert len(native['rejected_false_alternatives']) == 10
    assert all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']
    assert not native['external_complex_phase_input']
    assert not native['formal_proof_assistant_verified']

    wrong = list(command)
    wrong[-1] = '0' * 64
    result = subprocess.run(wrong, capture_output=True, text=True)
    assert result.returncode != 0 and 'input pin mismatch' in result.stderr
    # Explicitly exercise file-tamper rejection as well as checking current bytes.
    with tempfile.TemporaryDirectory(prefix='r18-tamper-') as tmp:
        root = Path(tmp)
        (root / 'altered.txt').write_text('altered source')
        rejected = False
        try:
            check_file(root, {'path': 'altered.txt', 'sha256': digest(b'original source'),
                              'git_blob_sha': '0' * 40})
        except ValueError:
            rejected = True
        assert rejected

    # Verify old sources again after executing the adapters, proving this run
    # did not refresh or rewrite the evidence it depended upon.
    for record in pins['preserved_evidence']:
        check_file(ROOT, record)
    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False) + '\n').encode()
    counters = {row['name']: {k: v for k, v in row.items() if k not in ('name', 'passed')}
                for row in native['exact_checks']}
    report = {
        'schema': 'extra-ideas.r18.verification.v1', 'date': '2026-10-01',
        'status': 'PASS_R18_CUT_ADDRESS_TRANSPORT',
        'written_results': 9, 'exact_check_groups': 10, 'native_symbolic_replays': 8,
        'counters': counters, 'rejected_false_alternatives': native['rejected_false_alternatives'],
        'altered_native_certificate_rejected': True, 'wrong_input_pin_rejected': True,
        'source_file_tamper_rejected': True,
        'previous_source_replayed': {'commit': pins['source_commit'], 'status': previous['status'],
            'native_certificate_sha256': previous['native_certificate_sha256'],
            'foundation': previous['foundation_reproduced']},
        'upstream_canonical_commit': pins['rkf']['commit'],
        'upstream_files_verified': len(pins['rkf']['files']),
        'comparison_source_commit': pins['publications_comparison_only']['commit'],
        'comparison_files_hash_verified': len(pins['publications_comparison_only']['files']),
        'comparison_physics_certificates_reexecuted': False,
        'preserved_prior_files': len(pins['preserved_evidence']),
        'source_pins_sha256': digest(pin_bytes), 'native_certificate_sha256': digest(native_bytes),
        'physical_metric_selection_proved': False, 'physical_c_or_alpha_derived': False,
        'continuum_convergence_proved': False, 'formal_proof_assistant_verified': False,
        'core_external_complex_or_wave_law_input': False,
        'address_and_aggregation_are_explicit_target_definitions': True,
        'known_hadamard_walk_lineage_acknowledged': True,
        'runtime': {'python': sys.version.split()[0], 'node': subprocess.run(
            ['node', '--version'], check=True, capture_output=True, text=True).stdout.strip()},
        'scope': 'Native cut-address construction; exact finite signed transport, wave stencil and formal jets. Physical protocol selection, large-scale convergence and rod/clock calibration remain open.'}
    report['canonical_report_sha256'] = digest(canonical(report))
    args.native_output.parent.mkdir(parents=True, exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(report['status'])
    print('9 written results; 10 exact check groups; 8 native replays; 10 false alternatives rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
