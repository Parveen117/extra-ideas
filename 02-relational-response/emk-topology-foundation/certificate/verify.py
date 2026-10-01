#!/usr/bin/env python3
"""Reproduce the scoped EMK constructive certificate against pinned RKF sources.

This verifies identities, finite witnesses, source identity and ledger coverage.
It does not type-check written proofs or certify the original universal claim.
"""
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


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def check_bytes(data, record):
    if sha256(data) != record['sha256']:
        raise ValueError('SHA-256 mismatch: '+record['path'])
    if 'git_blob_sha' in record:
        actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if actual != record['git_blob_sha']:
            raise ValueError('Git blob mismatch: '+record['path'])


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: certificate assertions must remain enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=HERE/'VERIFICATION.json')
    parser.add_argument('--native-output', type=Path, default=HERE/'NATIVE_REPLAY.json')
    args = parser.parse_args()
    pins_bytes = (HERE/'SOURCE_PINS.json').read_bytes()
    pins = json.loads(pins_bytes)
    for row in pins['local_files']:
        check_bytes((ROOT/row['path']).read_bytes(), row)
    for row in pins['rkf']['files']:
        check_bytes((args.rkf_root/row['path']).read_bytes(), row)

    source = (ROOT/'source/emk_topology.original.tex').read_bytes()
    ledger = json.loads((HERE/'ORIGINAL_CLAIM_LEDGER.json').read_text())
    assert sha256(source) == ledger['original_source_sha256']
    text = source.decode().replace('\r\n','\n').replace('\r','\n')
    pattern = r'\\begin\{(theorem|lemma|proposition|corollary)\}(?:\[([^\]]*)\])?(.*?)\\end\{\1\}'
    original = list(re.finditer(pattern, text, re.S))
    assert len(original) == len(ledger['claims']) == 48
    for index, (match, row) in enumerate(zip(original, ledger['claims']), 1):
        assert row['id'] == f'ET{index:02}'
        assert row['environment'] == match[1] and row['title'] == (match[2] or '')
        assert row['start_line'] == text[:match.start()].count('\n')+1
        assert row['end_line'] == text[:match.end()].count('\n')+1
        assert row['normalized_environment_sha256'] == sha256(match[0].encode())
        assert row['replacement_scope'] and not row['original_promoted_to_unconditional_theorem']
    assert ledger['whole_source_verdict'] == 'NOT_CERTIFIED_AS_WRITTEN'

    core = json.loads((HERE/'CORE_LEDGER.json').read_text())
    available = {x['id'] for x in core['premises']}
    for index, row in enumerate(core['theorems'], 1):
        assert row['id'] == f'C{index}'
        assert set(row['depends_on']) <= available
        assert row['status'] == 'WRITTEN_PROOF_IN_CORRECTED_CONSTRUCTION'
        assert not row['formal_proof_assistant_verified']
        available.add(row['id'])
    assert len(core['theorems']) == 14
    for row in core['analytic_citations']:
        assert set(row['input_theorems']) <= available
        available.add(row['id'])
    for row in ledger['claims']:
        assert set(row['replacement_ids']) <= available
    manuscript = (ROOT/'emk_topology_foundation.tex').read_text()
    found = re.findall(r'\\begin\{(?:theorem|proposition)\}\[(C\d+):', manuscript)
    assert found == [f'C{i}' for i in range(1,15)]
    assert manuscript.count(r'\begin{proof}') == manuscript.count(r'\end{proof}') == 14

    # Deliberately import only after both checker and verifier hashes pass.
    import checks
    finite = checks.build()
    native = json.loads(subprocess.run([
        'node', str(HERE/'native_bridge.cjs'), '--rkf-root', str(args.rkf_root),
        '--expected-input-sha256', pins['native_input_sha256'],
    ], check=True, capture_output=True, text=True).stdout)
    assert native['status'] == 'PASS_SCOPED_NATIVE_BRIDGE'
    replays = [r for group in native['groups'].values() for r in group['results']]
    assert len(replays) == 13 and all(r['replay']=='REPLAY_MATCH' for r in replays)
    assert native['full_algebra_dimension'] == 4
    assert native['all_coordinate_turn_components_zero']
    assert native['complete_cut_corner_basis_cases'] == 16
    assert native['negative_controls']['altered_rotation_square_rejected']

    analytic = []
    with tempfile.TemporaryDirectory(prefix='emk-native-analytic-') as tmp:
        for item in pins['analytic_replays']:
            output = Path(tmp)/(item['id']+'.json')
            subprocess.run([sys.executable, str(args.rkf_root/item['verifier']),
                            '--output', str(output)], check=True, capture_output=True, text=True)
            result = json.loads(output.read_text())
            expected_hash = result.pop('canonical_result_sha256')
            assert sha256(canonical(result)) == expected_hash == item['expected_canonical_result_sha256']
            assert result['status'] == item['passing_status']
            assert all(x['passed'] for x in result['checks'].values())
            assert not result['missing_obligations']
            analytic.append({'id':item['id'],'status':result['status'],
                             'canonical_result_sha256':expected_hash,
                             'checks':len(result['checks']),
                             'negative_controls':len(result['negative_controls']),
                             'scope':result['claim_boundary']})

    rejected_pin = False
    test_pin = pins['local_files'][0]
    try:
        check_bytes((ROOT/test_pin['path']).read_bytes()+b'\nMUTATED\n', test_pin)
    except ValueError:
        rejected_pin = True
    assert rejected_pin
    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False)+'\n').encode()
    report = {
        'schema':'emk-constructive-cut.verification.v1',
        'status':'PASS_SCOPED_CONSTRUCTIVE_CERTIFICATE',
        'written_core_statements':14,
        'written_proof_status':'Every new core statement has a proof; ledger coverage checked, proof text not mechanically type-checked.',
        'whole_original_manuscript_certified':False,
        'original_claims_accounted_for':48,
        'formal_proof_assistant_verified':False,
        'external_peer_review_claimed':False,
        'extra_work_integration_performed':False,
        'source_pins_sha256':sha256(pins_bytes),
        'local_pinned_files':len(pins['local_files']),
        'canonical_rkf_files_checked':len(pins['rkf']['files']),
        'rkf_commit':pins['rkf']['commit'],
        'finite_exact_checks':finite,
        'native_replay':{'proofs':len(replays),'complete_cut_corner_basis_cases':16,
                         'certificate_sha256':sha256(native_bytes),
                         'all_coordinate_turn_components_zero':True,
                         'tampered_certificate_rejected':True},
        'analytic_replays':analytic,
        'source_mutation_rejected':rejected_pin,
        'runtime':{'python':sys.version.split()[0],
                   'node':subprocess.run(['node','--version'],check=True,capture_output=True,text=True).stdout.strip()},
        'mathematical_scope':'Free typed records and explicitly constructed role, coefficient, recognition and completion objects; not uniqueness from the original prose or physical identification.',
    }
    report['canonical_report_sha256'] = sha256(canonical(report))
    args.native_output.parent.mkdir(parents=True,exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(report['status'])
    print('14 written statements; 48 original dispositions; 13 native replays; both native analytic audit hashes match.')
    print(report['canonical_report_sha256'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
