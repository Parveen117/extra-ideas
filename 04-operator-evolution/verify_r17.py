#!/usr/bin/env python3
"""Verify source-derived cut-history census, memory and record-erasure results."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FOUNDATION=ROOT/'02-relational-response/emk-topology-foundation'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()


def check_file(home,record):
    data=(home/record['path']).read_bytes()
    if digest(data)!=record['sha256']:
        raise ValueError('Source hash mismatch: '+record['path'])
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if actual!=record['git_blob_sha']:
        raise ValueError('Git blob mismatch: '+record['path'])


def main():
    if not __debug__:
        raise RuntimeError('Certificate assertions require Python without -O.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rkf-root',required=True,type=Path)
    ap.add_argument('--output',type=Path,default=HERE/'R17_VERIFICATION.json')
    ap.add_argument('--native-output',type=Path,default=HERE/'R17_NATIVE_CERTIFICATE.json')
    args=ap.parse_args()
    pin_bytes=(HERE/'R17_SOURCE_PINS.json').read_bytes()
    pins=json.loads(pin_bytes)
    for record in pins['local_inputs']+pins['preserved_evidence']:
        check_file(ROOT,record)
    for record in pins['rkf']['files']:
        check_file(args.rkf_root,record)

    # The R16 packet is a frozen source dependency. Reproduce it into temporary
    # outputs, never overwrite its old verification or source-pin files.
    with tempfile.TemporaryDirectory(prefix='r17-foundation-') as tmp:
        report_path=Path(tmp)/'foundation.json'
        replay_path=Path(tmp)/'foundation-native.json'
        subprocess.run([sys.executable,'-B',str(FOUNDATION/'certificate/verify.py'),
                        '--rkf-root',str(args.rkf_root),'--output',str(report_path),
                        '--native-output',str(replay_path)],check=True,capture_output=True,text=True)
        foundation=json.loads(report_path.read_text())
    assert foundation['status']=='PASS_SCOPED_CONSTRUCTIVE_CERTIFICATE'
    assert not foundation['whole_original_manuscript_certified']
    assert not foundation['formal_proof_assistant_verified']
    source_roles=foundation['finite_exact_checks']['signed_roles']
    native_input=json.loads((HERE/'R17_NATIVE_INPUT.json').read_text())
    assert native_input['constructed_maps']==source_roles['constructed_maps']
    assert native_input['ordered_roles']==source_roles['primitive_data']['ordered_roles']
    assert native_input['source_foundation']['pins_sha256']==digest((FOUNDATION/'certificate/SOURCE_PINS.json').read_bytes())
    assert native_input['source_foundation']['commit']==pins['foundation_commit']

    text=(ROOT/'02-relational-response/CUT_HISTORY_MEMORY_R17.md').read_text()
    ledger=json.loads((HERE/'R17_CLAIM_LEDGER.json').read_text())
    pattern=r'### R17\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{2,3} |\Z)'
    results=list(re.finditer(pattern,text,re.S))
    assert len(results)==len(ledger['claims'])==10
    for n,(match,claim) in enumerate(zip(results,ledger['claims']),1):
        assert claim['id']==f'R17.{n}' and int(match[1])==n
        assert claim['title']==match[2]
        assert claim['written_section_sha256']==digest(match[0].encode())
        assert '**Proof.**' in match[3]
        assert claim['evidence']=='WRITTEN_PROOF_WITH_SCOPED_EXACT_CERTIFICATE'
    assert ledger['native_source_count_is_empirical_detector_law'] is False

    command=['node',str(HERE/'cut_history_dynamics.cjs'),'--rkf-root',str(args.rkf_root),
             '--input',str(HERE/'R17_NATIVE_INPUT.json'),
             '--expected-input-sha256',pins['native_input_sha256']]
    native=json.loads(subprocess.run(command,check=True,capture_output=True,text=True).stdout)
    assert native['status']=='PASS_R17_NATIVE_CUT_HISTORY'
    assert native['check_count']==12 and all(x['passed'] for x in native['exact_checks'])
    assert native['symbolic_replay_count']==15
    assert all(x['replay']=='REPLAY_MATCH' for x in native['symbolic_replays'])
    assert len(native['rejected_false_alternatives'])==11
    assert all(native['rejected_false_alternatives'].values())
    assert native['altered_native_certificate_rejected']

    # The native call must reject an independently supplied wrong input pin.
    wrong=list(command);wrong[-1]='0'*64
    rejected=subprocess.run(wrong,capture_output=True,text=True)
    assert rejected.returncode!=0 and 'input pin mismatch' in rejected.stderr
    native_bytes=(json.dumps(native,indent=2,ensure_ascii=False)+'\n').encode()
    counters={row['name']:{k:v for k,v in row.items() if k not in ('name','passed')}
              for row in native['exact_checks']}
    report={
      'schema':'extra-ideas.r17.verification.v1','date':'2026-10-01',
      'status':'PASS_R17_SOURCE_COUNT_MEMORY',
      'new_written_theorems':10,'exact_check_groups':native['check_count'],
      'native_symbolic_replays':native['symbolic_replay_count'],
      'counters':counters,'rejected_false_alternatives':native['rejected_false_alternatives'],
      'altered_native_certificate_rejected':True,'wrong_native_input_pin_rejected':True,
      'foundation_reproduced':{'commit':pins['foundation_commit'],'status':foundation['status'],
        'written_statements':foundation['written_core_statements'],
        'native_replays':foundation['native_replay']['proofs'],
        'analytic_replays':foundation['analytic_replays'],
        'original_universal_claim_certified':False},
      'upstream_canonical_commit':pins['rkf']['commit'],
      'upstream_files_verified':len(pins['rkf']['files']),
      'preserved_prior_files':len(pins['preserved_evidence']),
      'source_pins_sha256':digest(pin_bytes),
      'native_certificate_sha256':digest(native_bytes),
      'foundation_connected_to_extra_ideas':True,
      'new_formal_proof_assistant_verified':False,
      'physical_frequency_selection_proved':False,
      'core_ordinary_complex_or_gaussian_input':False,
      'runtime':{'python':sys.version.split()[0],
         'node':subprocess.run(['node','--version'],check=True,capture_output=True,text=True).stdout.strip()},
      'scope':'Free typed cut-history census, derived signed role actions and seam aperture. Exact identities and written all-depth proofs; not a universal physical noise law or complete UGD-state replacement.'}
    report['canonical_report_sha256']=digest(canonical(report))
    args.native_output.parent.mkdir(parents=True,exist_ok=True)
    args.native_output.write_bytes(native_bytes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(report['status'])
    print('10 written results; 12 exact check groups; 15 symbolic replays; 11 false alternatives and 2 altered contracts rejected.')
    print(report['canonical_report_sha256'])
    return 0


if __name__=='__main__':
    raise SystemExit(main())
