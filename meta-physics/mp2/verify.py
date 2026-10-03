"""Replay MP-2 against byte-pinned Publications sources, Python 3.12 only."""
import argparse
import ast
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest
import evolving_response as e

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(b):
    return hashlib.sha256(b).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def json_value(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {k: json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    return value


def source_checks(publications):
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    homes = {'extra-ideas': ROOT, 'Publications': publications}
    for row in pins['consumed_unchanged_files']:
        home = homes[row['repository']]
        path = (home/row['path']).resolve()
        require(path.is_relative_to(home), 'Pinned path escapes checkout')
        b = path.read_bytes()
        require(sha(b) == row['sha256'], 'Source changed: '+str(path))
        gitsha = hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
        require(gitsha == row['git_blob_sha1'], 'Source Git digest changed')
    base = 'papers/yang-mills-certified-benchmark/certificates/'
    todo, seen = ['ym54_response_protocol'], set()
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        seen.add(name)
        syntax = ast.parse((publications/base/(name+'.py')).read_text())
        for node in ast.walk(syntax):
            names = ([x.name for x in node.names] if isinstance(node, ast.Import) else
                     [node.module] if isinstance(node, ast.ImportFrom) else [])
            for child in names:
                if child and (publications/base/(child+'.py')).is_file():
                    todo.append(child)
    require(sorted(seen) == pins['runtime_modules'], 'Runtime import closure changed')
    pinned = {(r['repository'], r['path']) for r in pins['consumed_unchanged_files']}
    require(all(('Publications', base+n+'.py') in pinned for n in seen), 'Unpinned runtime source')
    require(pins['idea_sources'] == ['033', '053', '093'], 'Changed idea-development scope')
    require(pins['raw_private_pdfs_published'] is False, 'Raw PDF promotion is not authorized here')
    require(not list(HERE.parent.rglob('*.pdf')), 'Raw source PDFs are not this public development')
    # MP-1's own frozen replay also checks the complete R1--R46 register.
    replay = subprocess.run([sys.executable, '-B', str(HERE.parent/'mp1/verify.py'), '--check'],
                            cwd=ROOT, capture_output=True, text=True, timeout=120)
    require(replay.returncode == 0, replay.stdout+replay.stderr)
    old = json.loads(replay.stdout)
    require(old['tests'] == 23, 'Prior test scope changed')
    require(old['certificate_sha256'] ==
            'a3198a0bfbbcf9be5e0a69667d917b0d6cc3784485ec8d716a0944f5890aaf7b',
            'The frozen MP-1 certificate changed')
    return dict(pinned_files=len(pins['consumed_unchanged_files']), runtime_modules=len(seen),
                frozen_r_stages=46, mp1_tests_replayed=23,
                mp1_certificate_sha256=old['certificate_sha256'])


def build_report(publications):
    require(sys.version_info[:2] == (3, 12), 'Python 3.12 only')
    integrity = source_checks(publications)
    e.configure(publications)
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(HERE), pattern='test_*.py')
    results = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    require(results.wasSuccessful(), stream.getvalue())
    require(results.testsRun == 22 and not results.skipped, 'Unexpected MP-2 test scope')
    proof = (HERE/'THEOREM.md').read_text()
    require(all(f'## T{i} ' in proof for i in range(1, 8)), 'Missing written theorem')
    inputs = ['meta-physics/mp2/'+name for name in
              ['THEOREM.md', 'evolving_response.py', 'test_evolving_response.py',
               'SOURCE_PINS.json', 'verify.py']]+['.github/workflows/meta-physics-mp2.yml']
    report = dict(
        schema='MP2-VERIFICATION-1', status='PASS_SCOPED_WRITTEN_PROOFS_AND_EXACT_CHECKS',
        python='3.12', written_results=7, tests_passed=results.testsRun,
        prior_mp1_tests_passed=23, cofactor_and_flow_grid_matrices=624,
        universal_claim_basis='Seven written proofs; finite tests do not prove the continuum claims',
        continuous_time_implemented_as_numeric_ode=False,
        chronological_contraction_test='Exact ordered resolvent controls; heat limit is proved in T4',
        native_engines_copied=0, integrity=integrity, examples=json_value(e.exact_examples()),
        idea_sources_developed_in_scope=['033', '053', '093'], complete_pdfs_certified=0,
        formal_proof_assistant_verified=False, independently_peer_reviewed=False,
        physical_law_or_empirical_effect_certified=False, clay_mass_gap_proved=False,
        selections=['radial-relaxation projection at fixed oriented area', 'initial response source',
                    'optional frame rotation', 'four equally counted compact turns', 'heat clock'],
        open_gates=['physical selection and units', 'state-dependent protocol',
                    'retarded memory/passivity realization (source 051)',
                    'moving interacting ground source', 'four-dimensional continuum Yang--Mills'],
        input_sha256={path: sha((ROOT/path).read_bytes()) for path in inputs})
    report['canonical_report_sha256'] = sha(canonical(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', required=True, type=Path)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    report = build_report(args.publications_root.resolve())
    data = json.dumps(report, indent=2, ensure_ascii=False)+'\n'
    path = HERE/'VERIFICATION.json'
    if args.write:
        path.write_text(data)
    else:
        require(path.read_text() == data, 'Certificate differs from replay; review changes')
    print(json.dumps(dict(status=report['status'], mp2_tests=report['tests_passed'],
                          prior_mp1_tests=23, frozen_r_stages=46,
                          certificate_sha256=sha(data.encode()))))


if __name__ == '__main__':
    main()
