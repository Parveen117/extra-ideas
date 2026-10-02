#!/usr/bin/env python3
"""Check the current register and replay frozen native applications without changing evidence.

Scope: R1--R45 recorded evidence integrity, the complete R16 foundation
checker, and R17--R45 native application/claim-graph replays. Historical
recursive verifier entry points also pin navigation at their own stage;
this runner uses their unchanged native command construction blocks instead.
Written proof bindings are checked, not mechanically type-checked.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def load(path):
    return json.loads(Path(path).read_bytes())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local_path(path):
    result = (ROOT / path).resolve()
    require(result.is_relative_to(ROOT), 'Path outside repository: ' + path)
    return result


def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def check_anchor(row):
    data = local_path(row['path']).read_bytes()
    require(blob(data) == row['git_blob_sha'], 'Changed registered file: ' + row['path'])


def check_file(home, row):
    path = (home / row['path']).resolve()
    require(path.is_relative_to(home), 'Source path outside checkout')
    data = path.read_bytes()
    require(digest(data) == row['sha256'], 'Source SHA-256 mismatch: ' + row['path'])
    require(blob(data) == row['git_blob_sha'], 'Source Git blob mismatch: ' + row['path'])


def verify_register(book):
    rows = book['developments']
    require([row['id'] for row in rows] == [f'R{n}' for n in range(1, 46)],
            'Current register must cover R1--R45 in order')
    actual = {int(re.search(r'R(\d+)_VERIFICATION', p.name)[1])
              for p in HERE.glob('R*_VERIFICATION.json')}
    require(actual | {16} == set(range(1, 46)), 'Unregistered development report')
    for row in rows:
        for key in ['proof', 'verification', 'source_pins', 'native_certificate', 'ledger', 'scope_correction']:
            if key in row:
                check_anchor(row[key])
        record = load(local_path(row['verification']['path']))
        require(record['status'] == row['verification']['status'],
                'Changed recorded status: ' + row['id'])
        if row['id'] == 'R15':
            require(record['status'] == 'SCOPED_CHECKS_PASS_ORIGINAL_NOT_CERTIFIED',
                    'R15 original manuscript promoted')
            require(record['whole_original_manuscript_certified'] is False,
                    'R15 original manuscript promoted')
        else:
            require(record['status'].startswith('PASS'), 'Failed record: ' + row['id'])
        if 'canonical_report_sha256' in record:
            value = dict(record)
            expected = value.pop('canonical_report_sha256')
            require(digest(canonical(value)) == expected, 'Report digest mismatch: ' + row['id'])
    require(book['latest_development'] == 'R45', 'Changed latest development')
    require(book['formal_proof_assistant_verified'] is False and
            book['physical_constants_selected'] is False, 'Physical/formal claim promotion')
    check_anchor(book['historical_manifest'])
    check_anchor(book['integrated_historical_note'])
    return len(rows)


def verify_current_inputs(book, args):
    pins = load(HERE / 'R45_SOURCE_PINS.json')
    record = load(HERE / 'R45_VERIFICATION.json')
    require(digest((HERE / 'R45_SOURCE_PINS.json').read_bytes()) == record['source_pins_sha256'],
            'R45 source pins changed')
    for row in pins['local_inputs'] + pins['preserved_evidence']:
        check_file(ROOT, row)
    for key, home in [('rkf', args.rkf_root),
                      ('publications_comparison_only', args.publications_root)]:
        for row in pins[key]['files']:
            check_file(home, row)
    require(book['latest_research_commit'] == book['research_snapshot']['commit'],
            'Research snapshot mismatch')
    return {'local_files': len(pins['local_inputs']) + len(pins['preserved_evidence']),
            'preserved_prior_files': len(pins['preserved_evidence']),
            'upstream_files': len(pins['rkf']['files']) +
                              len(pins['publications_comparison_only']['files'])}


def verify_history(book, skip):
    if skip:
        require(not __import__('os').environ.get('GITHUB_ACTIONS'),
                'CI may not skip branch history')
        return {'status': 'NOT_RUN_LOCAL_SNAPSHOT'}
    commits = [book['latest_research_commit']] + [r['commit'] for r in book['integrated_branches']]
    for commit in commits:
        result = subprocess.run(['git', 'merge-base', '--is-ancestor', commit, 'HEAD'], cwd=ROOT)
        require(result.returncode == 0, 'Branch/source not integrated: ' + commit)
    return {'status': 'PASS_ANCESTRY', 'commits_checked': len(commits)}


def run(command):
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    require(result.returncode == 0,
            'Native replay failed: ' + ' '.join(command[:2]) + '\n' + result.stderr[-8000:])
    return result


def foundation(args, output):
    packet = ROOT / '02-relational-response/emk-topology-foundation/certificate'
    report_path, native_path = output / 'r16.json', output / 'r16-native.json'
    run([sys.executable, '-B', str(packet / 'verify.py'), '--rkf-root', str(args.rkf_root),
         '--output', str(report_path), '--native-output', str(native_path)])
    fresh, frozen = load(report_path), load(packet / 'VERIFICATION.json')
    for key in ['status', 'written_core_statements', 'whole_original_manuscript_certified',
                'source_pins_sha256', 'finite_exact_checks', 'native_replay', 'analytic_replays']:
        require(fresh[key] == frozen[key], 'R16 replay mismatch: ' + key)
    require(digest(native_path.read_bytes()) == frozen['native_replay']['certificate_sha256'],
            'R16 native certificate mismatch')
    return {'id': 'R16', 'status': fresh['status'], 'written_results': fresh['written_core_statements'],
            'native_symbolic_replays': fresh['native_replay']['proofs']}


def native_command(number, pins, args):
    path = HERE / f'verify_r{number}.py'
    tree = ast.parse(path.read_text())
    candidates = [node for node in ast.walk(tree) if isinstance(node, ast.Assign) and
                  any(isinstance(target, ast.Name) and target.id == 'command' for target in node.targets)]
    require(len(candidates) == 1, 'Ambiguous frozen native command: ' + path.name)
    assignment = candidates[0]
    bodies = [node.body for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
              and assignment in node.body]
    require(len(bodies) == 1, 'Frozen command outside a unique function body')
    body = bodies[0]
    block = [assignment]
    following = body[body.index(assignment) + 1]
    if isinstance(following, ast.For):
        # R27 appends its two certificate dependencies in a literal finite loop.
        require(isinstance(following.target, ast.Name) and isinstance(following.iter, ast.List)
                and all(isinstance(item, ast.Constant) and isinstance(item.value, int)
                        for item in following.iter.elts)
                and not following.orelse and following.body
                and all(isinstance(item, ast.AugAssign) and isinstance(item.target, ast.Name)
                        and item.target.id == 'command' and isinstance(item.op, ast.Add)
                        for item in following.body), 'Unexpected frozen command construction')
        block.append(following)
    # This construction is read from the hash-checked verifier, not a new command recipe.
    scope = {'__builtins__': {'str': str}, 'HERE': HERE, 'args': args, 'pins': pins}
    exec(compile(ast.Module(body=block, type_ignores=[]), str(path), 'exec'), scope)
    command = scope['command']
    require(isinstance(command, list) and command[0] == 'node', 'Non-native command')
    require(Path(command[1]).resolve().is_relative_to(ROOT), 'Application outside repository')
    return command


def proof_and_graph(number, row, pins, frozen):
    ledger = load(local_path(row['ledger']['path']))
    proof = local_path(row['proof']['path']).read_text()
    pattern = rf'### R{number}\.(\d+) — ([^\n]+)\n(.*?)(?=\n#{{2,3}} |\Z)'
    sections = list(re.finditer(pattern, proof, re.S))
    require(len(sections) == len(ledger['claims']), 'Written claim count mismatch')
    for index, (section, claim) in enumerate(zip(sections, ledger['claims']), 1):
        require(claim['id'] == f'R{number}.{index}' and int(section[1]) == index,
                'Written claim order mismatch')
        require(section[2] == claim['title'] and '**Proof.**' in section[3],
                'Written claim binding mismatch')
        require(digest(section[0].encode()) == claim['written_section_sha256'],
                'Written proof section changed')
    spec = importlib.util.spec_from_file_location(f'ci_frozen_r{number}', HERE / f'verify_r{number}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mutation_count = 0
    if hasattr(module, 'derivation_gate'):
        graph = module.derivation_gate(ledger, pins)
        require(graph == frozen['derivation_graph'], 'Derivation graph changed')
        mutation_fn = getattr(module, 'mutation_checks', None) or getattr(module, 'mutated_graph_checks', None)
        require(callable(mutation_fn), 'Frozen graph mutation checker missing')
        mutations = mutation_fn(ledger, pins)
        require(mutations == frozen['derivation_graph_mutations_rejected'] and all(mutations.values()),
                'Graph mutation rejection changed')
        mutation_count = len(mutations)
    return len(sections), mutation_count


def application(number, row, args):
    pins = load(HERE / f'R{number}_SOURCE_PINS.json')
    frozen = load(HERE / f'R{number}_VERIFICATION.json')
    require(digest((HERE / f'R{number}_SOURCE_PINS.json').read_bytes()) == frozen['source_pins_sha256'],
            'Frozen pins changed')
    started = time.monotonic()
    command = native_command(number, pins, args)
    native = json.loads(run(command).stdout)
    native_bytes = (json.dumps(native, indent=2, ensure_ascii=False) + '\n').encode()
    require(digest(native_bytes) == frozen['native_certificate_sha256'],
            f'R{number} native certificate replay mismatch')
    require(native['exact_checks'] and all(c['passed'] for c in native['exact_checks']),
            'Failed native exact group')
    require(all(c['replay'] == 'REPLAY_MATCH' for c in native['symbolic_replays']),
            'Failed native word replay')
    require(all(native['rejected_false_alternatives'].values()) and
            native['altered_native_certificate_rejected'], 'Failed native negative control')
    if 'counters' in frozen:
        counters = {c['name']: {k: v for k, v in c.items() if k not in ('name', 'passed')}
                    for c in native['exact_checks']}
        require(counters == frozen['counters'], 'Frozen native counters changed')
    for option in ['--expected-input-sha256'] + [x for x in command if x.startswith('--expected-r')]:
        wrong = list(command)
        wrong[wrong.index(option) + 1] = '0' * 64
        result = subprocess.run(wrong, capture_output=True, text=True)
        require(result.returncode != 0 and 'pin mismatch' in result.stderr,
                'Wrong native input/parent pin accepted')
    proofs, mutations = proof_and_graph(number, row, pins, frozen)
    return {'id': f'R{number}', 'status': native['status'], 'written_results': proofs,
            'exact_check_groups': native['check_count'],
            'native_symbolic_replays': native['symbolic_replay_count'],
            'graph_mutations_rejected': mutations,
            'elapsed_seconds': round(time.monotonic() - started, 3)}


def main():
    require(__debug__, 'Run without -O; frozen assertions must remain enabled')
    require(sys.version_info[:2] in [(3, 11), (3, 12)], 'Use Python 3.11 or 3.12')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rkf-root', required=True, type=Path)
    parser.add_argument('--publications-root', required=True, type=Path)
    parser.add_argument('--skip-history', action='store_true', help='Local connector snapshot only; forbidden in CI')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    args.rkf_root, args.publications_root = args.rkf_root.resolve(), args.publications_root.resolve()
    if args.output:
        require(not args.output.resolve().is_relative_to(ROOT), 'Write CI output outside frozen evidence tree')
    book = load(ROOT / 'CURRENT_MANIFEST.json')
    records = verify_register(book)
    integrity = verify_current_inputs(book, args)
    history = verify_history(book, args.skip_history)
    print(f'PASS: {records} recorded stages; {integrity["local_files"]} local and {integrity["upstream_files"]} upstream pins.', flush=True)
    replays = []
    with tempfile.TemporaryDirectory(prefix='extra-ideas-ci-') as tmp:
        replays.append(foundation(args, Path(tmp)))
        print('PASS: R16 complete scoped foundation checker.', flush=True)
        for number in range(17, 46):
            replays.append(application(number, book['developments'][number - 1], args))
            print(f'PASS: R{number} native application, proof bindings and available graph controls.', flush=True)
    # Recheck source identity after runtime execution; applications must not refresh old evidence.
    verify_current_inputs(book, args)
    verify_register(book)
    result = {'schema': 'extra-ideas.current-ci.v1', 'status': 'PASS_CURRENT_REPOSITORY_CI',
              'recorded_stages': records, 'integrity': integrity, 'history': history,
              'runtime_replayed_stages': [row['id'] for row in replays], 'replays': replays,
              'native_symbolic_replays': sum(row['native_symbolic_replays'] for row in replays),
              'historical_recursive_entry_points_rerun': False,
              'earlier_R1_R15_runtime_reexecuted': False,
              'formal_proof_assistant_verified': False, 'physical_validation_claimed': False}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)
    print(f'30 stage replays; {result["native_symbolic_replays"]} native symbolic replays.', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
