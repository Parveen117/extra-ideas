"""Read-only, source-bound replay of CO1's declared symbolic model."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_files(home, files):
    for relative, expected in files.items():
        path = (home / relative).resolve()
        require(path.is_relative_to(home), 'Source outside checkout')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                'Source drift: ' + relative)


def verify(publications_root):
    require(sys.version_info[:2] == (3, 12), 'Use Python 3.12')
    require(__debug__, 'Optimization disables upstream assertions')
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    check_files(ROOT, pins['local'])
    check_files(publications_root, pins['publications']['files'])
    spec = importlib.util.spec_from_file_location('co1_replay', HERE / 'co1_expansion_as_fall.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.run()
    require(result == json.loads((HERE / 'CO1_RESULT.json').read_text()),
            'Replayed CO1 result differs from frozen evidence')
    tests = subprocess.run(
        [sys.executable, '-B', '-m', 'unittest', 'discover', '-s', str(HERE),
         '-p', 'test_co1.py', '-v'], capture_output=True, text=True)
    require(tests.returncode == 0, 'CO1 tests failed:\n' + tests.stdout + tests.stderr)
    require('Ran 4 tests' in tests.stderr and 'skipped=' not in tests.stderr,
            'Unexpected CO1 test scope')
    check_files(ROOT, pins['local'])
    check_files(publications_root, pins['publications']['files'])
    print(json.dumps({'status': 'PASS_CO1_SOURCE_BOUND_SYMBOLIC_REPLAY',
                      'tests_passed': 4, 'statement_groups': ['T1', 'T2', 'T3', 'T4', 'T5'],
                      'local_source_pins': len(pins['local']),
                      'publications_commit': pins['publications']['commit'],
                      'result_sha256': hashlib.sha256((HERE / 'CO1_RESULT.json').read_bytes()).hexdigest(),
                      'scope': pins['scope']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', required=True, type=Path)
    args = parser.parse_args()
    verify(args.publications_root.resolve())
