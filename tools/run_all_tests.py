"""Run every tool's tests (stdlib only; Python 3.12).  Usage: python run_all_tests.py"""
import pathlib
import subprocess
import sys

root = pathlib.Path(__file__).resolve().parent
total, failed = 0, []
for folder in sorted(p for p in root.iterdir() if p.is_dir() and list(p.glob('test_*.py'))):
    out = subprocess.run([sys.executable, '-B', '-m', 'unittest', 'discover', '-p', 'test_*.py'],
                         cwd=folder, capture_output=True, text=True)
    line = [l for l in out.stderr.splitlines() if l.startswith('Ran ')]
    n = int(line[0].split()[1]) if line else 0
    total += n
    ok = out.returncode == 0
    if not ok:
        failed.append(folder.name)
    print(f'{folder.name:5s} {n:3d} tests  {"OK" if ok else "FAILED"}')
print(f'total {total} tests in {sum(1 for _ in root.glob("*/test_*.py"))} files; failed: {failed or "none"}')
sys.exit(1 if failed else 0)
