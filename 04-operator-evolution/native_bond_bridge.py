"""R4 bridge to the existing RKF native return; no copied return solver.

The finite matrices are the stated faithful real calibration. General native
identities and their hypotheses are proved in NATIVE_BOND_BRIDGE_R4.md.
"""

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess

from aghora_return import add, identity, matrix, mul, scalar, scale, sub


NATIVE_R = matrix(((0, -1), (1, 0)))
NATIVE_K = matrix(((0, 1), (1, 0)))
NATIVE_L = mul(NATIVE_K, NATIVE_R)


def response_bridge(value):
    """B_x=(I+xL)/2 and its raw commutator; x is not normalized to one."""
    x = scalar(value)
    if x < 0:
        raise ValueError('The native response contract here is nonnegative')
    unit = identity(2)
    bond_candidate = scale(add(unit, scale(NATIVE_L, x)), Q(1, 2))
    commutator = sub(mul(NATIVE_K, bond_candidate), mul(bond_candidate, NATIVE_K))
    return {
        'x': x,
        'bond_candidate': bond_candidate,
        'commutator': commutator,
        'bond_defect': sub(mul(bond_candidate, bond_candidate), bond_candidate),
        'complex_defect': add(mul(commutator, commutator), unit),
        'exact_cut': x == 1,
    }


def closure_budget(lower, upper):
    """Transport a source interval to bond and raw complex-structure defects.

    A finite interval is an enclosure, not a proof of exact completed closure.
    The distance is in native coefficient mass (also the operator norm in the
    stated faithful calibration), not an undeclared physical norm.
    """
    lower, upper = scalar(lower), scalar(upper)
    if not 0 <= lower <= upper:
        raise ValueError('An ordered nonnegative source interval is required')
    return {
        'bond_defect_interval': ((lower**2 - 1)/4, (upper**2 - 1)/4),
        'complex_defect_interval': (1 - upper**2, 1 - lower**2),
        'distance_to_cut_bound': max(abs(lower - 1), abs(upper - 1))/2,
        'excludes_exact_cut': not lower <= 1 <= upper,
    }


def period_cut_condition(first, second):
    """Existing native SY1 criterion, used as a dependency, not a new theorem."""
    first, second = scalar(first), scalar(second)
    if first <= 0 or second <= 0:
        raise ValueError('The supplied two-cell profile must be positive')
    mismatch = first - second - 1
    return {'mismatch': mismatch, 'exact_completed_cut': mismatch == 0}


def check_runtime_pins(rkf_root, pins_path=None):
    here = Path(__file__).resolve().parent
    pins = json.loads(Path(pins_path or here / 'R4_SOURCE_PINS.json').read_text())
    runtime_pins = [entry for entry in pins['files'] if entry['role'] == 'runtime']
    if len(runtime_pins) != 4:
        raise ValueError('The R4 runtime contract requires four pinned source files')
    root, checked = Path(rkf_root).resolve(), {}
    for entry in runtime_pins:
        path = root / entry['path']
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise ValueError('Upstream source hash mismatch: ' + entry['path'])
        checked[entry['path']] = actual
    return checked


def load_native_probe(rkf_root):
    """Check source identity before executing its unchanged exact Node code."""
    checked = check_runtime_pins(rkf_root)
    probe = Path(__file__).with_name('native_bond_probe.cjs')
    result = subprocess.run(
        ['node', str(probe), str(Path(rkf_root).resolve())],
        check=True, capture_output=True, text=True, timeout=30)
    packet = json.loads(result.stdout)
    packet['checked_runtime_sha256'] = checked
    return packet
