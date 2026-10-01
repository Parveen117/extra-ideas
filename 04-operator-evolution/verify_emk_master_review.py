"""Reproduce the EMK master review against unchanged public source bytes."""

import argparse
import hashlib
import importlib.util
import inspect
import io
import json
from pathlib import Path
import sys
import unittest
from fractions import Fraction as F


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(root, entries):
    for item in entries:
        path = root / item['path']
        data = path.read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if blob != item['git_blob_sha'] or sha256(path) != item['sha256']:
            raise ValueError('Changed upstream source: ' + item['path'])
    return len(entries)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def cross_checks(root):
    g3 = load('emkg3_master_review', root / 'papers/emk-recognition-geometry/certificates/emkg3_helical_sheet_memory.py')
    t1 = load('emkt1_master_review', root / 'papers/emk-ugd-algebra/certificates/emkt1_master_tensor_and_time.py')
    g2 = g3.g2
    p = (F(0), F(0), (0, F(0), 0))
    moved = g3.deck(p, 1)
    before = g3.compensated(p[0], p[2][0], p[2][1])
    after = g3.compensated(moved[0], moved[2][0], moved[2][1])
    assert before[1] == 0 and after[1] == F(-4, 5)

    def repaired_compensation(u, phi, sigma):
        return ((phi + g3.P_PHI * u) % g3.KMOD, sigma + g3.P_SIGMA * u)

    checked = 0
    for u in (F(0), F(1, 3), F(-2)):
        for phi in (F(0), F(5, 4), F(-1)):
            for sigma in (F(0), F(3), F(-2, 7)):
                pt = (u, F(0), (phi, sigma, 2))
                for n in range(-3, 4):
                    d = g3.deck(pt, n)
                    assert repaired_compensation(d[0], d[2][0], d[2][1]) == repaired_compensation(pt[0], pt[2][0], pt[2][1])
                    checked += 1

    b1 = [[F(3), F(1)], [F(0), F(2)]]
    b2 = [[F(5), F(0)], [F(1), F(4)]]
    product = t1.mm(b1, b2)
    transport = [[F(2), F(0)], [F(0), F(2)]]
    assert t1.mm(transport, product) == t1.mm(product, transport)

    def coupling(t):
        return t1.det2(t1.mm(t, product)) / (t1.det2(t1.mm(t, b1)) * t1.det2(t1.mm(t, b2)))

    factor = coupling(transport)
    assert factor == F(1, 4)
    shear = [[F(1), F(1)], [F(0), F(1)]]
    assert t1.mm(shear, product) != t1.mm(product, shear)
    assert coupling(shear) == 1

    warp = g2.A_rational(F(-2))
    values = [g2.peval(warp, v) for v in (F(-1), F(0), F(1))]
    assert values == [F(-1), F(1), F(-1)]
    boundary = g2.holonomy_boundary(warp, F(0), F(2), F(-1), F(1))
    area = g2.holonomy_area(warp, F(0), F(2), F(-1), F(1))
    assert boundary == area == 16
    return {
        'EMKG3_T2': {
            'existing_certificate_reports': g3.certify_T2()['verdict'],
            'actual_deck_compensated_sigma_before': str(before[1]),
            'actual_deck_compensated_sigma_after': str(after[1]),
            'drift': str(after[1] - before[1]),
            'finding': 'Inverse deck and forward compensation use opposite signs.',
            'minimal_repair': 'Keep inverse deck; phi_hat=(phi+p_phi*u) mod K, sigma_hat=sigma+p_sigma*u.',
            'repair_exact_checks': checked,
        },
        'EMKT1_T6': {
            'commuting_transport': '2I', 'coupling_factor': str(factor),
            'noncommuting_determinant_one_transport_coupling_factor': str(coupling(shear)),
            'finding': 'The inverse transport determinant is not a noncommutativity diagnostic; its arithmetic is correct.',
        },
        'EMKG2_T5': {
            'c': '-2', 'v_interval': ['-1', '1'],
            'A_at_minus_one_zero_one': [str(x) for x in values],
            'formal_boundary_integral': str(boundary), 'formal_area_integral': str(area),
            'finding': 'A vanishes inside this rectangle; formal polynomial agreement does not establish Levi-Civita holonomy across a degenerate metric.',
            'minimal_repair': 'Restrict geometric holonomy checks to a declared nondegenerate region, here A>0.',
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publications-root', type=Path, required=True)
    parser.add_argument('--rkf-root', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if sys.version_info[:2] not in ((3, 11), (3, 12)):
        parser.error('Use Python 3.11 or 3.12.')
    here = Path(__file__).resolve().parent
    repo = here.parent
    pins = json.loads((here / 'EMK_MASTER_REVIEW_SOURCE_PINS.json').read_text())
    root = args.publications_root.resolve()
    count = check_sources(root, pins['publications']['files'])
    rkf_count = check_sources(args.rkf_root.resolve(), pins['rkf_context']['files']) if args.rkf_root else None
    suites = []
    test_paths = sorted(x['path'] for x in pins['publications']['files'] if '/tests/test_' in x['path'])
    assert len(test_paths) == 8
    for i, path in enumerate(test_paths):
        mod = load('emk_master_review_tests_' + str(i), root / path)
        functions = [fn for name, fn in inspect.getmembers(mod, inspect.isfunction) if name.startswith('test_')]
        assert all(not inspect.signature(fn).parameters for fn in functions)
        suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
        stream = io.StringIO()
        result = unittest.TextTestRunner(stream=stream).run(suite)
        if not result.wasSuccessful():
            raise AssertionError(path + '\n' + stream.getvalue())
        suites.append({'path': path, 'tests': result.testsRun, 'passed': True})
    assert sum(x['tests'] for x in suites) == 191
    certificate_pins = []
    for item in pins['publications']['files']:
        path = root / item['path']
        if path.name.startswith('EXPECTED_') and path.suffix == '.sha256':
            family = path.stem.removeprefix('EXPECTED_')
            expected = path.read_text().split()[0]
            actual = sha256(path.with_name(family + '_RESULT.json'))
            assert actual == expected
            certificate_pins.append({'family': family, 'sha256': actual, 'matches': True})
    assert len(certificate_pins) == 8
    findings = cross_checks(root)
    assert check_sources(root, pins['publications']['files']) == count
    preserved = {}
    for revision in range(1, 8):
        record = json.loads((here / ('R' + str(revision) + '_VERIFICATION.json')).read_text())
        assert record['status'].startswith('PASS')
        assert all(sha256(repo / path) == digest for path, digest in record['sha256'].items())
        preserved['R' + str(revision)] = {'files_checked': len(record['sha256']), 'all_recorded_hashes_match': True}
    report = {
        'status': 'PASS_REPRODUCIBLE_EMK_MASTER_REVIEW_WITH_THREE_FINDINGS',
        'review_date': '2026-10-01', 'publications_commit': pins['publications']['commit'],
        'unchanged_publications_files': count, 'checked_rkf_context_files': rkf_count,
        'runner': 'stdlib unittest.FunctionTestCase on unchanged no-argument upstream test functions',
        'suites': suites, 'total_upstream_tests': 191, 'certificate_pins': certificate_pins,
        'findings': findings, 'upstream_repairs_applied': False, 'preserved_revisions': preserved,
        'sha256': {p: sha256(repo / p) for p in (
            '02-relational-response/EMK_MASTER_TENSOR_REVIEW.md',
            '04-operator-evolution/verify_emk_master_review.py',
            '04-operator-evolution/EMK_MASTER_REVIEW_SOURCE_PINS.json')},
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        output = args.output.resolve()
        if output.is_relative_to(root) or (args.rkf_root and output.is_relative_to(args.rkf_root.resolve())):
            parser.error('Output must be outside the upstream source checkouts.')
        output.write_text(payload)
    print(payload, end='')


if __name__ == '__main__':
    main()
