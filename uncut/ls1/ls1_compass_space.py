"""Exact controls for the compass target's written geometric construction."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
I = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.diag(1, -1)
J = R*K
r, phi = sp.symbols('r phi', real=True)
SOURCES = (
    'uncut/up7/UP7_TURN_AND_CUT.md',
    'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
    'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
    'physics/yc23/YC23_PRIMITIVE_CUT_AND_AGITATION.md',
    '02-relational-response/OBSERVER_METRIC_FOUNDATION_R12.md',
    '02-relational-response/R14_SCOPE_CORRECTION.md',
)


def simplify(x):
    return sp.trigsimp(sp.simplify(sp.expand(x)))


def zero(x):
    return all(zero(v) for v in x) if isinstance(x, sp.MatrixBase) else simplify(x) == 0


def chart(depth=r, angle=phi):
    k = sp.cos(angle)*K+sp.sin(angle)*J
    return sp.sqrt(1+depth**2)*k+depth*R


def metric(depth=r):
    return sp.diag((1+2*depth**2)/(1+depth**2), 1+depth**2)


def gauss_curvature(depth=r):
    return -1/(1+2*depth**2)**2


def rational_cut(boost, cosine, sine):
    b, c, s = map(sp.Rational, (boost, cosine, sine))
    if b <= 0 or c*c+s*s != 1:
        raise ValueError('positive boost and unit angular pair required')
    depth, height = (b-1/b)/2, (b+1/b)/2
    return height*(c*K+s*J)+depth*R


def chord_squared(a, b):
    d = a-b
    return sp.simplify(sp.trace(d.T*d)/2)


def certificate():
    checks = {}
    def check(name, condition):
        checks[name] = bool(condition)
        if not checks[name]:
            raise AssertionError(name)

    u, v, w = sp.symbols('u v w', real=True)
    raw = u*K+v*J+w*R
    check('native quadratic classification', zero(raw*raw-(u*u+v*v-w*w)*I))
    check('smooth constraint has no critical point on its level', sp.expand((u*u+v*v+w*w)-(1+2*w*w)-(u*u+v*v-w*w-1)) == 0)
    kap = chart()
    check('global cylinder chart is an involution', zero(kap*kap-I))
    check('depth is the signed anticommutator', zero(kap*R+R*kap+2*r*I))
    check('zero cuts still do not commute with the turn', not zero(K*R-R*K))
    coords = (r, phi)
    g = sp.Matrix([[simplify(sp.trace(kap.diff(a).T*kap.diff(b))/2) for b in coords] for a in coords])
    check('induced positive metric formula', zero(g-metric()))
    alg = sp.Matrix([[simplify(sp.trace(kap.diff(a)*kap.diff(b))/2) for b in coords] for a in coords])
    check('algebra trace form has Lorentzian signature', zero(alg-sp.diag(-1/(1+r*r), 1+r*r)))
    check('metric determinant and regular zero seam', zero(g.det()-(1+2*r*r)) and g.subs(r, 0) == I)
    check('radial coefficient dominates cylinder metric', zero(g[0, 0]-1-r*r/(1+r*r)))

    k = sp.cos(phi)*K+sp.sin(phi)*J
    n = R*k
    root = sp.sqrt(1+r*r)*I-r*n
    root_inv = sp.sqrt(1+r*r)*I+r*n
    b = (1+2*r*r)*I-2*r*sp.sqrt(1+r*r)*n
    bi = (1+2*r*r)*I+2*r*sp.sqrt(1+r*r)*n
    check('native return and positive determinant-one root', zero((R*kap)**2-b) and zero(root*root-b) and zero(b.det()-1))
    check('return inverse formula', zero(b*bi-I) and zero(root*root_inv-I))
    check('entire zero circle has identity return', b.subs(r, 0) == I)
    check('opposite labelled cuts have the same return', zero((R*(-kap))**2-b))
    xs = [bi*b.diff(a) for a in coords]
    gb = sp.Matrix([[simplify(sp.trace(x*y)/8) for y in xs] for x in xs])
    check('return-space metric loses zero angle', zero(gb-sp.diag(1/(1+r*r), r*r*(1+r*r))) and zero(gb.det()-r*r))
    aa = [(root*x*root_inv/2-root.diff(a)*root_inv).applyfunc(simplify) for a, x in zip(coords, xs)]
    check('global native connection from full return', zero(aa[0]) and zero(aa[1]-r*r*R))
    curvature = aa[1].diff(r)-aa[0].diff(phi)+aa[0]*aa[1]-aa[1]*aa[0]
    check('native curvature linear in signed depth', zero(curvature-2*r*R))
    check('metric-normalized native curvature', zero(sp.trace(curvature.T*curvature)/2/g.det()-4*r*r/(1+2*r*r)))

    f = sp.sqrt(g[1, 1])
    kg = simplify(-sp.diff(sp.diff(f, r)/sp.sqrt(g[0, 0]), r)/(f*sp.sqrt(g[0, 0])))
    check('Gaussian curvature differs from native transport', zero(kg-gauss_curvature()) and kg.subs(r, 0) == -1)
    scale = sp.diag(2, sp.Rational(1, 2))
    tangent = sp.Matrix([[0, 1], [0, 0]])
    check('noncompact stabilizer expands a nonzero tangent', scale*K*scale.inv() == K and scale*tangent*scale.inv() == 4*tangent and zero(K*tangent+tangent*K))

    eigen = sp.Matrix([sp.cos(phi/2), sp.sin(phi/2)])
    check('oriented reset needs double cover', zero(k*eigen-eigen) and zero(eigen.subs(phi, phi+2*sp.pi)+eigen))
    p, q, ww = 3, 4, 1
    delta = p*p+q*q-ww*ww
    kc = (p*K+q*J+ww*R)/sp.sqrt(delta)
    mapped = sp.sqrt(1+sp.Rational(ww*ww, delta))*(sp.Rational(3, 5)*K+sp.Rational(4, 5)*J)+R/sp.sqrt(delta)
    check('actual response normalization maps into target', zero(kc-mapped) and zero(kc*kc-I))
    check('regular real depth one despite zero-series radius', metric(sp.S.One) == sp.diag(sp.Rational(3, 2), 2) and zero(gauss_curvature(sp.S.One)+sp.Rational(1, 9)))
    check('complex branch locations set chart-series radius', (1+r*r).subs(r, sp.I) == 0 and (1+r*r).subs(r, -sp.I) == 0)

    return dict(stage='LS1', exact_check_count=len(checks), checks=checks,
                source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
                target='all labelled nontrivial real two-mode cuts; R times S1',
                positive_metric='diag((1+2r^2)/(1+r^2), 1+r^2)',
                native_connection='r^2 R dphi', native_curvature='2r R dr wedge dphi',
                metric_gaussian_curvature='-1/(1+2r^2)^2',
                proof_scope='topology/completeness by written proof; exact algebra controls',
                claims=dict(zero_seam_retains_circle=True, complete_positive_metric=True,
                            global_native_connection=True, full_primitive_lambda_space_selected=False,
                            canonical_tower_interpolation=False, spacetime_replaced=False,
                            new_yang_mills_gap=False, continuum_mass_gap=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = certificate()
    if args.write:
        (HERE/'LS1_RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"LS1: {result['exact_check_count']} exact checks pass; zero seam retained.")
