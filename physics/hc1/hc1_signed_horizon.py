"""HC1: exact audit of causal signs and native/geometric curvature contractions.

The global horizon statement uses the written stationary-extension argument.
No Yang--Mills spectral bound is computed or altered by this packet.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import argparse
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def signed_reading(beta):
    beta = F(beta)
    if beta <= 0:
        raise ValueError('finite positive fall parameter required')
    return 2*beta/(1+beta*beta), (1-beta*beta)/(1+beta*beta)


def radial_speeds(beta):
    beta = F(beta)
    if beta <= 0:
        raise ValueError('finite positive fall parameter required')
    return 1-beta, -1-beta


def exterior_curvatures(r, rs=1):
    r, rs = F(r), F(rs)
    if not r > rs > 0:
        raise ValueError('the positive static response requires r > r_s > 0')
    radial = -rs**2/(r**2*(r-rs)**4)
    return {'radial_angular': radial, 'angular_angular': 4*radial,
            'native_diagnostic': -6*radial, 'geometric_square': 12*rs**2/r**6,
            'clock_squared': 1-rs/r}


@lru_cache(None)
def symbolic_checks():
    checks = {}
    b, r, rs, eps = sp.symbols('beta r r_s epsilon', positive=True)
    u = 2*b/(1+b*b)
    v = (1-b*b)/(1+b*b)
    zero = lambda x: sp.simplify(x) == 0
    checks['signed circle'] = zero(u*u+v*v-1)
    checks['inverse reconstructs beta'] = zero(u/(1+v)-b)
    checks['inverse reconstructs radius'] = zero((1+v)/(1-v)-1/b**2)
    checks['mirror preserves folded reading'] = zero(u.subs(b, 1/b)-u)
    checks['mirror reverses signed reading'] = zero(v.subs(b, 1/b)+v)
    checks['outgoing sign multiplier is positive for beta positive'] = zero(
        1-b-v*(1+b*b)/(1+b))
    g = sp.Matrix([[1-b*b, -b], [-b, -1]])
    checks['radial determinant stays minus one'] = g.det() == -1
    checks['time gradient has unit squared norm'] = zero(g.inv()[0, 0]-1)
    checks['both radial null roots'] = all(
        zero((sp.Matrix([1, speed]).T*g*sp.Matrix([1, speed]))[0])
        for speed in (1-b, -1-b))
    checks['falling clock has unit proper-time rate'] = zero(
        (sp.Matrix([1, -b]).T*g*sp.Matrix([1, -b]))[0]-1)
    f = (r+rs)**2/(4*r)
    checks['tower exterior fold'] = zero(f-rs-(r-rs)**2/(4*r))
    checks['tower radius mirror'] = zero(f.subs(r, rs**2/r)-f)
    checks['tower orientation reverses at horizon'] = zero(
        sp.diff(f, r)-(1-rs**2/r**2)/4)
    checks['escape deficit squares'] = zero(1-u-(1-b)**2/(1+b*b))
    checks['clock at next tower level'] = zero(1-u*u-v*v)

    # Independent 2x2 realization at a pole; angular derivatives are retained.
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    c, s = (1+b*b)/(1-b*b), 2*b/(1-b*b)
    dp = -b/(r*(1-b*b))
    inv = c*sp.eye(2)-s*sz
    X = [inv*s*sx/r, inv*s*sy/r, dp*sz]
    squares = []
    for i, j in ((0, 2), (1, 2), (0, 1)):
        curvature = -(X[i]*X[j]-X[j]*X[i])/4
        squares.append(sp.simplify(sp.trace(curvature*curvature)/2))
    expected = [-(dp*s)**2/(4*r*r)]*2+[-s**4/(4*r**4)]
    checks['matrix curvature preserves both transverse directions'] = all(
        zero(a-e) for a, e in zip(squares, expected))
    checks['fixed base squaring radial angular amplification'] = zero(
        -((2*dp)*(2*s*c))**2/(4*r*r)-16*c*c*expected[0])
    checks['fixed base squaring angular angular amplification'] = zero(
        -(2*s*c)**4/(4*r**4)-16*c**4*expected[2])
    checks['fixed base squaring changes native ratio'] = zero(
        r*(2*dp)/(2*s*c)+v/2)
    bp = sp.sqrt(rs/r)
    new_memory = 4*rs*r/(r+rs)**2
    checks['fixed base tower fails the unchanged vacuum profile'] = zero(
        sp.diff(r*new_memory, r)-8*rs**2*r/(r+rs)**3)
    checks['response rapidity derivative includes factor two'] = zero(
        2*sp.diff(bp, r)/(1-bp**2)-dp.subs(b, bp))
    ira = -rs**2/(r*r*(r-rs)**4)
    checks['native radial angular formula'] = zero(squares[0].subs(b, bp)-ira)
    checks['native angular angular formula'] = zero(squares[2].subs(b, bp)-4*ira)
    D, K, N2 = -6*ira, 12*rs**2/r**6, 1-rs/r
    checks['native geometric observer weight'] = zero(D*N2**4-K/2)
    checks['native fourth order horizon divergence'] = (
        sp.limit((r-rs)**4*D, r, rs, dir='+') == 6)
    checks['geometric horizon remains finite'] = zero(K.subs(r, rs)-12/rs**4)
    checks['geometric singular centre'] = sp.limit(K, r, 0, dir='+') == sp.oo
    checks['native monotonicity logarithmic derivative'] = zero(
        sp.diff(D, r)/D+2/r+4/(r-rs))
    primitive = -8*sp.pi*rs**2/(r-rs)**3
    checks['native diagnostic spatial measure primitive'] = zero(
        sp.diff(primitive, r)-4*sp.pi*r*r*D)
    checks['cutoff divergence coefficient'] = zero(
        primitive.subs(r, rs+eps)+8*sp.pi*rs**2/eps**3)
    # w>1 is required for the real logarithm; derivative identity is algebraic.
    w = sp.symbols('w', positive=True)
    tstar = rs*(w*w+2*w+2*sp.log(w-1))
    checks['outgoing travel-time primitive'] = zero(
        sp.diff(tstar, w)/(2*rs*w)-1/(1-1/w))
    speed = 1-bp
    checks['simple outgoing horizon zero'] = zero(sp.diff(speed, r).subs(r, rs)-1/(2*rs))
    checks['falling crossing time primitive'] = zero(
        sp.diff(2*r**sp.Rational(3, 2)/(3*sp.sqrt(rs)), r)-1/bp)
    checks['folded centre has two limiting sides'] = (
        sp.limit(u, b, 0) == sp.limit(u, b, sp.oo) == 0 and
        sp.limit(v, b, 0) == 1 and sp.limit(v, b, sp.oo) == -1)
    return checks


def source_pins():
    names = [
        'physics/mo1/MO1_MOTION.md',
        'physics/ma1/MA1_MO1_ASSUMPTION_FROM_THE_LAW.md',
        'physics/gr1/GR1_MEMORY_DICTIONARY.md',
        'physics/gb1/GB1_GRAVITY_IN_THE_BLOCK.md',
        'physics/cv1/CV1_FRAME_CURVATURE.md',
        'physics/nc1/NC1_NATIVE_CURVATURE_OF_THE_FIELD.md',
        'physics/nc1/nc1_native_curvature_of_the_field.py',
        'physics/lt1/LT1_LAMBDA_TOWER_INWARD.md',
        'physics/lt1/lt1_lambda_tower_inward.py',
        'physics/ts1/TS1_TWO_SHEET_CLOSURE.md',
        'physics/hm1/HM1_HORIZON_MEMORY_PRODUCT.md',
        'uncut/hx1/HX1_FALL_IS_A_TURN_OF_CUTS.md',
        'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
        'physics/hc1/HC1_SIGNED_HORIZON_AND_NATIVE_CURVATURE.md',
        'physics/hc1/hc1_signed_horizon.py',
        'physics/hc1/test_hc1.py',
    ]
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def encode(obj):
    if isinstance(obj, F): return str(obj)
    if isinstance(obj, dict): return {k: encode(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)): return [encode(v) for v in obj]
    return obj


def run():
    checks = dict(symbolic_checks())
    mirror = []
    for beta in (F(1, 2), F(1), F(2)):
        u, v = signed_reading(beta)
        mirror.append({'beta': beta, 'u': u, 'v': v, 'speeds': radial_speeds(beta)})
    checks['same folded reading opposite escape signs rational witness'] = (
        mirror[0]['u'] == mirror[2]['u'] == F(4, 5) and
        mirror[0]['speeds'][0] == F(1, 2) and mirror[2]['speeds'][0] == -1)
    checks['horizon has one zero and one inward null root'] = radial_speeds(1) == (0, -2)
    beta, orbit = F(1, 3), []
    for _ in range(4):
        u, v = signed_reading(beta)
        r = 1/beta**2
        orbit.append({'beta': beta, 'r_over_rs': r, 'escape_speed': 1-beta,
                      **exterior_curvatures(r)})
        if not (1-beta)**2/2 <= 1-u <= (1-beta)**2:
            raise AssertionError('escape-step inequality')
        beta = u
    checks['rational tower decreases escape speed and radius'] = all(
        a['r_over_rs'] > b['r_over_rs'] > 1 and a['escape_speed'] > b['escape_speed'] > 0
        for a, b in zip(orbit, orbit[1:]))
    checks['tower native diagnostic grows while geometric square stays below horizon'] = all(
        a['native_diagnostic'] < b['native_diagnostic'] and b['geometric_square'] < 12
        for a, b in zip(orbit, orbit[1:]))
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    return encode({'stage': 'HC1', 'base_commit': 'b8a4805', 'checks': checks,
        'signed_mirror_control': mirror, 'exterior_tower': orbit,
        'exact_relation': 'D = K/(2*N^8), exterior spherical profile only',
        'claim_boundary': {'conditional_spherical_horizon': True,
            'native_curvature_divergent_at_horizon': True,
            'spacetime_curvature_divergent_at_horizon': False,
            'infinite_physical_energy_proved': False,
            'tower_is_infalling_time_evolution': False,
            'positive_static_response_continued_inside': False,
            'fixed_base_squaring_preserves_spherical_vacuum_law': False,
            'ST_centre_identified_with_black_hole': False,
            'new_YM_gap_or_continuum_result': False,
            'finite_checks_are_formal_verification': False},
        'source_pins': source_pins()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'HC1_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"HC1: {len(result['checks'])} exact checks pass.")
    print('Signed causal reading and D=K/(2*N^8) verified on the declared profile; no new YM gap.')
