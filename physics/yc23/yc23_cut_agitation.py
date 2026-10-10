"""Exact controls for YC23's full-carrier written proofs, Python 3.12."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = (
    'physics/yc20/YC20_METRIC_RETAINED_RETURN.md',
    'physics/yc21/YC21_TWO_CUBE_BLOCK_AND_REMOVED_CHANNELS.md',
    'physics/yc22/YC22_GAUGE_PROTECTED_CHANNELS.md',
    'uncut/up8/UP8_DEFORMED_COMPASS_AND_SEAM_TRANSPORT.md',
    '02-relational-response/EMK_TOPOLOGY_R15.md',
)


def zero(x):
    if isinstance(x, sp.MatrixBase):
        return all(zero(v) for v in x)
    return sp.trigsimp(sp.simplify(x)) == 0


def comm(a, b):
    return a*b-b*a


def sharp_data(h, omega, cut):
    if not zero(h-h.conjugate().T) or not zero(h*omega):
        raise ValueError('Hermitian shifted Hamiltonian and null vacuum required')
    if not zero((omega.conjugate().T*omega)[0]-1):
        raise ValueError('normalized vacuum required')
    eye = sp.eye(h.rows)
    if not zero(cut-cut.conjugate().T) or not zero(cut*cut-eye):
        raise ValueError('self-adjoint involution required')
    rho = omega*omega.conjugate().T
    plus, minus = (eye+cut)/2, (eye-cut)/2
    post = plus*rho*plus+minus*rho*minus
    energy = sp.simplify(sp.trace(h*post))
    excite = sp.simplify(1-(omega.conjugate().T*post*omega)[0])
    direct = sp.simplify((omega.conjugate().T*cut*h*cut*omega)[0]/2)
    double = sp.simplify((omega.conjugate().T*comm(cut, comm(h, cut))*omega)[0]/4)
    second = sp.simplify((comm(h, cut)*omega).conjugate().dot(comm(h, cut)*omega))
    return energy, excite, direct, double, second


def quat(a, b):
    return (a[0]*b[0]-sum(a[i]*b[i] for i in range(1, 4)),
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


def conjugate(q):
    return (q[0], -q[1], -q[2], -q[3])


def plaquette():
    qs = [sp.symbols(f'q{e}_0:4', real=True) for e in range(4)]
    w = quat(quat(quat(qs[0], qs[1]), conjugate(qs[2])), conjugate(qs[3]))[0]
    return qs, sp.expand(w)


def sphere_laplacian(f, q):
    euler = sum(x*sp.diff(f, x) for x in q)
    return sp.expand(sum(sp.diff(f, x, 2) for x in q)
                     -sum(x*sp.diff(euler, x) for x in q)-2*euler)


def haar_moment(n):
    if type(n) is not int or n < 0:
        raise ValueError('nonnegative integer moment required')
    return sp.S.Zero if n % 2 else sp.catalan(n//2)/4**(n//2)


def haar_average(poly, c):
    return sp.expand(sum(coef*haar_moment(power[0])
                         for power, coef in sp.Poly(poly, c).terms()))


def free_energy_from_s(s):
    s = sp.Rational(s)
    if not 0 < s <= 1:
        raise ValueError('sqrt(1-eta²) must belong to (0,1]')
    return 1-2*s*s/(1+s)


def graph_ground():
    omega = sp.Matrix([sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(2, 3)])
    edges = {(0, 1): sp.Rational(1, 9), (1, 2): sp.Rational(2, 9),
             (0, 2): sp.Rational(3, 9)}
    h = sp.zeros(3)
    for (i, j), c in edges.items():
        h[i, j] = h[j, i] = -c/(omega[i]*omega[j])
        h[i, i] += c/omega[i]**2
        h[j, j] += c/omega[j]**2
    return omega, h, edges


def certificate():
    checks = {}
    def check(name, condition):
        checks[name] = bool(condition)
        if not checks[name]:
            raise AssertionError(name)

    h, om = sp.diag(0, 5), sp.Matrix([1, 0])
    k = sp.Matrix([[sp.Rational(3, 5), sp.Rational(4, 5)],
                   [sp.Rational(4, 5), -sp.Rational(3, 5)]])
    e, q, direct, double, second = sharp_data(h, om, k)
    check('projective instrument energy and double commutator', e == direct == double == sp.Rational(8, 5))
    check('vacuum escape probability and gap normalization', q == sp.Rational(8, 25) and e/q == 5)
    check('first and second energy moments are distinct', second == 16 and second != e)
    e0 = sharp_data(sp.diag(0, 1, 3), sp.Matrix([1, 0, 0]),
                    sp.diag(1, sp.Matrix([[0, 1], [1, 0]])))[0]
    check('noncommutation elsewhere does not agitate the vacuum', e0 == 0)
    check('rescaling the energy preserves escape probability', sharp_data(7*h, om, k)[:2] == (7*e, q))

    u = sp.symbols('u', real=True)
    mp, mm = sp.sqrt((1+u)/2), sp.sqrt((1-u)/2)
    metric = sp.diff(mp, u)**2+sp.diff(mm, u)**2
    fisher = sum(sp.diff(p, u)**2/p for p in ((1+u)/2, (1-u)/2))
    check('smooth branch normalization', zero(mp**2+mm**2-1))
    check('exact local cut energy density', zero(metric-1/(4*(1-u*u))))
    check('Fisher normalization is four times the metric', zero(fisher-4*metric))
    a = sp.symbols('a', real=True)
    w = sp.Matrix([sp.cos(a/2), sp.sin(a/2)])
    p = w*w.T
    kc = 2*p-sp.eye(2)
    dw = w.diff(a)
    check('lifted cut is an involution', zero(kc*kc-sp.eye(2)))
    check('fixed pointer compresses to the smooth reading', zero((w.T*sp.diag(1, -1)*w)[0]-sp.cos(a)))
    check('horizontal cut metric has quarter normalization', zero((dw.T*(sp.eye(2)-p)*dw)[0]-sp.Rational(1, 4)))
    check('cut derivative norm has eighth normalization', zero(sp.trace(kc.diff(a).T*kc.diff(a))/8-sp.Rational(1, 4)))
    check('projected connection is flat despite positive cut metric', zero(w.T*dw) and zero((dw.T*dw)[0]-sp.Rational(1, 4)))
    r = sp.Matrix([[0, -1], [1, 0]])
    check('native four-step return is identity for the same cut', zero((r*kc)**2-sp.eye(2)))

    eta, c = sp.symbols('eta c', real=True)
    mu2, mu3, mu4 = sp.symbols('mu2 mu3 mu4', real=True)
    avp = (1-mu2*eta**2/8+mu3*eta**3/16-5*mu4*eta**4/128)/sp.sqrt(2)
    avm = avp.subs(eta, -eta)
    qseries = sp.series(1-avp**2-avm**2, eta, 0, 5).removeO().expand()
    check('centered probe escape expansion', zero(qseries-mu2*eta**2/4-(5*mu4-mu2**2)*eta**4/64))
    fseries = sp.series(sp.sqrt((1+eta*c)/2), eta, 0, 5).removeO()
    av = haar_average(fseries, c)
    qhaar = sp.series(1-av**2-av.subs(eta, -eta)**2, eta, 0, 5).removeO()
    check('actual Haar escape expansion', zero(qhaar-eta**2/16-sp.Rational(9, 1024)*eta**4))

    qs, wp = plaquette()
    norms = [sum(x*x for x in qq) for qq in qs]
    check('four-link fundamental electric energy is twelve',
          zero(-sum(sphere_laplacian(wp, qq) for qq in qs)-12*wp))
    check('four exact single-link gradient identities', all(
        sp.expand(sum(sp.diff(wp, x)**2 for x in qs[i])-sp.prod(norms[j] for j in range(4) if j != i)) == 0
        for i in range(4)))
    check('Haar second and fourth moments', haar_moment(2) == sp.Rational(1, 4) and haar_moment(4) == sp.Rational(1, 8))
    check('plaquette weak-probe ratio is twelve', 4*(1-haar_moment(2))/haar_moment(2) == 12)
    t = sp.symbols('t')
    closed = 1-2*(1-t)/(1+sp.sqrt(1-t))
    series = sp.series(closed, t, 0, 8).removeO().expand()
    check('closed Haar energy agrees with seven exact moments',
          zero(series-sum((haar_moment(2*n)-haar_moment(2*n+2))*t**(n+1) for n in range(7))))
    check('energy weak expansion', zero(sp.series(closed, t, 0, 3).removeO()-3*t/4-t*t/8))
    check('finite Haar measurement is not an energy gap',
          all(0 < free_energy_from_s(s) <= 1-s*s for s in (sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(12, 13))))

    omega, hg, edges = graph_ground()
    fs = sp.symbols('f0:3', real=True)
    fv = sp.Matrix(fs)
    psi = sp.diag(*omega)*fv
    check('nonconstant true-ground form control', zero(hg*omega) and zero(
        (psi.T*hg*psi)[0]-sum(cc*(fs[i]-fs[j])**2 for (i, j), cc in edges.items())))
    return dict(stage='YC23', checks=checks, exact_check_count=len(checks),
                scope='written full finite-lattice proof; exact algebra and normalization controls',
                source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
                free_plaquette=dict(variance='1/4', dirichlet_energy='3', normalized_weak_probe_energy='12',
                                    energy_at_eta_3_5=str(free_energy_from_s(sp.Rational(4, 5)))),
                claims=dict(cut_energy_information_identity=True, calibrated_physical_gap_identity=True,
                            primitive_lambda_space_uniquely_derived=False, spontaneous_event_derived=False,
                            larger_gap_window=False, continuum_mass_gap=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = certificate()
    if args.write:
        (HERE/'YC23_RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC23: {result['exact_check_count']} exact checks pass; no enlarged gap window.")
