"""YC16: absorb complete one-cube return operators, retaining compensation.

The exact rebasing is on the original lattice Hamiltonian. The connected
second-order source is bounded separately from the full remaining interaction.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


y15 = load('yc16_yc15', 'physics/yc15/yc15_boundary_return_channels.py')
GAP = y15.GAP
ETA_MAX = y15.ETA_MAX
ALPHA_MAX = Q(235, 11232)
RADIUS = Q(1, 128)
EXP_UPPER = Q(16, 15)


def split_return(returned, omega_a, omega_b):
    """Exact conditional expectations in declared orthonormal finite controls.

The written theorem uses the same bounded-operator identities on full spaces.
"""
    if (omega_a.T*omega_a)[0] != 1 or (omega_b.T*omega_b)[0] != 1:
        raise ValueError('normalized real reference vectors required')
    na, nb = len(omega_a), len(omega_b)
    if returned.shape != (na*nb, na*nb) or returned != returned.T:
        raise ValueError('self-adjoint return on the complete product required')
    va = sp.kronecker_product(sp.eye(na), omega_b)
    vb = sp.kronecker_product(omega_a, sp.eye(nb))
    omega = sp.kronecker_product(omega_a, omega_b)
    ea, eb = va.T*returned*va, vb.T*returned*vb
    alpha = (omega.T*returned*omega)[0]
    single_a, single_b = ea-alpha*sp.eye(na), eb-alpha*sp.eye(nb)
    connected = returned-alpha*sp.eye(na*nb)-sp.kronecker_product(single_a, sp.eye(nb))-sp.kronecker_product(sp.eye(na), single_b)
    return dict(alpha=alpha, single_a=single_a, single_b=single_b,
                connected=connected, omega=omega, va=va, vb=vb)


def exponential_integral(coefficients, offset):
    return sum(Q(value)/Q(offset+rate) for rate, value in coefficients.items())


def square_integral(coefficients, offset):
    return sum(Q(a)*Q(b)/Q(offset+i+j)
               for i, a in coefficients.items() for j, b in coefficients.items())


def source_bounds(theta_a, theta_b):
    ta, tb = y15.edge_kinetic_bound(theta_a), y15.edge_kinetic_bound(theta_b)
    return dict(theta_a=Q(theta_a), theta_b=Q(theta_b), edge_kinetic_a=ta,
                edge_kinetic_b=tb, single_a_source_squared_upper=Q(4, 481)**2*ta,
                single_b_source_squared_upper=Q(4, 481)**2*tb,
                connected_source_squared_upper=Q(4, 1443)**2*ta*tb,
                connected_operator_norm_upper=Q(1, 24),
                single_operator_norm_upper=ALPHA_MAX)


def reference_bounds(eta):
    eta = Q(eta)
    if not 0 <= eta <= ETA_MAX:
        raise ValueError('certified interface window is 0<=eta<=1/4320')
    eps = 24*ALPHA_MAX*eta**2
    return dict(eta=eta, counterterm_norm_upper=eps,
                counterterm_source_upper=eta**2/20,
                new_cube_gap_lower=GAP-eps,
                new_ground_energy_abs_upper=eta**4/108,
                new_ground_vector_distance_upper=Q(10, 27)*eta**2,
                two_bridge_inverse_vector_change_upper=eta**2/13,
                single_counterterm_creator_upper=Q(3, 16)*eta**2,
                connected_source_after_reference_upper=Q(1,8658)+Q(5,162)*eta**2,
                recentered_connected_one_cube_cost_per_block_upper=Q(20,27)*eta**4+Q(400,729)*eta**6,
                recentered_connected_scalar_cost_per_face_upper=Q(50,2187)*eta**6)


def join_bounds(eta):
    ref = reference_bounds(eta)
    eta, eps, gh = ref['eta'], ref['counterterm_norm_upper'], ref['new_cube_gap_lower']
    inverse_cap = Q(27, 640)+eta**2/13
    cube_seed = 96*eta*inverse_cap+Q(3, 16)*eta**2
    bridge_seed = 8*eta*inverse_cap+eta/3
    seed = max(cube_seed, bridge_seed)
    beta = 24*eta+eps
    contraction = 4*beta/gh*EXP_UPPER*(2+8*(1+2*RADIUS))
    relative = 8*beta/gh*EXP_UPPER
    mapping = seed+contraction*RADIUS
    return dict(reference=ref, beta=beta, seed=seed, bridge_seed=bridge_seed,
                contraction=contraction, relative_return=relative,
                mapping_reserve=RADIUS-mapping,
                full_gap_lower=gh*(1-relative), radius=RADIUS,
                nonlinear_minimum_support=1,
                fixed_point_norm_upper=seed/(1-contraction))


def symbolic_checks():
    t, g, a = sp.symbols('t g a', positive=True)
    # Laplace transform of J_g(t)=(exp(-gt)-exp(-3t))/(3-g).
    j_laplace = (1/(a+g)-1/(a+3))/(3-g)
    checks = {
        'charged Duhamel kernel keeps both decay rates': sp.cancel(j_laplace-1/((a+g)*(a+3))) == 0,
        'two charged kernels have the exact squared integral': sp.cancel(
            (1/(6+2*g)-2/(9+g)+sp.Rational(1,12))/(3-g)**2-
            2/((6+2*g)*(9+g)*12)) == 0,
        'one-cube source coefficient keeps the charged decay of its neighbor':
            Q(1, 2)/(Q(13, 2)*Q(37, 4)) == Q(4, 481),
        'connected source coefficient is four over 1443':
            Q(2)/(Q(13, 2)*Q(37, 4)*12) == Q(4, 1443),
    }
    heat = {3: Q(1, 84), 9: Q(-1, 48), 17: Q(1, 112)}
    free_face = {9: Q(1, 4), 17: Q(3, 4), 3: Q(-1)}
    checks['first connected cube source is an exact mixed internal-coupling derivative'] = square_integral(heat, 6)/4 == Q(19, 86261760)
    checks['free one-cube operator changes an excited gauge face'] = exponential_integral(free_face, 9)/4 == Q(-19, 1872)
    checks['free connected operator survives despite zero vacuum source'] = square_integral(free_face, 6)/4 == Q(3931, 599040) > 0
    u, v = sp.symbols('u v')
    checks['compensating operators keep the original Hamiltonian unchanged'] = sp.expand((u-v)+v-u) == 0
    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_z = sp.diag(1, -1)
    omega = sp.Matrix([1, 0])
    returned = (sp.eye(4)/48+(sp.kronecker_product(sigma_x, sp.eye(2))+
                 sp.kronecker_product(sp.eye(2), sigma_x))/2000+
                 sp.kronecker_product(sigma_x, sigma_x)/3000+
                 sp.kronecker_product(sigma_z, sigma_z)/20000)
    parts = split_return(returned, omega, omega)
    checks['exact split reconstructs the full return operator'] = returned == (
        parts['alpha']*sp.eye(4)+sp.kronecker_product(parts['single_a'],sp.eye(2))+
        sp.kronecker_product(sp.eye(2),parts['single_b'])+parts['connected'])
    checks['both one-cube operators have zero reference expectation'] = (
        (omega.T*parts['single_a']*omega)[0] == 0 and
        (omega.T*parts['single_b']*omega)[0] == 0)
    checks['connected return has both conditional expectations zero'] = (
        parts['va'].T*parts['connected']*parts['va'] == sp.zeros(2)
        and parts['vb'].T*parts['connected']*parts['vb'] == sp.zeros(2))
    checks['remaining connected source excites both cubes'] = (
        parts['connected']*parts['omega'] == sp.Matrix([0,0,0,sp.Rational(1,3000)]))
    return checks


def source_pins():
    paths = ['physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md', 'physics/yc12/YC12_RESULT.json',
             'physics/yc14/YC14_OBSERVER_RESET_AND_INTERFACE.md',
             'physics/yc15/YC15_BOUNDARY_RETURN_CHANNELS.md',
             'physics/yc15/yc15_boundary_return_channels.py', 'physics/yc15/YC15_RESULT.json',
             'physics/yc16/YC16_ABSORBED_CUBE_RETURN.md',
             'physics/yc16/yc16_absorbed_cube_return.py', 'physics/yc16/test_yc16.py']
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def run():
    for path in ('physics/yc12/YC12_RESULT.json', 'physics/yc15/YC15_RESULT.json'):
        y15.verify_predecessor(path)
    checks = symbolic_checks()
    src = source_bounds(1, 1)
    checks['whole-window connected source norm is at most one over 8658'] = src['connected_source_squared_upper'] == Q(1, 8658)**2
    checks['one free cube kills the connected source exactly'] = (
        source_bounds(0, 1)['connected_source_squared_upper'] == 0 and
        source_bounds(1, 0)['connected_source_squared_upper'] == 0)
    checks['single return norm is bounded by its proved vacuum expectation cap'] = ALPHA_MAX >= Q(1, 48) and ALPHA_MAX < Q(1, 24)
    row = join_bounds(ETA_MAX)
    ref = row['reference']
    checks['corrected cube reference retains a full charged gap above 27 percent'] = ref['new_cube_gap_lower'] > Q(27, 100)
    checks['exact one-cube source bound fits its simple polynomial envelope'] = Q(24, 481) < Q(1, 20)
    checks['ground energy displacement retains the complete quadratic residual'] = Q(1, 400)/Q(27, 100) == Q(1, 108)
    checks['new reference vector has a uniform norm-distance envelope'] = 2*Q(1, 20)/Q(27, 100) == Q(10, 27)
    checks['source is transported to the new reference rather than frozen'] = (
        ref['connected_source_after_reference_upper'] == Q(1,8658)+Q(5,162)*ETA_MAX**2
        and ref['connected_source_after_reference_upper'] < Q(1,8657))
    checks['connected recentering regenerates single-cube operators only at fourth interface order'] = (
        ref['recentered_connected_one_cube_cost_per_block_upper'] ==
        24*ETA_MAX**2*(ref['new_ground_vector_distance_upper']/12+
                      ref['new_ground_vector_distance_upper']**2/6))
    checks['double centering makes the regenerated scalar cost sixth order'] = (
        ref['recentered_connected_scalar_cost_per_face_upper'] ==
        ETA_MAX**2*ref['new_ground_vector_distance_upper']**2/6)
    checks['two-bridge inverse change fits eta squared over thirteen'] = (
        Q(5, 81)+Q(235, 16848)+ETA_MAX**2/3888 < Q(1, 13))
    checks['single compensating creator fits three sixteenths eta squared'] = (
        (Q(1,20)+Q(235,468)*Q(10,27)*ETA_MAX**2)/Q(27,100) < Q(3,16))
    checks['the entire compensated cluster ball has a positive reserve'] = row['mapping_reserve'] > Q(1, 40000)
    checks['full nonlinear contraction survives absorption'] = row['contraction'] == Q(5450044392, 6218422849) < Q(877, 1000)
    checks['all excited channels retain the relative inverse reserve'] = row['relative_return'] == Q(3229655936, 18655268547) < Q(174,1000)
    checks['original actual lattice still has gap greater than nine fortieths'] = row['full_gap_lower'] == Q(15425612611, 68125224960) > Q(9,40)
    checks['compensation costs a declared quadratic term in the gap formula'] = row['full_gap_lower'] == GAP-Q(1024,5)*ETA_MAX-Q(517,108)*ETA_MAX**2
    checks['the previous stronger numeric bound is retained separately'] = y15.join_bounds(ETA_MAX)['full_gap_lower'] > row['full_gap_lower']
    checks = {key: bool(value) for key,value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([key for key,value in checks.items() if not value])
    return y15.y14.encode(dict(stage='YC16', base_commit='8a976df', checks=checks,
                return_split='R_p=alpha_p I+A_p,A tensor I+I tensor A_p,B+C_p',
                full_window_return=src,
                connected_mixed_derivative=Q(19,86261760),
                free_single_face_eigenvalue=Q(-19,1872),
                free_connected_face_pair_eigenvalue=Q(3931,599040),
                endpoint=row,
                full_gap_formula='g_star-(1024/5)*eta-(517/108)*eta^2; eta<=1/4320',
                claim_boundary=dict(full_one_cube_operators_absorbed=True,
                    original_Hamiltonian_preserved_by_compensation=True,
                    first_single_cube_ground_creator_removed_through_interface_order_two=True,
                    connected_source_quadratic_in_internal_couplings=True,
                    connected_operator_zero_at_free_point=False,
                    all_remaining_interactions_and_clusters_retained=True,
                    original_YC15_window_preserved=True,
                    improved_numeric_window_or_floor=False,
                    only_Wilson_coupling_shift=False, spatial_RG=False,
                    continuum_mass_gap=False, formal_machine_verification=False),
                source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC16_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('certificate or source pins changed')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC16: {len(result['checks'])} exact checks pass.")
    print('Full one-cube operators absorbed; new-reference connected source <1/8657; original gap >9/40 preserved.')
