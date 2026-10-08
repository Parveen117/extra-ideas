"""TD1: the thermodynamic response block read by the diagonal observer.  C_V/C_P is F = R - D; a gas of light has seen = lost.

NT-1: H = [[a, b], [b, c]] = second derivatives of U(S, V); Delta = ac - b^2; C_V = T/a, C_P = Tc/Delta, K_S = Vc, K_T = V Delta/a.
DO1: the diagonal observer reads every cut alike; for him F = R - D is the determinant.
Symbolic (sympy); fluid numbers from reference equations of state (CoolProp).  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def memory(U, S, V):
    a, b, c = sp.diff(U, S, 2), sp.diff(U, S, V), sp.diff(U, V, 2)
    return sp.simplify(b**2/(a*c))


def run():
    out, num = {}, {}
    a, b, c, T, V = sp.symbols('a b c T V', positive=True)
    Delta = a*c - b**2
    CV, CP, KS, KT = T/a, T*c/Delta, V*c, V*Delta/a                   # NT-1
    m = b**2/(a*c)
    # T1: the pure number of the block: C_V/C_P = K_T/K_S = 1 - m , m = b^2/(ac) unchanged by any choice of units for S and V
    out['T1_ratio'] = z(CV/CP - (1 - m)) and z(KT/KS - (1 - m))
    ls, lv = sp.symbols('l_s l_v', positive=True)
    out['T1_unit_free'] = z((b*ls*lv)**2/((a*ls**2)*(c*lv**2)) - m)
    # in units in which the two cuts read alike the block is sqrt(ac) [[1, t], [t, 1]], t^2 = m: the diagonal observer; F = R - D = 1 - m
    t = sp.Symbol('t', positive=True)
    Hd = sp.Matrix([[1, t], [t, 1]])
    out['T1_F_is_the_ratio'] = z(Hd.det() - (1 - t**2)) and z(Hd[0, 0]**2 - Hd[0, 1]**2 - (1 - t**2))
    # its own rates are 1 +- t: ratio e^l with tanh^2(l/2) = m  (RMG1's l; NC1's eta = l/2)
    ell = sp.log((1 + t)/(1 - t))
    out['T1_anisotropy'] = z(sp.tanh(ell/2).rewrite(sp.exp) - t)
    # T2: matter gas: counts fixed at fixed S, momentum ~ 1/L, energy ~ momentum^2: U = g(S) V^(-2/d); equal share: T ~ U
    S, d, k = sp.symbols('S d k', positive=True)
    U_m = sp.exp(2*S/(d*k))*V**(-2/d)
    out['T2_matter_gas'] = z(memory(U_m, S, V) - 2/(d + 2))
    # f readings in all (turning, etc.): U = exp(2S/(f k)) V^(-2/f): m = 2/(f+2), lost : seen = 2 : f
    f = sp.Symbol('f', positive=True)
    U_f = sp.exp(2*S/(f*k))*V**(-2/f)
    mf = memory(U_f, S, V)
    out['T2_counts'] = z(mf - 2/(f + 2)) and z(mf/(1 - mf) - 2/f)
    out['T2_two_cuts_is_the_diagonal'] = z(memory(U_m, S, V).subs(d, 2) - sp.Rational(1, 2))
    # T3: gas of light: energy ~ momentum, no count kept: U = k' S^((d+1)/d) V^(-1/d): m = 1 for every d: seen = lost, F = 0
    U_l = S**((d + 1)/d)*V**(-1/d)
    out['T3_light_gas_seen_equals_lost'] = z(memory(U_l, S, V) - 1)
    # because U is then homogeneous of degree one in (S, V): any such U has m = 1
    g = sp.Function('g')
    U_h = S*g(V/S)
    out['T3_degree_one_means_no_scale'] = z(memory(U_h, S, V) - 1)
    out['T3_pressure_and_entropy_numbers'] = z((-sp.diff(U_l, V)*V/U_l) - 1/d) and z((sp.diff(U_l, S)*S/U_l) - (d + 1)/d)
    # T4: van der Waals with f readings, on its critical isochore: gamma = 1 + (2/f)/(1 - Tc/T): the diagonal (gamma = 2) at T/Tc = f/(f-2)
    Tr = sp.Symbol('T_r', positive=True)
    gam = 1 + (2/f)/(1 - 1/Tr)
    out['T4_vdW_diagonal'] = sp.solve(sp.Eq(gam, 2), Tr) == [f/(f - 2)]
    # numbers from reference equations of state
    try:
        from CoolProp.CoolProp import PropsSI as P
        def mem(fl, T_, **kw):
            key, val = ('D', kw['D']) if 'D' in kw else ('P', kw['Pr'])
            return 1 - P('CVMASS', 'T', T_, key, val, fl)/P('CPMASS', 'T', T_, key, val, fl)
        low = {fl: round(mem(fl, 300, Pr=1000), 5) for fl in ('Helium', 'Neon', 'Argon', 'Krypton', 'Xenon', 'Nitrogen', 'Oxygen', 'CarbonMonoxide')}
        num['memory_at_low_density_300K'] = low
        diag = {}
        for fl in ('Neon', 'Argon', 'Krypton', 'Xenon', 'Nitrogen', 'Oxygen'):
            Tc, rc = P('Tcrit', fl), P('rhocrit', fl)
            lo, hi = Tc*1.0005, Tc*6
            fn = lambda T_: mem(fl, T_, D=rc) - 0.5
            for _ in range(60):
                mid = (lo + hi)/2
                if fn(lo)*fn(mid) <= 0: hi = mid
                else: lo = mid
            diag[fl] = round(lo/Tc, 4)
        num['diagonal_on_critical_isochore_T_over_Tc'] = diag
        num['near_critical_memory_at_1p001_Tc'] = {fl: round(mem(fl, P('Tcrit', fl)*1.001, D=P('rhocrit', fl)), 4) for fl in ('Argon', 'Krypton', 'Xenon', 'CarbonDioxide')}
        out['N1_noble_gases_two_fifths'] = all(abs(low[fl] - 0.4) < 2e-4 for fl in ('Helium', 'Neon', 'Argon', 'Krypton', 'Xenon'))
        out['N1_two_atom_gases_two_sevenths'] = all(abs(low[fl] - 2/7) < 3e-3 for fl in ('Nitrogen', 'Oxygen', 'CarbonMonoxide'))
        out['N2_diagonal_noble_between_2p6_and_2p9'] = all(2.6 < diag[fl] < 2.9 for fl in ('Neon', 'Argon', 'Krypton', 'Xenon'))
        out['N3_memory_goes_to_one_at_the_critical_point'] = all(v > 0.99 for v in num['near_critical_memory_at_1p001_Tc'].values())
    except ImportError:
        num['fluids'] = 'CoolProp not installed: fluid numbers skipped'
    out = {k_: bool(v) for k_, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'TD1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o); print(n)
