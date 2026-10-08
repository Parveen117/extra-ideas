"""DL1: the diagonal line of a fluid: where lost = seen (C_P = 2 C_V) on the diagram, for the noble fluids.

TD1: m = b^2/(ac) = (C_P - C_V)/C_P is the lost part for the diagonal observer; m = 1/2 is the diagonal.
Symbolic (sympy); fluid numbers from reference equations of state (CoolProp).  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    a, b, c = sp.symbols('a b c', positive=True)
    m = b**2/(a*c)
    # E1: the same line in five readings
    Delta = a*c - b**2
    T, V = sp.symbols('T V', positive=True)
    CV, CP, KS, KT = T/a, T*c/Delta, V*c, V*Delta/a
    half = {b: sp.sqrt(a*c/2)}
    out['E1_heat_and_stiffness'] = z((CP/CV).subs(half) - 2) and z((KS/KT).subs(half) - 2)
    # in the ensemble of QC2 the offsets (ds, dv) have covariance kappa H^-1: squared correlation = m ; same for (dT, dP) with kappa H
    H = sp.Matrix([[a, b], [b, c]]); Hi = H.inv()
    out['E1_correlation_of_fluctuations'] = z(Hi[0, 1]**2/(Hi[0, 0]*Hi[1, 1]) - m) and z(H[0, 1]**2/(H[0, 0]*H[1, 1]) - m)
    # own rates of the block on the diagonal: ratio (1 + sqrt m)/(1 - sqrt m) = 3 + 2 sqrt 2 at m = 1/2
    t = sp.sqrt(sp.Rational(1, 2))
    out['E1_rate_ratio'] = z((1 + t)/(1 - t) - (3 + 2*sp.sqrt(2)))
    # E2: van der Waals with f readings: C_P - C_V = Nk/(1 - X), X = Ts(v)/T with Ts the spinodal: diagonal T = Ts f/(f-2)
    f, v, Tr = sp.symbols('f v T_r', positive=True)
    Ts = (3*v - 1)**2/(4*v**3)                                       # spinodal, reduced
    Pr = 8*Tr/(3*v - 1) - 3/v**2
    out['E2_spinodal'] = z(sp.solve(sp.Eq(sp.diff(Pr, v), 0), Tr)[0] - Ts)
    gam = 1 + (2/f)/(1 - Ts/Tr)
    Td = sp.solve(sp.Eq(gam, 2), Tr)[0]
    out['E2_diagonal_is_spinodal_times_count'] = z(Td - Ts*f/(f - 2))
    out['E2_peak_at_critical_density'] = z(sp.diff(Ts, v).subs(v, 1)) and sp.diff(Ts, v, 2).subs(v, 1) < 0 and z(Td.subs(v, 1) - f/(f - 2))
    # numbers
    try:
        from CoolProp.CoolProp import PropsSI as P
        import numpy as np
        fluids = ('Neon', 'Argon', 'Krypton', 'Xenon')
        def mem(fl, T_, D): return 1 - P('CVMASS', 'T', T_, 'D', D, fl)/P('CPMASS', 'T', T_, 'D', D, fl)
        def bis(fn, lo, hi):
            if fn(lo)*fn(hi) > 0: return None
            for _ in range(60):
                mid = (lo + hi)/2
                if fn(lo)*fn(mid) <= 0: hi = mid
                else: lo = mid
            return lo
        def tsat(fl, D, Tc, rc):
            q = 1 if D < rc else 0
            return bis(lambda T_: P('D', 'T', T_, 'Q', q, fl) - D, P('Ttriple', fl)*1.001, Tc*0.99995)
        def diagT(fl, rr):
            Tc, rc = P('Tcrit', fl), P('rhocrit', fl)
            Ts_ = tsat(fl, rr*rc, Tc, rc)
            r_ = bis(lambda T_: mem(fl, T_, rr*rc) - 0.5, (Ts_ if Ts_ else Tc)*1.0003, Tc*8)
            return r_/Tc if r_ else None
        grid = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
        arch = {fl: [round(diagT(fl, rr), 4) for rr in grid] for fl in fluids}
        arch['van_der_Waals_f3'] = [round(float(3*(3/rr - 1)**2/(4/rr**3)), 4) for rr in grid]
        num['reduced_density_grid'] = grid
        num['diagonal_line_T_over_Tc'] = arch
        spread = lambda i, fls: (max(arch[fl][i] for fl in fls) - min(arch[fl][i] for fl in fls))/np.mean([arch[fl][i] for fl in fls])
        num['spread_four_fluids_gas_side'] = round(max(spread(i, fluids) for i in range(4)), 4)
        num['spread_krypton_xenon_whole_line'] = round(max(spread(i, ('Krypton', 'Xenon')) for i in range(len(grid))), 4)
        peaks, ends, triple, boyle = {}, {}, {}, {}
        for fl in fluids:
            Tc, rc = P('Tcrit', fl), P('rhocrit', fl)
            pk = max((diagT(fl, rr) or 0, round(float(rr), 2)) for rr in np.arange(0.80, 1.11, 0.02))
            peaks[fl] = dict(T_over_Tc=round(pk[0], 4), reduced_density=pk[1])
            Tv = bis(lambda T_: (1 - P('CVMASS', 'T', T_, 'Q', 1, fl)/P('CPMASS', 'T', T_, 'Q', 1, fl)) - 0.5, P('Ttriple', fl)*1.01, Tc*0.999)
            ends[fl] = dict(T_over_Tc=round(Tv/Tc, 4), reduced_density=round(P('D', 'T', Tv, 'Q', 1, fl)/rc, 4))
            Tt = P('Ttriple', fl)*1.0005
            triple[fl] = round(1 - P('CVMASS', 'T', Tt, 'Q', 0, fl)/P('CPMASS', 'T', Tt, 'Q', 0, fl), 4)
            boyle[fl] = round(bis(lambda T_: P('Bvirial', 'T', T_, 'D', 1e-6, fl), Tc*1.5, Tc*5)/Tc, 4)
        num['peak'] = peaks; num['vapour_end_on_coexistence'] = ends; num['lost_part_of_liquid_at_triple_point'] = triple
        num['Boyle_T_over_Tc'] = boyle
        num['peak_over_Boyle'] = {fl: round(peaks[fl]['T_over_Tc']/boyle[fl], 4) for fl in fluids}
        out['N1_gas_side_one_curve'] = num['spread_four_fluids_gas_side'] < 0.012
        out['N2_peak_near_critical_density'] = all(0.85 <= peaks[fl]['reduced_density'] <= 1.0 for fl in fluids)
        out['N3_triple_point_liquid_near_the_diagonal'] = all(abs(triple[fl] - 0.5) < 0.015 for fl in fluids)
        out['N4_vapour_end_near_three_quarters'] = all(abs(ends[fl]['T_over_Tc'] - 0.75) < 0.02 for fl in fluids)
        out['N5_count_form_is_ten_percent_high'] = all(0.88 < peaks[fl]['T_over_Tc']/3 < 0.95 for fl in fluids)
    except ImportError:
        num['fluids'] = 'CoolProp not installed: fluid numbers skipped'
    out = {k: bool(v_) for k, v_ in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'DL1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o)
    for k, v_ in n.items(): print(k, v_)
