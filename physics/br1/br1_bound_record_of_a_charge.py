"""BR1: the bound record of a charge around a charged centre.  Where the fine-structure number sits in a record, and the size it gives.

OB1-U6: the field changes the reading (n; r) by q(E.r ; nE + r x B).  PR2-T2: the boosted mass is (n; r) = gamma g (1; v).
OR1 / WQ1: a record is kept when the area of each pair is whole: L = l kappa, J_r = n_r kappa.
Units c = 1.  K = strength of the attraction between the two contents; alpha := K / kappa.
Symbolic (sympy) and mpmath.  Python 3.12."""
import json, os
import sympy as sp
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    g, K, L, E, u, ph = sp.symbols('g K L E u phi', positive=True)
    # B1: central attraction dp/dt = -K rhat / r^2 with p = gamma g v.  Conserved: E = gamma g - K/r and L = r x p.
    #     With u = 1/r:  gamma g = E + K u ,  p_r = -L du/dphi ,  and (gamma g)^2 = g^2 + p_r^2 + L^2 u^2 :
    U = sp.Function('U')(ph)
    shell = (E + K*U)**2 - g**2 - L**2*(sp.diff(U, ph)**2 + U**2)
    orbit = sp.diff(U, ph, 2) + (1 - K**2/L**2)*U - E*K/L**2
    dshell = sp.diff(shell, ph)
    out['B1_orbit_equation'] = z(dshell + 2*L**2*sp.diff(U, ph)*orbit)        # d(shell)/dphi = -2 L^2 U' x (orbit equation)
    # every bound history: in-out rate / round rate = sqrt(1 - K^2/L^2), exactly, at every amplitude
    w = sp.sqrt(1 - K**2/L**2)
    A, B = sp.symbols('A B', real=True)
    sol = E*K/(L**2 - K**2) + B*sp.cos(w*ph)
    out['B1_ratio_exact'] = z(orbit.subs(U, sol).doit())
    # B2: circles (B = 0): u0 = E K/(L^2 - K^2).  Speed = p/(gamma g) = L u0/(E + K u0) = K/L ; shell gives E = g sqrt(1 - K^2/L^2)
    u0 = E*K/(L**2 - K**2)
    out['B2_speed_is_K_over_L'] = z(L*u0/(E + K*u0) - K/L)
    Ec = g*sp.sqrt(1 - K**2/L**2)
    out['B2_energy_is_count_factor'] = z(((E + K*u0)**2 - g**2 - L**2*u0**2).subs(E, Ec))
    # so with L = l kappa:  speed = alpha/l ,  memory = alpha^2/l^2 ,  energy = g x count factor sqrt(1 - alpha^2/l^2)
    al, l, kap = sp.symbols('alpha l kappa', positive=True)
    r0 = (1/u0).subs(E, Ec).subs({K: al*kap, L: l*kap})
    out['B2_size'] = z(r0 - l**2*kap/(g*al)*sp.sqrt(1 - al**2/l**2))
    # B3: three lengths of one charge in ratio alpha: K/g (charge), kappa/g (turn; twice GM1's displacement a), kappa/(g alpha) (record)
    charge_len, turn_len, size = al*kap/g, kap/g, kap/(g*al)
    out['B3_lengths_in_ratio_alpha'] = z(charge_len/turn_len - al) and z(turn_len/size - al) and z(size - 2*(kap/(2*g))/al)
    # binding of the first circle: 1 - sqrt(1 - alpha^2) = alpha^2/2 + alpha^4/8 + ...
    ser = sp.series(1 - sp.sqrt(1 - al**2), al, 0, 6).removeO()
    out['B3_binding'] = z(ser - (al**2/2 + al**4/8))
    # B4: the in-out pair: J_r = 2 pi [ E K / sqrt(g^2 - E^2) - sqrt(L^2 - K^2) ]  (checked by quadrature below); whole: J_r = 2 pi n_r kappa
    nr = sp.Symbol('n_r', positive=True)
    Elev = g/sp.sqrt(1 + al**2/(nr + sp.sqrt(l**2 - al**2))**2)
    Jr = 2*sp.pi*(E*K/sp.sqrt(g**2 - E**2) - sp.sqrt(L**2 - K**2))
    out['B4_levels'] = z((Jr.subs({K: al*kap, L: l*kap}).subs(E, Elev) - 2*sp.pi*nr*kap)) or \
        abs(sp.N((Jr.subs({K: al*kap, L: l*kap}).subs(E, Elev) - 2*sp.pi*nr*kap).subs({al: sp.Rational(1, 7), l: 2, nr: 3, kap: 1, g: 1}))) < 1e-12
    out['B4_circles_are_n_r_zero'] = z(Elev.subs(nr, 0)**2 - (g**2*(1 - al**2/l**2)))
    # slopes of J_r: with respect to L it is -2 pi (round rate / in-out rate)
    out['B4_slope_in_L'] = z(sp.diff(Jr, L) + 2*sp.pi*L/sp.sqrt(L**2 - K**2))
    # quadrature of p_r dr between the turning points at sample values
    mp.mp.dps = 30
    gv, Kv, Lv, Ev = mp.mpf(1), mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('0.97')
    f = lambda r: (Ev + Kv/r)**2 - gv**2 - Lv**2/r**2
    a_, b_, c_ = Ev**2 - gv**2, 2*Ev*Kv, Kv**2 - Lv**2               # a r^2 + b r + c = 0
    r1 = (-b_ + mp.sqrt(b_**2 - 4*a_*c_))/(2*a_); r2 = (-b_ - mp.sqrt(b_**2 - 4*a_*c_))/(2*a_)
    quad = 2*mp.quad(lambda r: mp.sqrt(max(f(r), 0)), mp.linspace(min(r1, r2), max(r1, r2), 9))
    closed = 2*mp.pi*(Ev*Kv/mp.sqrt(gv**2 - Ev**2) - mp.sqrt(Lv**2 - Kv**2))
    out['B4_quadrature'] = abs(quad - closed) < mp.mpf(10)**(-10)
    # B5: second shell: (n_r, l) = (1, 1) and (0, 2): difference over the binding scale g alpha^2/2 is alpha^2/16 at leading order
    d2 = sp.series((Elev.subs({nr: 0, l: 2}) - Elev.subs({nr: 1, l: 1}))/(g*al**2/2), al, 0, 4).removeO()
    out['B5_split_of_second_shell'] = z(d2 - al**2/16)
    # B6: the turn left over per circuit: 2 pi (1/sqrt(1 - alpha^2/l^2) - 1) = pi alpha^2/l^2 + ...
    left = sp.series(2*sp.pi*(1/sp.sqrt(1 - al**2/l**2) - 1), al, 0, 4).removeO()
    out['B6_left_over_turn'] = z(left - sp.pi*al**2/l**2)
    # against the field of a mass at the same K/L: MO1-M5 ratio^2 = 1 - 3x with x = 2 (K/L)^2 at first order: six times the charge's
    x = sp.Symbol('x', positive=True)
    L2 = 1/(x*(2 - 3*x))                                             # MO1: L^2 = r_s r/(2 - 3x), r_s = 1, K = g r_s/2 -> (K/L)^2 = x(2-3x)/4 per g^2
    kl2 = x*(2 - 3*x)/4
    out['B6_mass_field_is_six_times'] = z(sp.limit((3*x)/kl2, x, 0) - 6)
    # numbers (constants as in FS1)
    hbar, c, me, e_, k_e = 1.054571817e-34, 299792458.0, 9.1093837015e-31, 1.602176634e-19, 8.9875517923e9
    alpha = k_e*e_**2/(hbar*c)
    a_disp = hbar/(2*me*c)
    num['alpha'] = alpha
    num['speed_of_first_circle_over_c'] = alpha
    num['memory_of_first_circle'] = alpha**2
    num['binding_over_rest_count'] = 1 - (1 - alpha**2)**0.5
    num['size_m'] = 2*a_disp/alpha*(1 - alpha**2)**0.5
    num['lengths_m'] = dict(charge=alpha*2*a_disp, turn=2*a_disp, record=2*a_disp/alpha)
    num['left_over_turn_first_circle_rad'] = float(2*mp.pi*(1/mp.sqrt(1 - alpha**2) - 1))
    num['split_of_second_shell_over_binding_scale'] = dict(line=alpha**2/16, measured_recalled=10969.0e6/3.2898419603e15)
    num['binding_scale_eV'] = 0.5*alpha**2*me*c**2/e_
    # ratios of the rates between circles (leading order): counts only
    rate = lambda m_, n_: sp.Rational(1, m_**2) - sp.Rational(1, n_**2)
    num['rate_ratios'] = {'(1,2):(2,3)': str(rate(1, 2)/rate(2, 3)), '(2,4):(2,3)': str(rate(2, 4)/rate(2, 3)), '(2,5):(2,4)': str(rate(2, 5)/rate(2, 4))}
    out['N1_rate_ratios_are_counts'] = rate(1, 2)/rate(2, 3) == sp.Rational(27, 5) and rate(2, 4)/rate(2, 3) == sp.Rational(27, 20)
    air = dict(Ha=656.28, Hb=486.13, Hg=434.05)                      # nm, recalled
    num['measured_wavelength_ratios_recalled'] = {'Ha/Hb': round(air['Ha']/air['Hb'], 5), 'Hb/Hg': round(air['Hb']/air['Hg'], 5)}
    out['N2_wavelength_ratios'] = abs(air['Ha']/air['Hb'] - 27/20) < 1e-4 and abs(air['Hb']/air['Hg'] - 28/25) < 1e-4
    out['N3_record_length_is_the_known_size'] = abs(num['lengths_m']['record'] - 5.29177e-11) < 1e-15 and num['size_m'] < num['lengths_m']['record']
    out['N4_split_within_half_percent'] = abs(num['split_of_second_shell_over_binding_scale']['measured_recalled']/(alpha**2/16) - 1) < 5e-3
    out = {k_: bool(v) for k_, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'BR1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o)
    for kk, v in n.items(): print(kk, v)
