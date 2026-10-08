"""CR1: chromium at its Neel transition, in the line's quantities. Inputs are published numbers (see the note);
everything else is arithmetic on them plus one exact identity (sympy)."""
import json, os, sympy as sp
out = {}
# ---- exact part: a first-order step of width d smeared uniformly: anomalous C and anomalous alpha_V
L, dV, V, d, TN, dTdP, Cb, ab = sp.symbols('L DeltaV V delta T_N dTdP C_bg alpha_bg', real=True)
C_an, a_an = L/d, dV/(V*d)
out['T1_ratio_width_free'] = sp.simplify(C_an/a_an - L*V/dV) == 0
cc = {dV: L*dTdP/TN}                                   # the step's own relation between L, DeltaV and dT_N/dP
out['T1_plateau'] = sp.simplify((C_an/(2*a_an)).subs(cc) - V*TN/(2*dTdP)) == 0      # w inside = V T_N / (2 dT_N/dP)
# alpha changes sign inside the window iff the width is below a critical width
dstar = sp.solve(sp.Eq(ab + dV/(V*d), 0), d)[0]
out['T2_critical_width'] = sp.simplify(dstar + dV/(V*ab)) == 0
# ---- numbers
inp = {'T_N_K': 311.4, 'latent_heat_J_per_mol_up': 1.10, 'latent_heat_J_per_mol_down': 0.97,
       'dTN_dP_K_per_Pa': -5.1e-8, 'molar_volume_m3_per_mol': 7.23e-6,
       'alpha_volume_bg_per_K': 1.47e-5, 'C_P_bg_J_per_molK': 23.3}
Lm = (inp['latent_heat_J_per_mol_up'] + inp['latent_heat_J_per_mol_down'])/2
DV = Lm*inp['dTN_dP_K_per_Pa']/inp['T_N_K']
rel = DV/inp['molar_volume_m3_per_mol']
num = {'mean_latent_heat': round(Lm, 3), 'DeltaV_m3_per_mol': DV, 'DeltaV_over_V': rel,
       'critical_width_K': round(-rel/inp['alpha_volume_bg_per_K'], 3),
       'w_background_J_per_mol': round(inp['C_P_bg_J_per_molK']/(2*inp['alpha_volume_bg_per_K'])),
       'w_plateau_J_per_mol': round(inp['molar_volume_m3_per_mol']*inp['T_N_K']/(2*inp['dTN_dP_K_per_Pa'])),
       'loop_heat_J_per_mol': round(inp['latent_heat_J_per_mol_up'] - inp['latent_heat_J_per_mol_down'], 3),
       'loop_count_J_per_molK': (inp['latent_heat_J_per_mol_up'] - inp['latent_heat_J_per_mol_down'])/inp['T_N_K']}
# profile of w for three widths
prof = {}
for width in (0.3, 1.0, 3.0):
    a_in = inp['alpha_volume_bg_per_K'] + rel/width
    c_in = inp['C_P_bg_J_per_molK'] + Lm/width
    prof[str(width)] = {'alpha_inside_per_K': a_in, 'C_inside': round(c_in, 3), 'w_inside_J_per_mol': round(c_in/(2*a_in)), 'sign_flips': a_in < 0}
out['T3_contracts_on_heating'] = rel < 0
out['T3_narrow_flips_wide_does_not'] = prof['0.3']['sign_flips'] and prof['1.0']['sign_flips'] and not prof['3.0']['sign_flips']
out['T3_plateau_negative'] = num['w_plateau_J_per_mol'] < 0 < num['w_background_J_per_mol']
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump({'checks': out, 'inputs': inp, 'numbers': num, 'profiles': prof}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CR1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out); print(num); print(prof)
