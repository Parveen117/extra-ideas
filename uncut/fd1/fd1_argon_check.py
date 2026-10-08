"""FD1 data check: the frame-defect number D of argon on the isochore 1 mol/l, from tabulated Cv, Cp, Joule-Thomson
(NIST Chemistry WebBook fluid tables, retrieved 8 Oct 2026), against the van der Waals form with argon's a, b."""
import csv, json, os
here = os.path.dirname(os.path.abspath(__file__))
rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(here, 'data', 'argon_isochore_1mol_per_l.csv')))]
R = 8.314462618
def Hh(r):                      # C_P/(alpha T), with alpha T = 1 + mu_JT * C_P * rho  (mu C_P = V (alpha T - 1))
    aT = 1 + r['JT_K_per_MPa']*r['Cp_J_per_molK']*r['rho_mol_per_l']*1e-3
    return r['Cp_J_per_molK']/aT
a_vdw, b_vdw = 0.1355, 3.201e-5     # Pa m^6/mol^2, m^3/mol (argon)
V = 1e-3                             # m^3/mol
out = []
for i in range(1, len(rows) - 1):
    lo, r, hi = rows[i-1], rows[i], rows[i+1]
    D = (r['T_K']/r['Cv_J_per_molK'])*(Hh(hi) - Hh(lo))/(hi['T_K'] - lo['T_K'])
    Dv = 2*a_vdw*(V - b_vdw)/(R*r['T_K']*V**2)
    out.append({'T': r['T_K'], 'D_data': round(D, 5), 'D_vdW': round(Dv, 5), 'ratio': round(D/Dv, 3)})
res = {'rows': out, 'all_positive': all(o['D_data'] > 0 for o in out),
       'falls_with_T': all(out[i]['D_data'] > out[i+1]['D_data'] for i in range(len(out)-1)),
       'ideal_value': 0}
json.dump(res, open(os.path.join(here, 'FD1_ARGON_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__':
    for o in out: print(o)
    print(res['all_positive'], res['falls_with_T'])
