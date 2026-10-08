"""FD1 step check: D of argon at 300 K from central differences with steps 20, 10, 5 K, a Richardson value,
and the effect of the printed digits (rounding noise) at each step."""
import csv, json, os
here = os.path.dirname(os.path.abspath(__file__))
rows = {float(r['T_K']): {k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(here, 'data', 'argon_isochore_1mol_per_l_fine.csv')))}
def Hh(r, dcp=0.0, dmu=0.0):
    cp, mu = r['Cp_J_per_molK'] + dcp, r['JT_K_per_MPa'] + dmu
    return cp/(1 + mu*cp*r['rho_mol_per_l']*1e-3)
def D(h, s=(0, 0, 0, 0)):
    lo, hi, mid = rows[300 - h], rows[300 + h], rows[300]
    return (300/mid['Cv_J_per_molK'])*(Hh(hi, s[0], s[1]) - Hh(lo, s[2], s[3]))/(2*h)
out = {f'D_step_{2*h}K': round(D(h), 6) for h in (20, 10, 5)}
out['richardson_10_5'] = round((4*D(5) - D(10))/3, 6)
# rounding: last printed digit of Cp is 1e-4 (+-5e-5), of JT 1e-5 (+-5e-6): worst case on the difference
for h in (20, 10, 5):
    worst = max(abs(D(h, (a*5e-5, b*5e-6, c*5e-5, d*5e-6)) - D(h)) for a in (-1, 1) for b in (-1, 1) for c in (-1, 1) for d in (-1, 1))
    out[f'rounding_bound_step_{2*h}K'] = round(worst, 6)
out['truncation_20K_vs_richardson'] = round(abs(D(10) - out['richardson_10_5']), 6)
out['stable_to_1_percent'] = abs(D(10) - out['richardson_10_5'])/out['richardson_10_5'] < 0.01
json.dump(out, open(os.path.join(here, 'FD1_STEP_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
