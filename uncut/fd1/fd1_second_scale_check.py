"""FD1 second-scale check. D' = (T/C_V) d[C_P/(alpha T) - C_V]/dT at fixed V: zero for an ideal gas with ANY heat capacity,
so internal motion is removed at the ideal level. Seven substances at rho/rho_c ~ 0.0746; excess over the monatomic
curve against the acentric factor (standard tables)."""
import csv, json, os
here = os.path.dirname(os.path.abspath(__file__))
subs = {  # file, T_c (K), rho_c (mol/l), acentric factor
 'argon': ('argon_isochore_1mol_per_l.csv', 150.687, 13.4074, 0.000),
 'krypton': ('krypton_isochore_0p808.csv', 209.48, 10.85, 0.000),
 'xenon': ('xenon_isochore_0p626.csv', 289.733, 8.4, 0.004),
 'oxygen': ('oxygen_isochore_1p017.csv', 154.581, 13.63, 0.022),
 'nitrogen': ('nitrogen_isochore_0p833.csv', 126.192, 11.1839, 0.037),
 'carbon_monoxide': ('carbon_monoxide_isochore_0p809.csv', 132.86, 10.85, 0.048),
 'carbon_dioxide': ('carbon_dioxide_isochore_0p791.csv', 304.128, 10.6249, 0.224),
}
def load(f): return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(here, 'data', f)))]
def G(r): return r['Cp_J_per_molK']/(1 + r['JT_K_per_MPa']*r['Cp_J_per_molK']*r['rho_mol_per_l']*1e-3) - r['Cv_J_per_molK']
pts = {}
for n, (f, Tc, rc, w) in subs.items():
    rows = load(f); pts[n] = []
    for i in range(1, len(rows) - 1):
        lo, r, hi = rows[i-1], rows[i], rows[i+1]
        pts[n].append((r['T_K']/Tc, (r['T_K']/r['Cv_J_per_molK'])*(G(hi) - G(lo))/(hi['T_K'] - lo['T_K']), r['Cv_J_per_molK']))
def at(n, Tr, k=1):
    p = pts[n]
    for a, b in zip(p, p[1:]):
        if a[0] <= Tr <= b[0]:
            t = (Tr - a[0])/(b[0] - a[0]); return a[k] + t*(b[k] - a[k])
grid = [1.5, 1.7, 1.9, 2.1]
tab = {Tr: {n: round(at(n, Tr), 4) for n in subs} for Tr in grid}
noble = ['argon', 'krypton', 'xenon']
rel = {Tr: {n: round(tab[Tr][n]/(sum(tab[Tr][m] for m in noble)/3) - 1, 4) for n in subs} for Tr in grid}
res = {'Dprime': tab, 'relative_to_monatomic_mean': rel, 'acentric': {n: subs[n][3] for n in subs},
       'noble_within_half_percent': all(abs(rel[Tr][n]) < 0.005 for Tr in grid for n in noble),
       'all_seven_within_7_percent': all(abs(v) < 0.07 for Tr in grid for v in rel[Tr].values()),
       'carbon_dioxide_within_3_percent': all(abs(rel[Tr]['carbon_dioxide']) < 0.03 for Tr in grid)}
json.dump(res, open(os.path.join(here, 'FD1_SECOND_SCALE_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__':
    for Tr in grid: print(Tr, tab[Tr]); print('   rel', rel[Tr])
    print(res['noble_within_half_percent'], res['all_seven_within_7_percent'], res['carbon_dioxide_within_3_percent'])
