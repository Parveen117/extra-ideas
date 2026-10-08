"""FD1 collapse check: the frame-defect number D at one reduced density (rho/rho_c ~ 0.0746) against reduced temperature,
for argon, krypton, xenon, nitrogen, and carbon dioxide. Tables: NIST Chemistry WebBook fluid data, retrieved 8 Oct 2026."""
import csv, json, os
here = os.path.dirname(os.path.abspath(__file__))
subs = {  # file, T_c (K), rho_c (mol/l)
 'argon': ('argon_isochore_1mol_per_l.csv', 150.687, 13.4074),
 'krypton': ('krypton_isochore_0p808.csv', 209.48, 10.85),
 'xenon': ('xenon_isochore_0p626.csv', 289.733, 8.4),
 'nitrogen': ('nitrogen_isochore_0p833.csv', 126.192, 11.1839),
 'carbon_dioxide': ('carbon_dioxide_isochore_0p791.csv', 304.128, 10.6249),
}
def load(f): return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(here, 'data', f)))]
def Hh(r): return r['Cp_J_per_molK']/(1 + r['JT_K_per_MPa']*r['Cp_J_per_molK']*r['rho_mol_per_l']*1e-3)
def curve(Tr, vr): return 2.25*(1 - 1/(3*vr))/(Tr*vr)        # FD1-T4, reduced
table = {}
for name, (f, Tc, rc) in subs.items():
    rows = load(f); pts = []
    for i in range(1, len(rows) - 1):
        lo, r, hi = rows[i-1], rows[i], rows[i+1]
        dT = hi['T_K'] - lo['T_K']
        D = (r['T_K']/r['Cv_J_per_molK'])*(Hh(hi) - Hh(lo))/dT
        Dint = (r['T_K']/r['Cv_J_per_molK'])*(hi['Cv_J_per_molK'] - lo['Cv_J_per_molK'])/dT   # part from a heat capacity that itself changes
        Tr = r['T_K']/Tc; vr = rc/r['rho_mol_per_l']
        pts.append({'Tr': round(Tr, 3), 'D': round(D, 4), 'D_minus_internal': round(D - Dint, 4), 'D_times_Tr': round(D*Tr, 4), 'curve': round(curve(Tr, vr), 4)})
    table[name] = pts
def interp(pts, Tr, key):
    for a, b in zip(pts, pts[1:]):
        if a['Tr'] <= Tr <= b['Tr']:
            w = (Tr - a['Tr'])/(b['Tr'] - a['Tr']); return a[key] + w*(b[key] - a[key])
grid = [1.4, 1.6, 1.8, 2.0, 2.2]
comp = {Tr: {n: (round(interp(p, Tr, 'D'), 4) if interp(p, Tr, 'D') is not None else None) for n, p in table.items()} for Tr in grid}
noble = ['argon', 'krypton', 'xenon']
spread = {Tr: round((max(comp[Tr][n] for n in noble) - min(comp[Tr][n] for n in noble))/(sum(comp[Tr][n] for n in noble)/3), 4) for Tr in grid}
n2_excess = {Tr: round(comp[Tr]['nitrogen']/(sum(comp[Tr][n] for n in noble)/3) - 1, 4) for Tr in grid}
res = {'points': table, 'D_at_common_Tr': comp, 'relative_spread_noble_three': spread, 'nitrogen_above_noble_mean': n2_excess,
       'noble_collapse_within_9_percent': all(s < 0.09 for s in spread.values()),
       'noble_collapse_within_5_percent_above_Tr_1p5': all(s < 0.05 for Tr, s in spread.items() if Tr > 1.5),
       'nitrogen_departs_15_to_33_percent': all(0.15 < v < 0.33 for v in n2_excess.values()),
       'carbon_dioxide_departs': all(abs(comp[Tr]['carbon_dioxide'] - comp[Tr]['argon'])/comp[Tr]['argon'] > 0.2 for Tr in grid if comp[Tr]['carbon_dioxide'] is not None)}
json.dump(res, open(os.path.join(here, 'FD1_COLLAPSE_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__':
    for Tr in grid: print(Tr, comp[Tr], 'spread', spread[Tr], 'curve', round(curve(Tr, 13.41), 4))
    print({n: [(p['Tr'], p['D_minus_internal']) for p in pts][::3] for n, pts in table.items() if n == 'carbon_dioxide'})
    print(spread, n2_excess); print(res['noble_collapse_within_9_percent'], res['noble_collapse_within_5_percent_above_Tr_1p5'], res['nitrogen_departs_15_to_33_percent'], res['carbon_dioxide_departs'])
