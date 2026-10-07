"""PH1 table 2: one reduced cycle, several fluids (reference EOS vs van der Waals)."""
import json
import ph1_real_fluid_holonomy as p

REDUCED = dict(T=(1.03, 1.30), rho=(0.40, 1.50))
FLUIDS = ['Argon', 'Krypton', 'Xenon', 'Methane', 'Nitrogen', 'CarbonDioxide', 'Water']


def run(n):
    rows = []
    for name in FLUIDS:
        H, crit, _ = p.real_fluid(name)
        T, rho = p.loop(n, REDUCED['T'][0]*crit['Tc'], REDUCED['T'][1]*crit['Tc'],
                        REDUCED['rho'][0]*crit['rhoc'], REDUCED['rho'][1]*crit['rhoc'])
        real = p.holonomy([H(float(t), float(r)) for t, r in zip(T, rho)])
        # van der Waals has rho_c = 1/(3b); use its own critical density for the same reduced cycle
        R, Tc, Pc = crit['R'], crit['Tc'], crit['Pc']
        rc_vdw = 1/(3*R*Tc/(8*Pc))
        Tv, rv = p.loop(n, REDUCED['T'][0]*Tc, REDUCED['T'][1]*Tc, REDUCED['rho'][0]*rc_vdw, REDUCED['rho'][1]*rc_vdw)
        Hv = p.model(name, 'vdw', 'cv0(T)')
        vdw = p.holonomy([Hv(float(t), float(r)) for t, r in zip(Tv, rv)])
        Hi = p.model(name, 'ideal', 'cv0(T)')
        ideal = p.holonomy([Hi(float(t), float(r)) for t, r in zip(T, rho)])
        rows.append(dict(fluid=name, theta_reference=real['theta_transport'], theta_loop=real['theta_loop'],
                         length_reference=real['shape_length'], theta_vdw=vdw['theta_transport'],
                         theta_ideal_cv0=ideal['theta_transport']))
    return rows


if __name__ == '__main__':
    import sys
    rows = run(int(sys.argv[1]) if len(sys.argv) > 1 else 4000)
    json.dump(dict(reduced_cycle=REDUCED, rows=rows), open('PH1_CORRESPONDING_STATES.json', 'w'), indent=1)
    for r in rows:
        print('%-14s Theta_ref=%+.5f (loop %+.5f)  L=%.3f  Theta_vdW=%+.5f  Theta_ideal=%+.6f  |Theta|/(L/2)=%.3f' % (
            r['fluid'], r['theta_reference'], r['theta_loop'], r['length_reference'], r['theta_vdw'],
            r['theta_ideal_cv0'], abs(r['theta_reference'])/(r['length_reference']/2)))
