"""PH1: response holonomy of a real fluid around a declared supercritical cycle.

H = D^2 U(s, v) per unit mass from measured-type response functions:
  a = T/c_v,  b = -(T/c_v)(dP/dT)_v,  c = rho^2 w^2   (w = speed of sound).
Routes: (A) loop formula of RMG1-T3, (B) direct transport dw = -(1/2) H^-1 dH w.
Reference equation of state through CoolProp (not part of the stdlib CI).
"""
import json
import math
import sys

import numpy as np
from CoolProp import CoolProp as CP


def loop(n, T0, T1, r0, r1):
    t = (np.arange(n)+0.5)*2*math.pi/n
    T = (T0+T1)/2+(T1-T0)/2*np.cos(t)
    rho = (r0+r1)/2+(r1-r0)/2*np.sin(t)
    return T, rho


def real_fluid(name):
    st = CP.AbstractState('HEOS', name)

    def H(T, rho):
        st.update(CP.DmassT_INPUTS, rho, T)
        cv, w = st.cvmass(), st.speed_sound()
        pT = st.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
        return np.array([[T/cv, -T*pT/cv], [-T*pT/cv, rho*rho*w*w]])
    R = st.gas_constant()/st.molar_mass()
    return H, dict(R=R, Tc=st.T_critical(), Pc=st.p_critical(), rhoc=st.rhomass_critical()), st


def model(name, kind, cv_mode):
    """kind: 'vdw' or 'ideal'; cv_mode: 'cv0(T)' (the fluid's ideal-gas c_v) or 'const'."""
    _, crit, st = real_fluid(name)
    R, Tc, Pc = crit['R'], crit['Tc'], crit['Pc']
    A = 27*R*R*Tc*Tc/(64*Pc) if kind == 'vdw' else 0.0
    B = R*Tc/(8*Pc) if kind == 'vdw' else 0.0

    def cv0(T):
        if cv_mode == 'const':
            T = Tc
        st.update(CP.DmassT_INPUTS, 1e-6, T)
        return st.cp0mass()-R

    def H(T, rho):
        v, cv = 1/rho, cv0(T)
        pT = R/(v-B)
        pv_T = -R*T/(v-B)**2+2*A/v**3
        return np.array([[T/cv, -T*pT/cv], [-T*pT/cv, -pv_T+T*pT*pT/cv]])
    return H


def holonomy(Hs, scale=(1.0, 1.0)):
    D = np.diag(scale)
    Hs = [D@h@D for h in Hs]
    n = len(Hs)
    a = np.array([h[0, 0] for h in Hs]); b = np.array([h[0, 1] for h in Hs]); c = np.array([h[1, 1] for h in Hs])
    u, v = (a-c)/(a+c), 2*b/(a+c)
    rho2 = u*u+v*v
    phi = np.unwrap(np.arctan2(v, u))
    dphi = np.diff(np.append(phi, phi[0]+round((phi[-1]-phi[0])/(2*math.pi))*2*math.pi))
    # periodic closure of the unwrapped angle
    dphi = np.array([((x+math.pi) % (2*math.pi))-math.pi for x in dphi])
    ch = 1/np.sqrt(1-rho2)
    ch_mid = (ch+np.roll(ch, -1))/2
    theta_loop = -0.5*float(np.sum((ch_mid-1)*dphi))
    w = np.eye(2)
    length = 0.0
    rap = 0.0
    for k in range(n):
        h0, h1 = Hs[k], Hs[(k+1) % n]
        hm = (h0+h1)/2
        X = np.linalg.solve(hm, h1-h0)
        g = np.trace(X@X)-0.5*np.trace(X)**2
        length += math.sqrt(max(g, 0.0))
        # second-order step of dw = -(1/2) H^-1 dH w
        step = np.eye(2)-0.5*X+0.125*X@X
        w = step@w
        ev = np.linalg.eigvalsh(h0)
        rap = max(rap, math.log(ev[1]/ev[0]))
    ev, O = np.linalg.eigh(Hs[0])
    root = O@np.diag(np.sqrt(ev))@O.T
    rot = root@w@np.linalg.inv(root)
    theta_tr = math.atan2(rot[1, 0]-rot[0, 1], rot[0, 0]+rot[1, 1])
    ortho = float(np.linalg.norm(rot@rot.T-np.eye(2)))
    return dict(theta_loop=theta_loop, theta_transport=theta_tr, orthogonality_defect=ortho,
                shape_length=length, max_rapidity=rap,
                bound=0.5*math.tanh(rap/4)*length, half_length=0.5*length)


def stage(name, T0, T1, r0, r1, n):
    T, rho = loop(n, T0, T1, r0, r1)
    Hreal, crit, _ = real_fluid(name)
    cases = {'reference_eos': Hreal,
             'van_der_waals_cv0(T)': model(name, 'vdw', 'cv0(T)'),
             'van_der_waals_const_cv': model(name, 'vdw', 'const'),
             'ideal_gas_cv0(T)': model(name, 'ideal', 'cv0(T)'),
             'ideal_gas_const_cv': model(name, 'ideal', 'const')}
    out = dict(fluid=name, cycle=dict(T_K=[T0, T1], rho_kg_m3=[r0, r1], shape='ellipse in (T, rho)', points=n),
               critical=crit, cases={})
    for key, H in cases.items():
        Hs = [H(float(t), float(r)) for t, r in zip(T, rho)]
        res = holonomy(Hs)
        res['unit_change_theta'] = holonomy(Hs, scale=(37.0, 0.004))['theta_transport']
        res['unit_change_length'] = holonomy(Hs, scale=(37.0, 0.004))['shape_length']
        out['cases'][key] = res
    return out


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    runs = [stage('CarbonDioxide', 310.0, 350.0, 200.0, 700.0, n),
            stage('Argon', 155.0, 200.0, 250.0, 800.0, n),
            stage('Nitrogen', 130.0, 170.0, 150.0, 480.0, n)]
    json.dump(runs, open('PH1_RESULT.json', 'w'), indent=1)
    for r in runs:
        print(r['fluid'], r['cycle']['T_K'], r['cycle']['rho_kg_m3'])
        for k, c in r['cases'].items():
            print('  %-24s Theta_loop=%+.6f Theta_tr=%+.6f unit=%+.6f L=%.4f (unit %.4f) r*=%.3f bound=%.4f' % (
                k, c['theta_loop'], c['theta_transport'], c['unit_change_theta'], c['shape_length'],
                c['unit_change_length'], c['max_rapidity'], c['bound']))
