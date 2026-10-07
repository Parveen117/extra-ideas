"""GR1: gravity as memory - the static dictionary and what G/c^2 means in it.

Static flow of MS1-T2 with a light-like (pure reading) inflow and accumulated flip G:
    n = cosh 2G,  j = 1,  sigma = sinh 2G,  N := v = j/n,  m := memory = 1 - N^2,  p = (1 + N)/2.
T1  local identities on any flip profile g(x):  d_x n = 2 g sigma,  d_x sigma = 2 g n,  d_x j = 0,
    hence  d_x N = -2 g N sqrt(m):  the change of N is flip rate x N x root of memory.
T2  dictionary (named physical input: m = r_s / r, r_s = 2 G_N M / c^2, the Newtonian 1/r law):
    N = sqrt(1 - r_s/r),  blueshift n = 1/N,  share p = (1 + N)/2,  horizon <=> p = 1/2,
    acceleration/c^2 = d_r N = r_s / (2 r^2 N)  with flip profile g(r) = beta / (4 r (1 - beta^2)), beta^2 = r_s/r.
T3  what is not derived: the 1/r law and the constant 2 G_N / c^2.
Exact rational arithmetic for T1, T2; plain floats only for the illustrative table.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../ms1')
import ms1_speed_and_mass as ms1


def local_control():
    """d_x psi = g K psi at a point: derivatives of the bilinears, exactly."""
    rows = []
    for psi, g in (((F(5, 4), F(3, 4)), F(2)), ((F(13, 12), F(5, 12)), F(1, 3)), ((F(3), F(1)), F(-1, 2)), ((F(1), F(0)), F(7))):
        n, j, sigma = ms1.state(psi)
        dpsi = (g*psi[1], g*psi[0])                               # g K psi
        dn = 2*(psi[0]*dpsi[0]+psi[1]*dpsi[1])
        dj = 2*(psi[0]*dpsi[0]-psi[1]*dpsi[1])
        dsigma = 2*(dpsi[0]*psi[1]+psi[0]*dpsi[1])
        if (dn, dj, dsigma) != (2*g*sigma, F(0), 2*g*n):
            raise ValueError('static derivative laws failed')
        N = j/n
        dN = -j*dn/(n*n)
        m = 1-N*N
        # dN = -2 g N sqrt(m), with sqrt(m) = sigma/n (sign of sigma kept)
        if dN != -2*g*N*(sigma/n) or (sigma/n)**2 != m:
            raise ValueError('dN = -2 g N sqrt(memory) failed')
        rows.append(dict(N=str(N), memory=str(m), flip_rate=str(g), dN=str(dN)))
    if rows[-1]['dN'] != '0' or rows[-1]['memory'] != '0':
        raise ValueError('with no memory nothing changes, whatever the flip rate')
    return rows


def dictionary_control():
    rows = []
    for beta in (F(3, 5), F(4, 5), F(5, 13), F(12, 13), F(8, 17), F(0)):
        m = beta*beta                                             # r_s / r
        N2 = 1-m
        # N rational for Pythagorean beta
        N = {F(3, 5): F(4, 5), F(4, 5): F(3, 5), F(5, 13): F(12, 13), F(12, 13): F(5, 13), F(8, 17): F(15, 17), F(0): F(1)}[beta]
        if N*N != N2:
            raise ValueError('N^2 + memory = 1 failed')
        n, j, sigma = 1/N, F(1), beta/N                           # cosh eta, 1, sinh eta
        if n*n != j*j+sigma*sigma:
            raise ValueError('state tensor is not null')
        p = (n+j)/(2*n)
        if p != (1+N)/2 or 4*p*(1-p) != m:
            raise ValueError('share and memory failed')
        for r in (F(1), F(7, 2), F(10)):
            r_s = m*r
            if beta:
                g = beta/(4*r*(1-m))                              # flip profile needed for m = r_s/r
                acc = 2*g*N*beta                                  # T1: |d_r N| = 2 g N sqrt(m)
                if acc != r_s/(2*r*r*N):
                    raise ValueError('acceleration is not r_s / (2 r^2 N)')
                # the profile is d(eta/2)/dr: with tanh(eta) = beta(r) = sqrt(r_s/r), d beta/dr = -beta/(2r)
                deta = (1/(1-m))*(beta/(2*r))
                if g != deta/2:
                    raise ValueError('flip profile is not half the rapidity gradient')
        rows.append(dict(r_s_over_r=str(m), N=str(N), blueshift=str(n), share=str(p), proper_density=str(sigma)))
    return rows


def horizon_control():
    """As memory -> 1 the share -> 1/2 from above and N -> 0; never reached at finite accumulated flip."""
    vals = []
    for b in (F(2), F(5), F(50), F(1000)):                        # e^G
        ch2 = ((b*b)+(1/(b*b)))/2                                 # cosh 2G
        N = 1/ch2
        p = (1+N)/2
        vals.append((p, 1-N*N))
    if any(p <= F(1, 2) or m >= 1 for p, m in vals):
        raise ValueError('horizon must not be reached at finite flip')
    if not all(a[0] > b[0] and a[1] < b[1] for a, b in zip(vals, vals[1:])):
        raise ValueError('share must fall to 1/2 and memory rise to 1')
    return dict(last_share=str(vals[-1][0]), last_memory_gap=str(1-vals[-1][1]))


def illustration():
    G_N, c = 6.67430e-11, 299792458.0
    out = dict(two_G_over_c2_m_per_kg=2*G_N/c**2)
    for name, M, Rr in (('Earth surface', 5.972e24, 6.371e6), ('Sun surface', 1.989e30, 6.957e8),
                        ('neutron star 1.4 Msun, 12 km', 1.4*1.989e30, 1.2e4)):
        m = 2*G_N*M/(Rr*c*c)
        N = (1-m)**0.5
        out[name] = dict(memory=m, N=N, share=(1+N)/2, acceleration_m_s2=c*c*m/(2*Rr*N))
    return out


def run():
    return dict(local=local_control(), dictionary=dictionary_control(), horizon=horizon_control(),
                illustration=illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('GR1_RESULT.json', 'w'), indent=1)
    for k in ('local', 'dictionary'):
        for r in res[k]:
            print(r)
    print(res['horizon'])
    for k, v in res['illustration'].items():
        print(k, v)
