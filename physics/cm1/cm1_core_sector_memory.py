"""CM1: turning-sector core readings and their exact missing-source Gram.

The physical carrier is the free 3xd matrix C of TC1/CR1, h=-Delta/2+e2.
Pairing is inherited from the declared Gaussian C-space realization, not a
primitive new framework axiom. Proof arithmetic uses Fraction only.
"""
from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cr = load(PHYSICS/'cr1/cr1_core_rates.py', 'cm1_cr1')
tc = cr.tc1
ZERO = (0, 0, 0)
E1, E2, E3 = ({(1, 0, 0): F(1)}, {(0, 1, 0): F(1)}, {(0, 0, 1): F(1)})


def inv_basis(top):
    return sorted([(a, b, c) for a in range(top+1) for b in range(top//2+1)
                   for c in range(top//3+1) if a+2*b+3*c <= top],
                  key=lambda e: (cr.degree(e), e))


def scalar_h(p, d, w):
    out = cr.padd({}, cr.gen(p, d), -1)
    out = cr.padd(out, {e: (2*w*cr.degree(e)+F(3*d)*w/2)*v for e, v in p.items()})
    out = cr.padd(out, cr.pmul(E1, p), -w*w/2)
    return cr.padd(out, cr.pmul(E2, p))


@lru_cache(None)
def power_sum(r):
    if r == 1:
        return E1
    if r == 2:
        return cr.padd(cr.pmul(E1, E1), E2, -2)
    if r == 3:
        return cr.padd(cr.padd(cr.pmul(E1, cr.pmul(E1, E1)), cr.pmul(E1, E2), -3), E3, 3)
    return cr.padd(cr.padd(cr.pmul(E1, power_sum(r-1)), cr.pmul(E2, power_sum(r-2)), -1),
                   cr.pmul(E3, power_sum(r-3)))


def channel_product(r, s, d):
    """Orientation-average of t_r t_s, t_r=(M^r)11-(M^r)22."""
    p = cr.padd(power_sum(r+s), cr.pmul(power_sum(r), power_sum(s)), -F(1, d))
    factor = F(4, (d-1)*(d+2))
    return {e: factor*v for e, v in p.items()}


def add_term(out, r, e, v, d):
    if not v or not r:
        return
    # For d=3, M^3=e1 M^2-e2 M+e3 I, and t0=0.
    # For d=4 use the rank-three M^4=e1 M^3-e2 M^2+e3 M.
    if (d == 3 and r >= 3) or r >= 4:
        for k, sign, shift in [(1, 1, (1,0,0)), (2, -1, (0,1,0)), (3, 1, (0,0,1))]:
            add_term(out, r-k, tuple(a+b for a,b in zip(e,shift)), sign*v, d)
    else:
        key = (r,)+e
        out[key] = out.get(key, F(0))+v
        if not out[key]:
            del out[key]


def spin_l(poly, d):
    """L=Delta_C/2 on sum_r t_r p_r(e), retaining all cross derivatives."""
    out = {}
    for (r,a,b,c), value in poly.items():
        e = (a,b,c)
        for f,v in cr.gen_mono(e,d).items():
            add_term(out,r,f,value*v,d)
        if r == 2:
            add_term(out,1,e,value*(d+8),d)
        elif r == 3:
            add_term(out,2,e,value*(2*d+15),d)
            add_term(out,1,(a+1,b,c),value*3,d)
        # Gamma(t_r,e1)=4r t_r.
        if a:
            add_term(out,r,(a-1,b,c),value*4*r*a,d)
        # Gamma(t_r,e2)=4r(e1 t_r-t_(r+1)).
        if b:
            add_term(out,r,(a+1,b-1,c),value*4*r*b,d)
            add_term(out,r+1,(a,b-1,c),-value*4*r*b,d)
        # Gamma(t_r,e3)=4r(e2 t_r-e1 t_(r+1)+t_(r+2)).
        if c:
            add_term(out,r,(a,b+1,c-1),value*4*r*c,d)
            add_term(out,r+1,(a+1,b,c-1),-value*4*r*c,d)
            add_term(out,r+2,(a,b,c-1),value*4*r*c,d)
    return out


def spin_h(poly, d, w):
    out = {}
    for key,v in spin_l(poly,d).items():
        add_term(out,key[0],key[1:],-v,d)
    for (r,a,b,c),v in poly.items():
        e=(a,b,c)
        add_term(out,r,e,v*(2*w*(r+cr.degree(e))+F(3*d)*w/2),d)
        add_term(out,r,(a+1,b,c),-v*w*w/2,d)
        add_term(out,r,(a,b+1,c),v,d)
    return out


def spin_inner(p, q, d, mean):
    out = {}
    for (r,*e),x in p.items():
        for (s,*f),y in q.items():
            base=tuple(a+b for a,b in zip(e,f))
            for z,v in channel_product(r,s,d).items():
                k=tuple(a+b for a,b in zip(base,z))
                out[k]=out.get(k,F(0))+x*y*v
    return mean(out)


def spin_basis(d, top):
    return sorted([(r,)+e for r in range(1,3 if d==3 else 4)
                   for e in inv_basis(top-r)],key=lambda k:(k[0]+cr.degree(k[1:]),k))


def spin_ladder(d, w, top):
    keys=spin_basis(d,top)
    polynomials=[{k:F(1)} for k in keys]
    mean=cr.free_means(d,w)
    hp=[spin_h(p,d,w) for p in polynomials]
    s=[[spin_inner(p,q,d,mean) for q in polynomials] for p in polynomials]
    h=[[spin_inner(p,q,d,mean) for q in hp] for p in polynomials]
    if s != transpose(s) or h != transpose(h):
        raise ValueError('pairing or generator is not symmetric')
    return keys,s,h


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def subtract(a,b):
    return [[x-y for x,y in zip(row,other)] for row,other in zip(a,b)]


def solve_spd(a, b):
    """Exact LDL solve for multiple RHS, refusing nonpositive pivots."""
    n=len(a)
    if a != transpose(a):
        raise ValueError('non-Hermitian Gram')
    l=[[F(int(i==j)) for j in range(n)] for i in range(n)]
    diagonal=[]
    for j in range(n):
        pivot=a[j][j]-sum(l[j][k]*l[j][k]*diagonal[k] for k in range(j))
        if pivot <= 0:
            raise ValueError('Gram is not strictly positive')
        diagonal.append(pivot)
        for i in range(j+1,n):
            l[i][j]=(a[i][j]-sum(l[i][k]*l[j][k]*diagonal[k] for k in range(j)))/pivot
    m=len(b[0])
    y=[[F(0)]*m for _ in range(n)]
    for i in range(n):
        y[i]=[b[i][j]-sum(l[i][k]*y[k][j] for k in range(i)) for j in range(m)]
    z=[[v/diagonal[i] for v in row] for i,row in enumerate(y)]
    x=[[F(0)]*m for _ in range(n)]
    for i in reversed(range(n)):
        x[i]=[z[i][j]-sum(l[k][i]*x[k][j] for k in range(i+1,n)) for j in range(m)]
    return x


def exact_rank(a):
    rows=[list(row) for row in a]
    rank=0
    for col in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        p=rows[rank][col]
        for i in range(rank+1,len(rows)):
            if rows[i][col]:
                factor=rows[i][col]/p
                rows[i]=[x-factor*y for x,y in zip(rows[i],rows[rank])]
        rank+=1
        if rank==len(rows):
            break
    return rank


def source_packet(d,w,top):
    """Actual Q h P columns, not an arbitrary finite matrix fixture."""
    basis,s,h=cr.ladder(d,w,top)
    ext=inv_basis(top+2)
    index={e:i for i,e in enumerate(ext)}
    n=len(basis)
    hp=[scalar_h({e:F(1)},d,w) for e in basis]
    mean=cr.free_means(d,w)
    second=[[mean(cr.pmul(p,q)) for q in hp] for p in hp]
    action=solve_spd(s,h)
    gram=subtract(second,matmul(h,action))
    k=[[F(0)]*n for _ in ext]
    for j,p in enumerate(hp):
        for e,v in p.items():
            k[index[e]][j]=v
        for i,e in enumerate(basis):
            k[index[e]][j]-=action[i][j]
    ext_s=[[mean({tuple(a+b for a,b in zip(e,f)):F(1)}) for f in ext] for e in ext]
    direct=matmul(transpose(k),matmul(ext_s,k))
    if direct != gram:
        raise ValueError('missing-source Gram does not reconstruct')
    if any(sum(ext_s[i][k0]*k[k0][j] for k0 in range(len(ext)))
           for i in range(n) for j in range(n)):
        raise ValueError('source columns are not orthogonal to retained cut')
    predicted=len(basis)-len(inv_basis(top-2)) if top>=2 else len(basis)
    tail_coeff=k[n:]
    rank=exact_rank(tail_coeff)
    if rank != predicted:
        raise ValueError('source-rank formula failed')
    # All columns from degrees <= D-2 are exactly zero, not just small.
    for j,e in enumerate(basis):
        if cr.degree(e)<=top-2 and any(row[j] for row in k):
            raise ValueError('low-degree source should vanish')
    return dict(basis=basis,s=s,h=h,gram=gram,columns=k,expanded_gram=ext_s,
                second=second,rank=rank,shell_dimension=len(ext)-n)


def leading_source_rank(d,w,top):
    """Actual above-cut action coefficients; no Gram solve needed for this rank."""
    basis=inv_basis(top)
    tail=[e for e in inv_basis(top+2) if cr.degree(e)>top]
    actions=[scalar_h({e:F(1)},d,w) for e in basis]
    return exact_rank([[p.get(e,F(0)) for p in actions] for e in tail])


def floor_count(packet,z,d0):
    if d0 <= z:
        raise ValueError('hidden floor must exceed energy threshold')
    a=[[h-z*s-g/(d0-z) for h,s,g in zip(hr,sr,gr)]
       for hr,sr,gr in zip(packet['h'],packet['s'],packet['gram'])]
    return cr.below([[F(int(i==j)) for j in range(len(a))] for i in range(len(a))],a,F(0))


def required_floor(packet,z,upper=F(100),width=F(1,1000)):
    """Conditional hidden-floor demand, never a proved physical floor."""
    lo=z
    hi=upper
    if floor_count(packet,z,hi)>1:
        raise ValueError('chosen upper floor does not suffice')
    while hi-lo>width:
        mid=(hi+lo)/2
        if floor_count(packet,z,mid)<=1:
            hi=mid
        else:
            lo=mid
    return lo,hi,floor_count(packet,z,hi)


def entry_m_powers(d):
    n=3*d
    def variable(k):
        e=[0]*n;e[k]=1
        return {tuple(e):F(1)}
    c=[[variable(a*d+i) for i in range(d)] for a in range(3)]
    m=[[{} for _ in range(d)] for _ in range(d)]
    for i in range(d):
        for j in range(d):
            for a in range(3):
                m[i][j]=tc.padd(m[i][j],tc.pmul(c[a][i],c[a][j]))
    powers=[m]
    for _ in range(2):
        last=powers[-1]
        out=[[{} for _ in range(d)] for _ in range(d)]
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    out[i][j]=tc.padd(out[i][j],tc.pmul(last[i][k],m[k][j]))
        powers.append(out)
    return [tc.padd(p[0][0],p[1][1],-1) for p in powers]


def spin_to_entries(poly,d):
    tr=entry_m_powers(d)
    gens=tc.entry_polys(d)
    out={}
    for (r,*e),value in poly.items():
        inv=cr.to_entries({tuple(e):value},gens)
        out=tc.padd(out,tc.pmul(tr[r-1],inv))
    return out


def source_pins():
    paths=[PHYSICS/'tc1/TC1_CORE_OF_TURNS.md',PHYSICS/'tc1/tc1_core_of_turns.py',
           PHYSICS/'cr1/CR1_CORE_RATES.md',PHYSICS/'cr1/cr1_core_rates.py',
           PHYSICS/'RH_YM_TRANSFER_MAP.md',HERE/'cm1_core_sector_memory.py',
           HERE/'CM1_CORE_SECTOR_MEMORY.md',HERE/'test_cm1.py']
    return {str(p.relative_to(PHYSICS.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run():
    checks={};numbers={}
    for d in (3,4):
        for r in (1,2,3) if d==4 else (1,2):
            # Independent entry-coordinate Laplacian route, including mixed derivatives.
            fixtures=[{(r,0,0,0):F(1)},{(r,1,0,0):F(1)},{(r,0,1,0):F(1)}]
            if r==1:fixtures.append({(r,0,0,1):F(1)})
            checks[f'generator entry identities d={d} r={r}']=all(
                cr.entry_laplacian_half(spin_to_entries(p,d))==spin_to_entries(spin_l(p,d),d)
                for p in fixtures)
        mean=cr.free_means(d,F(1))
        tr=entry_m_powers(d)
        checks[f'orientation average versus entry Wick pairing d={d}']=all(
            mean(channel_product(r,s,d))==tc.free_mean(tc.pmul(tr[r-1],tr[s-1]))
            for r,s in [(1,1),(1,2),(2,2),(1,3)])
        rows=[]
        w=cr.RATE[d]
        for top in (2,4,6,8):
            keys,s,h=spin_ladder(d,w,top)
            solve_spd(s,[[F(int(i==j)) for j in range(len(s))] for i in range(len(s))])
            lo,hi=cr.bracket(s,h,0,F(0),F(20),len(keys),grid=10000)
            rows.append(dict(degree=top,size=len(keys),ritz_lower=str(lo),ritz_upper=str(hi)))
            checks[f'spin-sector exact upper count d={d} D={top}']=cr.below(s,h,hi)>=1
            checks[f'spin-sector symmetry and positive Gram d={d} D={top}']=s==transpose(s) and h==transpose(h)
        numbers[f'turning_sector_d{d}']=rows
        checks[f'turning Ritz bound below CR1 second upper d={d}']=F(rows[-1]['ritz_upper']) < (F(8048,1000) if d==3 else F(11566,1000))
        checks[f'nested turning Ritz bounds d={d}']=all(F(b['ritz_upper'])<=F(a['ritz_upper']) for a,b in zip(rows,rows[1:]))
        packets=[]
        for top in (2,4,6):
            packet=source_packet(d,w,top)
            checks[f'exact source reconstruction and orthogonality d={d} D={top}']=True
            checks[f'actual source rank d={d} D={top}']=packet['rank']==len(inv_basis(top))-len(inv_basis(top-2))
            # A norm-free check that a coordinate dot product is not the Gram metric.
            naive=matmul(transpose(packet['columns']),packet['columns'])
            checks[f'Euclidean-source-Gram negative control d={d} D={top}']=naive!=packet['gram']
            z=F(6) if d==3 else F(9)
            low,high,count=required_floor(packet,z)
            packets.append(dict(degree=top,retained=len(packet['s']),shell_dimension=packet['shell_dimension'],
                                source_rank=packet['rank'],test_energy=str(z),
                                conditional_hidden_floor_lower=str(low),conditional_hidden_floor_upper=str(high),
                                lower_comparison_negative_count=count))
            checks[f'conditional floor is above CR1 global lower d={d} D={top}']=high>(F(414,100) if d==3 else F(632,100))
        numbers[f'scalar_source_d{d}']=packets
    ranks={str(d):leading_source_rank(d,cr.RATE[d],10) for d in (3,4)}
    checks['actual source rank 26 at D10, shell dimension 35']=all(
        rank==26 for rank in ranks.values()) and len(inv_basis(12))-len(inv_basis(10))==35
    numbers['scalar_source_D10']={'retained':len(inv_basis(10)),
        'next_shell_dimension':len(inv_basis(12))-len(inv_basis(10)),
        'actual_above_cut_source_rank_by_d':ranks,'full_source_gram':'NOT_COMPUTED_AT_D10'}
    return {'schema':'cm1_core_sector_memory_v1','arithmetic':'Fraction; exact polynomial means and inertia',
            'checks':checks,'numbers':numbers,
            'claim_boundary':{'turning_sector_ritz_upper':'PROVED_IN_DECLARED_SECTOR',
               'scalar_source_gram_and_rank':'EXACT',
               'conditional_schur_floor_demand':'EXACT_CONDITIONAL_NOT_AN_ACTUAL_HIDDEN_FLOOR',
               'excited_energy_lower':'OPEN','quantitative_core_gap_lower':'OPEN',
               'qualitative_core_gap':'WRITTEN_PROOF_WITH_STANDARD_ANALYTIC_INPUTS_NOT_A_FINITE_CERTIFICATE',
               'continuum_mass_gap':'OPEN'}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if not all(result['checks'].values()):
        raise SystemExit('FAIL: '+', '.join(k for k,v in result['checks'].items() if not v))
    result['source_sha256']=source_pins()
    path=HERE/'CM1_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result,indent=2)+'\n')
    if args.check and json.loads(path.read_text())!=result:
        raise SystemExit('FAIL: stale CM1 result/source pin')
    print('PASS',len(result['checks']),'exact checks')
    print(json.dumps(result['numbers'],indent=2))


if __name__=='__main__':
    main()
