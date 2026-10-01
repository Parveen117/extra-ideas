#!/usr/bin/env node
'use strict';
// Native application packet. Exact arithmetic and word proofs are supplied
// by the unchanged canonical RKF engine; no second engine is introduced.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const arg=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(arg('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs')),p=require(path.join(home,'core/paninian_operator.cjs')),w=require(path.join(home,'core/workbench.cjs'));
const ensure=(v,msg)=>{if(!v)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
const input=JSON.parse(fs.readFileSync(arg('--input'),'utf8')),inputHash=p.digest(input);
ensure(inputHash===arg('--expected-input-sha256'),'R31 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r30-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r30-sha256'),'R31 R30 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R30_NATIVE_DIRECTIONAL_LOOP_INTERACTION','R31 parent status');
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const I=o.identity(2),H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange),R=o.mul(K,H);
const P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),C0=o.add(H,K),A0=o.sub(H,K),E01=o.mul(P,K),E10=o.mul(Q,K);
const col=a=>o.matrix(a.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const key=(x,y)=>x+','+y,xy=k=>k.split(',').map(Number),unit=X=>new Map([['0,0',X]]);
const addAt=(a,k,M)=>{const v=o.add(a.get(k)||o.zeros(M.length,M[0].length),M);if(o.isZero(v))a.delete(k);else a.set(k,v);};
const plus=(a,b)=>{const c=new Map(a);for(const [k,M]of b)addAt(c,k,M);return c;};
const scale=(a,c)=>new Map([...a].map(([k,M])=>[k,o.scale(M,c)])),minus=(a,b)=>plus(a,scale(b,-1));
const shift=(a,dx,dy)=>new Map([...a].map(([k,M])=>{const [x,y]=xy(k);return [key(x+dx,y+dy),M];}));
const compose=(a,b)=>{const c=new Map();for(const [i,A]of a)for(const [j,B]of b){const [x,y]=xy(i),[u,v]=xy(j);addAt(c,key(x+u,y+v),o.mul(A,B));}return c;};
const dag=a=>new Map([...a].map(([k,M])=>{const [x,y]=xy(k);return [key(-x,-y),o.dagger(M)];}));
const eqFields=(a,b,r,c,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eq(a.get(k)||o.zeros(r,c),b.get(k)||o.zeros(r,c),msg+' at '+k);};
const power=(a,n,size=2)=>{let b=unit(o.identity(size));for(let i=0;i<n;i++)b=compose(a,b);return b;};
const energy=a=>[...a.values()].reduce((s,M)=>s.add(o.energy(M)),f(0));
const value=(a,x,y)=>a.get(key(x,y))||o.zeros(2,1);
const dot=(a,b)=>[...a].reduce((s,[k,v])=>s.add(o.mul(o.dagger(v),b.get(k)||o.zeros(v.length,1))[0][0].rad),f(0));
const lift=(A,a)=>new Map([...a].map(([k,M])=>[k,o.kron(A,M)]));
const at=(op,phi,x,y)=>[...op].reduce((s,[k,M])=>{const [dx,dy]=xy(k);return o.add(s,o.mul(M,phi(x-dx,y-dy)));},o.zeros(2,1));
const tx=shift(unit(I),1,0),ty=shift(unit(I),0,1),dx=minus(tx,unit(I)),dy=minus(ty,unit(I));
const nx=plus(unit(o.scale(o.mul(P,C0),'1/2')),shift(unit(o.scale(o.mul(Q,C0),'1/2')),-1,0));
const ny=plus(unit(o.scale(o.mul(P,C0),'1/2')),shift(unit(o.scale(o.mul(Q,C0),'1/2')),0,-1));
const vx=plus(unit(I),compose(nx,dx)),vy=plus(unit(I),compose(ny,dy)),vxD=dag(vx),vyD=dag(vy);
const forward=compose(vy,vx),reverse=compose(vx,vy),dual=compose(vyD,vxD);
const omega=minus(reverse,forward),g=compose(compose(dx,dy),minus(dag(ty),dag(tx)));
const loop=compose(compose(compose(vx,vy),vxD),vyD),loopDefect=minus(loop,unit(I));
const comm=(a,b)=>minus(compose(a,b),compose(b,a));
const lapx=compose(dag(dx),dx),lapy=compose(dag(dy),dy),sumlap=plus(lapx,lapy);
const skewx=minus(tx,dag(tx)),skewy=minus(ty,dag(ty));
const wave=minus(scale(sumlap,'1/2'),scale(compose(skewx,skewy),'1/4'));
const waveCoefficient=minus(unit(o.scale(I,2)),wave),increment=minus(forward,unit(I));
const scalarD=scale(sumlap,'1/4');
const skewQ=scale(plus(compose(minus(lapx,lapy),unit(R)),compose(plus(skewx,skewy),unit(A0))),'1/4');
const sx=compose(dx,vy),shear=plus(sx,compose(dy,vx));
const D=plus(unit(o.identity(4)),lift(E01,shear)),Di=minus(unit(o.identity(4)),lift(E01,shear));
const cx=plus(lift(I,vx),lift(E01,compose(dx,omega))),cy=minus(lift(I,vy),lift(E01,compose(dy,omega)));
const joint=compose(cy,cx),fullGram=compose(dag(D),D);
const stack=(a,b)=>{const out=new Map();for(const k of new Set([...a.keys(),...b.keys()]))out.set(k,[...(a.get(k)||o.zeros(2,1)),...(b.get(k)||o.zeros(2,1))]);return out;};
const split=a=>{const b=new Map(),phi=new Map();for(const [k,v]of a){addAt(b,k,v.slice(0,2));addAt(phi,k,v.slice(2));}return {b,phi};};
const probes=[unit(e0),unit(e1),new Map([['-1,0',col([2,-1])],['1,1',col([-1,3])],['0,-2',col([1,1])]])];
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

check('one_local_shear_unmixes_both_directional_field_updates',()=>{
 eqFields(comm(vx,shear),compose(dx,omega),2,2,'x drive is the shear commutator');
 eqFields(comm(vy,shear),scale(compose(dy,omega),-1),2,2,'y drive has reversed loop orientation');
 eqFields(compose(D,Di),unit(o.identity(4)),4,4,'local shear inverse');
 eqFields(compose(Di,D),unit(o.identity(4)),4,4,'local shear inverse other side');
 for(const [c,v]of [[cx,vx],[cy,vy]])eqFields(compose(D,c),compose(lift(I,v),D),4,4,'simultaneous source unmixing');
 eqFields(compose(D,joint),compose(lift(I,forward),D),4,4,'alternating block unmixing');
 eqFields(compose(compose(D,comm(cx,cy)),Di),lift(I,omega),4,4,'joint directional order residue is two source copies');
 const Dx=plus(unit(o.identity(4)),lift(E01,sx));
 eqFields(compose(Dx,cx),compose(lift(I,vx),Dx),4,4,'R30 x-only target has its smaller shear');
 return {complete_Laurent_intertwiner_and_inverse_identities:9,simultaneously_unmixed_directions:2,
  primitive_physical_interaction_inferred_from_off_diagonal_drive:false};
});

check('complete_source_energy_and_driven_solution_are_bounded',()=>{
 for(const c of [cx,cy,joint])eqFields(compose(compose(dag(c),fullGram),c),fullGram,4,4,'full source form is conserved');
 let identities=0,energies=0,bounds=0;
 for(let j=0;j<probes.length;j++){
  const phi0=probes[j],b0=j?probes[j-1]:new Map(),B0=plus(b0,compose(shear,phi0));
  let z=stack(b0,phi0),B=B0,phi=phi0,bfree=b0,preparationS=compose(shear,phi0);
  const norm=energy(B0).add(energy(phi0));
  for(let n=0;n<=7;n++){
   const parts=split(z),formula=minus(plus(bfree,preparationS),compose(shear,phi));
   eqFields(parts.b,formula,2,1,'exact driven solution');eqFields(parts.phi,phi,2,1,'retained source solution');
   eqFields(compose(D,z),stack(B,phi),4,1,'independent continued channels');identities+=3;
   ensure(energy(compose(D,z)).eq(norm),'complete positive energy');energies++;
   ensure(energy(minus(parts.b,bfree)).le(energy(phi0).mul(64)),'uniform block drive bound independent of duration');bounds++;
   if(n<7){z=compose(joint,z);B=compose(forward,B);phi=compose(forward,phi);bfree=compose(forward,bfree);preparationS=compose(forward,preparationS);}
  }
 }
 return {complete_Laurent_Gram_identities:3,exact_channel_and_solution_matches:identities,conserved_complete_energy_checks:energies,
  squared_drive_bound_checks:bounds,alternating_drive_norm_bound_coefficient:8,max_four_event_blocks:7,
  complete_source_form_is_distinct_from_R30_observable_form:true};
});

check('alternating_source_and_coupled_fields_obey_the_derived_wave',()=>{
 eqFields(plus(forward,dag(forward)),waveCoefficient,2,2,'exact two-direction scalar coefficient');
 eqFields(plus(minus(compose(forward,forward),compose(waveCoefficient,forward)),unit(I)),new Map(),2,2,'native alternating wave polynomial');
 eqFields(plus(minus(compose(joint,joint),compose(lift(I,waveCoefficient),joint)),unit(o.identity(4))),new Map(),4,4,'joint wave polynomial');
 eqFields(plus(dual,dag(dual)),waveCoefficient,2,2,'dual curvature wave has the same coefficient');
 let waves=0;
 for(const seed of probes){let old=seed,now=compose(forward,seed);for(let n=1;n<=6;n++){
  const next=compose(forward,now);eqFields(plus(minus(next,scale(now,2)),old),scale(compose(wave,now),-1),2,1,'finite native wave');waves++;
  old=now;now=next;
 }}
 return {complete_Laurent_wave_identities:4,finite_source_wave_checks:waves,mixed_difference_coefficient:'-1/4',
  classical_wave_equation_or_metric_used:false};
});

check('positive_native_squares_classify_stationary_source',()=>{
 eqFields(dag(skewQ),scale(skewQ,-1),2,2,'derived Q is skew under native pairing');
 eqFields(compose(unit(o.dagger(R)),minus(vx,vyD)),plus(scalarD,skewQ),2,2,'native rotated source difference');
 eqFields(comm(scalarD,skewQ),new Map(),2,2,'scalar D commutes with Q');
 eqFields(wave,plus(compose(scalarD,scalarD),compose(dag(skewQ),skewQ)),2,2,'positive native sum of squares');
 eqFields(compose(dag(increment),increment),wave,2,2,'increment Gram is the scalar wave response');
 const anti=scale(minus(forward,dag(forward)),'1/2');
 eqFields(minus(wave,scale(compose(wave,wave),'1/4')),compose(dag(anti),anti),2,2,'native upper bound without spectral theorem');
 let norms=0;for(const phi of probes){const a=energy(compose(increment,phi));
  ensure(a.eq(energy(compose(scalarD,phi)).add(energy(compose(skewQ,phi)))),'positive energy decomposition');
  ensure(a.le(energy(phi).mul(4)),'increment squared bound');norms+=2;
 }
 let affineChecks=0;for(const v of [e0,e1])for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++){
  eq(at(increment,(a,b)=>o.scale(v,a-b),x,y),o.zeros(2,1),'unbounded affine stationary witness outside finite-pairing domain');affineChecks++;
 }
 return {complete_Laurent_factorization_identities:6,finite_positive_form_checks:norms,unbounded_affine_boundary_checks:affineChecks,
  stationary_source_requires_both_address_differences_zero_on_finite_pairing_domain:true,spectral_theorem_used:false};
});

const reduce=(a,u,v)=>{const out=new Map();for(const [k,M]of a){const [x,y]=xy(k);addAt(out,key(u*x+v*y,0),M);}return out;};
const choose=(n,k)=>{let b=f(1);for(let j=0;j<k;j++)b=b.mul(n-j).div(j+1);return b;};
// Finite binomial coefficients solve multiplication by (1+p)^n modulo
// total degree three, including inverse shifts; no analytic Taylor input.
const jet=(a,degree)=>{const out=new Map();for(const [k,M]of a){const [x,y]=xy(k);for(let i=0;i<=degree;i++)for(let j=0;j+i<=degree;j++)addAt(out,key(i,j),o.scale(M,choose(x,i).mul(choose(y,j))));}return out;};
check('directional_count_quotients_expose_a_quartic_transverse_sector',()=>{
 eqFields(reduce(wave,1,0),scale(lapx,'1/2'),2,2,'single-direction wave');
 eqFields(reduce(wave,1,1),minus(scale(lapx,2),scale(compose(lapx,lapx),'1/4')),2,2,'equal-shift wave');
 eqFields(reduce(wave,1,-1),scale(compose(lapx,lapx),'1/4'),2,2,'opposite-shift quartic wave');
 eqFields(reduce(omega,1,1),new Map(),2,2,'equal-shift loop blindness');
 eqFields(reduce(omega,1,-1),scale(compose(compose(lapx,skewx),unit(A0)),'1/4'),2,2,'opposite-shift curvature');
 const quadratic=new Map([['2,0',o.scale(I,'-1/2')],['1,1',o.scale(I,-1)],['0,2',o.scale(I,'-1/2')]]);
 eqFields(jet(wave,2),quadratic,2,2,'complete degree-two count jet');
 eqFields(jet(increment,1),new Map([['1,0',o.scale(C0,'1/2')],['0,1',o.scale(C0,'1/2')]]),2,2,'one leading source direction');
 const form=o.matrix([['1/2','1/2'],['1/2','1/2']]);ensure(o.rank(form)===1,'native first-difference coefficient form rank');
 eq(o.mul(form,col([1,-1])),o.zeros(2,1),'transverse coefficient null direction');
 return {complete_reduction_and_jet_identities:7,native_quadratic_coefficient_rank:1,
  exact_transverse_wave_operator:'Delta^2/4',physical_two_dimensional_metric_selected:false};
});

const ring=(op,L,r=2,c=2)=>{const M=o.zeros(r*L*L,c*L*L);for(const [k,A]of op){const [dx,dy]=xy(k);
 for(let x=0;x<L;x++)for(let y=0;y<L;y++){const u=((x-dx)%L+L)%L,v=((y-dy)%L+L)%L,a=x*L+y,b=u*L+v;
  for(let i=0;i<r;i++)for(let j=0;j<c;j++)M[a*r+i][b*c+j]=M[a*r+i][b*c+j].add(A[i][j]);
 }}return M;};
const average=(L,a,b,weighted=false)=>{const op=new Map();for(let k=0;k<L;k++)addAt(op,key(k*a,k*b),o.scale(I,new o.F(weighted?k:1,L)));return op;};
// Each positive pivot is an exact completed native square; a zero pivot
// must have a zero remainder row. This checks finite witnesses, not a theorem.
const psd=M=>{eq(M,o.dagger(M),'self-dagger finite form');const a=M.map(row=>row.map(v=>{ensure(v.turn.zero(),'real signed role form');return v.rad;}));let positives=0;
 for(let k=0;k<a.length;k++){const pivot=a[k][k];ensure(f(0).le(pivot),'negative native square pivot');if(pivot.eq(0)){ensure(a[k].slice(k).every(x=>x.eq(0)),'nonzero zero-pivot row');continue;}positives++;
  for(let i=k+1;i<a.length;i++)for(let j=i;j<a.length;j++){a[i][j]=a[i][j].sub(a[i][k].mul(a[k][j]).div(pivot));a[j][i]=a[i][j];}
 }return positives;};
const cyclic=[];
check('finite_count_averages_give_an_explicit_propagation_gap_bound',()=>{
 let identities=0;const rows=[];
 for(let L=1;L<=5;L++){
  const n=2*L*L,id=o.identity(n),ax=average(L,1,0),ay=average(L,0,1),a0=compose(ax,ay),A=ring(a0,L),U=ring(forward,L),F=ring(increment,L),Lap=ring(wave,L),remain=o.sub(id,A);
  eq(o.mul(A,A),A,'uniform native average cut');eq(o.mul(F,A),o.zeros(n),'uniform stationary sector');identities+=2;
  let mu=null,pivots=0;
  if(L>1){const bx=average(L,1,0,true),by=average(L,0,1,true);
   eq(ring(compose(dx,bx),L),o.sub(id,ring(ax,L)),'x weighted count inverse');eq(ring(compose(dy,by),L),o.sub(id,ring(ay,L)),'y weighted count inverse');
   eq(o.add(ring(compose(dx,bx),L),ring(compose(ax,compose(dy,by)),L)),remain,'orthogonal two-direction mean removal');identities+=3;
   const count=f(L-1).pow(2).div(4);psd(o.sub(o.scale(ring(sumlap,L),count),remain));
   mu=f(1).div(f(L-1).pow(4));pivots=psd(o.sub(Lap,o.scale(remain,mu)));
  }
  rows.push({period:L,source_roles:n,mean_free_lower_bound:mu?mu.toString():null,positive_gap_witness_pivots:pivots});cyclic.push({L,n,id,A,U,F,Lap,remain,mu});
 }
 return {finite_average_and_inverse_identities:identities,count_quotients:rows,gap_lower_bound_formula:'1/(L-1)^4 for L>=2',
  bound_claimed_sharp:false,physical_particle_mass_identified:false};
});

check('stationary_and_loop_blind_sectors_are_distinct',()=>{
 const rows=[];let positiveJointForms=0;
 for(const {L,n,A,F,Lap}of cyclic){const Om=ring(omega,L),J=ring(joint,L,4,4);
  ensure(o.rank(F)===n-2&&o.rank(Lap)===n-2,'exact stationary source rank');
  ensure(o.rank(o.sub(J,o.identity(2*n)))===2*n-4,'joint stationary sector has four native roles');
  const blind=2*(3*L-2);ensure(o.rank(Om)===n-blind,'previous curvature kernel retained');
  if(L<=4){ensure(psd(ring(fullGram,L,4,4))===2*n,'complete joint form is positive before observable quotient');positiveJointForms++;}
  rows.push({period:L,stationary_source_roles:2,loop_blind_source_roles:blind,dynamical_loop_blind_roles:blind-2,joint_stationary_roles:4});
 }
 const {L,n,U}=cyclic[1],stripe=col(Array.from({length:n},(_,i)=>Math.floor(i/2/L)%2?-1:1));
 ensure(o.isZero(o.mul(ring(omega,L),stripe))&&!o.equal(o.mul(U,stripe),stripe),'period-two invisible source still evolves');
 return {count_quotients:rows,positive_complete_joint_form_rank_checks:positiveJointForms,period_two_blind_but_dynamical_witness:true};
});

check('one_event_block_increment_recovers_all_nonuniform_source_information',()=>{
 const rows=[];let exactInverses=0,identities=0,bounds=0;
 for(const {L,n,id,A,U,F,Lap,remain,mu}of cyclic){
  if(L>4)continue;
  const decoder=o.mul(o.inverse(o.add(Lap,A)),o.dagger(F));
  eq(o.mul(decoder,F),remain,'exact native increment decoder');eq(o.mul(decoder,A),o.zeros(n),'uniform increment blindness');exactInverses+=2;
  const step=o.sub(id,o.scale(Lap,'1/4')),q2=mu?f(1).sub(mu.div(4)):f(0),seed=col(Array.from({length:n},(_,i)=>(i%7)-3)),phi=o.mul(remain,seed),reading=o.mul(F,phi);
  let sum=o.zeros(n),term=id;const initial=o.energy(phi),readNorm=o.energy(reading);
  for(let N=0;N<=5;N++){
   sum=o.add(sum,term);term=o.mul(step,term);
   const finite=o.scale(o.mul(o.mul(sum,o.dagger(F)),remain),'1/4');
   eq(o.mul(finite,F),o.sub(remain,o.mul(term,remain)),'finite native decoder remainder identity');identities++;
   const err=o.energy(o.sub(o.mul(finite,reading),phi));
   ensure(err.le(q2.pow(N+1).mul(initial)),'derived finite decoder error bound');bounds++;
   if(mu){ensure(err.le(q2.pow(N+1).mul(readNorm).div(mu)),'error bound from observed increment alone');bounds++;}
  }
  rows.push({period:L,increment_observer_rank:n-2,unobservable_roles:2,finite_reconstruction_depths:6});
 }
 return {count_quotients:rows,exact_inverse_and_blindness_identities:exactInverses,complete_finite_decoder_remainder_identities:identities,
  certified_error_bound_checks:bounds,classical_spectral_or_Fourier_inverse_used:false};
});

check('native_address_propagation_has_a_finite_count_cone',()=>{
 let sourceChecks=0,jointChecks=0,extremes=0;
 let v=unit(e0),z=stack(new Map(),unit(e0));
 for(let n=0;n<=6;n++){
  ensure([...v.keys()].every(k=>xy(k).every(a=>Math.abs(a)<=n)),'source count cone');sourceChecks++;
  ensure([...z.keys()].every(k=>xy(k).every(a=>Math.abs(a)<=n+1)),'joint count cone with fixed shear collar');jointChecks++;
  ensure(!o.isZero(value(v,n,n)),'outer diagonal address is attained');extremes++;
  if(n<6){v=compose(forward,v);z=compose(joint,z);}
 }
 return {finite_source_support_checks:sourceChecks,joint_support_checks:jointChecks,attained_outer_addresses:extremes,
  source_radius_per_four_event_block:1,physical_speed_c_in_units_selected:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cp=add(hh,kk),cm=add(hh,sc(-1,kk));
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['source_census_square',mul(cp,cp),sc(2,ii)],['reverse_role_square',mul(cm,cm),sc(2,ii)],['native_turn_square',mul(rr,rr),sc(-1,ii)],
 ['turn_reverses_census',mul(rr,cp),sc(-1,cm)],['census_then_turn',mul(cp,rr),cm],
 ['turn_reverses_reflection',mul(rr,cm),cp],['reflection_then_turn',mul(cm,rr),sc(-1,cp)],
 ['plus_cut_square',mul(pp,pp),pp],['role_cuts_disjoint',mul(pp,qq),sc(0,ii)]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const termwiseZero=a=>[...a.values()].every(o.isZero),L2=cyclic[1];
const negatives={
 a_nonzero_off_diagonal_drive_forces_an_irreducible_interaction:termwiseZero(minus(compose(D,joint),compose(lift(I,forward),D))),
 one_shear_cannot_unmix_both_native_directions:termwiseZero(minus(compose(D,cy),compose(lift(I,vy),D))),
 y_direction_uses_the_same_signed_loop_drive_as_x:!termwiseZero(minus(comm(vy,shear),compose(dy,omega))),
 complete_joint_energy_requires_discarding_loop_blind_information:psd(ring(fullGram,2,4,4))===16,
 original_b_phi_norm_is_always_the_complete_conserved_form:(()=>{const a=stack(new Map(),unit(e0));return !energy(compose(cx,a)).eq(energy(a));})(),
 mixed_native_wave_term_can_be_omitted:!termwiseZero(minus(wave,scale(sumlap,'1/2'))),
 two_address_labels_force_a_rank_two_quadratic_geometry:o.rank(o.matrix([['1/2','1/2'],['1/2','1/2']]))===1,
 opposite_shift_wave_is_second_order:!termwiseZero(minus(reduce(wave,1,-1),scale(lapx,'1/2'))),
 equal_shift_curvature_is_nonzero:termwiseZero(reduce(omega,1,1)),
 zero_curvature_response_implies_stationarity:checks.find(r=>r.name==='stationary_and_loop_blind_sectors_are_distinct').period_two_blind_but_dynamical_witness,
 stationary_source_has_all_three_loop_blind_sectors:o.rank(cyclic[2].F)===cyclic[2].n-2&&2*(3*3-2)>2,
 finite_pairing_stationary_classification_covers_unbounded_affine_profiles:o.isZero(at(increment,(x,y)=>o.scale(e0,x-y),0,0))&&!o.isZero(at(dx,(x,y)=>o.scale(e0,x-y),0,0)),
 block_increment_observer_forgets_the_entire_curvature_kernel:o.rank(cyclic[2].F)>o.rank(ring(omega,3)),
 increment_reading_recovers_the_uniform_offset:o.isZero(o.mul(L2.F,L2.A)),
 finite_decoder_is_exact_after_a_generic_single_term:(()=>{const step=o.sub(cyclic[2].id,o.scale(cyclic[2].Lap,'1/4'));return !o.isZero(o.mul(step,cyclic[2].remain));})(),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r31.native-propagation-geometry.v1',status:'PASS_R31_NATIVE_PROPAGATION_GEOMETRY',
 input_sha256:inputHash,r30_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_coupling_unmixed_without_information_loss:true,native_two_direction_wave_and_quartic_sector_derived:true,
 native_positive_propagation_factorization_derived:true,complete_increment_observer_and_error_bound_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,
 physical_metric_c_alpha_derived:false,physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,
 primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native results: simultaneous local unmixing of loop-driven fields, full source energy and bounded drive, two-direction wave/count cone, positive factorization and stationary sector, directional quartic response and rank-one count jet, finite-count gap bound and complete block-increment reconstruction. Explicit mathematical targets; no physical force/EM/c/alpha/mass/metric selection.'},null,2)+'\n');
