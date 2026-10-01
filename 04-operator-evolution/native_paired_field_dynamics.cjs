#!/usr/bin/env node
'use strict';
// Application identities only. Native arithmetic and word replay remain
// in the unchanged canonical Recognition Kernel engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const arg=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(arg('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const ensure=(v,msg)=>{if(!v)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
const input=JSON.parse(fs.readFileSync(arg('--input'),'utf8')),inputHash=p.digest(input);
ensure(inputHash===arg('--expected-input-sha256'),'R29 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r28-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r28-sha256'),'R29 R28 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R28_NATIVE_CUT_CURVATURE_PROPAGATION','R29 parent status');
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const I=o.identity(2),H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange),R=o.mul(K,H);
const P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),C0=o.add(H,K),F=o.commutator(H,K),E01=o.mul(P,K),E10=o.mul(Q,K);
const col=a=>o.matrix(a.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const unit=X=>new Map([[0,X]]);
const addAt=(a,x,M)=>{const v=o.add(a.get(x)||o.zeros(M.length,M[0].length),M);if(o.isZero(v))a.delete(x);else a.set(x,v);};
const plus=(a,b)=>{const c=new Map(a);for(const [x,M]of b)addAt(c,x,M);return c;};
const scale=(a,c)=>new Map([...a].map(([x,M])=>[x,o.scale(M,c)]));
const minus=(a,b)=>plus(a,scale(b,-1)),shift=(a,d)=>new Map([...a].map(([x,M])=>[x+d,M]));
const compose=(a,b)=>{const c=new Map();for(const [i,A]of a)for(const [j,B]of b)addAt(c,i+j,o.mul(A,B));return c;};
const dag=a=>new Map([...a].map(([x,M])=>[-x,o.dagger(M)]));
const eqFields=(a,b,r,c,msg)=>{for(const x of new Set([...a.keys(),...b.keys()]))eq(a.get(x)||o.zeros(r,c),b.get(x)||o.zeros(r,c),msg+' at '+x);};
const energy=a=>[...a.values()].reduce((s,M)=>s.add(o.energy(M)),f(0));
const dot=(a,b)=>[...a].reduce((s,[x,v])=>s.add(o.mul(o.dagger(v),b.get(x)||o.zeros(v.length,1))[0][0].rad),f(0));
const sum=a=>[...a.values()].reduce((s,v)=>o.add(s,v),o.zeros(2,1));
const value=(a,y)=>a.get(y)||o.zeros(2,1);
const lift=(A,a)=>new Map([...a].map(([x,M])=>[x,o.kron(A,M)]));
const stack=(a,b)=>{const out=new Map();for(const x of new Set([...a.keys(),...b.keys()]))out.set(x,[...(a.get(x)||o.zeros(2)),...(b.get(x)||o.zeros(2))]);return out;};
const hstack=(a,b)=>{const out=new Map();for(const x of new Set([...a.keys(),...b.keys()])){const A=a.get(x)||o.zeros(2),B=b.get(x)||o.zeros(2);out.set(x,A.map((row,i)=>row.concat(B[i])));}return out;};
const grad=new Map([[1,I],[0,o.scale(I,-1)]]),gdag=dag(grad),lap=compose(gdag,grad);
const N=new Map([[0,o.scale(o.mul(P,C0),'1/2')],[-1,o.scale(o.mul(Q,C0),'1/2')]]);
const V=plus(unit(I),compose(N,grad)),Vdag=dag(V);
const raw=new Map([[1,o.mul(P,C0)],[-1,o.mul(Q,C0)]]);
const sourceV=new Map([...compose(raw,raw)].map(([x,M])=>[x/2,o.scale(M,'1/2')]));
const M=plus(plus(lift(P,minus(unit(I),scale(lap,'1/2'))),lift(E01,scale(gdag,'-1/2'))),plus(lift(E10,grad),lift(Q,unit(I))));
const Mi=plus(plus(lift(P,unit(I)),lift(E01,scale(gdag,'1/2'))),plus(lift(E10,scale(grad,-1)),lift(Q,minus(unit(I),scale(lap,'1/2')))));
const G=plus(plus(lift(P,unit(I)),lift(Q,unit(o.scale(I,'1/2')))),plus(lift(E01,scale(gdag,'1/4')),lift(E10,scale(grad,'1/4'))));
const include=stack(N,unit(I)),residue=hstack(unit(I),scale(N,-1));
const evolve=(e,b)=>{const bn=plus(b,compose(grad,e)),en=minus(e,scale(compose(gdag,bn),'1/2'));return {e:en,b:bn};};
const defect=(e,b)=>minus(e,compose(N,b));
const fieldEnergy=(e,b)=>energy(e).add(energy(b).div(2)).add(dot(b,compose(grad,e)).div(2));
const squares=(e,b)=>energy(plus(b,scale(compose(grad,e),'1/2'))).div(2).add(energy(e).div(2)).add(energy(plus(shift(e,1),e)).div(8));
const seeds=[unit(e0),unit(e1),new Map([[-2,col([2,-1])],[0,col([-3,4])],[3,col([1,2])]])];
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

check('source_transport_factors_into_a_native_constitutive_operator',()=>{
 eqFields(V,sourceV,2,2,'derive coarse V from two literal source steps');
 eqFields(minus(V,unit(I)),compose(N,grad),2,2,'source factorization');
 eqFields(compose(dag(N),N),unit(o.scale(I,'1/2')),2,2,'constitutive Gram');
 eqFields(compose(N,dag(N)),unit(o.scale(I,'1/2')),2,2,'constitutive co-Gram');
 eqFields(compose(Vdag,V),unit(I),2,2,'source pairing preservation');
 eqFields(plus(V,Vdag),minus(unit(o.scale(I,2)),scale(lap,'1/2')),2,2,'native wave coefficient');
 eqFields(compose(N,Vdag),plus(N,scale(gdag,'1/2')),2,2,'constitutive adjoint identity');
 let norms=0,preparations=0;
 for(const phi of seeds){const b=compose(grad,phi),e=minus(compose(V,phi),phi);
  eqFields(e,compose(N,b),2,1,'field relation from a source potential');ensure(energy(e).eq(energy(b).div(2)),'native norm ratio');norms++;
  eqFields(compose(unit(o.scale(F,'1/2')),compose(unit(R),phi)),phi,2,1,'native curvature preparation covers the source target');preparations++;
 }
 return {Laurent_operator_identities:7,native_norm_ratio_checks:norms,curvature_preparation_checks:preparations,
  constitutive_norm_ratio:'1/2',new_physical_coefficient_supplied:false};
});

check('paired_induction_is_local_reversible_and_obeys_the_native_wave',()=>{
 eqFields(compose(M,Mi),unit(o.identity(4)),4,4,'paired right inverse');eqFields(compose(Mi,M),unit(o.identity(4)),4,4,'paired left inverse');
 const coefficient=lift(I,minus(unit(o.scale(I,2)),scale(lap,'1/2')));
 eqFields(plus(minus(compose(M,M),compose(coefficient,M)),unit(o.identity(4))),new Map(),4,4,'complete paired wave polynomial');
 ensure([...M.keys()].every(x=>Math.abs(x)<=1)&&[...Mi.keys()].every(x=>Math.abs(x)<=1),'finite forward/backward range');
 let matched=0,waves=0;
 for(const phi0 of seeds){let phi=phi0,e=minus(compose(V,phi),phi),b=compose(grad,phi);const es=[e],bs=[b];
  for(let m=0;m<12;m++){
   const next=evolve(e,b);phi=compose(V,phi);
   eqFields(next.e,minus(compose(V,phi),phi),2,1,'direct-source event field');eqFields(next.b,compose(grad,phi),2,1,'direct-source address field');matched+=2;
   e=next.e;b=next.b;es.push(e);bs.push(b);
  }
  for(let m=1;m<12;m++)for(const states of [es,bs]){
   eqFields(plus(minus(states[m+1],scale(states[m],2)),states[m-1]),scale(compose(lap,states[m]),'-1/2'),2,1,'paired field wave');waves++;
  }
 }
 return {Laurent_inverse_and_wave_identities:3,source_field_matches:matched,finite_wave_checks:waves,maximum_pair_event:12,paired_event_range:1};
});

check('compatibility_residue_is_transported_and_selects_the_source_sector',()=>{
 eqFields(compose(M,include),compose(include,V),4,2,'native source graph intertwiner');
 eqFields(compose(residue,M),compose(Vdag,residue),2,4,'opposite-continuation residue intertwiner');
 let norms=0,forced=0,integrals=0;
 const pairs=[...seeds.map(phi=>({b:compose(grad,phi),e:minus(compose(V,phi),phi)})),
  {e:unit(e0),b:new Map()},{e:seeds[2],b:unit(e1)},{e:scale(compose(N,seeds[2]),-1),b:seeds[2]}];
 for(const pair of pairs){let {e,b}=pair,c=defect(e,b);const cenergy=energy(c),total=sum(b);
  for(let m=0;m<12;m++){
   const next=evolve(e,b),cn=defect(next.e,next.b);
   eqFields(cn,compose(Vdag,c),2,1,'residue recurrence');
   eqFields(next.b,plus(compose(V,b),compose(grad,c)),2,1,'explicit residue source');forced+=2;
   ensure(energy(cn).eq(cenergy),'residue norm is retained');norms++;eq(sum(next.b),total,'total address-field invariant');integrals++;
   e=next.e;b=next.b;c=cn;
  }
 }
 return {Laurent_intertwining_identities:2,exact_residue_and_forcing_checks:forced,residue_norm_checks:norms,total_field_checks:integrals,
  incompatible_pairs_included:true,spurious_wave_solutions_promoted_to_source:false};
});

check('paired_field_energy_has_an_exact_positive_conserved_form',()=>{
 eqFields(compose(compose(dag(M),G),M),G,4,4,'full native field Gram conservation');
 let norms=0,sourceNorms=0;
 const pairs=[{e:unit(e0),b:new Map()},{e:seeds[2],b:unit(e1)},...seeds.map(b=>({e:compose(N,b),b}))];
 for(const pair of pairs){let {e,b}=pair;const initial=fieldEnergy(e,b);
  for(let m=0;m<=12;m++){
   ensure(fieldEnergy(e,b).eq(initial),'paired conserved form');ensure(fieldEnergy(e,b).eq(squares(e,b)),'positive square decomposition');norms+=2;
   ensure(f(0).le(initial),'positive native matching form');
   if(m<12)({e,b}=evolve(e,b));
  }
 }
 for(const b0 of seeds){let b=b0;const initial=energy(b);
  for(let m=0;m<=12;m++){
   const e=compose(N,b),g=fieldEnergy(e,b);
   ensure(energy(e).add(energy(b).div(2)).eq(initial),'source balanced energy');
   ensure(g.eq(energy(b).sub(energy(compose(grad,b)).div(8))),'source restricted positive form');
   ensure(initial.div(2).le(g)&&g.le(initial),'source form bounds');sourceNorms+=3;
   b=compose(V,b);
  }
 }
 const off=evolve(unit(e0),new Map()),unbalanced=energy(off.e).add(energy(off.b).div(2));
 ensure(unbalanced.eq('3/2')&&!unbalanced.eq(1),'ordinary field sum fails off source graph');
 return {complete_Laurent_Gram_identities:1,conservation_and_positive_decomposition_checks:norms,source_form_and_bound_checks:sourceNorms,
  rejected_uncoupled_energy:{before:'1',after:unbalanced.toString()}};
});

const current=(b,y)=>{const a=value(b,y),d=value(b,y+1),h=a[0][0].rad.add(a[1][0].rad),k=d[0][0].rad.sub(d[1][0].rad);return h.pow(2).sub(k.pow(2)).div(4).add(h.mul(k).div(2));};
const row=xs=>o.matrix([xs]),outer=(a,b=a)=>o.mul(o.dagger(a),b);
const hp=row([1,1,0,0,0,0]),hc=row([0,0,1,1,0,0]),kc=row([0,0,1,-1,0,0]),kn=row([0,0,0,0,1,-1]);
const Jmatrix=(h,k)=>o.scale(o.add(o.sub(outer(h),outer(k)),o.add(outer(h,k),outer(k,h))),'1/4');
const Jright=Jmatrix(hc,kn),Jleft=Jmatrix(hp,kc);
check('local_field_current_is_a_complete_quadratic_identity',()=>{
 const out=[...o.scale(o.add(hp,kc),'1/2'),...o.scale(o.sub(hc,kn),'1/2')];
 const present=o.matrix([[0,0,1,0,0,0],[0,0,0,1,0,0]]);
 eq(o.sub(outer(out),outer(present)),o.sub(Jleft,Jright),'all coefficients of local quadratic conservation');
 let balances=0,sourced=0;
 for(const b0 of seeds){let b=b0;for(let m=0;m<10;m++){
  const bn=compose(V,b),lo=Math.min(...b.keys())-2,hi=Math.max(...b.keys())+2;
  for(let y=lo;y<=hi;y++){
   ensure(o.energy(value(bn,y)).sub(o.energy(value(b,y))).eq(current(b,y-1).sub(current(b,y))),'exact local source-field current');balances++;
  }b=bn;
 }}
 let e=seeds[2],b=unit(e1);
 for(let m=0;m<8;m++){
  const c=defect(e,b),vb=compose(V,b),dc=compose(grad,c),next=evolve(e,b);
  const lo=Math.min(...b.keys(),...dc.keys())-2,hi=Math.max(...b.keys(),...dc.keys())+2;
  for(let y=lo;y<=hi;y++){
   const cross=o.mul(o.dagger(value(vb,y)),value(dc,y))[0][0].rad.mul(2).add(o.energy(value(dc,y)));
   ensure(o.energy(value(next.b,y)).sub(o.energy(value(b,y))).eq(current(b,y-1).sub(current(b,y)).add(cross)),'local compatibility source');sourced++;
  }({e,b}=next);
 }
 return {complete_native_quadratic_coefficient_identity:true,quadratic_input_coordinates:6,local_balances:balances,local_defect_source_balances:sourced,
  interference_term_retained:true};
});

const recover=(b,qminus=o.zeros(2,1))=>y=>{let v=qminus;for(const [x,u]of b)if(x<=y)v=o.sub(v,u);return v;};
const at=(op,phi,y)=>[...op].reduce((s,[d,X])=>o.add(s,o.mul(X,phi(y-d))),o.zeros(2,1));
check('finite_difference_fields_reconstruct_their_source_without_extra_modes',()=>{
 let reconstructions=0,readouts=0;
 for(const phi of seeds){const b=compose(grad,phi),e=minus(compose(V,phi),phi),decoded=recover(b);
  eq(sum(b),o.zeros(2,1),'finite-potential total residue is zero');
  const lo=Math.min(...phi.keys(),...b.keys())-2,hi=Math.max(...phi.keys(),...b.keys())+2;
  for(let y=lo;y<=hi;y++){
   eq(decoded(y),value(phi,y),'finite potential recovered');reconstructions++;
   eq(at(V,decoded,y),o.add(decoded(y),value(e,y)),'reconstructed native event continuation');readouts++;
  }
 }
 const badE=scale(compose(N,seeds[2]),-1);
 ensure(energy(badE).eq(energy(seeds[2]).div(2))&&!energy(defect(badE,seeds[2])).eq(0),'norm ratio alone is not constitutive compatibility');
 return {finite_native_potential_reconstructions:reconstructions,reconstructed_event_readouts:readouts,
  field_norm_ratio_alone_rejected:true,zero_sum_and_constitutive_constraints_both_retained:true};
});

const ring=(op,L)=>{
 const first=[...op.values()][0],r=first.length,c=first[0].length,A=o.zeros(L*r,L*c);
 for(const [d,X]of op)for(let y=0;y<L;y++){const z=((y-d)%L+L)%L;
  for(let i=0;i<r;i++)for(let j=0;j<c;j++)A[y*r+i][z*c+j]=A[y*r+i][z*c+j].add(X[i][j]);
 }return A;
};
check('cyclic_zero_field_backgrounds_have_exactly_two_native_kernel_roles',()=>{
 const records=[];let recovered=0;
 for(let size=1;size<=8;size++){
  const Dr=ring(grad,size),Nr=ring(N,size),Vr=ring(V,size),Id=o.identity(2*size);
  const fieldMap=ring(stack(minus(V,unit(I)),grad),size),constraint=ring(residue,size);
  const total=o.zeros(2,4*size);for(let y=0;y<size;y++){total[0][4*y+2]=o.ONE;total[1][4*y+3]=o.ONE;}
  const allConstraints=constraint.concat(total);
  ensure(o.rank(fieldMap)===2*size-2,'exact native field-map rank');
  ensure(o.rank(o.sub(Vr,Id))===2*size-2,'all stationary potentials are uniform');
  ensure(o.rank(allConstraints)===2*size+2,'complete paired-field constraint rank');
  eq(o.mul(allConstraints,fieldMap),o.zeros(2*size+2,2*size),'all field images satisfy constraints');
  const phi=col(Array.from({length:2*size},(_,i)=>i%5-2)),b=o.mul(Dr,phi),decoded=[];let q=o.zeros(2,1);
  for(let y=0;y<size;y++){q=o.sub(q,b.slice(2*y,2*y+2));decoded.push(...q);}
  eq(o.mul(Dr,decoded),b,'cyclic finite-sum decoder');eq(o.mul(Nr,b),o.mul(o.sub(Vr,Id),decoded),'cyclic event readout decoder');recovered+=2;
  const background=col(Array.from({length:2*size},(_,i)=>i%2?-2:1));
  eq(o.mul(fieldMap,background),o.zeros(4*size,1),'uniform field blindness');eq(o.mul(Vr,background),background,'uniform paired-event stationarity');
  const curved=o.mul(o.kron(o.identity(size),F),background);ensure(o.energy(curved).eq(o.energy(background).mul(4)),'nonzero background curvature');
  records.push({cyclic_addresses:size,field_rank:2*size-2,uniform_kernel_roles:2,constraint_rank:2*size+2});
 }
 return {cyclic_count_quotients:records,finite_recovery_checks:recovered,physical_periodic_topology_selected:false,
  arbitrary_local_gauge_symmetry_claimed:false};
});

const stepInterface=unit(e0),zeroResidue=compose(grad,unit(e0));
check('boundary_interface_residue_and_native_width_bound_are_exact',()=>{
 const charge=col([2,-1]),profiles=[stepInterface,zeroResidue,...[1,2,3,5,8,13].map(M=>new Map(Array.from({length:M},(_,y)=>[y,o.scale(charge,new o.F(1,M))])))];
 let boundaries=0,continuations=0,bounds=0,norms=0;
 const widthRows=[];
 for(const b0 of profiles){let b=b0;const total=sum(b),initial=energy(b),qminus=col([3,-2]),qplus=o.sub(qminus,total);
  const M0=b.size;ensure(o.energy(total).div(M0).le(initial),'native width lower bound');bounds++;
  if(o.equal(total,charge)){ensure(initial.eq(o.energy(charge).div(M0)),'sharp width-bound family');widthRows.push({width:M0,energy:initial.toString()});}
  for(let m=0;m<=8;m++){
   const phi=recover(b,qminus),lo=Math.min(...b.keys())-2,hi=Math.max(...b.keys())+2,next=compose(V,b),e=compose(N,b);
   eq(phi(lo),qminus,'left end record');eq(phi(hi),qplus,'right end record');eq(sum(b),total,'conserved boundary residue');boundaries+=3;
   ensure(energy(b).eq(initial),'conserved interface field energy');norms++;
   for(let y=lo;y<=hi;y++){
    eq(o.sub(at(V,phi,y),phi(y)),value(e,y),'constant-tail source increment');
    eq(o.sub(at(V,phi,y-1),at(V,phi,y)),value(next,y),'independent potential-to-field continuation');continuations+=2;
   }b=next;
  }
 }
 eqFields(compose(V,stepInterface),new Map([[1,o.scale(e0,'1/2')],[0,o.scale(o.add(e0,e1),'1/2')],[-1,o.scale(e1,'-1/2')]]),2,1,'explicit interface spreading');
 ensure(energy(zeroResidue).eq(2)&&o.isZero(sum(zeroResidue)),'zero-total nonzero wave witness');
 return {interface_profiles:profiles.length,end_record_and_total_checks:boundaries,potential_field_continuation_checks:continuations,
  native_width_bound_checks:bounds,field_energy_checks:norms,sharp_width_family:widthRows,
  nonzero_residue_identified_as_electric_charge:false,universal_mass_from_boundary_residue:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh))),ff=add(hh,kk);
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical source presentation');
const tasks=[
 ['source_H_self_return',mul(hh,hh),ii],['source_K_self_return',mul(kk,kk),ii],
 ['native_source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['two_tag_normalization',mul(ff,ff),sc(2,ii)],['plus_cut_idempotent',mul(pp,pp),pp],
 ['minus_cut_idempotent',mul(qq,qq),qq],['source_cuts_disjoint',mul(pp,qq),sc(0,ii)],
 ['source_role_exchange',mul(kk,pp),mul(qq,kk)],
 ['curvature_preparation_inverse',mul(sc(-1,rr),rr),ii],
 ['constitutive_zero_shift_Gram',mul(sc('1/2',ff),sc('1/2',ff)),sc('1/2',ii)],
];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});
}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const bad=evolve(unit(e0),new Map()),NB=compose(N,stepInterface),currentWithoutCross=o.scale(o.sub(outer(hc),outer(kn)),'1/4');
const negatives={
 constitutive_relation_is_identity:!energy(minus(NB,stepInterface)).eq(0),
 field_norm_ratio_alone_selects_source:!energy(defect(scale(NB,-1),stepInterface)).eq(0),
 every_paired_wave_is_a_source_history:!energy(defect(unit(e0),new Map())).eq(0),
 compatibility_residue_disappears_after_one_step:energy(defect(bad.e,bad.b)).eq(1),
 ordinary_uncoupled_energy_is_conserved_for_all_pairs:energy(bad.e).add(energy(bad.b).div(2)).eq('3/2'),
 link_current_can_omit_interference:!o.equal(Jright,currentWithoutCross),
 zero_field_background_has_zero_local_native_curvature:!o.isZero(o.mul(F,e0)),
 arbitrary_local_potential_change_is_field_invisible:!energy(compose(grad,unit(e0))).eq(0),
 every_finite_address_field_has_zero_boundary_residue:!o.isZero(sum(stepInterface)),
 zero_total_residue_implies_zero_field:o.isZero(sum(zeroResidue))&&!energy(zeroResidue).eq(0),
 fixed_boundary_residue_fixes_a_universal_energy:(()=>{const rows=checks.find(x=>x.name==='boundary_interface_residue_and_native_width_bound_are_exact').sharp_width_family;return !f(rows.find(x=>x.width===1).energy).eq(rows.find(x=>x.width===8).energy);})(),
 native_constitutive_norm_ratio_is_pointwise:!o.energy(value(NB,0)).eq(o.energy(value(stepInterface,0)).div(2)),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r29.native-paired-field-dynamics.v1',status:'PASS_R29_NATIVE_PAIRED_FIELD_DYNAMICS',
 input_sha256:inputHash,r28_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_constitutive_operator_derived:true,paired_induction_and_source_constraint_derived:true,
 positive_field_form_and_local_current_derived:true,stationary_background_and_boundary_residue_classified:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,
 physical_metric_c_alpha_derived:false,physical_electric_charge_or_particle_mass_identified:false,formal_proof_assistant_verified:false,
 scope:'Seven written proofs: source-derived field relation, reversible induction, exact source-sector constraint, positive conserved form, local current and compatibility source, uniform stationary backgrounds and conserved boundary interfaces. The readout and boundary targets remain explicit; physical EM/vacuum/charge/c/alpha/mass are not identified.'},null,2)+'\n');
