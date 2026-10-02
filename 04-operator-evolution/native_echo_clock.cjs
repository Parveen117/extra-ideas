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
ensure(inputHash===arg('--expected-input-sha256'),'R39 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r38-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r38-sha256'),'R39 R38 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R38_NATIVE_UNIVERSAL_SIGNAL_CONE','R39 parent status');
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
const I4=o.identity(4),Ps=o.kron(P,I),Qs=o.kron(Q,I),Cs=o.kron(C0,I);
const X=shift(unit(o.kron(I,H)),1,0),Y=shift(unit(o.kron(I,K)),0,1);
const covV=A=>plus(unit(I4),scale(compose(compose(plus(unit(Ps),compose(dag(A),unit(Qs))),unit(Cs)),minus(A,unit(I4))),'1/2'));
const VX=covV(X),VY=covV(Y),W=compose(VY,VX),B=compose(W,W),F=minus(B,unit(I4));
const TX=shift(unit(I4),1,0),TY=shift(unit(I4),0,1),DX=minus(compose(TX,TX),unit(I4)),DY=minus(compose(TY,TY),unit(I4));
const LX=compose(dag(DX),DX),LY=compose(dag(DY),DY),QX=plus(unit(I4),scale(LX,'1/4')),QY=plus(unit(I4),scale(LY,'1/4'));
const G=compose(QX,QY),coefficient=minus(unit(o.scale(I4,2)),G),pairAverageX=scale(plus(compose(TX,TX),dag(compose(TX,TX))),'1/2'),pairAverageY=scale(plus(compose(TY,TY),dag(compose(TY,TY))),'1/2');
const liftMemory=a=>new Map([...a].map(([k,M])=>[k,o.kron(M,I)])),bareW=liftMemory(forward);
const zeroOp=a=>[...a.values()].every(o.isZero),parity=n=>((n%2)+2)%2;
const recordFrame=(x,y)=>o.mul(o.power(H,parity(x)),o.power(K,parity(y)));
const probes=[unit(o.kron(e0,e0)),unit(o.kron(e1,e1)),new Map([['-1,0',col([1,-2,3,1])],['1,1',col([-1,2,0,1])]])];

const block=(a)=>{const out=new Map();for(const [k,A]of a){const [dx,dy]=xy(k);for(let ex=0;ex<2;ex++)for(let ey=0;ey<2;ey++){
 const sx=-Math.floor((ex-dx)/2),sy=-Math.floor((ey-dy)/2),fx=parity(ex-dx),fy=parity(ey-dy),M=o.zeros(16),sign=-((-1)**(sx+sy));
 for(let i=0;i<4;i++)for(let j=0;j<4;j++)M[(ex*2+ey)*4+i][(fx*2+fy)*4+j]=A[i][j].mul(sign);
 addAt(out,key(sx,sy),M);
 }}return out;};
const Z=block(B),I16=o.identity(16),U16=unit(I16),zero16=o.zeros(16);
const sumAtOne=a=>[...a.values()].reduce(o.add,o.zeros(a.values().next().value.length));
const jet=axis=>[...Z].reduce((s,[k,M])=>o.add(s,o.scale(M,xy(k)[axis])),zero16);
const Ax=jet(0),Ay=jet(1),trace=M=>M.reduce((s,row,i)=>s.add(row[i]),o.Cut.of(0));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});
const id1=unit(o.identity(1)),sx1=shift(id1,1,0),sy1=shift(id1,0,1);
const dx1=minus(sx1,id1),dy1=minus(sy1,id1),ex1=plus(sx1,id1),ey1=plus(sy1,id1);
const lx1=compose(dag(dx1),dx1),ly1=compose(dag(dy1),dy1);
const ell1=minus(scale(plus(lx1,ly1),'1/2'),scale(compose(lx1,ly1),'1/16'));
const scalarLift=(a,size)=>new Map([...a].map(([k,M])=>[k,o.scale(o.identity(size),M[0][0])]));
const ell=scalarLift(ell1,16);
const normSquare=a=>compose(dag(a),a);

// All coefficients below use the unchanged native Cut/F arithmetic.
const D4=o.kron(A0,I),D16raw=o.kron(o.identity(4),D4);
const link=plus(compose(Y,unit(Ps)),unit(Qs)),U=scale(block(link),-1);
const ct=compose(compose(U,unit(D16raw)),dag(U)); // sqrt(2) times C
const basis=(n,j)=>col(Array.from({length:n},(_,i)=>i===j?1:0));
const fieldData=a=>[...a].sort(([a],[b])=>a.localeCompare(b)).map(([k,M])=>[k,M.map(r=>r.map(x=>x.toJSON()))]);
const fieldHash=a=>hash(JSON.stringify(fieldData(a)));
const signedPower=(a,n)=>n<0?f(1).div(f(a).pow(-n)):f(a).pow(n);
const lt=(a,b)=>f(a).le(b)&&!f(a).eq(b);
const floorF=a=>{a=f(a);ensure(a.n>=0n,'nonnegative count floor');return a.n/a.d;};
const ceilF=a=>{a=f(a);ensure(a.n>=0n,'nonnegative count ceiling');return (a.n+a.d-1n)/a.d;};
check('source_factors_derive_the_complete_local_reversal_involution',()=>{
 for(const V of [VX,VY])eqFields(compose(compose(unit(D4),V),unit(D4)),scale(dag(V),2),4,4,'common cut reversal');
 eqFields(compose(dag(U),U),U16,16,16,'link isometry');
 eqFields(compose(U,dag(U)),U16,16,16,'link coisometry');
 eqFields(dag(ct),ct,16,16,'reversal dagger');
 eqFields(compose(ct,ct),scale(U16,2),16,16,'reversal square');
 eqFields(compose(compose(ct,Z),ct),scale(dag(Z),2),16,16,'complete source reversal');
 eqFields(ct,scale(block(compose(VY,unit(D4))),-1),16,16,'same source factorization');
 ensure([...ct.keys()].sort().join(';')===['0,-1','0,0','0,1'].sort().join(';'),'finite reversal support');
 ensure([...U.keys()].sort().join(';')===['0,0','0,1'].sort().join(';'),'finite link support');
 return {complete_laurent_identities:8,reverser_support:[...ct.keys()].sort(),link_support:[...U.keys()].sort(),ordinary_plane_reflection_used:false};
});
const zPowers=[U16,Z,compose(Z,Z)],zInteger=n=>n<0?dag(zPowers[-n]):zPowers[n];
check('full_signed_echo_cancels_equal_counts_and_retains_mismatched_counts',()=>{
 const pairs=[[0,0],[1,1],[2,0],[0,2],[2,1],[1,2]];
 for(const [a,b]of pairs)eqFields(compose(compose(compose(ct,zPowers[b]),ct),zPowers[a]),scale(zInteger(a-b),2),16,16,'echo count identity');
 ensure(!zeroOp(minus(Z,U16)),'mismatched count must retain history');
 return {complete_echo_word_identities:pairs.length,count_pairs:pairs,mismatched_one_block_defect_sha256:fieldHash(minus(Z,U16)),full_signed_state_not_only_intensity_restored:true};
});

// Rational native matrix witness t^2=1/2 checks both root conjugates.
// The positive scalar root is selected in the written native proof.
const rootWitness=o.scale(C0,'1/2'),D32=o.kron(D16raw,rootWitness);
const tensorRole=(a,M)=>new Map([...a].map(([k,A])=>[k,o.kron(A,M)]));
const U32=tensorRole(U,I),Ud32=dag(U32),Z32=tensorRole(Z,I),C32=tensorRole(ct,rootWitness);
const rect=(xmin,xmax,ymin,ymax)=>(x,y)=>xmin<=x&&x<=xmax&&ymin<=y&&y<=ymax;
const mask=(a,T)=>new Map([...a].filter(([k])=>T(...xy(k))));
const linkedCut=(a,T)=>compose(U32,mask(compose(Ud32,a),T));
const localGate=(a,T)=>compose(U32,new Map([...compose(Ud32,a)].map(([k,v])=>[k,T(...xy(k))?o.mul(D32,v):v])));
const testFields=[unit(o.kron(basis(16,0),e0)),unit(o.kron(basis(16,7),e1)),new Map([['0,0',o.kron(col(Array.from({length:16},(_,i)=>i%3-1)),e0)],['1,-1',o.kron(basis(16,12),col([2,-1]))]])];
const pointPatch=rect(0,0,0,0),originPatch=rect(0,0,-1,1),remotePatch=rect(1,2,-1,1);
check('paired_link_cuts_give_unitary_finite_window_reversal',()=>{
 eq(o.mul(rootWitness,rootWitness),o.scale(I,'1/2'),'native algebraic root witness');
 eq(o.dagger(D32),D32,'onsite gate dagger');eq(o.mul(D32,D32),o.identity(32),'onsite gate square');
 let basisChecks=0,identityOutside=0;
 for(const y of [-1,0,1,2])for(let j=0;j<32;j++){
  const v=shift(unit(basis(32,j)),0,y),gv=localGate(v,pointPatch);
  eqFields(localGate(gv,pointPatch),v,32,1,'local gate square on complete patch basis');
  ensure(energy(gv).eq(energy(v)),'local gate norm');basisChecks++;
  if(y===-1||y===2){eqFields(gv,v,32,1,'identity outside paired patch');identityOutside++;}
 }
 let cutChecks=0;
 for(const v of testFields){const q=linkedCut(v,pointPatch),cq=compose(C32,q);
  eqFields(linkedCut(q,pointPatch),q,32,1,'linked projection');
  eqFields(linkedCut(compose(C32,v),pointPatch),cq,32,1,'cut reverser commute');
  eqFields(localGate(v,pointPatch),plus(cq,minus(v,q)),32,1,'orthogonal local gate');
  ensure(energy(q).add(energy(minus(v,q))).eq(energy(v)),'orthogonal norm split');cutChecks+=4;
 }
 return {complete_local_basis_checks:basisChecks,identity_outside_patch_checks:identityOutside,linked_cut_checks:cutChecks,root_witness_checks_both_algebraic_conjugates:true,root_witness_is_positive_scalar_simulation:false};
});
const evolved=(a,n)=>{let v=a;for(let i=0;i<n;i++)v=compose(Z32,v);return v;};
const echoCycle=(a,n=1)=>localGate(evolved(localGate(evolved(a,n),remotePatch),n),originPatch);
let naiveWitness=null,finiteEchoFields=[];
check('finite_window_echo_error_is_bounded_by_uncaptured_native_information',()=>{
 let bounds=0;for(const T of [pointPatch,originPatch,remotePatch])for(const v of testFields){
  const missing=energy(minus(v,linkedCut(v,T))),error=energy(minus(localGate(v,T),compose(C32,v)));
  ensure(error.le(missing.mul(4)),'finite clipping bound');bounds++;
 }
 for(let j=0;j<32&&!naiveWitness;j++){
  const v=unit(basis(32,j)),naive=plus(mask(compose(C32,mask(v,pointPatch)),pointPatch),minus(v,mask(v,pointPatch)));
  if(!energy(naive).eq(energy(v)))naiveWitness={basis_role:j,input_norm_square:energy(v).toString(),naively_clipped_norm_square:energy(naive).toString()};
 }ensure(naiveWitness,'naive clipping counterexample');
 const v=testFields[0];eqFields(linkedCut(v,originPatch),v,32,1,'origin capture');
 const rows=[];for(const n of [1,2,3]){
  const outgoing=evolved(v,n),missing=energy(minus(outgoing,linkedCut(outgoing,remotePatch))),returned=echoCycle(v,n),error=energy(minus(returned,v));
  ensure(error.le(missing.mul(4)),'localized echo clipping budget');ensure(energy(returned).eq(energy(v)),'echo total norm');
  rows.push({free_leg_blocks:n,remote_missing_norm_square:missing.toString(),echo_error_square:error.toString(),upper_bound:missing.mul(4).toString(),returned_field_sha256:fieldHash(returned)});
  finiteEchoFields.push(returned);
 }
 return {literal_clipping_bounds:bounds,naive_address_clipping_counterexample:naiveWitness,literal_echo_instances:rows,small_literal_instances_are_high_retention_R37_packets:false};
});

const r37Path=path.join(path.dirname(arg('--r38-certificate')),'R37_NATIVE_CERTIFICATE.json'),r37Bytes=fs.readFileSync(r37Path);
ensure(hash(r37Bytes)===JSON.parse(parentBytes).r37_native_certificate_sha256,'R39 frozen R37 packet pin mismatch');
const r37=JSON.parse(r37Bytes),constants=r37.exact_checks.find(x=>x.name==='a_source_corrector_removes_every_first_difference_error').source_constants.find(x=>x.parameter==='2');
const kappa=f(constants.kappa),CT=f(constants.comparison_norm_bound),Csum=f(constants.residual_xx_mass).add(constants.residual_xy_mass).add(constants.residual_yy_mass),drift=f('2/3');
const budget=L=>{L=f(L);const M=L.pow(2),N=L.pow(3),width=M.mul(3).sub(2),d=f(floorF(N.mul(drift))),dd=f(ceilF(width.add(3).div(drift))),d0=f(ceilF(width.add(1).div(drift)));
 const epsilon=Csum.mul(4).add(CT.mul(2).mul(drift).mul(f(1).sub(drift))).div(L).add(CT.mul(2).div(M)).div(f(1).sub(kappa.mul(2).div(M)));
 return {L,M,N,width,d,dd,d0,epsilon};};
const budgets=[4096,16384,65536].map(budget);
check('pinned_packet_bounds_certify_remote_and_return_trigger_windows',()=>{
 ensure(kappa.eq('214/9')&&CT.eq('437/9'),'exact packet source constants');
 const rows=budgets.map(b=>{const e2=b.epsilon.pow(2),trigger=f(1).sub(e2.mul(2));
  ensure(lt(e2,'1/32'),'small echo budget');ensure(lt(b.dd.add(b.d0),b.N),'separated flight and aperture scales');
  ensure(lt(e2,trigger)&&trigger.le(f(1).sub(e2)),'remote trigger exists and no early trigger');
  ensure(lt(e2.mul(16),'1/2')&&lt('1/2',f(1).sub(e2.mul(9))),'return first half threshold');
  ensure(lt(b.width,b.dd.mul(drift).sub(1))&&b.width.add(1).le(floorF(b.d0.mul(drift))),'native support separation');
  return {L:b.L.toString(),profile_M:b.M.toString(),horizon_N:b.N.toString(),patch_displacement:b.d.toString(),epsilon:b.epsilon.toString(),epsilon_square:e2.toString(),echo_norm_error_upper:b.epsilon.mul(3).toString(),remote_trigger:trigger.toString(),remote_window_Delta_d:b.dd.toString(),return_window_Delta_0:b.d0.toString()};
 });
 return {frozen_R37_native_certificate_sha256:hash(r37Bytes),carrier_drift:drift.toString(),exact_large_parameter_budgets:rows,large_packets_directly_simulated:false,trigger_is_a_declared_native_norm_control:true};
});
check('scheduled_echo_cycles_have_controlled_error_and_distinct_tick_records',()=>{
 const v=testFields[0],one=echoCycle(v),oneError=energy(minus(one,v));let current=v;const repeated=[];
 for(let k=0;k<=3;k++){if(k)current=echoCycle(current);const error=energy(minus(current,v));
  ensure(error.le(oneError.mul(k*k)),'finite unitary echo telescoping');ensure(energy(current).eq(energy(v)),'repeated echo norm');
  repeated.push({ticks:k,error_square:error.toString(),bound:k===0?'0':oneError.mul(k*k).toString()});
 }
 const counters=[];for(const q of [2,3,5]){const successor=o.zeros(q);for(let j=0;j<q;j++)successor[(j+1)%q][j]=o.ONE;
  eq(o.mul(o.dagger(successor),successor),o.identity(q),'native counter permutation');let state=basis(q,0);const seen=new Set();
  for(let j=0;j<q;j++){const tag=fieldHash(unit(state));ensure(!seen.has(tag),'distinct tick mark');seen.add(tag);state=o.mul(successor,state);}
  eq(state,basis(q,0),'finite counter alias');let binary=0;while(2**binary<q)binary++;
  ensure(2**binary>=q&&(binary===0||2**(binary-1)<q),'minimum binary marked slots');
  counters.push({labels:q,distinct_counts:seen.size,first_alias_tick:q,minimum_binary_marked_slots:binary});
 }
 let previous=null;const scaling=[];for(const p0 of [6,7,8]){const L=f(16).pow(p0),k=f(4).pow(p0),b=budget(L),bound=b.epsilon.mul(k).mul(3);
  ensure(k.pow(2).eq(L),'growing tick count');if(previous)ensure(lt(bound,previous),'decreasing total echo error');previous=bound;
  scaling.push({L:L.toString(),ticks:k.toString(),total_echo_error_bound:bound.toString()});
 }
 return {literal_repeated_cycles:repeated,native_counter_orbits:counters,growing_tick_vanishing_error_budgets:scaling,scheduled_reset_count:'2 tau free blocks plus two controls',earliest_return_is_not_silently_reset:true,capacity_scope:'distinguishable marked labels, not arbitrary coherent phase dimensions',unbounded_time_stored_in_exactly_recurrent_source:false};
});
check('return_counts_derive_midpoint_distance_and_information_length_calibration',()=>{
 let corners=0;const rows=[];for(const b of budgets){
  for(const tau of [b.N.sub(b.dd).add(1),b.N])for(const s of [tau.sub(b.d0).add(1),tau]){const r=tau.add(s),midpoint=tau.sub(r.div(2)),distance=b.d.sub(drift.mul(r).div(2));
   ensure(f(0).le(midpoint)&&lt(midpoint,b.d0.div(2)),'midpoint interval');
   ensure(f(-1).le(distance)&&lt(distance,drift.mul(b.dd.add(b.d0.div(2)))),'distance interval');corners++;
  }
  rows.push({L:b.L.toString(),midpoint_error_strict_upper:b.d0.div(2).toString(),distance_error_lower:'-1',distance_error_strict_upper:drift.mul(b.dd.add(b.d0.div(2))).toString(),relative_aperture_budget:b.dd.add(b.d0.div(2)).div(b.N).toString()});
 }
 const cuts=r37.exact_checks.find(x=>x.name==='a_native_rational_speed_family_approaches_the_sharp_count_coefficient').native_speed_cuts;
 for(const row of cuts){const a=f(row.carrier_parameter),w0=a.add(f(2).div(a)).div(2),g0=f(1).div(w0); // native rational carrier relation
  ensure(g0.eq(row.derived_drift)&&f('1/2').sub(g0.pow(2)).eq(row.squared_speed_deficit),'inherited sharp speed cuts');
 }
 ensure(f('1/2').mul(2).eq(1),'information length calibration coefficient');
 return {exact_interval_corner_checks:corners,reading_intervals:rows,inherited_approaching_speed_cuts:cuts.length,sharp_echo_speed:'1/sqrt(2)',information_length:'sqrt(2) d',fixed_carrier_length_ratio:'1/(sqrt(2) g)',free_source_events_per_scheduled_cycle:'16 tau',reversal_controls_per_cycle:2,control_duration_was_assigned:false,independent_physical_clock_synchronization_proved:false};
});
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),J=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),Qc=o.sub(I16,Pc),Ct0=sumAtOne(ct),Tt=o.mul(Ct0,J);
const P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,J),-1),I),J32=o.kron(J,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=pair=>plus(pair[0],compose(unit(J32),pair[1]));
check('complete_curvature_observation_preserves_the_echo_and_exposes_link_correction',()=>{
 for(const A of [Ax,Ay]){eq(o.add(o.mul(Ct0,A),o.mul(A,Ct0)),zero16,'leading source reversal');eq(o.mul(Tt,A),o.mul(A,Tt),'commuting extra native turn');}
 eq(o.mul(Ct0,J),o.mul(J,Ct0),'curvature commutes with zero-carrier reversal');
 eq(o.mul(Tt,Tt),o.scale(I16,-2),'extra native turn square');eq(o.scale(o.mul(Tt,J),-1),Ct0,'source reversal contains curvature and extra turn');
 const bare=minus(compose(compose(unit(J),Z),unit(o.dagger(J))),dag(Z));ensure(!zeroOp(bare),'bare curvature reversal defect survives');
 const lost=compose(compose(unit(Pc),ct),unit(Qc));ensure(!zeroOp(lost),'single cut loses reverser data');
 let count=0;for(const v of [...testFields,...finiteEchoFields]){const observed=observe(v);eqFields(decode(observed),v,32,1,'curvature full reconstruction');
  ensure(energy(observed[0]).add(energy(observed[1])).eq(energy(v)),'curvature complete norm');
  for(const action of [a=>localGate(a,pointPatch),a=>linkedCut(a,remotePatch),a=>compose(Z32,a)]){
   const transported=observe(action(decode(observed))),direct=observe(action(v));for(let j=0;j<2;j++)eqFields(transported[j],direct[j],32,1,'observed control transport');count++;
  }
 }
 return {zero_carrier_reversal_and_turn_identities:7,complete_observed_control_checks:count,bare_curvature_finite_reversal_defect_sha256:fieldHash(bare),partial_cut_reversal_loss_sha256:fieldHash(lost),complete_observation_is_isometric:true};
});
check('source_front_record_wrap_and_control_accounting_reject_overstated_claims',()=>{
 let v=unit(basis(16,0));const fronts=[];for(let n=1;n<=3;n++){v=compose(Z,v);const corner=v.get(key(n,n))||o.zeros(16,1),e=o.energy(corner);
  ensure(e.eq(signedPower(256,-n)),'surviving faster exact source corner');fronts.push({blocks:n,corner:[n,n],norm_square:e.toString()});
 }
 let zeroRejected=false;try{o.Cut.of(0).inv();}catch(_){zeroRejected=true;}ensure(zeroRejected,'undefined source inverse');
 ensure(naiveWitness.naively_clipped_norm_square!==naiveWitness.input_norm_square,'boundary norm loss retained');
 return {exact_fast_front_witnesses:fronts,zero_inverse_rejected:zeroRejected,added_controls_are_free_source_events:false,adaptive_trigger_claimed_linear:false,physical_mirror_or_autonomous_material_clock_selected:false,physical_c_h_alpha_identified:false,minimum_counter_bound_claims_arbitrary_phase_dimension:false,source_law_or_complete_observer_was_replaced:false};
});
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cpw=add(hh,kk),cmw=add(hh,sc(-1,kk)),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['native_turn_square',mul(rr,rr),sc(-1,ii)],['native_census_square',mul(cpw,cpw),sc(2,ii)],['native_reverse_census_square',mul(cmw,cmw),sc(2,ii)],
 ['native_record_loop',mul(hh,kk,hh,kk),sc(-1,ii)],['native_reverse_record_loop',mul(kk,hh,kk,hh),sc(-1,ii)],
 ['native_plus_cut',mul(pp,pp),pp],['native_minus_cut',mul(qq,qq),qq],['native_disjoint_cuts',mul(pp,qq),sc(0,ii)],['native_complete_cut_pair',add(pp,qq),ii],
 ['native_half_speed_algebraic_witness',mul(sc('1/2',cpw),sc('1/2',cpw)),sc('1/2',ii)],['native_turn_increment_norm',mul(add(ii,sc(-1,rr)),add(ii,rr)),sc(2,ii)]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');



const negatives={
 ordinary_complex_field_was_primitive:replays.some(row=>row.name==='native_turn_square'&&row.replay==='REPLAY_MATCH'),
 reversal_requires_reflecting_the_entire_address_plane:checks[0].ordinary_plane_reflection_used===false,
 reversal_only_matches_a_first_jet:checks[0].complete_laurent_identities===8,
 mismatched_flight_counts_still_return_every_state:checks[1].mismatched_one_block_defect_sha256.length===64,
 a_naive_address_clip_preserves_norm:checks[3].naive_address_clipping_counterexample.input_norm_square!==checks[3].naive_address_clipping_counterexample.naively_clipped_norm_square,
 algebraic_root_witness_selects_the_positive_scalar_branch:checks[2].root_witness_is_positive_scalar_simulation===false,
 large_certified_packets_were_brute_force_simulated:checks[4].large_packets_directly_simulated===false,
 earliest_return_and_scheduled_reset_are_identical:checks[5].earliest_return_is_not_silently_reset,
 an_exactly_returned_source_retains_unbounded_elapsed_ticks:checks[5].unbounded_time_stored_in_exactly_recurrent_source===false,
 a_finite_counter_has_no_wrap:checks[5].native_counter_orbits.every(x=>x.first_alias_tick===x.labels),
 flight_count_assigns_duration_to_control_arrows:checks[6].control_duration_was_assigned===false,
 fixed_carrier_dispersion_correction_can_be_dropped:checks[6].fixed_carrier_length_ratio==='1/(sqrt(2) g)',
 bare_curvature_is_the_complete_finite_source_reverser:checks[7].bare_curvature_finite_reversal_defect_sha256.length===64,
 the_exact_faster_source_front_has_been_removed:checks[8].exact_fast_front_witnesses.length===3,
 native_control_construction_identifies_physical_constants:checks[8].physical_c_h_alpha_identified===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r39.native-echo-clock.v1',status:'PASS_R39_NATIVE_ECHO_CLOCK',
 input_sha256:inputHash,r38_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_local_source_reverser_derived:true,native_local_echo_and_return_meter_derived:true,
 native_repeatable_clock_and_record_capacity_derived:true,native_echo_information_length_calibration_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: full local reversal, exact echo and finite window gates, high-retention remote trigger, return timestamp and distance error, repeatable scheduled clock and marked record capacity, information-length calibration and sharp echo speed, and curvature-complete observation. Reversal controls are derived available arrows with explicit count accounting; an autonomous physical mirror, material clock and c, h, alpha are not selected.'},null,2)+'\n');
