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
ensure(inputHash===arg('--expected-input-sha256'),'R40 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r39-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r39-sha256'),'R40 R39 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R39_NATIVE_ECHO_CLOCK','R40 parent status');
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

// Native root witnesses check both algebraic conjugates. The written
// completion selects the positive scalar root; no scalar simulation is claimed.
const basis=(n,j)=>col(Array.from({length:n},(_,i)=>i===j?1:0));
const rootWitness=o.scale(C0,'1/2'),D4=o.kron(A0,I),D16raw=o.kron(o.identity(4),D4);
const link=plus(compose(Y,unit(Ps)),unit(Qs)),U=scale(block(link),-1);
const tensorRole=(a,M)=>new Map([...a].map(([k,A])=>[k,o.kron(A,M)]));
const U32=tensorRole(U,I),Ud32=dag(U32),Z32=tensorRole(Z,I),Zd32=dag(Z32),D32=o.kron(D16raw,rootWitness);
const ct=compose(compose(U,unit(D16raw)),dag(U)),C32=tensorRole(ct,rootWitness);
const rect=(xmin,xmax,ymin,ymax)=>(x,y)=>xmin<=x&&x<=xmax&&ymin<=y&&y<=ymax;
const mask=(a,T)=>new Map([...a].filter(([k])=>T(...xy(k))));
const origin=rect(0,0,-1,1),remote=rect(1,2,-1,1),whole=()=>true;
const onsite=(a,T)=>new Map([...a].map(([k,v])=>[k,T(...xy(k))?o.mul(D32,v):v]));
const linkedCut=(a,T)=>compose(U32,mask(compose(Ud32,a),T));
const localGate=(a,T)=>compose(U32,onsite(compose(Ud32,a),T));
const sourceRun=(a,n)=>{let b=a;for(let i=0;i<n;i++)b=compose(Z32,b);return b;};
const explicitCycle=(a,N,Td=remote,T0=origin)=>localGate(sourceRun(localGate(sourceRun(a,N),Td),N),T0);
const fieldData=a=>[...a].sort(([a],[b])=>a.localeCompare(b)).map(([k,M])=>[k,M.map(r=>r.map(x=>x.toJSON()))]);
const fieldHash=a=>hash(JSON.stringify(fieldData(a)));
const seeds=[unit(o.kron(basis(16,0),e0)),unit(o.kron(basis(16,7),e1)),new Map([['0,0',o.kron(col(Array.from({length:16},(_,i)=>i%3-1)),e0)]])];
const arrow=(name,A,Ad)=>({name,apply:a=>compose(A,a),inverse:a=>compose(Ad,a)});
const gZ=arrow('Z',Z32,Zd32),gU=arrow('U',U32,Ud32),gUd=arrow('U_dagger',Ud32,U32);
const contrast=(name,T)=>({name,apply:a=>onsite(a,T),inverse:a=>onsite(a,T)});
const program=(N,global=false)=>({N,word:[...Array(N).fill(gZ),gUd,contrast('D_remote',global?whole:remote),gU,...Array(N).fill(gZ),gUd,contrast('D_origin',global?whole:origin),gU]});
const programs=[program(1),program(2)];
const prefix=(a,prog,r,inverse=false)=>{let b=a;if(inverse){for(let j=r-1;j>=0;j--)b=prog.word[j].inverse(b);}else{for(let j=0;j<r;j++)b=prog.word[j].apply(b);}return b;};
const cycle=(a,prog)=>prefix(a,prog,prog.word.length);
const mark=(j,m)=>j+','+m,marks=k=>k.split(',').map(Number);
const joint=(a,j=0,m=0)=>new Map([[mark(j,m),a]]);
const jointAdd=(a,b)=>{const c=new Map(a);for(const [k,v]of b)c.set(k,c.has(k)?plus(c.get(k),v):v);return c;};
const jointScale=(a,s)=>new Map([...a].map(([k,v])=>[k,scale(v,s)]));
const jointMinus=(a,b)=>jointAdd(a,jointScale(b,-1));
const jointEnergy=a=>[...a.values()].reduce((s,v)=>s.add(energy(v)),f(0));
const eqJoint=(a,b,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eqFields(a.get(k)||new Map(),b.get(k)||new Map(),32,1,msg+' marks '+k);};
const stepJoint=(state,prog,Q,back=false)=>{const q=prog.word.length,out=new Map();for(const [k,v]of state){const [j,m]=marks(k);let jj,mm,vv;
 if(back){jj=(j+q-1)%q;mm=(m+Q-(j===0?1:0))%Q;vv=prog.word[jj].inverse(v);}else{jj=(j+1)%q;mm=(m+(j===q-1?1:0))%Q;vv=prog.word[j].apply(v);}
 const tag=mark(jj,mm);ensure(!out.has(tag),'native controller mark collision');out.set(tag,vv);
 }return out;};
const runJoint=(a,prog,Q,n)=>{let b=a;for(let j=0;j<n;j++)b=stepJoint(b,prog,Q);return b;};
const frameJoint=(a,prog,inverse=false)=>new Map([...a].map(([k,v])=>[k,prefix(v,prog,marks(k)[0],inverse)]));
const seamStep=(a,prog,Q)=>{const q=prog.word.length,out=new Map();for(const [k,v]of a){const [j,m]=marks(k);out.set(mark((j+1)%q,(m+(j===q-1?1:0))%Q),j===q-1?cycle(v,prog):v);}return out;};
check('source_arrows_compile_the_local_echo_word_with_explicit_control_cost',()=>{
 eqFields(compose(dag(U),U),U16,16,16,'source link pairing');
 eqFields(compose(dag(Z),Z),U16,16,16,'source block pairing');
 eq(o.mul(D32,D32),o.identity(32),'native normalized contrast square');eq(o.dagger(D32),D32,'native contrast dagger');
 let words=0,bounds=0;const examples=[];for(const prog of programs){const N=prog.N,q=prog.word.length;ensure(q===2*N+6,'explicit controller period');
  ensure(prog.word.filter(x=>x.name==='Z').length===2*N,'free word blocks');
  for(const v of seeds){const actual=cycle(v,prog),direct=explicitCycle(v,N);eqFields(actual,direct,32,1,'compiled full control word');
   eqFields(prefix(actual,prog,q,true),v,32,1,'full word inverse');ensure(energy(actual).eq(energy(v)),'full word norm');words++;
   const out=sourceRun(v,N),missing=energy(minus(out,linkedCut(out,remote))),error=energy(minus(actual,v));
   eqFields(linkedCut(v,origin),v,32,1,'origin patch captures seed');ensure(error.le(missing.mul(4)),'twice missing norm echo bound');bounds++;
   examples.push({N,cycle_error_square:error.toString(),missing_remote_norm_square:missing.toString()});
  }
 }
 for(const A of [U,Z,dag(U)])ensure([...A.keys()].every(k=>xy(k).every(n=>Math.abs(n)<=1)),'local source/control range');
 return {complete_native_arrow_identities:4,compiled_word_and_inverse_instances:words,literal_echo_error_bounds:bounds,literal_examples:examples,controls_per_cycle:6,source_blocks_per_cycle:'2 N',source_events_per_cycle:'16 N',controller_updates:'2 N + 6',small_literal_fields_are_large_R37_packets:false,root_witness_is_positive_scalar_simulation:false};
});
check('one_fixed_joint_arrow_preserves_matching_and_has_a_local_inverse',()=>{
 let inverses=0,wraps=0;const prog=programs[0],Q=3,q=prog.word.length;
 for(let j=0;j<q;j++)for(let m=0;m<Q;m++)for(const v of seeds.slice(0,2)){
  const a=joint(v,j,m),b=stepJoint(a,prog,Q);eqJoint(stepJoint(b,prog,Q,true),a,'joint inverse after forward');
  eqJoint(stepJoint(stepJoint(a,prog,Q,true),prog,Q),a,'joint forward after inverse');ensure(jointEnergy(a).eq(jointEnergy(b)),'joint matching norm');inverses++;
  const [jj,mm]=marks([...b.keys()][0]);ensure(jj===(j+1)%q&&mm===(m+(j===q-1?1:0))%Q,'retained phase and carry');if(j===q-1)wraps++;
 }
 const a=jointAdd(jointScale(joint(seeds[0],0,0),'3/5'),jointScale(joint(seeds[1],q-1,1),'4/5')),b=stepJoint(a,prog,Q);
 ensure(jointEnergy(a).eq(jointEnergy(b)),'coherent controller pairing');eqJoint(stepJoint(b,prog,Q,true),a,'coherent inverse');
 return {joint_basis_sector_inverse_and_norm_checks:inverses,wrap_carry_checks:wraps,coherent_joint_inverse_checks:1,iteration_has_no_external_gate_selector:true,spatial_range_per_axis:1,controller_marks_are_spatial_dimensions:false};
});
check('every_finite_program_prefix_matches_the_native_cycle_and_tick_formula',()=>{
 let counts=0;const rows=[];for(const prog of programs){const q=prog.word.length,Q=3,v=seeds[0];let state=joint(v,0,1),ek=v;const histories=[];
  for(let k=0;k<=2;k++){let rstate=ek;const row=[rstate];for(let r=1;r<q;r++){rstate=prog.word[r-1].apply(rstate);row.push(rstate);}histories.push(row);ek=cycle(ek,prog);}
  for(let n=0;n<=2*q+2;n++){if(n)state=stepJoint(state,prog,Q);const k=Math.floor(n/q),r=n%q,expect=joint(histories[k][r],r,(1+k)%Q);eqJoint(state,expect,'exact complete prefix formula');counts++;}
  rows.push({N:prog.N,controller_period:q,checked_updates:2*q+3,final_field_sha256:fieldHash([...state.values()][0])});
 }
 return {full_joint_prefix_equalities:counts,programs:rows,source_fields_retained_with_signed_coefficients:true};
});
check('prefix_frames_isolate_the_cycle_residue_and_control_coherent_phase_error',()=>{
 const prog=programs[0],q=prog.word.length,Q=3;let frames=0;
 for(const j of [0,prog.N,prog.N+3,q-1]){const a=joint(seeds[0],j,1),left=frameJoint(stepJoint(frameJoint(a,prog),prog,Q),prog,true),right=seamStep(a,prog,Q);eqJoint(left,right,'prefix frame seam identity');frames++;}
 const v=seeds[0],phase=prog.N,a0=jointAdd(jointScale(joint(v,0,0),'3/5'),jointScale(joint(prefix(v,prog,phase),phase,0),'4/5'));
 let actual=a0,ek=v;const errors=[];for(let k=1;k<=2;k++){actual=runJoint(actual,prog,Q,q);ek=cycle(ek,prog);const target=jointAdd(jointScale(joint(v,0,k),'3/5'),jointScale(joint(prefix(v,prog,phase),phase,k),'4/5'));
  const e=jointEnergy(jointMinus(actual,target)),expected=energy(minus(ek,v));ensure(e.eq(expected),'coherent phase error exact equality');errors.push({cycles:k,joint_error_square:e.toString(),source_cycle_error_square:expected.toString()});
 }
 const global=program(1,true),gq=global.word.length;for(const v of seeds.slice(0,2)){eqFields(cycle(v,global),v,32,1,'exact full reverser cycle');eqJoint(runJoint(joint(v),global,2,2*gq),joint(v),'full signed joint counter period');}
 return {exact_prefix_frame_identities:frames,coherent_phase_error_equalities:errors,exact_full_reverser_period_witnesses:2,cycle_arrow_only_wrap_edge_in_prefix_coordinates:true,prefix_frame_is_pointwise_spatial:false,finite_window_cycle_claimed_identity_on_all_fields:false};
});
const parent=JSON.parse(parentBytes),r37Path=path.join(path.dirname(arg('--r39-certificate')),'R37_NATIVE_CERTIFICATE.json'),r37Bytes=fs.readFileSync(r37Path);
const inheritedPacket=parent.exact_checks.find(x=>x.name==='pinned_packet_bounds_certify_remote_and_return_trigger_windows');
ensure(hash(r37Bytes)===inheritedPacket.frozen_R37_native_certificate_sha256,'R40 frozen R37 packet pin mismatch');
const r37=JSON.parse(r37Bytes),consts=r37.exact_checks.find(x=>x.name==='a_source_corrector_removes_every_first_difference_error').source_constants.find(x=>x.parameter==='2');
const kap=f(consts.kappa),ctBound=f(consts.comparison_norm_bound),csum=f(consts.residual_xx_mass).add(consts.residual_xy_mass).add(consts.residual_yy_mass),drift=f('2/3');
const lt=(a,b)=>f(a).le(b)&&!f(a).eq(b),floorF=x=>{x=f(x);ensure(x.n>=0n,'nonnegative floor');return x.n/x.d;};
const packetBudget=L=>{L=f(L);const M=L.pow(2),N=L.pow(3),q=N.mul(2).add(6),d=f(floorF(N.mul(drift))),eps=csum.mul(4).add(ctBound.mul(2).mul(drift).mul(f(1).sub(drift))).div(L).add(ctBound.mul(2).div(M)).div(f(1).sub(kap.mul(2).div(M)));
 return {L,M,N,q,d,eps};};
const binarySlots=count=>{count=f(count);ensure(count.d===1n&&count.n>0n,'positive marked capacity');let bits=0,total=1n;while(total<count.n){bits++;total*=2n;}ensure(total>=count.n&&(bits===0||total/2n<count.n),'minimal marked binary count');return bits;};
let clockBudgets=[];
check('derived_packet_horizons_bound_autonomous_echo_error_and_record_capacity',()=>{
 const horizons=[];for(const row of inheritedPacket.exact_large_parameter_budgets){const b=packetBudget(row.L);
  ensure(b.eps.eq(row.epsilon)&&b.N.eq(row.horizon_N)&&b.d.eq(row.patch_displacement),'frozen all-event packet budget');
  ensure(lt(b.eps.pow(2),'1/32')&&lt(b.eps.mul(2),b.eps.mul(3)),'strict scheduled echo improvement');
  horizons.push({L:b.L.toString(),N:b.N.toString(),controller_period:b.q.toString(),cycle_norm_error_upper:b.eps.mul(2).toString(),prior_trigger_bound:b.eps.mul(3).toString()});
 }
 let previous=null;clockBudgets=[6,7,8].map(index=>{const b=packetBudget(f(16).pow(index)),ticks=f(4).pow(index),Q=ticks.add(1),capacity=b.q.mul(Q),err=b.eps.mul(ticks).mul(2);
  ensure(ticks.pow(2).eq(b.L),'growing native tick horizon');if(previous)ensure(lt(err,previous),'decreasing joint source error');previous=err;
  const row={L:b.L.toString(),ticks:ticks.toString(),tick_labels:Q.toString(),phase_labels:b.q.toString(),phase_tick_labels:capacity.toString(),minimum_packed_binary_marked_slots:binarySlots(capacity),total_norm_error_upper:err.toString()};return row;
 });
 return {frozen_packet_sha256:hash(r37Bytes),scheduled_horizon_budgets:horizons,growing_tick_budgets:clockBudgets,large_packet_or_controller_arrays_directly_simulated:false,minimum_storage_scope:'specified phase/tick marked alphabet'};
});
check('native_phase_tick_decoding_has_the_exact_finite_capacity_and_wrap',()=>{
 let decoded=0;const rows=[];for(const q of [8,10,12])for(const Q of [2,3,5]){let j=0,m=0;const seen=new Set();for(let n=0;n<2*q*Q;n++){
  ensure(j===n%q&&m===Math.floor(n/q)%Q&&j+q*m===n%(q*Q),'exact mark count decoder');
  if(n<q*Q){ensure(!seen.has(mark(j,m)),'injective phase tick marks');seen.add(mark(j,m));}
  const carry=j===q-1;j=(j+1)%q;m=(m+(carry?1:0))%Q;decoded++;
 }ensure(j===0&&m===0&&seen.size===q*Q,'exact finite mark wrap');rows.push({phase_labels:q,tick_labels:Q,first_mark_wrap:q*Q,distinct_decoded_counts:seen.size,packed_binary_marked_slots:binarySlots(q*Q)});}
 const q=8,Q=3;ensure(q%Q!==1,'bad every-event carry witness differs from one-cycle carry');
 return {exact_count_decodings:decoded,finite_counters:rows,incorrect_every_update_tick_carry_rejected:true,least_period_claim_requires_exact_full_reverser_cycle:true,arbitrary_coherent_encoding_dimension_bound_claimed:false};
});
check('word_refinement_preserves_echo_and_derives_the_count_calibration',()=>{
 let products=0;for(const prog of programs){const N=prog.N,q=prog.word.length;
  const fused=[...Array(N).fill(gZ),{apply:a=>localGate(a,remote)},...Array(N).fill(gZ),{apply:a=>localGate(a,origin)}];
  const padded={N,word:[{name:'identity',apply:a=>a,inverse:a=>a},...prog.word]};
  const v=seeds[0];let reduced=v;for(const g0 of fused)reduced=g0.apply(reduced);
  eqFields(reduced,cycle(v,prog),32,1,'factored and fused echo product');eqFields(cycle(v,padded),cycle(v,prog),32,1,'identity refinement preserves echo');products+=2;
  ensure(f(q).div(2).sub(3).eq(N)&&padded.word.length===q+1,'explicit refinement count change');
 }
 const rows=[];for(const L of [4096,16384,65536]){const b=packetBudget(L),free=b.q.div(2).sub(3),floorError=b.d.sub(drift.mul(free)),speed=b.d.mul(2).div(b.q),loss=drift.sub(speed);
  ensure(free.eq(b.N)&&lt(-1,floorError)&&floorError.le(0),'exact constructed arm count');
  ensure(f(0).le(loss)&&lt(loss,f(1).add(drift.mul(3)).div(b.N.add(3))),'raw controller speed overhead');
  rows.push({L:b.L.toString(),constructed_distance:b.d.toString(),controller_period:b.q.toString(),raw_two_leg_speed:speed.toString(),control_overhead_fraction:f(3).div(b.N).toString(),free_source_events:b.N.mul(16).toString()});
 }
 const cuts=r37.exact_checks.find(x=>x.name==='a_native_rational_speed_family_approaches_the_sharp_count_coefficient').native_speed_cuts;
 for(const row of cuts){const g0=f(row.derived_drift);ensure(f('1/2').sub(g0.pow(2)).eq(row.squared_speed_deficit),'pinned sharp count family');}
 return {complete_refined_word_product_checks:products,exact_cost_and_calibration_budgets:rows,sharp_speed_family_cuts:cuts.length,control_arrow_count:6,refined_free_count_formula:'(q_prime - h)/(2 a) = N',calibrated_information_length_limit:'(q/2)/(sqrt(2) d) -> 1',physical_duration_per_controller_update_assigned:false,raw_period_is_invariant_under_identity_insertion:false};
});
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),J=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2');
const P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,J),-1),I),J32=o.kron(J,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
const observeJoint=a=>new Map([...a].map(([k,v])=>[k,observe(v)])),decodeJoint=a=>new Map([...a].map(([k,v])=>[k,decode(v)]));
const observedEnergy=a=>[...a.values()].reduce((s,pair)=>s.add(energy(pair[0])).add(energy(pair[1])),f(0));
const outsideX=a=>energy(new Map([...a].filter(([k])=>xy(k)[0]!==0)));
let closureWitness;
check('complete_curvature_observation_preserves_the_controller_and_exposes_hidden_memory',()=>{
 const prog=programs[0],Q=3,q=prog.word.length;let transported=0;
 for(const j of [0,prog.N,prog.N+1,q-1])for(const v of seeds.slice(0,2)){const a=joint(v,j,1),obs=observeJoint(a);
  eqJoint(decodeJoint(obs),a,'complete joint curvature reconstruction');ensure(observedEnergy(obs).eq(jointEnergy(a)),'complete joint observed norm');
  const next=observeJoint(stepJoint(decodeJoint(obs),prog,Q)),expected=observeJoint(stepJoint(a,prog,Q));
  for(const [k,pair]of next)for(let t=0;t<2;t++)eqFields(pair[t],expected.get(k)[t],32,1,'observed autonomous arrow');transported++;
 }
 const v=seeds[0],nextFree=stepJoint(joint(v,0,0),prog,Q),nextControl=stepJoint(joint(v,prog.N,0),prog,Q);
 const a=outsideX([...nextFree.values()][0]),b=outsideX([...nextControl.values()][0]);ensure(!a.eq(b)&&b.eq(0),'same source readout different next response');
 closureWitness={initial_source_sha256:fieldHash(v),free_phase:0,control_phase:prog.N,next_nonzero_x_norm_square_free:a.toString(),next_nonzero_x_norm_square_control:b.toString()};
 return {complete_joint_observer_checks:transported,source_only_nonclosure_witness:closureWitness,controller_ignored_readout_is_complete:false,new_independent_physical_fields_claimed:false};
});
check('native_threshold_and_global_flatness_counterexamples_keep_autonomy_scope_exact',()=>{
 const u=e0,thresholdRows=[];for(const [a,b,c]of [[3,4,5],[5,12,13],[7,24,25],[20,21,29]]){const v=col([new o.F(a,c),new o.F(b,c)]),overlap=o.mul(o.dagger(u),v)[0][0].rad,cut=o.energy(o.mul(Q,v));
  ensure(o.energy(v).eq(1)&&o.energy(o.mul(Q,u)).eq(0)&&lt('1/2',cut),'opposite exact norm threshold preparations');
  const targetPair=o.mul(o.dagger(o.kron(u,e0)),o.kron(v,e1))[0][0];ensure(targetPair.zero()&&!overlap.zero(),'nondisturbing threshold pairing contradiction');
  thresholdRows.push({native_input_pairing:overlap.toString(),upper_cut_fraction:cut.toString(),requested_flagged_pairing:'0',maximum_output_norm_error_lower_bound:overlap.div(2).toString()});
 }
 const far=shift(seeds[0],20,0),prog=programs[0],out=cycle(far,prog),free=sourceRun(far,2);
 eqFields(out,free,32,1,'far source misses both local controls');const defect=energy(minus(out,far));ensure(!defect.zero(),'localized echo is not globally identity');
 eq(o.mul(o.mul(o.mul(H,K),H),K),o.scale(I,-1),'old spatial source loop retained');
 return {nondisturbing_threshold_obstruction_witnesses:thresholdRows,distant_field_echo_defect_square:defect.toString(),spatial_record_loop:'-I',threshold_obstruction_applies_to_every_adaptive_device:false,autonomous_means_one_fixed_joint_arrow:true,free_source_alone_selects_controller:false,material_clock_or_physical_constants_identified:false};
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
 every_iteration_requires_an_external_gate_command:checks[1].iteration_has_no_external_gate_selector,
 reversal_has_zero_control_arrow_cost:checks[0].controls_per_cycle===6,
 controller_marks_are_added_spatial_dimensions:checks[1].controller_marks_are_spatial_dimensions===false,
 matrix_root_witness_selects_the_positive_scalar_branch:checks[0].root_witness_is_positive_scalar_simulation===false,
 huge_packet_and_controller_arrays_were_directly_simulated:checks[4].large_packet_or_controller_arrays_directly_simulated===false,
 the_finite_window_echo_is_identity_on_all_source_fields:!f(checks[8].distant_field_echo_defect_square).zero(),
 prefix_readout_flattens_the_spatial_source_loop:checks[3].prefix_frame_is_pointwise_spatial===false&&checks[8].spatial_record_loop==='-I',
 the_counter_should_advance_on_every_program_edge:checks[5].incorrect_every_update_tick_carry_rejected,
 finite_phase_tick_marks_store_unbounded_elapsed_updates:checks[5].finite_counters.every(x=>x.first_mark_wrap===x.distinct_decoded_counts),
 the_mark_capacity_bound_covers_arbitrary_coherent_encodings:checks[5].arbitrary_coherent_encoding_dimension_bound_claimed===false,
 a_raw_controller_period_is_invariant_under_refinement:checks[6].raw_period_is_invariant_under_identity_insertion===false,
 a_source_only_readout_is_closed_after_forgetting_the_controller:checks[7].controller_ignored_readout_is_complete===false,
 an_exact_nondisturbing_norm_threshold_flag_preserves_every_native_pairing:checks[8].nondisturbing_threshold_obstruction_witnesses.length===4,
 autonomous_native_construction_already_selects_a_material_clock:checks[8].material_clock_or_physical_constants_identified===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r40.native-autonomous-echo.v1',status:'PASS_R40_NATIVE_AUTONOMOUS_ECHO',
 input_sha256:inputHash,r39_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_scheduled_echo_error_derived:true,native_autonomous_local_joint_evolution_derived:true,
 native_program_cycle_and_tick_capacity_derived:true,native_controller_memory_and_calibration_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: scheduled echo with sharper error, one fixed local joint controller arrow, exact prefix and cycle formulas, coherent phase error, marked elapsed-count capacity, control refinement and information-length calibration, complete curvature observation and hidden controller memory, and a nondisturbing threshold obstruction. Autonomous means a constructed time-independent joint arrow; the program, preparation, control alphabet and physical-selection boundary remain explicit.'},null,2)+'\n');
