#!/usr/bin/env node
'use strict';
// Application witnesses only: all scalar/matrix arithmetic and word proof
// replay use the unchanged canonical native RKF engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const arg=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(arg('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const ensure=(v,msg)=>{if(!v)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const input=JSON.parse(fs.readFileSync(arg('--input'),'utf8')),inputHash=p.digest(input);
ensure(inputHash===arg('--expected-input-sha256'),'R28 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r27-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r27-sha256'),'R28 R27 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R27_NATIVE_REPLICA_SELECTION','R28 parent status');
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const I=o.identity(2),Z=o.zeros(2),H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange);
const R=o.mul(K,H),F=o.commutator(H,K),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),C0=o.add(H,K);
const col=a=>o.matrix(a.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const even=(X,J=H)=>o.scale(o.add(X,o.mul(o.mul(J,X),J)),'1/2');
const odd=(X,J=H)=>o.sub(X,even(X,J));
const addAt=(a,x,M)=>{const v=o.add(a.get(x)||o.zeros(M.length,M[0].length),M);if(o.isZero(v))a.delete(x);else a.set(x,v);};
const plus=(a,b)=>{const c=new Map(a);for(const [x,M]of b)addAt(c,x,M);return c;};
const scale=(a,c)=>new Map([...a].map(([x,M])=>[x,o.scale(M,c)]));
const minus=(a,b)=>plus(a,scale(b,-1));
const shift=(a,d)=>new Map([...a].map(([x,M])=>[x+d,M]));
const left=(X,a)=>new Map([...a].map(([x,M])=>[x,o.mul(X,M)]));
const compose=(a,b)=>{const c=new Map();for(const [i,A]of a)for(const [j,B]of b)addAt(c,i+j,o.mul(A,B));return c;};
const dag=a=>new Map([...a].map(([x,M])=>[-x,o.dagger(M)]));
const eqFields=(a,b,r,c,msg)=>{for(const x of new Set([...a.keys(),...b.keys()]))eq(a.get(x)||o.zeros(r,c),b.get(x)||o.zeros(r,c),msg+' at '+x);};
const unit=X=>new Map([[0,X]]),fieldEnergy=a=>[...a.values()].reduce((s,M)=>s.add(o.energy(M)),f(0));
const power=(a,n)=>{let b=unit(I);for(let i=0;i<n;i++)b=compose(a,b);return b;};
const W=new Map([[1,o.mul(P,C0)],[-1,o.mul(Q,C0)]]),D=new Map([[1,P],[-1,Q]]);
const block=(X,Y)=>compose(compose(unit(X),W),unit(Y));
const A=block(P,P),B=block(P,Q),L=block(Q,P),M=block(Q,Q);
const step=a=>compose(W,a),den=n=>new o.F(1,1n<<BigInt(n));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});
const probes=[unit(e0),unit(e1),new Map([[-2,col([2,-1])],[0,col([-3,4])],[3,col([1,2])]])];

check('native_cut_curvature_is_hidden_but_not_flat',()=>{
 eq(F,o.scale(R,-2),'native signed order curvature');eq(even(F),Z,'odd curvature cut reading');
 eq(o.dagger(F),o.scale(F,-1),'curvature native dagger');eq(o.mul(o.dagger(F),F),o.scale(I,4),'curvature square');
 const hol=o.mul(o.mul(o.mul(H,K),H),K);eq(hol,o.scale(I,-1),'four-letter holonomy');
 eq(o.sub(hol,I),o.mul(o.mul(F,H),K),'order residue versus loop displacement');
 eq(even(o.mul(o.dagger(F),F)),o.scale(I,4),'visible quadratic response');
 eq(o.mul(P,o.mul(F,e0)),o.zeros(2,1),'zero state readout for oriented preparation');
 return {exact_identities:8,curvature:'-2 KH',cut_even_curvature:'0',quadratic_response:'4 I',loop_displacement:'-2 I'};
});

check('propagation_cut_differential_equals_native_curvature',()=>{
 const delta=minus(compose(unit(P),W),compose(W,unit(P)));
 eqFields(delta,scale(compose(D,unit(F)),'1/2'),2,2,'unnormalized cut differential');
 eqFields(compose(dag(delta),delta),unit(I),2,2,'unnormalized exchange square');
 eqFields(A,new Map([[1,P]]),2,2,'visible block');eqFields(B,new Map([[1,o.mul(o.mul(P,K),Q)]]),2,2,'repair block');
 eqFields(L,new Map([[-1,o.mul(o.mul(Q,K),P)]]),2,2,'creation block');eqFields(M,new Map([[-1,o.scale(Q,-1)]]),2,2,'hidden block');
 eqFields(compose(dag(B),B),unit(Q),2,2,'hidden to visible Gram');eqFields(compose(dag(L),L),unit(P),2,2,'visible to hidden Gram');
 eqFields(compose(dag(W),W),unit(o.scale(I,2)),2,2,'raw transport Gram');
 return {Laurent_operator_identities:9,normalization_squared:'1/2',normalized_transition_differential_Gram:'I/2',extra_coupling_parameter:false};
});

check('eliminated_visible_memory_kernel_matches_full_words',()=>{
 let kernels=0,recurrences=0,series=0;
 for(let k=0;k<=16;k++){
  const actual=compose(compose(B,power(M,k)),L),expected=new Map([[-k,o.scale(P,k%2?-1:1)]]);
  eqFields(actual,expected,2,2,'raw memory kernel');kernels++;
  if(k){const previous=new Map([[1-k,o.scale(P,(k-1)%2?-1:1)]]);
   eqFields(plus(actual,shift(previous,-1)),new Map(),2,2,'formal inverse positive coefficient');series++;
  }
 }
 for(const seed of probes){const states=[seed];for(let n=0;n<14;n++)states.push(step(states[n]));
  for(let n=0;n<14;n++){
   let rhs=compose(A,left(P,states[n]));
   for(let k=0;k<n;k++)rhs=plus(rhs,compose(compose(compose(B,power(M,k)),L),left(P,states[n-1-k])));
   rhs=plus(rhs,compose(compose(B,power(M,n)),left(Q,seed)));
   eqFields(left(P,states[n+1]),rhs,2,1,'finite visible memory recurrence');recurrences++;
   const injection=new Map();for(const [x,v]of left(Q,seed))addAt(injection,x+1-n,o.scale(o.mul(K,v),n%2?-1:1));
   eqFields(compose(compose(B,power(M,n)),left(Q,seed)),injection,2,1,'hidden initial response term');recurrences++;
  }
 }
 return {kernel_coefficients:kernels,finite_memory_and_initial_term_checks:recurrences,formal_inverse_cancellations:series,
  maximum_memory_power:16,coefficient_arithmetic:'Exact unnormalized W numerators; U coefficient at word length j is numerator / sqrt(2)^j.'};
});

const hiddenSeed=unit(o.scale(e1,-1)),hiddenOne=step(hiddenSeed);
check('zero_visible_input_can_reappear_from_retained_hidden_energy',()=>{
 eqFields(left(P,hiddenSeed),new Map(),2,1,'zero initial visible reading');
 eqFields(hiddenSeed,unit(o.scale(o.mul(F,e0),'1/2')),2,1,'curvature preparation');
 eqFields(hiddenOne,new Map([[1,o.scale(e0,-1)],[-1,e1]]),2,1,'first native response');
 ensure(fieldEnergy(left(P,hiddenOne)).mul('1/2').eq('1/2'),'visible half-energy');
 ensure(fieldEnergy(left(Q,hiddenOne)).mul('1/2').eq('1/2'),'hidden half-energy');
 eqFields(step(unit(o.zeros(2,1))),new Map(),2,1,'actual zero state remains zero');
 return {initial_visible_energy:'0',initial_retained_energy:'1',next_visible_energy:'1/2',next_hidden_energy:'1/2',external_injection_after_preparation:false};
});

const histories=[];
check('curvature_response_obeys_finite_wave_and_conservation',()=>{
 const W2=power(W,2),W4=power(W,4);
 eqFields(minus(minus(W2,minus(shift(W,1),shift(W,-1))),unit(o.scale(I,2))),new Map(),2,2,'one-event polynomial');
 eqFields(plus(minus(W4,plus(plus(shift(W2,2),scale(W2,2)),shift(W2,-2))),unit(o.scale(I,4))),new Map(),2,2,'two-event wave polynomial');
 let waves=0,norms=0,cones=0,orders=0;
 for(const v of probes){const state=[left(o.scale(F,'1/2'),v)],hk=[left(o.mul(H,K),v)],kh=[left(o.mul(K,H),v)];
  const lo=Math.min(...v.keys()),hi=Math.max(...v.keys()),norm=fieldEnergy(v);
  for(let n=0;n<=16;n++){
   ensure(fieldEnergy(state[n]).mul(den(n)).eq(norm),'curvature norm conservation');norms++;
   ensure([...state[n].keys()].every(x=>x>=lo-n&&x<=hi+n),'native address cone');cones++;
   eqFields(state[n],scale(minus(hk[n],kh[n]),'1/2'),2,1,'same continuation ordered response');orders++;
   if(n<16){state.push(step(state[n]));hk.push(step(hk[n]));kh.push(step(kh[n]));}
  }
  for(let n=2;n<=14;n++){
   const lhs=plus(minus(state[n+2],scale(state[n],4)),scale(state[n-2],4));
   const rhs=plus(minus(shift(state[n],2),scale(state[n],2)),shift(state[n],-2));
   eqFields(lhs,rhs,2,1,'raw signed-response wave');waves++;
  }
  histories.push(state);
 }
 return {Laurent_polynomial_identities:2,signed_wave_checks:waves,norm_checks:norms,support_checks:cones,ordered_preparation_checks:orders,maximum_event:16};
});

check('literal_cut_words_reproduce_curvature_preparation_and_propagation',()=>{
 let prefixes=0,comparisons=0;
 // Start in either native role. Every H/K word remains a signed single role,
 // so its next address is reconstructed from the source role, not a stencil.
 for(const role of [e0,e1]){let words=[{x:0,v:role}],direct=unit(role);
  for(let n=1;n<=10;n++){
   const next=[];for(const old of words)for(const G of [H,K]){const v=o.mul(G,old.v),dir=o.isZero(o.mul(P,v))?-1:1;next.push({x:old.x+dir,v});}
   words=next;const collected=new Map();for(const row of words)addAt(collected,row.x,row.v);
   direct=step(direct);eqFields(direct,collected,2,1,'literal H/K source census');prefixes+=words.length;comparisons++;
  }
 }
 const states=histories[0],rho=(n,x)=>o.energy(states[n].get(x)||o.zeros(2,1)).mul(den(n));
 const lhs=rho(4,0).sub(rho(2,0).mul(2)).add(rho(0,0));
 const rhs=rho(2,-2).sub(rho(2,0).mul(2)).add(rho(2,2)).mul('1/2');
 ensure(lhs.eq('1/8')&&rhs.eq('-1/4')&&!lhs.eq(rhs),'squared response is not a linear wave');
 return {literal_source_prefixes:prefixes,independent_field_comparisons:comparisons,maximum_literal_event:10,
  intensity_wave_rejection:{left:lhs.toString(),right:rhs.toString()}};
});

const J=o.kron(H,I),X=[o.kron(P,o.mul(P,K)),o.kron(o.mul(P,K),o.mul(Q,K)),o.kron(o.mul(Q,K),P)];
const cyc=(xs,fn)=>xs.reduce((acc,_v,i)=>o.add(acc,fn(i,(i+1)%3,(i+2)%3)),o.zeros(xs[0].length));
const projectedCycle=cyc(X,(a,b,c)=>o.commutator(even(X[a],J),even(o.commutator(X[b],X[c]),J)));
check('cyclic_native_curvature_has_exact_hidden_source',()=>{
 let identities=0;
 const cases=[{J:H,X:[H,K,R]},{J,X},...[-2,-1,1,2,3].map(t=>({J,X:X.map((A,i)=>o.add(o.scale(A,t+i),o.scale(o.identity(4),i-1)))}))];
 for(const {J:cut,X:xs}of cases){const full=cyc(xs,(a,b,c)=>o.commutator(xs[a],o.commutator(xs[b],xs[c])));
  const vis=cyc(xs,(a,b,c)=>o.commutator(even(xs[a],cut),even(o.commutator(xs[b],xs[c]),cut)));
  const hid=cyc(xs,(a,b,c)=>o.commutator(odd(xs[a],cut),odd(o.commutator(xs[b],xs[c]),cut)));
  eq(full,o.zeros(cut.length),'finite cyclic identity');eq(o.add(vis,hid),o.zeros(cut.length),'projected cyclic source');identities+=2;
  for(let a=0;a<3;a++)for(let b=0;b<3;b++){
   eq(even(o.commutator(xs[a],xs[b]),cut),o.add(o.commutator(even(xs[a],cut),even(xs[b],cut)),o.commutator(odd(xs[a],cut),odd(xs[b],cut))),'even curvature retains odd pair');identities++;
  }
 }
 eq(projectedCycle,o.kron(P,H),'nonzero visible cyclic witness');
 ensure(!o.isZero(projectedCycle),'hidden source is load-bearing');
 return {native_operator_families:cases.length,cyclic_and_grading_identities:identities,nonzero_visible_source:'P tensor H',imported_coordinate_derivatives:false};
});

check('same_native_curvature_allows_distinct_pairing_preserving_continuations',()=>{
 eq(o.mul(o.dagger(R),R),I,'local turn preserves pairing');
 const local=left(R,hiddenSeed);ensure([...local.keys()].every(x=>x===0),'local turn does not spread');
 ensure([...hiddenOne.keys()].sort().join(',')==='-1,1','native wave target spreads');
 ensure(fieldEnergy(local).eq(fieldEnergy(hiddenSeed))&&fieldEnergy(hiddenOne).mul('1/2').eq(fieldEnergy(hiddenSeed)),'same census conservation, different propagation');
 return {native_continuations:['identity','local KH','R18 address transport'],same_source_curvature:true,unique_dynamics_from_algebra_and_pairing:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const sub=(a,b)=>add(a,sc(-1,b)),comm=(a,b)=>sub(mul(a,b),mul(b,a));
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),ff=comm(hh,kk);
const pp=sc('1/2',add(ii,hh)),qq=sub(ii,pp),ee=X=>sc('1/2',add(X,mul(hh,X,hh)));
const aa=word('A'),bb=word('B'),cc=word('C'),jj=word('J');
const ev=X=>sc('1/2',add(X,mul(jj,X,jj))),od=X=>sub(X,ev(X));
const cycle=(fn)=>add(fn(aa,bb,cc),fn(bb,cc,aa),fn(cc,aa,bb));
const banks=[{name:'canonical_native_source',spec:JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,tasks:[
 ['source_curvature',ff,sc(-2,rr)],['curvature_cut_odd',ee(ff),sc(0,ii)],
 ['curvature_square',mul(ff,ff),sc(-4,ii)],['four_arrow_loop',mul(hh,kk,hh,kk),sc(-1,ii)],
 ['cut_transition_source',comm(pp,add(hh,kk)),sc('1/2',ff)],
 ['hidden_excursion',mul(pp,kk,qq,kk,pp),pp],['native_balanced_Gram',mul(add(hh,kk),add(hh,kk)),sc(2,ii)],
 ]},{name:'native_associative_arrows_and_involution',spec:{tokens:['A','B','C','J'],rules:[{id:'JJ',lhs:['J','J'],rhs:[[[],1]],source:'Native cut self-return only; no arrow commutativity'}]},tasks:[
 ['cyclic_word_cancellation',cycle((a,b,c)=>comm(a,comm(b,c))),sc(0,ii)],
 ['cut_even_curvature_ledger',ev(comm(aa,bb)),add(comm(ev(aa),ev(bb)),comm(od(aa),od(bb)))],
 ['cut_projected_cyclic_source',cycle((a,b,c)=>comm(ev(a),ev(comm(b,c)))),sc(-1,cycle((a,b,c)=>comm(od(a),od(comm(b,c)))))],
 ]}];
const replays=[],presentations=[];let firstSystem;
for(const bank of banks){const system=w.presentation(bank.spec),audit=system.audit();ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','native presentation audit');
 presentations.push({name:bank.name,specification:bank.spec,audit});if(!firstSystem)firstSystem=system;
 for(const [name,left,right]of bank.tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
  ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native symbolic replay '+name);
  replays.push({name,presentation:bank.name,result,replay:'REPLAY_MATCH'});
 }
}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{firstSystem.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const negatives={
 zero_cut_reading_implies_zero_curvature:o.isZero(even(F))&&!o.isZero(F),
 zero_cut_reading_implies_zero_quadratic_response:!o.isZero(even(o.mul(o.dagger(F),F))),
 order_residue_equals_loop_displacement:!o.equal(F,o.scale(I,-2)),
 cut_projection_preserves_products:!o.equal(even(o.mul(F,F)),o.mul(even(F),even(F))),
 zero_visible_input_implies_zero_future_visible_response:fieldEnergy(left(P,hiddenSeed)).eq(0)&&!fieldEnergy(left(P,hiddenOne)).eq(0),
 actual_zero_state_can_supply_response:fieldEnergy(step(unit(o.zeros(2,1)))).eq(0),
 intensity_obeys_signed_response_wave:checks.find(x=>x.name==='literal_cut_words_reproduce_curvature_preparation_and_propagation').intensity_wave_rejection.left!=='-1/4',
 visible_cyclic_curvature_always_closes_without_memory:!o.isZero(projectedCycle),
 nonzero_native_curvature_alone_forces_address_spreading:[...left(R,hiddenSeed).keys()].every(x=>x===0),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r28.native-cut-curvature-propagation.v1',status:'PASS_R28_NATIVE_CUT_CURVATURE_PROPAGATION',
 input_sha256:inputHash,r27_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:presentations,symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_curvature_exchange_derived:true,hidden_memory_kernel_derived:true,signed_curvature_response_wave_derived:true,
 native_chart_free_cyclic_ledger_derived:true,ordinary_complex_or_classical_field_premise:false,
 physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,formal_proof_assistant_verified:false,
 scope:'Native cut curvature, exact transition coupling, eliminated hidden memory, signed response propagation and finite cyclic curvature ledger. The R18 wave target is reused explicitly. No physical vacuum, EM field, c or alpha identification.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
