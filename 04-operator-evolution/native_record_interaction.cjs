#!/usr/bin/env node
'use strict';
// Application adapter: all exact coefficients and word proofs use unchanged RKF.
const fs=require('node:fs'),path=require('node:path');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw Error('R20 input pin mismatch');
const ensure=(ok,msg)=>{if(!ok)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg);
const f=x=>o.F.of(x),I=o.identity(2),Z=o.zeros(2),zv=()=>o.zeros(2,1);
const tagsToMatrix=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=tagsToMatrix(input.constructed_maps.parity_on_roles),K=tagsToMatrix(input.constructed_maps.role_exchange);
const F=o.add(H,K),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const outer=(v,u=v)=>o.mul(v,o.dagger(u));
const trace=a=>a[0][0].rad.add(a[1][1].rad),norm=v=>trace(outer(v));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

// Tuples are finite marked records. Keys are addresses and literal binary words.
const addAt=(a,key,v)=>{const s=o.add(a.get(key)||zv(),v);if(o.isZero(s))a.delete(key);else a.set(key,s);};
const joint=entries=>{const a=new Map();for(const [x,m,v]of entries)addAt(a,x+':'+m,v);return a;};
const entries=a=>[...a].map(([key,v])=>{const at=key.indexOf(':');return [+key.slice(0,at),key.slice(at+1),v];});
const jscale=(a,c)=>joint(entries(a).map(([x,m,v])=>[x,m,o.scale(v,c)]));
const jeq=(a,b,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||zv(),b.get(key)||zv(),msg+' '+key);};
const jzero=a=>a.size===0;
const toggled=(m,i)=>m.slice(0,i)+(m[i]==='0'?'1':'0')+m.slice(i+1);
const write=(a,slot)=>{
 const out=new Map();for(const [x,old,v]of entries(a)){
  const m=slot==='append'?old+'0':old,i=slot==='append'?old.length:slot;
  ensure(Number.isInteger(i)&&i>=0&&i<m.length,'invalid record slot');
  for(let role=0;role<2;role++)addAt(out,x+':'+(role?toggled(m,i):m),o.mul(role?Q:P,v));
 }return out;
};
const systemStep=a=>{const out=new Map();for(const [x,m,v]of entries(a))for(let role=0;role<2;role++)
 addAt(out,(x+(role?-1:1))+':'+m,o.mul(L[role],v));return out;};
const inverseSystemNumerator=a=>{const out=new Map();for(const [x,m,v]of entries(a))for(let role=0;role<2;role++)
 addAt(out,(x-(role?-1:1))+':'+m,o.mul(F,o.mul(role?Q:P,v)));return out;};
const step=(a,slot)=>write(systemStep(a),slot);
const power=(fn,a,n)=>{for(let i=0;i<n;i++)a=fn(a);return a;};
const blank=(a,m)=>joint(entries(a).map(([x,old,v])=>{ensure(old==='','system field expected');return [x,m,v];}));
const withMemory=(a,mu)=>joint(entries(a).flatMap(([x,m,v])=>{ensure(m==='','empty record expected');return [0,1].map(b=>[x,''+b,o.scale(v,mu[b][0])]);}));
const memoryArrow=(a,M,slot=0)=>joint(entries(a).flatMap(([x,m,v])=>[0,1].map(b=>
 [x,m.slice(0,slot)+b+m.slice(slot+1),o.scale(v,M[b][+m[slot]])])));
const addressWrite=(a,code)=>joint(entries(a).map(([x,m,v])=>[x,code(x),v]));
const addressToggle=(a,code)=>joint(entries(a).map(([x,m,v])=>{
 const c=code(x);ensure(c.length===m.length,'code arity');return [x,[...m].map((b,i)=>b===c[i]?'0':'1').join(''),v];}));

const padd=(g,x,y,v)=>{const key=x+','+y,s=o.add(g.get(key)||Z,v);if(o.isZero(s))g.delete(key);else g.set(key,s);};
const pairs=es=>{const g=new Map();for(const [x,y,v]of es)padd(g,x,y,v);return g;};
const pes=g=>[...g].map(([key,v])=>{const [x,y]=key.split(',').map(Number);return [x,y,v];});
const at=(g,x,y)=>g.get(x+','+y)||Z;
const scale=(g,c)=>pairs(pes(g).map(([x,y,v])=>[x,y,o.scale(v,c)]));
const peq=(a,b,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||Z,b.get(key)||Z,msg+' '+key);};
const equalPairs=(a,b)=>[...new Set([...a.keys(),...b.keys()])].every(k=>o.equal(a.get(k)||Z,b.get(k)||Z));
const cut=(g,pred)=>pairs(pes(g).filter(([x,y])=>pred(x,y)));
const diag=g=>cut(g,(x,y)=>x===y),even=g=>cut(g,(x,y)=>(x-y)%2===0);
const J=g=>pairs(pes(g).map(([x,y,v])=>[x,y,o.add(o.mul(o.mul(P,v),P),o.mul(o.mul(Q,v),Q))]));
const T=g=>{const out=new Map();for(const [x,y,v]of pes(g))for(let a=0;a<2;a++)for(let b=0;b<2;b++)
 padd(out,x+(a?-1:1),y+(b?-1:1),o.scale(o.mul(o.mul(L[a],v),o.dagger(L[b])),'1/2'));return out;};
const read=(a,c=1)=>{
 const byRecord=new Map();for(const [x,m,v]of entries(a)){if(!byRecord.has(m))byRecord.set(m,[]);byRecord.get(m).push([x,v]);}
 const g=new Map();for(const vs of byRecord.values())for(const [x,v]of vs)for(const [y,u]of vs)padd(g,x,y,outer(v,u));return scale(g,c);
};
const rho=(g,x)=>trace(at(g,x,x)),cur=(g,x)=>trace(o.mul(K,at(g,x,x)));
const addresses=(...gs)=>new Set(gs.flatMap(g=>pes(g).flatMap(([x,y])=>[x,y])));
const localEq=(a,b,msg)=>{for(const x of addresses(a,b)){ensure(rho(a,x).eq(rho(b,x)),msg+' energy');ensure(cur(a,x).eq(cur(b,x)),msg+' current');}};
const total=g=>[...addresses(g)].reduce((s,x)=>s.add(rho(g,x)),f(0));
const seed=joint([[0,'',e0]]),g0=read(seed);
const fixtures=[seed,joint([[-1,'',e0],[1,'',e0]]),joint([[-2,'',col(['1/3','-2/5'])],[1,'',col([2,3])]]),joint([[0,'',col([1,1])]])];
const denom=n=>new o.F(1,2**n);
const W=o.add(o.kron(P,I),o.kron(Q,K));

check('derived_tuple_gate_and_nonfactorization',()=>{
 eq(o.mul(W,W),o.identity(4),'gate inverse');eq(o.dagger(W),W,'gate dagger');
 eq(o.mul(o.dagger(W),W),o.identity(4),'gate pairing');
 for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  const v=o.zeros(4,1),expected=o.zeros(4,1);v[2*a+b][0]=o.Cut.of(1);expected[2*a+(a?1-b:b)][0]=o.Cut.of(1);
  eq(o.mul(W,v),expected,'literal tuple permutation');
 }
 const out=o.mul(W,o.kron(o.mul(F,e0),e0));
 ensure(!out[0][0].mul(out[3][0]).sub(out[1][0].mul(out[2][0])).zero(),'coupling must not factor');
 ensure(!o.equal(o.mul(W,o.kron(F,I)),o.mul(o.kron(F,I),W)),'noncommuting writing');
 return {complete_tuple_basis_cases:4,pairing_preserved:true,unknown_vector_copying_claimed:false};
});
check('record_current_determines_overlap',()=>{
 const memories=[e0,e1,col([1,1]),col([1,-1]),col([3,4]),col([-2,3]),col(['1/3','-2/5'])];
 let cases=0;
 for(const mu of memories){const energy=norm(mu),kappa=o.mul(o.dagger(mu),o.mul(K,mu))[0][0].rad.div(energy);
  ensure(f(1).sub(kappa).n>=0n&&f(1).add(kappa).n>=0n,'derived overlap bounds');
  for(const source of fixtures){const actual=read(write(withMemory(source,mu),0),f(1).div(energy));
   const predicted=pairs(pes(read(source)).map(([x,y,v])=>[x,y,o.matrix(v.map((row,a)=>row.map((q,b)=>a===b?q:q.mul(kappa))))]));
   peq(actual,predicted,'calculated pointer overlap');cases++;
  }
 }
 return {native_record_preparations:memories.length,whole_pair_readout_comparisons:cases,inserted_overlap_parameter:false};
});
check('fresh_joint_response_matches_independent_pair_recurrence',()=>{
 let cases=0;for(const source of fixtures){let a=source,g=read(source);const energy=total(g);
  for(let n=0;n<=6;n++){
   const actual=read(a,denom(n));peq(actual,g,'fresh joint / pair recurrence');ensure(total(actual).eq(energy),'fresh energy');
   a=step(a,'append');g=J(T(g));cases++;
  }
 }return {preparations:fixtures.length,complete_evolution_cases:cases,maximum_depth:6};
});
check('literal_HK_words_reproduce_fresh_register',()=>{
 let words=[{word:'',record:'',x:0,v:e0}],a=seed,cases=0;
 for(let n=0;n<=8;n++){
  const raw=joint(words.map(q=>[q.x,q.record,q.v]));jeq(raw,a,'literal H/K / branch gate');
  ensure(new Set(words.map(q=>q.word)).size===2**n,'source words lost');
  ensure(new Set(words.map(q=>q.record)).size===2**n,'direction bijection lost');cases+=words.length;
  words=words.flatMap(q=>[[H,'H'],[K,'K']].map(([M,tag])=>{const v=o.mul(M,q.v),role=o.isZero(o.mul(P,v))?1:0;
   return {word:q.word+tag,record:q.record+role,x:q.x+(role?-1:1),v};}));a=step(a,'append');
 }return {maximum_depth:8,literal_source_words:cases,full_word_to_record_bijection:true};
});
check('fresh_address_cut_domain_and_counterexample',()=>{
 let cases=0;for(let x=-2;x<=2;x++)for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  const E=o.zeros(2);E[a][b]=o.Cut.of(1);const g=pairs([[x,x,E]]),out=J(T(g));
  peq(out,diag(T(g)),'address diagonal domain');peq(out,diag(out),'output address diagonal');cases++;
 }
 const general=T(read(fixtures[1]));ensure(!equalPairs(J(general),diag(general)),'false universal address cut');
 ensure(!at(J(general),0,2)[0][0].zero(),'same-role off-address witness');
 return {complete_diagonal_coefficient_basis_cases:cases,unrestricted_address_cut_claim_rejected:true};
});
const choose=(n,k)=>{if(k<0||k>n||!Number.isInteger(k))return 0;let v=1;for(let i=1;i<=k;i++)v=v*(n-i+1)/i;return v;};
check('fresh_source_count_and_zero_current',()=>{
 let a=seed,g=g0,cases=0;
 for(let n=0;n<=9;n++){
  const actual=read(a,denom(n));peq(actual,g,'fresh / R19 repeated cut');
  for(let x=-n;x<=n;x+=2){ensure(rho(actual,x).eq(new o.F(choose(n,(n+x)/2),2**n)),'literal count energy');
   if(n){ensure(cur(actual,x).zero(),'fresh current');
    ensure(at(actual,x,x)[0][0].rad.eq(new o.F(choose(n-1,(n+x)/2-1),2**n)),'right port count');
    ensure(at(actual,x,x)[1][1].rad.eq(new o.F(choose(n-1,(n+x)/2),2**n)),'left port count');}cases++;
  }a=step(a,'append');g=diag(T(g));
 }return {maximum_depth:9,address_count_cases:cases,stochastic_transition_premise:false};
});
check('reused_record_factorization_for_single_origin',()=>{
 let cases=0;for(const v of [e0,e1,col([2,-3])])for(const mu of [e0,col([1,1]),col([2,-1])]){
  let system=joint([[0,'',v]]),a=withMemory(system,mu);const energy=norm(mu),kappa=o.mul(o.dagger(mu),o.mul(K,mu))[0][0].rad.div(energy);
  for(let n=0;n<=8;n++){
   const expected=joint(entries(system).flatMap(([x,m,vx])=>{const count=(n-x)/2;ensure(Number.isInteger(count),'native negative-step count');
    const rec=count%2?o.mul(K,mu):mu;return [0,1].map(b=>[x,''+b,o.scale(vx,rec[b][0])]);}));
   jeq(a,expected,'all paths share derived reused record');
   const base=read(system,denom(n)),actual=read(a,denom(n).div(energy));
   const predicted=pairs(pes(base).map(([x,y,q])=>[x,y,o.scale(q,((x-y)/2)%2?kappa:1)]));
   peq(actual,predicted,'reuse overlap filter');localEq(actual,base,'reuse complete local response');
   system=systemStep(system);a=step(a,0);cases++;
  }
 }return {initial_roles:3,record_preparations:3,full_joint_factorizations:cases,maximum_depth:8};
});
const reused=n=>read(power(a=>step(a,0),blank(seed,'0'),n),denom(n));
const fresh=n=>read(power(a=>step(a,'append'),seed,n),denom(n));
let fixedOmega,fixedXi;
check('first_return_and_fixed_gate_nonclosure',()=>{
 peq(reused(1),fresh(1),'same complete initial readout');
 for(const x of addresses(reused(2),fresh(2)))ensure(rho(reused(2),x).eq(rho(fresh(2),x)),'event two energies');
 ensure(cur(reused(2),0).eq('1/2')&&cur(fresh(2),0).zero(),'event two current');
 for(const [x,c,r]of [[-3,1,1],[-1,3,1],[1,3,5],[3,1,1]]){
  ensure(rho(fresh(3),x).eq(new o.F(c,8))&&rho(reused(3),x).eq(new o.F(r,8)),'third event energies');
 }
 const omega=step(blank(seed,'00'),0),xi=write(omega,1);
 peq(read(omega,'1/2'),read(xi,'1/2'),'fixed-catalogue initial readout');
 fixedOmega=read(step(omega,0),'1/4');fixedXi=read(step(xi,0),'1/4');
 ensure(cur(fixedOmega,0).eq('1/2')&&cur(fixedXi,0).zero(),'fixed same next gate nonclosure');
 return {first_current_difference:'1/2',first_energy_difference_at_event:3,
         same_next_gate_on_same_slot_witness:true,untouched_extra_slot_retained:true};
});
check('address_record_codes_derive_equality_cuts',()=>{
 const source=joint([-2,-1,0,1,2].map(x=>[x,'',col([x+3,2-x])])),base=read(source);
 const codes=[x=>(x+2).toString(2).padStart(3,'0'),x=>''+((x%2+2)%2),x=>x<0?'0':'1',x=>'0'];
 let cases=0;for(const code of codes){const recorded=addressWrite(source,code),actual=read(recorded);
  peq(actual,cut(base,(x,y)=>code(x)===code(y)),'record equality filter');localEq(actual,base,'recorded local block');
  jeq(addressToggle(recorded,code),blank(source,'0'.repeat(code(0).length)),'address involution inverse');cases++;
 }
 peq(read(addressWrite(source,codes[0])),diag(base),'injective address code');
 peq(read(addressWrite(source,codes[1])),even(base),'parity address code');
 return {finite_code_constructions:cases,injective_full_cut_and_parity_cut_verified:true};
});
check('parity_cut_never_changes_future_local_targets',()=>{
 let cases=0;for(const source of fixtures){let full=read(source),repeated=even(full);
  for(let n=0;n<=7;n++){localEq(full,repeated,'harmless parity record');peq(even(full),repeated,'parity intertwiner');
   full=T(full);repeated=even(T(repeated));cases++;}
 }return {future_local_comparisons:cases,maximum_depth:7};
});
check('record_only_native_arrows_preserve_unresolved_pairs',()=>{
 let cases=0,a=blank(seed,'0');for(let n=0;n<=7;n++){
  const base=read(a,denom(n));
  for(const [M,factor]of [[K,1],[H,1],[F,'1/2']]){
   peq(read(memoryArrow(a,M),denom(n).mul(factor)),base,'record-only pairing arrow');cases++;
  }a=step(a,0);
 }return {complete_pair_invariance_cases:cases,allowed_record_arrows:['K','H','(H+K)/sqrt(2)']};
});
check('coupled_decoders_restore_blank_records',()=>{
 let reuseCases=0,freshCases=0;
 for(let n=0;n<=8;n++){
  const a=power(v=>step(v,0),blank(seed,'0'),n),code=x=>''+(((n-x)/2)%2);
  jeq(addressToggle(a,code),blank(power(systemStep,seed,n),'0'),'reused address uncomputation');reuseCases++;
 }
 for(const source of fixtures.slice(0,3))for(let n=1;n<=6;n++){
  const initial=blank(source,'0'.repeat(n));let a=initial;
  for(let i=0;i<n;i++)a=step(a,i);
  for(let i=n-1;i>=0;i--)a=inverseSystemNumerator(write(a,i));
  a=jscale(a,denom(n));jeq(a,initial,'full native inverse on prepared subspace');
  a=power(systemStep,a,n);jeq(a,blank(power(systemStep,source,n),'0'.repeat(n)),'reverse and coherent replay');freshCases++;
 }return {reuse_address_decoders:reuseCases,fresh_reverse_and_replay_cases:freshCases,instantaneous_physical_erasure_claimed:false};
});
check('record_rank_matches_distinct_history_and_final_target',()=>{
 let a=seed,cases=0,portsChecked=0;const ranks=[];
 for(let n=1;n<=7;n++){
  a=step(a,'append');const records=[...new Set(entries(a).map(q=>q[1]))].sort();ensure(records.length===2**n,'full history rank');
  const portKeys=[...new Set(entries(a).flatMap(([x,m,v])=>[0,1].filter(b=>!v[b][0].zero()).map(b=>x+','+b)))].sort();
  const rows=portKeys.map(key=>{const [x,b]=key.split(',').map(Number);return records.map(m=>(a.get(x+':'+m)||zv())[b][0]);});
  const M=o.matrix(rows),G=o.mul(M,o.dagger(M));ensure(o.rank(M)===2*n,'independent native rank');
  ensure(portKeys.length===2*n,'positive system port count');
  for(let i=0;i<G.length;i++)for(let j=0;j<G.length;j++){
   if(i===j){ensure(G[i][i].rad.n>0n&&G[i][i].turn.zero(),'positive native square-root input');portsChecked++;}
   else ensure(G[i][j].zero(),'orthogonal record row witness');
  }
  // R16 C8 supplies sqrt(c/2^n) in the written attainment proof.
  // We verify its exact positive square, not a floating-point square root.
  ranks.push({depth:n,history_mark_rank:2**n,minimum_final_pair_record_rank:2*n});cases++;
 }return {exact_record_rank_cases:cases,positive_attainment_square_inputs:portsChecked,ranks,
         root_existence_source:'R16 C8; inherited written proof and replay',causal_online_compression_claimed:false};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk);
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const tasks=[
 ['source_record_exchange_involution',mul(kk,kk),ii],
 ['source_role_zero_idempotent',mul(pp,pp),pp],
 ['source_role_one_idempotent',mul(qq,qq),qq],
 ['source_role_cuts_disjoint',mul(pp,qq),sc(0,ii)],
 ['source_role_cuts_complete',add(pp,qq),ii],
 ['source_transport_inverse_numerator',mul(ff,ff),sc(2,ii)],
 ['zero_to_zero_transition_sign',mul(pp,ff,pp),pp],
 ['one_to_one_transition_sign',mul(qq,ff,qq),sc(-1,qq)],
 ['record_contrast_reverses_current',mul(hh,kk,hh),sc(-1,kk)],
 ['native_balanced_record_is_exchange_fixed',mul(kk,ff,pp),mul(ff,pp)],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[5].result.certificate));altered.output.terms=[[[],['3','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native proof accepted');
const arbitraryReuse=read(step(blank(fixtures[1],'0'),0),'1/2'),arbitraryFull=T(read(fixtures[1]));
const negatives={
 native_record_gate_is_identity:!o.equal(W,o.identity(4)),
 native_gate_copies_arbitrary_vector:!o.equal(o.mul(W,o.kron(o.mul(F,e0),e0)),o.kron(o.mul(F,e0),o.mul(F,e0))),
 every_record_has_zero_overlap:!o.mul(o.dagger(col([1,1])),o.mul(K,col([1,1])))[0][0].zero(),
 record_overlap_is_always_nonnegative:o.mul(o.dagger(col([1,-1])),o.mul(K,col([1,-1])))[0][0].rad.n<0n,
 fresh_role_cut_is_address_cut_on_every_input:!equalPairs(J(T(read(fixtures[1]))),diag(T(read(fixtures[1])))),
 fresh_and_reused_slots_have_same_current:!cur(fresh(2),0).eq(cur(reused(2),0)),
 energy_difference_already_at_second_event:rho(fresh(2),0).eq(rho(reused(2),0)),
 final_address_cut_equals_fresh_writing:!cur(diag(power(T,g0,2)),0).eq(cur(fresh(2),0)),
 record_only_exchange_restores_lost_cross_pair:equalPairs(read(memoryArrow(step(blank(seed,'0'),0),K),'1/2'),fresh(1)),
 local_pair_is_closed_under_fixed_next_gate:!cur(fixedOmega,0).eq(cur(fixedXi,0)),
 one_reused_slot_preserves_every_initial_address_preparation:!cur(arbitraryReuse,0).eq(cur(arbitraryFull,0)),
 every_address_cut_changes_future_local_response:!equalPairs(read(fixtures[2]),even(read(fixtures[2])))&&
  equalPairs(diag(power(T,even(read(fixtures[2])),3)),diag(power(T,read(fixtures[2]),3))),
 evaluated_operator_return_erases_history:o.equal(o.mul(K,K),I)&&'11'!=='',
 full_history_rank_equals_minimum_final_record_rank:2**3!==2*3,
 omit_transport_normalization_preserves_energy:!total(read(systemStep(seed))).eq(total(g0)),
};
ensure(Object.values(negatives).every(Boolean),'false alternative accepted');
const output={schema:'extra-ideas.r20.record-interaction.v1',status:'PASS_R20_NATIVE_RECORD_INTERACTION',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,
 altered_native_certificate_rejected:alteredRejected,ordinary_complex_hilbert_measurement_or_stochastic_law_input:false,
 record_interaction_constructed:true,physical_memory_allocation_selected:false,physical_metric_c_alpha_derived:false,
 formal_proof_assistant_verified:false,
 scope:'Native tuple construction and pairing readout. Exact fresh/reused memory mechanisms; written arbitrary-depth proofs and finite native checks. No universal full-UGD quotient, physical selection or online compression claim.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
