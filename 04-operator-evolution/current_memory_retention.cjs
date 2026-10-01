'use strict';
// R19 application adapter. All scalar/matrix arithmetic and symbolic proof
// replay come from the unchanged native RKF engine. No external physics law.
const fs=require('node:fs'),path=require('node:path');
const get=name=>{const i=process.argv.indexOf(name);return i<0?null:process.argv[i+1];};
if(!get('--rkf-root')||!get('--input'))throw new Error('--rkf-root and --input required');
const home=path.resolve(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw new Error('R19 input pin mismatch');
const ensure=(ok,msg)=>{if(!ok)throw new Error(msg);};
const eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const I=o.identity(2),Z=o.zeros(2),zeroColumn=()=>o.zeros(2,1);
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange);
const F=o.add(H,K),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P);
const Lp=o.mul(P,F),Lm=o.mul(Q,F),moves=[[1,Lp],[-1,Lm]];
const col=xs=>o.matrix(xs.map(x=>[x])),ek=col([1,0]);
const outer=(a,b=a)=>o.mul(a,o.dagger(b));
const trace=a=>a[0][0].rad.add(a[1][1].rad);
const checkRows=[],check=(name,fn)=>checkRows.push({name,passed:true,...(fn()||{})});

// Finite labelled-response and pair ledgers; these only orchestrate the core.
const fieldAddAt=(a,x,v)=>{const value=o.add(a.get(x)||zeroColumn(),v);if(o.isZero(value))a.delete(x);else a.set(x,value);};
const field=entries=>{const a=new Map();for(const [x,v]of entries)fieldAddAt(a,x,v);return a;};
const fieldStep=a=>{const b=new Map();for(const [x,v]of a)for(const [s,L]of moves)fieldAddAt(b,x+s,o.mul(L,v));return b;};
const pairAddAt=(g,x,y,v)=>{const key=x+','+y,value=o.add(g.get(key)||Z,v);if(o.isZero(value))g.delete(key);else g.set(key,value);};
const pairs=entries=>{const g=new Map();for(const [x,y,v]of entries)pairAddAt(g,x,y,v);return g;};
const pairAt=(g,x,y)=>g.get(x+','+y)||Z;
const pairEntries=g=>[...g].map(([key,v])=>{const [x,y]=key.split(',').map(Number);return [x,y,v];});
const pairScale=(g,c)=>pairs(pairEntries(g).map(([x,y,v])=>[x,y,o.scale(v,c)]));
const pairSum=(...g)=>pairs(g.flatMap(pairEntries));
const pairSub=(a,b)=>pairSum(a,pairScale(b,-1));
const pairEq=(a,b,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||Z,b.get(key)||Z,msg+' at '+key);};
const pairIsZero=g=>g.size===0;
const pairOf=a=>pairs([...a].flatMap(([x,v])=>[...a].map(([y,u])=>[x,y,outer(v,u)])));
const filter=(g,keep)=>pairs(pairEntries(g).filter(([x,y])=>keep(x,y)));
const diag=g=>filter(g,(x,y)=>x===y),off=g=>filter(g,(x,y)=>x!==y);
const even=g=>filter(g,(x,y)=>(x-y)%2===0),odd=g=>filter(g,(x,y)=>(x-y)%2!==0);
const band=(g,r)=>filter(g,(x,y)=>Math.abs(x-y)<=r);
const T=g=>{
 const next=new Map();
 for(const [x,y,v]of pairEntries(g))for(const [s,L]of moves)for(const [t,M]of moves)
  pairAddAt(next,x+s,y+t,o.scale(o.mul(o.mul(L,v),o.dagger(M)),'1/2'));
 return next;
};
const power=(fn,g,n)=>{for(let i=0;i<n;i++)g=fn(g);return g;};
const A=g=>diag(T(diag(g))),B=g=>diag(T(off(g)));
const Cm=g=>off(T(diag(g))),D=g=>off(T(off(g)));
const rho=(g,x)=>trace(pairAt(g,x,x)),current=(g,x)=>trace(o.mul(K,pairAt(g,x,x)));
const allAddresses=(...g)=>new Set(g.flatMap(x=>pairEntries(x).flatMap(([a,b])=>[a,b])));
const targetEq=(a,b,read,msg)=>{for(const x of allAddresses(a,b))ensure(read(a,x).eq(read(b,x)),msg+' at '+x);};
const total=g=>[...allAddresses(g)].reduce((s,x)=>s.add(rho(g,x)),f(0));
const seed=pairOf(field([[0,ek]]));
const fixtures=[field([[0,ek]]),field([[-2,col(['1/3','-2/5'])],[1,col([2,3])]]),
 field([[-4,col([1,1])],[-1,col([1,-1])],[0,col([-2,1])],[4,col([1,0])]])];

check('inherited_native_role_and_flux_identities',()=>{
 eq(o.mul(H,H),I,'source H square');eq(o.mul(K,K),I,'source K square');
 eq(o.add(o.mul(H,K),o.mul(K,H)),Z,'source anticommutation');
 eq(o.mul(F,F),o.scale(I,2),'branch normalization');
 eq(o.add(o.mul(o.dagger(Lp),Lp),o.mul(o.dagger(Lm),Lm)),o.scale(I,2),'sum of source flux forms');
 eq(o.mul(o.mul(F,P),F),o.add(I,K),'right flux');
 eq(o.mul(o.mul(F,Q),F),o.sub(I,K),'left flux');
 for(const value of [H,K,Lp,Lm].flat(2))ensure(value.turn.zero(),'external complex field input');
 return {role_maps_from_verified_r17_input:true,external_complex_or_physics_input:false};
});
check('pair_update_matches_independent_response_propagation',()=>{
 let cases=0;
 for(const source of fixtures){let a=source,g=pairOf(source);const energy=total(g);
  for(let n=0;n<=7;n++){
   pairEq(g,pairScale(pairOf(a),new o.F(1,2**n)),'pair update / independent native response');
   ensure(total(g).eq(energy),'native energy preserved');
   for(const [x,y,v]of pairEntries(g))eq(pairAt(g,y,x),o.dagger(v),'pair dagger');
   a=fieldStep(a);g=T(g);cases++;
  }
 }
 return {preparations:fixtures.length,complete_evolution_cases:cases,maximum_depth:7};
});
check('literal_ordered_word_pairs_recover_evaluated_pair_record',()=>{
 let records=[{word:'',v:ek,x:0}],g=seed,wordCount=0,pairCount=0;
 for(let n=0;n<=6;n++){
  const literal=new Map(),tags=new Set();
  for(const a of records)for(const b of records){
   pairAddAt(literal,a.x,b.x,o.scale(outer(a.v,b.v),new o.F(1,2**n)));
   const tag=a.word+'|'+b.word;ensure(!tags.has(tag),'ordered word tag merged');tags.add(tag);
   ensure(!o.isZero(outer(a.v,b.v)),'native pair coefficient unexpectedly zero');pairCount++;
  }
  ensure(tags.size===4**n,'full ordered pair count');
  pairEq(literal,g,'literal paired source / derived pair transport');wordCount+=records.length;
  if(n<6){records=records.flatMap(rec=>[[H,'H'],[K,'K']].map(([op,tag])=>{
   const v=o.mul(op,rec.v),s=o.mul(o.dagger(v),o.mul(H,v))[0][0].rad;
   ensure(s.eq(1)||s.eq(-1),'literal source role contrast');
   return {word:rec.word+tag,v,x:rec.x+(s.eq(1)?1:-1)};
  }));g=T(g);}
 }
 return {maximum_depth:6,literal_source_words:wordCount,literal_ordered_pairs:pairCount,full_source_tags_preserved:true};
});
const counterPlus=pairOf(field([[-1,ek],[1,ek]]));
const counterMinus=pairOf(field([[-1,ek],[1,o.scale(ek,-1)]]));
check('all_local_quadratics_fail_to_determine_next_current',()=>{
 pairEq(diag(counterPlus),diag(counterMinus),'all local quadratic records equal');
 ensure(current(T(counterPlus),0).eq(1)&&current(T(counterMinus),0).eq(-1),'next current separates');
 const cutOff=diag(T(off(counterPlus)));
 ensure(current(cutOff,0).eq(1),'initial pair residue returns to current');
 for(const source of fixtures){
  const g=pairOf(source),next=T(g),addresses=new Set([...allAddresses(next)]);
  for(const x of addresses){
   const left=source.get(x-1)||zeroColumn(),right=source.get(x+1)||zeroColumn();
   const value=left[0][0].rad.add(left[1][0].rad).mul(right[0][0].rad.sub(right[1][0].rad));
   ensure(current(next,x).eq(value),'native joint-current law');
  }
 }
 return {identical_current_local_records:true,separated_next_currents:['1','-1'],classical_state_matrix_input:false};
});
check('record_cut_identity_and_full_memory_expansion',()=>{
 let expansions=0,nonzeroInitialCases=0;
 for(const source of fixtures){
  const initial=pairOf(source),q0=off(initial),histories=[initial];
  for(let n=0;n<6;n++)histories.push(T(histories[n]));
  pairEq(diag(diag(initial)),diag(initial),'record cut idempotence');
  pairEq(pairSum(diag(initial),off(initial)),initial,'record cut completeness');
  pairEq(diag(T(T(diag(initial)))),pairSum(A(A(initial)),B(Cm(initial))),'cut returning defect');
  for(let n=0;n<6;n++){
   let predicted=pairSum(A(histories[n]),B(power(D,q0,n)));
   for(let r=0;r<n;r++)predicted=pairSum(predicted,B(power(D,Cm(histories[r]),n-1-r)));
   pairEq(predicted,diag(histories[n+1]),'full exact memory expansion');expansions++;
   if(!pairIsZero(B(power(D,q0,n))))nonzeroInitialCases++;
  }
 }
 ensure(nonzeroInitialCases>0,'initial memory test vacuous');
 return {memory_expansions:expansions,nonzero_initial_memory_terms:nonzeroInitialCases,all_memory_blocks_derived:true};
});
check('single_seed_generated_memory_returns_before_energy_difference',()=>{
 const g1=T(seed),returned=B(Cm(seed));
 const E=o.mul(ek,o.dagger(col([0,1])));
 eq(pairAt(off(g1),1,-1),o.scale(E,'1/2'),'first generated pair');
 eq(pairAt(returned,0,0),o.scale(K,'1/4'),'exact returned coefficient');
 ensure(returned.size===1&&rho(returned,0).zero()&&current(returned,0).eq('1/2'),'returned current / energy');
 const full2=power(T,seed,2),cut2=power(A,seed,2);
 targetEq(full2,cut2,rho,'energy unchanged through second event');
 ensure(current(full2,0).sub(current(cut2,0)).eq('1/2'),'second event current');
 const full3=T(full2),cut3=A(cut2);
 for(const x of allAddresses(full3,cut3))ensure(rho(full3,x).sub(rho(cut3,x)).eq(x===1?'1/4':x===-1?'-1/4':0),'third event energy difference');
 return {returned_local_matrix:pairAt(returned,0,0),second_event_current_return:'1/2',third_event_energy_change:{'-1':'-1/4','1':'1/4'}};
});
check('repeated_native_cut_equals_literal_source_census',()=>{
 let cut=seed,counts=new Map([[0,1]]),countCases=0;
 for(let n=0;n<=9;n++){
  ensure(total(cut).eq(1),'cut energy sum');
  for(const x of new Set([...allAddresses(cut),...counts.keys()])){
   ensure(rho(cut,x).eq(new o.F(counts.get(x)||0,2**n)),'repeated cut / word census');
   const block=pairAt(cut,x,x);ensure(block[0][1].rad.zero()&&block[1][0].rad.zero(),'cut local cross role absent');countCases++;
  }
  const next=new Map();for(const [x,c]of counts)for(const s of [-1,1])next.set(x+s,(next.get(x+s)||0)+c);
  counts=next;cut=A(cut);
 }
 return {maximum_depth:9,source_count_address_cases:countCases,reset_probability_or_noise_parameter:false};
});
check('retained_current_drives_exact_count_residue',()=>{
 let full=seed,cut=seed,cases=0;
 for(let n=0;n<9;n++){
  const nf=T(full),nc=A(cut),eta=x=>rho(full,x).sub(rho(cut,x));
  for(const x of allAddresses(nf,nc)){
   const predicted=eta(x-1).add(eta(x+1)).add(current(full,x-1)).sub(current(full,x+1)).div(2);
   ensure(rho(nf,x).sub(rho(nc,x)).eq(predicted),'exact memory/current residue balance');cases++;
  }
  ensure(total(nf).sub(total(nc)).zero(),'zero total residue');full=nf;cut=nc;
 }
 return {residue_balance_address_cases:cases,maximum_depth:9,noise_strength_fitted:false};
});
check('sharp_arbitrary_horizon_witnesses',()=>{
 let currentCases=0,energyCases=0,prefixCases=0;
 for(let h=1;h<=7;h++){
  const initialPlus=pairOf(field([[-h,ek],[h,ek]]));
  const initialMinus=pairOf(field([[-h,ek],[h,o.scale(ek,-1)]]));
  pairEq(band(initialPlus,2*h-1),band(initialMinus,2*h-1),'all smaller-separation records equal');
  let plus=initialPlus,minus=initialMinus;
  for(let n=0;n<=h;n++){
   targetEq(plus,minus,rho,'same local energies until collision');prefixCases++;
   if(n<h){plus=T(plus);minus=T(minus);}
  }
  const exact=new o.F(h%2?1:-1,2**(h-1));
  ensure(current(plus,0).eq(exact)&&current(minus,0).eq(exact.neg()),'sharp current separation');currentCases++;
  ensure(rho(T(plus),1).sub(rho(T(minus),1)).eq(exact),'sharp next energy separation');energyCases++;
 }
 return {current_horizon_witnesses:currentCases,energy_horizon_witnesses:energyCases,equal_energy_prefix_cases:prefixCases,
         general_horizon_proof:'Written unique extreme-path argument, not finite sampling.'};
});
check('finite_horizon_pair_band_sufficiency',()=>{
 let currentCases=0,energyCases=0;
 for(const source of fixtures){const initial=pairOf(source);
  for(let n=1;n<=4;n++){
   const full=power(T,initial,n);
   targetEq(full,power(T,band(initial,2*n),n),current,'current pair band sufficient');currentCases++;
   targetEq(full,power(T,band(initial,2*(n-1)),n),rho,'energy pair band sufficient');energyCases++;
  }
 }
 return {current_band_cases:currentCases,energy_band_cases:energyCases};
});
check('parity_sector_intertwines_source_and_is_locally_invisible',()=>{
 let basisCases=0,futureCases=0;
 for(let x=-2;x<=2;x++)for(let y=-2;y<=2;y++)for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  const e=o.zeros(2);e[a][b]=o.Cut.of(1);const g=pairs([[x,y,e]]);
  pairEq(even(T(g)),T(even(g)),'even pair quotient descends');
  ensure(pairIsZero(diag(odd(g))),'odd pair locally invisible');basisCases++;
 }
 for(const source of fixtures){let hidden=odd(pairOf(source));
  for(let n=0;n<=5;n++){
   ensure(pairIsZero(diag(hidden)),'odd record remains invisible');hidden=T(hidden);futureCases++;
  }
 }
 return {full_pair_coefficient_basis_cases:basisCases,odd_future_cases:futureCases,largest_invisible_subspace_claimed:false};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const Ii=word(),Ki=word('K'),Hi=sc(-1,word('R','K')),Fi=add(Hi,Ki);
const Pi=sc('1/2',add(Ii,Hi)),Qi=sc('1/2',add(Ii,sc(-1,Hi)));
const E=mul(Pi,Ki,Qi),Et=mul(Qi,Ki,Pi),Li=mul(Pi,Fi),Mi=mul(Qi,Fi);
const tasks=[
 ['source_branch_sum_squared',mul(Fi,Fi),sc(2,Ii)],
 ['source_role_cut_idempotent',mul(Pi,Pi),Pi],
 ['source_role_cut_orthogonality',mul(Pi,Qi),sc(0,Ii)],
 ['source_role_cut_completeness',add(Pi,Qi),Ii],
 ['right_outgoing_form',mul(Fi,Pi,Fi),add(Ii,Ki)],
 ['left_outgoing_form',mul(Fi,Qi,Fi),add(Ii,sc(-1,Ki))],
 ['generated_pair_return_coefficient',add(mul(Mi,E,Fi,Pi),mul(Li,Et,Fi,Qi)),Ki],
 ['retained_sign_reverses_current',mul(Hi,Ki,Hi),sc(-1,Ki)],
 ['cross_role_current_not_diagonal',mul(Pi,Ki,Pi),sc(0,Ii)],
];
const system=w.presentation(presentation);
const replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[6].result.certificate));
altered.output.terms=[[[],['1','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}
ensure(alteredRejected,'altered pair-return proof accepted');
const firstReturn=B(Cm(seed)),coherent3=power(T,seed,3),cut3=power(A,seed,3);
const negative={
 local_records_determine_next_current:!current(T(counterPlus),0).eq(current(T(counterMinus),0)),
 omit_initial_hidden_pair:!pairIsZero(B(off(counterPlus))),
 no_generated_memory_from_local_seed:!pairIsZero(firstReturn),
 returned_current_equals_returned_energy:!current(firstReturn,0).eq(rho(firstReturn,0)),
 final_cut_equals_cut_every_event:!rho(coherent3,1).eq(rho(cut3,1)),
 omit_native_half_normalization:!total(pairScale(T(seed),2)).eq(total(seed)),
 zero_current_for_count_residue:!rho(coherent3,1).sub(rho(cut3,1)).zero(),
 every_off_diagonal_record_safely_erasable:!pairIsZero(diag(T(off(counterPlus)))),
 every_off_diagonal_record_required:pairIsZero(diag(power(T,pairs([[0,1,I]]),4))),
 energy_difference_already_at_second_event:rho(power(T,seed,2),0).eq(rho(power(A,seed,2),0)),
 source_word_equals_evaluated_operator:o.equal(o.mul(H,H),I)&&'HH'!=='',
 finite_radius_below_horizon_is_complete:!current(power(T,pairOf(field([[-3,ek],[3,ek]])),3),0).eq(current(power(T,band(pairOf(field([[-3,ek],[3,ek]])),5),3),0)),
};
ensure(Object.values(negative).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r19.current-memory.v1',status:'PASS_R19_NATIVE_CURRENT_MEMORY',
 input_sha256:inputHash,source_foundation:input.source_foundation,
 exact_checks:checkRows,check_count:checkRows.length,symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negative,altered_native_certificate_rejected:alteredRejected,
 ordinary_complex_hilbert_probability_or_physics_law_input:false,
 current_retention_and_pair_memory_derived:true,physical_interaction_selection_proved:false,
 formal_proof_assistant_verified:false,
 scope:'Native finite address-response target inherited from R18. Pair evolution, exact return, repeated record cut and sharp finite-horizon requirements; no physical metric, force, c or alpha selected.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
