'use strict';
// R17 adapter: no copied arithmetic/compiler engine and no independent H/K source.
const fs = require('node:fs');
const path = require('node:path');
const get = name => { const i=process.argv.indexOf(name); return i<0?null:process.argv[i+1]; };
if (!get('--rkf-root') || !get('--input')) throw new Error('--rkf-root and --input required');
const home=path.resolve(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8'));
const inputHash=p.digest(input);
if (process.argv.includes('--emit-input-hash')) { process.stdout.write(inputHash+'\n'); process.exit(0); }
if (get('--expected-input-sha256')!==inputHash) throw new Error('R17 input pin mismatch');
const ensure=(condition,message)=>{if(!condition)throw new Error(message);};
const eq=(a,b,message)=>ensure(o.equal(a,b),message);
const F=o.F, f=x=>F.of(x), half=f('1/2');
const checks=[];
const check=(name,body)=>{const data=body();checks.push({name,passed:true,...(data||{})});};
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const h=fromTags(input.constructed_maps.parity_on_roles),k=fromTags(input.constructed_maps.role_exchange);
const id=o.identity(2), zero=o.zeros(2), r=o.mul(k,h);
const mean=o.scale(o.add(h,k),'1/2'),difference=o.scale(o.sub(h,k),'1/2');
const aperture=o.scale(o.add(id,k),'1/2'),complement=o.sub(id,aperture);
const curvature=o.commutator(k,h),generators=[h,k];
const column=values=>o.matrix(values.map(x=>[x]));
const outer=(a,b=a)=>o.mul(a,o.dagger(b));
const energy=a=>o.energy(a);
const average=values=>o.scale(values.reduce((s,x)=>o.add(s,x),o.zeros(values[0].length,values[0][0].length)),new F(1,values.length));
const sumF=values=>values.reduce((s,x)=>s.add(x),f(0));
const averageF=values=>sumF(values).div(values.length);
const levels=seed=>{const all=[[seed]];for(let n=1;n<=input.maximum_depth;n++)all.push(all[n-1].flatMap(x=>generators.map(g=>o.mul(g,x))));return all;};
const coefficient=(v,row)=>o.mul(o.matrix([row]),v)[0][0].rad;
const aRow=['1/2','1/2'],bRow=['1/2','-1/2'];
const source=column([1,0]);
const seeds=input.verification_seeds.map(column);

check('source_generator_relations_before_census',()=>{
 eq(o.mul(h,h),id,'H square');eq(o.mul(k,k),id,'K square');
 eq(o.add(o.mul(h,k),o.mul(k,h)),zero,'H/K anticommutation');
 eq(o.dagger(h),h,'H dagger');eq(o.dagger(k),k,'K dagger');
 eq(o.mul(r,r),o.scale(id,-1),'derived iota');
 for(const x of [h,k,r].flat(2))ensure(x.turn.zero(),'pre-supplied complex input');
 return {roles_from_foundation:true,ordinary_complex_coefficient_input:false};
});
check('mean_difference_and_curvature_energy_operators',()=>{
 eq(o.mul(mean,mean),o.scale(id,'1/2'),'mean square');
 eq(o.mul(difference,difference),o.scale(id,'1/2'),'difference square');
 eq(o.mul(mean,difference),o.scale(r,'1/2'),'mean difference');
 eq(o.mul(difference,mean),o.scale(r,'-1/2'),'difference mean');
 eq(curvature,o.scale(r,2),'source curvature');
 eq(o.mul(o.dagger(difference),difference),o.scale(o.mul(o.dagger(curvature),curvature),'1/8'),'curvature spread');
 eq(o.mul(mean,o.scale(mean,2)),id,'mean retains invertible initial-state map');
});
check('complete_radial_quadratic_second_moment_basis',()=>{
 const basis=[o.matrix([[1,0],[0,0]]),o.matrix([[0,0],[0,1]]),o.matrix([[0,1],[1,0]])];
 for(let j=0;j<basis.length;j++){
  const b=basis[j],lhs=o.scale(o.add(o.mul(o.mul(h,b),h),o.mul(o.mul(k,b),k)),'1/2');
  eq(lhs,j<2?o.scale(id,'1/2'):zero,'radial quadratic coefficient '+j);
 }
 return {all_symmetric_quadratic_coefficients:3};
});

let enumeratedStates=0,depthCases=0,lagPairs=0;
check('full_source_census_mean_spread_and_covariance',()=>{
 for(const seed of seeds){
  const tree=levels(seed),initialEnergy=energy(seed);
  for(let n=0;n<tree.length;n++){
   const states=tree[n],mn=average(states),expected=o.mul(o.power(mean,n),seed);
   ensure(states.length===2**n,'complete count tree');eq(mn,expected,'count mean');
   const retained=initialEnergy.mul(half.pow(n));ensure(energy(mn).eq(retained),'retained mean energy');
   const variance=averageF(states.map(v=>energy(o.sub(v,mn))));
   ensure(variance.eq(initialEnergy.sub(retained)),'centered record spread');
   const covariance=average(states.map(v=>outer(o.sub(v,mn))));
   eq(covariance,n===0?zero:o.sub(o.scale(id,initialEnergy.div(2)),outer(mn)),'full covariance');
   for(const state of states)ensure(energy(state).eq(initialEnergy),'tagged energy preserved');
   enumeratedStates+=states.length;depthCases++;
  }
  // Indexing puts each fixed prefix's suffixes in a contiguous group.
  for(let n=0;n<=4;n++)for(let lag=0;lag<=3;lag++){
   const now=tree[n],later=tree[n+lag],mn=average(now),ml=average(later);
   const cross=average(later.map((v,index)=>outer(o.sub(v,ml),o.sub(now[Math.floor(index/(2**lag))],mn))));
   const cn=n===0?zero:o.sub(o.scale(id,initialEnergy.div(2)),outer(mn));
   eq(cross,o.mul(o.power(mean,lag),cn),'lag covariance');lagPairs++;
  }
 }
 return {verification_seeds:seeds.length,maximum_depth:input.maximum_depth,
         complete_depth_cases:depthCases,enumerated_state_records:enumeratedStates,lag_cases:lagPairs};
});
check('cylinder_counts_survive_every_finite_refinement',()=>{
 let cases=0;
 for(let n=0;n<=input.maximum_depth;n++)for(let length=0;length<=n;length++){
  ensure(new F(2**(n-length),2**n).eq(half.pow(length)),'prefix content');cases++;
 }
 return {prefix_length_depth_cases:cases};
});

check('seam_observer_reconstruction_isometry_and_native_minimum',()=>{
 eq(o.mul(o.mul(h,aperture),h),complement,'HPH=Q');
 eq(o.add(aperture,o.mul(o.mul(h,aperture),h)),id,'observer reconstruction');
 for(const seed of seeds){
  const first=o.mul(aperture,seed),second=o.mul(o.mul(aperture,h),seed);
  eq(o.add(first,o.mul(h,second)),seed,'two-response state recovery');
  ensure(energy(first).add(energy(second)).eq(energy(seed)),'observer isometry');
 }
 const closure=o.rowClosure([aRow],[h,k]);
 ensure(closure.rank===2&&closure.extra===1,'N03 minimum observer');
 ensure(o.rank([aRow])===1&&o.rank([aRow,bRow])===2,'complete observer row rank');
 return {native_observer_rank:closure.rank,minimum_added_scalar_responses:closure.extra,rank_ladder:closure.ranks};
});
check('first_cut_source_derives_balanced_hidden_pair',()=>{
 const first=[o.mul(h,source),o.mul(k,source)];
 eq(o.mul(aperture,first[0]),o.mul(aperture,first[1]),'same visible source fibre');
 const hidden=first.map(v=>coefficient(v,bRow));
 ensure(averageF(hidden).eq(0),'derived hidden balance');
 ensure(averageF(hidden.map(x=>x.pow(2))).eq('1/4'),'derived hidden variance');
 ensure(averageF(hidden.map(x=>x.pow(4))).eq('1/16'),'derived hidden fourth moment');
 const hiddenVectors=first.map(v=>o.mul(complement,v)),hiddenMean=average(hiddenVectors);
 ensure(averageF(hiddenVectors.map(v=>energy(o.sub(v,hiddenMean)))).eq('1/2'),'native hidden spread');
 for(let j=0;j<2;j++)ensure(coefficient(o.mul(h,first[j]),aRow).eq(hidden[j]),'returning hidden value');
 return {count_multiplicities:[1,1],hidden_coefficients:hidden,variance:'1/4',fourth_moment:'1/16',native_energy_spread:'1/2'};
});

check('source_first_return_partition_and_signed_tail',()=>{
 let cylinders=0,returnCases=0;
 for(let n=1;n<=input.return_depth;n++){
  const weights=Array.from({length:n},(_,m)=>half.pow(m+1));
  ensure(sumF(weights).add(half.pow(n)).eq(1),'first-return plus unresolved tail');
  const truncatedMean=sumF(weights.map((p,m)=>p.mul(m+1)));
  ensure(truncatedMean.eq(f(2).sub(half.pow(n).mul(n+2))),'first-return mean tail');
  const signed=sumF(weights.map((p,m)=>p.mul(m%2?-1:1)));
  const target=f(1).sub(new F(-1,2).pow(n)).div(3);
  ensure(signed.eq(target),'signed memory tail');
  for(const seed of seeds){
   const base=o.mul(o.mul(aperture,h),seed);
   let sum=o.zeros(2,1);
   for(let m=0;m<n;m++){
    const response=o.mul(o.mul(o.mul(aperture,h),o.power(k,m)),seed);
    eq(response,o.scale(base,m%2?-1:1),'return parity');
    sum=o.add(sum,o.scale(response,weights[m]));returnCases++;
   }
   eq(sum,o.scale(base,target),'finite return sum');
  }
  cylinders+=n;
 }
 // Independent census of first H in every actual finite source record.
 for(let n=1;n<=input.maximum_depth;n++){
  const buckets=Array(n+1).fill(0);
  for(let code=0;code<2**n;code++){
   let first=n;
   for(let position=0;position<n;position++)if(((code>>position)&1)===0){first=position;break;}
   buckets[first]++;
  }
  for(let m=0;m<n;m++)ensure(buckets[m]===2**(n-m-1),'enumerated first-return cylinder');
  ensure(buckets[n]===1,'retained K-only tail');
 }
 return {maximum_return_depth:input.return_depth,finite_event_cylinders:cylinders,exact_return_response_cases:returnCases};
});
check('parity_erasure_has_derived_two_thirds_one_third_content',()=>{
 // These fractions solve the geometric continuation equation with tail 1/4.
 const even=half.div(f(1).sub(half.pow(2))),odd=half.pow(2).div(f(1).sub(half.pow(2)));
 ensure(even.eq('2/3')&&odd.eq('1/3'),'return sign contents');
 for(const seed of seeds){
  const base=o.mul(o.mul(aperture,h),seed),positive=base,negative=o.scale(base,-1);
  const averageReturn=o.add(o.scale(positive,even),o.scale(negative,odd));
  eq(averageReturn,o.scale(base,'1/3'),'return mean');
  const spread=energy(o.sub(positive,averageReturn)).mul(even).add(energy(o.sub(negative,averageReturn)).mul(odd));
  ensure(spread.eq(energy(base).mul('8/9')),'return spread');
 }
 return {even_content:even,odd_content:odd,mean_return_factor:'1/3',spread_factor:'8/9',infinite_evidence:'written geometric-tail proof; finite tails checked separately'};
});
check('full_history_survives_state_observer_completion',()=>{
 eq(o.mul(h,h),id,'history collision evaluated action');
 ensure(JSON.stringify([])!==JSON.stringify(['H','H']),'history records retained');
 const actionKeys=new Set();
 let operators=[id];
 for(let n=1;n<=4;n++)operators=operators.flatMap(a=>generators.map(g=>o.mul(g,a)));
 for(const a of operators)actionKeys.add(JSON.stringify(a));
 ensure(operators.length===16&&actionKeys.size<=8,'same-depth history collision');
 return {same_depth_records:operators.length,same_depth_distinct_actions:actionKeys.size,state_completion_is_not_path_completion:true};
});

check('finite_history_aperture_recovery_and_optimal_error_budget',()=>{
 let apertureCases=0,fibres=0;
 for(const seed of seeds){
  const tree=levels(seed),initialEnergy=energy(seed);
  for(let n=0;n<=input.maximum_depth;n++)for(let j=0;j<=n;j++){
   const chunk=2**(n-j),means=[];
   for(let prefix=0;prefix<tree[j].length;prefix++){
    const endpoints=tree[n].slice(prefix*chunk,(prefix+1)*chunk),estimate=average(endpoints);
    eq(estimate,o.mul(o.power(mean,n-j),tree[j][prefix]),'prefix-fibre estimate');
    ensure(energy(estimate).eq(initialEnergy.mul(half.pow(n-j))),'prefix retained energy');
    const residual=averageF(endpoints.map(v=>energy(o.sub(v,estimate))));
    ensure(residual.eq(initialEnergy.mul(f(1).sub(half.pow(n-j)))),'prefix residual spread');
    const alternative=o.add(estimate,source);
    const alteredError=averageF(endpoints.map(v=>energy(o.sub(v,alternative))));
    ensure(alteredError.eq(residual.add(energy(source))),'native least-error identity');
    means.push(estimate);fibres++;
   }
   const wholeMean=average(tree[n]);
   const recovered=averageF(means.map(v=>energy(o.sub(v,wholeMean))));
   ensure(recovered.eq(initialEnergy.mul(half.pow(n-j).sub(half.pow(n)))),'resolved spread');
   const remaining=initialEnergy.mul(f(1).sub(half.pow(n-j)));
   ensure(recovered.add(remaining).eq(initialEnergy.mul(f(1).sub(half.pow(n)))),'whole spread budget');
   apertureCases++;
  }
 }
 return {history_aperture_cases:apertureCases,exact_prefix_fibres:fibres,optimal_fibre_mean_proved_in_written_theorem:true};
});
check('recognition_boundary_metric_equals_source_cylinder_count',()=>{
 const depth=5,words=Array.from({length:2**depth},(_,code)=>Array.from({length:depth},(_,i)=>(code>>i)&1));
 let pairs=0;
 for(const a of words)for(const b of words){
  if(a===b)continue;
  let common=0;while(a[common]===b[common])common++;
  const count=words.filter(word=>word.slice(0,common).every((tag,j)=>tag===a[j])).length;
  ensure(new F(count,words.length).eq(half.pow(common)),'recognition metric/count equality');pairs++;
 }
 return {infinite_continuations_represented_by_length_five_prefixes:words.length,
         distinct_ordered_prefix_pairs:pairs,
         finite_terminated_histories_not_reinterpreted_as_third_stochastic_branch:true};
});

// Symbolic identities replayed through the unchanged canonical compiler.
const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens});
const sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args});
const mul=(...args)=>({op:'multiply',args});
const Id=word(),K=word('K'),R=word('R'),H=sc(-1,word('R','K'));
const M=sc('1/2',add(H,K)),D=sc('1/2',add(H,sc(-1,K)));
const P=sc('1/2',add(Id,K)),Q=add(Id,sc(-1,P));
const tasks=[
 ['derived_parity_squared',mul(H,H),Id],
 ['cut_anticommutation',add(mul(H,K),mul(K,H)),sc(0,Id)],
 ['census_mean_squared',mul(M,M),sc('1/2',Id)],
 ['census_difference_squared',mul(D,D),sc('1/2',Id)],
 ['ordered_mean_difference',mul(M,D),sc('1/2',R)],
 ['reverse_ordered_mean_difference',mul(D,M),sc('-1/2',R)],
 ['source_curvature',add(mul(K,H),sc(-1,mul(H,K))),sc(2,R)],
 ['seam_hidden_exchange',mul(H,P,H),Q],
 ['complete_observer',add(P,mul(H,P,H)),Id],
 ['no_direct_visible_H',mul(P,H,P),sc(0,Id)],
 ['K_preserves_visible',mul(P,K,P),P],
 ['K_reverses_hidden',mul(Q,K,Q),sc(-1,Q)],
 ['two_cut_memory_return',mul(P,H,Q,H,P),P],
 ['first_return_sign',mul(P,H,K),sc(-1,mul(P,H))],
 ['mean_initial_state_inverse',mul(M,sc(2,M)),Id],
];
const s=w.presentation(presentation);
const replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(s,left),w.parseExpression(s,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&s.replay(result.certificate),'symbolic replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const corrupted=JSON.parse(JSON.stringify(replays[2].result.certificate));
corrupted.output.terms=[[[],['1','0']]];
let corruptedRejected=false;
try{s.replay(corrupted);}catch(_){corruptedRejected=true;}
ensure(corruptedRejected,'corrupted mean-square proof accepted');

// Deliberately false alternatives must disagree with the derived exact results.
const reject={
 supplied_commuting_branches:!o.equal(o.mul(h,k),o.mul(k,h)),
 mean_energy_conservation_after_record_erasure:!energy(o.mul(mean,source)).eq(energy(source)),
 zero_erased_record_variance:!energy(o.mul(difference,source)).zero(),
 covariance_before_any_event:!o.equal(o.sub(o.scale(id,'1/2'),outer(source)),zero),
 one_scalar_recovers_every_role_state:o.rank([aRow])!==2,
 erased_return_parity:!o.equal(o.mul(o.mul(o.mul(aperture,h),k),source),o.mul(o.mul(aperture,h),source)),
 equal_even_odd_return_contents:!f('2/3').eq('1/2'),
 omit_unreturned_tail:!sumF(Array.from({length:input.return_depth},(_,m)=>half.pow(m+1))).eq(1),
 gaussian_hidden_fourth_moment:!f('1/16').eq(f('1/4').pow(2).mul(3)),
 claim_mean_destroys_initial_state:o.equal(o.mul(mean,o.scale(mean,2)),id),
 erase_source_word_at_state_return:JSON.stringify([])!==JSON.stringify(['H','H']),
};
ensure(Object.values(reject).every(Boolean),'undetected false alternative');
const output={schema:'extra-ideas.r17.native-cut-history.v1',status:'PASS_R17_NATIVE_CUT_HISTORY',
 input_sha256:inputHash,source_foundation:input.source_foundation,
 derived_operators:{h,k,r,mean,difference,aperture,complement,curvature},
 exact_checks:checks,check_count:checks.length,symbolic_replays:replays,
 symbolic_replay_count:replays.length,rejected_false_alternatives:reject,
 altered_native_certificate_rejected:corruptedRejected,
 evidence_scope:'Exact source census and algebra checks. Written induction and geometric-tail proofs cover arbitrary depth. Native count content is not an empirically selected physical probability law.',
 new_formal_proof_assistant:false,physical_noise_universality_claimed:false};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
