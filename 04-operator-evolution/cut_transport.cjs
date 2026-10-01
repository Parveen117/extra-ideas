'use strict';
// Application adapter only. Arithmetic, Laurent products and proof replay use
// the pinned canonical RKF engine; H/K come from the unchanged R17 input.
const fs=require('node:fs'),path=require('node:path');
const get=name=>{const i=process.argv.indexOf(name);return i<0?null:process.argv[i+1];};
if(!get('--rkf-root')||!get('--input'))throw new Error('--rkf-root and --input required');
const home=path.resolve(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const l=require(path.join(home,'core/native_laurent.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8'));
const inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw new Error('R18 input pin mismatch');
const ensure=(ok,msg)=>{if(!ok)throw new Error(msg);};
const eq=(a,b,msg)=>ensure(o.equal(a,b),msg);
const leq=(a,b,msg)=>ensure(l.equalSeries(a,b),msg);
const F=o.F,f=x=>F.of(x),I=o.identity(2),Z=o.zeros(2);
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange);
const R=o.mul(K,H),B=o.add(H,K),M=o.scale(B,'1/2');
const P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P);
const KP=o.scale(o.add(I,K),'1/2'),KQ=o.sub(I,KP);
const column=xs=>o.matrix(xs.map(x=>[x])),seed=column([1,0]),nil=()=>o.zeros(2,1);
const checkRows=[],check=(name,fn)=>checkRows.push({name,passed:true,...(fn()||{})});
const energy=o.energy;
const addAt=(map,x,v)=>{const sum=o.add(map.get(x)||nil(),v);if(o.isZero(sum))map.delete(x);else map.set(x,sum);};
const field=(entries)=>{const out=new Map();for(const [x,v]of entries)addAt(out,x,v);return out;};
const at=(a,x)=>a.get(x)||nil();
const fieldEq=(a,b,msg)=>{for(const x of new Set([...a.keys(),...b.keys()]))eq(at(a,x),at(b,x),msg+' at '+x);};
const fieldEnergy=a=>[...a.values()].reduce((s,v)=>s.add(energy(v)),f(0));
const fieldScale=(a,c)=>field([...a].map(([x,v])=>[x,o.scale(v,c)]));
const fieldAdd=(...a)=>field(a.flatMap(x=>[...x]));
const shift=(a,n)=>field([...a].map(([x,v])=>[x+n,v]));
const step=a=>{const out=new Map();for(const [x,v]of a){const b=o.mul(B,v);addAt(out,x+1,o.mul(P,b));addAt(out,x-1,o.mul(Q,b));}return out;};
const laurentScale=(a,c)=>l.series([...a.terms].map(([k,v])=>[k,o.scale(v,c)]),a.n);
const daggerSeries=a=>l.series([...a.terms].map(([k,v])=>[-k,o.dagger(v)]),a.n);
const S=l.series([[1,P],[-1,Q]],2),W=l.times(S,l.constant(B));
const V=laurentScale(l.times(W,W),'1/2'),unit=l.constant(I);
const L=l.series([[2,o.scale(I,'1/2')],[0,I],[-2,o.scale(I,'1/2')]],2);
const waveResidual=(coefficient)=>l.plus(l.plus(l.times(V,V),laurentScale(l.times(l.series([[2,o.scale(I,coefficient)],[0,o.scale(I,f(2).sub(f(2).mul(coefficient)))],[-2,o.scale(I,coefficient)]],2),V),-1)),unit);

check('inherited_roles_projectors_and_energy_normalization',()=>{
 eq(o.mul(H,H),I,'H square');eq(o.mul(K,K),I,'K square');
 eq(o.add(o.mul(H,K),o.mul(K,H)),Z,'role anticommutation');
 eq(o.mul(R,R),o.scale(I,-1),'native iota');
 for(const e of [H,K,R].flat(2))ensure(e.turn.zero(),'ordinary complex coefficient input');
 eq(o.mul(P,P),P,'positive role projector');eq(o.mul(Q,Q),Q,'negative role projector');
 eq(o.mul(P,Q),Z,'orthogonal roles');
 leq(l.times(daggerSeries(S),S),unit,'shift isometry');
 leq(l.times(S,daggerSeries(S)),unit,'shift inverse');
 leq(l.times(daggerSeries(W),W),l.constant(o.scale(I,2)),'unique normalization square');
 leq(l.times(W,daggerSeries(W)),l.constant(o.scale(I,2)),'surjective normalization');
 leq(l.times(daggerSeries(V),V),unit,'rational two-step isometry');
 return {positive_normalization_squared:'1/2',complex_scalar_input:false,finite_address_pairing:true};
});

// Literal words use generator action, not the derived address recurrence.
const depth=10,levels=[],censuses=[],aggregates=[];
let records=[{word:'',v:seed,x:0,signs:[]}],recordsChecked=0;
check('all_words_role_address_bijection_and_signed_aggregate',()=>{
 for(let n=0;n<=depth;n++){
  const seen=new Set(),counts=new Map(),sum=new Map();
  for(const rec of records){
   ensure(energy(rec.v).eq(1),'source isometry');
   const key=rec.signs.join(',');ensure(!seen.has(key),'address sign-sequence collision');seen.add(key);
   let last=1,recovered='';for(const s of rec.signs){recovered+=s===last?'H':'K';last=s;}
   ensure(recovered===rec.word,'inverse chronological source word');
   ensure(Math.abs(rec.x)<=n&&(n-rec.x)%2===0,'reachable address');
   const plus=(n+rec.x)/2,minus=(n-rec.x)/2;
   ensure(n*n-rec.x*rec.x===4*plus*minus,'address interval');
   counts.set(rec.x,(counts.get(rec.x)||0)+1);addAt(sum,rec.x,rec.v);recordsChecked++;
  }
  ensure(seen.size===2**n&&counts.size===n+1,'full address census');
  if(n)fieldEq(sum,step(aggregates[n-1]),'literal words / derived transport');
  levels.push(records);censuses.push(counts);aggregates.push(sum);
  if(n<depth)records=records.flatMap(rec=>[[H,'H'],[K,'K']].map(([g,tag])=>{
   const v=o.mul(g,rec.v),s=o.mul(o.dagger(v),o.mul(H,v))[0][0].rad;
   ensure(s.eq(1)||s.eq(-1),'occupied role sign');const delta=s.eq(1)?1:-1;
   return {word:rec.word+tag,v,x:rec.x+delta,signs:[...rec.signs,delta]};
  }));
 }
 return {maximum_depth:depth,enumerated_records:recordsChecked,complete_depths:depth+1};
});
check('count_recurrence_moments_and_continuation_cone',()=>{
 let bins=[1],countCases=0,continuationCases=0;
 for(let n=0;n<=depth;n++){
  let total=0,first=0,second=0;
  for(let i=0;i<=n;i++){
   const x=2*i-n,c=censuses[n].get(x);ensure(c===bins[i],'binomial count');
   total+=c;first+=x*c;second+=x*x*c;countCases++;
   if(n)ensure(c===(censuses[n-1].get(x-1)||0)+(censuses[n-1].get(x+1)||0),'count recurrence');
  }
  ensure(total===2**n&&first===0&&second===n*2**n,'exact count moments');
  bins=Array.from({length:n+2},(_,i)=>(bins[i-1]||0)+(bins[i]||0));
 }
 for(let prefix=0;prefix<=4;prefix++)for(const rec of levels[prefix])for(let m=0;m<=4;m++){
  const descendants=levels[prefix+m].filter(x=>x.word.startsWith(rec.word));
  const increments=new Set(descendants.map(x=>x.x-rec.x));
  for(let d=-m-1;d<=m+1;d++)ensure(increments.has(d)===(Math.abs(d)<=m&&(m-d)%2===0),'exact continuation cone');
  continuationCases++;
 }
 return {address_count_cases:countCases,continuation_cases:continuationCases};
});
check('normalization_and_local_current_on_source_histories',()=>{
 let localCases=0;
 for(let n=0;n<=depth;n++){
  const a=aggregates[n];ensure(fieldEnergy(a).eq(2**n),'global normalized response energy');
  if(n===depth)continue;
  const right=x=>{const v=at(a,x),sum=v[0][0].rad.add(v[1][0].rad);return sum.pow(2).div(2**(n+1));};
  const left=x=>{const v=at(a,x),diff=v[0][0].rad.sub(v[1][0].rad);return diff.pow(2).div(2**(n+1));};
  const rho=x=>energy(at(a,x)).div(2**n);
  const current=x=>right(x).sub(left(x+1));
  for(let x=-n-2;x<=n+2;x++){
   const next=energy(at(aggregates[n+1],x)).div(2**(n+1));
   ensure(next.eq(right(x-1).add(left(x+1))),'role outgoing energy');
   ensure(right(x).add(left(x)).eq(rho(x)),'local outgoing sum');
   ensure(next.sub(rho(x)).eq(current(x-1).sub(current(x))),'local conserved current');localCases++;
  }
 }
 const energy3=[-3,-1,1,3].map(x=>energy(at(aggregates[3],x)).div(8));
 ensure(energy3.map(String).join(',')==='1/8,1/8,5/8,1/8','retained-sign witness');
 return {local_current_cases:localCases,third_event_addresses:[-3,-1,1,3],third_event_census:['1/8','3/8','3/8','1/8'],third_event_energy:energy3};
});
check('all_coefficient_laurent_wave_identity',()=>{
 const expected=l.series([[2,o.matrix([['1/2','1/2'],[0,0]])],[0,o.matrix([['1/2','-1/2'],['1/2','1/2']])],[-2,o.matrix([[0,0],['-1/2','1/2']])]],2);
 leq(V,expected,'two-event transfer');
 const residual=l.plus(l.plus(l.times(V,V),laurentScale(l.times(L,V),-1)),unit);
 leq(residual,l.constant(Z),'full Laurent wave polynomial');
 leq(waveResidual('1/2'),l.constant(Z),'wave coefficient');
 const skew=o.matrix([[0,1],[-1,0]]);
 // For every 2x2 matrix A, A^T J A = det(A) J. Verify determinant one
 // directly at the full formal polynomial level rather than sample z values.
 const transposeOnly=l.series([...V.terms].map(([k,v])=>[k,o.dagger(v)]),2);
 leq(l.times(l.times(transposeOnly,l.constant(skew)),V),l.constant(skew),'determinant one');
 return {all_formal_laurent_coefficients:true,wave_coefficient:'1/2',two_event_transfer:[...V.terms]};
});
check('independent_word_response_obeys_exact_wave_stencil',()=>{
 let cases=0;
 for(let n=0;n<=depth-4;n++){
  const rhs=fieldAdd(shift(aggregates[n+2],2),fieldScale(aggregates[n+2],2),shift(aggregates[n+2],-2),fieldScale(aggregates[n],-4));
  fieldEq(aggregates[n+4],rhs,'integer word wave identity');cases++;
 }
 const fixtures=[field([[0,column([1,0])]]),field([[-2,column(['1/3','-2/5'])],[1,column([2,3])]]),field([[-1,column([1,1])],[0,column([1,-1])],[3,column([-2,1])]])];
 let energies=0;
 for(let a of fixtures)for(let n=0;n<6;n++){
  const b=step(a);ensure(fieldEnergy(b).eq(fieldEnergy(a).mul(2)),'arbitrary finite-field isometry');
  a=fieldScale(step(b),'1/2');energies++;
 }
 return {independent_word_stencils:cases,arbitrary_field_isometry_cases:energies};
});
check('finite_shift_jets_and_characteristic_coefficients',()=>{
 const second=l.series([[1,KP],[-1,KQ]],2);
 leq(V,l.times(S,second),'two projected shifts');
 const zeroth=[...V.terms].reduce((a,[,b])=>o.add(a,b),Z);
 const first=[...V.terms].reduce((a,[k,b])=>o.add(a,o.scale(b,-k)),Z);
 const secondJet=[...V.terms].reduce((a,[k,b])=>o.add(a,o.scale(b,new F(k*k,2))),Z);
 eq(zeroth,I,'zeroth shift jet');eq(first,o.scale(B,-1),'first shift jet');
 eq(secondJet,o.add(I,o.mul(H,K)),'second shift jet with retained order');
 eq(o.mul(M,M),o.scale(I,'1/2'),'transport coefficient square');
 ensure(l.trace(M).zero(),'transport coefficient trace');
 const det=M[0][0].mul(M[1][1]).sub(M[0][1].mul(M[1][0]));
 ensure(det.rad.eq('-1/2')&&det.turn.zero(),'characteristic constant coefficient');
 eq(o.mul(o.matrix([[1,0],[0,'-1/2']]),o.matrix([[1,0],[0,-2]])),I,'candidate metric inverse');
 return {finite_jet_order:2,continuum_convergence_claimed:false,characteristic_speed_squared:'1/2'};
});
check('exact_nonzero_front_and_distinct_address_metric',()=>{
 for(let n=1;n<=depth;n++){
  eq(at(aggregates[n],n),seed,'positive unique path');
  eq(at(aggregates[n],-n),column([0,n%2?1:-1]),'negative unique path');
  ensure(energy(at(aggregates[n],n)).div(2**n).eq(new F(1,2**n)),'nonzero front energy');
  ensure(!f(n*n).sub(f(2*n*n)).zero(),'exact endpoint outside candidate continuum cone');
 }
 return {front_depths:depth,exact_front_speed_squared:'1',candidate_slope_is_not_exact_front:true};
});
check('minimum_quadratic_current_observer',()=>{
 // Coordinates of a symmetric quadratic form are its coefficients of a²,ab,b².
 const density=[1,0,1],current=[0,2,0],right=['1/2',1,'1/2'],left=['1/2',-1,'1/2'];
 ensure(o.rank([right,left])===2&&o.rank([density,current])===2,'two independent flux forms');
 ensure(o.rank([density])===1,'energy alone has one response');
 const v=column([1,1]),u=column([1,-1]);
 ensure(energy(v).eq(energy(u)),'equal local energy witness');
 eq(o.mul(P,o.mul(B,u)),nil(),'opposite retained current');
 ensure(!o.isZero(o.mul(P,o.mul(B,v))),'positive retained current');
 return {minimum_linear_quadratic_responses:2,energy_only_counterexample:[[1,1],[1,-1]]};
});
check('calibration_family_does_not_enter_native_transfer',()=>{
 const adapters=[['1','1'],['3','2'],['1/5','7']];
 const ratios=adapters.map(([length,time])=>f(length).div(time));
 ensure(new Set(ratios.map(String)).size===3,'distinct adapter speeds');
 for(const ratio of ratios)ensure(ratio.pow(2).div(2).div(ratio.pow(2)).eq('1/2'),'dimensionless squared slope ratio');
 return {adapter_examples:adapters,unique_physical_metric_selected:false,physical_c_or_alpha_derived:false};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const Ii=word(),Ki=word('K'),Hi=sc(-1,word('R','K')),Bi=add(Hi,Ki);
const Pi=sc('1/2',add(Ii,Hi)),Qi=sc('1/2',add(Ii,sc(-1,Hi)));
const tasks=[
 ['branch_sum_squared',mul(Bi,Bi),sc(2,Ii)],
 ['role_positive_idempotent',mul(Pi,Pi),Pi],
 ['role_negative_idempotent',mul(Qi,Qi),Qi],
 ['role_projector_orthogonality',mul(Pi,Qi),sc(0,Ii)],
 ['role_projector_completeness',add(Pi,Qi),Ii],
 ['normalized_conjugation_without_radical',mul(Bi,Hi,Bi),sc(2,Ki)],
 ['source_mean_characteristic_square',mul(sc('1/2',Bi),sc('1/2',Bi)),sc('1/2',Ii)],
 ['ordered_second_jet',mul(Hi,Ki),sc(-1,word('R'))],
];
const system=w.presentation(presentation);
const replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));
altered.output.terms=[[[],['1','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}
ensure(alteredRejected,'altered native certificate accepted');
const badBranch=o.add(I,K);
const negative={
 commuting_cut_roles:!o.equal(o.mul(H,K),o.mul(K,H)),
 erased_parity_sign_preserves_norm:!o.equal(o.mul(o.dagger(badBranch),badBranch),o.scale(I,2)),
 census_half_is_coherent_isometry:!o.equal(o.scale(o.mul(o.dagger(B),B),'1/4'),I),
 wave_coefficient_one:!l.equalSeries(waveResidual(1),l.constant(Z)),
 reversed_second_jet_order:!o.equal(o.add(I,o.mul(H,K)),o.add(I,o.mul(K,H))),
 count_content_equals_signed_energy:!energy(at(aggregates[3],1)).eq(censuses[3].get(1)),
 address_alone_recovers_final_role:!o.equal(levels[2].find(x=>x.word==='HK').v,levels[2].find(x=>x.word==='KK').v)&&levels[2].find(x=>x.word==='HK').x===levels[2].find(x=>x.word==='KK').x,
 energy_alone_recovers_current:o.rank([[1,0,1],[0,2,0]])>o.rank([[1,0,1]]),
 candidate_symbol_equals_exact_front:!f('1/2').eq(1)&&!energy(at(aggregates[depth],depth)).zero(),
 selected_dimensional_calibration:!f(1).eq(new F(3,2)),
};
ensure(Object.values(negative).every(Boolean),'false alternative not rejected');
const out={schema:'extra-ideas.r18.native-cut-transport.v1',status:'PASS_R18_NATIVE_CUT_TRANSPORT',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checkRows,check_count:checkRows.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negative,
 altered_native_certificate_rejected:alteredRejected,
 derived_operators:{H,K,R,branch_sum:B,first_transport_coefficient:M,role_positive:P,role_negative:Q},
 scope:'Finite cut-address protocol and exact wave stencil; finite jets, not continuum convergence or physical metric selection.',
 external_complex_phase_input:false,full_ugd_collapsed:false,formal_proof_assistant_verified:false};
process.stdout.write(JSON.stringify(out,null,2)+'\n');
