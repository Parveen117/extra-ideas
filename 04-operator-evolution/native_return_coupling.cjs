#!/usr/bin/env node
'use strict';
// Application adapter only. Exact cut arithmetic and proof replay stay upstream.
const fs=require('node:fs'),path=require('node:path');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw Error('R24 input pin mismatch');
const ensure=(ok,msg)=>{if(!ok)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const I=o.identity(2),Z=o.zeros(2),fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange);
const R=o.mul(K,H),B=o.add(I,R),F=o.add(H,K),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]),outer=v=>o.mul(v,o.dagger(v));
const norm=v=>o.energy(v),moment=(v,M)=>o.mul(o.mul(o.dagger(v),M),v)[0][0].rad;
const pow2=n=>1n<<BigInt(n),den=n=>new o.F(1,pow2(n));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});
const addAt=(a,x,M)=>{const z=o.add(a.get(x)||o.zeros(M.length,M[0].length),M);if(o.isZero(z))a.delete(x);else a.set(x,z);};
const step=a=>{const b=new Map();for(const [x,M]of a)for(let j=0;j<2;j++)addAt(b,x+1-2*j,o.mul(L[j],M));return b;};
const stopped=[new Map([[0,I]])],returns=[Z],unrestricted=[I];
let live=stopped[0],free=new Map([[0,I]]);
for(let n=1;n<=64;n++){
 const b=step(live);returns.push(b.get(0)||Z);b.delete(0);live=b;stopped.push(live);
 free=step(free);unrestricted.push(free.get(0)||Z);
}
const gram=a=>[...a.values()].reduce((g,M)=>o.add(g,o.mul(o.dagger(M),M)),Z);
const normReturn=n=>o.scale(o.mul(o.dagger(returns[n]),returns[n]),den(n));
const cumulative=n=>{let G=Z;for(let j=1;j<=n;j++)G=o.add(G,normReturn(j));return G;};

check('stopped_native_norm_balance_and_source_preserving_renewal',()=>{
 let balances=0,renewals=0,translations=0,G=Z;
 for(let n=1;n<=64;n++){
  G=o.add(G,normReturn(n));eq(o.add(G,o.scale(gram(stopped[n]),den(n))),I,'stopped native isometry');balances++;
  let t=Z;for(let j=1;j<=n;j++)t=o.add(t,o.mul(unrestricted[n-j],returns[j]));
  eq(t,unrestricted[n],'ordered first-return renewal');renewals++;
 }
 for(const origin of [-5,3]){
  let a=new Map([[origin,I]]);
  for(let n=1;n<=16;n++){a=step(a);eq(a.get(origin)||Z,returns[n],'translated first-return map');a.delete(origin);translations++;}
 }
 return {exact_pairing_balances:balances,ordered_renewal_coefficients:renewals,translated_return_maps:translations,max_event:64};
});

const choose=(n,k)=>{let x=1n;for(let j=1;j<=k;j++)x=x*BigInt(n-k+j)/BigInt(j);return x;};
const catalan=n=>choose(2*n,n)/BigInt(n+1);
const C=[1n];for(let n=1;n<=128;n++){let z=0n;for(let j=0;j<n;j++)z+=C[j]*C[n-1-j];C.push(z);}
const t=[0n,1n];for(let n=2;n<=128;n++){let z=-2n*t[n-1];for(let j=1;j<n;j++)z+=t[j]*t[n-j];t.push(z);}
const a=n=>new o.F(choose(2*n,n),pow2(2*n)),b=n=>a(n).pow(2);
check('signed_excursion_grammar_matches_native_return_propagation',()=>{
 let matrices=0,series=0,vanishing=0;
 for(let m=1;m<=128;m++){
  const expected=m===1?1n:m%2===0?(m/2%2?-1n:1n)*C[m/2-1]:0n;
  ensure(t[m]===expected,'grammar / signed word-count formula');series++;
 }
 for(let n=1;n<=64;n++){
  const expected=n%2===0?o.scale(B,t[n/2]):Z;
  eq(returns[n],expected,'native stopped propagation / grammar');matrices++;
  if(o.isZero(expected))vanishing++;
 }
 return {formal_coefficients:series,independent_native_return_matrices:matrices,zero_return_matrices:vanishing,
  first_nonzero_events:[2,4,8,12,16],event_six_return_zero:o.isZero(returns[6])};
});

check('native_unsigned_counts_and_rational_tail_inequalities',()=>{
 let countCases=0,wordCases=0,inequalities=0;
 for(let k=0;k<=128;k++){
  ensure(C[k]===catalan(k),'recursive count / reflected-prefix binomial difference');
  ensure(C[k]===choose(2*k,k)-(k?choose(2*k,k-1):0n),'bad-prefix subtraction');countCases++;
  ensure(b(k).le(new o.F(1,k+1)),'squared central count bound');inequalities++;
  ensure(a(k+1).eq(a(k).mul(new o.F(2*k+1,2*k+2))),'native count ratio');
  const lhs=new o.F(1,BigInt(k+2)**3n),rhs=new o.F(1,2).mul(new o.F(1,(k+1)**2).sub(new o.F(1,(k+2)**2)));
  ensure(lhs.le(rhs),'elementary squared-tail telescoping bound');inequalities++;
 }
 for(let k=0;k<=7;k++){
  let good=0;for(let mask=0;mask<2**(2*k);mask++){
   let x=0,ok=true;for(let j=0;j<2*k;j++){x+=((mask>>j)&1)?-1:1;if(x<0)ok=false;}
   if(ok&&x===0)good++;wordCases++;
  }ensure(BigInt(good)===C[k],'literal unsigned positive words');
 }
 return {independent_count_formulas:countCases,literal_direction_words:wordCases,exact_tail_inequalities:inequalities,max_count_index:128};
});

check('literal_HK_first_return_histories_retain_signs_and_cancellation',()=>{
 let visited=0,returned=0,cancelledWords=0;
 for(const seed of [e0,e1]){
  let histories=[{x:0,v:seed,word:''}];
  for(let n=1;n<=12;n++){
   const next=[],arrivals=[];
   for(const h of histories)for(const [M,letter]of [[H,'H'],[K,'K']]){
    const v=o.mul(M,h.v),role=o.isZero(o.mul(P,v))?1:0,x=h.x+1-2*role;
    const row={x,v,word:h.word+letter};visited++;
    if(x===0){arrivals.push(row);returned++;}else next.push(row);
   }
   const sum=arrivals.reduce((v,h)=>o.add(v,h.v),o.zeros(2,1));
   eq(sum,o.mul(returns[n],seed),'literal H/K first-return sum');
   ensure(new Set(arrivals.map(h=>h.word)).size===arrivals.length,'source tags stay distinct');
   if(n===6){ensure(arrivals.length===4&&o.isZero(sum),'nonempty six-event history cancellation');cancelledWords+=arrivals.length;}
   histories=next;
  }
 }
 return {literal_HK_prefixes:visited,distinct_first_return_histories:returned,six_event_histories_with_zero_aggregate:cancelledWords,max_event:12};
});

const floorScaled=(x,digits)=>{const s=10n**BigInt(digits),z=x.n*s;let q=z/x.d;if(z<0n&&z%x.d!==0n)q--;return q;};
const decimal=(x,digits,upper=false)=>{
 let q=floorScaled(x,digits);const scale=10n**BigInt(digits);if(upper&&!new o.F(q,scale).eq(x))q++;
 const sign=q<0n?'-':'';if(q<0n)q=-q;const str=q.toString().padStart(digits+1,'0');
 return sign+str.slice(0,-digits)+'.'+str.slice(-digits);
};
const enclose=(lo,hi,digits=10)=>({lower:decimal(lo,digits),upper:decimal(hi,digits,true),
 lower_rational:new o.F(floorScaled(lo,digits),10n**BigInt(digits)).toString(),
 upper_rational:new o.F(floorScaled(hi,digits)+1n,10n**BigInt(digits)).toString()});
let sqrtLo=f(1),sqrtHi=f(2);
for(let j=0;j<96;j++){const mid=sqrtLo.add(sqrtHi).div(2);if(mid.pow(2).le(2))sqrtLo=mid;else sqrtHi=mid;}
const gLo=sqrtLo.sub(1),gHi=sqrtHi.sub(1),fLo=f(1).sub(f(1).div(sqrtLo)),fHi=f(1).sub(f(1).div(sqrtHi));
let coherentData;
check('coherent_completion_branch_and_certified_rational_enclosures',()=>{
 ensure(sqrtLo.pow(2).le(2)&&f(2).le(sqrtHi.pow(2)),'C8 square-root bracket');
 ensure(gLo.pow(2).add(gLo.mul(2)).le(1)&&f(1).le(gHi.pow(2).add(gHi.mul(2))),'positive return-scale quadratic');
 let sum=f('1/2'),bounds=0;
 for(let k=0;k<=128;k++){
  sum=sum.add(new o.F((k%2?1n:-1n)*C[k],pow2(2*k+2)));
  const next=new o.F((k%2?-1n:1n)*catalan(k+1),pow2(2*k+4));
  const lower=sum.le(sum.add(next))?sum:sum.add(next),upper=sum.le(sum.add(next))?sum.add(next):sum;
  ensure(lower.le(fLo)&&fHi.le(upper),'alternating return sum brackets the selected completed root');bounds++;
 }
 coherentData={g:enclose(gLo,gHi,15),g_squared:enclose(gLo.pow(2),gHi.pow(2),15),
  rational_sqrt_two_bracket:[sqrtLo.toString(),sqrtHi.toString()],bisections:96};
 return {alternating_error_enclosures:bounds,...coherentData,no_fitted_or_empirical_input:true};
});

const etaClosed=M=>f(2*M).add('1/2').add(new o.F(1,8n*BigInt(M+1)**2n)).mul(b(M));
let etaData,etaLo,etaHi;
check('retained_arrival_energy_telescopes_with_certified_completion_tail',()=>{
 let partial=f('1/2'),identities=0,balances=0;
 for(let k=0;k<=128;k++){
  partial=partial.add(new o.F(C[k]*C[k],pow2(4*k+3)));
  ensure(partial.eq(etaClosed(k)),'finite squared-count telescoping identity');identities++;
  if(4*(k+1)<=64){eq(cumulative(4*(k+1)),o.scale(I,partial),'propagated / counted arrival energy');balances++;}
 }
 const M=4096;etaLo=etaClosed(M);const tail=new o.F(1,16n*BigInt(M+1)**2n);etaHi=etaLo.add(tail);
 ensure(f('5/8').le(etaLo)&&etaHi.le('11/16'),'strict nonunit return-energy interval');
 etaData={count_index:M,event_horizon:4*(M+1),tail_upper:tail.toString(),eta:enclose(etaLo,etaHi),
  balanced_record_retention:enclose(f(1).sub(etaHi),f(1).sub(etaLo)),
  negative_role_record_retention:enclose(f(1).sub(etaHi.mul(2)),f(1).sub(etaLo.mul(2))),
  exact_partial_sum_sha256:p.digest(etaLo.toString())};
 return {finite_telescope_identities:identities,independent_propagated_energy_sums:balances,...etaData,
  imported_integral_period_or_special_function_identity:false};
});

const preparations=[e0,e1,col([1,1]),col([1,-1]),col([1,2]),col([-2,3]),col(['1/3','2/5'])];
check('normalized_return_overlap_is_source_parity_not_return_scale',()=>{
 eq(o.mul(o.dagger(B),B),o.scale(I,2),'return normalization');
 eq(o.mul(o.mul(o.dagger(B),K),B),o.scale(H,2),'return exchange is input parity');
 let cases=0,retained=0;
 for(const v of preparations){
  const h=moment(v,H).div(norm(v));let energy=f(0),current=f(0);
  for(let n=1;n<=64;n++){
   const q=o.mul(returns[n],v);energy=energy.add(norm(q).mul(den(n)));current=current.add(moment(q,K).mul(den(n)));
   if(!o.isZero(q)){ensure(moment(q,K).div(norm(q)).eq(h),'every normalized returned record');cases++;}
   if(!energy.zero()){ensure(current.div(energy).eq(h),'normalized arrival-tagged record');retained++;}
  }
 }
 return {nonzero_return_preparations:cases,arrival_retained_normalizations:retained,source_preparations:preparations.length,
  canonical_source_retention:['1','-1','0']};
});

const W=o.add(o.kron(P,I),o.kron(Q,K));
const unresolved=M=>o.matrix([0,1].map(a=>[0,1].map(b=>M[a*2][b*2].add(M[a*2+1][b*2+1]))));
const pairBasis=[];for(let a=0;a<2;a++)for(let b=0;b<2;b++){const G=o.zeros(2);G[a][b]=o.ONE;pairBasis.push(G);}
check('full_joint_return_controlled_gate_derives_completed_coupling',()=>{
 let coefficients=0,branches=0,errorBounds=0;
 for(const v of preparations)for(const M of [0,1,2,3]){
  const N=4*(M+1),E=norm(v),h=moment(v,H).div(E),eta=etaClosed(M),kappa=f(1).sub(eta).add(eta.mul(h));
  const entries=[];for(let n=1;n<=N;n++)if(!o.isZero(returns[n]))entries.push({v:o.mul(returns[n],v),scale:den(n).div(E),returned:true});
  for(const q of stopped[N].values())entries.push({v:o.mul(q,v),scale:den(N).div(E),returned:false});
  branches+=entries.length;
  for(const G of pairBasis){
   let actual=Z;
   for(const entry of entries){const joint=o.scale(o.kron(G,outer(entry.v)),entry.scale);
    actual=o.add(actual,unresolved(entry.returned?o.mul(o.mul(W,joint),o.dagger(W)):joint));}
   const expected=o.matrix(G.map((row,a)=>row.map((z,b)=>a===b?z:z.mul(kappa))));
   eq(actual,expected,'explicit joint branch gate / derived cross coefficient');coefficients+=4;
  }
  const lower=f(1).sub(etaHi.mul(f(1).sub(h))),upper=f(1).sub(etaLo.mul(f(1).sub(h)));
  const tail=new o.F(1,16*(M+1)**2).mul(f(1).sub(h));
  ensure(lower.le(upper)&&kappa.sub(lower).le(tail),'completed coupling within certified tail');errorBounds++;
 }
 return {explicit_joint_pair_coefficients:coefficients,retained_live_and_arrival_branches:branches,
  completed_coupling_error_bounds:errorBounds,physical_interface_selection_assumed:false};
});

check('full_word_memory_changes_return_energy_and_completed_interaction',()=>{
 let words=0,norms=0,currents=0,gates=0;
 for(const v of preparations.slice(0,5)){
  let total=f(0),arrivedCurrent=f(0);
  // Here each entry retains its own full direction word; no equal-address aggregation.
  let rows=[{x:0,v,word:''}];
  for(let n=1;n<=12;n++){
   const next=[],arrived=[];
   for(const h of rows)for(let d=0;d<2;d++){
    const row={x:h.x+1-2*d,v:o.mul(L[d],h.v),word:h.word+d};words++;
    if(row.x===0)arrived.push(row);else next.push(row);
   }
   let hit=f(0);for(const row of arrived){hit=hit.add(norm(row.v).mul(den(n)).div(norm(v)));arrivedCurrent=arrivedCurrent.add(moment(row.v,K).mul(den(n)).div(norm(v)));}
   total=total.add(hit);ensure(arrivedCurrent.zero(),'word-record returned role current zero');currents++;
   if(n%2===0){const m=n/2;ensure(hit.eq(new o.F(2n*catalan(m-1),pow2(2*m))),'literal full-word first-return energy');
    ensure(total.eq(f(1).sub(a(m))),'native first-return word count telescopes');norms++;
    let surviving=f(0);for(const row of next)surviving=surviving.add(norm(row.v).mul(den(n)).div(norm(v)));
    ensure(surviving.add(total).eq(1)&&surviving.add(arrivedCurrent).eq(a(m)),'word-memory returned-controlled gate');gates++;
   }
   ensure(new Set(arrived.map(h=>h.word)).size===arrived.length,'word memory labels unique');rows=next;
  }
 }
 ensure(cumulative(6)[0][0].rad.eq('5/8')&&f(1).sub(a(3)).eq('11/16'),'sharp six-event protocol contrast');
 return {literal_recorded_direction_prefixes:words,exact_return_count_sums:norms,zero_current_checks:currents,word_memory_gate_checks:gates,
  first_energy_difference_event:6,arrival_only_energy_at_six:'5/8',full_word_energy_at_six:'11/16'};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk);
const bb=add(ii,rr),bd=add(ii,sc(-1,rr)),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const tasks=[
 ['role_parity_involution',mul(hh,hh),ii],
 ['role_exchange_involution',mul(kk,kk),ii],
 ['derived_return_orientation',mul(kk,hh),rr],
 ['native_return_orientation_square',mul(rr,rr),sc(-1,ii)],
 ['native_transport_normalization',mul(ff,ff),sc(2,ii)],
 ['return_map_pairing_scale',mul(bd,bb),sc(2,ii)],
 ['returned_exchange_is_source_parity',mul(bd,kk,bb),sc(2,hh)],
 ['returned_parity_is_negative_source_exchange',mul(bd,hh,bb),sc(-2,kk)],
 ['native_role_cuts_disjoint',mul(pp,qq),sc(0,ii)],
 ['native_role_cuts_complete',add(pp,qq),ii],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[5].result.certificate));altered.output.terms=[[[],['3','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native return proof accepted');
const negatives={
 every_possible_first_return_has_nonzero_coherent_signal:o.isZero(returns[6])&&catalan(2)===2n,
 all_even_first_return_amplitudes_are_positive:returns[4][0][0].rad.eq(-1),
 first_return_energy_is_total_return_energy:!o.equal(returns[4],unrestricted[4]),
 stopping_destroys_the_unreturned_source_norm:!o.isZero(gram(stopped[8])),
 all_source_energy_eventually_returns_with_arrival_only_memory:etaHi.le('11/16'),
 coherent_scale_squared_equals_arrival_retained_energy:gHi.pow(2).le('1/5')&&f('5/8').le(etaLo),
 normalized_return_overlap_equals_the_coherent_scale:moment(o.mul(B,e0),K).div(norm(o.mul(B,e0))).eq(1)&&gHi.le('1/2'),
 all_native_source_preparations_have_the_same_record_overlap:moment(e0,H).eq(1)&&moment(e1,H).eq(-1),
 conditional_record_normalization_keeps_return_energy:moment(o.mul(B,col([1,1])),K).zero()&&f(1).sub(etaHi).le('1/2')&&f('1/4').le(f(1).sub(etaHi)),
 cross_retention_must_be_nonnegative:f(1).sub(etaLo.mul(2)).le('-1/4'),
 full_word_recording_preserves_the_arrival_only_return_response:!cumulative(6)[0][0].rad.eq(f(1).sub(a(3))),
 native_gates_alone_select_one_physical_coupling:f('1/4').le(f(1).sub(etaHi))&&a(128).le('1/10'),
 omitted_positive_tail_can_be_set_to_zero:new o.F(C[2]*C[2],pow2(11)).n>0n,
 other_quadratic_branch_is_the_completed_return:f('1/4').le(fLo)&&fHi.le('3/8')&&f(1).le(f(1).add(f(1).div(sqrtHi))),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r24.native-return-coupling.v1',status:'PASS_R24_NATIVE_RETURN_COUPLING',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,
 altered_native_certificate_rejected:alteredRejected,native_first_return_coefficients_derived:true,native_return_controlled_coupling_derived:true,
 ordinary_complex_probability_measurement_or_classical_return_premise:false,
 physical_interaction_interface_selected:false,physical_metric_c_alpha_derived:false,formal_proof_assistant_verified:false,
 scope:'Native signed/refined real role target. Eight written all-depth/completion results, finite exact native regressions and rational error enclosures. Fixed coefficients of explicitly constructed transport, return cuts and record interfaces; no physical source/interface or electromagnetic identification.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
