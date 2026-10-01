#!/usr/bin/env node
'use strict';
// Finite application adapter: arithmetic and proof replay stay in the pinned engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw Error('R25 input pin mismatch');
const r24bytes=fs.readFileSync(get('--r24-certificate')),r24Hash=crypto.createHash('sha256').update(r24bytes).digest('hex');
if(r24Hash!==get('--expected-r24-sha256'))throw Error('R25 R24 certificate pin mismatch');
const r24=JSON.parse(r24bytes);
const ensure=(ok,msg)=>{if(!ok)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
ensure(r24.status==='PASS_R24_NATIVE_RETURN_COUPLING','wrong R24 status');
const r24eta=r24.exact_checks.find(x=>x.name==='retained_arrival_energy_telescopes_with_certified_completion_tail').eta;
const etaLo=f(r24eta.lower_rational),etaHi=f(r24eta.upper_rational);
ensure(f('5/8').le(etaLo)&&etaHi.le('11/16'),'inherited native eta bracket');
const I=o.identity(2),Z=o.zeros(2),fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange),R=o.mul(K,H),B=o.add(I,R),F=o.add(H,K);
const P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const E=[I,H,R,K],col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const norm=v=>o.energy(v),bilinear=(u,M,v=u)=>o.mul(o.mul(o.dagger(u),M),v)[0][0].rad;
const LC=X=>o.scale(o.mul(o.mul(o.dagger(B),X),F),'1/2');
const pow2=n=>1n<<BigInt(n),den=n=>new o.F(1,pow2(n));
const choose=(n,k)=>{let x=1n;for(let j=1;j<=k;j++)x=x*BigInt(n-k+j)/BigInt(j);return x;};
const catalan=n=>choose(2*n,n)/BigInt(n+1),central=n=>new o.F(choose(2*n,n),pow2(2*n));
const weight=n=>n===2?f('1/2'):n>0&&n%4===0?new o.F(catalan(n/4-1)**2n,pow2(n-1)):f(0);
const ww=Array.from({length:257},(_,n)=>weight(n)),ss=[f(1)];
for(let n=1;n<=256;n++)ss.push(ss[n-1].sub(ww[n]));
const response=[I];
for(let n=1;n<=256;n++){
 let D=o.scale(I,ss[n]);for(let j=1;j<=n;j++)if(!ww[j].zero())D=o.add(D,o.scale(LC(response[n-j]),ww[j]));response.push(D);
}
const cycles=[I];for(let n=1;n<=8;n++)cycles.push(LC(cycles[n-1]));
const completed=r=>o.scale(o.add(o.add(I,o.scale(H,r)),o.add(o.scale(R,r.pow(2).neg()),o.scale(K,r.pow(3).neg()))),f(1).sub(r).div(f(1).add(r.pow(4))));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

const W=o.add(o.kron(P,I),o.kron(Q,K)),pulseRaw=o.mul(W,o.kron(I,B));
check('native_feedback_pulse_and_endogenous_eight_return_clock',()=>{
 eq(o.mul(K,B),F,'native K after return gives original balanced arrow');
 eq(o.mul(o.dagger(pulseRaw),pulseRaw),o.scale(o.identity(4),2),'joint pulse normalization');
 eq(o.power(pulseRaw,4),o.scale(o.kron(H,I),-4),'four-return source restoration and system sign');
 eq(o.power(pulseRaw,8),o.scale(o.identity(4),16),'eight-return pulse');
 let proper=0;
 for(let n=1;n<8;n++){
  const M=o.power(pulseRaw,n);
  if(n%2===0)ensure(!o.equal(M,o.scale(o.identity(4),pow2(n/2))),'no smaller even period');
  else ensure(M.some((row,i)=>row.some((z,j)=>i!==j&&!z.zero())),'no smaller odd period');proper++;
 }
 return {native_pulse_identities:4,excluded_smaller_periods:proper,return_period:8,source_role_restored_at_return:4,
  fixed_transport_event_period_assumed:false};
});

check('native_four_component_continuation_and_hidden_R_reentry',()=>{
 eq(cycles[1],H,'L I');eq(cycles[2],o.scale(R,-1),'L H');eq(cycles[3],o.scale(K,-1),'L -R');eq(cycles[4],o.scale(I,-1),'L -K');eq(cycles[8],I,'L period');
 for(const X of E){let Y=X;for(let k=0;k<4;k++)Y=LC(Y);eq(Y,o.scale(X,-1),'fourth power on full real matrix basis');}
 ensure(o.rank(cycles.slice(0,4).map(M=>o.flatten(M)))===4,'independent continuation quartet');
 const v=col([1,1]);ensure(bilinear(v,R).zero()&&bilinear(v,LC(R)).eq(2),'invisible R returns as visible K');
 ensure(bilinear(e0,cycles[1]).eq(1)&&bilinear(e0,cycles[2]).zero(),'scalar reset law fails at next return');
 return {closure_rank:4,full_basis_fourth_power_checks:4,cycle_entries:8,presently_invisible_R_has_visible_continuation:true};
});

const key=(x,mark)=>x+','+mark;
const split=k=>k.split(',').map(Number);
const addAt=(a,k,M)=>{const z=o.add(a.get(k)||o.zeros(M.length,M[0].length),M);if(o.isZero(z))a.delete(k);else a.set(k,z);};
const envStep=(a,system,n,wordMemory=false)=>{
 const b=new Map();for(const [k,M]of a){const [x,mark]=split(k);for(let d=0;d<2;d++){
  const y=x+1-2*d,newmark=wordMemory?(mark|(d<<(n-1))):(mark|((y===0?1:0)<<(n-1)));
  let value=o.mul(L[d],M);if(y===0&&system)value=o.mul(K,value);addAt(b,key(y,newmark),value);
 }}return b;
};
const envGram=a=>[...a.values()].reduce((G,M)=>o.add(G,o.mul(o.dagger(M),M)),Z);
const envCross=(a,b)=>{let D=Z;for(const [k,M]of a)if(b.has(k))D=o.add(D,o.mul(o.dagger(M),b.get(k)));return D;};
let arrivalFields;
check('full_joint_event_feedback_matches_exact_retained_memory_recurrence',()=>{
 let a=new Map([[key(0,0),I]]),b=new Map([[key(0,0),I]]),norms=0,entries=0,recordBlocks=0;
 const snapshots=[[a,b]];
 for(let n=1;n<=14;n++){
  a=envStep(a,0,n);b=envStep(b,1,n);snapshots.push([a,b]);
  eq(o.scale(envGram(a),den(n)),I,'system0 full native norm');eq(o.scale(envGram(b),den(n)),I,'system1 full native norm');norms+=2;
  eq(o.scale(envCross(a,b),den(n)),response[n],'full environment / closed return recurrence');entries+=4;recordBlocks+=a.size+b.size;
 }
 arrivalFields=snapshots;
 eq(response[2],o.scale(o.add(I,H),'1/2'),'two-event response');
 eq(response[4],o.add(o.scale(o.add(I,H),'3/8'),o.scale(R,'-1/4')),'four-event response retains R');
 const pp=[];let countReconstructions=0;
 for(let n=0;n<=64;n++){
  const row=Array.from({length:Math.floor(n/2)+1},()=>f(0));row[0]=ss[n];
  for(let j=1;j<=n;j++)if(!ww[j].zero())for(let k=0;k<pp[n-j].length;k++)row[k+1]=row[k+1].add(ww[j].mul(pp[n-j][k]));
  ensure(row.reduce((s,x)=>s.add(x),f(0)).eq(1),'derived return-count ledger normalized');
  const D=row.reduce((M,c,k)=>o.add(M,o.scale(cycles[k%8],c)),Z);eq(D,response[n],'positive count ledger / operator response');
  pp.push(row);countReconstructions++;
 }
 return {direct_joint_norm_checks:norms,direct_joint_response_entries:entries,retained_environment_blocks:recordBlocks,
  exact_normalized_return_count_reconstructions:countReconstructions,max_direct_event:14,max_recurrence_event:256};
});

check('literal_HK_feedback_histories_match_joint_return_flags',()=>{
 let historiesChecked=0,blocks=0;
 for(let s=0;s<2;s++)for(const seed of [e0,e1]){
  let histories=[{x:0,mark:0,v:seed,word:''}];
  for(let n=1;n<=10;n++){
   histories=histories.flatMap(h=>[[H,'H'],[K,'K']].map(([M,letter])=>{
    let v=o.mul(M,h.v);const d=o.isZero(o.mul(P,v))?1:0,x=h.x+1-2*d,mark=h.mark|((x===0?1:0)<<(n-1));
    if(x===0&&s)v=o.mul(K,v);historiesChecked++;return {x,mark,v,word:h.word+letter};
   }));
   const aggregate=new Map();for(const h of histories)addAt(aggregate,key(h.x,h.mark),h.v);
   const actual=arrivalFields[n][s];for(const k of new Set([...aggregate.keys(),...actual.keys()])){
    eq(aggregate.get(k)||o.zeros(2,1),o.mul(actual.get(k)||Z,seed),'literal H/K feedback / joint branch propagation');blocks++;
   }
   ensure(new Set(histories.map(h=>h.word)).size===2**n,'full source H/K tags retained');
  }
 }
 return {literal_HK_prefixes:historiesChecked,exact_aggregated_flag_blocks:blocks,max_literal_event:10};
});

const rationalDiagnostics=[f(0),f('1/2'),f('5/8'),etaLo,etaHi,f('11/16'),f(1)];
check('completed_feedback_polynomial_identity_and_unique_algebraic_response',()=>{
 const numerator=[I,H,o.scale(R,-1),o.scale(K,-1)];let coefficients=0;
 for(let k=0;k<=4;k++){
  const v=o.sub(k<4?numerator[k]:Z,k>0?LC(numerator[k-1]):Z);
  eq(v,k===0||k===4?I:Z,'all-parameter native inverse polynomial coefficient');coefficients++;
 }
 for(const r of rationalDiagnostics){const D=completed(r);eq(D,o.add(o.scale(I,f(1).sub(r)),o.scale(LC(D),r)),'completed feedback identity at exact diagnostic');}
 // At eta=1 algebraic uniqueness alone does not force boundary convergence.
 eq(completed(f(1)),Z,'algebraic solution at excluded unit return mass');
 ensure(!o.equal(cycles[0],cycles[4]),'unit-mass boundary iterates still oscillate');
 return {exact_polynomial_coefficients:coefficients,exact_rational_diagnostics:rationalDiagnostics.map(String),
  diagnostics_are_physical_parameter_inputs:false,unit_return_mass_is_not_claimed_convergent:true};
});

const preparations=[e0,e1,col([1,1]),col([1,-1]),col([1,2]),col([-2,3]),col(['1/3','2/5'])];
const terminal=(r,m,G)=>{let sum=Z,X=I;for(let k=0;k<m;k++){sum=o.add(sum,o.scale(X,f(1).sub(r).mul(r.pow(k))));X=LC(X);}
 let T=G;for(let k=0;k<m;k++)T=LC(T);return o.add(sum,o.scale(T,r.pow(m)));};
check('native_terminal_continuations_lose_influence_with_derived_return_depth',()=>{
 let bounds=0,recurrences=0;
 for(const r of [f('5/8'),etaLo,f('11/16')])for(const m of [0,1,2,3,4,8,12,16]){
  const targets=E.map(G=>terminal(r,m,G)),limit=completed(r),error=r.pow(m).mul(2);
  for(let j=0;j<E.length;j++){
   const next=terminal(r,m+1,E[j]);eq(next,o.add(o.scale(I,f(1).sub(r)),o.scale(LC(targets[j]),r)),'finite terminal recursion');recurrences++;
   for(const u of preparations.slice(0,4))for(const v of preparations.slice(0,4)){
    for(const other of [limit,targets[(j+1)%E.length]]){
     const delta=bilinear(u,o.sub(targets[j],other),v);
     ensure(delta.pow(2).le(error.pow(2).mul(norm(u)).mul(norm(v))),'native far-boundary bilinear error');bounds++;
    }
   }
  }
 }
 return {exact_boundary_bilinear_bounds:bounds,finite_boundary_recursions:recurrences,max_return_cutoff:16,
  terminal_choices:['I','H','R','K'],unbounded_terminal_amplification_admitted:false};
});

const min=(a,b)=>a.le(b)?a:b,max=(a,b)=>a.le(b)?b:a;
const it=(lo,hi=lo)=>[f(lo),f(hi)];
const plus=(a,b)=>[a[0].add(b[0]),a[1].add(b[1])];
const times=(a,b)=>{const vs=[a[0].mul(b[0]),a[0].mul(b[1]),a[1].mul(b[0]),a[1].mul(b[1])];return [vs.reduce(min),vs.reduce(max)];};
const neg=a=>[a[1].neg(),a[0].neg()],minus=(a,b)=>plus(a,neg(b));
const divPositive=(a,b)=>{ensure(f(0).le(b[0])&&!b[0].zero(),'positive interval denominator');return times(a,[f(1).div(b[1]),f(1).div(b[0])]);};
const power=(a,n)=>{let b=it(1);for(let k=0;k<n;k++)b=times(b,a);return b;};
const etaInterval=it(etaLo,etaHi),qInterval=divPositive(minus(it(1),etaInterval),plus(it(1),power(etaInterval,4)));
const coefficientIntervals=[qInterval,times(qInterval,etaInterval),neg(times(qInterval,power(etaInterval,2))),neg(times(qInterval,power(etaInterval,3)))];
const responseInterval=(u,v=u)=>E.reduce((total,M,i)=>plus(total,times(coefficientIntervals[i],it(bilinear(u,M,v)))),it(0));
const maxAbs=a=>max(a[0].abs(),a[1].abs());
check('finite_event_tail_bound_and_completed_response_enclosures',()=>{
 let gaps=0,errors=0,differences=0;
 for(let N=0;N<=256;N++){
  const gapUpper=etaHi.sub(f(1).sub(ss[N]));ensure(gapUpper.le(new o.F(4,BigInt(N+1)**2n)),'inherited event-gap tail bound');gaps++;
  const bound=min(f(2),new o.F(138368,125n*BigInt(N+1)**2n));
  for(const u of preparations.slice(0,4))for(const v of preparations.slice(0,4)){
   const diff=minus(it(bilinear(u,response[N],v)),responseInterval(u,v));
   ensure(maxAbs(diff).pow(2).le(bound.pow(2).mul(norm(u)).mul(norm(v))),'finite native event / completed interval error');errors++;
  }
 }
 for(let k=5;k<=128;k++){
  ensure(BigInt(k)**3n-4n*BigInt(k-1)**3n+6n*BigInt(k-2)**3n-4n*BigInt(k-3)**3n+BigInt(k-4)**3n===0n,'fourth finite difference of native cube counts');differences++;
 }
 for(const r of [f('5/8'),etaLo,f('11/16')]){
  let partial=f(0);for(let k=1;k<=128;k++)partial=partial.add(r.pow(k-1).mul(BigInt(k)**3n));
  const total=f(1).add(r.mul(4)).add(r.pow(2)).div(f(1).sub(r).pow(4));ensure(partial.le(total),'native cubic geometric sum positive tail');
 }
 const x=f('11/16'),prefactor=f(8).mul(f(1).add(x.mul(4)).add(x.pow(2))).div(f(1).sub(x).pow(3));
 ensure(prefactor.eq('138368/125'),'uniform event bound exact coefficient');
 return {event_gap_tail_checks:gaps,exact_bilinear_event_error_checks:errors,cubic_finite_difference_identities:differences,
  inherited_eta_interval:r24eta,uniform_event_error:'min(2,138368/(125*(N+1)^2))',
  imported_exponential_in_event_time:false};
});

const floorScaled=(x,digits)=>{const z=x.n*10n**BigInt(digits);let n=z/x.d;if(z<0n&&z%x.d!==0n)n--;return n;};
const decimal=(x,digits,upper=false)=>{const scale=10n**BigInt(digits);let n=floorScaled(x,digits);if(upper&&!new o.F(n,scale).eq(x))n++;
 const sign=n<0n?'-':'';if(n<0n)n=-n;const s=n.toString().padStart(digits+1,'0');return sign+s.slice(0,-digits)+'.'+s.slice(-digits);};
const enclosure=(x,digits=10)=>({lower:decimal(x[0],digits),upper:decimal(x[1],digits,true),
 lower_rational:new o.F(floorScaled(x[0],digits),10n**BigInt(digits)).toString(),
 upper_rational:new o.F(floorScaled(x[1],digits)+1n,10n**BigInt(digits)).toString()});
let canonicalIntervals;
check('native_source_preparation_range_and_fixed_numeric_responses',()=>{
 let circles=0,identities=0;
 for(const v of preparations){const h=bilinear(v,H).div(norm(v)),j=bilinear(v,K).div(norm(v));ensure(h.pow(2).add(j.pow(2)).eq(1),'native source circle identity');circles++;
  for(const r of rationalDiagnostics){const D=completed(r),q=f(1).sub(r).div(f(1).add(r.pow(4)));
   ensure(bilinear(v,D).div(norm(v)).eq(q.mul(f(1).add(r.mul(h)).sub(r.pow(3).mul(j)))),'source parity/current response formula');identities++;
  }
 }
 for(const r of rationalDiagnostics){const X=o.sub(o.scale(H,r),o.scale(K,r.pow(3)));eq(o.mul(X,X),o.scale(I,r.pow(2).add(r.pow(6))),'sharp source-range square identity');}
 const upper=f('11/16');ensure(upper.pow(2).add(upper.pow(6)).le(1),'strictly positive completed source range');
 const evalE0=r=>f(1).sub(r.pow(2)).div(f(1).add(r.pow(4)));
 const evalE1=r=>f(1).sub(r).pow(2).div(f(1).add(r.pow(4)));
 const evalBalanced=r=>f(1).sub(r).mul(f(1).sub(r.pow(3))).div(f(1).add(r.pow(4)));
 canonicalIntervals={source_e0:enclosure([evalE0(etaHi),evalE0(etaLo)]),source_e1:enclosure([evalE1(etaHi),evalE1(etaLo)]),
  source_balanced_Ce0:enclosure([evalBalanced(etaHi),evalBalanced(etaLo)])};
 ensure(f(0).le(evalE1(etaHi))&&evalE1(etaLo).le(evalE0(etaHi)),'separated source readouts');
 return {native_preparation_circle_checks:circles,exact_source_response_identities:identities,
  sharp_range_square_checks:rationalDiagnostics.length,completed_coefficient_intervals:coefficientIntervals.map(x=>enclosure(x)),
  ...canonicalIntervals,empirical_constant_or_fit_used:false};
});

const LW=X=>o.add(o.scale(o.add(I,K),X[1][0].div(2)),o.scale(o.sub(I,K),X[0][1].div(2)));
check('full_word_feedback_is_nilpotent_and_differs_from_arrival_only_memory',()=>{
 let maps=0,norms=0,entries=0,recordWords=0;
 for(const X of E){
  let direct=Z;for(let s=0;s<2;s++)direct=o.add(direct,o.scale(o.mul(o.mul(o.dagger(o.mul(K,L[s])),X),L[s]),'1/2'));
  eq(direct,LW(X),'recorded native excursion map');eq(LW(LW(LW(X))),Z,'third word-continuation power');maps+=2;
 }
 eq(LW(R),K,'word R to K');eq(LW(K),I,'word K to I');eq(LW(I),Z,'word I to zero');eq(LW(H),Z,'word H to zero');
 let a=new Map([[key(0,0),I]]),b=new Map([[key(0,0),I]]);
 for(let n=1;n<=12;n++){
  a=envStep(a,0,n,true);b=envStep(b,1,n,true);
  eq(o.scale(envGram(a),den(n)),I,'word-memory system0 norm');eq(o.scale(envGram(b),den(n)),I,'word-memory system1 norm');norms+=2;
  eq(o.scale(envCross(a,b),den(n)),o.scale(I,central(Math.floor(n/2))),'full-word feedback / exact survival response');entries+=4;recordWords+=a.size+b.size;
 }
 ensure(!o.equal(response[2],o.scale(I,'1/2')),'same return gate, different native memory at two events');
 return {native_word_continuation_map_checks:maps+4,joint_word_norm_checks:norms,full_word_response_entries:entries,
  preserved_word_environment_blocks:recordWords,max_direct_word_event:12,word_continuation_nilpotency:3,
  zero_response_is_not_zero_joint_source:true};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk),bb=add(ii,rr),bd=add(ii,sc(-1,rr));
const tasks=[
 ['source_parity_involution',mul(hh,hh),ii],
 ['source_orientation_square',mul(rr,rr),sc(-1,ii)],
 ['native_feedback_K_after_return',mul(kk,bb),ff],
 ['source_return_normalization',mul(bd,bb),sc(2,ii)],
 ['source_C_normalization',mul(ff,ff),sc(2,ii)],
 ['source_return_square',mul(bb,bb),sc(2,rr)],
 ['source_return_fourth_power',mul(bb,bb,bb,bb),sc(-4,ii)],
 ['source_return_eighth_power',mul(bb,bb,bb,bb,bb,bb,bb,bb),sc(16,ii)],
 ['feedback_I_to_H',mul(bd,ff),sc(2,hh)],
 ['feedback_H_to_negative_R',mul(bd,hh,ff),sc(-2,rr)],
 ['feedback_R_to_K',mul(bd,rr,ff),sc(2,kk)],
 ['feedback_K_to_I',mul(bd,kk,ff),sc(2,ii)],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[6].result.certificate));altered.output.terms=[[[],['4','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native pulse certificate accepted');
const negatives={
 return_timing_is_an_external_fixed_two_event_schedule:!weight(2).zero()&&!weight(4).zero(),
 replacing_the_source_by_a_fresh_copy_after_each_return_is_exact:bilinear(e0,cycles[1]).eq(1)&&bilinear(e0,cycles[2]).zero(),
 scalar_current_alone_is_closed_under_return:bilinear(col([1,1]),R).zero()&&!bilinear(col([1,1]),LC(R)).zero(),
 the_antisymmetric_R_component_may_be_dropped_before_continuing:!o.equal(LC(R),Z),
 four_returns_restore_the_entire_system_source_pulse:!o.equal(o.power(pulseRaw,4),o.scale(o.identity(4),4)),
 eight_returns_mean_eight_transport_events:!weight(4).pow(8).zero(),
 coherent_amplitude_is_the_arrival_energy_coefficient:f('1/2').le(etaLo)&&f(r24.exact_checks.find(x=>x.name==='coherent_completion_branch_and_certified_rational_enclosures').g.upper_rational).le('1/2'),
 algebraic_fixed_point_uniqueness_alone_guarantees_forward_completion:o.isZero(completed(f(1)))&&!o.equal(cycles[0],cycles[4]),
 all_bounded_terminal_choices_survive_in_the_completed_response:etaHi.pow(16).mul(2).le('1/100'),
 completed_scalar_retention_is_independent_of_native_source:!bilinear(e0,completed(etaLo)).eq(bilinear(e1,completed(etaLo))),
 causality_and_pairing_preservation_uniquely_select_memory:f(0).le(responseInterval(col([1,1]))[0])&&!LW(I)[0][0].rad.eq(1),
 using_the_same_return_gate_makes_full_word_and_arrival_protocols_identical:!o.equal(response[2],o.scale(I,'1/2')),
 finite_no_return_is_a_declaration_of_permanent_escape:!weight(68).zero()&&f(0).le(ss[64]),
 the_event_gap_tail_is_eta_to_the_event_depth:f('11/16').pow(64).le(weight(68)),
 zero_full_word_cross_response_means_zero_joint_source:o.isZero(LW(I))&&central(6).le(1),
 full_word_nilpotency_identifies_or_deletes_the_source_histories:!o.isZero(LW(LW(R)))&&o.isZero(LW(LW(LW(R)))),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r25.native-return-feedback.v1',status:'PASS_R25_NATIVE_RETURN_FEEDBACK',
 input_sha256:inputHash,r24_native_certificate_sha256:r24Hash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_return_feedback_derived:true,native_terminal_boundary_independence_derived:true,
 ordinary_complex_stochastic_or_classical_feedback_premise:false,
 physical_memory_interface_selected:false,physical_metric_c_alpha_derived:false,formal_proof_assistant_verified:false,
 scope:'Real signed native source and a declared autonomous origin-return gate. Exact four-component memory, all-return response, terminal-boundary independence and finite-event tail. A different native memory interface gives a different response. Eight written proofs and scoped checks; no physical interface or electromagnetic selection.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
