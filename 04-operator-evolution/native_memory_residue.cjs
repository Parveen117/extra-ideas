#!/usr/bin/env node
'use strict';
// R26 application adapter. All scalar/matrix arithmetic and symbolic replay
// use the unchanged, pinned RKF engine; this file is not a second engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw Error('R26 input pin mismatch');
const ensure=(yes,msg)=>{if(!yes)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const I=o.identity(2),I4=o.identity(4),Z=o.zeros(2),Z4=o.zeros(4);
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange),R=o.mul(K,H);
const F=o.add(H,K),B=o.add(I,R),D=o.sub(K,H),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const E=[I,H,R,K],col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const preps=[e0,e1,col([1,1]),col([1,-2]),col(['1/3','2/5'])];
const pow2=n=>1n<<BigInt(n),den=n=>new o.F(1,pow2(n));
const bilinear=(u,M,v=u)=>o.mul(o.mul(o.dagger(u),M),v)[0][0].rad;
const addAt=(field,k,M)=>{const value=o.add(field.get(k)||o.zeros(M.length,M[0].length),M);if(o.isZero(value))field.delete(k);else field.set(k,value);};
const equalFields=(a,b,rows,cols,label)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||o.zeros(rows,cols),b.get(key)||o.zeros(rows,cols),label+' at '+key);};
const fieldEnergy=a=>[...a.values()].reduce((z,M)=>z.add(o.energy(M)),f(0));
const fieldScale=(a,c)=>new Map([...a].map(([x,M])=>[x,o.scale(M,c)]));
const fieldSub=(a,b)=>{const c=new Map(a);for(const [x,M]of b)addAt(c,x,o.scale(M,-1));return c;};
const shift=(a,d)=>new Map([...a].map(([x,M])=>[x+d,M]));
const jointStep=(field,arrow)=>{const next=new Map();for(const [x,M]of field)for(let d=0;d<2;d++)addAt(next,x+1-2*d,o.mul(o.kron(L[d],arrow(x,d)),M));return next;};
const bareStep=field=>{const next=new Map();for(const [x,M]of field)for(let d=0;d<2;d++)addAt(next,x+1-2*d,o.mul(L[d],M));return next;};
const hk=(_x,d)=>d?K:H,commuting=(_x,d)=>d?o.scale(H,-1):H;
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

check('flat_reciprocal_and_commuting_records_cancel_in_source_and_feedback',()=>{
 const G=x=>o.mul(o.power(H,((x%2)+2)%2),o.power(K,((Math.floor(x/2)%2)+2)%2));
 const reciprocal=(x,d)=>o.mul(G(x+1-2*d),o.dagger(G(x)));
 let bare=new Map([[0,I]]),rec=new Map([[0,I4]]),com=new Map([[0,I4]]),blocks=0;
 for(let n=1;n<=18;n++){
  bare=bareStep(bare);rec=jointStep(rec,reciprocal);com=jointStep(com,commuting);
  const er=new Map(),ec=new Map();for(const [x,M]of bare){er.set(x,o.kron(M,G(x)));ec.set(x,o.kron(M,o.mul(o.power(H,(n+x)/2),o.power(o.scale(H,-1),(n-x)/2))));}
  equalFields(rec,er,4,4,'reciprocal record endpoint factor');equalFields(com,ec,4,4,'commuting record endpoint factor');blocks+=rec.size+com.size;
 }
 let flagBlocks=0;
 for(const arrow of [reciprocal,commuting])for(let controller=0;controller<2;controller++){
  let a=new Map([['0,0',I]]),b=new Map([['0,0',I4]]);
  for(let n=1;n<=8;n++){
   const aa=new Map(),bb=new Map();
   for(const [key,M]of a){const [x,mark]=key.split(',').map(Number);for(let d=0;d<2;d++){
    const y=x+1-2*d,k=y+','+(mark|((y===0?1:0)<<(n-1)));let v=o.mul(L[d],M);if(controller&&y===0)v=o.mul(K,v);addAt(aa,k,v);
   }}
   for(const [key,M]of b){const [x,mark]=key.split(',').map(Number);for(let d=0;d<2;d++){
    const y=x+1-2*d,k=y+','+(mark|((y===0?1:0)<<(n-1)));let v=o.mul(o.kron(L[d],arrow(x,d)),M);if(controller&&y===0)v=o.mul(o.kron(K,I),v);addAt(bb,k,v);
   }}
   a=aa;b=bb;const expected=new Map();for(const [key,M]of a){const x=Number(key.split(',')[0]);const g=arrow===reciprocal?G(x):o.mul(o.power(H,(n+x)/2),o.power(o.scale(H,-1),(n-x)/2));expected.set(key,o.kron(M,g));}
   equalFields(b,expected,4,4,'flat record with retained feedback flags');flagBlocks+=b.size;
  }
 }
 return {endpoint_operator_blocks:blocks,retained_feedback_blocks:flagBlocks,max_unflagged_event:18,max_feedback_event:8,
  reciprocal_edge_arrows_need_not_commute:true,multiple_address_preparation_equivalence_claimed:false};
});

check('two_event_commutator_square_and_nonremovable_HK_readout',()=>{
 const rotation=o.matrix([['3/5','-4/5'],['4/5','3/5']]);let bounds=0,frames=0;
 for(const A of [I,H,K,R,rotation])for(const BB of [I,H,K,R,rotation])for(const mu of preps){
  const arrow=(_x,d)=>d?BB:A,seed=o.kron(e0,mu);let field=new Map([[0,seed]]);field=jointStep(jointStep(field,arrow),arrow);
  const v=field.get(0)||o.zeros(4,1),energy=o.energy(mu),delta=o.energy(o.mul(o.sub(o.mul(A,BB),o.mul(BB,A)),mu)).div(energy);
  ensure(f(0).le(delta)&&delta.le(4),'native commutator-square range');
  ensure(o.energy(v).div(energy.mul(4)).eq('1/2'),'two-event return energy');
  ensure(bilinear(v,o.kron(K,I)).div(energy.mul(4)).eq(f('1/2').sub(delta.div(4))),'commutator current identity');bounds++;
  const diff=o.sub(o.mul(A,BB),o.mul(BB,A)),framed=o.mul(o.mul(rotation,diff),o.dagger(R));
  ensure(o.energy(o.mul(framed,o.mul(R,mu))).eq(o.energy(o.mul(diff,mu))),'native endpoint-frame invariance');frames++;
 }
 eq(o.add(o.mul(H,K),o.mul(K,H)),Z,'record anticommutation');eq(o.mul(o.dagger(o.mul(H,K)),o.mul(K,H)),o.scale(I,-1),'relative record arrow');
 let field=new Map([[0,o.kron(e0,e0)]]);for(let n=1;n<=3;n++)field=jointStep(field,hk);
 const energies=[-3,-1,1,3].map(x=>o.energy(field.get(x)).div(8));
 ensure(energies.every((value,i)=>value.eq(['1/8','5/8','1/8','1/8'][i])),'third-event observable residue');
 return {commutator_response_checks:bounds,frame_invariance_checks:frames,HK_defect:'4',HK_return_current:'-1/2',
  third_event_addresses:[-3,-1,1,3],third_event_energies:energies.map(String)};
});

check('normal_ordered_memory_matches_full_joint_signed_stencil',()=>{
 let joint=new Map([[0,I4]]),phi=new Map([[0,I]]),blocks=0;
 for(let n=0;n<32;n++){
  joint=jointStep(joint,hk);const next=new Map();for(const [x,M]of phi){addAt(next,x+1,o.mul(L[0],M));addAt(next,x-1,o.scale(o.mul(L[1],M),((n+x)/2)%2?-1:1));}phi=next;
  const expected=new Map();for(const [x,M]of phi)expected.set(x,o.kron(M,o.mul(o.power(H,(n+1+x)/2),o.power(K,(n+1-x)/2))));
  equalFields(joint,expected,4,4,'normal-ordered signed source stencil');blocks+=joint.size;
 }
 return {exact_joint_operator_blocks:blocks,max_event:32,record_preparation_parameter_required:false};
});

check('native_quartic_shift_identity_and_sharp_two_event_gap',()=>{
 const M=new Map([[1,o.kron(L[0],H)],[-1,o.kron(L[1],K)]]);
 // Finite convolution of native matrices is a coefficient check, not a new scalar engine.
 const compose=(a,b)=>{const c=new Map();for(const [i,A]of a)for(const [j,BB]of b)addAt(c,i+j,o.mul(A,BB));return c;};
 const M2=compose(M,M),M4=compose(M2,M2),right=new Map();for(const [x,A]of M2){addAt(right,x+2,A);addAt(right,x-2,A);}addAt(right,0,o.scale(I4,-4));
 equalFields(M4,right,4,4,'all Laurent coefficients of native quartic');
 let norms=0,sharp=0;
 const seeds=[o.kron(e0,e0),o.kron(col([1,2]),col([2,-1])),col([1,2,-3,4])];
 for(const seed of seeds)for(const m of [1,2,3,8,16,32])for(const alternate of [false,true]){
  const packet=new Map();for(let j=1;j<=m;j++)packet.set(2*j,o.scale(seed,alternate&&j%2?-1:1));
  const V2=fieldScale(jointStep(jointStep(packet,hk),hk),'1/2'),lhs=fieldEnergy(fieldSub(V2,packet)),energy=fieldEnergy(packet);
  const delta=fieldEnergy(fieldSub(shift(packet,2),packet));ensure(lhs.eq(energy.add(delta.div(2))),'native defect norm identity');norms++;
  ensure(lhs.div(energy).eq(alternate?f(3).sub(new o.F(1,m)):f(1).add(new o.F(1,m))),'sharp packet sequence');sharp++;
 }
 for(const m of [1,2,3,8,16,32]){
  const packet=new Map();for(let j=1;j<=m;j++)packet.set(2*j,e0);
  const U2=fieldScale(bareStep(bareStep(packet)),'1/2');ensure(fieldEnergy(fieldSub(U2,packet)).div(fieldEnergy(packet)).eq(new o.F(1,m)),'bare native gap is zero');
 }
 return {quartic_coefficient_addresses:[...new Set([...M4.keys(),...right.keys()])].sort((a,b)=>a-b),
  finite_field_gap_equalities:norms,sharp_packet_checks:sharp,bare_packet_checks:6,sharp_infimum:'1',sharp_supremum:'3',
  defect:'squared norm of (V^2-I) divided by source norm',physical_mass_interpretation_proved:false};
});

// Independent integer recurrence from formal coefficient comparison in R26.6.
const a=[0n,1n,1n];
for(let k=3;k<=4096;k++){
 const n=BigInt(k),value=4n*(4n*n-3n)*a[k-1]-16n*(4n*n-9n)*a[k-2]+64n*(2n*n-6n)*a[k-3];
 ensure(value%(2n*n)===0n,'integral derived coefficient');a.push(value/(2n*n));
}
const t=m=>m===1?1n:m%2===0?a[m/2]:-2n*a[(m-1)/2];
const weight=n=>n>0&&n%2===0?new o.F(t(n/2)**2n,pow2(n-1)):f(0);
const ww=Array.from({length:129},(_,n)=>weight(n)),ss=[f(1)];
for(let n=1;n<=128;n++)ss.push(ss[n-1].sub(ww[n]));

check('first_return_joint_operators_match_new_signed_grammar',()=>{
 let live=new Map([[0,I4]]),returns=Z4,operatorEntries=0,grams=0;const rows=[];
 for(let n=1;n<=64;n++){
  live=jointStep(live,hk);const ret=live.get(0)||Z4;live.delete(0);
  const expected=n%2?Z4:o.scale(o.kron((n/2)%2?D:B,o.power(R,n/2)),t(n/2));
  eq(ret,expected,'native first-return raw operator');operatorEntries+=16;
  const gram=o.scale(o.mul(o.dagger(ret),ret),den(n));eq(gram,o.scale(I4,ww[n]),'derived return norm is scalar');grams++;
  returns=o.add(returns,gram);const liveGram=[...live.values()].reduce((sum,A)=>o.add(sum,o.mul(o.dagger(A),A)),Z4);
  eq(o.add(returns,o.scale(liveGram,den(n))),I4,'first-return survivor telescope');grams++;
  if(n%2===0)rows.push({event:n,signed_coefficient:t(n/2).toString(),weight:String(ww[n])});
 }
 return {joint_return_operator_entries:operatorEntries,return_and_survivor_gram_checks:grams,max_event:64,first_returns:rows};
});

let pLo,pHi,returnBracket;
const floorScaled=(x,digits)=>{const n=x.n*10n**BigInt(digits);let q=n/x.d;if(n<0n&&n%x.d)q--;return q;};
const outward=(lo,hi,digits=12)=>{
 const scale=10n**BigInt(digits),low=floorScaled(lo,digits),high=floorScaled(hi,digits)+1n;
 const decimal=n=>{const sign=n<0n?'-':'';if(n<0n)n=-n;const s=n.toString().padStart(digits+1,'0');return sign+s.slice(0,-digits)+'.'+s.slice(-digits);};
 return {lower:decimal(low),upper:decimal(high),lower_rational:String(new o.F(low,scale)),upper_rational:String(new o.F(high,scale))};
};
check('native_coefficient_completion_and_exact_return_tail_witness',()=>{
 let quadratic=0,grammar=0;
 for(let k=1;k<=128;k++){
  let rhs=0n;for(let i=1;i<k;i++)rhs+=a[i]*a[k-i];for(let i=1;i<k-1;i++)rhs-=4n*a[i]*a[k-1-i];if(k===1)rhs=1n;
  ensure(a[k]===rhs,'linear versus independent quadratic coefficient recurrence');quadratic++;
 }
 const tv=Array.from({length:129},(_,m)=>m?t(m):0n),square=Array(129).fill(0n);
 for(let i=1;i<=128;i++)for(let j=1;i+j<=128;j++)square[i+j]+=tv[i]*tv[j];
 for(let n=0;n<=128;n++){
  const coefficient=(square[n]||0n)+2n*(square[n-1]||0n)-(tv[n]||0n)-2n*(tv[n-1]||0n)-4n*(tv[n-2]||0n)+(n===1?1n:n===2?2n:0n);
  ensure(coefficient===0n,'new signed first-return polynomial grammar');grammar++;
 }
 let Zcount=0n,denom=1n,Z128=0n,pairs=0;
 for(let k=1;k<=4096;k++){
  denom*=16n;Zcount=16n*Zcount+a[k]*a[k];if(k===128)Z128=Zcount;
  ensure(a[k]**2n*BigInt(k+1)**3n<=64n*denom,'finite coefficient bound witness');
  if(k<=31){ensure(weight(4*k).eq(weight(4*k+2)),'paired even/odd return weights');pairs++;}
 }
 ensure(2000n*Z128<141n*(16n**128n),'exact strict upper return-mass witness');
 ensure(new o.F(128,129n**2n).le('9/1000'),'M128 tail gap to four fifths');
 const partial=f('1/2').add(new o.F(4n*Zcount,denom)),tail=new o.F(128,4097n**2n);
 returnBracket=outward(partial,partial.add(tail));pLo=f(returnBracket.lower_rational);pHi=f(returnBracket.upper_rational);
 ensure(f('3/4').le(pLo)&&!pLo.eq('3/4')&&pHi.le('4/5'),'derived strict return bracket');
 return {independent_quadratic_coefficients:quadratic,first_return_grammar_coefficients:grammar,finite_coefficient_bound_checks:4096,
  equal_return_weight_pairs:pairs,Z128:Z128.toString(),integer_witness:'2000 Z128 < 141 * 16^128',
  return_mass:{...returnBracket,partial_count:4096,tail_upper_rational:String(tail),partial_rational_sha256:hash(String(partial)),
   all_depth_tail:'0 <= p-p_M <= 128/(M+1)^2',written_native_tail_proof:'R26.6'},
  floating_point_or_fitted_input_used:false};
});

check('coherent_return_map_has_native_sqrt_three_normalization',()=>{
 eq(o.mul(o.dagger(B),B),o.scale(I,2),'even pulse norm');eq(o.mul(D,D),o.scale(I,2),'odd pulse norm');
 eq(o.mul(o.dagger(B),D),o.scale(K,2),'first cross pulse');eq(o.mul(D,B),o.scale(K,2),'second cross pulse');
 const X=o.kron(B,I),Y=o.kron(D,R);
 eq(o.add(o.mul(o.dagger(X),Y),o.mul(o.dagger(Y),X)),Z4,'coherent cross cancellation');
 // Write E=(1-r)/4, O=(1+r)/4. All coefficients are checked before r^2=3.
 const A0=o.scale(o.add(X,Y),'1/4'),A1=o.scale(o.sub(Y,X),'1/4');
 eq(o.add(o.mul(o.dagger(A0),A0),o.scale(o.mul(o.dagger(A1),A1),3)),I4,'native sqrt3 constant coefficient');
 eq(o.add(o.mul(o.dagger(A0),A1),o.mul(o.dagger(A1),A0)),Z4,'native sqrt3 linear coefficient');
 ensure(pHi.le('4/5'),'coherent isometry is not unit retained-arrival energy');
 // Lambda=(1+iota*sqrt3)/2: norm and geometric denominator use only iota^2=-1.
 ensure(new o.Cut(0,1).norm2().eq(1)&&o.IOTA.mul(o.IOTA).eq(-1),'derived native scalar turn');
 ensure(f('1/4').add(f('3/4')).eq(1),'lambda and one-minus-lambda norms');
 return {exact_source_and_joint_identities:7,native_radical_relation:'r^2=3, r>0',
  E:'(1-sqrt(3))/4',O:'(1+sqrt(3))/4',coherent_fold_is_norm_preserving_erasure_of_arrival_flags:false};
});

const Le=X=>o.scale(o.mul(o.mul(o.dagger(B),X),F),'1/2');
const Lo=X=>o.scale(o.mul(o.mul(D,X),o.dagger(B)),'1/2');
const averaged=(X,r)=>o.add(o.scale(Le(X),r.sub('1/2').div(2)),o.scale(Lo(X),r.add('1/2').div(2)));
const response=[I];
for(let n=1;n<=128;n++){
 let value=o.scale(I,ss[n]);for(let j=2;j<=n;j+=2)value=o.add(value,o.scale((j%4?Lo:Le)(response[n-j]),ww[j]));response.push(value);
}
check('return_parity_fixes_feedback_imbalance_and_fourth_power',()=>{
 const even=[H,o.scale(R,-1),K,I],odd=[o.scale(H,-1),R,K,I];
 for(let i=0;i<4;i++){eq(Le(E[i]),even[i],'even return cross map');eq(Lo(E[i]),odd[i],'odd return cross map');}
 eq(o.mul(K,B),F,'even controlled pulse');eq(o.mul(K,D),o.dagger(B),'odd controlled pulse');
 // M=M0+p*M1: polynomial coefficient identity for all p, not fitted samples.
 const map0=X=>o.scale(o.sub(Lo(X),Le(X)),'1/4'),map1=X=>o.scale(o.add(Le(X),Lo(X)),'1/2');let coeffs=0;
 for(const X of E){let pol=[X];for(let k=0;k<4;k++){const next=Array.from({length:pol.length+1},()=>Z);for(let j=0;j<pol.length;j++){next[j]=o.add(next[j],map0(pol[j]));next[j+1]=o.add(next[j+1],map1(pol[j]));}pol=next;}
  for(let j=0;j<5;j++){eq(pol[j],j===2?o.scale(X,'-1/4'):Z,'all-parameter feedback fourth-power coefficient');coeffs++;}}
 for(let m=1;m<=31;m++){
  let evenWeight=f(0),oddWeight=f('1/2');for(let k=1;k<=m;k++){evenWeight=evenWeight.add(weight(4*k));oddWeight=oddWeight.add(weight(4*k+2));}
  ensure(evenWeight.sub(oddWeight).eq('-1/2'),'paired-tail cancellation fixes imbalance');
 }
 return {pulse_and_basis_identities:10,all_parameter_polynomial_coefficients:coeffs,exact_finite_imbalance_checks:31,
  completed_even_minus_odd_weight:'-1/2',fourth_power:'-p^2/4 times identity'};
});

const envStep=(field,controller,n)=>{
 const next=new Map();for(const [key,M]of field){const [x,mark]=key.split(',').map(Number);for(let d=0;d<2;d++){
  const y=x+1-2*d,k=y+','+(mark|((y===0?1:0)<<(n-1)));let v=o.mul(o.kron(L[d],d?K:H),M);if(controller&&y===0)v=o.mul(o.kron(K,I),v);addAt(next,k,v);
 }}return next;
};
const envGram=field=>[...field.values()].reduce((sum,M)=>o.add(sum,o.mul(o.dagger(M),M)),Z4);
const envCross=(a,b)=>{let X=Z4;for(const [key,M]of a)if(b.has(key))X=o.add(X,o.mul(o.dagger(M),b.get(key)));return X;};
let snapshots;
check('full_joint_feedback_matches_parity_recurrence_with_record_retained',()=>{
 let a0=new Map([['0,0',I4]]),a1=new Map([['0,0',I4]]),blocks=0;snapshots=[[a0,a1]];
 for(let n=1;n<=12;n++){
  a0=envStep(a0,0,n);a1=envStep(a1,1,n);snapshots.push([a0,a1]);
  eq(o.scale(envGram(a0),den(n)),I4,'system0 full environment norm');eq(o.scale(envGram(a1),den(n)),I4,'system1 full environment norm');
  eq(o.scale(envCross(a0,a1),den(n)),o.kron(response[n],I),'full joint feedback versus parity renewal');blocks+=a0.size+a1.size;
 }
 eq(response[2],o.scale(o.sub(I,H),'1/2'),'first-return feedback signs');
 return {joint_norm_equalities:24,joint_response_entries:192,retained_environment_blocks:blocks,max_direct_event:12,max_recurrence_event:128,
  auxiliary_record_is_not_reset:true,joint_relative_response_has_record_identity_factor:true};
});

check('literal_source_letters_match_direction_memory_and_feedback',()=>{
 let prefixes=0,blocks=0;
 for(let controller=0;controller<2;controller++)for(const source of [e0,e1])for(const record of [e0,e1]){
  const seed=o.kron(source,record);let histories=[{x:0,mark:0,v:seed,word:''}];
  for(let n=1;n<=8;n++){
   histories=histories.flatMap(h=>[[H,'H'],[K,'K']].map(([letter,tag])=>{
    let v=o.mul(o.kron(letter,I),h.v);const d=o.isZero(o.mul(o.kron(P,I),v))?1:0;
    v=o.mul(o.kron(I,d?K:H),v);const x=h.x+1-2*d,mark=h.mark|((x===0?1:0)<<(n-1));if(controller&&x===0)v=o.mul(o.kron(K,I),v);
    prefixes++;return {x,mark,v,word:h.word+tag};
   }));
   ensure(new Set(histories.map(h=>h.word)).size===2**n,'literal source tags retain all words');
   const aggregate=new Map();for(const h of histories)addAt(aggregate,h.x+','+h.mark,h.v);
   const expected=new Map([...snapshots[n][controller]].map(([key,M])=>[key,o.mul(M,seed)]));equalFields(aggregate,expected,4,1,'literal source/record histories');blocks+=aggregate.size;
  }
 }
 return {literal_HK_prefixes:prefixes,exact_aggregate_blocks:blocks,max_literal_event:8,
  source_letter_tag_is_not_the_record_direction_control:true};
});

const completed=r=>o.scale(o.add(o.add(I,o.scale(H,'-1/2')),o.add(o.scale(R,'-1/4'),o.scale(K,r.div(-4)))),f(1).sub(r).div(f(1).add(r.pow(2).div(4))));
const normalizedResponses=r=>{
 const q=f(1).sub(r).div(f(1).add(r.pow(2).div(4)));return [q.div(2),q.mul('3/2'),q.mul(f(1).sub(r.div(4)))];
};
check('completed_feedback_identity_ratio_three_and_certified_enclosures',()=>{
 const N0=o.sub(o.sub(I,o.scale(H,'1/2')),o.scale(R,'1/4')),N1=o.scale(K,'-1/4');
 const map0=X=>o.scale(o.sub(Lo(X),Le(X)),'1/4'),map1=X=>o.scale(o.add(Le(X),Lo(X)),'1/2');
 eq(o.sub(N0,map0(N0)),I,'inverse polynomial constant');eq(o.sub(o.sub(N1,map0(N1)),map1(N0)),Z,'inverse polynomial linear');eq(o.scale(map1(N1),-1),o.scale(I,'1/4'),'inverse polynomial quadratic');
 let identities=0,sourceChecks=0;
 for(const r of [f('3/4'),pLo,pHi,f('4/5')]){
  const value=completed(r);eq(value,o.add(o.scale(I,f(1).sub(r)),averaged(value,r)),'completed fixed point');identities++;
  for(const v of preps){const norm=o.energy(v),h=bilinear(v,H).div(norm),j=bilinear(v,K).div(norm),q=f(1).sub(r).div(f(1).add(r.pow(2).div(4)));
   ensure(bilinear(v,value).div(norm).eq(q.mul(f(1).sub(h.div(2)).sub(r.mul(j).div(4)))),'source readout');sourceChecks++;}
  ensure(bilinear(e1,value).eq(bilinear(e0,value).mul(3))&&!bilinear(e0,value).zero(),'exact channel ratio three');
 }
 const lo=normalizedResponses(pHi),hi=normalizedResponses(pLo),enclosures={};
 for(const [i,key]of ['source_e0','source_e1','source_balanced_Ce0'].entries())enclosures[key]=outward(lo[i],hi[i]);
 // Norm bound for arbitrary terminal actions is written in R26.8. These
 // exact native-arrow diagnostics independently exercise its finite recurrence.
 let terminalBounds=0;
 for(const r of [f('3/4'),f('4/5')])for(const G of E){let boundary=G;
  for(let k=0;k<=16;k++){
   for(const u of preps.slice(0,3))for(const v of preps.slice(0,3)){
    const delta=bilinear(u,o.sub(boundary,completed(r)),v);ensure(delta.pow(2).le(r.pow(2*k).mul(4).mul(o.energy(u)).mul(o.energy(v))),'terminal-depth native pairing bound');terminalBounds++;
   }boundary=o.add(o.scale(I,f(1).sub(r)),averaged(boundary,r));
  }
 }
 const k=64,M=4096,N=k*(4*M+2),bound=pHi.pow(k+1).add(new o.F(128,BigInt(M+1)**2n).div(f(1).sub(pHi))).mul(2);
 ensure(bound.le('1/10000'),'finite event completion witness');
 return {all_parameter_fixed_point_coefficients:3,exact_fixed_point_diagnostics:identities,source_formula_checks:sourceChecks,
  terminal_bilinear_bound_checks:terminalBounds,...enclosures,exact_channel_ratio:'3',
  event_horizon_witness:{k,M,N,error_upper:outward(f(0),bound),formula:'2*(p^(k+1)+128/((M+1)^2*(1-p)))'},
  event_witness_is_tail_proof_evaluation_not_million_event_enumeration:true,physical_charge_or_alpha_identification_proved:false};
});

check('literal_letter_copying_requires_zero_cut_real_flag_overlap',()=>{
 const HH=o.kron(H,H),KK=o.kron(K,K),RR=o.kron(R,R),raw=o.add(HH,KK),kernel=col([1,0,0,-1]);
 eq(o.mul(HH,KK),RR,'doubled letters product');eq(o.mul(KK,HH),RR,'doubled letters commute');
 eq(o.scale(o.mul(o.dagger(raw),raw),'1/2'),o.add(I4,RR),'unflagged lift Gram');eq(o.mul(raw,kernel),o.zeros(4,1),'nonzero lift kernel');
 const flags=[[e0,e0],[e0,e1],[e0,o.scale(e0,o.IOTA)],[e0,col(['3/5','4/5'])],
  [e0,col([new o.Cut(0,'3/5'),'4/5'])]];let norms=0,passing=0;
 for(const [h,k]of flags){ensure(o.energy(h).eq(1)&&o.energy(k).eq(1),'normalized native flags');
  const lift=o.add(o.kron(HH,h),o.kron(KK,k)),overlap=o.mul(o.dagger(h),k)[0][0];
  const gram=o.scale(o.mul(o.dagger(lift),lift),'1/2');eq(gram,o.add(I4,o.scale(RR,overlap.rad)),'full flag-pairing expansion');norms++;
  ensure(o.equal(gram,I4)===overlap.rad.zero(),'iff zero cut-real flag overlap');if(overlap.rad.zero())passing++;
 }
 const turnLift=o.add(HH,o.scale(KK,o.IOTA));
 eq(o.scale(o.mul(o.dagger(turnLift),turnLift),'1/2'),I4,'one-direction native-iota repair');
 return {doubled_letter_and_kernel_identities:4,exact_flag_gram_checks:norms,pairing_preserving_flag_examples:passing,
  real_flag_rule:'orthogonal',turn_valued_flag_rule:'zero cut-real overlap',one_scalar_direction_iota_flags_distinguish_letters:false};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk),bb=add(ii,rr),bd=add(ii,sc(-1,rr)),dd=add(kk,sc(-1,hh));
const tasks=[
 ['source_H_involution',mul(hh,hh),ii],['source_K_involution',mul(kk,kk),ii],
 ['derived_R_square',mul(rr,rr),sc(-1,ii)],['native_HK_anticommutator',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['even_return_norm',mul(bd,bb),sc(2,ii)],['odd_return_norm',mul(dd,dd),sc(2,ii)],
 ['even_odd_cross',mul(bd,dd),sc(2,kk)],['odd_controlled_pulse',mul(kk,dd),bd],
 ['even_I_to_H',mul(bd,ff),sc(2,hh)],['even_H_to_negative_R',mul(bd,hh,ff),sc(-2,rr)],
 ['even_R_to_K',mul(bd,rr,ff),sc(2,kk)],['even_K_to_I',mul(bd,kk,ff),sc(2,ii)],
 ['odd_I_to_negative_H',mul(dd,bd),sc(-2,hh)],['odd_H_to_R',mul(dd,hh,bd),sc(2,rr)],
 ['odd_R_to_K',mul(dd,rr,bd),sc(2,kk)],['odd_K_to_I',mul(dd,kk,bd),sc(2,ii)],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[2].result.certificate));altered.output.terms=[[[],['1','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const negatives={
 every_reused_record_changes_source_readouts:checks[0].endpoint_operator_blocks>0,
 HK_direction_memory_is_reciprocal_flat_memory:!o.equal(o.mul(H,K),o.mul(K,H)),
 memory_only_changes_unobserved_labels:checks[1].third_event_energies[1]==='5/8',
 common_endpoint_frame_can_remove_the_relative_minus_sign:!o.equal(o.scale(I,-1),I),
 retained_HK_record_preparation_supplies_a_free_coupling:checks[2].record_preparation_parameter_required===false,
 old_R24_eta_can_be_reused_as_the_new_return_mass:f('3/4').le(pLo),
 coherent_unit_return_map_means_unit_retained_arrival_energy:pHi.le('4/5'),
 all_return_parities_have_the_same_feedback_pulse:!o.equal(Le(I),Lo(I)),
 dropping_R_before_continuation_preserves_future_response:!o.isZero(Le(R)),
 feedback_channel_ratio_requires_fitting_the_return_tail:normalizedResponses(pLo)[1].eq(normalizedResponses(pLo)[0].mul(3)),
 unflagged_coherent_literal_HK_copy_preserves_pairing:!o.isZero(o.kron(R,R)),
 real_flag_orthogonality_is_forced_in_every_native_turn_sector:o.IOTA.rad.zero()&&o.IOTA.norm2().eq(1),
 direction_control_is_identical_to_original_source_letter_control:!o.isZero(o.mul(P,o.mul(H,e0)))&&!o.isZero(o.mul(P,o.mul(K,e1))),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r26.native-memory-residue.v1',status:'PASS_R26_NATIVE_MEMORY_RESIDUE',
 input_sha256:inputHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_commutator_residue_derived:true,native_two_event_gap_derived:true,native_feedback_ratio_three_derived:true,
 ordinary_complex_noise_mass_or_classical_feedback_premise:false,
 physical_memory_interface_selected:false,physical_metric_c_alpha_derived:false,physical_mass_gap_identified:false,
 formal_proof_assistant_verified:false,
 scope:'Nine written native results with scoped exact checks: flat memory, commutator readout, signed H/K memory, two-event propagation defect, new return grammar/tail, parity feedback, exact channel ratio and literal-letter flag constraint. Interfaces remain explicit constructions; no physical mass, charge or alpha selection.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
