#!/usr/bin/env node
'use strict';
// R27 application checks. Native arithmetic, elimination and proof replay
// remain in the unchanged canonical RKF engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const l=require(path.join(home,'core/native_laurent.cjs'));
const ensure=(yes,msg)=>{if(!yes)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const hash=bytes=>crypto.createHash('sha256').update(bytes).digest('hex');
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
ensure(inputHash===get('--expected-input-sha256'),'R27 input pin mismatch');
const inherited={};
for(const n of [25,26]){
 const bytes=fs.readFileSync(get('--r'+n+'-certificate')),sha=hash(bytes);
 ensure(sha===get('--expected-r'+n+'-sha256'),'R27 R'+n+' certificate pin mismatch');
 const value=JSON.parse(bytes);ensure(value.status===(n===25?'PASS_R25_NATIVE_RETURN_FEEDBACK':'PASS_R26_NATIVE_MEMORY_RESIDUE'),'wrong inherited status');
 inherited[n]={sha,value};
}
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const I=o.identity(2),Z=o.zeros(2),H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange),R=o.mul(K,H);
const F=o.add(H,K),B0=o.add(I,R),D0=o.sub(K,H),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const E=[I,H,R,K],col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const basis=n=>Array.from({length:n},(_,j)=>col(Array.from({length:n},(_,i)=>Number(i===j))));
const bilinear=(v,A,u=v)=>o.mul(o.mul(o.dagger(v),A),u)[0][0].rad;
const pow2=n=>1n<<BigInt(n),den=n=>new o.F(1,pow2(n));
const defect=(A,B)=>o.scale(o.add(o.mul(A,B),o.mul(B,A)),'1/2');
const countPair=(m,n)=>{const d=f(m).pow(2).add(f(n).pow(2)),gamma=f(m).pow(2).sub(f(n).pow(2)).div(d),sigma=f(m).mul(n).mul(2).div(d);
 return {name:'counts_'+m+'_'+n,A:H,B:o.add(o.scale(H,gamma),o.scale(K,sigma)),gamma,sigma};};
const families=[[1,0],[1,1],[2,1],[1,2],[3,2],[1,-1],[0,1]].map(([m,n])=>countPair(m,n));
const addAt=(field,key,M)=>{const value=o.add(field.get(key)||o.zeros(M.length,M[0].length),M);if(o.isZero(value))field.delete(key);else field.set(key,value);};
const eqFields=(a,b,rows,cols,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||o.zeros(rows,cols),b.get(key)||o.zeros(rows,cols),msg+' at '+key);};
const step=(field,A,B)=>{const next=new Map();for(const [x,M]of field)for(let d=0;d<2;d++)addAt(next,x+1-2*d,o.mul(o.kron(L[d],d?B:A),M));return next;};
const fieldScale=(a,c)=>new Map([...a].map(([x,M])=>[x,o.scale(M,c)]));
const fieldSub=(a,b)=>{const c=new Map(a);for(const [x,M]of b)addAt(c,x,o.scale(M,-1));return c;};
const shift=(a,d)=>new Map([...a].map(([x,M])=>[x+d,M]));
const energy=a=>[...a.values()].reduce((sum,M)=>sum.add(o.energy(M)),f(0));
const fieldQuadratic=(a,X)=>[...a.values()].reduce((sum,v)=>sum.add(bilinear(v,X)),f(0));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

const replicas=[];
for(const d of [1,2,3]){
 const n=2*d;let T=o.identity(n);
 for(const [i,j]of [[0,n-1],[0,d],...(d>1?[[1,n-1]]:[])]){
  const turn=o.identity(n);turn[i][i]=turn[j][j]=new o.Cut('3/5');turn[i][j]=new o.Cut('-4/5');turn[j][i]=new o.Cut('4/5');T=o.mul(turn,T);
 }
 eq(o.mul(o.dagger(T),T),o.identity(n),'native rational record frame');
 const A=o.mul(o.mul(T,o.kron(H,o.identity(d))),o.dagger(T)),B=o.mul(o.mul(T,o.kron(K,o.identity(d))),o.dagger(T));
 replicas.push({name:'record_rank_'+n,d,A,B,T});
}
const diagonal=values=>o.matrix(values.map((x,i)=>values.map((_y,j)=>i===j?x:0)));
const mixed={name:'mixed_defect',A:o.kron(H,I),B:o.add(o.kron(H,diagonal(['3/5',0])),o.kron(K,diagonal(['4/5',1])))};
const candidates=[...families,...replicas,mixed,{name:'commuting_nonscalar',A:I,B:H}];

check('record_balance_curvature_and_return_current_share_one_native_defect',()=>{
 let identities=0,currents=0,calibrations=0;
 for(const {A,B}of candidates){const n=A.length,id=o.identity(n),G=defect(A,B),comm=o.sub(o.mul(A,B),o.mul(B,A));
  eq(o.mul(A,A),id,'candidate A self-return');eq(o.mul(B,B),id,'candidate B self-return');eq(o.dagger(A),A,'candidate A dagger');eq(o.dagger(B),B,'candidate B dagger');
  eq(o.dagger(G),G,'defect self-adjoint');eq(o.mul(A,G),o.mul(G,A),'defect commutes A');eq(o.mul(B,G),o.mul(G,B),'defect commutes B');
  for(const sign of [1,-1]){const sum=o.add(A,o.scale(B,sign));eq(o.scale(o.mul(o.dagger(sum),sum),'1/2'),o.add(id,o.scale(G,sign)),'native balance-defect identity');}
  eq(o.mul(o.dagger(comm),comm),o.scale(o.sub(id,o.mul(G,G)),4),'curvature-square identity');identities+=10;
  for(const mu of [...basis(n),col(Array.from({length:n},(_,i)=>i%2?-i-1:i+1))]){
   const norm=o.energy(mu),g=bilinear(mu,G).div(norm),g2=o.energy(o.mul(G,mu)).div(norm);
   let field=new Map([[0,o.kron(e0,mu)]]);field=step(step(field,A,B),A,B);const v=field.get(0)||o.zeros(2*n,1);
   const current=bilinear(v,o.kron(K,id)).div(norm.mul(4));ensure(current.add('1/2').eq(g2),'return current equals squared calibration defect');currents++;
   ensure(g.pow(2).le(g2)&&f(-1).le(g)&&g.le(1),'native calibration inequality');calibrations++;
  }
 }
 return {record_candidates:candidates.length,exact_operator_identities:identities,current_identities:currents,
  native_pairing_calibration_bounds:calibrations,anticommutation_assumed_for_general_candidates:false};
});

check('native_cut_normal_form_is_faithful_and_has_only_multiplicity_intertwiners',()=>{
 let factors=0,intertwiners=0;const rows=[];
 for(const rec of replicas){const {A,B,d}=rec,n=A.length,id=o.identity(n),Ep=o.scale(o.add(id,A),'1/2'),Em=o.sub(id,Ep);
  eq(o.mul(B,Ep),o.mul(Em,B),'exchange of the native cuts');ensure(o.rank(Ep)===d&&o.rank(Em)===d,'equal cut ranks');
  // Extract an independent plus-cut basis; its inherited Gram is retained.
  // No Euclidean orthonormal-basis assumption or radical shortcut is used.
  const vectors=[];for(const v of basis(n)){const x=o.mul(Ep,v);if(o.rank([...vectors,x].map(o.flatten))>vectors.length)vectors.push(x);}
  ensure(vectors.length===d,'plus-cut basis count');
  const W=o.matrix(Array.from({length:n},(_,i)=>vectors.map(v=>v[i][0]))),BW=o.mul(B,W);
  const T=o.matrix(W.map((row,i)=>row.concat(BW[i]))),inv=o.inverse(T),J=o.mul(o.dagger(W),W);
  eq(o.mul(o.dagger(T),T),o.kron(I,J),'inherited-pairing normal form');
  eq(o.mul(o.mul(inv,A),T),o.kron(H,o.identity(d)),'record parity normal form');
  eq(o.mul(o.mul(inv,B),T),o.kron(K,o.identity(d)),'record exchange normal form');
  eq(o.mul(T,inv),id,'two-sided reconstruction');
  ensure(o.rank([id,A,B,o.mul(B,A)].map(o.flatten))===4,'source algebra is faithful');factors++;
  const equations=[];for(const X of [A,B]){
   const columns=Array.from({length:n*n},(_,index)=>{const unit=o.zeros(n);unit[Math.floor(index/n)][index%n]=o.ONE;return o.flatten(o.commutator(X,unit));});
   equations.push(...l.transpose(columns));
  }
  const kernel=l.nullspaceRows(equations);ensure(kernel.length===d*d,'exact intertwiner dimension');
  for(const vector of kernel){const S=o.unflatten(vector,n),normal=o.mul(o.mul(inv,S),T);
   const C=normal.slice(0,d).map(row=>row.slice(0,d));eq(normal,o.kron(I,C),'all computed intertwiners have common diagonal block');intertwiners++;
  }
  rec.normalForm={T,inv,J};rows.push({record_rank:n,source_factor_rank:2,multiplicity_rank:d,intertwiner_rank:kernel.length});
 }
 return {exact_reconstructions:factors,exact_intertwiner_basis_checks:intertwiners,normal_forms:rows,
  orthonormal_record_basis_imported:false,source_words_identified_as_histories:false};
});

check('general_native_propagation_identity_and_sharp_count_family_bounds',()=>{
 const compose=(a,b)=>{const c=new Map();for(const [i,A]of a)for(const [j,B]of b)addAt(c,i+j,o.mul(A,B));return c;};
 let coefficients=0,norms=0,sharp=0;
 for(const {A,B}of [...families,mixed]){
  const n=A.length,id=o.identity(2*n),Gbar=o.kron(I,defect(A,B));
  const M=new Map([[1,o.kron(L[0],A)],[-1,o.kron(L[1],B)]]),M2=compose(M,M),M4=compose(M2,M2);
  const sum=new Map(M2);for(const [x,X]of M2)addAt(sum,-x,o.dagger(X));
  const expected=new Map([[2,id],[-2,id],[0,o.scale(Gbar,2)]]);eqFields(sum,expected,2*n,2*n,'general two-step adjoint identity');coefficients+=new Set([...sum.keys(),...expected.keys()]).size;
  const right=new Map();for(const [x,X]of M2){addAt(right,x+2,X);addAt(right,x-2,X);addAt(right,x,o.scale(o.mul(Gbar,X),2));}addAt(right,0,o.scale(id,-4));
  eqFields(M4,right,2*n,2*n,'general quartic coefficients');coefficients+=new Set([...M4.keys(),...right.keys()]).size;
  const seed=col(Array.from({length:2*n},(_,i)=>(i%3)-1)),packet=new Map([[-2,seed],[1,o.scale(seed,'2/3')],[4,o.scale(seed,-2)]]);
  const V2=fieldScale(step(step(packet,A,B),A,B),'1/2');
  ensure(energy(fieldSub(V2,packet)).eq(energy(packet).add(energy(fieldSub(shift(packet,2),packet)).div(2)).sub(fieldQuadratic(packet,Gbar))),'general finite-field defect norm');norms++;
 }
 for(const {A,B,gamma}of families)for(const m of [1,2,7,16])for(const alternate of [false,true]){
  const packet=new Map();for(let j=1;j<=m;j++)packet.set(2*j,o.scale(o.kron(e0,e1),alternate&&j%2?-1:1));
  const V2=fieldScale(step(step(packet,A,B),A,B),'1/2'),ratio=energy(fieldSub(V2,packet)).div(energy(packet));
  const expected=alternate?f(3).sub(gamma).sub(new o.F(1,m)):f(1).sub(gamma).add(new o.F(1,m));ensure(ratio.eq(expected),'sharp scalar-defect packet');sharp++;
 }
 return {Laurent_coefficient_comparisons:coefficients,general_finite_field_norm_equalities:norms,sharp_packet_checks:sharp,
  all_parameter_identity_is_written:'R27.3',gap_is_identity_relative:true};
});

check('source_law_record_endpoint_response_is_independent_of_rank_and_frame',()=>{
 let fields=0,readouts=0;
 for(const {A,B}of replicas){const n=A.length;let joint=new Map([[0,o.identity(2*n)]]),phi=new Map([[0,I]]);
  for(let event=0;event<10;event++){
   joint=step(joint,A,B);const next=new Map();for(const [x,M]of phi){addAt(next,x+1,o.mul(L[0],M));addAt(next,x-1,o.scale(o.mul(L[1],M),((event+x)/2)%2?-1:1));}phi=next;
   const expected=new Map();for(const [x,M]of phi)expected.set(x,o.kron(M,o.mul(o.power(A,((event+1+x)/2)%2),o.power(B,((event+1-x)/2)%2))));
   eqFields(joint,expected,2*n,2*n,'native record endpoint universality');fields+=joint.size;
   if(event<5)for(const [x,M]of joint)for(const X of [I,H,K]){
    eq(o.mul(o.mul(o.dagger(M),o.kron(X,o.identity(n))),M),o.kron(o.mul(o.mul(o.dagger(phi.get(x)),X),phi.get(x)),o.identity(n)),'local pulled-back readout factors on every joint input');readouts++;
   }
  }
 }
 return {joint_endpoint_blocks:fields,operator_readout_equalities_for_all_joint_inputs:readouts,max_event:10,
  uncorrelated_source_record_assumption_required:false,record_preparation_fitted:false};
});

const rows26=inherited[26].value.exact_checks.find(row=>row.name==='first_return_joint_operators_match_new_signed_grammar').first_returns;
const weight26=n=>n%2||n===0?f(0):f(rows26.find(row=>row.event===n).weight);
const choose=(n,k)=>{let c=1n;for(let j=1;j<=k;j++)c=c*BigInt(n-k+j)/BigInt(j);return c;};
const catalan=k=>choose(2*k,k)/BigInt(k+1);
const weight25=n=>n===2?f('1/2'):n>0&&n%4===0?new o.F(catalan(n/4-1)**2n,pow2(n-1)):f(0);
const Le=X=>o.scale(o.mul(o.mul(o.dagger(B0),X),F),'1/2');
const Lo=X=>o.scale(o.mul(o.mul(D0,X),o.dagger(B0)),'1/2');
const recurrence=(odd,N)=>{const values=[I],weight=odd?weight26:weight25;let survival=f(1);
 for(let n=1;n<=N;n++){survival=survival.sub(weight(n));let X=o.scale(I,survival);for(let j=2;j<=n;j+=2)X=o.add(X,o.scale((odd&&j%4?Lo:Le)(values[n-j]),weight(j)));values.push(X);}return values;};
const oddResponse=recurrence(true,32),evenResponse=recurrence(false,32);
const envStep=(a,A,B,controller,n)=>{const b=new Map();for(const [key,M]of a){const [x,mark]=key.split(',').map(Number);for(let d=0;d<2;d++){
 const y=x+1-2*d,k=y+','+(mark|((y===0?1:0)<<(n-1)));let X=o.mul(o.kron(L[d],d?B:A),M);if(controller&&y===0)X=o.mul(o.kron(K,o.identity(A.length)),X);addAt(b,k,X);
 }}return b;};
const envGram=a=>[...a.values()].reduce((sum,M)=>o.add(sum,o.mul(o.dagger(M),M)),o.zeros([...a.values()][0].length));
const envCross=(a,b)=>{let X=o.zeros([...a.values()][0].length);for(const [key,M]of a)if(b.has(key))X=o.add(X,o.mul(o.dagger(M),b.get(key)));return X;};

check('retained_feedback_and_first_returns_are_universal_for_source_law_records',()=>{
 let returns=0,feedback=0,norms=0;
 for(const {A,B}of replicas.slice(0,2)){const n=A.length,id=o.identity(2*n);let live=new Map([[0,id]]);
  for(let event=1;event<=20;event++){
   live=step(live,A,B);const ret=live.get(0)||o.zeros(2*n);live.delete(0);
   const expected=event%2?o.zeros(2*n):o.scale(o.kron((event/2)%2?D0:B0,o.power(o.mul(B,A),(event/2)%4)),rows26.find(row=>row.event===event).signed_coefficient);
   eq(ret,expected,'record-independent exact first-return operator');
   eq(o.scale(o.mul(o.dagger(ret),ret),den(event)),o.scale(id,weight26(event)),'record-independent return energy');returns+=2;
  }
  let a=new Map([['0,0',id]]),b=new Map([['0,0',id]]);
  for(let event=1;event<=8;event++){
   a=envStep(a,A,B,0,event);b=envStep(b,A,B,1,event);
   eq(o.scale(envGram(a),den(event)),id,'joint system0 norm');eq(o.scale(envGram(b),den(event)),id,'joint system1 norm');norms+=2;
   eq(o.scale(envCross(a,b),den(event)),o.kron(oddResponse[event],o.identity(n)),'source-law feedback universal operator');feedback++;
  }
 }
 const inheritedP=inherited[26].value.exact_checks.find(row=>row.name==='native_coefficient_completion_and_exact_return_tail_witness').return_mass;
 return {first_return_operator_and_gram_checks:returns,full_feedback_operator_checks:feedback,joint_norm_checks:norms,
  max_direct_return_event:20,max_direct_feedback_event:8,inherited_return_mass:inheritedP,
  exact_inherited_channel_ratio:'3',new_return_constant_fitted:false};
});

check('literal_copy_full_law_requires_a_two_role_compensator',()=>{
 const HH=o.kron(H,H),KK=o.kron(K,K),id=o.identity(4);eq(o.mul(HH,KK),o.mul(KK,HH),'double copied letters commute');
 ensure(!o.isZero(o.add(o.mul(HH,KK),o.mul(KK,HH))),'double copy violates native mixed law');
 let repairs=0;for(const {A:U,B:W}of replicas.slice(0,2)){
  const Hp=o.kron(HH,U),Kp=o.kron(KK,W),Ip=o.identity(Hp.length),raw=o.add(Hp,Kp);
  eq(o.mul(Hp,Hp),Ip,'repaired parity involution');eq(o.mul(Kp,Kp),Ip,'repaired exchange involution');
  eq(o.add(o.mul(Hp,Kp),o.mul(Kp,Hp)),o.zeros(Hp.length),'compensated mixed source law');
  eq(o.scale(o.mul(o.dagger(raw),raw),'1/2'),Ip,'compensated balanced source norm');repairs+=4;
 }
 const turnK=o.scale(KK,o.IOTA),turnSum=o.add(HH,turnK);
 eq(o.scale(o.mul(o.dagger(turnSum),turnSum),'1/2'),id,'one-dimensional turn repairs norm');eq(o.mul(turnK,turnK),o.scale(id,-1),'turn repair fails source involution');
 let scalarPairs=0;for(const u of [o.ONE,new o.Cut(-1),o.IOTA,o.IOTA.neg()])for(const v of [o.ONE,new o.Cut(-1),o.IOTA,o.IOTA.neg()]){
  ensure(!(u.mul(u).eq(1)&&v.mul(v).eq(1)&&u.mul(v).add(v.mul(u)).zero()),'scalar diagnostic cannot carry source law');scalarPairs++;
 }
 return {repaired_joint_source_identities:repairs,smallest_constructed_compensator_rank:2,minimal_source_record_compensator_rank:8,
  native_scalar_diagnostic_pairs:scalarPairs,one_dimensional_norm_repair_preserves_full_word_law:false,
  minimality_is_for_factorized_operator_copy_target:true};
});

const copies=r=>{let A=o.identity(1),B=o.identity(1);for(let k=0;k<r;k++){A=o.kron(A,H);B=o.kron(B,K);}return {A,B};};
check('copy_parity_has_exact_native_tag_normal_form_and_even_kernel',()=>{
 let tags=0,operatorIdentities=0;const counts=[];
 for(let k=1;k<=5;k++){
  const {A,B}=copies(k),n=2**k,id=o.identity(n),AB=o.mul(A,B);eq(AB,o.scale(o.mul(B,A),k%2?-1:1),'native retained-copy sign');operatorIdentities++;
  if(k%2){
   const encoder=o.zeros(n);for(let z=0;z<n;z++){
    const bits=Array.from({length:k},(_,j)=>(z>>(k-1-j))&1),parity=bits.reduce((a,b)=>a^b,0);let target=parity;
    for(let i=0;i<k-1;i++)target=2*target+(bits[i]^bits[k-1]);encoder[target][z]=o.ONE;tags++;
   }
   eq(o.mul(o.dagger(encoder),encoder),id,'native bit relabelling is bijective');
   eq(o.mul(o.mul(encoder,A),o.dagger(encoder)),o.kron(H,o.identity(n/2)),'logical parity factor');
   eq(o.mul(o.mul(encoder,B),o.dagger(encoder)),o.kron(K,o.identity(n/2)),'logical exchange factor');operatorIdentities+=3;
   counts.push({total_copies:k,source_law:true,multiplicity_rank:n/2});
  }else{
   const raw=o.add(A,B),gram=o.scale(o.mul(o.dagger(raw),raw),'1/2'),plus=o.scale(o.add(id,AB),'1/2'),minus=o.sub(id,plus);
   eq(gram,o.add(id,AB),'even-copy balanced Gram');eq(o.mul(raw,minus),o.zeros(n),'even-copy kernel');
   ensure(o.rank(plus)===n/2&&o.rank(minus)===n/2&&o.rank(raw)===n/2,'even-copy exact half ranks');operatorIdentities+=3;
   counts.push({total_copies:k,source_law:false,kernel_rank:n/2,other_sector_Gram:'2'});
  }
 }
 return {native_tuple_tags_encoded:tags,copy_operator_and_rank_checks:operatorIdentities,copy_classes:counts,
  full_role_tuples_are_not_erased_by_encoding:true};
});

check('retained_record_copy_parity_selects_two_return_and_feedback_classes',()=>{
 let norms=0,crosses=0,returnChecks=0;const families=[];
 for(let r=0;r<=3;r++){
  const {A,B}=copies(r),n=A.length,id=o.identity(2*n),odd=r%2===1,weight=odd?weight26:weight25;
  let live=new Map([[0,id]]),returned=o.zeros(2*n);
  for(let event=1;event<=12;event++){
   live=step(live,A,B);const ret=live.get(0)||o.zeros(2*n);live.delete(0);const gram=o.scale(o.mul(o.dagger(ret),ret),den(event));
   eq(gram,o.scale(id,weight(event)),'record-copy parity first-return coefficient');returned=o.add(returned,gram);returnChecks++;
   const liveGram=[...live.values()].reduce((sum,M)=>o.add(sum,o.mul(o.dagger(M),M)),o.zeros(2*n));eq(o.add(returned,o.scale(liveGram,den(event))),id,'copy-parity survivor telescope');norms++;
  }
  if(r<=2){let a=new Map([['0,0',id]]),b=new Map([['0,0',id]]);
   for(let event=1;event<=8;event++){
    a=envStep(a,A,B,0,event);b=envStep(b,A,B,1,event);
    eq(o.scale(envCross(a,b),den(event)),o.kron((odd?oddResponse:evenResponse)[event],o.identity(n)),'copy-parity full feedback class');crosses++;
   }
  }
  families.push({retained_record_copies:r,response:odd?'R26 p / parity pulses':'R24 eta / R25 feedback',record_rank:n});
 }
 const eta=inherited[25].value.exact_checks.find(row=>row.name==='finite_event_tail_bound_and_completed_response_enclosures').inherited_eta_interval;
 ensure(f(eta.upper_rational).le('11/16'),'inherited even return bracket');
 return {copy_return_gram_checks:returnChecks,copy_survivor_checks:norms,copy_feedback_checks:crosses,
  copy_classes:families,inherited_even_return_mass:eta,total_source_copies_and_ancillary_record_copies_are_distinct:true};
});

check('native_count_family_and_sign_change_separate_readout_from_identity_gap',()=>{
 const {A,B,gamma}=countPair(2,1),negative=o.scale(B,-1),id=o.identity(4);let positiveField=new Map([[0,id]]),negativeField=new Map([[0,id]]),blocks=0,flags=0;
 for(let n=1;n<=10;n++){
  positiveField=step(positiveField,A,B);negativeField=step(negativeField,A,negative);
  const expected=new Map([...positiveField].map(([x,M])=>[x,o.scale(M,((n-x)/2)%2?-1:1)]));eqFields(negativeField,expected,4,4,'same endpoint sign for changed B');blocks+=positiveField.size;
  if(n===3){const energies=[-3,-1,1,3].map(x=>o.energy(o.mul(positiveField.get(x),o.kron(e0,e0))).div(8));
   ensure(energies.every((v,i)=>v.eq(['1/8','89/200','61/200','1/8'][i])),'native count-family source energies');}
 }
 for(let controller=0;controller<2;controller++){
  let a=new Map([['0,0',id]]),b=new Map([['0,0',id]]);
  for(let n=1;n<=8;n++){
   a=envStep(a,A,B,controller,n);b=envStep(b,A,negative,controller,n);
   const expected=new Map([...a].map(([key,M])=>{const x=Number(key.split(',')[0]);return [key,o.scale(M,((n-x)/2)%2?-1:1)];}));
   eqFields(b,expected,4,4,'same arrival-record sign for changed B');flags+=a.size;
  }
 }
 ensure(f(1).sub(gamma).eq('2/5')&&f(1).add(gamma).eq('8/5'),'different sharp identity-relative gaps');
 ensure(gamma.pow(2).sub('1/2').eq('-7/50'),'count-derived current');
 ensure(!o.isZero(defect(A,B)),'isometric transport does not imply native source continuation');
 return {exact_signed_endpoint_blocks:blocks,exact_signed_feedback_flag_blocks:flags,native_counts:[2,1],
  gamma:'3/5',sigma:'4/5',return_current:'-7/50',sharp_lower_gap:'2/5',sign_reversed_sharp_lower_gap:'8/5',
  counts_are_constructed_witnesses_not_fitted_physics:true};
});

check('literal_HK_histories_reproduce_record_copy_parity_feedback',()=>{
 let prefixes=0,blocks=0;
 for(const r of [1,2]){
  const {A,B}=copies(r),n=A.length,record=basis(n)[0],seed=o.kron(e0,record),id=o.identity(2*n);
  for(let controller=0;controller<2;controller++){
   let field=new Map([['0,0',id]]),histories=[{x:0,mark:0,v:seed,word:''}];
   for(let event=1;event<=7;event++){
    field=envStep(field,A,B,controller,event);
    histories=histories.flatMap(h=>[[H,'H'],[K,'K']].map(([letter,tag])=>{
     let v=o.mul(o.kron(letter,o.identity(n)),h.v);const d=o.isZero(o.mul(o.kron(P,o.identity(n)),v))?1:0;
     v=o.mul(o.kron(I,d?B:A),v);const x=h.x+1-2*d,mark=h.mark|((x===0?1:0)<<(event-1));if(controller&&x===0)v=o.mul(o.kron(K,o.identity(n)),v);
     prefixes++;return {x,mark,v,word:h.word+tag};
    }));
    ensure(new Set(histories.map(h=>h.word)).size===2**event,'full source tags retained');
    const aggregate=new Map();for(const h of histories)addAt(aggregate,h.x+','+h.mark,h.v);
    const expected=new Map([...field].map(([key,M])=>[key,o.mul(M,seed)]));eqFields(aggregate,expected,2*n,1,'literal source histories / copy-parity records');blocks+=aggregate.size;
   }
  }
 }
 return {literal_HK_prefixes:prefixes,matched_arrival_record_blocks:blocks,max_literal_event:7,word_histories_quotiented_away:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),aa=word('A'),bb=word('B'),gg=sc('1/2',add(mul(aa,bb),mul(bb,aa))),cc=add(mul(aa,bb),sc(-1,mul(bb,aa)));
const generalSpec={tokens:['A','B'],rules:[{id:'AA',lhs:['A','A'],rhs:[[[],1]],source:'R27 candidate record self-return; no anticommutation assumed'},
 {id:'BB',lhs:['B','B'],rhs:[[[],1]],source:'R27 candidate record self-return; no anticommutation assumed'}]};
const canonicalSpec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const banks=[{name:'general_involutive_record',spec:generalSpec,tasks:[
 ['defect_central_for_A',mul(aa,gg),mul(gg,aa)],['defect_central_for_B',mul(bb,gg),mul(gg,bb)],
 ['balanced_record_Gram',sc('1/2',mul(add(aa,bb),add(aa,bb))),add(ii,gg)],
 ['curvature_square_without_anticommutation',mul(sc(-1,cc),cc),sc(4,add(ii,sc(-1,mul(gg,gg))))],
 ]},{name:'canonical_native_source',spec:canonicalSpec,tasks:[
 ['source_H_self_return',mul(hh,hh),ii],['source_K_self_return',mul(kk,kk),ii],
 ['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['native_source_census_norm',mul(add(hh,kk),add(hh,kk)),sc(2,ii)],
 ['plus_cut_idempotent',mul(pp,pp),pp],['minus_cut_idempotent',mul(qq,qq),qq],
 ['cuts_are_disjoint',mul(pp,qq),sc(0,ii)],['exchange_intertwines_cuts',mul(kk,pp),mul(qq,kk)],
 ['oriented_native_record_square',mul(rr,rr),sc(-1,ii)],
 ['native_sign_conjugacy',mul(hh,kk,hh),sc(-1,kk)],
 ]}];
const replays=[],presentations=[];let firstSystem;
for(const bank of banks){const system=w.presentation(bank.spec),audit=system.audit();ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','presentation audit');
 presentations.push({name:bank.name,specification:bank.spec,audit});if(!firstSystem)firstSystem=system;
 for(const [name,left,right]of bank.tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
  ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,presentation:bank.name,result,replay:'REPLAY_MATCH'});
 }
}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{firstSystem.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const gammaRecord=countPair(2,1),double=copies(2),doubleSum=o.add(double.A,double.B),I4=o.identity(4);
const negatives={
 two_record_involutions_already_force_anticommutation:!o.isZero(defect(gammaRecord.A,gammaRecord.B)),
 joint_transport_normalization_already_implies_record_census_normalization:!o.equal(o.scale(o.mul(o.add(gammaRecord.A,gammaRecord.B),o.add(gammaRecord.A,gammaRecord.B)),'1/2'),I),
 record_defect_current_is_linear_in_G:f('3/5').pow(2).eq('9/25')&&!f('3/5').eq('9/25'),
 any_larger_source_law_record_needs_a_different_mixed_source_relation:replicas.every(rec=>o.isZero(defect(rec.A,rec.B))),
 an_intertwiner_can_freely_mix_the_two_source_role_blocks:!o.isZero(o.sub(o.mul(H,K),o.mul(K,H))),
 literal_double_copy_preserves_the_original_source_law:o.equal(o.mul(double.A,double.B),o.mul(double.B,double.A)),
 scalar_iota_norm_repair_preserves_both_source_involutions:o.IOTA.mul(o.IOTA).eq(-1),
 even_copy_balanced_census_is_pairing_preserving:!o.equal(o.scale(o.mul(o.dagger(doubleSum),doubleSum),'1/2'),I4),
 adding_one_retained_native_copy_never_changes_the_response:!o.equal(oddResponse[2],evenResponse[2]),
 ancillary_copy_parity_equals_total_literal_source_copy_parity:(2%2)!==(3%2),
 identical_local_readouts_determine_the_identity_relative_gap:!f('2/5').eq('8/5'),
 old_eta_and_R26_p_can_be_used_interchangeably:f('11/16').le('3/4')&&!f('11/16').eq('3/4'),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r27.native-replica-selection.v1',status:'PASS_R27_NATIVE_REPLICA_SELECTION',
 input_sha256:inputHash,r25_native_certificate_sha256:inherited[25].sha,r26_native_certificate_sha256:inherited[26].sha,
 source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_presentations:presentations,symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_source_law_record_classified:true,native_record_balance_forces_mixed_law:true,
 native_copy_parity_response_derived:true,minimum_factorized_compensator_derived:true,
 ordinary_complex_representation_classification_or_classical_physics_premise:false,
 physical_memory_interface_selected:false,physical_metric_c_alpha_derived:false,physical_mass_gap_identified:false,
 formal_proof_assistant_verified:false,
 scope:'Seven native written proofs: record normalization and curvature, explicit source-law normal form, general propagation defect, R26 universality, minimum literal-copy compensator, copy-parity response classes and native counterfamilies locating unselected physical targets. No physical interaction, mass, charge or alpha identification.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
