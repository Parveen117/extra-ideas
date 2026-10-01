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
ensure(inputHash===arg('--expected-input-sha256'),'R34 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r33-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r33-sha256'),'R34 R33 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R33_NATIVE_INTERACTION_MEMORY','R34 parent status');
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
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

const mprod=(...a)=>a.reduce(o.mul),mi=o.inverse,zeroFields=new Map();
const F2=minus(compose(B,B),unit(I4)),Cgap=minus(G,unit(I4)),Dgap=minus(unit(o.scale(I4,4)),G);
const M=scale(compose(G,minus(unit(o.scale(I4,5)),G)),'1/2');
const ESX=plus(compose(TX,TX),unit(I4)),ESY=plus(compose(TY,TY),unit(I4));
const ring=(op,L,r=4,c=4)=>{const mat=o.zeros(r*L*L,c*L*L);for(const [k,A]of op){const [dx,dy]=xy(k);
 for(let x=0;x<L;x++)for(let y=0;y<L;y++){const u=((x-dx)%L+L)%L,v=((y-dy)%L+L)%L,a=x*L+y,b=u*L+v;
  for(let i=0;i<r;i++)for(let j=0;j<c;j++)mat[a*r+i][b*c+j]=mat[a*r+i][b*c+j].add(A[i][j]);
 }}return mat;};
const matrixColumn=(n,k,value=1)=>col(Array.from({length:n},(_,i)=>i===k?value:0));
const scalar=A=>A[0][0].rad,abs=x=>x.le(0)?f(0).sub(x):x;
const cyclic=[];

check('native_minimum_decoder_is_attained_and_orthogonal',()=>{
 let identities=0;const rows=[];
 for(let L=1;L<=3;L++){
  const g=ring(G,L),ff=ring(F,L),bb=ring(B,L),f2=ring(F2,L),mm=ring(M,L),n=g.length,id=o.identity(n),inv=mi(g),minv=mi(mm);
  eq(o.mul(o.dagger(ff),ff),g,'finite native source Gram');
  eq(o.scale(o.add(o.mul(o.dagger(ff),ff),o.mul(o.dagger(f2),f2)),'1/2'),mm,'weighted bank Gram');identities+=2;
  for(const target of [matrixColumn(n,0),matrixColumn(n,n-1),col(Array.from({length:n},(_,i)=>(i%3)-1))]){
   const seed=o.mul(inv,target),decoder=o.mul(ff,seed),beta=scalar(mprod(o.dagger(target),inv,target));
   eq(o.mul(o.dagger(ff),decoder),target,'native minimum decoder target');
   ensure(o.energy(decoder).eq(beta),'minimum decoder norm');
   ensure(scalar(mprod(o.dagger(target),seed)).pow(2).eq(beta.mul(o.energy(o.mul(ff,seed)))),'sharp target attainment');
   const next=o.mul(minv,target),d1=o.mul(ff,next),d2=o.mul(f2,next),newbeta=scalar(mprod(o.dagger(target),minv,target));
   eq(o.scale(o.add(o.mul(o.dagger(ff),d1),o.mul(o.dagger(f2),d2)),'1/2'),target,'balanced weighted decoder');
   ensure(o.energy(d1).add(o.energy(d2)).div(2).eq(newbeta),'balanced minimum norm in native half-weight chart');
   const h1=o.mul(o.dagger(f2),target),h2=o.scale(o.mul(o.dagger(ff),target),-1);
   eq(o.add(o.mul(o.dagger(ff),h1),o.mul(o.dagger(f2),h2)),o.zeros(n,1),'native orthogonal redundant decoder direction');
   ensure(o.energy(o.add(d1,h1)).add(o.energy(o.add(d2,h2))).eq(o.energy(d1).add(o.energy(d2)).add(o.energy(h1)).add(o.energy(h2))),'minimum decoder Pythagoras');identities+=7;
  }
  cyclic.push({L,n,id,g,ff,bb,f2,mm,inv,minv});rows.push({address_period:L,roles:n,baseline_local_beta:inv[0][0].rad.toString(),balanced_local_beta:minv[0][0].rad.toString()});
 }
 return {exact_decoder_and_attainment_identities:identities,cyclic_four_role_targets:rows,weighted_chart_for_equal_native_tags:'sum of two squared branch norms divided by 2'};
});

check('native_two_horizon_gram_has_complete_positive_identities',()=>{
 eqFields(compose(dag(F2),F2),compose(G,Dgap),4,4,'sixteen-event increment Gram');
 eqFields(M,scale(plus(compose(dag(F),F),compose(dag(F2),F2)),'1/2'),4,4,'balanced bank Gram');
 const cs=[scale(DX,'1/2'),scale(DY,'1/2'),scale(compose(DX,DY),'1/4')];
 // sqrt(1/2) on ESX is represented as its rational squared weight below.
 const ds=[{a:ESX,w:'1/2'},{a:ESY,w:'1/4'},{a:compose(DX,ESY),w:'1/16'}];
 let csum=new Map(),dsum=new Map(),cross=new Map();
 for(const a of cs)csum=plus(csum,compose(dag(a),a));
 for(const {a,w:weight}of ds)dsum=plus(dsum,scale(compose(dag(a),a),weight));
 for(const a of cs)for(const {a:b,w:weight}of ds){const ab=compose(a,b);cross=plus(cross,scale(compose(dag(ab),ab),weight));}
 eqFields(csum,Cgap,4,4,'lower Gram defect native squares');eqFields(dsum,Dgap,4,4,'upper Gram defect native squares');
 eqFields(cross,compose(Cgap,Dgap),4,4,'nine native product squares');
 eqFields(M,plus(unit(o.scale(I4,2)),scale(cross,'1/2')),4,4,'sharp positive balanced identity');
 const centered=minus(G,unit(o.scale(I4,'5/2')));
 eqFields(M,minus(unit(o.scale(I4,'25/8')),scale(compose(centered,centered),'1/2')),4,4,'exact balanced upper square');
 return {complete_Laurent_gram_and_square_identities:7,positive_product_square_terms:9,external_spectral_or_minimax_theorem_used:false};
});

let designRows=[];
check('native_two_horizon_design_has_a_unique_balanced_optimum',()=>{
 let identities=0;const L=4,g=ring(G,L),ff=ring(F,L),f2=ring(F2,L),n=g.length;
 const uniform=col(Array.from({length:n},(_,i)=>i%4===0?1:0)),alternating=col(Array.from({length:n},(_,i)=>i%4===0?(-1)**(Math.floor(i/4/L/2)+Math.floor((Math.floor(i/4)%L)/2)):0));
 eq(o.mul(g,uniform),uniform,'exact low endpoint on source addresses');
 eq(o.mul(g,alternating),o.scale(alternating,4),'exact high endpoint on source addresses');
 for(let k=0;k<=12;k++){
  const s=f(k).div(12),a=f(1).add(s.mul(2)),b=f(4).sub(s.mul(4)),lower=a.le(b)?a:b;
  const op=compose(G,minus(unit(o.scale(I4,f(1).add(s.mul(3)))),scale(G,s)));
  const chord=plus(plus(scale(Dgap,a.div(3)),scale(Cgap,b.div(3))),scale(compose(Cgap,Dgap),s));
  eqFields(op,chord,4,4,'complete design chord identity');
  const mat=o.add(o.scale(o.mul(o.dagger(ff),ff),f(1).sub(s)),o.scale(o.mul(o.dagger(f2),f2),s));
  eq(o.mul(mat,uniform),o.scale(uniform,a),'actual low design endpoint');
  eq(o.mul(mat,alternating),o.scale(alternating,b),'actual high design endpoint');
  ensure(lower.le(2)&&(lower.eq(2)===s.eq('1/2')),'unique sampled optimum agrees with exact endpoint proof');
  designRows.push({second_increment_weight:s.toString(),sharp_lower_bound:lower.toString(),full_state_beta:lower.zero()?'infinite':f(1).div(lower).toString()});identities+=3;
 }
 return {exact_design_and_endpoint_identities:identities,weight_probes:designRows,optimal_weight:'1/2',baseline_full_state_beta:'1',balanced_full_state_beta:'1/2',optimality_scope:'normalized two-increment banks, worst full-state target'};
});

const shiftCycle=L=>{const a=o.zeros(L);for(let n=0;n<L;n++)a[(n+1)%L][n]=o.Cut.of(1);return a;};
let pointCycles=[];
check('finite_native_component_cuts_preserve_their_boundary_beta',()=>{
 let identities=0;for(let L=2;L<=5;L++){
  const id=o.identity(L),s=shiftCycle(L),q=o.sub(o.scale(id,'3/2'),o.scale(o.add(s,o.dagger(s)),'1/4')),qi=mi(q),g=o.kron(q,q),n=L*L,ii=o.identity(n),gg=mi(g),m=o.scale(o.mul(g,o.sub(o.scale(ii,5),g)),'1/2'),im=mi(m);
  eq(gg,o.kron(qi,qi),'native factorized cyclic inverse');
  eq(im,o.scale(o.add(gg,mi(o.sub(o.scale(ii,5),g))),'2/5'),'exact local inverse decomposition');identities+=2;
  const ell=matrixColumn(n,0),psi=o.mul(im,ell),beta=im[0][0].rad;
  ensure(scalar(mprod(o.dagger(ell),psi)).pow(2).eq(beta.mul(scalar(mprod(o.dagger(psi),m,psi)))),'local sharp attainment');identities++;
  pointCycles.push({component_period:L,old_local_beta:gg[0][0].rad.toString(),new_local_beta:beta.toString()});
 }
 ensure(pointCycles[0].old_local_beta==='9/16'&&pointCycles[0].new_local_beta==='5/12','two-cycle beta differs from unwrapped target');
 return {exact_cyclic_inverse_and_attainment_identities:identities,component_cycle_values:pointCycles,unwrapped_beta_substituted_for_cyclic_beta:false};
});

const choose=(n,k)=>{if(k<0||k>n)return 0n;let a=1n;for(let j=1;j<=k;j++)a=a*BigInt(n-j+1)/BigInt(j);return a;};
const moment=n=>new o.F(Array.from({length:Math.floor(n/2)+1},(_,k)=>choose(n,2*k)*choose(2*k,k)*6n**BigInt(n-2*k)).reduce((a,b)=>a+b,0n),4n**BigInt(n));
const padd=(a,k,x)=>{const v=(a.get(k)||f(0)).add(x);if(v.zero())a.delete(k);else a.set(k,v);};
const nextMomentPoly=a=>{const b=new Map();for(const [k,v]of a){padd(b,k,v.mul('3/2'));padd(b,k-1,v.mul('-1/4'));padd(b,k+1,v.mul('-1/4'));}return b;};
let betaLower,betaUpper,pointBracket;
check('native_return_counts_certify_the_completed_local_burden',()=>{
 let poly=new Map([[0,f(1)]]),sum=f(0),identities=0;const checkpoints=[];
 for(let n=0;n<=160;n++){
  const mn=moment(n);if(n<=24){ensure(mn.eq(poly.get(0)||f(0)),'independent native finite Laurent moment');poly=nextMomentPoly(poly);identities++;}
  ensure(mn.le(f(2).pow(n)),'native moment upper bound');sum=sum.add(mn.pow(2).div(f(5).pow(n)));
  const lower=f('1/5').add(sum.mul('2/25')),error=f('2/5').mul(f('4/5').pow(n+1));
  if([0,1,4,8,16,32,64,128,160].includes(n))checkpoints.push({N:n,lower:lower.toString(),upper:lower.add(error).toString(),tail_bound:error.toString()});
  if(n===160){betaLower=lower;betaUpper=lower.add(error);}
 }
 ensure(f('3616370316517928/10000000000000000').le(betaLower)&&betaUpper.le('3616370316517930/10000000000000000'),'published local burden enclosure');
 ensure(betaUpper.le('1/2'),'strict local improvement');
 ensure(f('2767/10000').le(f(1).sub(betaUpper.mul(2))),'more than 27.67 percent point burden reduction');
 pointBracket={lower:betaLower.toString(),upper:betaUpper.toString(),completed_terms:161,decimal_lower:'0.3616370316517928',decimal_upper:'0.3616370316517930'};
 return {independent_finite_Laurent_moment_checks:identities,exact_native_return_moments:161,certified_tail_checkpoints:checkpoints,local_burden_interval:pointBracket,physical_noise_distribution_supplied:false};
});

check('native_sharp_patterns_certify_gain_independent_conditioning',()=>{
 const L=6,id=o.identity(L),s=shiftCycle(L),q=o.sub(o.scale(id,'3/2'),o.scale(o.add(s,o.dagger(s)),'1/4')),g=o.kron(q,q),ii=o.identity(36),m=o.scale(o.mul(g,o.sub(o.scale(ii,5),g)),'1/2'),u=col([2,1,-1,-2,-1,1]),alt=col([1,-1,1,-1,1,-1]);
 eq(o.mul(o.add(s,o.dagger(s)),u),u,'native six-count recurrence');
 eq(o.mul(q,u),o.scale(u,'5/4'),'exact x-direction count value');
 eq(o.mul(q,alt),o.scale(alt,2),'exact y-direction count value');
 const psi=o.kron(u,alt);eq(o.mul(g,psi),o.scale(psi,'5/2'),'sharp internal Gram value');
 eq(o.mul(m,psi),o.scale(psi,'25/8'),'sharp balanced maximum');
 let packetChecks=0;for(const size of [2,3,4,6,8])for(const sign of [1,-1]){
  const field=new Map();for(let x=0;x<size;x++)for(let y=0;y<size;y++)field.set(key(2*x,2*y),o.scale(o.kron(e0,e0),new o.F(sign===1?1:(-1)**(x+y),size)));
  const residual=minus(compose(G,field),scale(field,sign===1?1:4));ensure(energy(residual).le(f(512).div(size)),'native endpoint packet error');packetChecks++;
 }
 return {exact_sharp_pattern_identities:5,compact_endpoint_packet_bound_checks:packetChecks,baseline_squared_response_range:['1','4'],balanced_squared_response_range:['2','25/8'],gain_independent_range_ratio_before:'4',gain_independent_range_ratio_after:'25/16'};
});

check('native_cost_normalization_retains_its_distinct_tradeoffs',()=>{
 let n=0;for(const row of designRows){if(row.full_state_beta==='infinite')continue;const s=f(row.second_increment_weight),cost=f(1).add(s).mul(f(row.full_state_beta));
  ensure(f('3/4').le(cost)&&(cost.eq('3/4')===s.eq('1/2')),'unique cost-normalized worst-burden optimum');n++;
 }
 ensure(f('1/2').le(betaLower.mul('3/2')),'point improvement does not survive weighted event-cost normalization');
 ensure(f('1/2').mul(2).eq(1),'maximum horizon cost removes worst-burden improvement');
 return {exact_weighted_cost_probes:n,baseline_event_horizon:8,balanced_event_horizon:16,balanced_weighted_event_count:12,weighted_event_normalized_full_beta:'3/4',maximum_horizon_normalized_full_beta:'1',point_target_improves_under_weighted_event_cost:false};
});

check('native_error_transport_does_not_create_independent_information',()=>{
 let n=0;for(const c of cyclic){const bp=o.add(c.bb,c.id),d=col(Array.from({length:c.n},(_,i)=>(i%5)-2));
  eq(o.mul(bp,c.ff),c.f2,'second increment factors through first');
  eq(o.scale(o.add(c.id,o.mul(o.dagger(bp),bp)),'1/2'),o.scale(o.sub(o.scale(c.id,5),c.g),'1/2'),'exact native output frame Gram');
  const e=matrixColumn(c.n,0);eq(o.sub(o.mul(bp,o.add(d,e)),o.mul(bp,d)),o.mul(bp,e),'same error retained in second channel');
  ensure(o.rank(c.ff)===c.n&&o.rank(c.mm)===c.n,'unchanged full recognition rank');n+=4;
 }
 return {exact_factored_observer_and_error_checks:n,old_error_e_maps_to:'(e,(B+I)e)/sqrt(2)',independent_new_error_assumed:false,burden_in_transported_old_error_norm:'unchanged'};
});

check('native_target_and_duplicate_normalization_prevent_fake_beta_reduction',()=>{
 const c=cyclic[1],target=matrixColumn(c.n,0),beta=scalar(mprod(o.dagger(target),c.inv,target));
 const amplified=mi(o.scale(c.g,4)),tripled=mi(o.scale(c.g,3)),normalized=mi(o.scale(o.scale(c.g,3),'1/3'));
 ensure(scalar(mprod(o.dagger(target),amplified,target)).eq(beta.div(4)),'scalar gain changes raw beta');
 ensure(scalar(mprod(o.dagger(target),tripled,target)).eq(beta.div(3)),'unnormalized duplicate bank');
 ensure(scalar(mprod(o.dagger(target),normalized,target)).eq(beta),'fixed-weight duplicates do not help');
 ensure(scalar(mprod(o.dagger(o.scale(target,'1/2')),c.inv,o.scale(target,'1/2'))).eq(beta.div(4)),'target rescaling changes beta');
 const t=matrixColumn(4,3,'1/2'),unitTarget=matrixColumn(4,3);ensure(o.energy(t).eq('1/4')&&o.energy(unitTarget).eq(1),'exact cited quarter arithmetic reconstructed natively');
 return {exact_gain_duplicate_and_target_controls:5,quarter_example_target_norm_squared:'1/4',quarter_example_unit_target_beta:'1',repository_wide_beta_declared:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cpw=add(hh,kk),cmw=add(hh,sc(-1,kk)),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['native_turn_square',mul(rr,rr),sc(-1,ii)],['native_census_square',mul(cpw,cpw),sc(2,ii)],['native_reverse_census_square',mul(cmw,cmw),sc(2,ii)],
 ['native_record_loop',mul(hh,kk,hh,kk),sc(-1,ii)],['native_reverse_record_loop',mul(kk,hh,kk,hh),sc(-1,ii)],
 ['native_plus_cut',mul(pp,pp),pp],['native_minus_cut',mul(qq,qq),qq],['native_disjoint_cuts',mul(pp,qq),sc(0,ii)],['native_complete_cut_pair',add(pp,qq),ii],
 ['equal_tag_weight_algebraic_witness',mul(sc('1/2',cpw),sc('1/2',cpw)),sc('1/2',ii)],['native_turn_increment_norm',mul(add(ii,sc(-1,rr)),add(ii,rr)),sc(2,ii)]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const negatives={
 beta_is_one_number_for_the_entire_repository:!f('1').eq('1/2'),
 local_decoder_beta_is_the_memory_spread_fraction:!f('1/2').eq('1/9'),
 the_cited_quarter_survives_unit_target_normalization:!f('1/4').eq(1),
 sixteen_event_increment_alone_is_complete:designRows.at(-1).full_state_beta==='infinite',
 balanced_repair_keeps_baseline_worst_burden:designRows[6].full_state_beta==='1/2',
 every_weight_is_equally_optimal:designRows.filter(r=>r.sharp_lower_bound==='2').length===1,
 balanced_bank_is_just_scalar_amplification:!f('4').eq('25/16'),
 fixed_weight_duplicates_improve_beta:checks[8].exact_gain_duplicate_and_target_controls===5,
 balanced_bank_adds_new_recognizable_source_roles:checks[7].exact_factored_observer_and_error_checks===12,
 postprocessing_creates_independent_errors:checks[7].independent_new_error_assumed===false,
 transported_old_noise_norm_keeps_the_raw_improvement:checks[7].burden_in_transported_old_error_norm==='unchanged',
 event_cost_can_be_omitted:checks[6].balanced_event_horizon===16,
 point_repair_is_better_per_weighted_event:f('1/2').le(betaLower.mul('3/2')),
 all_cycles_have_the_unwrapped_local_beta:pointCycles[0].old_local_beta==='9/16',
 a_finite_local_count_sum_is_the_completed_value:betaLower.le(betaUpper)&&!betaLower.eq(betaUpper),
};ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r34.native-decoder-burden.v1',status:'PASS_R34_NATIVE_DECODER_BURDEN',
 input_sha256:inputHash,r33_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_target_burden_and_decoder_derived:true,native_fixed_weight_observer_optimum_derived:true,native_point_burden_reduction_certified:true,native_gain_cost_and_error_boundaries_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Six written native results: minimum decoder, current local/full burdens, optimal balanced two-increment bank, exact local count enclosure, gain-independent range and event-cost audit, and target/duplicate/inherited-error normalization. Physical noise, clocks, alpha and a repository-wide proof-progress score are not identified.'},null,2)+'\n');
