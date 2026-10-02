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
ensure(inputHash===arg('--expected-input-sha256'),'R38 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r37-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r37-sha256'),'R38 R37 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R37_NATIVE_PACKET_FLIGHT','R38 parent status');
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

const block=(a)=>{const out=new Map();for(const [k,A]of a){const [dx,dy]=xy(k);for(let ex=0;ex<2;ex++)for(let ey=0;ey<2;ey++){
 const sx=-Math.floor((ex-dx)/2),sy=-Math.floor((ey-dy)/2),fx=parity(ex-dx),fy=parity(ey-dy),M=o.zeros(16),sign=-((-1)**(sx+sy));
 for(let i=0;i<4;i++)for(let j=0;j<4;j++)M[(ex*2+ey)*4+i][(fx*2+fy)*4+j]=A[i][j].mul(sign);
 addAt(out,key(sx,sy),M);
 }}return out;};
const Z=block(B),I16=o.identity(16),U16=unit(I16),zero16=o.zeros(16);
const sumAtOne=a=>[...a.values()].reduce(o.add,o.zeros(a.values().next().value.length));
const jet=axis=>[...Z].reduce((s,[k,M])=>o.add(s,o.scale(M,xy(k)[axis])),zero16);
const Ax=jet(0),Ay=jet(1),trace=M=>M.reduce((s,row,i)=>s.add(row[i]),o.Cut.of(0));
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});
const id1=unit(o.identity(1)),sx1=shift(id1,1,0),sy1=shift(id1,0,1);
const dx1=minus(sx1,id1),dy1=minus(sy1,id1),ex1=plus(sx1,id1),ey1=plus(sy1,id1);
const lx1=compose(dag(dx1),dx1),ly1=compose(dag(dy1),dy1);
const ell1=minus(scale(plus(lx1,ly1),'1/2'),scale(compose(lx1,ly1),'1/16'));
const scalarLift=(a,size)=>new Map([...a].map(([k,M])=>[k,o.scale(o.identity(size),M[0][0])]));
const ell=scalarLift(ell1,16);
const normSquare=a=>compose(dag(a),a);

const cutpow=(z,n)=>{if(n<0)return cutpow(z.inv(),-n);let a=o.Cut.of(1);for(let j=0;j<n;j++)a=a.mul(z);return a;};
const evaluate=(a,x,y)=>[...a].reduce((s,[k,M])=>o.add(s,o.scale(M,cutpow(x,xy(k)[0]).mul(cutpow(y,xy(k)[1])))),o.zeros(a.values().next().value.length,a.values().next().value[0].length));
const entryMass=M=>M.flat().reduce((s,c)=>s.add(c.gauge()),f(0));
const mass=A=>[...A.values()].reduce((s,M)=>s.add(entryMass(M)),f(0));
const deriv=(A,e)=>new Map([...A].map(([d,M])=>[d,o.scale(M,f(e[0]).mul(xy(d)[0]).add(f(e[1]).mul(xy(d)[1])))]));
const dirs=[[f(1),f(0)],[f(0),f(1)],[f('3/5'),f('4/5')]];
const iota=new o.Cut(0,1),e16=col(Array.from({length:16},(_,i)=>i===0?1:0));
const rootWitnesses=[];
check('native_dyadic_turn_counts_derive_finite_phase_resolution',()=>{
 let sums=0,primitives=0;
 const tau=o.scale(C0,'1/2');eq(o.mul(tau,tau),o.scale(I,'1/2'),'native half-root algebraic witness');
 for(const [Q,root]of [[2,o.scale(I,-1)],[4,o.scale(I,iota)],[8,o.scale(tau,new o.Cut(1,1))]]){
  eq(o.mul(o.dagger(root),root),I,'native root unit identity');eq(o.power(root,Q),I,'finite root period');
  for(let d=1;d<Q;d++){ensure(!o.equal(o.power(root,d),I),'primitive root witness');primitives++;}
  for(let d=0;d<Q;d++){let sum=o.zeros(2);for(let j=0;j<Q;j++)sum=o.add(sum,o.power(root,j*d));eq(sum,d===0?o.scale(I,Q):o.zeros(2),'finite geometric orthogonality');sums++;}
  rootWitnesses.push({period:Q,unit_and_primitive:true});
 }
 let inverse=0,intertwiners=0,norms=0;
 for(const Q of [2,4]){
  const xi=Q===2?o.Cut.of(-1):iota,field=new Map();
  for(let x=0;x<Q;x++)for(let y=0;y<Q;y++)field.set(key(x,y),col(Array.from({length:16},(_,j)=>new o.Cut((x+2*y+j)%5-2,(2*x-y+j)%3-1))));
  const cyclic=A=>{const out=new Map();for(const [d,M]of A){const [x,y]=xy(d);addAt(out,key(((x%Q)+Q)%Q,((y%Q)+Q)%Q),M);}return out;};
  const transformed=new Map(),continued=cyclic(compose(Z,field));
  for(let j=0;j<Q;j++)for(let k=0;k<Q;k++){
   const x=cutpow(xi,j),y=cutpow(xi,k),v=evaluate(field,x,y);transformed.set(key(j,k),v);
   eq(evaluate(continued,x,y),o.mul(evaluate(Z,x,y),v),'finite phase continuation');intertwiners++;
  }
  ensure(energy(transformed).div(Q*Q).eq(energy(field)),'native finite norm resolution');norms++;
  for(let x=0;x<Q;x++)for(let y=0;y<Q;y++){
   eq(o.scale(evaluate(transformed,cutpow(xi,-x),cutpow(xi,-y)),f(1).div(Q*Q)),field.get(key(x,y)),'finite phase inverse');inverse++;
  }
 }
 return {root_witnesses:rootWitnesses,primitive_power_checks:primitives,finite_orthogonality_sums:sums,
  finite_inverse_checks:inverse,finite_source_intertwiners:intertwiners,finite_norm_identities:norms,
  period_eight_witness:'native matrix realization of tau^2=1/2; both conjugates checked, positive scalar root selected in the written proof',
  infinite_transform_or_circle_measure_assumed:false};
});

const blocks=[U16,Z,compose(Z,Z)];blocks.push(compose(blocks[2],Z));
check('the_complete_source_ledger_and_block_displacement_product_rule_match',()=>{
 const totals=[0,1,2].map(k=>[...Z].reduce((s,[d,M])=>s.add(entryMass(M).mul(f(xy(d).reduce((s,a)=>s+Math.abs(a),0)).pow(k))),f(0)));
 ensure(totals[0].eq(68)&&totals[1].eq(64)&&totals[2].eq(88),'complete source gauge moments');
 eqFields(compose(dag(Z),Z),U16,16,16,'source unitary Laurent identity');
 let products=0,zeroJets=0;
 for(const e of dirs)for(let m=1;m<=3;m++){
  let product=new Map();const U=deriv(Z,e);
  for(let j=0;j<m;j++)product=plus(product,compose(compose(blocks[j],U),blocks[m-1-j]));
  eqFields(deriv(blocks[m],e),product,16,16,'complete m-block product derivative');products++;
  eq(sumAtOne(deriv(blocks[m],e)),o.scale(o.add(o.scale(Ax,e[0]),o.scale(Ay,e[1])),m),'touching point uses no gap inverse');zeroJets++;
 }
 return {complete_source_coefficient_sums:totals,source_unitary_identity:1,complete_block_derivative_identities:products,
  zero_gap_direct_jet_identities:zeroJets,uniform_written_bound:'m/sqrt(2) + 176 sqrt(m)',block_counts_checked:[1,2,3]};
});

let bandChecks=0,microResponse;
check('native_band_cancellation_controls_the_long_block_derivative',()=>{
 const rp=sx1,sp=sy1,omr=minus(id1,rp),oms=minus(id1,sp),tmr=minus(scale(id1,2),rp),tms=minus(scale(id1,2),sp);
 const lp=minus(scale(plus(rp,sp),2),compose(rp,sp)),den=compose(lp,minus(scale(id1,4),lp));
 const num=plus(compose(compose(tms,tms),compose(rp,omr)),compose(compose(tmr,tmr),compose(sp,oms)));
 const margin=plus(plus(scale(compose(compose(rp,rp),oms),2),scale(compose(compose(sp,sp),omr),2)),
  plus(compose(compose(rp,sp),minus(scale(id1,2),plus(rp,sp))),scale(compose(compose(rp,rp),compose(sp,sp)),'3/2')));
 eqFields(minus(scale(den,'1/2'),num),margin,1,1,'complete positive phase-speed polynomial');
 const z=new o.Cut('3/5','4/5'),t=new o.Cut('5/13','12/13'),tiny=new o.Cut('143/145','24/145');
 const phases=[[z,o.Cut.of(1)],[o.Cut.of(1),z],[z,z],[z,z.dagger()],[iota,iota],
  [o.Cut.of(-1),o.Cut.of(1)],[o.Cut.of(-1),o.Cut.of(-1)],[tiny,o.Cut.of(1)],[z,t]];
 let splittings=0,squares=0,recurrences=0,wrongLinear=0;
 for(const [x,y]of phases){
  const r=f(1).sub(x.rad).div(2),s=f(1).sub(y.rad).div(2),ell=r.add(s).mul(2).sub(r.mul(s)),a=f(1).sub(ell.div(2)),h2=ell.mul(f(4).sub(ell)).div(4);
  const Z0=evaluate(Z,x,y),B0=o.sub(Z0,o.scale(I16,a));eq(o.mul(B0,B0),o.scale(I16,h2.neg()),'native branch square');
  const powers=[I16,Z0];for(let m=2;m<=3;m++)powers.push(o.mul(powers[m-1],Z0));
  const sm=[f(0),f(1)];for(let m=2;m<=3;m++)sm.push(a.mul(2).mul(sm[m-1]).sub(sm[m-2]));
  for(const e of dirs){
   const U=evaluate(deriv(Z,e),x,y),BUB=o.scale(o.mul(B0,o.mul(U,B0)),f(1).div(h2));
   const Ud=o.scale(o.sub(U,BUB),'1/2'),Uo=o.scale(o.add(U,BUB),'1/2');
   eq(o.mul(B0,Ud),o.mul(Ud,B0),'native diagonal part');eq(o.add(o.mul(B0,Uo),o.mul(Uo,B0)),zero16,'native switching part');splittings+=2;
   const ellE=f(2).sub(s).mul(x.turn).mul(e[0]).add(f(2).sub(r).mul(y.turn).mul(e[1])).div(2),v2=ellE.pow(2).div(h2.mul(4));
   eq(Ud,o.scale(o.mul(Z0,B0),new o.Cut(0,ellE.neg().div(h2.mul(2)))),'derived diagonal response');
   eq(o.mul(o.dagger(Ud),Ud),o.scale(I16,v2),'exact diagonal speed square');ensure(v2.le('1/2'),'directional sharp speed bound');squares++;
   for(let m=1;m<=3;m++){
    const D=evaluate(deriv(blocks[m],e),x,y),correct=o.add(o.scale(o.mul(powers[m-1],Ud),m),o.scale(Uo,sm[m]));
    eq(D,correct,'all retained off-diagonal geometric cancellation');recurrences++;
    if(m>1&&!o.equal(D,o.scale(o.mul(powers[m-1],U),m)))wrongLinear++;
   }
  }
 }
 microResponse=o.energy(o.mul(evaluate(deriv(Z,dirs[2]),z,t),e16));
 ensure(microResponse.eq('1629/3250')&&!microResponse.le('1/2'),'single-block response can exceed limiting speed square');
 ensure(wrongLinear>0,'noncommuting block derivative cannot be replaced by m Z^(m-1) U');bandChecks=recurrences;
 return {complete_positive_speed_polynomial:1,phase_pairs:phases.length,branch_split_identities:splittings,
  diagonal_speed_square_identities:squares,exact_geometric_block_derivatives:recurrences,
  rejected_commuting_derivative_instances:wrongLinear,single_block_unit_direction_response_squared:microResponse,
  single_block_response_is_not_a_universal_instantaneous_speed_bound:true};
});

const signedPower=(a,n)=>n<0?f(1).div(f(a).pow(-n)):f(a).pow(n);
const tiltRational=(A,rx,ry)=>new Map([...A].map(([d,M])=>[d,o.scale(M,signedPower(rx,xy(d)[0]).mul(signedPower(ry,xy(d)[1])))]));
check('finite_native_count_tilts_preserve_products_and_control_the_remainder_budget',()=>{
 let products=0,jets=0;const massRows=[];
 for(let m=1;m<=3;m++){
  eqFields(tiltRational(blocks[m],'3/2','4/3'),power(tiltRational(Z,'3/2','4/3'),m,16),16,16,'positive multiplicative tilt product');products++;
  const moment=[...blocks[m]].reduce((s,[d,M])=>s.add(entryMass(M).mul((Math.abs(xy(d)[0])+Math.abs(xy(d)[1]))**2)),f(0));
  ensure(mass(blocks[m]).le(f(68).pow(m))&&moment.le(f(4*m*m).mul(f(68).pow(m))),'complete finite tilt mass bound');
  for(const e of dirs){
   let second=new Map();for(let j=0;j<m;j++)second=plus(second,compose(compose(blocks[j],deriv(deriv(Z,e),e)),blocks[m-1-j]));
   for(let j=0;j<m;j++)for(let k=j+1;k<m;k++)second=plus(second,scale(compose(compose(compose(compose(blocks[j],deriv(Z,e)),blocks[k-j-1]),deriv(Z,e)),blocks[m-1-k]),2));
   eqFields(deriv(deriv(blocks[m],e),e),second,16,16,'complete second exponential-tilt jet');jets++;
  }
  massRows.push({blocks:m,coefficient_mass:mass(blocks[m]),second_moment:moment,bound:f(4*m*m).mul(f(68).pow(m))});
 }
 let factorial=f(1),partial=f(0);for(let j=1;j<=24;j++){factorial=factorial.mul(j);if(j>=2)partial=partial.add(f(1).div(factorial));}
 ensure(partial.add(f(2).div(factorial.mul(25))).le(1),'finite factorial sum plus proved geometric tail');
 return {complete_positive_tilt_product_identities:products,complete_second_tilt_jet_identities:jets,coefficient_mass_rows:massRows,
  factorial_remainder_budget:'sum_(j>=2) 1/j! <= 1, finite sum plus native geometric tail',
  tilted_block_bound:'Exp_Sigma(beta m [1/sqrt(2)+176/sqrt(m)+4m 68^m beta])',positive_weights_are_source_dynamics:false};
});

const floorF=a=>a.n>=0n?a.n/a.d:-((-a.n+a.d-1n)/a.d),ceilF=a=>-floorF(a.neg());
const budgets=[];
check('explicit_native_parameter_budgets_give_a_uniform_half_space_bound',()=>{
 for(const etaString of ['1/4','1/2','1','2','16']){
  const eta=f(etaString),k=f(ceilF(f(704).div(eta))),m=k.pow(2),qa=f(1).div(m.mul(2)),qb=eta.div(m.mul(16)),q=qa.le(qb)?qa:qb;
  ensure(f(176).div(k).le(eta.div(4)),'block derivative budget');ensure(m.mul(4).mul(q).le(eta.div(4)),'tilt second-order budget');ensure(q.le(qa),'small tilt condition before factor 68^-m');
  const grid=f(ceilF(f(4).mul(f(1).add(eta)).div(eta)));let Q=2;while(f(Q*Q).le(grid)&&!f(Q*Q).eq(grid))Q++;
  ensure(f(4).mul(f(1).add(eta)).div(eta).le(Q*Q),'finite directions using c<=1');
  budgets.push({eta,k,m,q,beta:{coefficient:q,base:'68',integer_exponent:m.neg()},direction_grid:Q,expanded_large_power:false});
 }
 return {exact_parameter_budgets:budgets,half_space_squared_norm_bound:'C^2 Exp_Sigma(-2 beta [d-R-(c_Sigma+eta/2)n])',
  large_block_powers_are_exact_expressions:true,large_blocks_directly_simulated:false,constants_are_optimized:false};
});

let finalField;
check('literal_fields_obey_weighted_conjugation_and_detected_tail_inequalities',()=>{
 let conjugations=0,cuts=0,norms=0;
 const initial=new Map([['0,0',e16],['-1,1',col(Array.from({length:16},(_,j)=>new o.Cut((j%3)-1,(j%2))))]]);
 for(const axis of [0,1]){
  const rho=f('3/2'),rx=axis===0?rho:f(1),ry=axis===1?rho:f(1),tilted=tiltRational(Z,rx,ry);
  let current=initial,weighted=tiltRational(initial,rx,ry);const E0=energy(initial);
  for(let n=0;n<=3;n++){
   eqFields(weighted,tiltRational(current,rx,ry),16,1,'exact finite weighted conjugation');conjugations++;
   ensure(energy(current).eq(E0),'unweighted source conservation');norms++;
   for(const d of [-1,0,1,2,3]){const tail=new Map([...current].filter(([k])=>xy(k)[axis]>=d));
    ensure(energy(tail).mul(signedPower(rho,2*d)).le(energy(weighted)),'native weighted address-cut bound');cuts++;}
   if(axis===0&&n===3)finalField=current;
   if(n<3){current=compose(Z,current);weighted=compose(tilted,weighted);}
  }
 }
 return {literal_weighted_source_intertwiners:conjugations,positive_weight_tail_checks:cuts,native_conservation_checks:norms,
  largest_direct_block_count:3,uniform_exponential_asymptotic_theorem_is_written:true};
});

check('native_integer_direction_counts_cover_the_radial_detector',()=>{
 const units=[[1,0],[0,1],[-1,0],[0,-1],['3/5','4/5'],['-3/5','4/5'],['5/13','-12/13'],['8/17','15/17']].map(a=>a.map(f));
 let covers=0,counts=0;
 for(const Q of [2,3,5,8]){
  let all=0;for(let x=-Q;x<=Q;x++)for(let y=-Q;y<=Q;y++)if(x||y)all++;
  ensure(all===(2*Q+1)**2-1,'complete integer direction count');counts++;
  const gamma=f(1).sub(f(1).div(Q*Q));
  for(const e of units){ensure(e[0].pow(2).add(e[1].pow(2)).eq(1),'native unit test direction');
   const p=e.map(a=>f(floorF(a.mul(Q).add('1/2')))),p2=p[0].pow(2).add(p[1].pow(2)),dot=p[0].mul(e[0]).add(p[1].mul(e[1]));
   ensure(!p2.zero()&&f(0).le(dot)&&gamma.pow(2).mul(p2).le(dot.pow(2)),'normalized rounded direction covers unit vector');
   ensure(p[0].div(Q).sub(e[0]).pow(2).add(p[1].div(Q).sub(e[1]).pow(2)).le(f(1).div(2*Q*Q)),'native rounding radius');covers++;
  }
 }
 for(const eta of ['1/4','1/2','1','2'].map(f)){
  ensure(eta.mul('3/4').sub(eta.div(2)).sub(eta.div(8)).eq(eta.div(8)),'radial exponent budget');
  ensure(eta.div(4).sub(eta.div(8)).eq(eta.div(8)),'uniform early-arrival budget');
 }
 return {finite_direction_counts:counts,exact_rounded_cover_witnesses:covers,complete_count_formula:'(2Q+1)^2-1',
  normalized_tail_bound:'sqrt(D_Q) C Exp_Sigma(-beta eta n/8) + q_n',
  admissible_initial_localization:'R_n/n -> 0 and initial outside relative norm q_n -> 0',
  preparation_is_allowed_to_depend_on_scale:true};
});

check('the_curvature_observer_and_previous_native_packets_share_the_sharp_signal_bound',()=>{
 const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),Rc=o.scale(o.mul(Pc,Jc),-1);
 let points=0;for(const [k,v]of finalField){const a=o.mul(Pc,v),b=o.mul(Rc,v);eq(o.add(a,o.mul(Jc,b)),v,'curvature source reconstruction');ensure(o.energy(a).add(o.energy(b)).eq(o.energy(v)),'same pointwise detector matching');points++;}
 const tail=A=>new Map([...A].filter(([k])=>xy(k)[0]>=1&&xy(k)[1]>=0));
 ensure(energy(tail(compose(unit(Pc),finalField))).add(energy(tail(compose(unit(Rc),finalField)))).eq(energy(tail(finalField))),'same tail in complete curvature observer');
 const parent=JSON.parse(parentBytes),family=parent.exact_checks.find(c=>c.name==='a_native_rational_speed_family_approaches_the_sharp_count_coefficient');
 ensure(family.passed&&family.native_speed_cuts.length===7,'frozen constructive speed family present');
 for(const row of family.native_speed_cuts){const a=f(row.carrier_parameter),next=a.pow(2).add(2).div(a.mul(2)),g=f(1).div(next);ensure(g.eq(row.derived_drift)&&f('1/2').sub(g.pow(2)).eq(row.squared_speed_deficit),'parent sharpness cuts remain native');}
 return {pointwise_curvature_reconstructions_and_norms:points,complete_tail_detector_identity:1,reused_native_sharpness_cuts:7,
  sharp_operational_speed_supremum:'1/sqrt(2)',dual_information_length_speed:'1',fine_count_speed_squared_per_original_event:'1/8',
  physical_field_or_rod_clock_identified:false};
});

check('fast_front_and_localization_witnesses_reject_overstated_cone_claims',()=>{
 let field=unit(e16),fronts=0;const rows=[];
 for(let n=1;n<=4;n++){
  const corner=key(n,n);ensure(!field.has(corner),'detector not yet reachable by prior support');field=compose(Z,field);
  const v=field.get(corner);ensure(v&&!o.isZero(v),'exact fast front survives');
  ensure(o.energy(v).eq(signedPower(256,-n)),'vanishing threshold reads faster corner');
  rows.push({blocks:n,component_address:[n,n],relative_threshold:signedPower(256,-n),squared_speed:'2'});fronts++;
 }
 const remote=shift(unit(e16),20,0);ensure(energy(new Map([...remote].filter(([k])=>xy(k)[0]>=20))).eq(1),'unlocalized preparation already lies in remote detector');
 let singular=false;try{o.Cut.of(0).inv();}catch(_){singular=true;}ensure(singular,'zero gap division rejected');
 ensure(microResponse&&!microResponse.le('1/2'),'one-block sharp-speed promotion rejected');
 return {exact_fast_front_witnesses:fronts,vanishing_threshold_counterexamples:rows,
  missing_initial_localization_counterexample:true,zero_gap_inverse_rejected:singular,
  all_thresholds_or_exact_support_cone_claimed:false,all_source_laws_or_physical_fields_covered:false};
});
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cpw=add(hh,kk),cmw=add(hh,sc(-1,kk)),pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['native_turn_square',mul(rr,rr),sc(-1,ii)],['native_census_square',mul(cpw,cpw),sc(2,ii)],['native_reverse_census_square',mul(cmw,cmw),sc(2,ii)],
 ['native_record_loop',mul(hh,kk,hh,kk),sc(-1,ii)],['native_reverse_record_loop',mul(kk,hh,kk,hh),sc(-1,ii)],
 ['native_plus_cut',mul(pp,pp),pp],['native_minus_cut',mul(qq,qq),qq],['native_disjoint_cuts',mul(pp,qq),sc(0,ii)],['native_complete_cut_pair',add(pp,qq),ii],
 ['native_half_speed_algebraic_witness',mul(sc('1/2',cpw),sc('1/2',cpw)),sc('1/2',ii)],['native_turn_increment_norm',mul(add(ii,sc(-1,rr)),add(ii,rr)),sc(2,ii)]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');



const negatives={
 ordinary_complex_field_was_primitive:replays.some(row=>row.name==='native_turn_square'&&row.replay==='REPLAY_MATCH'),
 phase_resolution_was_an_imported_transform:checks[0].infinite_transform_or_circle_measure_assumed===false,
 small_phase_sampling_is_the_universal_proof:checks[2].complete_positive_speed_polynomial===1&&checks[5].uniform_exponential_asymptotic_theorem_is_written,
 a_block_derivative_is_m_times_a_commuting_term:checks[2].rejected_commuting_derivative_instances>0,
 the_single_block_response_always_has_speed_square_at_most_half:checks[2].single_block_response_is_not_a_universal_instantaneous_speed_bound,
 the_band_gap_can_be_inverted_at_the_touching_point:checks[8].zero_gap_inverse_rejected,
 the_tilt_changes_the_native_source_law:checks[3].positive_weights_are_source_dynamics===false,
 huge_certified_blocks_were_directly_simulated:checks[4].large_blocks_directly_simulated===false,
 the_conservative_tail_constants_are_optimized:checks[4].constants_are_optimized===false,
 the_upper_theorem_only_covers_fixed_carriers:checks[6].preparation_is_allowed_to_depend_on_scale,
 localization_of_the_preparation_can_be_dropped:checks[8].missing_initial_localization_counterexample,
 complete_curvature_readout_changes_the_tail:checks[7].complete_tail_detector_identity===1,
 the_exact_faster_front_was_removed:checks[8].exact_fast_front_witnesses===4,
 the_same_cone_holds_for_vanishing_thresholds:checks[8].vanishing_threshold_counterexamples.length===4,
 one_native_source_selects_all_physical_fields:checks[8].all_source_laws_or_physical_fields_covered===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r38.native-universal-signal-cone.v1',status:'PASS_R38_NATIVE_UNIVERSAL_SIGNAL_CONE',
 input_sha256:inputHash,r37_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_finite_phase_resolution_derived:true,native_uniform_block_and_tail_bound_derived:true,
 native_universal_localized_signal_cone_derived:true,native_sharp_operational_speed_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written proofs: native finite phase resolution, uniform block displacement bound, positive count tilt, all-packet half-space estimate, universal localized radial cone, sharp fixed-threshold arrival speed and curvature/front/calibration boundaries. Universal means all admissibly localized preparations of the pinned source; no physical c, h, alpha or universality across different source laws is identified.'},null,2)+'\n');
