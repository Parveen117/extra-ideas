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
ensure(inputHash===arg('--expected-input-sha256'),'R33 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r32-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r32-sha256'),'R33 R32 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R32_NATIVE_RETAINED_LOOP_GAP','R33 parent status');
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

// The following 2x2 matrices witness scalar polynomial identities in sqrt(2).
// C0 itself is not declared positive. Rational source cuts isolate the positive
// scalar branch separately. Physical role-frame probes have their own factor.
const rho=o.sub(o.scale(I,3),o.scale(C0,2)),alpha=o.mul(rho,rho);
const zero=o.zeros(2),mprod=(...a)=>a.reduce(o.mul),mi=o.inverse;
const blocks=rows=>rows.flatMap(row=>[0,1].map(i=>row.flatMap(M=>M[i])));
const part=(M,ri,ci)=>ri.map(i=>ci.map(j=>M[i][j]));
const cp=r=>blocks([[o.scale(I,'1/2'),o.scale(r,'1/2')],[o.scale(r,'1/2'),o.scale(I,'1/2')]]);
const kp=r=>{const t=o.scale(mi(o.sub(I,o.mul(r,r))),2);return blocks([[t,o.scale(o.mul(t,r),-1)],[o.scale(o.mul(t,r),-1),t]]);};
const abs=x=>x.le(0)?f(0).sub(x):x;
let rootBracket,rootLo,rootHi,alphaLo,alphaHi;

check('native_increment_gram_supplies_the_cut_quadratic',()=>{
 eqFields(compose(dag(F),F),G,4,4,'source-derived increment Gram');
 eqFields(G,plus(unit(I4),plus(scale(plus(LX,LY),'1/4'),scale(compose(LX,LY),'1/16'))),4,4,'native positive square expansion');
 let count=0;for(const psi of probes){
  const e=energy(psi),rhs=e.add(energy(compose(DX,psi)).div(4)).add(energy(compose(DY,psi)).div(4)).add(energy(compose(DX,compose(DY,psi))).div(16));
  ensure(energy(compose(F,psi)).eq(rhs),'source increment norm agrees with derived squares');
  ensure(e.le(rhs)&&rhs.le(e.mul(4)),'native lower and upper Gram bounds');count+=2;
 }
 return {complete_Laurent_identities:2,finite_norm_checks:count,physical_quadratic_action_supplied:false};
});

check('native_two_cut_elimination_fixes_roundtrip_gain',()=>{
 let n=0;for(let ell=1;ell<=6;ell++){
  const r=o.power(rho,ell),c=cp(r),k=kp(r),t=o.scale(mi(o.sub(I,o.mul(r,r))),2);
  eq(o.mul(k,c),o.identity(4),'two-cut completed inverse');
  const J=blocks([[I,zero],[o.scale(r,-1),I]]),D=blocks([[o.scale(I,2),zero],[zero,t]]);
  eq(mprod(o.dagger(J),D,J),k,'two-cut square completion');
  const release=blocks([[I],[r]]);eq(o.mul(k,release),blocks([[o.scale(I,2)],[zero]]),'unique released second cut');
  eq(o.mul(r,r),o.power(alpha,ell),'closed coefficient is the square of the one-way transfer');n+=4;
 }
 let lo=f(1),hi=f('3/2');for(let k=0;k<64;k++){const mid=lo.add(hi).div(2);if(mid.pow(2).le(2))lo=mid;else hi=mid;}
 rootLo=f(3).sub(hi.mul(2));rootHi=f(3).sub(lo.mul(2));alphaLo=rootLo.pow(2);alphaHi=rootHi.pow(2);
 ensure(f(0).le(rootLo)&&!rootLo.zero()&&rootHi.le('1/5'),'positive decaying root');
 ensure(f('2943725152285/100000000000000').le(alphaLo)&&alphaHi.le('2943725152287/100000000000000'),'published alpha_rt decimal enclosure');
 rootBracket={sqrt2_lower:lo.toString(),sqrt2_upper:hi.toString(),rho_lower:rootLo.toString(),rho_upper:rootHi.toString(),alpha_rt_lower:alphaLo.toString(),alpha_rt_upper:alphaHi.toString()};
 return {exact_cut_and_release_identities:n,component_separations_checked:6,positive_root_bisections:64,
  positive_native_root_interval:rootBracket,elementary_roundtrip_amplitude:'17 - 12 sqrt(2)',one_way_squared_norm_fraction:'alpha_rt',roundtrip_squared_norm_fraction:'alpha_rt^2'};
});

check('native_closed_feedback_has_exact_geometric_remainder',()=>{
 let n=0;for(let ell=1;ell<=3;ell++){
  const a=o.power(alpha,ell),den=o.sub(I,a),inv=mi(den);let sum=o.zeros(2),power=I;
  for(let N=0;N<=7;N++){
   sum=o.add(sum,power);power=o.mul(a,power);
   eq(o.mul(den,sum),o.sub(I,power),'complete finite return cancellation');
   eq(o.sub(inv,sum),o.mul(power,inv),'exact completed return remainder');n+=2;
  }
 }
 return {exact_feedback_and_remainder_identities:n,maximum_return_cutoff:7,physical_time_evolution_of_relaxation_claimed:false};
});

check('native_collinear_coarsening_multiplies_memory',()=>{
 let n=0;for(const [l,m]of [[1,1],[1,2],[2,3],[3,4]]){
  const a=o.power(rho,l),b=o.power(rho,m),ab=o.mul(a,b),da=mi(o.sub(I,o.mul(a,a))),db=mi(o.sub(I,o.mul(b,b)));
  const c=o.scale(blocks([[I,a,ab],[a,I,b],[ab,b,I]]),'1/2');
  const k=o.scale(blocks([[da,o.scale(o.mul(a,da),-1),zero],[o.scale(o.mul(a,da),-1),o.add(da,mprod(b,b,db)),o.scale(o.mul(b,db),-1)],[zero,o.scale(o.mul(b,db),-1),db]]),2);
  eq(o.mul(k,c),o.identity(6),'three-cut inverse from ordered count kernel');
  const ends=[0,1,4,5],middle=[2,3],schur=o.sub(part(k,ends,ends),mprod(part(k,ends,middle),mi(part(k,middle,middle)),part(k,middle,ends)));
  eq(schur,kp(ab),'release middle cut equals direct elimination');
  eq(o.mul(ab,ab),o.power(alpha,l+m),'coarsened roundtrip is multiplied');n+=3;
 }
 ensure(!o.equal(rho,o.power(rho,3)),'backtracking does not add shortest endpoint length');
 return {exact_three_cut_and_coarsening_identities:n,backtracking_counterexample_checked:true,no_nontrivial_resolution_fixed_point_selected:true};
});

check('native_localized_load_has_exact_spread_fraction',()=>{
 const s=mprod(o.add(I,alpha),mi(o.sub(I,alpha))),total=o.scale(o.mul(s,s),'1/4');
 eq(o.mul(s,s),o.scale(I,'9/8'),'exact completed squared-kernel count');
 eq(total,o.scale(I,'9/32'),'completed point response norm');
 eq(o.sub(total,o.scale(I,'1/4')),o.scale(total,'1/9'),'off-origin response fraction');
 let n=0;for(let N=0;N<=7;N++){
  let finite=I;for(let k=1;k<=N;k++)finite=o.add(finite,o.scale(o.power(alpha,k),2));
  const tail=o.scale(o.mul(o.power(alpha,N+1),mi(o.sub(I,alpha))),2);
  eq(o.sub(s,finite),tail,'one-direction squared-count tail');
  eq(o.sub(o.mul(s,s),o.mul(finite,finite)),o.mul(tail,o.add(s,finite)),'two-direction exact omitted norm');n+=2;
 }
 return {completed_kernel_norm_identities:3,finite_squared_count_tail_identities:n,total_response_norm_coefficient:'9/32',loaded_cut_norm_coefficient:'1/4',increment_norm_coefficient:'1/2',nonlocal_response_fraction:'1/9'};
});

// Scalar Laurent polynomials in the even-address component shift, using only
// the canonical native rationals. No floating-point series or Fourier step.
const fa=(a,b)=>{const c=new Map(a);for(const [k,v]of b){const n=(c.get(k)||f(0)).add(v);if(n.zero())c.delete(k);else c.set(k,n);}return c;};
const fscl=(a,b)=>new Map([...a].map(([k,v])=>[k,v.mul(b)]));
const fm=(a,b)=>{let c=new Map();for(const [i,u]of a)for(const [j,v]of b)c=fa(c,new Map([[i+j,u.mul(v)]]));return c;};
const feq=(a,b,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))ensure((a.get(k)||f(0)).eq(b.get(k)||f(0)),msg);};
const one=new Map([[0,f(1)]]),average=new Map([[-1,f('1/6')],[1,f('1/6')]]),qpoly=new Map([[0,f('3/2')],[-1,f('-1/4')],[1,f('-1/4')]]);
let finiteRows=[];
check('finite_native_cut_extraction_encloses_completed_coefficient',()=>{
 let sum=new Map(),power=one,identities=0,bounds=0;
 for(let N=0;N<=12;N++){
  sum=fa(sum,power);power=fm(average,power);const pn=fscl(sum,'2/3');
  feq(fm(qpoly,pn),fa(one,fscl(power,-1)),'all Laurent coefficients of finite inverse remainder');identities++;
  if(N===0)continue;
  const a=(pn.get(0)||f(0)).pow(2),b=(pn.get(0)||f(0)).mul(pn.get(1)||f(0)),rn=b.div(a),an=rn.pow(2),tail=f(1).div(f(3).pow(N+1)),eps=tail.mul(2).add(tail.pow(2));
  const d=eps.mul('6/5').div(f('1/2').sub(eps)),alphaErr=d.mul(f('2/5').add(d)),invErr=eps.div(f('2/5').mul(f('2/5').sub(eps)));
  ensure(abs(a.sub('1/2')).le(eps),'cut diagonal enclosure');
  for(const r of [rootLo,rootHi]){
   ensure(abs(b.sub(r.div(2))).le(eps),'cut cross enclosure');
   ensure(abs(rn.sub(r)).le(d),'one-way coefficient enclosure');
   ensure(abs(an.sub(r.pow(2))).le(alphaErr),'roundtrip coefficient enclosure');
   for(const sign of [-1,1])ensure(abs(f(1).div(a.add(b.mul(sign))).sub(f(2).div(f(1).add(r.mul(sign))))).le(invErr),'inverse plus/minus cut enclosure');
   bounds+=5;
  }bounds++;
  finiteRows.push({N,one_way_rational:rn.toString(),roundtrip_rational:an.toString(),cut_matrix_error_bound:eps.toString(),roundtrip_error_bound:alphaErr.toString()});
 }
 return {complete_finite_Laurent_remainder_identities:identities,certified_cut_inverse_and_coefficient_bound_checks:bounds,finite_native_extractions:finiteRows};
});

const cyclicShift=L=>{const s=o.zeros(L);for(let n=0;n<L;n++)s[(n+1)%L][n]=o.Cut.of(1);return s;};
let cycleRows=[];
check('finite_cyclic_elimination_matches_wrapped_native_kernel',()=>{
 let wraps=0,squares=0;
 for(let L=2;L<=5;L++){
  const id=o.identity(L),s=cyclicShift(L),q=o.sub(o.scale(id,'3/2'),o.scale(o.add(s,o.dagger(s)),'1/4')),qi=mi(q),rL=o.power(rho,L),norm=o.mul(o.scale(C0,'1/2'),mi(o.sub(I,rL)));
  for(let n=0;n<L;n++){
   const kernel=o.mul(norm,o.add(o.power(rho,n),o.power(rho,L-n)));
   eq(kernel,o.scale(I,qi[n][0]),'native wrapped inverse coefficient');wraps++;
  }
  const gg=o.kron(q,q),ginv=mi(gg),size=L*L,e=o.zeros(size,2);e[0][0]=o.Cut.of(1);e[1][1]=o.Cut.of(1);
  const c=mprod(o.dagger(e),ginv,e),k=mi(c),u=col([2,-1]),phi=mprod(ginv,e,k,u),z=col(Array.from({length:size},(_,i)=>i<2?0:(i%5)-2));
  eq(o.mul(o.dagger(e),phi),u,'finite released interior preserves both cuts');
  eq(mprod(o.dagger(o.add(phi,z)),gg,o.add(phi,z)),o.add(mprod(o.dagger(u),k,u),mprod(o.dagger(z),gg,z)),'exact finite minimizing-square identity');squares+=2;
  const r=c[0][1].rad.div(c[0][0].rad);cycleRows.push({component_period:L,normalized_neighbor_transfer:r.toString(),roundtrip_coefficient:r.pow(2).toString()});
 }
 ensure(cycleRows[0].normalized_neighbor_transfer==='1/3','two-cycle boundary changes native coefficient');
 return {wrapped_inverse_coefficient_identities:wraps,finite_interior_elimination_identities:squares,cyclic_response_coefficients:cycleRows,unwrapped_coefficients_reused_on_a_cycle:false};
});

check('native_frame_changes_preserve_full_role_roundtrip',()=>{
 const rotation=o.add(o.scale(I,'3/5'),o.scale(R,'4/5'));
 const frames=[I4,o.kron(H,K),o.kron(rotation,H),o.kron(K,rotation)];let n=0;
 for(let k=0;k<frames.length;k++){
  const u=frames[k],v=frames[(k+1)%frames.length],t=o.kron(rho,o.mul(o.dagger(v),u));
  eq(o.mul(o.dagger(u),u),I4,'native frame isometry');
  eq(o.mul(o.dagger(t),t),o.kron(alpha,I4),'preparation-independent forward Gram and closed amplitude');
  eq(mprod(o.dagger(t),t,o.dagger(t),t),o.kron(o.mul(alpha,alpha),I4),'complete roundtrip energy has another square');n+=3;
 }
 const [u,v,z]=frames,tuv=o.kron(rho,o.mul(o.dagger(v),u)),tvz=o.kron(rho,o.mul(o.dagger(z),v)),tzu=o.kron(rho,o.mul(o.dagger(u),z));
 eq(mprod(tzu,tvz,tuv),o.kron(o.power(rho,3),I4),'response-frame triangle telescopes');
 return {exact_native_frame_and_roundtrip_identities:n,closed_response_triangle_identities:1,full_four_role_preparation_independent:true,signed_record_holonomy_recovered_from_gram:false};
});

let averagedTransfer;
check('native_readout_and_aperture_choices_change_scalar_memory',()=>{
 const c=o.scale(blocks(Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>o.power(rho,Math.abs(i-j))))),'1/2'),norm=o.scale(C0,'1/2');
 const aperture=blocks([[norm,zero],[norm,zero],[zero,norm],[zero,norm]]),ca=mprod(o.dagger(aperture),c,aperture);
 eq(o.mul(o.dagger(aperture),aperture),o.identity(4),'disjoint normalized aperture cuts');
 const diagonal=o.scale(o.add(I,rho),'1/2'),cross=o.scale(mprod(rho,o.add(I,rho),o.add(I,rho)),'1/4');
 eq(ca,blocks([[diagonal,cross],[cross,diagonal]]),'actual four-site aperture compression');
 averagedTransfer=o.mul(mi(diagonal),cross);eq(averagedTransfer,o.scale(o.mul(rho,o.add(I,rho)),'1/2'),'changed normalized aperture transfer');
 ensure(!o.equal(averagedTransfer,rho),'aperture does not preserve elementary transfer');
 eq(mprod(o.dagger(e0),e1),o.zeros(1),'single role channel can miss complete transfer');
 eq(mprod(o.dagger(e0),e0),o.identity(1),'aligned role channel retains full transfer');
 ensure([...G.keys()].every(k=>xy(k).every(n=>n%2===0)),'Gram separates address parity components');
 eq(mprod(H,K,H,K),o.scale(I,-1),'retained record holonomy remains minus identity');
 return {exact_aperture_and_channel_identities:6,parity_separation_checked:true,aperture_transfer:'rho(1+rho)/2',single_channel_response_range:'0 through alpha_rt',physical_alpha_selected:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cpw=add(hh,kk),rw=add(sc(3,ii),sc(-2,cpw)),aw=mul(rw,rw),oneMinus=add(ii,sc(-1,aw)),onePlus=add(ii,aw);
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[
 ['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['native_census_square',mul(cpw,cpw),sc(2,ii)],['native_turn_square',mul(rr,rr),sc(-1,ii)],['native_record_holonomy',mul(hh,kk,hh,kk),sc(-1,ii)],
 ['inverse_kernel_root',add(aw,sc(-6,rw),ii),sc(0,ii)],['inverse_kernel_normalization',mul(sc('1/2',cpw),add(ii,rw)),add(ii,sc(-1,rw))],
 ['reciprocal_native_roots',mul(rw,add(sc(3,ii),sc(2,cpw))),ii],['roundtrip_root_polynomial',add(mul(aw,aw),sc(-34,aw),ii),sc(0,ii)],
 ['squared_count_numerator',onePlus,sc(6,rw)],['squared_count_denominator',oneMinus,sc(4,mul(cpw,rw))],
 ['complete_squared_count_norm',sc(8,mul(onePlus,onePlus)),sc(9,mul(oneMinus,oneMinus))],
 ['roundtrip_exact_native_expression',aw,add(sc(17,ii),sc(-12,cpw))]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const negatives={
 roundtrip_amplitude_equals_one_way_amplitude:!o.equal(alpha,rho),
 roundtrip_energy_fraction_equals_roundtrip_amplitude:!o.equal(o.mul(alpha,alpha),alpha),
 full_role_fraction_depends_on_preparation:checks[7].full_four_role_preparation_independent,
 local_native_frames_change_the_closed_scalar:checks[7].exact_native_frame_and_roundtrip_identities===12,
 one_role_sensor_is_a_complete_readout:o.isZero(o.mul(o.dagger(e0),e1)),
 all_cut_separations_have_one_coupling:!o.equal(alpha,o.power(alpha,2)),
 this_nonzero_coefficient_is_resolution_invariant:!o.equal(alpha,o.power(alpha,3)),
 normalized_averaged_apertures_preserve_point_coupling:!o.equal(averagedTransfer,rho),
 total_spatial_memory_fraction_is_the_two_cut_alpha:!o.equal(alpha,o.scale(I,'1/9')),
 inverse_gram_transport_retains_minus_identity_holonomy:!o.equal(o.power(rho,4),o.scale(I,-1)),
 cyclic_boundary_has_the_unwrapped_elementary_coefficient:alphaHi.le('1/9')&&!alphaHi.eq('1/9'),
 a_generic_finite_inverse_series_is_the_completed_coefficient:f(finiteRows[0].roundtrip_rational).le(alphaLo)&&!f(finiteRows[0].roundtrip_rational).eq(alphaLo),
 growing_native_root_is_a_localized_tail:f(3).add(f(2).mul(f(rootBracket.sqrt2_lower))).le(1)===false,
 different_address_parities_have_nonzero_gram_transfer:[...G.keys()].every(k=>xy(k).every(n=>n%2===0)),
 feedback_denominator_uses_one_way_amplitude:!o.equal(o.sub(I,alpha),o.sub(I,rho)),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r33.native-interaction-memory.v1',status:'PASS_R33_NATIVE_INTERACTION_MEMORY',
 input_sha256:inputHash,r32_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_cut_elimination_and_roundtrip_invariant_derived:true,native_memory_coarsening_and_feedback_derived:true,
 native_complete_spread_fraction_derived:true,finite_extraction_and_selection_counterexamples_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,
 physical_metric_c_alpha_derived:false,physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,
 primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native results: increment-derived cut elimination, two-cut roundtrip invariant, exact closed feedback, coarsening by releasing cuts, complete response spread fraction, finite extraction and cyclic boundaries, and role/aperture/holonomy selection counterexamples. No physical alpha identification, electromagnetic normalization or empirical constant is supplied.'},null,2)+'\n');
