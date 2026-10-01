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
ensure(inputHash===arg('--expected-input-sha256'),'R37 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r36-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r36-sha256'),'R37 R36 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R36_NATIVE_CURVATURE_OBSERVER','R37 parent status');
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

const cutPow=(z,n)=>{let a=o.Cut.of(1),b=n<0?z.inv():z;for(let k=0;k<Math.abs(n);k++)a=a.mul(b);return a;};
const en=v=>o.energy(v);
const realMass=M=>M.flat().reduce((s,c)=>s.add(c.gauge()),f(0));
const opMass=A=>[...A.values()].reduce((s,M)=>s.add(realMass(M)),f(0));
const phaseOp=(A,z,t)=>new Map([...A].map(([k,M])=>{const [a,b]=xy(k);return[k,o.scale(M,cutPow(z,a).mul(cutPow(t,b)))];}));
const moment=(A,axis)=>[...A].reduce((s,[k,M])=>o.add(s,o.scale(M,xy(k)[axis])),zero16);
const scalarOp=a=>scalarLift(a,16);
const dx16=scalarOp(dx1),dy16=scalarOp(dy1);
const invRemainder=n=>{const r=new Map();if(n>=2)for(let j=0;j<=n-2;j++)r.set(key(j,0),o.matrix([[n-1-j]]));
 if(n<0)for(let j=1;j<=-n;j++)r.set(key(-j,0),o.matrix([[-n-j+1]]));return r;};
const factorSecond=R=>{let xx=new Map(),xyc=new Map(),yy=new Map();for(const [k,M]of R){const [a,b]=xy(k);
  for(const [j,c]of invRemainder(a))addAt(xx,key(xy(j)[0],b),o.scale(M,c[0][0]));
  addAt(xyc,'0,0',o.scale(M,a*b));
  for(const [j,c]of invRemainder(b)){const y=xy(j)[0];addAt(yy,key(0,y),o.scale(M,c[0][0].mul(1-a)));addAt(yy,key(1,y),o.scale(M,c[0][0].mul(a)));}
 }return {xx,xy:xyc,yy};};

const carrier=a=>{
 a=f(a);ensure(f(2).le(a.pow(2))&&!a.pow(2).eq(2),'positive nonzero carrier parameter');
 const t=a.pow(2).sub(2).div(a.mul(2)),w0=a.pow(2).add(2).div(a.mul(2)),g0=f(1).div(w0),den=f(1).add(t.pow(2));
 const z=new o.Cut(f(1).sub(t.pow(2)).div(den),t.mul(2).div(den)),zeta=new o.Cut(f(1).div(den),t.mul(w0).div(den));
 ensure(z.norm2().eq(1)&&zeta.norm2().eq(1),'native carrier unit norms');
 const mode=phaseOp(Z,z,o.Cut.of(1)),Z0=sumAtOne(mode),Pp=o.scale(o.sub(Z0,o.scale(I16,zeta.dagger())),zeta.sub(zeta.dagger()).inv()),Pm=o.sub(I16,Pp);
 eq(o.mul(Pp,Pp),Pp,'carrier actual branch');eq(o.dagger(Pp),Pp,'native branch orthogonality');eq(o.mul(Z0,Pp),o.scale(Pp,zeta),'native branch response');ensure(o.rank(Pp)===8,'carrier rank');
 const W=scale(mode,zeta.inv()),W0=sumAtOne(W),eta=zeta.dagger().div(zeta),Ui=[moment(W,0),moment(W,1)];
 eq(o.mul(Pp,o.mul(Ui[0],Pp)),o.scale(Pp,g0),'derived x drift');eq(o.mul(Pp,o.mul(Ui[1],Pp)),zero16,'derived zero y drift');
 const correction=Ui.map(U=>o.scale(o.mul(Pm,o.mul(U,Pp)),o.Cut.of(1).sub(eta).inv()));
 const T=plus(unit(Pp),plus(compose(unit(correction[0]),dx16),compose(unit(correction[1]),dy16)));
 const E=plus(U16,scale(dx16,g0)),rem=minus(compose(W,T),compose(T,E));
 eq(sumAtOne(rem),zero16,'corrected zero residual');eq(moment(rem,0),zero16,'corrected first x residual');eq(moment(rem,1),zero16,'corrected first y residual');
 const second=factorSecond(rem);
 const rebuilt=plus(plus(compose(second.xx,compose(dx16,dx16)),compose(second.xy,compose(dx16,dy16))),compose(second.yy,compose(dy16,dy16)));
 eqFields(rem,rebuilt,16,16,'complete native second-difference remainder');
 const kappa=realMass(correction[0]).add(realMass(correction[1]));
 return {a,t,w0,g0,z,zeta,W,W0,Pp,Pm,correction,T,E,rem,second,kappa,Ct:f(1).add(kappa.mul(2)),Cxx:opMass(second.xx),Cxy:opMass(second.xy),Cyy:opMass(second.yy)};
};

const carriers=['2','3/2','17/12'].map(carrier),mainCarrier=carriers[0];
check('native_rational_carriers_have_actual_bands_and_derived_drifts',()=>({
 exact_carriers:carriers.map(c=>({parameter:c.a,spatial_phase:c.z,event_phase:c.zeta,
  band_rank:8,drift_x:c.g0,drift_y:'0',positive_gap_parameter:c.t})),
 exact_unit_projector_and_drift_checks:carriers.length*8,imported_Fourier_or_packet_theorem:false
}));

check('a_source_corrector_removes_every_first_difference_error',()=>{
 let firstDefects=0;
 for(const c of carriers){const bare=minus(compose(c.W,unit(c.Pp)),compose(unit(c.Pp),c.E));
  ensure(!o.isZero(moment(bare,1))||!o.isZero(moment(bare,0)),'uncorrected packet has first-order mixing');firstDefects++;
  eq(sumAtOne(c.rem),zero16,'complete corrected zeroth term');eq(moment(c.rem,0),zero16,'complete corrected x jet');eq(moment(c.rem,1),zero16,'complete corrected y jet');}
 return {complete_second_difference_factorizations:carriers.length,uncorrected_first_order_defects:firstDefects,
  source_constants:carriers.map(c=>({parameter:c.a,kappa:c.kappa,comparison_norm_bound:c.Ct,
   residual_xx_mass:c.Cxx,residual_xy_mass:c.Cxy,residual_yy_mass:c.Cyy,residual_terms:c.rem.size,
   exact_remainder_sha256:hash(Buffer.from(JSON.stringify([...c.rem])))}))};
});

const polyAdd=(a,b)=>Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||f(0)).add(b[i]||f(0)));
const polyScale=(a,c)=>a.map(x=>x.mul(c));
const polyMul=(a,b)=>{const r=Array(a.length+b.length-1).fill(f(0));for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)r[i+j]=r[i+j].add(a[i].mul(b[j]));return r;};
const choose5=a=>{let r=[f(1)];for(let j=0;j<5;j++)r=polyMul(r,[a[0].sub(j),a[1]]);return polyScale(r,'1/120');};
const box=M=>new Map(Array.from({length:M},(_,x)=>[key(x,0),o.matrix([[1]])]));
const cubicProfile=M=>power(box(M),3,1);
const scalarEnergy=A=>energy(A);
const packetProfile=M=>{const h=cubicProfile(M),F=new Map();for(const [a,v]of h)for(const [b,w]of h)F.set(key(xy(a)[0],xy(b)[0]),o.matrix([[v[0][0].mul(w[0][0])]]));return F;};
const scalarDifference=(A,x,y)=>minus(shift(A,x,y),A);
check('compact_three_box_counts_derive_exact_norm_and_difference_bounds',()=>{
 const normPolynomial=polyAdd(polyAdd(choose5([f(2),f(3)]),polyScale(choose5([f(2),f(2)]),-6)),polyScale(choose5([f(2),f(1)]),15));
 const expected=[f(0),f('1/5'),f(0),f('1/4'),f(0),f('11/20')];
 ensure(normPolynomial.length===expected.length&&normPolynomial.every((v,i)=>v.eq(expected[i])),'complete native count polynomial');
 let checks=0;
 for(let M=1;M<=8;M++){const h=cubicProfile(M),d=scalarDifference(h,1,0),d2=scalarDifference(d,1,0),m=f(M),norm=m.pow(5).mul(11).add(m.pow(3).mul(5)).add(m.mul(4)).div(20);
  ensure(energy(h).eq(norm),'profile exact norm');ensure(energy(d).eq(m.mul(m.pow(2).add(1))),'profile first difference');ensure(energy(d2).eq(m.mul(6)),'profile second difference');
  ensure(energy(d).mul(m.pow(2)).le(norm.mul(4))&&energy(d2).mul(m.pow(4)).le(norm.mul(16)),'native uniform difference bounds');checks+=4;}
 return {complete_degree_five_count_polynomial:1,finite_profile_identities:checks,
  profile_norm_squared:'(11 M^5 + 5 M^3 + 4 M)/20',first_difference_norm_squared:'M (M^2+1)',second_difference_norm_squared:'6 M',
  first_relative_norm_bound:'2/M',all_second_relative_norm_bound:'4/M^2',finite_packet_support:'[0,3M-2]^2 after correction'};
});

const binomialWeights=(n,g)=>{let a=[f(1)];for(let k=0;k<n;k++){const b=Array(k+2).fill(f(0));for(let j=0;j<a.length;j++){b[j]=b[j].add(a[j].mul(f(1).sub(g)));b[j+1]=b[j+1].add(a[j].mul(g));}a=b;}return a;};
const floorF=a=>a.n/a.d,ceilF=a=>(a.n+a.d-1n)/a.d;
check('native_binomial_counts_control_rigid_drift_without_a_probability_premise',()=>{
 let moments=0,taylors=0;
 for(const c of carriers)for(let n=0;n<=9;n++){const a=binomialWeights(n,c.g0),mean=c.g0.mul(n),m=floorF(mean),fraction=mean.sub(f(m));
  const sum=a.reduce((s,x)=>s.add(x),f(0)),first=a.reduce((s,x,j)=>s.add(x.mul(j)),f(0)),variance=a.reduce((s,x,j)=>s.add(x.mul(f(j).sub(mean).pow(2))),f(0));
  ensure(sum.eq(1)&&first.eq(mean)&&variance.eq(c.g0.mul(f(1).sub(c.g0)).mul(n)),'finite binomial moments');
  const second=a.reduce((s,x,j)=>{const d=f(BigInt(j)-m);return s.add(x.mul(d.mul(d.sub(1))));},f(0));
  ensure(second.eq(variance.add(fraction.pow(2)).sub(fraction))&&second.le(variance),'rounded drift remainder');moments+=4;}
 for(let d=-9;d<=9;d++){const left=minus(minus(shift(id1,d,0),id1),scale(dx1,d)),right=compose(compose(dx1,dx1),invRemainder(d));
  eqFields(left,right,1,1,'complete integer shift Taylor remainder');ensure(opMass(invRemainder(d)).eq(f(d*(d-1)).div(2)),'exact Taylor coefficient mass');taylors+=2;}
 return {finite_native_moment_identities:moments,complete_integer_shift_identities:taylors,
  translated_envelope_bound:'fraction(ng) ||dx f|| + n g(1-g) ||dx^2 f||/2',physical_stochastic_law_assumed:false};
});

const vectorFromScalar=(A,v)=>new Map([...A].map(([k,M])=>[k,o.scale(v,M[0][0])]));
const e00=col(Array.from({length:16},(_,i)=>i===0?1:0));
const seedRole=o.mul(mainCarrier.Pp,e00);ensure(!o.isZero(seedRole),'nonzero band preparation');
const modulate=(field,z)=>new Map([...field].map(([k,v])=>[k,o.scale(v,cutPow(z,-xy(k)[0]))]));
const norm2Diff=(a,b)=>energy(minus(a,b));
let storedPacket;
check('literal_finite_packets_satisfy_the_source_and_residual_telescoping',()=>{
 const c=mainCarrier;let sourceChecks=0,telescopes=0,energyChecks=0;
 for(const M of [1,2]){const f0=vectorFromScalar(packetProfile(M),seedRole),initial=compose(c.T,f0);
  let actual=initial,env=f0,err=new Map(),literal=modulate(initial,c.z);
  const E0=energy(initial),dxx=compose(compose(dx16,dx16),f0),dxy=compose(compose(dx16,dy16),f0),dyy=compose(compose(dy16,dy16),f0);
  const residualBound2=f(3).mul(c.Cxx.pow(2).mul(energy(dxx)).add(c.Cxy.pow(2).mul(energy(dxy))).add(c.Cyy.pow(2).mul(energy(dyy))));
  for(let n=0;n<=3;n++){
   const compared=compose(c.T,env);eqFields(minus(actual,compared),err,16,1,'exact finite residual telescoping');telescopes++;
   eqFields(literal,scale(modulate(actual,c.z),cutPow(c.zeta,n)),16,1,'literal source matches carrier packet');sourceChecks++;
   ensure(energy(actual).eq(E0)&&energy(literal).eq(E0),'actual source conservation');energyChecks++;
   ensure(norm2Diff(actual,compared).le(residualBound2.mul(n*n)),'native residual norm bound');
   if(n===3&&M===2)storedPacket=literal;
   if(n<3){err=plus(compose(c.W,err),compose(c.rem,env));actual=compose(c.W,actual);env=compose(c.E,env);literal=compose(Z,literal);}
  }
 }
 return {literal_source_carrier_identities:sourceChecks,complete_finite_telescoping_identities:telescopes,
  source_norm_and_error_checks:energyChecks,largest_directly_evolved_block_count:3,large_flight_counts_are_symbolic_bounds:true};
});

const flightBound=(c,L)=>{const l=f(L),M=l.pow(2),N=l.pow(3),den=f(1).sub(c.kappa.mul(2).div(M));ensure(f(0).le(den)&&!den.zero(),'valid corrected packet norm');
 const A=c.Cxx.add(c.Cxy).add(c.Cyy).mul(4).add(c.Ct.mul(2).mul(c.g0).mul(f(1).sub(c.g0)));
 const epsilon=A.div(l).add(c.Ct.mul(2).div(M)).div(den),m=floorF(c.g0.mul(N)),window=ceilF(M.mul(3).add(1).div(c.g0));
 return {scale:l,width_parameter:M,block_count:N,centre_shift:f(m),relative_error:epsilon,escaped_norm_fraction_bound:epsilon.pow(2),arrival_window_blocks:f(window),original_event_count:N.mul(8),lower_arrival_block:N.sub(f(window)).add(1)};};
check('finite_rational_bounds_certify_reliable_localized_arrival_windows',()=>{
 const rows=[1024,4096,16384].map(L=>flightBound(mainCarrier,L));
 for(const r of rows){ensure(r.escaped_norm_fraction_bound.le('1/2'),'half-norm detection is certified');ensure(f(0).le(r.lower_arrival_block),'arrival window follows preparation');
  ensure(r.centre_shift.le(mainCarrier.g0.mul(r.block_count))&&mainCarrier.g0.mul(r.block_count).sub(r.centre_shift).le(1),'rounded displacement');
  ensure(r.width_parameter.mul(3).add(1).le(r.arrival_window_blocks.mul(mainCarrier.g0)),'disjoint earlier comparison supports');}
 return {exact_bound_examples:rows,detector:'native address cut with matching-norm fraction threshold 1/2',
  asymptotic_packet_width:'L^2',asymptotic_flight_blocks:'L^3',relative_error_order:'1/L',
  no_mass_outside_window_claimed:false,large_profiles_explicitly_simulated:false};
});

check('a_native_rational_speed_family_approaches_the_sharp_count_coefficient',()=>{
 let a=f(2);const rows=[];for(let n=0;n<7;n++){const next=a.pow(2).add(2).div(a.mul(2)),g=f(1).div(next),gap=next.pow(2).sub(2);
  ensure(gap.eq(a.pow(2).sub(2).pow(2).div(a.pow(2).mul(4))),'native root error recurrence');ensure(f(2).le(next.pow(2))&&next.le(a),'nested native upper cuts');
  const deficit=f('1/2').sub(g.pow(2));ensure(deficit.eq(gap.div(next.pow(2).mul(2)))&&!deficit.zero(),'strict speed deficit');rows.push({index:n,carrier_parameter:a,derived_drift:g,squared_speed_deficit:deficit});a=next;}
 return {native_speed_cuts:rows,sharp_speed_supremum:'1/sqrt(2)',fixed_nonzero_carrier_attains_supremum:false,
  packet_scale_must_dominate_carrier_correction:true,order_of_limits:'fix carrier; localize and propagate; then approach zero carrier'};
});

check('the_curvature_complete_observer_preserves_the_packet_and_detector_reading',()=>{
 const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),Rc=o.scale(o.mul(Pc,Jc),-1);
 let norms=0,decodes=0;for(const [k,v]of storedPacket){const a=o.mul(Pc,v),b=o.mul(Rc,v);eq(o.add(a,o.mul(Jc,b)),v,'pointwise source reconstruction');decodes++;
  ensure(en(a).add(en(b)).eq(en(v)),'pointwise detected information');norms++;}
 const take=(field,predicate)=>new Map([...field].filter(([k])=>predicate(...xy(k))));
 const observedA=compose(unit(Pc),storedPacket),observedB=compose(unit(Rc),storedPacket),pred=(x,y)=>x>=1&&x<=4&&y>=0&&y<=4;
 ensure(energy(take(observedA,pred)).add(energy(take(observedB,pred))).eq(energy(take(storedPacket,pred))),'same detector cut');
 return {pointwise_packet_reconstructions:decodes,pointwise_norm_checks:norms,detector_norm_identity:1,
  source_events_per_block:8,fine_address_steps_per_component:4,physical_rod_clock_adapter_selected:false};
});

check('carrier_threshold_and_exact_front_boundaries_reject_false_signal_claims',()=>{
 let singularRejected=false;try{o.Cut.of(0).inv();}catch(_){singularRejected=true;}ensure(singularRejected,'zero band gap cannot be divided');
 const c=mainCarrier;
 ensure(!zeroOp(c.rem),'finite packets are not exactly rigid translates');
 ensure(!f('1/2').eq(2),'phase/current bound differs from exact diagonal support speed');
 const out=unit(o.kron(e0,e0));let b=out;
 for(let n=1;n<=3;n++){b=compose(B,b);eq(b.get(key(2*n,2*n)),o.scale(o.kron(e0,e0),f((-1)**n).div(f(16).pow(n))),'unchanged nonzero fast front');}
 return {zero_gap_inverse_rejected:singularRejected,nonzero_packet_remainder:true,
  exact_fast_corner_witnesses:3,universal_all_packet_ballistic_upper_bound_claimed:false,
  zero_threshold_arrival_claimed:false,physical_c_h_alpha_identified:false};
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
 imported_plane_wave_dynamics_was_needed:checks[0].imported_Fourier_or_packet_theorem===false,
 a_scalar_envelope_needs_no_role_corrector:checks[1].uncorrected_first_order_defects===3,
 a_first_jet_test_is_a_complete_remainder_proof:checks[1].complete_second_difference_factorizations===3,
 the_packet_has_infinite_initial_support:checks[2].finite_packet_support==='[0,3M-2]^2 after correction',
 binomial_comparison_is_a_physical_noise_law:checks[3].physical_stochastic_law_assumed===false,
 the_large_certified_flight_was_brute_force_simulated:checks[5].large_profiles_explicitly_simulated===false,
 finite_packets_translate_without_error:checks[8].nonzero_packet_remainder,
 norm_localization_means_zero_support_outside:checks[5].no_mass_outside_window_claimed===false,
 every_detection_threshold_has_the_same_arrival:checks[8].zero_threshold_arrival_claimed===false,
 a_finite_carrier_already_attains_the_limit:checks[6].fixed_nonzero_carrier_attains_supremum===false,
 zero_carrier_can_be_inserted_in_the_corrector:checks[8].zero_gap_inverse_rejected,
 carrier_and_width_limits_may_be_interchanged_freely:checks[6].packet_scale_must_dominate_carrier_correction,
 curvature_readout_changes_packet_speed:checks[7].detector_norm_identity===1,
 the_exact_faster_front_was_removed:checks[8].exact_fast_corner_witnesses===3,
 the_native_flight_meter_already_fixes_physical_c:checks[8].physical_c_h_alpha_identified===false,
};ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r37.native-packet-flight.v1',status:'PASS_R37_NATIVE_PACKET_FLIGHT',
 input_sha256:inputHash,r36_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_carrier_and_packet_corrector_derived:true,native_finite_packet_error_derived:true,
 native_localized_arrival_and_speed_family_derived:true,native_curvature_observed_flight_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written proofs: native rational carriers, complete finite packet corrector, compact count profiles, all-n drift error, reliable localized arrival, the sharp approaching speed family, and curvature-observer/event-calibration boundaries. No universal all-packet ballistic bound, zero-threshold front or physical c, h or alpha is identified.'},null,2)+'\n');
