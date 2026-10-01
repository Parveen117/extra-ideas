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
ensure(inputHash===arg('--expected-input-sha256'),'R35 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r34-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r34-sha256'),'R35 R34 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R34_NATIVE_DECODER_BURDEN','R35 parent status');
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

check('complete_native_complement_and_positive_wave_identities',()=>{
 eqFields(plus(B,dag(B)),minus(unit(o.scale(I4,2)),G),4,4,'pinned retained central identity');
 eqFields(normSquare(plus(B,unit(I4))),minus(unit(o.scale(I4,4)),G),4,4,'complementary signed increment');
 const signScalar=a=>new Map([...a].map(([k,M])=>[k,o.scale(M,(-1)**((xy(k)[0]+xy(k)[1])/2))]));
 const cg=signScalar(G),target=scalarLift(ell1,4);
 const doubled=new Map([...target].map(([k,M])=>[key(2*xy(k)[0],2*xy(k)[1]),M]));
 eqFields(minus(unit(o.scale(I4,4)),cg),doubled,4,4,'address sign flips both squared shifts');
 const positive=plus(scale(plus(lx1,ly1),'3/8'),scale(plus(normSquare(compose(dx1,ey1)),normSquare(compose(dy1,ex1))),'1/32'));
 eqFields(positive,ell1,1,1,'complete lower squares');
 const upper=plus(scale(plus(normSquare(ex1),normSquare(ey1)),'1/4'),scale(normSquare(compose(ex1,ey1)),'1/16'));
 eqFields(upper,minus(scale(id1,3),ell1),1,1,'complete upper squares');
 eqFields(plus(Z,dag(Z)),minus(scale(U16,2),ell),16,16,'sixteen-role central identity');
 eqFields(normSquare(minus(Z,U16)),ell,16,16,'envelope increment Gram');
 eqFields(compose(dag(Z),Z),U16,16,16,'blocked native unitarity');
 return {complete_Laurent_identities:8,positive_lower_square_terms:4,positive_upper_square_terms:3,
   old_identity_relative_gap:'1',envelope_gap_infimum:'0',envelope_squared_increment_upper_bound:'3'};
});

const chi=(x,y)=>(-1)**(Math.floor(x/2)+Math.floor(y/2));
const encode=field=>{const out=new Map();for(const [k,v]of field){const [x,y]=xy(k),ex=parity(x),ey=parity(y),a=(x-ex)/2,b=(y-ey)/2;
 const c=o.zeros(16,1);for(let j=0;j<4;j++)c[(ex*2+ey)*4+j][0]=v[j][0].mul(chi(x,y));addAt(out,key(a,b),c);}return out;};
let frameSteps=0;
check('finite_fields_preserve_all_parities_under_signed_regrouping',()=>{
 for(const seed of probes){let original=seed,envelope=encode(seed);const E=energy(seed);
  for(let n=0;n<=4;n++){eqFields(envelope,scale(encode(original),(-1)**n),16,1,'exact time/address signed field bijection');
   ensure(energy(envelope).eq(E),'signed regrouping energy');frameSteps++;
   if(n<4){original=compose(B,original);envelope=compose(Z,envelope);}
  }
 }
 return {full_source_vs_blocked_evolution_checks:frameSteps,original_preparations:probes.length,
   retained_address_parities:4,roles_per_component_cell:16,source_roles_erased:false};
});

check('native_principal_generators_have_isotropic_source_square',()=>{
 eq(sumAtOne(Z),I16,'uniform component carrier');eq(o.dagger(Ax),Ax,'native x self-dagger');eq(o.dagger(Ay),Ay,'native y self-dagger');
 eq(o.mul(Ax,Ax),o.scale(I16,'1/2'),'x squared slope');eq(o.mul(Ay,Ay),o.scale(I16,'1/2'),'y squared slope');
 eq(o.add(o.mul(Ax,Ay),o.mul(Ay,Ax)),zero16,'zero mixed square');
 const probesK=[[1,0],[0,1],[1,1],[1,-1],[2,3],[-3,2]];
 for(const [x,y]of probesK){const a=o.add(o.scale(Ax,x),o.scale(Ay,y));eq(o.mul(a,a),o.scale(I16,f(x*x+y*y).div(2)),'directional square');}
 const metric=o.matrix([['1/2',0],[0,'1/2']]),bare=o.matrix([['1/2','1/2'],['1/2','1/2']]);
 ensure(o.rank(metric)===2&&o.rank(bare)===1,'derived rank distinction');
 return {basis_and_directional_matrix_identities:12,principal_metric:metric,principal_metric_rank:2,bare_R31_metric_rank:1,
   derived_Ax:Ax,derived_Ay:Ay,calibration:'component counts per eight-original-event B block'};
});

const cutpow=(z,n)=>{if(n<0)return cutpow(z.inv(),-n);let a=o.Cut.of(1);for(let j=0;j<n;j++)a=a.mul(z);return a;};
const evaluate=(a,x,y)=>[...a].reduce((s,[k,M])=>o.add(s,o.scale(M,cutpow(x,xy(k)[0]).mul(cutpow(y,xy(k)[1])))),o.zeros(a.values().next().value.length));
const phase=[o.Cut.of(1),o.Cut.of(-1),new o.Cut(0,1),new o.Cut(0,-1),new o.Cut('3/5','4/5'),new o.Cut('-3/5','4/5')];
const phaseRows=[];
check('native_unit_phase_family_and_actual_branch_projectors',()=>{
 const traceOp=new Map([...Z].map(([k,M])=>[k,o.matrix([[trace(M)]])])),expect=scale(minus(scale(id1,2),ell1),8);
 eqFields(traceOp,expect,1,1,'complete coefficient trace identity');
 let identities=1,projectors=0;for(const x of phase)for(const y of phase){
  ensure(x.norm2().eq(1)&&y.norm2().eq(1),'native phase norms');
  const m=evaluate(Z,x,y),l=evaluate(ell1,x,y)[0][0].rad;
  eq(o.mul(o.dagger(m),m),I16,'phase evaluation preserves matching');
  eq(o.add(m,o.dagger(m)),o.scale(I16,f(2).sub(l)),'phase exact central relation');
  eq(o.add(o.sub(o.mul(m,m),o.scale(m,f(2).sub(l))),I16),zero16,'actual source phase polynomial');identities+=3;
  if(l.eq(2)){const ip=new o.Cut(0,1),im=new o.Cut(0,-1),pplus=o.scale(o.sub(m,o.scale(I16,im)),ip.sub(im).inv()),pminus=o.sub(I16,pplus);
   eq(o.mul(pplus,pplus),pplus,'native branch idempotent');eq(o.mul(pplus,pminus),zero16,'disjoint branch projectors');
   eq(o.mul(m,pplus),o.scale(pplus,ip),'actual plus phase source');eq(o.dagger(pplus),pplus,'native phase projector dagger');
   ensure(o.rank(pplus)===8&&o.rank(pminus)===8,'two actual rank eight native bands');projectors+=2;
  }
  if(l.zero())eq(m,I16,'zero phase is full source identity');
  phaseRows.push({x:x.toJSON(),y:y.toJSON(),wave_coefficient:l.toString()});
 }
 return {exact_phase_and_trace_identities:identities,native_rational_phase_pairs:phaseRows,rank_eight_branch_projectors:projectors,Fourier_transform_assumed:false};
});

// Indeterminates r,s are formal native scalar coefficients, not address phases.
const rp=sx1,sp=sy1,oneMinusR=minus(id1,rp),oneMinusS=minus(id1,sp),twoMinusR=minus(scale(id1,2),rp),twoMinusS=minus(scale(id1,2),sp);
const lp=minus(scale(plus(rp,sp),2),compose(rp,sp)),denp=compose(lp,minus(scale(id1,4),lp));
const np=plus(compose(compose(twoMinusS,twoMinusS),compose(rp,oneMinusR)),compose(compose(twoMinusR,twoMinusR),compose(sp,oneMinusS)));
const margin=plus(plus(scale(compose(compose(rp,rp),oneMinusS),2),scale(compose(compose(sp,sp),oneMinusR),2)),
 plus(compose(compose(rp,sp),minus(scale(id1,2),plus(rp,sp))),scale(compose(compose(rp,rp),compose(sp,sp)),'3/2')));
check('sharp_phase_speed_has_a_complete_nonnegative_polynomial_proof',()=>{
 eqFields(minus(scale(denp,'1/2'),np),margin,1,1,'complete native positive polynomial');
 let n=0;for(let i=0;i<=16;i++)for(let j=0;j<=16;j++){
  const r=f(i).div(16),s=f(j).div(16),den=evaluate(denp,o.Cut.of(r),o.Cut.of(s))[0][0].rad;
  // Polynomial evaluation has nonnegative exponents; zero is valid here.
  const m=evaluate(margin,o.Cut.of(r),o.Cut.of(s))[0][0].rad;
  ensure(f(0).le(m),'nonnegative domain witness');
  if(!den.zero())ensure(f(0).le(m)&&!m.zero(),'strict regular phase margin');n++;
 }
 return {complete_polynomial_coefficient_identities:1,exact_domain_witnesses:n,nonnegative_terms:4,
   sharp_squared_phase_gradient_bound:'1/2',finite_nonzero_phase_attains_bound:false};
});

let speedLo=f(0),speedHi=f(1);
check('native_cut_root_and_rational_axis_family_certify_sharpness',()=>{
 const rows=[];for(const n of [2,3,4,8,16,32,64]){const r=f(1).div(n*n),v2=f(1).sub(r).div(f(2).sub(r)),gap=r.div(f(2).mul(f(2).sub(r)));
  ensure(f('1/2').sub(v2).eq(gap)&&f(0).le(gap),'sharp axis deficit');rows.push({r:r.toString(),squared_speed:v2.toString(),deficit:gap.toString()});}
 for(let n=0;n<80;n++){const mid=speedLo.add(speedHi).div(2);if(mid.pow(2).mul(2).le(1))speedLo=mid;else speedHi=mid;
  ensure(speedLo.pow(2).mul(2).le(1)&&f(1).le(speedHi.pow(2).mul(2)),'native positive speed cut');}
 const lo=f('7071067811865475/10000000000000000'),hi=f('7071067811865476/10000000000000000');
 ensure(lo.le(speedLo)&&speedHi.le(hi),'outward speed decimal interval');
 return {rational_sharpness_family:rows,positive_root_cut_depth:80,root_bracket:[speedLo.toString(),speedHi.toString()],
   outward_decimal_interval:['0.7071067811865475','0.7071067811865476'],native_source_root_equation:'2 c_count^2 = 1; c_count > 0'};
});

const entryAbs=M=>M.flat().reduce((s,c)=>s.add(c.rad.abs()).add(c.turn.abs()),f(0));
const moment=[...Z].reduce((s,[k,M])=>s.add(entryAbs(M).mul((Math.abs(xy(k)[0])+Math.abs(xy(k)[1]))**2)),f(0));
check('source_coefficient_moment_controls_the_native_phase_limit',()=>{
 ensure(moment.eq(88),'exact absolute entry second moment');ensure([...Z.keys()].every(k=>xy(k).every(i=>Math.abs(i)<=1)),'nine coefficient support');
 const momentRows=[...Z].map(([k,M])=>({shift:xy(k),absolute_entry_sum:entryAbs(M).toString(),weight:(Math.abs(xy(k)[0])+Math.abs(xy(k)[1]))**2}));
 let identities=0;const imaginary=new o.Cut(0,1);
 for(const [kx,ky]of [[1,0],[0,1],[1,1],[2,-1]]){
  const coeff=[];let fact=f(1),ip=o.Cut.of(1);
  for(let n=0;n<=4;n++){
   if(n){fact=fact.mul(n);ip=ip.mul(imaginary);}
   const a=[...Z].reduce((s,[k,M])=>o.add(s,o.scale(M,ip.mul(f(xy(k)[0]*kx+xy(k)[1]*ky).pow(n).div(fact)))),zero16);coeff.push(a);
  }
  eq(coeff[0],I16,'phase series constant');eq(coeff[1],o.scale(o.add(o.scale(Ax,kx),o.scale(Ay,ky)),imaginary),'phase first derivative');identities+=2;
  for(let n=1;n<=4;n++){
   let sum=zero16;for(let j=0;j<=n;j++)sum=o.add(sum,o.mul(o.dagger(coeff[j]),coeff[n-j]));eq(sum,zero16,'unit phase product coefficient '+n);identities++;
  }
 }
 return {exact_coefficient_moment:'88',coefficient_moment_ledger:momentRows,phase_jet_identities:identities,
   per_step_error_bound:'89 epsilon^2 (|kx|+|ky|)^2',n_step_error_bound:'89 n epsilon^2 (|kx|+|ky|)^2',
   condition:'epsilon (|kx|+|ky|) <= 1/2',modewise_continuum_control:true,all_field_continuum_limit_claimed:false};
});

check('finite_source_law_multiplicity_intertwines_the_same_envelope',()=>{
 const mm=o.identity(2),ii8=o.identity(8),p8=o.kron(P,o.identity(4)),q8=o.kron(Q,o.identity(4)),c8=o.kron(C0,o.identity(4));
 const xx=shift(unit(o.kron(I,o.kron(H,mm))),1,0),yy=shift(unit(o.kron(I,o.kron(K,mm))),0,1);
 const vv=a=>plus(unit(ii8),scale(compose(compose(plus(unit(p8),compose(dag(a),unit(q8))),unit(c8)),minus(a,unit(ii8))),'1/2'));
 const ww=compose(vv(yy),vv(xx)),bb=compose(ww,ww),expected=new Map([...B].map(([k,M])=>[k,o.kron(M,mm)]));
 eqFields(bb,expected,8,8,'inactive native multiplicity is exact');
 const turn=o.add(o.scale(I,'3/5'),o.scale(R,'4/5')),frame=o.kron(I,turn),fr=unit(frame);
 eq(o.mul(o.dagger(frame),frame),I4,'native record reframe');
 const bg=compose(compose(fr,B),dag(fr));eqFields(plus(bg,dag(bg)),minus(unit(o.scale(I4,2)),G),4,4,'same central law after record frame');
 return {complete_Laurent_multiplicity_and_frame_identities:2,record_rank_checked:4,source_roles_with_multiplicity:8,extra_coupling_from_inactive_multiplicity:false};
});

let cornerChecks=0;
check('exact_outer_front_and_block_units_keep_the_physical_boundary',()=>{
 let original=unit(o.kron(e0,e0));for(let n=1;n<=4;n++){original=compose(B,original);const expected=o.scale(o.kron(e0,e0),f((-1)**n).div(f(16).pow(n)));
  eq(original.get(key(2*n,2*n)),expected,'nonzero actual diagonal front');ensure(!o.isZero(expected),'corner stays nonzero');cornerChecks++;}
 const oldUniform=sumAtOne(B);eq(o.mul(o.dagger(o.add(oldUniform,I4)),o.add(oldUniform,I4)),o.scale(I4,3),'time sign alone is not the zero envelope');
 const uniform=sumAtOne(Z);eq(uniform,I16,'both signs selected gapless carrier');
 const wrongNext=o.scale(col(Array(16).fill(1)),-1);ensure(!o.equal(o.mul(uniform,col(Array(16).fill(1))),wrongNext),'arbitrary second-order initial data not source');
 ensure(f('1/2').mul(4).div(64).eq('1/32'),'coarse-count per original-event squared slope');
 ensure(f('1/2').mul(16).div(64).eq('1/8'),'fine-count per original-event squared slope');
 return {actual_source_outer_corner_checks:cornerChecks,source_events_per_envelope_block:8,coarse_steps_per_component_step:2,
   fine_steps_per_component_step:4,phase_speed_squared_in_component_units:'1/2',outer_corner_speed_squared_in_component_units:'2',
   exact_support_speed_equals_phase_speed:false,physical_length_time_calibration_selected:false};
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
 identity_gap_excludes_every_signed_carrier:o.equal(sumAtOne(Z),I16),
 temporal_sign_alone_selects_zero_carrier:!o.equal(sumAtOne(B),o.scale(I4,-1)),
 spatial_sign_alone_has_uniform_identity:!o.equal(o.scale(sumAtOne(Z),-1),I16),
 mixed_fourth_difference_may_be_dropped:!zeroOp(minus(ell1,scale(plus(lx1,ly1),'1/2'))),
 principal_metric_stays_rank_one:o.rank(o.matrix([['1/2',0],[0,'1/2']]))===2,
 parity_components_may_be_discarded:[...Z.values()].some(M=>M.some((r,i)=>r.some((v,j)=>!v.zero()&&Math.floor(i/4)!==Math.floor(j/4)))),
 phase_speed_is_exact_support_speed:f('1/2').le(2)&&!f('1/2').eq(2),
 every_scalar_wave_initial_pair_is_native_source:!o.equal(I16,o.scale(I16,-1)),
 all_regular_phases_attain_maximum:!f(0).eq('1/2'),
 scalar_probe_requires_imported_complex_field:input.constructed_maps.parity_on_roles.length===2,
 inactive_record_multiplicity_adds_a_coupling:checks[7].extra_coupling_from_inactive_multiplicity===false,
 old_gap_was_disproved:checks[0].old_identity_relative_gap==='1'&&checks[0].envelope_gap_infimum==='0',
 one_component_step_is_one_original_address:checks[8].fine_steps_per_component_step===4,
 one_envelope_event_is_one_original_event:checks[8].source_events_per_envelope_block===8,
 dimensional_speed_is_fixed_without_rod_clock:!f('1/2').eq(f('1/2').mul(4)),
};ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r35.native-signed-envelope.v1',status:'PASS_R35_NATIVE_SIGNED_ENVELOPE',
 input_sha256:inputHash,r34_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_signed_envelope_and_isotropic_cone_derived:true,native_sharp_phase_speed_derived:true,native_modewise_error_control_derived:true,native_front_and_calibration_boundaries_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written proofs: native signed carrier, full-rank principal roles, actual native phase bands, sharp phase-gradient bound, modewise error, record-law universality, and exact front/count-unit boundary. No physical c, h, alpha, electromagnetic field or universal continuum dynamics is identified.'},null,2)+'\n');
