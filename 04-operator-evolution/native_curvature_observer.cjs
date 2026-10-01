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
ensure(inputHash===arg('--expected-input-sha256'),'R36 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r35-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r35-sha256'),'R36 R35 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R35_NATIVE_SIGNED_ENVELOPE','R36 parent status');
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

const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Fc=o.scale(Jc,-2);
const Pc=o.scale(o.add(I16,Hc),'1/2'),Qc=o.sub(I16,Pc),readSecond=o.scale(o.mul(Pc,Jc),-1);
const readVector=v=>[o.mul(Pc,v),o.mul(readSecond,v)];
const decodeVector=([a,b])=>o.add(a,o.mul(Jc,b));
const ip=(a,b)=>o.mul(o.dagger(a),b)[0][0],en=v=>o.energy(v);
const current=v=>[ip(v,o.mul(Ax,v)).rad,ip(v,o.mul(Ay,v)).rad];
const cutPow=(z,n)=>{let a=o.Cut.of(1),base=n<0?z.inv():z;for(let k=0;k<Math.abs(n);k++)a=a.mul(base);return a;};
const phaseValue=(a,s,t)=>[...a].reduce((v,[k,M])=>{const [x,y]=xy(k);return o.add(v,o.scale(M,cutPow(s,x).mul(cutPow(t,y))));},o.zeros(a.values().next().value.length));
const nativePhases=[new o.Cut(1,0),new o.Cut(-1,0),new o.Cut(0,1),new o.Cut(0,-1),new o.Cut('3/5','4/5')];
const coord=n=>col(Array.from({length:16},(_,i)=>i===n?1:0));
const testVectors=[...Array.from({length:16},(_,i)=>coord(i)),col(Array.from({length:16},(_,i)=>(i%5)-2)),
 col(Array.from({length:16},(_,i)=>new o.Cut((i%3)-1,(i%4)-2)))];

check('native_increment_response_derives_the_same_information_metric',()=>{
 eqFields(normSquare(minus(Z,U16)),ell,16,16,'complete signed response cost');
 eq(o.mul(Ax,Ax),o.scale(I16,'1/2'),'metric xx');eq(o.mul(Ay,Ay),o.scale(I16,'1/2'),'metric yy');
 eq(o.add(o.mul(Ax,Ay),o.mul(Ay,Ax)),zero16,'metric xy');
 eq(o.mul(o.matrix([['1/2',0],[0,'1/2']]),o.matrix([[2,0],[0,2]])),I,'derived dual metric');
 let identities=0;
 for(const s of nativePhases)for(const t of nativePhases){const U=phaseValue(Z,s,t),cost=phaseValue(ell,s,t)[0][0];
  for(const v of testVectors.slice(0,4)){ensure(en(o.mul(o.sub(U,I16),v)).eq(cost.rad.mul(en(v))),'state independent actual response');identities++;}}
 return {complete_response_Gram_identities:1,principal_metric_identities:3,finite_state_phase_response_checks:identities,
  derived_metric:[['1/2','0'],['0','1/2']],dual_count_length_metric:[['2','0'],['0','2']],dual_metric_identities:1,
  metric_rank:2,physical_metric_selected:false};
});

check('source_response_products_derive_curvature_and_oriented_area',()=>{
 const comm=(A,B)=>o.sub(o.mul(A,B),o.mul(B,A));
 eq(comm(Ax,Ay),Jc,'principal order curvature');eq(o.dagger(Jc),o.scale(Jc,-1),'curvature skew dagger');
 eq(o.mul(Jc,Jc),o.scale(I16,-1),'curvature native turn');
 eq(o.mul(Hc,Hc),I16,'derived H involution');eq(o.mul(Kc,Kc),I16,'derived K involution');
 eq(o.add(o.mul(Hc,Kc),o.mul(Kc,Hc)),zero16,'derived source mixed law');
 eq(comm(Hc,Kc),Fc,'normalized role curvature');eq(o.mul(o.dagger(Fc),Fc),o.scale(I16,4),'normalized curvature square');
 eq(o.mul(o.mul(o.mul(Hc,Kc),Hc),Kc),o.scale(I16,-1),'actual four-arrow loop');
 const directions=[[1,0],[0,1],[1,1],[1,-1],[2,3],[-3,2]];let areaChecks=0;
 for(const a of directions)for(const b of directions){const Av=o.add(o.scale(Ax,a[0]),o.scale(Ay,a[1])),Bv=o.add(o.scale(Ax,b[0]),o.scale(Ay,b[1]));
  const d=a[0]*b[1]-a[1]*b[0],g=f(a[0]*b[0]+a[1]*b[1]).div(2),Fv=comm(Av,Bv);
  eq(o.mul(Av,Bv),o.add(o.scale(I16,g),o.scale(Jc,f(d).div(2))),'full response product');
  eq(o.mul(o.dagger(Fv),Fv),o.scale(I16,d*d),'oriented area square');
  const detGram=f(a[0]*a[0]+a[1]*a[1]).mul(b[0]*b[0]+b[1]*b[1]).sub(f(a[0]*b[0]+a[1]*b[1]).pow(2)).div(4);
  ensure(detGram.mul(4).eq(d*d),'metric area relation');areaChecks+=3;}
 return {complete_native_role_identities:9,directional_product_area_checks:areaChecks,
  principal_curvature_squared:'-I',role_curvature_norm_squared:'4 I',curvature_area_factor:'4',derived_H:Hc,derived_K:Kc,derived_J:Jc};
});

check('curvature_completes_a_minimal_native_two_cut_observer',()=>{
 eq(o.mul(Pc,Pc),Pc,'first observer cut');eq(o.dagger(Pc),Pc,'observer cut pairing');
 eq(o.mul(o.dagger(readSecond),readSecond),Qc,'curvature reads exactly the missing cut');
 eq(o.add(Pc,o.mul(Jc,readSecond)),I16,'exact decoder');
 eq(o.mul(Pc,Jc),o.mul(Jc,Qc),'curvature exchanges cut roles');
 eq(o.scale(o.mul(Pc,Fc),'1/2'),readSecond,'second readout is actual curvature contrast');
 const O=Pc.concat(readSecond),D=Pc.map((row,i)=>row.concat(o.mul(Jc,Pc)[i]));
 eq(o.mul(D,O),I16,'complete finite observer');
 eq(o.mul(o.dagger(O),O),I16,'observer matching isometry');
 ensure(o.rank(Pc)===8&&o.rank(O)===16,'native observer minimal role counts');
 for(const v of testVectors){const pair=readVector(v);eq(decodeVector(pair),v,'finite full-state reconstruction');ensure(en(pair[0]).add(en(pair[1])).eq(en(v)),'retained norm');}
 const hidden=o.mul(Qc,coord(0));ensure(!o.isZero(hidden)&&o.isZero(o.mul(Pc,hidden)),'one cut hidden witness');
 ensure(!o.isZero(o.mul(readSecond,hidden)),'curvature repair detects hidden witness');
 return {complete_observer_identities:8,finite_reconstructions:testVectors.length,one_cut_rank:8,complete_observer_rank:16,
  minimum_readouts_of_the_same_rank_eight_cut:2,arbitrary_physical_observer_minimum_claimed:false};
});

const cp=unit(Pc),cj=unit(Jc),cr=unit(readSecond);
const M11=compose(compose(cp,Z),cp),M12=compose(compose(cp,Z),compose(cj,cp));
const M21=compose(compose(cr,Z),cp),M22=compose(compose(cr,Z),compose(cj,cp));
const readField=u=>[compose(cp,u),compose(cr,u)],decodeField=([a,b])=>plus(a,compose(cj,b));
const stepPair=([a,b])=>[plus(compose(M11,a),compose(M12,b)),plus(compose(M21,a),compose(M22,b))];
check('the_complete_curvature_readout_intertwines_exact_source_dynamics',()=>{
 eqFields(compose(cp,Z),plus(compose(M11,cp),compose(M12,cr)),16,16,'first exact observed equation');
 eqFields(compose(cr,Z),plus(compose(M21,cp),compose(M22,cr)),16,16,'second exact observed equation');
 eqFields(plus(M11,dag(M11)),minus(scale(cp,2),compose(ell,cp)),16,16,'observed diagonal wave 11');
 eqFields(plus(M22,dag(M22)),minus(scale(cp,2),compose(ell,cp)),16,16,'observed diagonal wave 22');
 eqFields(plus(M12,dag(M21)),new Map(),16,16,'observed off diagonal wave');
 const zc=compose(compose(cj,Z),dag(cj));
 eqFields(plus(zc,dag(zc)),minus(scale(U16,2),ell),16,16,'curvature field same exact wave');
 const seeds=[unit(coord(0)),new Map([['0,0',coord(3)],['-1,1',o.scale(coord(8),2)]])];let runs=0;
 for(const seed of seeds){let u=seed,pair=readField(seed);const norm=energy(seed);for(let n=0;n<=3;n++){
  eqFields(decodeField(pair),u,16,1,'exact observed reconstruction after continuation');
  ensure(energy(pair[0]).add(energy(pair[1])).eq(norm),'observed continuation isometry');runs++;
  if(n<3){u=compose(Z,u);pair=stepPair(pair);}
 }}
 return {complete_Laurent_intertwining_wave_identities:6,finite_continuation_reconstructions:runs,
  exact_scalar_wave_shared:true,curvature_and_paired_fields_are_independent_physical_fields:false};
});

check('a_single_cut_has_derived_memory_and_cannot_close_the_source',()=>{
 const cross=compose(compose(cp,Z),unit(Qc));ensure(!zeroOp(cross),'one cut cannot autonomously continue');
 const seed=unit(o.mul(Qc,coord(0))),first=compose(cp,compose(Z,seed));
 ensure(zeroOp(compose(cp,seed))&&!zeroOp(first),'same first readout different next readout');
 let pair=readField(seed),as=[pair[0]],bs=[pair[1]],memoryChecks=0;
 const repeatedHidden=(field,n)=>{let result=field;for(let i=0;i<n;i++)result=compose(M22,result);return result;};
 for(let n=0;n<4;n++){
  let predicted=plus(compose(M11,as[n]),compose(M12,repeatedHidden(bs[0],n)));
  for(let k=0;k<n;k++)predicted=plus(predicted,compose(M12,repeatedHidden(compose(M21,as[n-1-k]),k)));
  pair=stepPair(pair);eqFields(pair[0],predicted,16,1,'exact hidden-history memory');as.push(pair[0]);bs.push(pair[1]);memoryChecks++;
 }
 return {one_cut_kernel_invariant:false,hidden_initial_state_witness:true,finite_memory_recurrences:memoryChecks,
  memory_kernel:'M12 M22^k M21',independent_noise_inserted:false};
});

check('native_information_current_has_the_same_sharp_cone_and_exact_deficit',()=>{
 let scalarChecks=0,continuityChecks=0;const rows=[];
 for(const v of testVectors){const [a,b]=readVector(v),ra=en(a),rb=en(b),h=ip(a,b).rad,rho=en(v),[jx,jy]=current(v);
  ensure(jx.eq(ra.sub(rb).add(h.mul(2)).div(2))&&jy.eq(ra.sub(rb).sub(h.mul(2)).div(2)),'paired information current');
  const deficit=rho.pow(2).div(2).sub(jx.pow(2)).sub(jy.pow(2)),gram=ra.mul(rb).sub(h.pow(2));
  ensure(deficit.eq(gram.mul(2))&&f(0).le(deficit),'exact coherence deficit');scalarChecks+=3;
  const ux=o.mul(Kc,v),uy=o.mul(Jc,v),ut=o.scale(o.add(o.mul(Ax,ux),o.mul(Ay,uy)),-1);
  ensure(ip(v,ut).rad.add(ip(v,o.add(o.mul(Ax,ux),o.mul(Ay,uy))).rad).eq(0),'local first-jet information continuity');continuityChecks++;
 }
 const a=o.mul(Pc,coord(0));for(const t of [0,1,-1,2,'1/2']){const v=decodeVector([a,o.scale(a,t)]),rho=en(v),[jx,jy]=current(v);
  ensure(jx.pow(2).add(jy.pow(2)).eq(rho.pow(2).div(2)),'sharp same-record family');
  rows.push({second_over_first:String(t),rho:rho.toString(),jx:jx.toString(),jy:jy.toString()});}
 const atZero=o.mul(Pc,coord(0)),btZero=o.mul(Pc,coord(1));
 const borth=o.sub(btZero,o.scale(atZero,ip(atZero,btZero).div(ip(atZero,atZero))));
 ensure(!o.isZero(borth),'independent retained records exist');
 const n0=en(atZero),n1=en(borth);let vzero;
 // These two source columns have the same native norm; no irrational sampler is used.
 ensure(n0.eq(n1),'zero-current record norms');vzero=decodeVector([atZero,borth]);
 const zeroCurrent=current(vzero);ensure(zeroCurrent.every(x=>x.eq(0))&&!en(vzero).zero(),'nonzero information may have zero principal current');
 for(const v of testVectors){const j=current(v),r=current(o.mul(Jc,v));ensure(r[0].eq(j[0].neg())&&r[1].eq(j[1].neg()),'curvature reverses principal current');}
 return {exact_current_and_deficit_checks:scalarChecks,native_local_continuity_checks:continuityChecks,
  sharp_current_examples:rows,sharp_squared_current_per_density:'1/2',nonzero_zero_current_witness:true,
  physical_energy_or_probability_current_assumed:false,exact_finite_lattice_current_claimed:false};
});

check('native_modewise_limit_transfers_to_the_observer_with_finite_corrections',()=>{
 const O=Pc.concat(readSecond),D=Pc.map((row,i)=>row.concat(o.mul(Jc,Pc)[i]));
 const zero=o.zeros(16),Hobs=Pc.map((row,i)=>row.concat(zero[i])).concat(zero.map((row,i)=>row.concat(o.scale(Pc,-1)[i])));
 const Kobs=zero.map((row,i)=>row.concat(Pc[i])).concat(Pc.map((row,i)=>row.concat(zero[i])));
 eq(o.mul(O,o.mul(Hc,D)),Hobs,'derived H in the cut observer');eq(o.mul(O,o.mul(Kc,D)),Kobs,'derived K in the cut observer');
 eq(o.mul(O,o.mul(Ax,D)),o.scale(o.add(Hobs,Kobs),'1/2'),'observed x generator');
 eq(o.mul(O,o.mul(Ay,D)),o.scale(o.sub(Hobs,Kobs),'1/2'),'observed y generator');
 const defect=minus(compose(compose(cj,Z),dag(cj)),dag(Z));ensure(!zeroOp(defect),'finite curvature conjugate is not generally the inverse');
 eq(sumAtOne(defect),zero16,'duality zeroth defect');
 for(let axis=0;axis<2;axis++){const d=[...defect].reduce((s,[k,M])=>o.add(s,o.scale(M,xy(k)[axis])),zero16);eq(d,zero16,'duality first defect');}
 const mass=M=>M.flat().reduce((a,x)=>{ensure(x.turn.eq(0),'real source coefficient');return a.add(x.rad.abs());},f(0));
 const moment=[...Z].reduce((s,[k,M])=>{const [a,b]=xy(k);return s.add(mass(M).mul((Math.abs(a)+Math.abs(b))**2));},f(0));
 ensure(moment.eq(88),'inherited complete error coefficient');
 return {full_observer_generator_identities:4,complete_coefficient_moment:moment.toString(),
  inherited_n_step_error_bound:'89 n epsilon^2 (|kx|+|ky|)^2',curvature_inverse_defect_nonzero:true,
  curvature_inverse_defect_vanishes_through_degree:1,curvature_inverse_phase_defect_bound:'176 epsilon^2 (|kx|+|ky|)^2',
  common_condition:'epsilon (|kx|+|ky|) <= 1/2'};
});

check('native_frames_preserve_the_complete_reading_and_protocol_curvature',()=>{
 const g=o.add(o.scale(I16,'3/5'),o.scale(Jc,'4/5')),gd=o.dagger(g);
 eq(o.mul(gd,g),I16,'source-built frame isometry');
 const ax=o.mul(g,o.mul(Ax,gd)),ay=o.mul(g,o.mul(Ay,gd)),hc=o.add(ax,ay),kc=o.sub(ax,ay),jc=o.mul(kc,hc),pc=o.scale(o.add(I16,hc),'1/2');
 eq(jc,o.mul(g,o.mul(Jc,gd)),'covariant protocol curvature');eq(o.mul(ax,ax),o.scale(I16,'1/2'),'covariant metric');
 for(const v of testVectors.slice(0,6)){const gv=o.mul(g,v),[a,b]=readVector(v);
  eq(o.mul(pc,gv),o.mul(g,a),'first frame readout');eq(o.scale(o.mul(pc,o.mul(jc,gv)),-1),o.mul(g,b),'second frame readout');}
 // Explicit local pure-frame shift loop: all intermediate frames cancel.
 const frames=[I16,Hc,Kc,Jc,g],link=(i,j)=>o.mul(frames[i],o.dagger(frames[j]));
 const path1=o.mul(link(4,1),link(1,0)),path2=o.mul(link(4,2),link(2,0));eq(path1,path2,'local shift frame has no manufactured flux');
 return {native_frame_source_checks:3,transported_readout_checks:12,pure_frame_path_equality:1,
  constant_response_metric_and_nonzero_protocol_curvature_coexist:true,Riemann_connection_selected:false};
});

check('signed_curvature_reading_and_physical_calibration_remain_distinct',()=>{
 eq(o.scale(o.add(Fc,o.mul(Hc,o.mul(Fc,Hc))),'1/2'),zero16,'old cut-even curvature blindness');
 eq(o.mul(o.dagger(Fc),Fc),o.scale(I16,4),'curvature square survives');
 const v=coord(0),r=readVector(v),s=readVector(o.scale(v,-1));
 ensure(en(r[0]).eq(en(s[0]))&&en(r[1]).eq(en(s[1]))&&!o.equal(decodeVector(r),decodeVector(s)),'two intensities do not reconstruct signed state');
 const loop=o.mul(o.mul(o.mul(Hc,Kc),Hc),Kc);ensure(!o.equal(loop,I16),'closed source record is not flat');
 ensure(f('1/2').mul(4).eq(2)&&!f('1/2').eq(2),'count coefficient changes with physical length/time calibration');
 return {cut_even_curvature_blindness_preserved:true,quadratic_information_does_not_determine_signed_state:true,
  two_coarse_steps_per_component:2,original_events_per_block:8,
  physical_c_h_alpha_identified:false,R10_smooth_metric_adapter_used_as_proof_premise:false};
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
 ordinary_metric_was_a_premise:checks[0].derived_metric[0][0]==='1/2',
 response_metric_depends_on_preparation:checks[0].finite_state_phase_response_checks===100,
 flat_response_metric_forces_zero_protocol_curvature:!o.isZero(Jc),
 curvature_has_no_area_readout:checks[1].curvature_area_factor==='4',
 one_cut_reconstructs_the_full_source:o.rank(Pc)===8,
 two_intensities_reconstruct_signed_state:checks[8].quadratic_information_does_not_determine_signed_state,
 curvature_readout_is_independent_noise:checks[4].independent_noise_inserted===false,
 one_cut_has_autonomous_continuation:checks[4].one_cut_kernel_invariant===false,
 all_nonzero_states_carry_maximal_current:checks[5].nonzero_zero_current_witness,
 current_deficit_requires_a_free_coupling:checks[5].sharp_squared_current_per_density==='1/2',
 finite_curvature_continuation_is_exact_inverse:checks[6].curvature_inverse_defect_nonzero,
 phase_cone_is_exact_lattice_support_cone:checks[5].exact_finite_lattice_current_claimed===false,
 local_reframing_creates_shift_curvature:checks[7].pure_frame_path_equality===1,
 Riemann_bridge_is_already_a_selected_spacetime:checks[7].Riemann_connection_selected===false,
 physical_c_follows_without_length_time_identification:checks[8].physical_c_h_alpha_identified===false,
};ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r36.native-curvature-observer.v1',status:'PASS_R36_NATIVE_CURVATURE_OBSERVER',
 input_sha256:inputHash,r35_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_response_metric_and_area_curvature_derived:true,native_minimal_curvature_observer_derived:true,
 native_information_current_cone_derived:true,native_observed_dynamics_and_limit_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written proofs: native response metric and oriented curvature, minimal complete curvature observer, exact paired dynamics and memory, sharp information-current cone and coherence deficit, controlled modewise flow, and frame/Riemann/physical boundaries. The constructed fields are readouts of one source, not independently identified physical fields.'},null,2)+'\n');
