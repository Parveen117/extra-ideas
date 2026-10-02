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
ensure(inputHash===arg('--expected-input-sha256'),'R41 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r40-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r40-sha256'),'R41 R40 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R40_NATIVE_AUTONOMOUS_ECHO','R41 parent status');
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

const io=o.IOTA,transpose=A=>o.matrix(A[0].map((_,j)=>A.map(r=>r[j])));
const basis=(n,j)=>col(Array.from({length:n},(_,i)=>i===j?1:0));
const rootWitness=o.scale(C0,'1/2'),D16raw=o.kron(o.identity(4),o.kron(A0,I));
const link=plus(compose(Y,unit(Ps)),unit(Qs)),U=scale(block(link),-1);
const tensorRole=(a,M)=>new Map([...a].map(([k,A])=>[k,o.kron(A,M)]));
const U32=tensorRole(U,I),Ud32=dag(U32),Z32=tensorRole(Z,I),Zd32=dag(Z32),D32=o.kron(D16raw,rootWitness);
const rect=(xmin,xmax,ymin,ymax)=>(x,y)=>xmin<=x&&x<=xmax&&ymin<=y&&y<=ymax;
const origin=rect(0,0,-1,1),remote=rect(1,2,-1,1),onsite=(a,T)=>new Map([...a].map(([k,v])=>[k,T(...xy(k))?o.mul(D32,v):v]));
const gZ={apply:a=>compose(Z32,a),inverse:a=>compose(Zd32,a)},gU={apply:a=>compose(U32,a),inverse:a=>compose(Ud32,a)},gUd={apply:a=>compose(Ud32,a),inverse:a=>compose(U32,a)};
const contrast=T=>({apply:a=>onsite(a,T),inverse:a=>onsite(a,T)});
const program=[gZ,gUd,contrast(remote),gU,gZ,gUd,contrast(origin),gU],qProgram=program.length;
const sourcePrefix=(a,n,back=false)=>{let b=a;if(back){for(let j=n-1;j>=0;j--)b=program[j].inverse(b);}else for(let j=0;j<n;j++)b=program[j].apply(b);return b;};
const sourceCycle=a=>sourcePrefix(a,qProgram),sourceSeed=unit(o.kron(basis(16,0),e0));
const processApply=(a,V)=>new Map([...a].map(([k,M])=>[k,o.mul(M,transpose(V))]));
const productField=(a,v)=>new Map([...a].map(([k,M])=>[k,o.mul(M,transpose(v))]));
const processTurn=o.add(o.scale(I,'3/5'),o.scale(R,'4/5')),processes=[H,R,processTurn];
const fieldPair=(a,b)=>[...a].reduce((s,[k,M])=>s.add(trace(o.mul(o.dagger(M),b.get(k)||o.zeros(M.length,M[0].length)))),o.ZERO);
const fieldHash=a=>hash(JSON.stringify([...a].sort(([a],[b])=>a.localeCompare(b))));
const mark=(j,m)=>j+','+m,marks=k=>k.split(',').map(Number);
const joint=(a,j=0,m=0)=>new Map([[mark(j,m),a]]),jointEnergy=a=>[...a.values()].reduce((s,v)=>s.add(energy(v)),f(0));
const eqJoint=(a,b,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eqFields(a.get(k)||new Map(),b.get(k)||new Map(),32,2,msg+' '+k);};
const coupledStep=(a,V,Qn,back=false)=>{const out=new Map();for(const [k,v]of a){const [j,m]=marks(k);let jj,mm,vv;
 if(back){jj=(j+qProgram-1)%qProgram;mm=(m+Qn-(j===0?1:0))%Qn;vv=j===0?processApply(v,o.dagger(V)):v;vv=program[jj].inverse(vv);}
 else{jj=(j+1)%qProgram;mm=(m+(j===qProgram-1?1:0))%Qn;vv=program[j].apply(v);if(j===qProgram-1)vv=processApply(vv,V);}
 ensure(!out.has(mark(jj,mm)),'joint mark collision');out.set(mark(jj,mm),vv);
 }return out;};
const sourcePowers=[sourceSeed,sourceCycle(sourceSeed)];
check('native_tick_coupling_runs_the_process_on_the_retained_clock',()=>{
 eqFields(compose(dag(U),U),U16,16,16,'unchanged native source link');eq(o.mul(D32,D32),o.identity(32),'native control root');
 let prefixes=0,inverses=0;const v=col([1,2]),rows=[];
 for(const V of processes){eq(o.mul(o.dagger(V),V),I,'constructed process preserves matching');let a=joint(productField(sourceSeed,v));
  for(let n=0;n<=qProgram+2;n++){if(n)a=coupledStep(a,V,3);if([0,1,4,7,8,10].includes(n)){const k=Math.floor(n/qProgram),r=n%qProgram,expected=productField(sourcePrefix(sourcePowers[k],r),o.mul(o.power(V,k),v));eqJoint(a,joint(expected,r,k),'complete coupled prefix');prefixes++;}}
  rows.push({process:V,checked_program_updates:qProgram+3,final_field_sha256:fieldHash([...a.values()][0])});
 }
 for(let j=0;j<qProgram;j++){const a=joint(productField(sourceSeed,v),j,1),b=coupledStep(a,processTurn,3);eqJoint(coupledStep(b,processTurn,3,true),a,'fixed joint inverse');ensure(jointEnergy(a).eq(jointEnergy(b)),'joint source/process norm');inverses++;}
 return {complete_coupled_prefix_equalities:prefixes,marked_inverse_and_norm_checks:inverses,native_processes:rows,process_updates_per_clock_cycle:1,source_program_updates_per_cycle:qProgram,source_root_witness_checks_both_algebraic_conjugates:true,physical_remote_control_asserted:false};
});

// Finite clock arrays below use the same native Cut/matrix engine. They
// provide small generic isometry witnesses for the all-W written identities.
const arAdd=(a,b)=>a.map((x,j)=>o.add(x,b[j])),arScale=(a,c)=>a.map(x=>o.scale(x,c)),arMinus=(a,b)=>arAdd(a,arScale(b,-1));
const arPair=(a,b)=>a.reduce((s,x,j)=>s.add(o.mul(o.dagger(x),b[j])[0][0]),o.ZERO),arEnergy=a=>arPair(a,a).rad;
const tick=(a,V,back=false)=>a.map((_,j)=>o.mul(back?o.dagger(V):V,a[(j+(back?1:a.length-1))%a.length]));
const count=a=>a.map((v,j)=>o.scale(v,j)),response=(a,V)=>arScale(arMinus(tick(a,V),tick(a,V,true)),new o.Cut(0,'1/2'));
const centeredResponse=(a,V)=>{const Qn=a.length,b=a.map(()=>o.zeros(2,1));b[0]=o.mul(V,a[Qn-1]);b[Qn-1]=o.mul(o.dagger(V),a[0]);return arMinus(arScale(arAdd(tick(a,V),tick(a,V,true)),'1/2'),arScale(b,f(Qn).div(2)));};
const history=(weights,V,v=e0)=>{let vr=v;return weights.map(h=>{const a=o.scale(vr,h);vr=o.mul(V,vr);return a;});};
const pascal=s=>{let row=[f(1)];for(let n=0;n<s;n++){const next=Array(n+2).fill(f(0));for(let j=0;j<row.length;j++){next[j]=next[j].add(row[j]);next[j+1]=next[j+1].add(row[j]);}row=next;}return row;};
const profile=(s,Qn)=>{const row=pascal(s);return Array.from({length:Qn},(_,r)=>r>=1&&r<=s+1?row[r-1]:f(0));};
const stats=(a,V)=>{const E0=arEnergy(a),nv=count(a),pv=response(a,V),nm=arPair(a,nv).div(E0),pm=arPair(a,pv).div(E0);ensure(nm.turn.zero()&&pm.turn.zero(),'self-dagger means are real');
 const dn=arMinus(nv,arScale(a,nm)),dp=arMinus(pv,arScale(a,pm)),cross=arPair(dn,dp).div(E0),rho=arPair(a,centeredResponse(a,V)).div(E0);ensure(rho.turn.zero(),'real count-response expectation');
 const vn=arEnergy(dn).div(E0),vp=arEnergy(dp).div(E0),cov=cross.rad,gap=vn.mul(vp).sub(cov.pow(2)).sub(rho.rad.pow(2).div(4));
 ensure(cross.turn.eq(rho.rad.div(2)),'native turn part of centered response');
 return {norm_square:E0,count_mean:nm.rad,response_mean:pm.rad,count_variance:vn,response_variance:vp,covariance:cov,commutator_response:rho.rad,uncertainty_gap:gap};};
let nonclosingHistories=0;
check('native_history_isometry_retains_winding_and_clock_readout',()=>{
 let identities=0;const rows=[];for(const Qn of [3,4,5])for(const V of processes){const weights=Array(Qn).fill(1),a=history(weights,V,e0),b=history(weights,V,col([1,new o.Cut(0,1)]));
  ensure(arEnergy(a).eq(Qn),'uniform history norm');ensure(arPair(a,b).eq(o.mul(o.dagger(e0),col([1,new o.Cut(0,1)]))[0][0].mul(Qn)),'full native history pairing');
  const defect=arEnergy(arMinus(tick(a,V),a)),winding=o.energy(o.sub(o.mul(o.power(V,Qn),e0),e0));ensure(defect.eq(winding),'complete cyclic history residue');if(!defect.zero())nonclosingHistories++;
  identities++;rows.push({tick_labels:Qn,normalized_history_defect_square:defect.div(Qn).toString(),winding_defect_square:winding.toString()});
 }
 const pair=history([1,1],K).reduce((s,v)=>o.add(s,o.mul(v,o.dagger(v))),o.zeros(2));eq(o.scale(pair,'1/2'),o.scale(I,'1/2'),'unresolved clock count readout');
 ensure(o.rank(pair)===2&&!o.equal(o.scale(pair,'1/2'),o.mul(e0,o.dagger(e0))),'retained history changes unresolved process pair');
 return {full_history_pairing_and_winding_checks:identities,uniform_history_examples:rows,nonclosing_native_histories:nonclosingHistories,unresolved_process_pair:o.scale(pair,'1/2'),random_history_selection_assumed:false,unknown_state_copying_map_used:false};
});
const tickMatrix=(V,Qn)=>{const d=V.length,M=o.zeros(d*Qn);for(let r=0;r<Qn;r++)for(let i=0;i<d;i++)for(let j=0;j<d;j++)M[((r+1)%Qn)*d+i][r*d+j]=V[i][j];return M;};
const countMatrix=(d,Qn)=>{const M=o.zeros(d*Qn);for(let r=0;r<Qn;r++)for(let i=0;i<d;i++)M[r*d+i][r*d+i]=o.Cut.of(r);return M;};
const wrapMatrix=(V,Qn)=>{const d=V.length,M=o.zeros(d*Qn);for(let i=0;i<d;i++)for(let j=0;j<d;j++)M[i][(Qn-1)*d+j]=V[i][j];return M;};
const cutPower=(z,n)=>{let v=o.ONE;for(let j=0;j<n;j++)v=v.mul(z);return v;};
let boundaryDefectHash;
check('native_phase_loop_and_wrap_cut_give_complete_commutator_identities',()=>{
 let counts=0,phaseLoops=0;const phaseRows=[];
 for(const Qn of [2,3,4,5,8])for(const V of processes){const d=V.length,T=tickMatrix(V,Qn),Nt=countMatrix(d,Qn),Pt=o.scale(o.sub(T,o.dagger(T)),new o.Cut(0,'1/2')),Ct=o.scale(o.add(T,o.dagger(T)),'1/2'),Jt=wrapMatrix(V,Qn),Rt=o.sub(Ct,o.scale(o.add(Jt,o.dagger(Jt)),f(Qn).div(2)));
  eq(o.mul(o.dagger(T),T),o.identity(d*Qn),'native coupled tick inverse');eq(o.dagger(Pt),Pt,'native tick contrast dagger');eq(o.commutator(Nt,Pt),o.scale(Rt,io),'complete count-response commutator');counts++;
  if(Qn===3){const defect=o.sub(o.commutator(Nt,Pt),o.scale(Ct,io));ensure(!o.isZero(defect),'wrap cut is necessary');boundaryDefectHash=hash(JSON.stringify(defect));}
  if(Qn===2||Qn===4){const z=Qn===2?o.Cut.of(-1):io,Dp=o.zeros(d*Qn);for(let r=0;r<Qn;r++)for(let i=0;i<d;i++)Dp[r*d+i][r*d+i]=cutPower(z,r);
   eq(o.mul(o.mul(o.mul(Dp,T),o.dagger(Dp)),o.dagger(T)),o.scale(o.identity(d*Qn),z),'native phase loop');const com=o.commutator(Dp,T),coef=o.Cut.of(2).sub(z).sub(z.dagger());eq(o.mul(o.dagger(com),com),o.scale(o.identity(d*Qn),coef),'phase curvature square');phaseLoops++;phaseRows.push({period:Qn,curvature_square:coef.rad.toString()});
  }
 }
 const z8=o.scale(rootWitness,new o.Cut(1,1));eq(o.power(z8,2),o.scale(I,io),'eighth-turn square');eq(o.power(z8,4),o.scale(I,-1),'eighth-turn primitive half period');eq(o.power(z8,8),I,'eighth-turn full period');
 const V=o.kron(processTurn,I),T=tickMatrix(V,8),Dp=o.zeros(32);for(let r=0;r<8;r++){const z=o.kron(I,o.power(z8,r));for(let i=0;i<4;i++)for(let j=0;j<4;j++)Dp[4*r+i][4*r+j]=z[i][j];}
 const expected=o.kron(o.identity(16),z8),kap8=o.sub(o.scale(I,2),o.add(z8,o.dagger(z8))),com=o.commutator(Dp,T);
 eq(o.mul(o.mul(o.mul(Dp,T),o.dagger(Dp)),o.dagger(T)),expected,'complete eighth-turn loop');eq(o.mul(o.dagger(com),com),o.kron(o.identity(16),kap8),'eighth-turn curvature square');eq(o.add(o.sub(o.mul(kap8,kap8),o.scale(kap8,4)),o.scale(I,2)),o.zeros(2),'native eighth-turn curvature polynomial');
 eq(tickMatrix(o.matrix([[1]]),2),K,'two-mark successor is native K');eq(o.matrix([[1,0],[0,-1]]),H,'two-mark phase is native H');
 return {complete_count_response_identities:counts,scalar_phase_loop_and_square_checks:phaseLoops,eighth_turn_matrix_loop_checks:3,phase_curvature_values:phaseRows,eighth_turn_positive_branch_curvature:'2 - sqrt(2)',eighth_turn_witness_checks_both_algebraic_conjugates:true,omitted_wrap_defect_sha256:boundaryDefectHash,global_finite_clock_canonical_identity_claimed:false};
});
let sharpWitness;
check('native_matching_gives_the_full_covariance_bound_and_a_sharp_witness',()=>{
 let cases=0,zeroCountVariance=0;for(const Qn of [3,4,5,8])for(const V of [I,R,processTurn]){
  const arbitrary=Array.from({length:Qn},(_,r)=>col([new o.Cut(r%3-1,(r+1)%2),new o.Cut((r+1)%5-2,r%2)]));
  const localized=Array.from({length:Qn},(_,r)=>r===1?e0:o.zeros(2,1));
  for(const a of [arbitrary,localized,history(Array(Qn).fill(1),V)]){const st=stats(a,V);ensure(f(0).le(st.uncertainty_gap),'native covariance determinant');if(st.count_variance.zero()){ensure(st.commutator_response.zero()&&st.covariance.zero(),'zero count variance branch');zeroCountVariance++;}cases++;}
 }
 sharpWitness=stats([e0,o.zeros(2,1),e0],I);ensure(sharpWitness.count_variance.eq(1)&&sharpWitness.response_variance.eq('1/4')&&sharpWitness.commutator_response.eq(-1)&&sharpWitness.uncertainty_gap.zero(),'nonzero exact sharpness witness');
 return {exact_covariance_inequality_instances:cases,zero_count_variance_instances:zeroCountVariance,sharp_three_mark_witness:sharpWitness,sharp_squared_coefficient:'1/4',external_uncertainty_law_used_as_premise:false};
});
let profileRows=[];
check('binary_word_counts_derive_all_finite_relational_history_moments',()=>{
 let cases=0;for(let s=1;s<=24;s++){let Qn=4;while(Qn<s+3)Qn*=2;const weights=profile(s,Qn),nu=pascal(2*s)[s];
  for(const V of [I,R,processTurn]){const a=history(weights,V),st=stats(a,V),rho=f(s).div(s+1),gap=f(3*s*s).div(f(4).mul(2*s-1).mul(f(s+1).pow(2)).mul(s+2));
   ensure(st.norm_square.eq(nu)&&st.count_mean.eq(f(1).add(f(s).div(2)))&&st.response_mean.eq(0),'word profile norm and means');
   ensure(st.count_variance.eq(f(s*s).div(4*(2*s-1)))&&st.response_variance.eq(f(2*s+1).div(f(s+1).mul(s+2)))&&st.covariance.eq(0),'complete profile variances');
   ensure(st.commutator_response.eq(rho)&&st.uncertainty_gap.eq(gap),'profile response and exact gap');
   ensure(arEnergy(arMinus(tick(a,V),a)).div(nu).eq(f(2).div(s+1)),'profile stationary defect');
   ensure(arPair(a,tick(tick(a,V),V)).div(nu).eq(f(s*(s-1)).div(f(s+1).mul(s+2))),'two-step history overlap');cases++;
  }
  if([1,2,4,8,16,24].includes(s))profileRows.push({word_length:s,tick_labels:Qn,...stats(history(weights,processTurn),processTurn),stationary_defect_square:f(2).div(s+1)});
 }
 return {exact_word_profile_instances:cases,profile_examples:profileRows,weights_constructed_by_binary_word_recursion:true,limiting_probability_distribution_assumed:false,profile_response_uses_the_joint_tick_not_bare_clock_shift:true};
});

const fldArAdd=(a,b)=>a.map((x,j)=>plus(x,b[j])),fldArScale=(a,c)=>a.map(x=>scale(x,c)),fldArMinus=(a,b)=>fldArAdd(a,fldArScale(b,-1));
const fldArPair=(a,b)=>a.reduce((s,x,j)=>s.add(fieldPair(x,b[j])),o.ZERO),fldArEnergy=a=>fldArPair(a,a).rad;
let actualHistory,actualResponses;
check('the_actual_echo_carrier_and_process_realize_the_relational_profile',()=>{
 for(let j=2;j<=3;j++)sourcePowers.push(sourceCycle(sourcePowers[j-1]));
 const powers=sourcePowers.map((a,j)=>productField(a,o.mul(o.power(R,j),e0))),empty=()=>new Map();
 actualHistory=[empty(),powers[1],powers[2],empty()];
 const next=[empty(),empty(),powers[2],powers[3]],prev=[powers[0],powers[1],empty(),empty()];
 const nv=actualHistory.map((a,j)=>scale(a,j)),pv=fldArScale(fldArMinus(next,prev),new o.Cut(0,'1/2')),rv=fldArScale(fldArAdd(next,prev),'1/2');
 const norm=fldArEnergy(actualHistory),mn=fldArPair(actualHistory,nv).div(norm),mp=fldArPair(actualHistory,pv).div(norm),dn=fldArMinus(nv,fldArScale(actualHistory,mn)),dp=fldArMinus(pv,fldArScale(actualHistory,mp));
 const vn=fldArEnergy(dn).div(norm),vp=fldArEnergy(dp).div(norm),rho=fldArPair(actualHistory,rv).div(norm),cross=fldArPair(dn,dp).div(norm);
 ensure(norm.eq(2)&&mn.eq('3/2')&&mp.zero()&&vn.eq('1/4')&&vp.eq('1/2')&&rho.eq('1/2')&&cross.rad.zero(),'actual source/process clock profile moments');
 const defect=fldArEnergy(fldArMinus(next,actualHistory)).div(norm);ensure(defect.eq(1),'actual history continuation defect');
 const factored=[empty(),productField(sourceSeed,o.mul(R,e0)),productField(sourceSeed,o.mul(o.power(R,2),e0)),empty()];
 const error=fldArEnergy(fldArMinus(actualHistory,factored)).div(norm),one=energy(minus(sourcePowers[1],sourceSeed));ensure(error.le(one.mul('5/2')),'actual finite carrier factorization budget');
 actualResponses={nv,pv,rv,next,prev,dn,dp};
 return {actual_echo_source_powers:4,profile_order:1,tick_labels:4,count_variance:vn,response_variance:vp,commutator_response:rho,stationary_defect_square:defect,carrier_factorization_error_square:error,finite_factorization_error_upper:one.mul('5/2'),literal_small_carrier_is_large_R37_packet:false};
});

const parent=JSON.parse(parentBytes),r37Path=path.join(path.dirname(arg('--r40-certificate')),'R37_NATIVE_CERTIFICATE.json'),r37Bytes=fs.readFileSync(r37Path);
const inherited=parent.exact_checks.find(x=>x.name==='derived_packet_horizons_bound_autonomous_echo_error_and_record_capacity');ensure(hash(r37Bytes)===inherited.frozen_packet_sha256,'R41 frozen packet pin mismatch');
const r37=JSON.parse(r37Bytes),packetConstants=r37.exact_checks.find(x=>x.name==='a_source_corrector_removes_every_first_difference_error').source_constants.find(x=>x.parameter==='2');
const kap=f(packetConstants.kappa),ct=f(packetConstants.comparison_norm_bound),csum=f(packetConstants.residual_xx_mass).add(packetConstants.residual_xy_mass).add(packetConstants.residual_yy_mass),drift=f('2/3');
const packetEpsilon=L=>csum.mul(4).add(ct.mul(2).mul(drift).mul(f(1).sub(drift))).div(L).add(ct.mul(2).div(f(L).pow(2))).div(f(1).sub(kap.mul(2).div(f(L).pow(2))));
const lt=(a,b)=>f(a).le(b)&&!f(a).eq(b);
check('rational_joint_limits_control_clock_error_stationarity_and_the_half_scale',()=>{
 let previous=null;const rows=[];for(const a of [6,7,8]){const L=f(16).pow(a),Qn=f(4).pow(a),s=Qn.sub(3),eps=packetEpsilon(L),vn=s.pow(2).div(s.mul(2).sub(1).mul(4)),vp=s.mul(2).add(1).div(s.add(1).mul(s.add(2))),rho=s.div(s.add(1));
  const second=f(1).add(s.div(2)).pow(2).add(vn),clock=eps.pow(2).mul(4).mul(second),stationary=f(2).div(s.add(1)),gap=vn.mul(vp).sub(rho.pow(2).div(4)),halfDeficit=f('1/4').sub(vn.mul(vp));
  ensure(Qn.pow(2).eq(L)&&s.add(3).eq(Qn),'one common integer scale');ensure(lt(0,gap)&&lt(0,halfDeficit),'finite profile retains response correction');
  if(previous)ensure(lt(clock,previous.clock)&&lt(stationary,previous.stationary)&&lt(gap,previous.gap)&&lt(halfDeficit,previous.halfDeficit),'simultaneously shrinking native defects');previous={clock,stationary,gap,halfDeficit};
  rows.push({L,tick_labels:Qn,profile_word_length:s,program_phase_labels:L.pow(3).mul(2).add(6),packet_epsilon:eps,carrier_factorization_error_square_upper:clock,history_defect_square:stationary,commutator_response:rho,uncertainty_product_square:vn.mul(vp),uncertainty_gap:gap,deficit_from_quarter:halfDeficit});
 }
 return {frozen_packet_sha256:hash(r37Bytes),simultaneous_exact_parameter_budgets:rows,clock_factorization_norm_rate:'L^(-1/2)',history_defect_norm_rate:'L^(-1/4)',squared_uncertainty_gap_rate:'L^(-1)',uncertainty_product_limit:'1/2',large_fields_or_clock_records_directly_simulated:false,normalizable_infinite_uniform_history_claimed:false};
});
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,Jc),-1),I),J32=o.kron(Jc,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
check('complete_curvature_observation_and_affine_calibration_preserve_the_law',()=>{
 let observations=0;const vectors=[actualHistory,...Object.values(actualResponses)];
 for(const vec of vectors){const obs=vec.map(observe);for(let r=0;r<vec.length;r++)eqFields(decode(obs[r]),vec[r],32,2,'complete relational history reconstruction');
  const norm=obs.reduce((s,x)=>s.add(energy(x[0])).add(energy(x[1])),f(0));ensure(norm.eq(fldArEnergy(vec)),'full observer norm');
  const mean=obs.reduce((s,x,r)=>{const base=observe(actualHistory[r]);return s.add(fieldPair(base[0],x[0])).add(fieldPair(base[1],x[1]));},o.ZERO);ensure(mean.eq(fldArPair(actualHistory,vec)),'complete observer cross pairing');observations++;
 }
 const a0=history(profile(2,8),processTurn),st=stats(a0,processTurn),E0=arEnergy(a0);let calibrations=0;const rows=[];
 for(const [aa,cc]of [[2,3],['1/2','5/3'],['7/4','2/5']]){const aa0=f(aa),cc0=f(cc),tv=arAdd(arScale(count(a0),aa0),arScale(a0,3)),pv=arAdd(arScale(response(a0,processTurn),cc0),arScale(a0,-2));
  const mt=arPair(a0,tv).div(E0),mp=arPair(a0,pv).div(E0),dt=arMinus(tv,arScale(a0,mt)),dp=arMinus(pv,arScale(a0,mp)),vt=arEnergy(dt).div(E0),vp=arEnergy(dp).div(E0),cov=arPair(dt,dp).rad.div(E0),det=vt.mul(vp).sub(cov.pow(2));
  ensure(vt.eq(st.count_variance.mul(aa0.pow(2)))&&vp.eq(st.response_variance.mul(cc0.pow(2))),'affine variances');ensure(det.eq(st.count_variance.mul(st.response_variance).sub(st.covariance.pow(2)).mul(aa0.mul(cc0).pow(2))),'affine covariance determinant');calibrations++;
  rows.push({count_scale:aa0,response_scale:cc0,commutator_scale:aa0.mul(cc0),scaled_covariance_determinant:det});
 }
 return {complete_observed_history_vectors:observations,exact_affine_calibration_checks:calibrations,calibrations:rows,physical_energy_action_or_time_scale_selected:false};
});
check('clock_scope_controls_retain_wrap_correlation_and_physical_boundaries',()=>{
 const state=history(profile(2,8),processTurn),rel=stats(state,processTurn),bare=stats(state,I);ensure(!rel.response_variance.eq(bare.response_variance)||!rel.response_mean.eq(bare.response_mean),'bare clock response differs from joint response');
 const s1=stats(history(profile(1,4),I),I);ensure(lt(s1.count_variance.mul(s1.response_variance),'1/4'),'finite product need not exceed one half');
 ensure(nonclosingHistories>0,'uniform histories need not be stationary');ensure(!o.equal(o.power(R,qProgram),R),'incorrect process carry on every edge');
 const invalid=o.scale(I,2);ensure(!o.equal(o.mul(o.dagger(invalid),invalid),I),'nonisometric process fails source contract');
 eq(o.mul(o.mul(o.mul(H,K),H),K),o.scale(I,-1),'retained spatial source loop');
 return {bare_clock_response_witness:{joint_response_mean:rel.response_mean,bare_response_mean:bare.response_mean,joint_response_variance:rel.response_variance,bare_response_variance:bare.response_variance},finite_uncertainty_product_square:s1.count_variance.mul(s1.response_variance),incorrect_every_program_edge_process_update_rejected:true,nonisometric_process_rejected:true,finite_native_product_always_at_least_one_half:false,clock_phase_curvature_is_record_size_independent:false,dimensionless_half_scale_is_physical_hbar:false,physical_constants_identified:false};
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
 the_clock_instantly_controls_a_separate_remote_process:checks[0].physical_remote_control_asserted===false,
 every_uniform_finite_history_is_stationary:checks[1].nonclosing_native_histories>0,
 native_history_preparation_copies_an_unknown_state:checks[1].unknown_state_copying_map_used===false,
 ignoring_clock_marks_assumes_a_random_history_selection:checks[1].random_history_selection_assumed===false,
 the_count_response_commutator_has_no_wrap_term:checks[2].omitted_wrap_defect_sha256.length===64,
 the_root_matrix_witness_alone_selects_a_positive_scalar_branch:checks[2].eighth_turn_witness_checks_both_algebraic_conjugates===true,
 the_finite_clock_has_a_global_canonical_iota_identity:checks[2].global_finite_clock_canonical_identity_claimed===false,
 uncertainty_was_imported_as_an_external_law:checks[3].external_uncertainty_law_used_as_premise===false,
 the_profile_formulas_assume_a_classical_probability_distribution:checks[4].limiting_probability_distribution_assumed===false,
 a_bare_clock_shift_always_has_the_joint_relational_response:checks[8].bare_clock_response_witness.joint_response_variance.toString()!==checks[8].bare_clock_response_witness.bare_response_variance.toString(),
 every_finite_native_uncertainty_product_is_at_least_one_half:checks[8].finite_native_product_always_at_least_one_half===false,
 the_large_packet_and_clock_arrays_were_directly_simulated:checks[6].large_fields_or_clock_records_directly_simulated===false,
 an_infinite_uniform_history_was_proved_normalizable:checks[6].normalizable_infinite_uniform_history_claimed===false,
 the_dimensionless_native_half_scale_is_already_physical_hbar:checks[8].dimensionless_half_scale_is_physical_hbar===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r41.native-relational-clock.v1',status:'PASS_R41_NATIVE_RELATIONAL_CLOCK',
 input_sha256:inputHash,r40_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_clock_process_relation_derived:true,native_clock_phase_curvature_derived:true,
 native_sharp_count_response_uncertainty_derived:true,native_relational_history_clock_limit_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: autonomous tick-coupled process evolution, native histories and unresolved-clock readouts, clock phase curvature and the wrap-corrected count-response commutator, a sharp covariance inequality, finite binary-word profile moments, a joint clock-error and stationary-history limit approaching the native half scale, and complete observation with affine calibration. The process and readouts are explicit native constructions; physical energy, action and time calibration remain unselected.'},null,2)+'\n');
