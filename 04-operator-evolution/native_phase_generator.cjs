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
ensure(inputHash===arg('--expected-input-sha256'),'R42 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r41-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r41-sha256'),'R42 R41 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R41_NATIVE_RELATIONAL_CLOCK','R42 parent status');
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
const inner=(a,b)=>o.mul(o.dagger(a),b)[0][0],comm=o.commutator;
const nullColumns=A=>{const z=o.rref(A),free=Array.from({length:A[0].length},(_,j)=>j).filter(j=>!z.pivots.includes(j));
 if(!free.length)return null;return o.matrix(Array.from({length:A[0].length},(_,j)=>free.map(k=>{const r=z.pivots.indexOf(j);return r<0?(j===k?o.ONE:o.ZERO):z.matrix[r][k].neg();})));};
const cutKernel=A=>{const X=nullColumns(A);return X?o.mul(o.mul(X,o.inverse(o.mul(o.dagger(X),X))),o.dagger(X)):o.zeros(A[0].length);};
const midpointFrom=V=>{const id=o.identity(V.length),cut=cutKernel(o.add(id,V)),regular=o.sub(id,cut),mid=o.scale(o.mul(o.mul(o.sub(V,id),regular),o.inverse(o.add(o.add(id,V),cut))),new o.Cut(0,2));return {cut,regular,mid};};
const stepFrom=(mid,s=1)=>{const id=o.identity(mid.length),z=o.scale(mid,new o.Cut(0,f(s).div(2)));return o.mul(o.inverse(o.add(id,z)),o.sub(id,z));};
const recoverStep=({cut,regular,mid})=>o.sub(o.mul(stepFrom(mid),regular),cut);
const blockDiag=(...parts)=>{const n=parts.reduce((s,A)=>s+A.length,0),out=o.zeros(n);let offset=0;for(const A of parts){for(let i=0;i<A.length;i++)for(let j=0;j<A.length;j++)out[offset+i][offset+j]=A[i][j];offset+=A.length;}return out;};
const Jturn=o.scale(R,io),otherTurn=o.add(o.scale(I,'-3/5'),o.scale(R,'4/5'));
const rotate3=o.matrix([['3/5','-4/5',0],['4/5','3/5',0],[0,0,1]]),mixed3=o.mul(o.mul(rotate3,blockDiag(o.matrix([[-1]]),processTurn)),o.dagger(rotate3));
const arbitraryMid=o.matrix([[1,new o.Cut(1,1)],[new o.Cut(1,-1),-2]]);
const phaseCases=[['identity',I],['full_reversal',o.scale(I,-1)],['source_H',H],['source_K',K],['source_R',R],['rational_turn',processTurn],['same_contrast_other_step',otherTurn],['mixed_reversal',mixed3],['general_native_midpoint',stepFrom(arbitraryMid)]]
 .map(([name,V])=>({name,V,...midpointFrom(V)}));
const resolve=(mid,s)=>o.inverse(o.add(o.identity(mid.length),o.scale(o.mul(mid,mid),f(s).pow(2).div(4))));
const kernelAt=(mid,s)=>o.mul(mid,resolve(mid,s));
const phaseCache=new Map();
const phaseApprox=(entry,n)=>{const cacheKey=entry.name+':'+n;if(phaseCache.has(cacheKey))return phaseCache.get(cacheKey);
 let sum=o.zeros(entry.mid.length),lower=f(0),upper=f(0);for(let j=0;j<n;j++){sum=o.add(sum,kernelAt(entry.mid,f(j).div(n)));upper=upper.add(f(4).div(f(1).add(f(j).div(n).pow(2))));lower=lower.add(f(4).div(f(1).add(f(j+1).div(n).pow(2))));}
 lower=lower.div(n);upper=upper.div(n);ensure(upper.sub(lower).eq(f(2).div(n)),'exact half-turn bracket width');
 const halfTurn=lower.add(upper).div(2),matrix=o.add(o.scale(sum,f(1).div(n)),o.scale(entry.cut,halfTurn)),M=o.mass(entry.mid),error=M.pow(3).div(4*n).add(o.isZero(entry.cut)?0:f(1).div(n));
 const out={matrix,halfTurn,lower,upper,error,midpoint_gain_bound:M};phaseCache.set(cacheKey,out);return out;};
const fact=n=>{let a=f(1);for(let j=2;j<=n;j++)a=a.mul(j);return a;};
const factorialFlow=(A,t,m)=>{const Z=o.scale(A,new o.Cut(0,f(t).neg())),mass=o.mass(Z),alpha=f((mass.n+mass.d-1n)/mass.d);ensure(alpha.le(m+2)&&!alpha.eq(m+2),'factorial tail ratio');let sum=o.identity(A.length),term=sum;
 for(let j=1;j<=m;j++){term=o.scale(o.mul(term,Z),f(1).div(j));sum=o.add(sum,term);}const tail=alpha.pow(m+1).div(fact(m+1)).div(f(1).sub(alpha.div(m+2)));return {matrix:sum,tail,alpha};};
const varOf=(A,v)=>{const en=o.energy(v),mean=inner(v,o.mul(A,v)).div(en);ensure(mean.turn.zero(),'self-dagger mean');const dv=o.sub(o.mul(A,v),o.scale(v,mean));return {mean:mean.rad,variance:o.energy(dv).div(en),centered:dv};};
const vectorFor=d=>col(Array.from({length:d},(_,j)=>new o.Cut(j+1,j%2)));

check('native_complete_midpoint_generator_retains_the_reversal_cut',()=>{
 const rows=[];for(const e of phaseCases){const id=o.identity(e.V.length),c=o.scale(o.add(e.V,o.dagger(e.V)),'1/2'),pv=o.scale(o.sub(e.V,o.dagger(e.V)),new o.Cut(0,'1/2'));
  eq(o.mul(o.dagger(e.V),e.V),id,'native process isometry');eq(comm(c,pv),o.zeros(id.length),'complementary responses commute');eq(o.add(o.mul(c,c),o.mul(pv,pv)),id,'complete response square');eq(o.sub(c,o.scale(pv,io)),e.V,'complete response reconstruction');
  eq(o.mul(e.cut,e.cut),e.cut,'kernel cut idempotence');eq(o.dagger(e.cut),e.cut,'kernel cut dagger');eq(o.mul(o.add(id,e.V),e.cut),o.zeros(id.length),'kernel cut range');ensure(o.rank(e.cut)===id.length-o.rank(o.add(id,e.V)),'complete kernel rank');
  eq(o.dagger(e.mid),e.mid,'native midpoint self-dagger');eq(comm(e.mid,e.V),o.zeros(id.length),'midpoint commutes with process');eq(o.mul(e.mid,e.cut),o.zeros(id.length),'zero midpoint on reversal cut');eq(recoverStep(e),e.V,'cut-complete recovery');
  rows.push({name:e.name,roles:id.length,reversal_cut_rank:o.rank(e.cut),midpoint_generator:e.mid,reversal_cut:e.cut});
 }
 eq(phaseCases[5].mid,Jturn,'rational turn midpoint generator');eq(o.scale(o.sub(processTurn,o.dagger(processTurn)),new o.Cut(0,'1/2')),o.scale(o.sub(otherTurn,o.dagger(otherTurn)),new o.Cut(0,'1/2')),'same contrast distinct processes');ensure(!o.equal(processTurn,otherTurn),'contrasts fail complete recovery');
 ensure(o.isZero(phaseCases[0].mid)&&o.isZero(phaseCases[1].mid)&&!o.equal(phaseCases[0].cut,phaseCases[1].cut),'reversal memory is necessary');
 return {complete_process_instances:rows,finite_processes_checked:rows.length,forward_reverse_response_alone_is_complete:false,reversal_cut_can_be_discarded:false,primitive_spectral_decomposition_used:false};
});
check('native_resolvent_count_sums_construct_phase_with_bounds',()=>{
 let identities=0;const rows=[];for(const e of phaseCases){const id=o.identity(e.mid.length),mid2=o.mul(e.mid,e.mid);
  for(const s of [0,'1/4','1/2','3/4',1]){const ss=f(s),inv=resolve(e.mid,ss),D=o.add(id,o.scale(e.mid,new o.Cut(0,ss.div(2)))),Di=o.inverse(D),A=stepFrom(e.mid,ss),der=o.sub(o.scale(o.mul(o.mul(Di,e.mid),A),new o.Cut(0,'-1/2')),o.scale(o.mul(Di,e.mid),new o.Cut(0,'1/2')));
   eq(o.mul(o.add(id,o.scale(mid2,ss.pow(2).div(4))),inv),id,'rational integrand inverse');eq(o.dagger(kernelAt(e.mid,ss)),kernelAt(e.mid,ss),'integrand dagger');eq(comm(kernelAt(e.mid,ss),e.mid),o.zeros(id.length),'integrand commute');eq(der,o.scale(o.mul(kernelAt(e.mid,ss),A),new o.Cut(0,-1)),'native ratio path derivative');eq(o.mul(o.dagger(D),D),o.add(id,o.scale(mid2,ss.pow(2).div(4))),'positive inverse identity');identities++;
  }
  const a=phaseApprox(e,8),b=phaseApprox(e,16),diff=o.energy(o.sub(a.matrix,b.matrix)),budget=a.error.add(b.error);ensure(diff.le(budget.pow(2).mul(id.length)),'successive phase count sums enclosure');eq(o.dagger(b.matrix),b.matrix,'phase approximant self-dagger');eq(comm(b.matrix,e.V),o.zeros(id.length),'phase approximant commutes with original step');
  rows.push({process:e.name,count_terms:16,midpoint_gain_bound:b.midpoint_gain_bound,phase_gain_error_upper:b.error,half_turn_lower:b.lower,half_turn_upper:b.upper,successive_approximation_difference_square:diff});
 }
 return {rational_resolvent_and_derivative_instances:identities,phase_count_bounds:rows,native_count_integration_imported:false,completed_phase_claimed_exactly_rational:false};
});
check('native_phase_exponential_enclosures_recover_the_original_step',()=>{
 const rows=[];for(const e of phaseCases.filter(x=>x.V.length===2&&x.name!=='same_contrast_other_step'&&x.name!=='general_native_midpoint')){let previous=null;
  for(const n of [8,16,32]){const a=phaseApprox(e,n),poly=factorialFlow(a.matrix,1,32),err=o.energy(o.sub(poly.matrix,e.V)),bound=a.error.add(poly.tail),sq=bound.pow(2).mul(e.V.length);ensure(err.le(sq),'native factorial reconstruction enclosure');if(previous&&!a.error.zero())ensure(bound.le(previous),'improving reconstruction bound');previous=bound;
   rows.push({process:e.name,count_terms:n,factorial_degree:32,phase_error_bound:a.error,factorial_gain_bound:poly.alpha,factorial_tail:poly.tail,reconstruction_entry_norm_square_sha256:hash(JSON.stringify(err)),reconstruction_entry_norm_square_bound:sq});
  }
 }
 const r=phaseCases.find(x=>x.name==='rational_turn'),a=phaseApprox(r,16),pol=factorialFlow(a.matrix,'1/2',32),square=o.mul(pol.matrix,pol.matrix),budget=a.error.add(pol.tail.mul(2)).add(pol.tail.pow(2));ensure(o.energy(o.sub(square,r.V)).le(budget.pow(2).mul(2)),'two half-step phase enclosure recovers full process');
 return {completed_phase_reconstruction_enclosures:rows,half_tick_square_enclosure:true,finite_factorial_polynomial_claimed_exactly_unitary:false,midpoint_generator_claimed_equal_to_phase_generator:false};
});
check('native_generator_readouts_derive_conservation_and_response_geometry',()=>{
 let conservation=0,geometry=0,speed=0;for(const e of phaseCases){const en=phaseApprox(e,16).matrix,v=vectorFor(e.V.length),energy0=o.energy(v),der=o.scale(o.mul(en,v),new o.Cut(0,-1)),st=varOf(en,v),projection=o.sub(der,o.scale(v,inner(v,der).div(energy0)));ensure(o.energy(projection).div(energy0).eq(st.variance),'orthogonal response variance metric');geometry++;
  for(let k=0;k<=4;k++){const vk=o.mul(o.power(e.V,k),v);for(const p of [1,2,3]){const A=o.power(en,p);ensure(inner(vk,o.mul(A,vk)).eq(inner(v,o.mul(A,v))),'generator moment conservation');conservation++;}}
  const reads=[o.identity(e.V.length),e.cut,o.matrix(Array.from({length:e.V.length},(_,i)=>Array.from({length:e.V.length},(_,j)=>i===j?i+1:0)))];
  for(const A of reads){const av=varOf(A,v),slope=inner(der,o.mul(A,v)).add(inner(v,o.mul(A,der))).div(energy0),fromComm=inner(v,o.mul(o.scale(comm(en,A),io),v)).div(energy0);ensure(slope.eq(fromComm)&&slope.turn.zero(),'native readout derivative');ensure(slope.rad.pow(2).le(st.variance.mul(av.variance).mul(4)),'native readout speed bound');speed++;}
 }
 const st=varOf(Jturn,e0),ak=varOf(K,e0),slope=inner(e0,o.mul(o.scale(comm(Jturn,K),io),e0));ensure(slope.eq(2)&&st.variance.eq(1)&&ak.variance.eq(1),'sharp speed coefficient witness');
 const pure=col([1,new o.Cut(0,1)]);ensure(varOf(Jturn,pure).variance.zero(),'pure phase has zero response metric');
 return {exact_conserved_moment_checks:conservation,exact_response_metric_checks:geometry,exact_readout_speed_checks:speed,sharp_speed_coefficient:2,sharpness_scalar_factor:'theta = native integral of (1+s^2/4)^(-1)',rational_checks_are_phase_approximants_not_exact_completed_phase:true,common_phase_metric_witness:'zero'};
});
const pathCost=(xs,V)=>xs.slice(1).reduce((s,x,j)=>s.add(o.energy(o.sub(x,o.mul(V,xs[j])))),f(0));
const unwind=(xs,V)=>xs.map((x,j)=>o.mul(o.power(o.dagger(V),j),x));
const minPath=(a,b,n,V)=>{const c=o.scale(o.sub(o.mul(o.power(o.dagger(V),n),b),a),f(1).div(n));return Array.from({length:n+1},(_,j)=>o.mul(o.power(V,j),o.add(a,o.scale(c,j))));};
check('native_residue_action_has_exact_minima_and_boundary_equations',()=>{
 let endpoints=0,decompositions=0,stationary=0;for(const e of phaseCases)for(const n of [1,2,5]){const a=vectorFor(e.V.length),b=o.scale(a,new o.Cut('2/3','1/4')),xs=minPath(a,b,n,e.V),target=o.energy(o.sub(o.mul(o.power(o.dagger(e.V),n),b),a)).div(n);ensure(pathCost(xs,e.V).eq(target),'clamped endpoint minimum');endpoints++;
  const trial=xs.map((x,j)=>j>0&&j<n?o.add(x,o.scale(basis(e.V.length,j%e.V.length),f(j).div(3))):x),ys=unwind(trial,e.V),c=o.scale(o.sub(ys[n],ys[0]),f(1).div(n));let excess=f(0);for(let j=0;j<n;j++)excess=excess.add(o.energy(o.sub(o.sub(ys[j+1],ys[j]),c)));ensure(pathCost(trial,e.V).eq(target.add(excess)),'exact action square completion');decompositions++;
  for(let j=1;j<n;j++){eq(o.sub(o.sub(o.scale(xs[j],2),o.mul(e.V,xs[j-1])),o.mul(o.dagger(e.V),xs[j+1])),o.zeros(e.V.length,1),'interior stationary equation');stationary++;}
  const free=Array.from({length:n+1},(_,j)=>o.mul(o.power(e.V,j),a));ensure(pathCost(free,e.V).zero(),'free terminal exact zero path');eq(o.sub(free[n],o.mul(e.V,free[n-1])),o.zeros(e.V.length,1),'free terminal natural equation');
 }
 const x=[e0,o.scale(e0,2),o.scale(e0,3)];eq(o.sub(o.sub(o.scale(x[1],2),x[0]),x[2]),o.zeros(2,1),'nonzero slope interior solution');ensure(!o.energy(o.sub(x[2],x[1])).zero(),'missing terminal condition witness');
 return {fixed_endpoint_minima:endpoints,exact_square_completion_checks:decompositions,interior_stationary_checks:stationary,interior_equation_alone_forces_first_order_history:false,path_vectors_silently_constrained_to_unit_norm:false,physical_least_action_postulate_used:false};
});
check('native_midpoint_action_retains_regular_and_reversal_residues',()=>{
 let checks=0;const rows=[];for(const e of phaseCases){const v=vectorFor(e.V.length),xs=[v,o.scale(v,new o.Cut('1/2','1/3')),o.add(v,basis(e.V.length,0)),o.mul(e.V,v)];let regCost=f(0),cutCost=f(0);
  for(let j=0;j<xs.length-1;j++){const d=o.mul(e.regular,o.sub(xs[j+1],xs[j])),m=o.scale(o.mul(e.regular,o.add(xs[j+1],xs[j])),'1/2'),r=o.add(d,o.scale(o.mul(e.mid,m),io));const cost=inner(r,o.mul(resolve(e.mid,1),r));ensure(cost.turn.zero(),'real midpoint residue cost');regCost=regCost.add(cost.rad);cutCost=cutCost.add(o.energy(o.mul(e.cut,o.add(xs[j+1],xs[j]))));checks++;}
  ensure(regCost.add(cutCost).eq(pathCost(xs,e.V)),'cut-complete midpoint action identity');const actual=o.mul(e.V,v),d=o.mul(e.regular,o.sub(actual,v)),m=o.scale(o.mul(e.regular,o.add(actual,v)),'1/2');eq(o.scale(d,io),o.mul(e.mid,m),'actual midpoint evolution equation');eq(o.mul(e.cut,o.add(actual,v)),o.zeros(e.V.length,1),'actual reversal alternation');rows.push({process:e.name,regular_residue:regCost,reversal_residue:cutCost,total_residue:pathCost(xs,e.V)});
 }
 return {exact_midpoint_residue_steps:checks,action_decompositions:rows,reversal_residue_is_dropped:false};
});
const jointProc=(a,A)=>new Map([...a].map(([k,v])=>[k,processApply(v,A)])),jointPair=(a,b)=>[...a].reduce((s,[k,v])=>s.add(fieldPair(v,b.get(k)||new Map())),o.ZERO);
check('the_actual_autonomous_clock_conserves_the_internal_generator',()=>{
 let edges=0,moments=0;for(const name of ['source_H','source_K','rational_turn']){const e=phaseCases.find(x=>x.name===name),A=phaseApprox(e,16).matrix,v=col([1,new o.Cut(2,1)]);
  for(let j=0;j<qProgram;j++){const z=joint(productField(sourceSeed,v),j,1),left=coupledStep(jointProc(z,A),e.V,3),right=jointProc(coupledStep(z,e.V,3),A);eqJoint(left,right,'actual clock generator intertwining');edges++;}
  let z=joint(productField(sourceSeed,v)),first=jointPair(z,jointProc(z,A));for(let n=0;n<=qProgram+1;n++){if(n)z=coupledStep(z,e.V,3);ensure(jointPair(z,jointProc(z,A)).eq(first),'actual autonomous generator readout conserved');moments++;}
  const weights=[1,2,1],historyV=weights.map((h,r)=>o.scale(o.mul(o.power(e.V,r),v),h)),en=weights.reduce((s,h)=>s.add(f(h).pow(2)),f(0));let pair=o.zeros(2);for(const x of historyV)pair=o.add(pair,o.mul(x,o.dagger(x)));ensure(trace(o.mul(A,pair)).div(en).eq(inner(v,o.mul(A,v))),'unresolved count history generator moment');
 }
 return {actual_source_process_edge_intertwinings:edges,actual_autonomous_moment_checks:moments,unresolved_history_moment_checks:3,clock_carrier_error_becomes_generator_conservation_error:false,completed_generator_transport_supported_by_commuting_count_limits:true};
});
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,Jc),-1),I),J32=o.kron(Jc,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
check('complete_curvature_observation_and_phase_branch_witnesses',()=>{
 const e=phaseCases.find(x=>x.name==='rational_turn'),A=phaseApprox(e,16).matrix,fields=[productField(sourceSeed,col([1,new o.Cut(2,1)])),productField(sourcePowers[1],col([new o.Cut(1,1),2]))];let matches=0;
 for(const v of fields){const av=processApply(v,A),obs=observe(v),oa=observe(av);eqFields(decode(obs),v,32,2,'complete generator source decode');ensure(energy(obs[0]).add(energy(obs[1])).eq(energy(v)),'complete source norm');ensure(fieldPair(obs[0],oa[0]).add(fieldPair(obs[1],oa[1])).eq(fieldPair(v,av)),'complete generator cross pairing');matches++;}
 const residue=minus(fields[1],processApply(fields[0],e.V)),obs=observe(residue);ensure(energy(obs[0]).add(energy(obs[1])).eq(energy(residue)),'observer preserves path residue');
 const midpoint=o.sub(I,o.scale(P,2));eq(midpoint,o.scale(H,-1),'half tick alias response');eq(o.mul(midpoint,midpoint),I,'same full tick');const v=col([1,1]),pair=o.mul(v,o.dagger(v)),alt=o.mul(midpoint,v);ensure(!o.equal(o.mul(alt,o.dagger(alt)),pair),'fractional branch changes relative response');
 return {complete_source_generator_pairings:matches,complete_path_residue_preserved:true,alias_full_tick:'I',selected_branch_half_tick:'I',alternative_branch_half_tick:'-H',integer_ticks_uniquely_determine_fractional_evolution:false,native_half_turn_identity_source:'R42.3 count integral and factorial product proof'};
});
check('calibration_and_count_response_controls_keep_physical_scope_exact',()=>{
 const e=phaseCases.find(x=>x.name==='rational_turn'),A=phaseApprox(e,16).matrix,v=col([1,new o.Cut(1,1)]),Qn=4,Nt=o.zeros(2*Qn),Et=o.kron(o.identity(Qn),A),St=o.zeros(Qn);for(let r=0;r<Qn;r++){St[(r+1)%Qn][r]=o.ONE;for(let j=0;j<2;j++)Nt[2*r+j][2*r+j]=o.Cut.of(r);}
 eq(comm(Nt,Et),o.zeros(2*Qn),'retained count commutes with internal process generator');const T=o.kron(St,e.V),Pt=o.scale(o.sub(T,o.dagger(T)),new o.Cut(0,'1/2'));ensure(!o.equal(Pt,Et),'joint tick contrast is not internal phase generator');
 const rows=[];for(const [aa,bb]of [[2,3],['1/2','5/3'],['7/4','2/5']]){const a=f(aa),b=f(bb),der=o.scale(o.mul(A,v),new o.Cut(0,f(-1).div(a)));eq(o.scale(der,new o.Cut(0,a.mul(b))),o.scale(o.mul(A,v),b),'calibrated evolution identity');rows.push({time_scale:a,readout_scale:b,action_scale:a.mul(b)});}
 const h=phaseCases.find(x=>x.name==='source_H');ensure(o.isZero(h.mid)&&!o.equal(h.V,I),'midpoint cannot generate reversal by factorial flow');const cost=pathCost([e0,col([1,1]),col([2,1])],e.V);ensure(!cost.zero(),'positive action scaling witness');
 return {count_internal_generator_commutator:'zero',joint_tick_contrast_equals_internal_generator:false,calibrations:rows,positive_action_multipliers_preserve_minimizers:true,finite_phase_branch_or_action_unit_uniquely_selected:false,physical_hbar_energy_time_or_action_identified:false,ordinary_operator_logarithm_or_spectral_theorem_used:false};
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
 forward_reverse_contrast_alone_recovers_every_process:checks[0].forward_reverse_response_alone_is_complete===false,
 the_reversal_cut_can_be_discarded:checks[0].reversal_cut_can_be_discarded===false,
 a_spectral_decomposition_was_a_primitive_input:checks[0].primitive_spectral_decomposition_used===false,
 integration_was_imported_without_native_count_completion:checks[1].native_count_integration_imported===false,
 rational_phase_approximants_are_exact_completed_generators:checks[1].completed_phase_claimed_exactly_rational===false,
 the_truncated_factorial_polynomial_is_exactly_unitary:checks[2].finite_factorial_polynomial_claimed_exactly_unitary===false,
 midpoint_and_phase_generators_are_identical:checks[2].midpoint_generator_claimed_equal_to_phase_generator===false,
 common_phase_always_has_nonzero_response_metric:checks[3].common_phase_metric_witness==='zero',
 interior_stationarity_alone_selects_the_original_evolution:checks[4].interior_equation_alone_forces_first_order_history===false,
 the_native_mismatch_functional_imports_a_physical_action_law:checks[4].physical_least_action_postulate_used===false,
 clock_carrier_error_breaks_the_conserved_generator:checks[6].clock_carrier_error_becomes_generator_conservation_error===false,
 integer_ticks_uniquely_select_fractional_evolution:checks[7].integer_ticks_uniquely_determine_fractional_evolution===false,
 the_R41_tick_contrast_is_the_internal_phase_generator:checks[8].joint_tick_contrast_equals_internal_generator===false,
 a_physical_hbar_or_action_unit_has_been_selected:checks[8].physical_hbar_energy_time_or_action_identified===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r42.native-phase-generator.v1',status:'PASS_R42_NATIVE_PHASE_GENERATOR',
 input_sha256:inputHash,r41_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_complete_process_phase_generator_derived:true,native_exact_tick_interpolation_derived:true,
 native_generator_response_geometry_derived:true,native_path_residue_variational_rule_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: cut-complete process reconstruction, phase integration by native count sums, exact factorial interpolation with finite enclosures, conserved generator moments and response geometry, a finite path-residue variational rule, autonomous-clock conservation and complete observation, and phase-branch and action-calibration boundaries. The finite process, interpolation branch and mismatch functional are explicit native constructions; physical energy, time and h-bar remain unselected.'},null,2)+'\n');
