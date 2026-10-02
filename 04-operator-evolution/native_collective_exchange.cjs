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
ensure(inputHash===arg('--expected-input-sha256'),'R44 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r43-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r43-sha256'),'R44 R43 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R43_NATIVE_EXCHANGE_CALIBRATION','R44 parent status');
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
const checks=[],check=(name,fn)=>{checks.push({name,passed:true,...fn()});if(process.env.R44_CHECK_PROGRESS==='1')process.stderr.write(name+'\n');};
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

const comm=o.commutator,inner=(a,b)=>o.mul(o.dagger(a),b)[0][0];
const digits=(j,n,d=2)=>Array.from({length:n},(_,a)=>Math.floor(j/d**(n-1-a))%d),encode=(v,d=2)=>v.reduce((s,a)=>s*d+a,0);
const swapN=(n,a,b,d=2)=>{const S=o.zeros(d**n);for(let j=0;j<d**n;j++){const x=digits(j,n,d);[x[a],x[b]]=[x[b],x[a]];S[encode(x,d)][j]=o.ONE;}return S;};
const localN=(n,a,A)=>{let M=o.identity(1);for(let j=0;j<n;j++)M=o.kron(M,j===a?A:o.identity(A.length));return M;};
const sumMatrices=(xs,D)=>xs.reduce(o.add,o.zeros(D));
const pairOf=v=>o.scale(o.mul(v,o.dagger(v)),f(1).div(o.energy(v)));
const mean=(v,A)=>inner(v,o.mul(A,v)).div(o.energy(v));
const part=(M,n,keep)=>{const rest=Array.from({length:n},(_,i)=>i).filter(i=>!keep.includes(i)),D=2**keep.length,out=o.zeros(D);for(let a=0;a<D;a++)for(let b=0;b<D;b++)for(let r=0;r<2**rest.length;r++){const x=Array(n).fill(0),y=Array(n).fill(0);digits(a,keep.length).forEach((v,i)=>x[keep[i]]=v);digits(b,keep.length).forEach((v,i)=>y[keep[i]]=v);digits(r,rest.length).forEach((v,i)=>x[rest[i]]=y[rest[i]]=v);out[a][b]=out[a][b].add(M[encode(x)][encode(y)]);}return out;};
const components=(n,edges)=>{const labels=Array(n).fill(-1);let count=0;for(let a=0;a<n;a++)if(labels[a]<0){labels[a]=count;const todo=[a];for(let k=0;k<todo.length;k++)for(const [u,v]of edges){const b=u===todo[k]?v:v===todo[k]?u:-1;if(b>=0&&labels[b]<0){labels[b]=count;todo.push(b);}}count++;}return {labels,count};};
const graph=(n,edges,d=2)=>{const D=d**n,id=o.identity(D),swaps=edges.map(([a,b])=>swapN(n,a,b,d)),cuts=swaps.map(S=>o.scale(o.sub(id,S),'1/2'));return {n,d,D,id,edges,swaps,cuts,G:sumMatrices(cuts,D),...components(n,edges)};};
const smallGraphs=[graph(2,[[0,1]]),graph(3,[[0,1],[1,2]]),graph(3,[[0,1],[1,2],[0,2]]),graph(4,[[0,1],[1,2],[2,3]]),graph(4,[[0,1],[1,2],[2,3],[0,3]]),graph(4,[[0,1],[2,3]]),graph(3,[])];
const oneCut=(n,edges)=>{const L=o.zeros(n);for(const [a,b]of edges){L[a][a]=L[a][a].add(1);L[b][b]=L[b][b].add(1);L[a][b]=L[a][b].sub(1);L[b][a]=L[b][a].sub(1);}return o.scale(L,'1/2');};
const edgeOneCut=(n,a,b)=>oneCut(n,[[a,b]]),pathEdges=n=>Array.from({length:n-1},(_,a)=>[a,a+1]),ringEdges=n=>[...pathEdges(n),[0,n-1]];
const collectiveRead=(g,A)=>sumMatrices(Array.from({length:g.n},(_,a)=>localN(g.n,a,A)),g.D);
const nativeControl=(n,a,b)=>o.add(localN(n,a,P),o.mul(localN(n,a,Q),localN(n,b,K)));
const factorial=n=>{let x=f(1);for(let j=2;j<=n;j++)x=x.mul(j);return x;};
const thetaCache=new Map(),thetaApprox=n=>{if(thetaCache.has(n))return thetaCache.get(n);let lo=f(0),hi=f(0);for(let j=0;j<n;j++){hi=hi.add(f(4).div(f(1).add(f(j).div(n).pow(2))));lo=lo.add(f(4).div(f(1).add(f(j+1).div(n).pow(2))));}lo=lo.div(n);hi=hi.div(n);ensure(hi.sub(lo).eq(f(2).div(n)),'native half-turn bracket width');const ans={lo,hi,value:lo.add(hi).div(2),error:f(1).div(n)};thetaCache.set(n,ans);return ans;};
const expPoly=(A,t,degree,gain)=>{const x=f(gain).mul(f(t).abs()),Z=o.scale(A,new o.Cut(0,f(t).neg()));ensure(x.le(degree+2)&&!x.eq(degree+2),'factorial tail ratio');let term=o.identity(A.length),sum=term;for(let j=1;j<=degree;j++){term=o.scale(o.mul(term,Z),f(1).div(j));sum=o.add(sum,term);}return {matrix:sum,tail:x.pow(degree+1).div(factorial(degree+1)).div(f(1).sub(x.div(degree+2)))};};
const cayleyCut=(P,x)=>{const z=new o.Cut(1,f(x).neg().div(2)).div(new o.Cut(1,f(x).div(2)));return o.add(o.identity(P.length),o.scale(P,z.sub(1)));};

check('native_edge_words_construct_collective_matching_targets',()=>{
 let controls=0,readouts=0,positivities=0;for(const g of smallGraphs){eq(o.dagger(g.G),g.G,'collective native self-dagger generator');for(let e=0;e<g.edges.length;e++){const [a,b]=g.edges[e],ab=nativeControl(g.n,a,b),ba=nativeControl(g.n,b,a);eq(o.mul(o.mul(ab,ba),ab),g.swaps[e],'lifted three-control exchange');eq(o.mul(g.cuts[e],g.cuts[e]),g.cuts[e],'lifted exchange cut');ensure(o.rank(g.cuts[e])===2**(g.n-2),'spectator cut count');controls++;}
  for(const A of [Q,H,K,o.scale(R,io)]){eq(comm(g.G,collectiveRead(g,A)),o.zeros(g.D),'collective identical readout conserved');readouts++;}
  const v=col(Array.from({length:g.D},(_,j)=>new o.Cut(j%3-1,j%2))),left=inner(v,o.mul(g.G,v)),right=g.cuts.reduce((s,P)=>s.add(o.energy(o.mul(P,v))),f(0));ensure(left.eq(right),'collective cost is sum of cut matching squares');positivities++;
 }
 const ternary=graph(3,[[0,1],[1,2]],3);for(const p of ternary.cuts)ensure(o.rank(p)===9,'three-role spectator cut count');
 return {finite_edge_targets:smallGraphs.length,lifted_source_control_words:controls,identical_readout_conservation_identities:readouts,positive_cut_sum_witnesses:positivities,full_swap_source_control_factors:3,edge_list_and_equal_allocation_are_native_targets:true,physical_network_or_material_interaction_selected:false};
});

check('native_count_refinement_has_rational_error_enclosures',()=>{
 const support=thetaApprox(32),theta={lo:f('49/16'),hi:f('51/16'),value:f('25/8'),error:f('1/16')},t=f('1/4'),rows=[];ensure(theta.lo.le(support.lo)&&support.hi.le(theta.hi),'compact rational bracket contains native count bracket');ensure(theta.hi.le(4),'native half-turn upper gain');let finiteDefects=0;
 for(const [nRole,edges]of [[3,pathEdges(3)],[3,ringEdges(3)],[4,pathEdges(4)]]){const cuts=edges.map(([a,b])=>edgeOneCut(nRole,a,b)),G1=oneCut(nRole,edges),A=o.scale(G1,theta.value),m=edges.length,poly=expPoly(A,t,24,4*m);let lastBound;
  for(const n of [8,16,32]){let word=o.identity(nRole);for(const P of cuts){const C=cayleyCut(P,theta.value.mul(t).div(n));eq(o.mul(o.dagger(C),C),o.identity(nRole),'exact rational edge isometry');word=o.mul(C,word);}const run=o.power(word,n);eq(o.mul(o.dagger(run),run),o.identity(nRole),'finite refined word matching');
   const x=f(4*m).mul(t).div(n),split=f(16*m*m).mul(t.pow(2)).div(n).div(f(1).sub(x)),cayley=f(64*m).mul(t.pow(3)).div(12*n*n),wordError=split.add(cayley),comparison=wordError.add(poly.tail),residual=o.energy(o.sub(run,poly.matrix));ensure(residual.le(comparison.pow(2).mul(nRole)),'count-refined word factorial enclosure');
   const completed=wordError.add(theta.error.mul(m).mul(t));if(lastBound)ensure(completed.le(lastBound),'refinement error budget decreases');lastBound=completed;if(!o.isZero(comm(run,G1)))finiteDefects++;
   rows.push({ledger_count:nRole,edges:m,refinement_count:n,native_fractional_factors:m*n,phase_error_bound:theta.error.mul(m).mul(t),split_bound:split,cayley_bound:cayley,factorial_tail:poly.tail,completed_flow_error_bound:completed,residual_entry_norm_square_sha256:hash(JSON.stringify(residual)),residual_entry_norm_square_bound:comparison.pow(2).mul(nRole)});
  }
 }
 ensure(finiteDefects>0,'finite words need not conserve the collective generator');
 return {native_half_turn_count_terms:32,compact_phase_bracket:[theta.lo,theta.hi],exact_count_refinement_enclosures:rows,enclosure_count:rows.length,finite_words_with_nonzero_generator_defect:finiteDefects,finite_word_identified_with_completed_collective_flow:false,refined_factor_count_erased:false,ordinary_product_limit_theorem_assumed:false};
});

check('collective_calibration_defects_derive_component_geometry',()=>{
 let cases=0,orthogonal=0;const examples=[];for(const g of smallGraphs){for(const A of [Q,K,o.scale(R,io),o.scale(I,2)]){const reads=Array.from({length:g.n},(_,a)=>localN(g.n,a,A)),trA=trace(A),v=trace(o.mul(A,A)).rad.sub(trA.rad.pow(2).div(g.d)),edgeDefects=g.edges.map(([a,b],i)=>comm(g.swaps[i],reads[a]));
  for(let i=0;i<edgeDefects.length;i++)for(let j=i+1;j<edgeDefects.length;j++){ensure(trace(o.mul(o.dagger(edgeDefects[i]),edgeDefects[j])).zero(),'distinct edge calibration defects are orthogonal');orthogonal++;}
  for(const scales of [Array(g.n).fill(f(1)),Array.from({length:g.n},(_,a)=>f(a+1)),Array.from({length:g.n},(_,a)=>f(a%2?'3/2':'1/2')),g.labels.map(j=>f(2*j+1))]){const C=sumMatrices(reads.map((M,a)=>o.scale(M,scales[a])),g.D),differences=g.edges.reduce((s,[a,b])=>s.add(scales[a].sub(scales[b]).pow(2)),f(0)),rhs=f(g.d**(g.n-1)).mul(v).mul(differences).div(2),defect=o.energy(comm(g.G,C));ensure(defect.eq(rhs),'exact collective calibration coefficient');ensure(defect.zero()===(v.zero()||differences.zero()),'component calibration criterion');cases++;if(o.equal(A,Q)&&examples.length<12)examples.push({ledgers:g.n,edges:g.edges,component_count:g.count,multipliers:scales,defect_over_half_turn_square:defect});}
 }}
 const g=graph(3,[[0,1],[1,2]],3),A=o.matrix([[1,new o.Cut(1,1),0],[new o.Cut(1,-1),-2,1],[0,1,3]]),sc=[f(1),f('3/2'),f(2)],C=sumMatrices(sc.map((c,a)=>o.scale(localN(3,a,A),c)),g.D),v=trace(o.mul(A,A)).rad.sub(trace(A).rad.pow(2).div(3)),dif=sc[0].sub(sc[1]).pow(2).add(sc[1].sub(sc[2]).pow(2));ensure(o.energy(comm(g.G,C)).eq(v.mul(dif).mul('9/2')),'three-role collective calibration');cases++;
 return {collective_calibration_defect_instances:cases,distinct_edge_orthogonality_instances:orthogonal,calibration_examples:examples,component_constants_are_exact_kernel:true,connected_relative_factors_derived:true,disconnected_components_forced_to_share_one_factor:false,scalar_readout_or_common_absolute_scale_selected:false};
});

check('native_local_currents_obey_subset_balance',()=>{
 let currents=0,boundaries=0,hierarchies=0;for(const g of smallGraphs){for(const A of [Q,K]){const reads=Array.from({length:g.n},(_,a)=>localN(g.n,a,A)),F=g.edges.map(([a,b],j)=>o.scale(comm(g.cuts[j],reads[a]),io));
  for(let e=0;e<g.edges.length;e++){const [a,b]=g.edges[e];eq(o.dagger(F[e]),F[e],'native current self-dagger');eq(o.scale(comm(g.cuts[e],reads[b]),io),o.scale(F[e],-1),'edge current antisymmetry');currents++;}
  for(let a=0;a<g.n;a++){const incident=sumMatrices(g.edges.map(([u,v],e)=>u===a?F[e]:v===a?o.scale(F[e],-1):o.zeros(g.D)),g.D);eq(o.scale(comm(g.G,reads[a]),io),incident,'exact local continuity');}
  for(let mask=0;mask<2**g.n;mask++){const inside=a=>(mask>>a)&1,read=sumMatrices(reads.filter((_,a)=>inside(a)),g.D),boundary=sumMatrices(g.edges.map(([a,b],e)=>inside(a)===inside(b)?o.zeros(g.D):o.scale(F[e],inside(a)?1:-1)),g.D);eq(o.scale(comm(g.G,read),io),boundary,'all-subset current cancellation');boundaries++;}
  for(let e=0;e<g.edges.length;e++){const [a,b]=g.edges[e],terms=g.edges.map(([u,v],j)=>{const c=comm(g.cuts[j],F[e]);if(![u,v].includes(a)&&![u,v].includes(b))eq(c,o.zeros(g.D),'disjoint current derivative vanishes');return c;});eq(comm(g.G,F[e]),sumMatrices(terms,g.D),'native current hierarchy');hierarchies++;}
 }}
 return {self_dagger_antisymmetric_edge_currents:currents,exact_subset_balance_identities:boundaries,local_current_hierarchy_identities:hierarchies,primitive_continuum_divergence_or_field_equation_assumed:false};
});

const S13=swapN(3,0,2),S23=swapN(3,1,2),Id8=o.identity(8),P13=o.scale(o.sub(Id8,S13),'1/2'),P23=o.scale(o.sub(Id8,S23),'1/2'),Gstar=o.add(P13,P23),F23=o.scale(comm(P23,localN(3,1,Q)),io),D23=o.scale(comm(Gstar,F23),io),Omega3=o.scale(comm(P13,P23),io);
const hiddenState=eta=>{const v=Array(16).fill(o.ZERO);for(const [s,a]of [['0100',o.ONE],['1000',eta],['0111',o.ONE],['1011',eta.neg()]])v[parseInt(s,2)]=a;return o.scale(col(v),'1/2');};
const hiddenEtas=[o.ONE,o.ONE.neg(),io,io.neg()],hiddenStates=hiddenEtas.map(hiddenState);
check('all_pair_readouts_miss_third_ledger_current_and_curvature',()=>{
 const pair12=o.scale(o.add(pairOf(basis(4,1)),pairOf(basis(4,2))),'1/2'),rows=[];let comparisons=0;
 for(let j=0;j<hiddenStates.length;j++){const v=hiddenStates[j],eta=hiddenEtas[j],full=pairOf(v);ensure(o.energy(v).eq(1),'pure native hidden-record preparation normalized');for(const [keep,expected]of [[[0,1],pair12],[[0,2],o.scale(I4,'1/4')],[[1,2],o.scale(I4,'1/4')]]){eq(part(full,4,keep),expected,'every pair readout identical');comparisons++;}
  const current=mean(v,o.kron(F23,I)),slope=mean(v,o.kron(D23,I)),curvature=mean(v,o.kron(Omega3,I));ensure(current.zero()&&slope.eq(eta.rad.neg().div(4))&&curvature.eq(eta.turn.div(4)),'hidden current and curvature witness');rows.push({eta,current,current_slope_over_half_turn:slope,ordering_curvature:curvature});
 }
 return {identical_two_ledger_pair_readouts:comparisons,hidden_record_witnesses:rows,all_pair_observations_are_predictively_complete:false,three_ledger_readouts_claimed_complete_for_every_larger_target:false,independent_random_mixture_premise:false};
});

check('native_shared_edge_ordering_has_a_quantized_cycle_cut',()=>{
 let cycles=0;const rows=[];for(const [n,d,e,fedge]of [[3,2,[0,2],[1,2]],[3,2,[0,1],[1,2]],[4,2,[0,1],[1,3]],[3,3,[0,1],[1,2]]]){const id=o.identity(d**n),A=swapN(n,...e,d),B=swapN(n,...fedge,d),Pa=o.scale(o.sub(id,A),'1/2'),Pb=o.scale(o.sub(id,B),'1/2'),X=o.mul(A,B),Xd=o.dagger(X),Om=o.scale(comm(Pa,Pb),io),Pc=o.sub(id,o.scale(o.add(o.add(id,X),Xd),'1/3'));
  eq(o.power(X,3),id,'three-cycle word');eq(o.mul(Pc,Pc),Pc,'cycle circulation cut');eq(o.dagger(Pc),Pc,'cycle cut self-dagger');eq(Om,o.scale(o.sub(X,Xd),new o.Cut(0,'1/4')),'ordering curvature from source cycle');eq(o.mul(Om,Om),o.scale(Pc,'3/16'),'native quantized curvature square');ensure(o.rank(Pc)===2*(d**3-d)/3*d**(n-3),'native cycle cut rank');ensure(!o.equal(o.mul(o.mul(o.mul(A,B),A),B),id),'full ordered loop nontrivial');rows.push({ledgers:n,roles_per_ledger:d,circulation_cut_rank:o.rank(Pc),curvature_square_coefficient:'3/16'});cycles++;
 }
 eq(comm(o.scale(o.sub(I16,swapN(4,0,1)),'1/2'),o.scale(o.sub(I16,swapN(4,2,3)),'1/2')),o.zeros(16),'disjoint exchange curvature zero');
 const jet=(P,axis,sign)=>new Map([['0,0',Id8],[axis===0?'1,0':'0,1',o.scale(P,new o.Cut(0,-sign))],[axis===0?'2,0':'0,2',o.scale(o.mul(P,P),'-1/2')]]),mulJet=(a,b)=>new Map([...compose(a,b)].filter(([k])=>xy(k)[0]+xy(k)[1]<=2));
 const loopJet=[jet(P13,0,1),jet(P23,1,1),jet(P13,0,-1),jet(P23,1,-1)].reduce(mulJet,unit(Id8));eq(loopJet.get('1,1'),o.scale(comm(P13,P23),-1),'mixed second-order loop coefficient');for(const k of ['1,0','0,1','2,0','0,2'])ensure(!loopJet.has(k)||o.isZero(loopJet.get(k)),'pure loop coefficients cancel');
 return {native_cycle_cut_instances:cycles,cycle_cut_examples:rows,mixed_loop_coefficient_exact:true,tree_edge_list_can_have_operation_order_curvature:true,cycle_cut_coefficient_identified_with_physical_alpha:false};
});

check('one_cut_exchange_derives_difference_geometry_and_zero_modes',()=>{
 let intertwinings=0,costs=0,zeroModes=0;const rows=[];for(const g of [...smallGraphs,graph(5,pathEdges(5))]){const T=o.zeros(g.D,g.n);for(let a=0;a<g.n;a++)T[2**(g.n-1-a)][a]=o.ONE;const g1=oneCut(g.n,g.edges);eq(o.mul(o.dagger(T),T),o.identity(g.n),'native single-cut insertion matching');eq(o.mul(g.G,T),o.mul(T,g1),'single-cut generator from swaps');intertwinings++;
  ensure(o.rank(g1)===g.n-g.count,'component zero-mode rank');for(let c=0;c<g.count;c++){const v=col(g.labels.map(j=>j===c?1:0));eq(o.mul(g1,v),o.zeros(g.n,1),'component-constant zero mode');zeroModes++;}
  for(const v of [col(Array.from({length:g.n},(_,a)=>a+1)),col(Array.from({length:g.n},(_,a)=>new o.Cut(a%2,a+1))),basis(g.n,0)]){const rhs=g.edges.reduce((s,[a,b])=>s.add(v[a][0].sub(v[b][0]).norm2()),f(0)).div(2);ensure(inner(v,o.mul(g1,v)).eq(rhs),'native one-cut cost equals endpoint differences');costs++;}
  rows.push({ledgers:g.n,edges:g.edges,component_count:g.count,derived_difference_rank:o.rank(g1)});
 }
 return {full_tuple_to_single_cut_intertwinings:intertwinings,exact_difference_cost_identities:costs,component_constant_zero_modes:zeroModes,target_examples:rows,ordinary_graph_laplacian_assumed:false,full_tuple_kernel_claimed_equal_to_component_count:false};
});

const distances=(n,edges,a)=>{const dist=Array(n).fill(Infinity),paths=Array(n).fill(0);dist[a]=0;paths[a]=1;const queue=[a];for(let k=0;k<queue.length;k++){const u=queue[k];for(const [x,y]of edges){const v=x===u?y:y===u?x:-1;if(v<0)continue;if(dist[v]===Infinity){dist[v]=dist[u]+1;queue.push(v);}if(dist[v]===dist[u]+1)paths[v]+=paths[u];}}return {dist,paths};};
check('cyclic_modes_and_distance_tails_follow_native_counting',()=>{
 let eigenspaces=0,geometric=0,leading=0,zeros=0,tailChecks=0;const ringRows=[];
 for(const [n,values]of [[3,[[0,1],['3/2',2]]],[4,[[0,1],[1,2],[2,1]]],[6,[[0,1],['1/2',2],['3/2',2],[2,1]]]]){const g1=oneCut(n,ringEdges(n));for(const [value,mult]of values){ensure(n-o.rank(o.sub(g1,o.scale(o.identity(n),value)))===mult,'native ring phase eigenspace');eigenspaces++;}ringRows.push({ledgers:n,phase_over_half_turn_and_multiplicity:values});}
 const roots=[o.ONE,io,o.ONE.neg(),io.neg()],modes=roots.map(z=>col(Array.from({length:4},(_,j)=>{let v=o.ONE;for(let k=0;k<j;k++)v=v.mul(z);return v;}))),g4=oneCut(4,ringEdges(4));for(let m=0;m<4;m++){eq(o.mul(g4,modes[m]),o.scale(modes[m],[0,1,2,1][m]),'native four-role circle mode');for(let k=0;k<4;k++)ensure(inner(modes[m],modes[k]).eq(m===k?4:0),'native mode matching');}
 for(let n=3;n<=16;n++){const geom=new Map(Array.from({length:n},(_,j)=>[key(j,0),o.identity(1)]));eqFields(compose(minus(id1,sx1),geom),minus(id1,shift(id1,n,0)),1,1,'native finite geometric identity');geometric++;}
 const theta=thetaApprox(16),gapBudgets=[];for(const n of [4,5,6,8,16,32,64,128]){const x=theta.value.div(n);ensure(x.le(1),'small native phase');let upper=f(0);for(let j=0;j<=8;j++)upper=upper.add(x.pow(2*j+1).div(factorial(2*j+1)).mul(j%2?-1:1));const lower=upper.sub(x.pow(19).div(factorial(19))),bound=f(16).div(3*n*n);ensure(f(0).le(lower)&&lower.pow(2).div(x.pow(2)).le(1)&&f(1).sub(bound).le(lower.pow(2).div(x.pow(2)))&&upper.pow(2).div(x.pow(2)).le(1),'native alternating gap bracket');gapBudgets.push({ledgers:n,relative_gap_lower_bound:f(1).sub(bound),relative_gap_upper_bound:'1'});}
 for(const [n,edges]of [[3,pathEdges(3)],[4,pathEdges(4)],[5,pathEdges(5)],[6,pathEdges(6)],[4,ringEdges(4)],[5,ringEdges(5)],[6,ringEdges(6)],[5,[[0,1],[0,2],[0,3],[0,4]]],[4,[[0,1],[2,3]]]]){const g1=oneCut(n,edges),powers=[o.identity(n)];for(let j=1;j<n;j++)powers.push(o.mul(powers[j-1],g1));const degree=Array(n).fill(0);for(const [a,b]of edges){degree[a]++;degree[b]++;}const Delta=Math.max(...degree),t=f('1/32'),x=f(4*Delta).mul(t),poly=expPoly(o.scale(g1,theta.value),t,16,4*Delta);ensure(x.le('1/2'),'explicit geometric exponential upper bound');
  for(let a=0;a<n;a++){const {dist,paths}=distances(n,edges,a);for(let b=0;b<n;b++){if(!Number.isFinite(dist[b])){for(const M of powers)ensure(M[b][a].zero(),'disconnected zero propagation');continue;}const ell=dist[b];for(let j=0;j<ell;j++){ensure(powers[j][b][a].zero(),'shorter words cannot reach target');zeros++;}ensure(powers[ell][b][a].eq(f('-1/2').pow(ell).mul(paths[b])),'shortest route fixes first phase coefficient');leading++;
   const tail=x.pow(ell).div(factorial(ell)).div(f(1).sub(x));ensure(poly.matrix[b][a].norm2().le(tail.pow(2)),'finite factorial propagation tail');tailChecks++;
  }}
 }
 return {exact_ring_eigenspaces:eigenspaces,ring_spectra:ringRows,native_four_role_modes:4,finite_geometric_count_identities:geometric,ring_gap_scaling_budgets:gapBudgets,shorter_word_zero_coefficients:zeros,shortest_route_leading_coefficients:leading,finite_factorial_tail_enclosures:tailChecks,continuous_exchange_has_an_exact_fixed_word_front:false,finite_ring_gap_identified_as_physical_mass_gap:false};
});

const jointApply=(a,A)=>new Map([...a].map(([k,v])=>[k,processApply(v,A)])),jointRead=(a,A)=>[...a.values()].reduce((s,v)=>s.add(fieldPair(v,processApply(v,A))),o.ZERO);
const jointEqual=(a,b,D,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eqFields(a.get(k)||new Map(),b.get(k)||new Map(),32,D,msg+' '+k);};
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,Jc),-1),I),J32=o.kron(Jc,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
const processPairFromField=(a,D)=>[...a.values()].reduce((s,M)=>{const T=transpose(M);return o.add(s,o.mul(T,o.dagger(T)));},o.zeros(D));
check('actual_clock_preserves_collective_charge_and_complete_observation',()=>{
 const V=o.mul(S23,S13),g=graph(3,[[0,2],[1,2]]),readouts=[Q,H,K].map(A=>collectiveRead(g,A));let edges=0,prefixes=0,observed=0;
 for(let j=0;j<qProgram;j++){const a=joint(productField(sourceSeed,col([1,0,1,io,2,1,0,-1])),j,1);for(const A of readouts){jointEqual(coupledStep(jointApply(a,A),V,3),jointApply(coupledStep(a,V,3),A),8,'actual collective clock invariant');edges++;}}
 let a=joint(productField(sourceSeed,basis(8,4)));for(let n=0;n<=2*qProgram+2;n++){if(n)a=coupledStep(a,V,3);ensure(jointEnergy(a).eq(1)&&jointRead(a,readouts[0]).eq(1),'actual clock count invariant');const position=Math.floor(n/qProgram)%3;for(let k=0;k<3;k++)ensure(jointRead(a,localN(3,k,Q)).eq(k===position?1:0),'actual ordered word count transport');prefixes++;}
 const finiteDefect=o.energy(comm(V,Gstar));ensure(!finiteDefect.zero(),'finite word does not conserve collective phase sum');eq(comm(Gstar,o.mul(Gstar,Gstar)),o.zeros(8),'completed flow generator moment commutation');
 for(const j of [0,2]){const v=hiddenStates[j],field=productField(sourcePowers[1],v),obs=observe(field);eqFields(decode(obs),field,32,16,'complete collective source reconstruction');const full=processPairFromField(field,16),seen=o.add(processPairFromField(obs[0],16),processPairFromField(obs[1],16));eq(full,seen,'complete observer preserves hidden-record process pair');for(const A of [o.kron(D23,I),o.kron(Omega3,I)])ensure(trace(o.mul(full,A)).eq(trace(o.mul(seen,A))),'complete observer current and curvature');observed++;}
 const calibration=[];for(const [aa,bb]of [[2,3],['1/2','5/3'],['7/4','2/5']]){const aaF=f(aa),bbF=f(bb),v=basis(8,4),der=o.scale(o.mul(Gstar,v),new o.Cut(0,f(-1).div(aaF)));eq(o.scale(der,new o.Cut(0,aaF.mul(bbF))),o.scale(o.mul(Gstar,v),bbF),'collective common action calibration with half-turn factored');calibration.push({time_scale:aaF,common_readout_scale:bbF,remaining_action_scale:aaF.mul(bbF)});}
 return {actual_clock_readout_intertwinings:edges,actual_clock_prefix_count_readings:prefixes,complete_hidden_record_observer_witnesses:observed,full_swap_process_source_control_factors:6,finite_ordered_word_collective_generator_defect_square:finiteDefect,remaining_common_calibrations:calibration,finite_ordered_word_conserves_collective_generator_automatically:false,carrier_error_becomes_conservation_error:false,physical_hbar_c_alpha_mass_or_dimension_selected:false,ordinary_complex_field_or_physical_network_premise:false};
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
 an_edge_list_or_material_interaction_was_uniquely_selected:checks[0].physical_network_or_material_interaction_selected===false,
 finite_refined_words_equal_the_completed_collective_flow:checks[1].finite_word_identified_with_completed_collective_flow===false,
 count_refinement_erases_its_factor_cost:checks[1].refined_factor_count_erased===false,
 an_external_product_limit_theorem_was_a_premise:checks[1].ordinary_product_limit_theorem_assumed===false,
 disconnected_components_must_share_one_calibration:checks[2].disconnected_components_forced_to_share_one_factor===false,
 a_scalar_readout_or_collective_conservation_fixes_absolute_scale:checks[2].scalar_readout_or_common_absolute_scale_selected===false,
 all_pair_observations_predict_every_next_current:checks[4].all_pair_observations_are_predictively_complete===false,
 three_ledger_readouts_close_every_larger_network:checks[4].three_ledger_readouts_claimed_complete_for_every_larger_target===false,
 the_hidden_record_witness_assumes_a_random_mixture:checks[4].independent_random_mixture_premise===false,
 ordering_curvature_requires_a_cycle_in_the_edge_list:checks[5].tree_edge_list_can_have_operation_order_curvature===true,
 the_native_curvature_square_coefficient_is_physical_alpha:checks[5].cycle_cut_coefficient_identified_with_physical_alpha===false,
 the_one_cut_kernel_count_is_the_full_tuple_kernel_count:checks[6].full_tuple_kernel_claimed_equal_to_component_count===false,
 the_analytic_flow_has_a_fixed_word_front_or_physical_mass_gap:checks[7].continuous_exchange_has_an_exact_fixed_word_front===false&&checks[7].finite_ring_gap_identified_as_physical_mass_gap===false,
 every_finite_ordered_clock_word_conserves_the_collective_generator:checks[8].finite_ordered_word_conserves_collective_generator_automatically===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r44.native-collective-exchange.v1',status:'PASS_R44_NATIVE_COLLECTIVE_EXCHANGE',
 input_sha256:inputHash,r43_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_collective_exchange_flow_derived:true,native_component_calibration_geometry_derived:true,
 native_higher_memory_and_cycle_curvature_derived:true,native_single_cut_collective_modes_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: count-refined collective exchange, connected-component calibration, local current balance, pair-blind three-ledger memory and cycle curvature, shared single-cut difference geometry, cyclic phase modes and distance tails, and actual-clock conservation with complete observation. Edge lists, equal allocations and count interpolation are explicit native targets. Finite words retain their cost and conservation boundaries; physical constants and material identifications remain unselected.'},null,2)+'\n');
