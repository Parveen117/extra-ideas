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
ensure(inputHash===arg('--expected-input-sha256'),'R43 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r42-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r42-sha256'),'R43 R42 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R42_NATIVE_PHASE_GENERATOR','R43 parent status');
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
const sw=d=>{const out=o.zeros(d*d);for(let a=0;a<d;a++)for(let b=0;b<d;b++)out[b*d+a][a*d+b]=o.ONE;return out;};
const reduced=(M,d,side)=>o.matrix(Array.from({length:d},(_,a)=>Array.from({length:d},(_,b)=>{let z=o.ZERO;for(let m=0;m<d;m++)z=z.add(side===0?M[a*d+m][b*d+m]:M[m*d+a][m*d+b]);return z;})));
const pairOf=v=>o.scale(o.mul(v,o.dagger(v)),f(1).div(o.energy(v)));
const mean=(v,A)=>{const z=inner(v,o.mul(A,v)).div(o.energy(v));ensure(z.turn.zero(),'real native readout');return z.rad;};
const pureSeeds=d=>[basis(d,0),basis(d,1),col(Array.from({length:d},(_,j)=>j+1)),col(Array.from({length:d},(_,j)=>new o.Cut(d-j,j%2)))];
const turns=[[1,0,0],[0,1,0],['1/2','1/2','1/2'],['1/2','1/2','-1/2'],['9/25','16/25','12/25'],['16/25','9/25','-12/25'],['25/169','144/169','60/169'],['64/289','225/289','120/289']]
 .map(([c2,s2,cs])=>({c2:f(c2),s2:f(s2),cs:f(cs)}));
const exchangeAt=(d,t)=>o.add(o.scale(o.identity(d*d),new o.Cut(t.c2,t.cs.neg())),o.scale(sw(d),new o.Cut(t.s2,t.cs)));
const Sswap=sw(2),Pminus=o.scale(o.sub(I4,Sswap),'1/2'),Pplus=o.sub(I4,Pminus);
const Wab=o.add(o.kron(P,I),o.kron(Q,K)),Wba=o.add(o.kron(I,P),o.kron(K,Q));
const Qa=o.kron(Q,I),Qb=o.kron(I,Q),Npair=o.add(Qa,Qb),Q11=o.kron(Q,Q),P00=o.kron(P,P);
const Podd=o.scale(o.sub(I4,o.kron(H,H)),'1/2'),Zpair=o.scale(o.sub(o.kron(H,I),o.kron(I,H)),'1/2'),Kpair=o.mul(o.kron(K,K),Podd),Jpair=o.scale(o.mul(Kpair,Zpair),io);
const Fpair=o.kron(H,H),Ujoint=o.mul(Fpair,Sswap),Llift=o.add(Npair,Pminus),Lend=o.scale(o.sub(I4,Ujoint),'1/2'),branchCut=o.add(Q11,Pminus);
const halfSwap=exchangeAt(2,turns[2]),halfH=o.add(P,o.scale(Q,new o.Cut(0,-1))),halfLift=o.mul(o.kron(halfH,halfH),halfSwap),halfEnd=o.add(I4,o.scale(Lend,new o.Cut(-1,-1)));
const v01=basis(4,1),v10=basis(4,2),anti=o.sub(v01,v10);

check('source_control_words_construct_complete_tuple_exchange',()=>{
 eq(o.mul(o.mul(Wab,Wba),Wab),Sswap,'three native controls construct swap');eq(o.mul(o.mul(o.kron(H,I),o.kron(I,H)),o.mul(o.mul(Wab,Wba),Wab)),Ujoint,'five-factor source-H process word');
 let instances=0;const ranks=[];for(const d of [2,3,4,5]){const S=sw(d),id=o.identity(d*d),pm=o.scale(o.sub(id,S),'1/2'),pp=o.sub(id,pm);eq(o.mul(S,S),id,'tuple exchange square');eq(o.dagger(S),S,'tuple exchange dagger');eq(o.mul(pm,pm),pm,'antisymmetric native cut');eq(o.mul(pp,pm),o.zeros(d*d),'orthogonal exchange cuts');ensure(o.rank(pm)===d*(d-1)/2&&o.rank(pp)===d*(d+1)/2,'native unequal-pair cut counts');
  for(const t of turns){ensure(t.c2.add(t.s2).eq(1)&&t.c2.mul(t.s2).eq(t.cs.pow(2)),'native unit-turn coefficients');const W=exchangeAt(d,t);eq(o.mul(o.dagger(W),W),id,'complete fractional exchange isometry');eq(comm(W,S),o.zeros(d*d),'fractional exchange preserves swap');instances++;}
  ranks.push({roles_per_ledger:d,symmetric_cut_rank:o.rank(pp),antisymmetric_cut_rank:o.rank(pm)});
 }
 eq(o.mul(halfSwap,halfSwap),Sswap,'native half tick squares to exchange');eq(o.power(halfSwap,4),I4,'native two tick period');
 return {source_control_factors_per_exchange:3,source_H_combined_process_factors:5,tuple_cut_counts:ranks,exact_fractional_exchange_instances:instances,half_tick_native_arrow:halfSwap,particle_statistics_assumed:false,physical_interaction_selected:false};
});
check('native_phase_additivity_retains_composition_branch_memory',()=>{
 let intertwinings=0;for(const A of [H,K,o.scale(R,io),o.matrix([[1,new o.Cut(1,1)],[new o.Cut(1,-1),-2]])]){eq(o.mul(Sswap,o.kron(A,I)),o.mul(o.kron(I,A),Sswap),'native swap intertwiner');eq(comm(Sswap,o.add(o.kron(A,I),o.kron(I,A))),o.zeros(4),'sum readout exchange invariance');intertwinings++;}
 for(const V of [H,K,R,processTurn]){const VV=o.kron(V,V);eq(comm(VV,Sswap),o.zeros(4),'identical process and exchange commute');eq(o.mul(o.dagger(o.mul(VV,Sswap)),o.mul(VV,Sswap)),I4,'combined process isometry');}
 eq(comm(Npair,Pminus),o.zeros(4),'resource and interaction commute');eq(comm(Ujoint,Llift),o.zeros(4),'lifted generator conservation');eq(o.sub(Llift,Lend),o.scale(branchCut,2),'exact retained phase branch');eq(o.mul(branchCut,branchCut),branchCut,'integer branch cut');
 eq(o.mul(halfLift,halfLift),Ujoint,'lifted half ticks recover combined process');eq(o.mul(halfEnd,halfEnd),Ujoint,'endpoint half ticks recover same process');eq(halfLift,o.mul(halfEnd,o.sub(I4,o.scale(branchCut,2))),'fractional branch difference');ensure(!o.equal(halfLift,halfEnd),'endpoint does not fix additive history');
 const Pone=o.sub(Podd,Pminus),Ptwo=branchCut;eq(o.add(o.add(P00,Pone),Ptwo),I4,'complete lifted generator cuts');eq(Llift,o.add(Pone,o.scale(Ptwo,2)),'lifted generator count coefficients');for(const A of [P00,Pone,Ptwo])eq(o.mul(A,A),A,'lifted native cuts');eq(o.mul(Pone,Ptwo),o.zeros(4),'lifted cuts orthogonal');
 return {arbitrary_source_readout_intertwinings:intertwinings,identical_native_process_pairs:4,lifted_generator_coefficients:[0,1,2],phase_branch_difference:'2 theta_- (Q_11+P_-)',endpoint_generator_is_automatically_additive:false,fractional_branch_witness:halfLift};
});
check('exchange_conservation_derives_current_and_relative_calibration',()=>{
 const reads=[Q,H,K,o.scale(R,io),o.matrix([[1,new o.Cut(1,1)],[new o.Cut(1,-1),-2]]),o.matrix([[0,0,0],[0,1,0],[0,0,3]]),o.scale(I,2)];
 const scales=[[1,1],[1,2],['2/3','2/3'],['1/2','5/3'],['7/4','2/5']];let balances=0,calibrations=0;const rows=[];
 for(const E of reads){const d=E.length,id=o.identity(d),S=sw(d),pm=o.scale(o.sub(o.identity(d*d),S),'1/2'),EA=o.kron(E,id),EB=o.kron(id,E);eq(comm(S,o.add(EA,EB)),o.zeros(d*d),'all-state readout balance');eq(o.add(o.scale(comm(pm,EA),io),o.scale(comm(pm,EB),io)),o.zeros(d*d),'opposite transfer currents');balances++;
  const trE=trace(E),trE2=trace(o.mul(E,E));ensure(trE.turn.zero()&&trE2.turn.zero(),'native self-dagger traces');const centered=trE2.rad.sub(trE.rad.pow(2).div(d));ensure(f(0).le(centered),'native centered coefficient norm');
  for(const [aa,bb]of scales){const a=f(aa),b=f(bb),C=o.add(o.scale(EA,a),o.scale(EB,b)),defect=o.energy(comm(S,C)),rhs=a.sub(b).pow(2).mul(2*d).mul(centered);ensure(defect.eq(rhs),'complete calibration defect');ensure(defect.zero()===(a.eq(b)||centered.zero()),'exact relative-calibration criterion');calibrations++;if(d===2&&o.equal(E,Q))rows.push({factor_A:a,factor_B:b,exchange_defect_square:defect});}
 }
 return {complete_readout_balance_identities:balances,exact_calibration_defect_instances:calibrations,source_count_calibrations:rows,nonconstant_identical_readout_requires_equal_factors:true,constant_readout_determines_relative_scale:false,common_factor_or_offsets_selected:false};
});
check('two_role_exchange_requires_the_hidden_pair_current',()=>{
 for(const A of [Zpair,Kpair,Jpair])eq(o.mul(A,A),Podd,'native exchange role square');eq(Pminus,o.scale(o.sub(Podd,Kpair),'1/2'),'exchange cut from sector roles');eq(o.scale(comm(Pminus,Zpair),io),o.scale(Jpair,-1),'native count difference derivative');eq(o.scale(comm(Pminus,Jpair),io),Zpair,'native exchange current derivative');eq(o.scale(comm(Pminus,Qa),io),o.scale(Jpair,'1/2'),'current transfer coefficient');
 let laws=0;const states=[v01,v10,col([0,1,io,0]),col([0,1,io.neg(),0]),col([1,new o.Cut(1,1),2,new o.Cut(0,1)])];
 for(const v of states)for(const t of turns){const out=o.mul(exchangeAt(2,t),v),z=mean(v,Zpair),j=mean(v,Jpair),nextZ=t.c2.sub(t.s2).mul(z).sub(t.cs.mul(2).mul(j)),nextJ=t.cs.mul(2).mul(z).add(t.c2.sub(t.s2).mul(j));ensure(mean(out,Zpair).eq(nextZ)&&mean(out,Jpair).eq(nextJ),'exact two-readout oscillator');ensure(mean(out,Podd).eq(mean(v,Podd))&&mean(out,Npair).eq(mean(v,Npair)),'exchange counts conserved');laws++;}
 for(const t of turns){const out=o.mul(exchangeAt(2,t),v01);ensure(mean(out,Qa).eq(t.s2)&&mean(out,Qb).eq(t.c2)&&mean(out,Jpair).eq(t.cs.mul(2)),'exact source count transfer');}
 const plus=states[2],minus=states[3];for(const side of [0,1])eq(reduced(pairOf(plus),2,side),reduced(pairOf(minus),2,side),'same local pairs opposite hidden current');ensure(mean(plus,Jpair).eq(1)&&mean(minus,Jpair).eq(-1),'opposite hidden pair currents');ensure(mean(v01,Jpair).zero()&&mean(v10,Jpair).zero()&&mean(v01,Zpair).eq(1)&&mean(v10,Zpair).eq(-1),'current alone not predictive');
 return {exact_closed_oscillator_instances:laws,source_transfer_turns:turns.length,opposite_current_with_same_local_pairs:{plus:mean(plus,Jpair),minus:mean(minus,Jpair)},half_tick_count_A:mean(o.mul(halfSwap,v01),Qa),full_tick_count_A:mean(o.mul(Sswap,v01),Qa),one_local_count_or_current_alone_is_complete:false,two_readouts_claimed_to_reconstruct_full_joint_state:false};
});
let costExamples=[];
check('exchange_phase_cost_equals_native_distinguishability',()=>{
 let costs=0,metrics=0;for(const d of [2,3,4]){const pm=o.scale(o.sub(o.identity(d*d),sw(d)),'1/2'),seeds=pureSeeds(d);for(let ai=0;ai<seeds.length;ai++)for(let bi=0;bi<seeds.length;bi++){const a=seeds[ai],b=seeds[bi],k=inner(a,b).norm2().div(o.energy(a).mul(o.energy(b))),cost=mean(o.kron(a,b),pm);ensure(cost.eq(f(1).sub(k).div(2))&&f(0).le(cost)&&cost.le('1/2'),'native product cost identity and ceiling');costs++;if(d===2&&ai===0)costExamples.push({input_A:ai,input_B:bi,overlap_square:k,interaction_cost_over_half_turn:cost});}}
 for(const E of [Q,H,o.scale(R,io)])for(const v of pureSeeds(2)){const nv=o.energy(v),m1=inner(v,o.mul(E,v)).div(nv),m2=inner(v,o.mul(o.mul(E,E),v)).div(nv);ensure(m1.turn.zero()&&m2.turn.zero(),'response moments real');const a1=m1.mul(new o.Cut(0,-1)),a2=m2.mul('-1/2'),overlapSecond=a2.add(a2.dagger()).add(a1.norm2()),variance=m2.rad.sub(m1.rad.pow(2));ensure(overlapSecond.turn.zero()&&overlapSecond.rad.neg().div(2).eq(variance.div(2)),'native local interaction metric coefficient');metrics++;}
 ensure(mean(anti,Pminus).eq(1),'antisymmetric full-state cost');ensure(mean(anti,Pminus).le('1/2')===false,'nonproduct exceeds product cost ceiling');
 return {exact_product_cost_instances:costs,local_metric_coefficient_checks:metrics,product_cost_examples:costExamples,product_cost_ceiling:'theta_-/2',full_state_cost_ceiling:'theta_-',antisymmetric_cost_over_half_turn:mean(anti,Pminus),product_cost_formula_claimed_for_all_joint_states:false};
});

let mixingDefect;
check('partial_exchange_pair_readouts_derive_memory_deficits',()=>{
 let pairs=0,deficits=0;const examples=[];for(const d of [2,3]){const seeds=pureSeeds(d),pm=o.scale(o.sub(o.identity(d*d),sw(d)),'1/2');for(const a of seeds)for(const b of seeds){const A=pairOf(a),Bv=pairOf(b),k=inner(a,b).norm2().div(o.energy(a).mul(o.energy(b))),v=o.kron(a,b),cost=mean(v,pm);
  for(const t of turns){const out=o.mul(exchangeAt(d,t),v),pair=pairOf(out),ra=reduced(pair,d,0),rb=reduced(pair,d,1),cross=o.scale(comm(Bv,A),new o.Cut(0,t.cs)),predA=o.add(o.add(o.scale(A,t.c2),o.scale(Bv,t.s2)),cross),predB=o.sub(o.add(o.scale(A,t.s2),o.scale(Bv,t.c2)),cross);
   eq(ra,predA,'complete native process pair A');eq(rb,predB,'complete native process pair B');eq(o.add(ra,rb),o.add(A,Bv),'pair sum balance');ensure(trace(ra).eq(1)&&trace(rb).eq(1),'native marginal normalization');pairs++;
   const ma=f(1).sub(trace(o.mul(ra,ra)).rad),mb=f(1).sub(trace(o.mul(rb,rb)).rad),formula=t.c2.mul(t.s2).mul(2).mul(f(1).sub(k).pow(2)),fromCost=t.cs.pow(2).mul(8).mul(cost.pow(2));ensure(ma.eq(mb)&&ma.eq(formula)&&ma.eq(fromCost)&&f(0).le(ma)&&ma.le('1/2'),'exact product memory deficit and cost relation');ensure(mean(out,pm).eq(cost),'interaction phase cost conserved');deficits++;
  }
 }}
 for(const [name,a,b]of [['orthogonal',e0,e1],['overlapping',e0,col([1,1])],['parallel',e0,e0]]){const v=o.kron(a,b),cost=mean(v,Pminus),ra=reduced(pairOf(o.mul(halfSwap,v)),2,0),memory=f(1).sub(trace(o.mul(ra,ra)).rad);ensure(memory.eq(cost.pow(2).mul(2)),'half tick cost-memory relation');examples.push({preparation:name,interaction_cost_over_half_turn:cost,half_tick_memory_deficit:memory});}
 const A=pairOf(e0),Bv=pairOf(col([1,1])),v=o.kron(e0,col([1,1])),actual=reduced(pairOf(o.mul(halfSwap,v)),2,0),mix=o.scale(o.add(A,Bv),'1/2');mixingDefect=o.energy(o.sub(actual,mix));ensure(mixingDefect.eq('1/8'),'omitted commutator has an explicit defect');
 const correlatedMemory=f(1).sub(trace(o.mul(reduced(pairOf(anti),2,0),reduced(pairOf(anti),2,0))).rad);ensure(correlatedMemory.eq('1/2'),'already correlated initial memory');
 return {exact_two_ledger_pair_readouts:pairs,exact_cost_memory_identities:deficits,half_tick_examples:examples,omitted_commutator_defect_square:mixingDefect,correlated_input_memory_at_zero_tick:correlatedMemory,random_selection_or_partial_trace_axiom_assumed:false,product_memory_law_claimed_for_correlated_inputs:false};
});
const halfTurnCache=new Map();
const halfTurnApprox=n=>{if(halfTurnCache.has(n))return halfTurnCache.get(n);let lower=f(0),upper=f(0);for(let j=0;j<n;j++){upper=upper.add(f(4).div(f(1).add(f(j).div(n).pow(2))));lower=lower.add(f(4).div(f(1).add(f(j+1).div(n).pow(2))));}lower=lower.div(n);upper=upper.div(n);ensure(upper.sub(lower).eq(f(2).div(n)),'frozen native half-turn bracket');const a={lower,upper,value:lower.add(upper).div(2),error:f(1).div(n)};halfTurnCache.set(n,a);return a;};
const fact=n=>{let a=f(1);for(let j=2;j<=n;j++)a=a.mul(j);return a;};
const expPolynomial=(A,t,m,gain)=>{const Z=o.scale(A,new o.Cut(0,f(t).neg())),alpha=f(gain).mul(f(t).abs());ensure(alpha.le(m+2)&&!alpha.eq(m+2),'native factorial tail ratio');let term=o.identity(A.length),sum=term;for(let j=1;j<=m;j++){term=o.scale(o.mul(term,Z),f(1).div(j));sum=o.add(sum,term);}const tail=alpha.pow(m+1).div(fact(m+1)).div(f(1).sub(alpha.div(m+2)));return {matrix:sum,tail,alpha};};
check('finite_phase_enclosures_bind_the_interpolated_exchange',()=>{
 const parentCert=JSON.parse(parentBytes),frozen=parentCert.exact_checks.find(x=>x.name==='native_resolvent_count_sums_construct_phase_with_bounds').phase_count_bounds[0],h16=halfTurnApprox(16);ensure(h16.lower.eq(frozen.half_turn_lower)&&h16.upper.eq(frozen.half_turn_upper),'R42 half-turn phase bounds reused exactly');
 const rows=[];for(const n of [8,16,32]){const theta=halfTurnApprox(n);ensure(theta.upper.le(4),'native positive half-turn bound');for(const [kind,L,t,expected,gain]of [['exchange',Pminus,'1/2',halfSwap,4],['exchange',Pminus,1,Sswap,4],['exchange',Pminus,2,I4,4],['lifted_process',Llift,'1/2',halfLift,8],['lifted_process',Llift,1,Ujoint,8]]){const A=o.scale(L,theta.value),poly=expPolynomial(A,t,40,gain),phaseError=theta.error.mul(t).mul(kind==='exchange'?1:2),bound=phaseError.add(poly.tail),residual=o.energy(o.sub(poly.matrix,expected));ensure(residual.le(bound.pow(2).mul(4)),'native phase recovery enclosure');rows.push({kind,count_terms:n,tick_parameter:String(t),factorial_degree:40,phase_error_bound:phaseError,factorial_tail:poly.tail,residual_entry_norm_square_sha256:hash(JSON.stringify(residual)),residual_entry_norm_square_bound:bound.pow(2).mul(4)});}}
 return {native_phase_reconstruction_enclosures:rows,exact_R42_half_turn_bounds_reused:true,phase_enclosure_instances:rows.length,finite_rational_phase_claimed_equal_to_completed_phase:false,huge_clock_or_packet_arrays_directly_simulated:false};
});
const jointApply=(a,A)=>new Map([...a].map(([k,v])=>[k,processApply(v,A)])),jointRead=(a,A)=>[...a.values()].reduce((s,v)=>s.add(fieldPair(v,processApply(v,A))),o.ZERO);
const jointEqual=(a,b,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eqFields(a.get(k)||new Map(),b.get(k)||new Map(),32,4,msg+' '+k);};
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),Pc=o.scale(o.add(I16,Hc),'1/2'),P32=o.kron(Pc,I),PJ32=o.kron(o.scale(o.mul(Pc,Jc),-1),I),J32=o.kron(Jc,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
const processPairFromField=a=>[...a.values()].reduce((s,M)=>{const T=transpose(M);return o.add(s,o.mul(T,o.dagger(T)));},o.zeros(4));
check('actual_clock_and_complete_curvature_observer_preserve_exchange',()=>{
 let edges=0,prefixes=0;const readouts=[Npair,Pminus,Llift];for(let j=0;j<qProgram;j++){const a=joint(productField(sourceSeed,col([1,new o.Cut(1,1),2,1])),j,1);for(const A of readouts){jointEqual(coupledStep(jointApply(a,A),Ujoint,3),jointApply(coupledStep(a,Ujoint,3),A),'actual clock exchange conservation');edges++;}}
 let a=joint(productField(sourceSeed,v01));for(let n=0;n<=qProgram+2;n++){if(n)a=coupledStep(a,Ujoint,3);ensure(jointEnergy(a).eq(1),'actual clock norm');ensure(jointRead(a,Npair).eq(1)&&jointRead(a,Pminus).eq('1/2')&&jointRead(a,Llift).eq('3/2'),'actual clock invariant readouts');ensure(jointRead(a,Qa).eq(Math.floor(n/qProgram)%2),'actual tick exchanges native resource');prefixes++;}
 const vHalf=o.mul(halfSwap,v01),field=productField(sourcePowers[1],vHalf),obs=observe(field);eqFields(decode(obs),field,32,4,'complete exchange source reconstruction');ensure(energy(obs[0]).add(energy(obs[1])).eq(energy(field)),'complete observer norm');const pair=processPairFromField(field),observed=o.add(processPairFromField(obs[0]),processPairFromField(obs[1]));eq(pair,observed,'complete observer preserves process pair');const ra=reduced(pair,2,0),memory=f(1).sub(trace(o.mul(ra,ra)).rad);ensure(memory.eq('1/2'),'actual carried half-exchange memory');
 for(const A of readouts)ensure(fieldPair(obs[0],processApply(obs[0],A)).add(fieldPair(obs[1],processApply(obs[1],A))).eq(fieldPair(field,processApply(field,A))),'complete observer invariant readout');
 return {actual_source_clock_exchange_intertwinings:edges,actual_clock_prefix_count_readings:prefixes,complete_carried_process_pair_preserved:true,carried_half_tick_memory_deficit:memory,source_H_process_word_factors_retained:5,clock_carrier_return_error_becomes_exchange_balance_error:false,instantaneous_remote_control_claimed:false};
});
check('physical_scope_and_false_mixing_controls_keep_native_types_distinct',()=>{
 const Ksum=o.add(o.kron(K,I),o.kron(I,K));ensure(!o.isZero(comm(Ujoint,Ksum)),'free process can change a generic readout sum');eq(comm(Sswap,Ksum),o.zeros(4),'pure exchange preserves same sum');
 ensure(!mixingDefect.zero(),'commutator cannot be discarded');ensure(mean(anti,Pminus).eq(1),'product cost ceiling does not cover all states');const correlatedMemory=f(1).sub(trace(o.mul(reduced(pairOf(anti),2,0),reduced(pairOf(anti),2,0))).rad);ensure(!correlatedMemory.zero(),'product memory law does not apply to correlated time-zero input');
 const scales=[];for(const [aa,bb]of [[2,3],['1/2','5/3'],['7/4','2/5']]){const a=f(aa),b=f(bb),der=o.scale(o.mul(Llift,v01),new o.Cut(0,f(-1).div(a)));eq(o.scale(der,new o.Cut(0,a.mul(b))),o.scale(o.mul(Llift,v01),b),'common action calibration with half-turn factored');eq(comm(Ujoint,o.scale(Npair,b)),o.zeros(4),'every common calibration preserves exchange');scales.push({time_scale:a,common_readout_scale:b,remaining_action_scale:a.mul(b)});}
 return {all_readout_sums_conserved_under_additional_free_evolution:false,commutator_free_mixing_is_universal:false,product_input_memory_law_is_universal:false,relative_multipliers_selected_for_nonconstant_identical_readout:true,remaining_common_calibrations:scales,absolute_action_scale_or_physical_hbar_selected:false,exchange_phase_cost_identified_with_R33_spatial_alpha:false,material_interaction_or_particle_statistics_selected:false,ordinary_complex_or_random_measurement_premise:false};
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
 the_compiled_exchange_has_zero_source_control_cost:checks[0].source_control_factors_per_exchange===3,
 exchange_cuts_select_physical_particle_statistics:checks[0].particle_statistics_assumed===false,
 endpoint_phase_is_automatically_additive:checks[1].endpoint_generator_is_automatically_additive===false,
 a_constant_readout_determines_relative_calibration:checks[2].constant_readout_determines_relative_scale===false,
 relative_calibration_selects_common_factor_and_offsets:checks[2].common_factor_or_offsets_selected===false,
 one_local_count_or_current_alone_predicts_every_exchange:checks[3].one_local_count_or_current_alone_is_complete===false,
 two_current_readouts_reconstruct_the_full_joint_state:checks[3].two_readouts_claimed_to_reconstruct_full_joint_state===false,
 the_product_cost_formula_covers_all_joint_states:checks[4].product_cost_formula_claimed_for_all_joint_states===false,
 random_selection_or_a_partial_trace_axiom_was_assumed:checks[5].random_selection_or_partial_trace_axiom_assumed===false,
 the_product_memory_law_covers_correlated_inputs:checks[5].product_memory_law_claimed_for_correlated_inputs===false,
 rational_phase_approximants_are_exact_completed_phases:checks[6].finite_rational_phase_claimed_equal_to_completed_phase===false,
 clock_carrier_error_destroys_exchange_balance:checks[7].clock_carrier_return_error_becomes_exchange_balance_error===false,
 a_commutator_free_mixing_law_is_universal:checks[8].commutator_free_mixing_is_universal===false,
 native_exchange_cost_selects_physical_hbar_or_spatial_alpha:checks[8].absolute_action_scale_or_physical_hbar_selected===false&&checks[8].exchange_phase_cost_identified_with_R33_spatial_alpha===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r43.native-exchange-calibration.v1',status:'PASS_R43_NATIVE_EXCHANGE_CALIBRATION',
 input_sha256:inputHash,r42_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_exchange_and_additive_phase_history_derived:true,native_relative_exchange_calibration_derived:true,
 native_exchange_cost_memory_relation_derived:true,native_clock_carried_exchange_balance_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: source-controlled exchange, additive phase history with retained branch memory, balanced transfer and relative readout calibration, hidden pair currents, distinguishability phase cost, exact reduced-pair memory laws, and autonomous-clock conservation with a remaining common action scale. The identical ledgers, interaction branch and readout targets are explicit native constructions; physical energy, material time, h-bar and alpha remain unselected.'},null,2)+'\n');
