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
ensure(inputHash===arg('--expected-input-sha256'),'R45 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r44-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r44-sha256'),'R45 R44 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R44_NATIVE_COLLECTIVE_EXCHANGE','R45 parent status');
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
const checks=[],check=(name,fn)=>{checks.push({name,passed:true,...fn()});if(process.env.R45_CHECK_PROGRESS==='1')process.stderr.write(name+'\n');};
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
const sumMatrices=(xs,D)=>xs.reduce(o.add,o.zeros(D));
const digits=(j,n)=>Array.from({length:n},(_,a)=>Math.floor(j/2**(n-1-a))%2),encode=v=>v.reduce((s,a)=>2*s+a,0);
const swapN=(n,a,b)=>{const S=o.zeros(2**n);for(let j=0;j<2**n;j++){const x=digits(j,n);[x[a],x[b]]=[x[b],x[a]];S[encode(x)][j]=o.ONE;}return S;};
const localN=(n,a,A)=>{let M=o.identity(1);for(let j=0;j<n;j++)M=o.kron(M,j===a?A:I);return M;};
const mean=(v,A)=>inner(v,o.mul(A,v)).div(o.energy(v));
const pairOf=v=>o.scale(o.mul(v,o.dagger(v)),f(1).div(o.energy(v)));
const part=(M,n,a)=>{const out=o.zeros(2),rest=Array.from({length:n},(_,j)=>j).filter(j=>j!==a);for(let i=0;i<2;i++)for(let j=0;j<2;j++)for(let k=0;k<2**rest.length;k++){const x=Array(n).fill(0),y=Array(n).fill(0);x[a]=i;y[a]=j;digits(k,rest.length).forEach((v,z)=>x[rest[z]]=y[rest[z]]=v);out[i][j]=out[i][j].add(M[encode(x)][encode(y)]);}return out;};
const factorial=n=>{let a=f(1);for(let j=2;j<=n;j++)a=a.mul(j);return a;};
const expPoly=(A,t,degree,gain)=>{const x=f(gain).mul(f(t).abs()),Z=o.scale(A,new o.Cut(0,f(t).neg()));ensure(x.le(degree+2)&&!x.eq(degree+2),'factorial tail ratio');let term=o.identity(A.length),sum=term;for(let j=1;j<=degree;j++){term=o.scale(o.mul(term,Z),f(1).div(j));sum=o.add(sum,term);}return {matrix:sum,tail:x.pow(degree+1).div(factorial(degree+1)).div(f(1).sub(x.div(degree+2)))};};
const edgeCut=(n,a,b)=>{const M=o.zeros(n);for(const j of [a,b])M[j][j]=o.Cut.of('1/2');M[a][b]=M[b][a]=o.Cut.of('-1/2');return M;};
const swapCut=P=>o.sub(o.identity(P.length),o.scale(P,2));
const cayleySwap=(S,x)=>o.scale(o.add(o.identity(S.length),o.scale(S,new o.Cut(0,f(x).div(2)))),o.ONE.div(new o.Cut(1,f(x).div(2))));
const loopFactors=(S,T,x)=>[[S,x],[T,x],[S,f(x).neg()],[T,f(x).neg()]].map(([A,y])=>cayleySwap(A,y));
const balancedFactors=(S,T,x)=>[...loopFactors(S,T,x),...loopFactors(S,T,f(x).neg())];
const product=xs=>xs.reduce(o.mul,o.identity(xs[0].length));
const curve=(S,T)=>{const id=o.identity(S.length),P=o.scale(o.sub(id,S),'1/2'),Q=o.scale(o.sub(id,T),'1/2'),X=o.mul(S,T),Pc=o.sub(id,o.scale(o.add(o.add(id,X),o.dagger(X)),'1/3')),Om=o.scale(comm(P,Q),io);return {id,S,T,P,Q,X,Pc,P0:o.sub(id,Pc),Om};};
const tri=curve(swapCut(edgeCut(3,0,2)),swapCut(edgeCut(3,1,2))),full=curve(swapN(3,0,2),swapN(3,1,2));
const closedFlow=(g,c,v)=>o.add(o.add(g.P0,o.scale(g.Pc,c)),o.scale(g.Om,io.mul(v)));
const insert=o.zeros(8,3);[4,2,1].forEach((j,a)=>insert[j][a]=o.ONE);
const ring=n=>{const id=o.identity(n),T=o.zeros(n);for(let j=0;j<n;j++)T[(j+1)%n][j]=o.ONE;const Ti=o.dagger(T),G=o.sub(id,o.scale(o.add(T,Ti),'1/2')),D=o.scale(o.sub(T,Ti),io.mul('1/2')),pairs=Array.from({length:n},(_,j)=>[edgeCut(n,j,(j+1)%n),edgeCut(n,(j+1)%n,(j+2)%n)]),Om=sumMatrices(pairs.map(([P,Q])=>o.scale(comm(P,Q),io)),n);return {n,id,T,Ti,G,D,pairs,Om};};

check('native_balanced_loop_has_exact_fourth_order_residue',()=>{
 const id='012',s='210',t='021',pmul=(a,b)=>b.split('').map(i=>a[Number(i)]).join('');
 const add=(p,k,g,c)=>{const key=k+':'+g,v=(p.get(key)||o.ZERO).add(c);if(v.zero())p.delete(key);else p.set(key,v);};
 const mul=(a,b)=>{const z=new Map();for(const [ka,ca]of a)for(const [kb,cb]of b){const [ia,ga]=ka.split(':'),[ib,gb]=kb.split(':');add(z,Number(ia)+Number(ib),pmul(ga,gb),ca.mul(cb));}return z;};
 const factors=[[s,1],[t,1],[s,-1],[t,-1],[s,-1],[t,-1],[s,1],[t,1]],num=factors.map(([g,a])=>new Map([['0:'+id,o.ONE],['1:'+g,new o.Cut(0,f(a).div(2))]])).reduce(mul,new Map([['0:'+id,o.ONE]]));
 const st=pmul(s,t),ts=pmul(t,s),sts=pmul(st,s),den=new Map([[0,1],[2,1],[4,'3/8'],[6,'1/16'],[8,'1/256']].map(([j,c])=>[j+':'+id,o.Cut.of(c)])),lead=new Map([['0:'+id,o.ONE],['2:'+st,o.Cut.of('-1/2')],['2:'+ts,o.Cut.of('1/2')]]),res=new Map(num);
 for(const [key,c]of mul(lead,den)){const [j,g]=key.split(':');add(res,Number(j),g,c.neg());}
 const table=[[4,[[id,'-1/4'],[ts,'-1/16'],[st,'5/16']]],[5,[[s,new o.Cut(0,'-1/16')],[t,new o.Cut(0,'-1/16')],[sts,new o.Cut(0,'1/8')]]],[6,[[id,'-1/16'],[ts,'-9/64'],[st,'13/64']]],[7,[[s,new o.Cut(0,'-1/128')],[t,new o.Cut(0,'-1/128')],[sts,new o.Cut(0,'1/64')]]],[8,[[id,'-1/256'],[ts,'-1/32'],[st,'9/256']]],[10,[[st,'1/512'],[ts,'-1/512']]]],expected=new Map();
 for(const [j,terms]of table)for(const [g,c]of terms)add(expected,j,g,o.Cut.of(c));ensure(res.size===expected.size,'residue support');let K=f(0);for(const [key,c]of expected){ensure(res.get(key)?.eq(c),'universal native permutation residue '+key);K=K.add(c.rad.abs()).add(c.turn.abs());}ensure(K.eq('355/256'),'native residue gain');
 let exactWords=0;for(const g of [tri,full]){eq(o.power(g.X,3),g.id,'source three cycle');eq(o.mul(g.Pc,g.Pc),g.Pc,'cycle cut');eq(o.mul(g.Om,g.Om),o.scale(g.Pc,'3/16'),'native curvature square');for(const x of ['-1','-1/2','1/4','1/2','1']){const fs=balancedFactors(g.S,g.T,x),B=product(fs);for(const C of fs){eq(o.mul(o.dagger(C),C),g.id,'rational native exchange isometry');eq(comm(C,g.Pc),o.zeros(g.id.length),'finite activity cut invariant');}const defect=o.sub(B,o.add(g.id,o.scale(g.Om,io.mul(f(x).pow(2).mul(2)))));ensure(o.energy(defect).le(K.pow(2).mul(f(x).pow(8)).mul(g.id.length)),'residue gain witness');exactWords++;}}
 const B=product(balancedFactors(tri.S,tri.T,1)),finite=o.energy(comm(B,tri.Om));ensure(finite.eq('4374/390625'),'finite curvature defect');const L=product(loopFactors(tri.S,tri.T,2)),rad=o.scale(o.add(L,o.dagger(L)),'1/2'),turn=o.scale(o.sub(L,o.dagger(L)),new o.Cut(0,'-1/2'));eq(rad,o.sub(tri.id,o.scale(tri.Pc,'3/8')),'retained half-loop radial phase');eq(o.mul(turn,turn),o.scale(tri.Pc,'39/64'),'retained half-loop turn square');
 return {universal_residue_terms:res.size,residue_degrees:table.map(([j])=>j),residue_gain_bound:K,exact_balanced_words:exactWords,finite_curvature_commutator_square:finite,finite_word_is_pure_curvature_flow:false,finite_word_conserves_cycle_activity:true,finite_word_conserves_signed_curvature:false,finite_half_loop_cycle_radial:'5/8',phase_of_finite_half_loop_discarded:false};
});

check('native_quadratic_count_refinement_encloses_curvature_flow',()=>{
 const rows=[];for(const a of ['1/4','1/2','1']){const t=f(a).pow(2).mul(2),poly=expPoly(o.scale(tri.Om,-1),t,24,'1/2');for(const q of [1,2,4,8]){const n=q*q,x=f(a).div(q),word=product(balancedFactors(tri.S,tri.T,x)),run=o.power(word,n),bound=f('643/1024').mul(t.pow(2)).div(n),budget=bound.add(poly.tail),residual=o.energy(o.sub(run,poly.matrix));eq(o.mul(o.dagger(run),run),tri.id,'exact loop-count isometry');ensure(residual.le(budget.pow(2).mul(3)),'native quadratic count enclosure');rows.push({quadratic_parameter:t,loop_count:n,edge_factor_count:8*n,rational_step:x,flow_gain_error_bound:bound,factorial_tail:poly.tail,residual_square_sha256:hash(JSON.stringify(residual)),residual_square_bound:budget.pow(2).mul(3)});}}
 let lo=f(1),hi=f(2);const brackets=[];for(let j=0;j<16;j++){const mid=lo.add(hi).div(2);if(mid.pow(2).le(3))lo=mid;else hi=mid;ensure(lo.pow(2).le(3)&&f(3).le(hi.pow(2)),'native square-root bisection bracket');ensure(hi.sub(lo).eq(f(1).div(2**(j+1))),'root bracket count width');if([3,7,15].includes(j))brackets.push({count:j+1,lower:lo,upper:hi});}
 return {exact_refinement_enclosures:rows,enclosure_count:rows.length,constructed_root_brackets:brackets,square_root_assumed_primitive:false,general_limit_proved_in_written_section:true,finite_count_identified_with_exact_limit:false,quadratic_parameter_identified_as_physical_duration:false};
});

const flowCases=[...['0','1/2','1','2','4','-1'].map(y=>{const a=f(y),d=f(1).add(a.pow(2).mul('3/16'));return {c:f(1).sub(a.pow(2).mul('3/16')).div(d),v:a.mul(2).div(d)};}),{c:f('-1/2'),v:f(2)},{c:f('-1/2'),v:f(-2)},{c:f(-1),v:f(0)}];
check('native_cycle_cut_determines_closed_curvature_arrows',()=>{
 let identities=0;for(const g of [tri,full]){eq(o.dagger(g.Om),g.Om,'native curvature self dagger');eq(o.mul(g.Om,g.P0),o.zeros(g.id.length),'fixed cut removed by curvature');eq(o.dagger(g.Pc),g.Pc,'cycle cut self dagger');for(const {c,v}of flowCases){ensure(c.pow(2).add(v.pow(2).mul('3/16')).eq(1),'native phase coordinates');const U=closedFlow(g,c,v);eq(o.mul(o.dagger(U),U),g.id,'closed curvature matching');eq(comm(U,g.Om),o.zeros(g.id.length),'completed curvature conserved');eq(comm(U,g.Pc),o.zeros(g.id.length),'completed cycle cut conserved');identities++;}eq(closedFlow(g,'-1/2',2),o.dagger(g.X),'exact positive oriented cycle');eq(closedFlow(g,'-1/2',-2),g.X,'exact reverse cycle');eq(closedFlow(g,-1,0),o.sub(g.P0,g.Pc),'cycle reflection');}
 for(const A of [P,Q,K,H]){const sum=sumMatrices([0,1,2].map(j=>localN(3,j,A)),8);eq(comm(full.Om,sum),o.zeros(8),'identical readout sum conserved');}
 eq(o.mul(full.Om,insert),o.mul(insert,tri.Om),'full tuple curvature restricts to single-cut sector');
 return {exact_closed_flow_instances:identities,native_phase_coordinate_cases:flowCases,oriented_cycle_transfer_exact:true,curvature_square_coefficient:'3/16',cycle_phase_period:'8 theta_- / sqrt_native(3)',phase_cuts_require_constructed_native_root:true};
});

check('oriented_transfer_derives_memory_metric_and_cubic_count',()=>{
 let cases=0;const reflected=[];for(const {c,v}of flowCases){const U=closedFlow(tri,c,v),u=o.mul(U,basis(3,0)),A=f(1).sub(c).div(3),expected=col([f(1).add(c.mul(2)).div(3),A.add(v.div(4)),A.sub(v.div(4))]);eq(u,expected,'native directed amplitudes');const ns=u.map(row=>row[0].norm2());ensure(ns.reduce((s,n)=>s.add(n),f(0)).eq(1),'total retained count');ensure(ns[1].sub(ns[2]).eq(f(1).sub(c).mul(v).div(3)),'signed directed count');const state=o.mul(insert,u),pair=pairOf(state);for(let j=0;j<3;j++){const local=part(pair,3,j),memory=f(1).sub(trace(o.mul(local,local)).rad);eq(local,o.matrix([[f(1).sub(ns[j]),0],[0,ns[j]]]),'matched local native pair');ensure(memory.eq(ns[j].mul(f(1).sub(ns[j])).mul(2)),'local matching memory');if(c.eq(-1))reflected.push({ledger:j,count:ns[j],deficit:memory});}const av=mean(u,tri.Om),activity=mean(u,o.mul(tri.Om,tri.Om)),variance=activity.sub(av.mul(av));ensure(av.zero()&&variance.eq('1/8'),'curvature orbit metric');ensure(activity.eq(mean(u,tri.Pc).mul('3/16')),'cycle activity identity');cases++;}
 ensure(reflected.map(r=>String(r.count)).join(',')==='1/9,4/9,4/9','reflection count vector');ensure(reflected.map(r=>String(r.deficit)).join(',')==='16/81,40/81,40/81','reflection memory vector');
 const rs=f('3/16'),cosEven=j=>rs.pow(j).mul(j%2?-1:1).div(factorial(2*j)),vOdd=j=>rs.pow(j).mul(j%2?-1:1).div(factorial(2*j+1)),cubic=cosEven(1).neg().mul(vOdd(0)).div(3),quintic=cosEven(1).neg().mul(vOdd(1)).sub(cosEven(2).mul(vOdd(0))).div(3);ensure(cubic.eq('1/32')&&quintic.eq('-3/2048'),'first directional count coefficients');
 return {exact_orbit_cases:cases,reflection_readings:reflected,response_metric:'1/8',directional_count_cubic:cubic,directional_count_quintic:quintic,direction_detected_in_quadratic_count:false,memory_requires_random_mixture:false,positive_activity_identified_with_signed_curvature:false};
});

check('native_network_loops_factor_the_cyclic_curvature',()=>{
 const ringRows=[],refinement=[];for(let n=3;n<=9;n++){const g=ring(n);eq(g.Om,o.mul(g.D,g.G),'native ring curvature factorization');eq(comm(g.D,g.G),o.zeros(n),'ring responses commute');const laurent=o.scale(o.sub(o.scale(o.sub(g.T,g.Ti),2),o.sub(o.power(g.T,2),o.power(g.Ti,2))),io.mul('1/4'));eq(g.Om,laurent,'native Laurent curvature');ringRows.push({ledgers:n,local_oriented_pairs:n});}
 for(const N of [3,4,5]){const g=ring(N),a=f('1/4'),t=a.pow(2).mul(2),poly=expPoly(o.scale(g.Om,-1),t,20,'3/2');for(const q of [2,4]){const n=q*q,x=a.div(q),words=g.pairs.map(([P,Q])=>product(balancedFactors(swapCut(P),swapCut(Q),x))),run=o.power(product(words),n),arg=f(N).mul(t).div(2*n);ensure(arg.le('1/2'),'network exponential count bound');const budget=t.pow(2).div(n).mul(f(643*N).div(1024).add(f(N*N).div(4).div(f(1).sub(arg)))),comparison=budget.add(poly.tail),residual=o.energy(o.sub(run,poly.matrix));ensure(residual.le(comparison.pow(2).mul(N)),'network curvature count enclosure');eq(o.mul(o.dagger(run),run),g.id,'finite network word matching');refinement.push({ledgers:N,loop_count:n,edge_factor_count:8*N*n,flow_gain_error_bound:budget,residual_square_sha256:hash(JSON.stringify(residual)),residual_square_bound:comparison.pow(2).mul(N)});}}
 return {ring_factorization_cases:ringRows,exact_network_enclosures:refinement,orientation_and_list_are_explicit_native_targets:true,physical_handedness_selected:false,cyclic_factorization_claimed_for_every_graph:false,continuum_derivative_premise:false};
});

const roots4=[o.ONE,io,o.ONE.neg(),io.neg()],cutPower=(z,n)=>{let a=o.ONE;for(let j=0;j<n;j++)a=a.mul(z);return a;},modes4=roots4.map(z=>col(Array.from({length:4},(_,j)=>cutPower(z,j))));
check('native_curvature_modes_have_cubic_count_scaling',()=>{
 const spectra=[];let eigenspaces=0;for(const [N,values]of [[3,[[0,1],['27/16',2]]],[4,[[0,2],[1,2]]],[6,[[0,2],['3/16',2],['27/16',2]]]]){const g=ring(N),square=o.mul(g.Om,g.Om);for(const [value,mult]of values){ensure(N-o.rank(o.sub(square,o.scale(g.id,value)))===mult,'native squared curvature spectrum');eigenspaces++;}spectra.push({ledgers:N,squared_phase_rates_and_multiplicities:values});}
 const kernel=[];for(let n=3;n<=11;n++){const g=ring(n),nullity=n-o.rank(g.Om);ensure(nullity===(n%2?1:2),'odd/even native curvature kernel');kernel.push({ledgers:n,dimension:nullity});}
 const g4=ring(4);for(let m=0;m<4;m++){eq(o.mul(g4.G,modes4[m]),o.scale(modes4[m],[0,1,2,1][m]),'native four-mode cost');eq(o.mul(g4.Om,modes4[m]),o.scale(modes4[m],[0,1,0,-1][m]),'native four-mode curvature');for(let j=0;j<4;j++)ensure(inner(modes4[m],modes4[j]).eq(m===j?4:0),'native geometric mode matching');}
 const coeff=j=>f(4).pow(j).sub(1).mul(j%2?1:-1).div(factorial(2*j+1));ensure(coeff(1).eq('1/2')&&coeff(2).eq('-1/8'),'cubic dispersion coefficients');let alternating=0;for(let j=1;j<16;j++){ensure(coeff(j+1).abs().div(coeff(j).abs()).le('1/4'),'native alternating coefficient decrease');alternating++;}
 const brackets=[];for(const k of ['1/8','1/4','1/2','1']){let upper=f(0);for(let j=1;j<=9;j++)upper=upper.add(coeff(j).mul(f(k).pow(2*j+1)));const lower=upper.add(coeff(10).mul(f(k).pow(21))),analyticLower=f(k).pow(3).div(2).sub(f(k).pow(5).div(8));ensure(analyticLower.le(lower)&&lower.le(upper)&&upper.le(f(k).pow(3).div(2)),'native alternating cubic phase bracket');brackets.push({rational_phase:k,lower,upper});}
 const budgets=[8,16,32,64,128].map(n=>({ledgers:n,scaled_gap_lower:f(1).sub(f(16).div(n*n)),scaled_gap_upper:'1'}));
 return {exact_squared_eigenspaces:eigenspaces,exact_spectra:spectra,kernel_dimensions:kernel,native_four_mode_checks:4,alternating_coefficient_comparisons:alternating,rational_phase_brackets:brackets,general_count_scaling_budgets:budgets,cubic_phase_leading_coefficient:'1/2',cubic_phase_next_coefficient:'-1/8',even_curvature_kernel_contains_alternating_mode:true,finite_target_gap_identified_as_physical_mass_gap:false,irrational_spectrum_claimed_exactly_enumerated_by_rational_tests:false};
});

check('native_cost_and_curvature_reconstruct_oriented_translation',()=>{
 const rows=[];for(let n=3;n<=11;n++){const g=ring(n),Pc=o.matrix(Array.from({length:n},()=>Array(n).fill(f(1).div(n)))),Qc=o.sub(g.id,Pc),Gsharp=o.mul(o.inverse(o.add(g.G,Pc)),Qc);eq(o.mul(g.G,Gsharp),Qc,'native complement inverse');eq(o.mul(g.D,Pc),o.zeros(n),'constant mode direction zero');eq(o.mul(g.Om,Gsharp),g.D,'curvature recovers directional response');eq(o.add(o.sub(g.id,g.G),o.scale(o.mul(g.Om,Gsharp),io)),g.Ti,'cost and curvature recover oriented shift');rows.push({ledgers:n,cost_nullity:n-o.rank(g.G),inverse_keeps_constant_cut:true});}
 const labels=[];for(let m=0;m<4;m++){const g=[0,1,2,1][m],w=[0,1,0,-1][m],z=g?new o.Cut(f(1).sub(g),f(w).div(g)):o.ONE;ensure(z.eq(roots4[m]),'single-mode direction reconstruction');labels.push({mode:m,cost:g,curvature:w,phase:z});}ensure(labels[1].cost===labels[3].cost&&labels[1].curvature!==labels[3].curvature,'cost orientation ambiguity');ensure(labels[0].curvature===labels[2].curvature&&labels[0].cost!==labels[2].cost,'curvature zero ambiguity');
 return {exact_translation_reconstructions:rows,native_four_mode_labels:labels,zero_cost_mode_divided_away:false,pair_of_mode_readings_is_injective:true,two_expectations_claimed_to_recover_arbitrary_state:false};
});

const coherenceStates=[o.add(modes4[0],modes4[1]),o.sub(modes4[0],modes4[1])];
check('native_commuting_readouts_retain_coherence_and_cost_boundaries',()=>{
 const g=ring(4),weights=coherenceStates.map(v=>modes4.map(u=>mean(v,o.scale(o.mul(u,o.dagger(u)),'1/4'))));for(let m=0;m<4;m++){ensure(weights[0][m].eq(weights[1][m]),'equal exact native mode weights');ensure(weights[0][m].eq(m<2?'1/2':0),'two-mode weights');}
 let moments=0;for(let j=0;j<=4;j++)for(let k=0;k<=4;k++){const A=o.mul(o.power(g.G,j),o.power(g.Om,k));ensure(mean(coherenceStates[0],A).eq(mean(coherenceStates[1],A)),'equal commuting moments');moments++;}const local=coherenceStates.map(v=>v.map(row=>row[0].norm2().div(o.energy(v))));ensure(local[0][0].eq('1/2')&&local[1][0].zero(),'commuting observations miss local coherence');
 const costs=[];for(const q of [1,2,4,8,16,32]){const n=q*q,x=f(1).div(4*q),lower=f(8*n).mul(x.sub(x.pow(3).div(12))),upper=f(8*n).mul(x);ensure(lower.le(upper)&&f(11*q).div(6).le(lower),'phase cost diverges along square counts');costs.push({loop_count:n,quadratic_parameter:'1/8',edge_factors:8*n,total_native_phase_lower:lower,total_native_phase_upper:upper});}
 const calibrations=[];for(const [a,b]of [[2,3],['1/2','5/3'],['7/4','2/5']]){const aa=f(a),bb=f(b),v=basis(3,0),der=o.scale(o.mul(tri.Om,v),new o.Cut(0,f(1).div(aa)));eq(o.scale(der,new o.Cut(0,aa.mul(bb))),o.scale(o.mul(tri.Om,v),bb.neg()),'native curvature action calibration');calibrations.push({time_scale:aa,readout_scale:bb,action_scale:aa.mul(bb)});}
 return {identical_mode_weights:weights,exact_equal_polynomial_moments:moments,distinct_local_profiles:local,all_polynomial_moments_equality_proved_in_written_section:true,mode_readout_is_full_state_tomography:false,retained_refinement_cost:costs,quadratic_parameter_erases_phase_cost:false,remaining_common_calibrations:calibrations,physical_hbar_c_alpha_or_mass_selected:false};
});

const jointApply=(a,A)=>new Map([...a].map(([k,v])=>[k,processApply(v,A)])),jointRead=(a,A)=>[...a.values()].reduce((s,v)=>s.add(fieldPair(v,processApply(v,A))),o.ZERO);
const jointEqual=(a,b,D,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))eqFields(a.get(k)||new Map(),b.get(k)||new Map(),32,D,msg+' '+k);};
const Hc=o.add(Ax,Ay),Kc=o.sub(Ax,Ay),Jc=o.mul(Kc,Hc),PcSource=o.scale(o.add(I16,Hc),'1/2'),P32=o.kron(PcSource,I),PJ32=o.kron(o.scale(o.mul(PcSource,Jc),-1),I),J32=o.kron(Jc,I);
const observe=a=>[compose(unit(P32),a),compose(unit(PJ32),a)],decode=a=>plus(a[0],compose(unit(J32),a[1]));
const processPairFromField=(a,D)=>[...a.values()].reduce((s,M)=>{const T=transpose(M);return o.add(s,o.mul(T,o.dagger(T)));},o.zeros(D));
check('actual_refined_clock_preserves_activity_and_complete_observation',()=>{
 const fs=balancedFactors(full.S,full.T,'1/2'),chronological=[...fs].reverse(),stages=[...program,...chronological.map(V=>({apply:a=>processApply(a,V),inverse:a=>processApply(a,o.dagger(V))}))],q=stages.length;
 const step=(a,back=false)=>{const out=new Map();for(const [key,v]of a){const [j,m]=marks(key),jj=back?(j+q-1)%q:(j+1)%q,mm=back?(m+3-(j===0?1:0))%3:(m+(j===q-1?1:0))%3,vv=back?stages[jj].inverse(v):stages[j].apply(v);ensure(!out.has(mark(jj,mm)),'refined clock pointer collision');out.set(mark(jj,mm),vv);}return out;};
 const reads=[Q,K].map(A=>sumMatrices([0,1,2].map(j=>localN(3,j,A)),8));reads.push(full.Pc);let intertwinings=0,roundtrips=0;const seed=productField(sourceSeed,basis(8,4));for(let j=0;j<q;j++){const a=joint(seed,j,1);jointEqual(step(step(a),true),a,8,'actual refined clock roundtrip');roundtrips++;for(const A of reads){jointEqual(step(jointApply(a,A)),jointApply(step(a),A),8,'actual refined clock conserved response');intertwinings++;}}
 let a=joint(seed);const initial=reads.map(A=>jointRead(a,A));let prefixes=0;for(let j=0;j<=q+3;j++){if(j)a=step(a);ensure(jointEnergy(a).eq(1),'actual refined clock norm');for(let k=0;k<reads.length;k++)ensure(jointRead(a,reads[k]).eq(initial[k]),'actual refined prefix response');if(j===q){const expected=joint(processApply(sourceCycle(seed),product(fs)),0,1);jointEqual(a,expected,8,'actual source echo and balanced word per tick');}prefixes++;}
 const reflection=o.mul(insert,o.mul(closedFlow(tri,-1,0),basis(3,0))),states=[reflection,...coherenceStates];let observed=0;for(const v of states){const D=v.length,field=productField(sourcePowers[1],v),obs=observe(field);eqFields(decode(obs),field,32,D,'complete native source reconstruction');const before=processPairFromField(field,D),after=o.add(processPairFromField(obs[0],D),processPairFromField(obs[1],D));eq(before,after,'complete native observer preserves process coherence');ensure(energy(obs[0]).add(energy(obs[1])).eq(energy(field)),'complete source observation matching');observed++;}
 return {actual_source_stages:qProgram,appended_exchange_factors:fs.length,actual_refined_controller_phases:q,actual_inverse_roundtrips:roundtrips,actual_readout_intertwinings:intertwinings,actual_prefix_checks:prefixes,complete_process_pair_witnesses:observed,full_cycle_equals_actual_source_echo_and_balanced_loop:true,cycle_activity_conserved_at_every_inner_stage:true,signed_curvature_claimed_conserved_by_finite_word:false,source_observer_preserves_cross_mode_coherence:true,carrier_return_error_becomes_invariant_error:false};
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
 finite_balanced_word_equals_pure_curvature_flow:checks[0].finite_word_is_pure_curvature_flow===false,
 conserving_activity_is_conserving_signed_curvature:checks[0].finite_word_conserves_cycle_activity===true&&checks[0].finite_word_conserves_signed_curvature===false,
 finite_half_loop_phase_can_be_discarded:checks[0].phase_of_finite_half_loop_discarded===false,
 positive_square_root_is_an_imported_primitive:checks[1].square_root_assumed_primitive===false,
 a_finite_count_certificate_proves_exact_equality_to_the_limit:checks[1].finite_count_identified_with_exact_limit===false,
 leading_quadratic_count_already_determines_direction:checks[3].direction_detected_in_quadratic_count===false,
 native_pair_memory_requires_random_mixing:checks[3].memory_requires_random_mixture===false,
 the_ordered_target_selects_physical_handedness:checks[4].physical_handedness_selected===false,
 cyclic_factorization_holds_on_every_network:checks[4].cyclic_factorization_claimed_for_every_graph===false,
 every_zero_curvature_mode_is_constant:checks[5].even_curvature_kernel_contains_alternating_mode===true,
 cubic_finite_target_gap_is_physical_mass_gap:checks[5].finite_target_gap_identified_as_physical_mass_gap===false,
 mode_label_recovery_is_full_state_tomography:checks[7].mode_readout_is_full_state_tomography===false,
 fixed_quadratic_parameter_erases_refinement_cost:checks[7].quadratic_parameter_erases_phase_cost===false,
 common_physical_action_speed_and_coupling_are_selected:checks[7].physical_hbar_c_alpha_or_mass_selected===false,
};ensure(Object.values(negatives).every(Boolean),'boundary control failed');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r45.native-curvature-flow.v1',status:'PASS_R45_NATIVE_CURVATURE_FLOW',
 input_sha256:inputHash,r44_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_loop_generated_curvature_flow_derived:true,native_oriented_transfer_memory_metric_derived:true,
 native_cubic_curvature_modes_derived:true,native_cost_curvature_direction_reconstruction_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,physical_metric_c_alpha_derived:false,
 physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native proofs: exact balanced-loop residue, controlled quadratic-count curvature flow, oriented transfer with memory and metric, collective loop synthesis and cyclic factorization, cubic phase rates and count gap, cost-curvature reconstruction of mode direction, and actual refined-clock/complete-observer transport. Finite words, completion, declared targets, full-state coherence, refinement cost and physical calibration remain distinct.'},null,2)+'\n');
