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
ensure(inputHash===arg('--expected-input-sha256'),'R30 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r29-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r29-sha256'),'R30 R29 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R29_NATIVE_PAIRED_FIELD_DYNAMICS','R30 parent status');
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
const lap=compose(dag(dx),dx);
const paired=plus(plus(lift(P,minus(unit(I),scale(lap,'1/2'))),lift(E01,scale(dag(dx),'-1/2'))),plus(lift(E10,dx),lift(Q,unit(I))));
const pairGram=plus(plus(lift(P,unit(I)),lift(Q,unit(o.scale(I,'1/2')))),plus(lift(E01,scale(dag(dx),'1/4')),lift(E10,scale(dx,'1/4'))));
// Input pair is (b, phi); its observable pair is (e=N_x b+Omega phi, b).
const observable=plus(plus(lift(P,nx),lift(E01,omega)),lift(E10,unit(I)));
const coupled=plus(lift(I,vx),lift(E01,compose(dx,omega)));
const coupledInverse=minus(lift(I,vxD),lift(E01,compose(dx,omega)));
const jointGram=compose(compose(dag(observable),pairGram),observable);
const efield=(b,phi)=>plus(compose(nx,b),compose(omega,phi));
const fieldEnergy=(e,b)=>energy(e).add(energy(b).div(2)).add(dot(b,compose(dx,e)).div(2));
const coupledStep=(b,phi)=>({b:plus(compose(vx,b),compose(dx,compose(omega,phi))),phi:compose(vx,phi)});
const probes=[unit(e0),unit(e1),new Map([['-1,0',col([2,-1])],['1,1',col([-1,3])],['0,-2',col([1,1])]])];
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

check('two_address_source_histories_reproduce_shared_role_transport',()=>{
 let prefixes=0,comparisons=0;
 for(const role of [e0,e1])for(const swapped of [false,true]){let histories=[{x:0,y:0,v:role}],direct=unit(role);
  for(let n=1;n<=8;n++){
   const direction=(Math.floor((n-1)/2)%2)^(swapped?1:0),next=[];
   for(const old of histories)for(const letter of [H,K]){const v=o.mul(letter,old.v),sign=o.isZero(o.mul(P,v))?-1:1;
    next.push({x:old.x+(direction?0:sign),y:old.y+(direction?sign:0),v});}
   histories=next;prefixes+=histories.length;
   if(n%4===0){const collected=new Map();for(const h of histories){ensure(h.x%2===0&&h.y%2===0,'coarse native counts');addAt(collected,key(h.x/2,h.y/2),o.scale(h.v,new o.F(1,1n<<BigInt(n/2))));}
    direct=compose(swapped?reverse:forward,direct);eqFields(direct,collected,2,1,'literal native direction-block census');comparisons++;
   }
  }
 }
 const independentX=new Map([...vx].map(([k,M])=>[k,o.kron(M,I)])),independentY=new Map([...vy].map(([k,M])=>[k,o.kron(I,M)]));
 eqFields(compose(independentX,independentY),compose(independentY,independentX),4,4,'independent role copies commute');
 ensure(!energy(compose(omega,unit(e0))).eq(0),'shared role target is order-sensitive');
 return {literal_source_prefixes:prefixes,independent_endpoint_comparisons:comparisons,max_literal_event:8,
  address_labels_come_from_event_block_counter:true,shared_source_role_is_load_bearing:true};
});

check('directional_order_residue_and_closed_loop_have_exact_native_factorization',()=>{
 eqFields(omega,scale(compose(g,unit(A0)),'1/4'),2,2,'exact directional curvature factor');
 eqFields(compose(dag(omega),omega),scale(compose(dag(g),g),'1/8'),2,2,'central curvature square');
 eqFields(loopDefect,compose(compose(omega,vxD),vyD),2,2,'actual closed protocol loop');
 eqFields(compose(dag(loopDefect),loopDefect),compose(dag(omega),omega),2,2,'same-input order/loop quadratic response');
 eqFields(compose(dag(loop),loop),unit(I),2,2,'closed protocol pairing');
 for(const v of [vx,vy])eqFields(compose(dag(v),v),unit(I),2,2,'directional source pairing');
 ensure(g.size===6,'six-term native difference residue');
 const norm=energy(compose(omega,unit(e0)));ensure(norm.eq('3/4'),'localized unit preparation curvature energy');
 return {complete_Laurent_identities:7,scalar_difference_terms:g.size,localized_unit_response_energy:norm.toString(),
  zero_response_for_all_compact_preparations:false,closed_protocol_is_not_assumed_a_spacetime_plaquette:true};
});

check('curvature_intertwines_inverse_directional_transport',()=>{
 let operators=0,fields=0,norms=0;
 for(const v of [vx,vy]){
  eqFields(scale(compose(compose(unit(A0),v),unit(A0)),'1/2'),dag(v),2,2,'native role reflection reverses directional transport');
  eqFields(compose(omega,v),compose(dag(v),omega),2,2,'curvature reverse-continuation intertwiner');operators+=2;
 }
 eqFields(compose(omega,forward),compose(dual,omega),2,2,'ordered two-direction dual continuation');
 eqFields(loopDefect,compose(reverse,omega),2,2,'closed-loop response has a native continuation decoder');operators+=2;
 for(const seed of probes){let phi=seed,kappa=compose(omega,seed);const initial=energy(kappa);
  for(let n=0;n<=5;n++){
   eqFields(compose(omega,phi),kappa,2,1,'curvature field follows the dual transport');fields++;
   ensure(energy(kappa).eq(initial),'native curvature response energy invariant');
   ensure(energy(compose(loopDefect,phi)).eq(initial),'actual loop response energy invariant');norms+=2;
   if(n<5){phi=compose(forward,phi);kappa=compose(dual,kappa);}
  }
 }
 return {complete_Laurent_intertwiners:operators,dual_field_evolution_checks:fields,curvature_and_loop_norm_checks:norms,
  max_four_event_block:5,dual_is_inverse_of_reversed_block_order:true};
});

check('intersecting_native_interfaces_localize_curvature',()=>{
 let cornerChecks=0,blindChecks=0,loopChecks=0;
 for(const v of [e0,e1,col([2,-1])]){
  const corner=(x,y)=>x<0&&y<0?v:o.zeros(2,1),response=new Map([['0,-1',o.scale(o.mul(A0,v),'1/4')],['-1,0',o.scale(o.mul(A0,v),'-1/4')]]);
  for(let x=-3;x<=3;x++)for(let y=-3;y<=3;y++){
   eq(at(omega,corner,x,y),value(response,x,y),'localized corner formula');cornerChecks++;
   eq(at(omega,(a,b)=>a>=0&&b<0?v:o.zeros(2,1),x,y),o.scale(value(response,x,y),-1),'corner orientation reversal');cornerChecks++;
   for(const blind of [(a,b)=>a<0?v:o.zeros(2,1),(a,b)=>b<0?v:o.zeros(2,1),(a,b)=>a+b<0?v:o.zeros(2,1),()=>v]){
    eq(at(omega,blind,x,y),o.zeros(2,1),'straight/diagonal/uniform loop blindness');blindChecks++;
   }
  }
  ensure(energy(response).eq(o.energy(v).div(4)),'corner response energy');
  const loopResponse=compose(reverse,response);
  for(let x=-4;x<=4;x++)for(let y=-4;y<=4;y++){eq(at(loopDefect,corner,x,y),value(loopResponse,x,y),'actual corner loop response');loopChecks++;}
 }
 return {oriented_corner_readings:cornerChecks,straight_diagonal_and_uniform_readings:blindChecks,actual_loop_readings:loopChecks,
  corner_response_energy_ratio:'1/4',infinite_background_norm_used:false};
});

check('loop_readout_closes_a_reversible_paired_field_interaction',()=>{
 eqFields(compose(observable,coupled),compose(paired,observable),4,4,'loop-sourced paired-field closure');
 eqFields(compose(coupled,coupledInverse),unit(o.identity(4)),4,4,'coupled inverse');
 eqFields(compose(coupledInverse,coupled),unit(o.identity(4)),4,4,'coupled inverse other side');
 eqFields(compose(compose(dag(coupled),jointGram),coupled),jointGram,4,4,'derived joint energy conservation');
 let matches=0,norms=0;const cases=probes.map((phi,i)=>({phi,b:i?probes[i-1]:new Map()}));
 for(const initial of cases){let {b,phi}=initial;const start=fieldEnergy(efield(b,phi),b),kenergy=energy(compose(omega,phi));
  for(let n=0;n<10;n++){
   const e=efield(b,phi),kappa=compose(omega,phi),next=coupledStep(b,phi),en=efield(next.b,next.phi);
   eqFields(next.b,plus(b,compose(dx,e)),2,1,'first paired induction equation');
   eqFields(en,minus(e,scale(compose(dag(dx),next.b),'1/2')),2,1,'second paired induction equation');
   eqFields(minus(en,compose(nx,next.b)),compose(vxD,kappa),2,1,'curvature is the continued compatibility residue');matches+=3;
   ensure(fieldEnergy(en,next.b).eq(start),'joint native field form is conserved');
   ensure(energy(compose(omega,next.phi)).eq(kenergy),'loop driver response norm is conserved');norms+=2;
   ({b,phi}=next);
  }
 }
 const first=coupledStep(new Map(),unit(e0));ensure(!energy(first.b).eq(0),'loop creates a nonzero coupled-field response');
 return {complete_Laurent_closure_inverse_and_Gram_identities:4,field_closure_matches:matches,conserved_form_checks:norms,
  first_induced_b_energy:energy(first.b).toString(),interaction_readout_c_is_Omega_phi:true,physical_interaction_selected:false};
});

const ring=(op,L,r=2,c=2)=>{const M=o.zeros(r*L*L,c*L*L);for(const [k,A]of op){const [dx,dy]=xy(k);
 for(let x=0;x<L;x++)for(let y=0;y<L;y++){const u=((x-dx)%L+L)%L,v=((y-dy)%L+L)%L,a=x*L+y,b=u*L+v;
  for(let i=0;i<r;i++)for(let j=0;j<c;j++)M[a*r+i][b*c+j]=M[a*r+i][b*c+j].add(A[i][j]);
 }}return M;};
const turn=(a,b)=>{const den=f(a).pow(2).add(f(b).pow(2));return o.add(o.scale(I,f(a).pow(2).sub(f(b).pow(2)).div(den)),o.scale(R,f(a).mul(b).mul(2).div(den)));};
check('local_native_frames_preserve_loop_response_and_do_not_create_curvature',()=>{
 const L=3,n=2*L*L,id=o.identity(n),G=o.zeros(n),options=[[1,0],[1,1],[2,1],[3,2]];
 for(let z=0;z<L*L;z++){let g=turn(...options[z%4]);if(z%3===0)g=o.mul(g,H);eq(o.mul(o.dagger(g),g),I,'native count frame norm');
  for(let i=0;i<2;i++)for(let j=0;j<2;j++)G[2*z+i][2*z+j]=g[i][j];}
 const GD=o.dagger(G),change=A=>o.mul(o.mul(G,A),GD),Tx=ring(tx,L),Ty=ring(ty,L),X=ring(vx,L),Y=ring(vy,L),Omega=ring(omega,L);
 eq(o.mul(GD,G),id,'local native frame isometry');
 eq(o.commutator(change(Tx),change(Ty)),o.zeros(n),'pure frame links remain flat');
 eq(o.commutator(change(X),change(Y)),change(Omega),'curvature transforms by the same native frame');
 const lp=o.sub(o.mul(o.mul(o.mul(change(X),change(Y)),o.dagger(change(X))),o.dagger(change(Y))),id);
 eq(lp,change(ring(loopDefect,L)),'actual loop response is covariant');
 const Pg=change(o.kron(o.identity(L*L),P)),Qg=o.sub(id,Pg),Hg=change(o.kron(o.identity(L*L),H)),Kg=change(o.kron(o.identity(L*L),K));
 const Ng=o.scale(o.mul(o.add(Pg,o.mul(o.dagger(change(Tx)),Qg)),o.add(Hg,Kg)),'1/2');
 eq(o.add(id,o.mul(Ng,o.sub(change(Tx),id))),change(X),'transform role cuts and shifts together');
 const q=col(Array.from({length:n},(_,i)=>i%2?-1:2)),qg=o.mul(G,q);
 eq(o.mul(o.sub(change(Tx),id),qg),o.zeros(n,1),'covariant background difference is zero');
 ensure(!o.isZero(o.mul(o.sub(Tx,id),qg)),'untransformed difference produces a false apparent field');
 const v=col(Array.from({length:n},(_,i)=>(i%7)-3));
 ensure(o.energy(o.mul(change(Omega),o.mul(G,v))).eq(o.energy(o.mul(Omega,v))),'native loop energy is frame invariant');
 return {cyclic_addresses:L*L,local_native_count_frames:L*L,exact_covariance_and_flatness_identities:8,
  position_dependent_frame_is_not_a_physical_gauge_field:true,ordinary_complex_phase_input:false};
});

const average=(L,a,b,weighted=false)=>{const op=new Map();for(let k=0;k<L;k++)addAt(op,key(k*a,k*b),o.scale(I,new o.F(weighted?k:1,L)));return op;};
const cyclic=[];
check('finite_cyclic_curvature_blindness_is_exactly_three_average_sectors',()=>{
 const rows=[];let identities=0;
 for(let L=1;L<=5;L++){
  const size=2*L*L,id=o.identity(size),ax=average(L,1,0),ay=average(L,0,1),ad=average(L,1,-1),a0=compose(ax,ay);
  const pi=minus(plus(plus(ax,ay),ad),scale(a0,2)),Pi=ring(pi,L),Omega=ring(omega,L),Ax=ring(ax,L),Ay=ring(ay,L),Ad=ring(ad,L),Amean=ring(a0,L);
  for(const A of [Ax,Ay,Ad]){eq(o.mul(A,A),A,'native average is a cut');eq(o.dagger(A),A,'native average is self-dagger');identities+=2;}
  for(const [A,B]of [[Ax,Ay],[Ax,Ad],[Ay,Ad]]){eq(o.mul(A,B),Amean,'distinct average cuts meet in uniform sector');identities++;}
  eq(o.mul(Pi,Pi),Pi,'complete blind-sector projector');eq(o.mul(Omega,Pi),o.zeros(size),'all three sectors are curvature blind');identities+=2;
  const bx=average(L,1,0,true),by=average(L,0,1,true),bd=average(L,1,-1,true),dd=minus(compose(tx,dag(ty)),unit(I));
  for(const [d,b,a]of [[dx,bx,ax],[dy,by,ay],[dd,bd,ad]]){eq(ring(compose(d,b),L),o.sub(id,ring(a,L)),'native weighted-count difference inverse');identities++;}
  const reconstruct=scale(compose(compose(compose(compose(tx,bx),by),bd),unit(A0)),2);
  eq(o.mul(ring(reconstruct,L),Omega),o.sub(id,Pi),'curvature reconstructs exactly the retained sector');identities++;
  const blind=2*(3*L-2),rank=2*(L-1)*(L-2);
  ensure(o.rank(Pi)===blind&&o.rank(Omega)===rank,'complete native curvature ranks');
  rows.push({period:L,source_roles:size,blind_roles:blind,curvature_rank:rank});cyclic.push({L,size,pi,Pi,Omega,rank,blind});
 }
 ensure(o.isZero(cyclic[1].Omega)&&!o.isZero(cyclic[2].Omega),'two-period aliasing versus nonzero three-period response');
 return {count_quotients:rows,complete_average_inverse_and_reconstruction_identities:identities,
  Fourier_or_classical_spectral_theorem_used:false,period_two_is_completely_curvature_blind:true};
});

check('curvature_reconstruction_and_coupled_energy_quotient_are_exact',()=>{
 const rows=[];let invariants=0,gramRanks=0;
 for(const {L,size,pi,Pi,blind}of cyclic){
  const X=ring(vx,L),Y=ring(vy,L),id=o.identity(size);eq(o.mul(Pi,X),o.mul(X,Pi),'blind sector x-invariant');eq(o.mul(Pi,Y),o.mul(Y,Pi),'blind sector y-invariant');invariants+=2;
  const J=ring(observable,L,4,4),C=ring(coupled,L,4,4),lost=ring(lift(Q,pi),L,4,4),G=ring(jointGram,L,4,4),n=2*size;
  eq(o.mul(J,lost),o.zeros(n),'only curvature-blind phi is forgotten');eq(o.mul(C,lost),o.mul(lost,C),'blind interaction subspace is invariant');invariants+=2;
  const wanted=2*size-blind;ensure(o.rank(J)===wanted,'exact interaction readout rank');
  eq(o.mul(G,lost),o.zeros(n),'derived joint energy vanishes on the exact blind kernel');invariants++;
  if(L<=4){ensure(o.rank(G)===wanted,'no extra null direction in the joint field form');gramRanks++;}
  rows.push({period:L,joint_input_roles:n,forgotten_roles:blind,closed_interaction_quotient_roles:wanted});
 }
 return {quotients:rows,exact_invariance_and_kernel_identities:invariants,full_joint_Gram_rank_checks:gramRanks,
  joint_form_positive_on_its_observable_quotient:true,blind_source_information_erased_or_declared_nonexistent:false};
});

const rowTotals=b=>{const sums=new Map();for(const [k,v]of b){const [,y]=xy(k);addAt(sums,''+y,v);}return sums;};
const eqTotals=(a,b)=>{for(const k of new Set([...a.keys(),...b.keys()]))eq(a.get(k)||o.zeros(2,1),b.get(k)||o.zeros(2,1),'boundary residue on native row '+k);};
const current=(b,x,y)=>{const a=value(b,x,y),d=value(b,x+1,y),h=a[0][0].rad.add(a[1][0].rad),k=d[0][0].rad.sub(d[1][0].rad);return h.pow(2).sub(k.pow(2)).div(4).add(h.mul(k).div(2));};
check('loop_driven_local_field_balance_retains_boundary_residue',()=>{
 let balances=0,totals=0;
 for(const phi0 of probes){let phi=phi0,b=unit(e0);const initial=rowTotals(b);
  for(let n=0;n<6;n++){
   const kappa=compose(omega,phi),source=compose(dx,kappa),vb=compose(vx,b),next=coupledStep(b,phi),keys=[...b.keys(),...next.b.keys(),...source.keys()].map(xy);
   const xmin=Math.min(...keys.map(x=>x[0]))-1,xmax=Math.max(...keys.map(x=>x[0]))+1,ymin=Math.min(...keys.map(x=>x[1]))-1,ymax=Math.max(...keys.map(x=>x[1]))+1;
   for(let x=xmin;x<=xmax;x++)for(let y=ymin;y<=ymax;y++){
    const s=value(source,x,y),v=value(vb,x,y),injection=o.mul(o.dagger(v),s)[0][0].rad.mul(2).add(o.energy(s));
    ensure(o.energy(value(next.b,x,y)).sub(o.energy(value(b,x,y))).eq(current(b,x-1,y).sub(current(b,x,y)).add(injection)),'loop-driven local native balance');balances++;
   }
   eqTotals(rowTotals(next.b),initial);totals++;({b,phi}=next);
  }
 }
 return {local_current_and_loop_source_balances:balances,retained_row_boundary_totals:totals,
  net_boundary_residue_created_by_loop_drive:false,physical_electric_current_identified:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cp=add(hh,kk),cm=add(hh,sc(-1,kk));
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['source_census_square',mul(cp,cp),sc(2,ii)],['reverse_role_square',mul(cm,cm),sc(2,ii)],
 ['source_and_reverse_role_anticommute',add(mul(cp,cm),mul(cm,cp)),sc(0,ii)],
 ['plus_cut_square',mul(pp,pp),pp],['minus_cut_square',mul(qq,qq),qq],['role_cuts_disjoint',mul(pp,qq),sc(0,ii)],
 ['role_exchange',mul(kk,pp),mul(qq,kk)],['native_turn_square',mul(rr,rr),sc(-1,ii)],
 ['native_turn_reversal',mul(hh,rr,hh),sc(-1,rr)]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const induced=coupledStep(new Map(),unit(e0));
const negatives={
 independent_address_shifts_force_directional_transports_to_commute:!energy(compose(omega,unit(e0))).eq(0),
 separate_source_role_copies_have_the_same_nonzero_commutator:(()=>{const x=new Map([...vx].map(([k,A])=>[k,o.kron(A,I)])),y=new Map([...vy].map(([k,A])=>[k,o.kron(I,A)]));return minus(compose(x,y),compose(y,x)).size===0;})(),
 order_residue_is_the_identical_signed_closed_loop_response:minus(omega,loopDefect).size!==0,
 curvature_follows_the_same_forward_directional_update:minus(compose(omega,vx),compose(vx,omega)).size!==0,
 dual_of_two_direction_order_is_always_its_plain_inverse:minus(dual,dag(forward)).size!==0,
 straight_interface_residue_alone_forces_loop_curvature:o.isZero(at(omega,(x,y)=>x<0?e0:o.zeros(2,1),0,0)),
 intersecting_zero_field_backgrounds_cannot_have_local_loop_response:!o.isZero(at(omega,(x,y)=>x<0&&y<0?e0:o.zeros(2,1),0,-1)),
 arbitrary_native_frame_change_creates_curvature_of_pure_shifts:checks.find(x=>x.name==='local_native_frames_preserve_loop_response_and_do_not_create_curvature').position_dependent_frame_is_not_a_physical_gauge_field,
 period_two_flat_readout_proves_source_is_flat:o.isZero(cyclic[1].Omega)&&!energy(compose(omega,unit(e0))).eq(0),
 periodic_curvature_blindness_contains_only_uniform_backgrounds:cyclic[2].blind>2,
 full_loop_driven_joint_form_is_positive_before_removing_blind_periodic_modes:cyclic[2].blind>0,
 loop_residue_cannot_drive_a_nonzero_native_field:!energy(induced.b).eq(0),
 ordinary_joint_b_phi_norm_is_the_conserved_interaction_energy:!energy(induced.b).add(energy(induced.phi)).eq(1),
 loop_drive_creates_net_row_boundary_residue:(()=>{const sums=rowTotals(induced.b);return [...sums.values()].every(o.isZero);})(),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r30.native-directional-loop-interaction.v1',status:'PASS_R30_NATIVE_DIRECTIONAL_LOOP_INTERACTION',
 input_sha256:inputHash,r29_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_shared_role_loop_curvature_derived:true,loop_driven_paired_field_closure_derived:true,
 native_local_frame_covariance_derived:true,complete_cyclic_curvature_quotient_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,
 physical_metric_c_alpha_derived:false,physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,
 formal_proof_assistant_verified:false,
 scope:'Eight written native results: internal direction-block addresses, exact protocol-loop curvature, dual response transport, corner localization, closed loop-driven field interaction, native local-frame covariance, complete cyclic blind-sector reconstruction and its positive interaction quotient. No physical EM/vacuum/charge/c/alpha/mass/dimension/gauge-group identification.'},null,2)+'\n');
