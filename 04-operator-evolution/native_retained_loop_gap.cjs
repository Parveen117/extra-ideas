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
ensure(inputHash===arg('--expected-input-sha256'),'R32 input pin mismatch');
const parentBytes=fs.readFileSync(arg('--r31-certificate')),parentHash=hash(parentBytes);
ensure(parentHash===arg('--expected-r31-sha256'),'R32 R31 certificate pin mismatch');
ensure(JSON.parse(parentBytes).status==='PASS_R31_NATIVE_PROPAGATION_GEOMETRY','R32 parent status');
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

check('native_source_law_record_has_a_nonflat_reciprocal_loop',()=>{
 eqFields(plus(compose(X,Y),compose(Y,X)),new Map(),4,4,'native record makes link shifts anticommute');
 const closed=compose(compose(compose(X,Y),dag(X)),dag(Y));
 eqFields(closed,unit(o.scale(I4,-1)),4,4,'record plaquette is minus identity');
 const order=minus(compose(X,Y),compose(Y,X));eqFields(compose(dag(order),order),unit(o.scale(I4,4)),4,4,'all-preparation native loop residue');
 for(const A of [X,Y,VX,VY,W,B])eqFields(compose(dag(A),A),unit(I4),4,4,'native link/transport pairing');
 eq(o.mul(o.dagger(C0),C0),o.scale(I,2),'same source-law census on memory');
 return {complete_Laurent_loop_and_pairing_identities:9,source_census_normalization_identity:1,
  record_loop:'-I',loop_order_residue_squared:'4 I',new_record_law_or_angle_supplied:false};
});

check('literal_fine_source_histories_reproduce_coarse_retained_transport',()=>{
 let prefixes=0,endpoints=0;
 for(const source of [e0,e1])for(const memory of [e0,e1])for(const swapped of [false,true]){
  let histories=[{x:0,y:0,v:source,mu:memory}],direct=unit(o.kron(source,memory));
  for(let n=1;n<=8;n++){
   const axis=(Math.floor((n-1)/2)%2)^(swapped?1:0),next=[];
   for(const old of histories)for(const letter of [H,K]){
    const v=o.mul(letter,old.v),sign=o.isZero(o.mul(P,v))?-1:1,coordinate=axis?old.y:old.x;
    const arrow=parity(Math.min(coordinate,coordinate+sign))?(axis?K:H):I;
    next.push({x:old.x+(axis?0:sign),y:old.y+(axis?sign:0),v,mu:o.mul(arrow,old.mu)});
   }
   histories=next;prefixes+=histories.length;
   if(n%4===0){const census=new Map();for(const h of histories){ensure(h.x%2===0&&h.y%2===0,'even endpoint address');addAt(census,key(h.x/2,h.y/2),o.scale(o.kron(h.v,h.mu),new o.F(1,1n<<BigInt(n/2))));}
    direct=compose(swapped?compose(VX,VY):W,direct);eqFields(direct,census,4,1,'literal source/edge-memory census');endpoints++;
   }
  }
 }
 return {literal_source_prefixes:prefixes,independent_endpoint_comparisons:endpoints,max_literal_event:8,
  memory_written_on_alternating_fine_edges:true,coarse_covariant_transport_is_source_derived:true};
});

const staggeredY=a=>{const out=new Map();for(const [k,v]of a){const [x,y]=xy(k);for(const [j,M]of vy){const [u,d]=xy(j);addAt(out,key(x+u,y+d),o.scale(o.mul(M,v),d&&parity(x)?-1:1));}}return out;};
check('area_parity_normal_form_changes_an_actual_source_reading',()=>{
 let joint=unit(I4),signed=unit(I),blocks=0,energies=0;const preps=[e0,e1,col([1,2])];
 for(let n=0;n<=8;n++){
  const expected=new Map();for(const [k,M]of signed){const [x,y]=xy(k);expected.set(k,o.kron(M,recordFrame(x,y)));}
  eqFields(joint,expected,4,4,'source word normal form with retained endpoint record');blocks+=joint.size;
  for(const mu of preps){const full=compose(joint,unit(o.kron(e0,mu)));ensure(energy(full).eq(o.energy(mu)),'complete memory source norm');energies++;}
  if(n<8){joint=compose(W,joint);signed=staggeredY(compose(vx,signed));}
 }
 const seed=o.kron(e0,e0),retained=compose(B,unit(seed)),flat=compose(compose(bareW,bareW),unit(seed));
 const a=o.energy(retained.get('0,0')),b=o.energy(flat.get('0,0'));ensure(a.eq('1/64')&&b.eq('9/64'),'two-block source return contrast');
 let rectangles=0;for(let a=1;a<=5;a++)for(let b=1;b<=5;b++){
  const loop=o.mul(o.mul(o.mul(o.power(H,a),o.power(K,b)),o.power(H,a)),o.power(K,b));eq(loop,o.scale(I,(a*b)%2?-1:1),'rectangular native area parity');rectangles++;
 }
 return {normal_ordered_endpoint_operator_blocks:blocks,complete_source_norm_checks:energies,max_four_event_blocks:8,
  rectangular_loop_parity_checks:rectangles,origin_response_after_two_blocks:{retained:'1/64',flat:'9/64'},
  record_preparation_independence_scope:'one-origin product input on the unwrapped count plane',record_multiplicity_can_be_factored_with_staggered_signs:true};
});

check('retained_loop_derives_a_gapped_paired_wave_and_sharp_bounds',()=>{
 eqFields(plus(B,dag(B)),coefficient,4,4,'complete paired wave coefficient');
 eqFields(plus(minus(compose(B,B),compose(coefficient,B)),unit(I4)),new Map(),4,4,'fourth-degree W identity');
 eqFields(compose(dag(F),F),G,4,4,'factorized positive native defect');
 eqFields(G,plus(plus(unit(I4),scale(plus(LX,LY),'1/4')),scale(compose(LX,LY),'1/16')),4,4,'complete difference-square ledger');
 let squares=0,waves=0,sharp=0;
 for(const seed of probes){let old=seed,now=compose(B,seed);
  const right=energy(seed).add(energy(compose(DX,seed)).div(4)).add(energy(compose(DY,seed)).div(4)).add(energy(compose(compose(DX,DY),seed)).div(16));
  ensure(energy(compose(F,seed)).eq(right),'positive defect sum');squares++;
  for(let n=1;n<=5;n++){const next=compose(B,now);eqFields(plus(minus(next,scale(now,2)),old),scale(compose(G,now),-1),4,1,'native eight-event paired wave');old=now;now=next;waves++;}
 }
 for(const M of [1,2,4,8])for(const alternating of [false,true])for(const v of [o.kron(e0,e0),col([1,-2,3,1])]){
  const packet=new Map();for(let a=1;a<=M;a++)for(let b=1;b<=M;b++)packet.set(key(2*a,2*b),o.scale(v,alternating&&(a+b)%2?-1:1));
  const ratio=energy(compose(F,packet)).div(energy(packet)),target=(alternating?f(2).sub(new o.F(1,2*M)):f(1).add(new o.F(1,2*M))).pow(2);
  ensure(ratio.eq(target)&&f(1).le(ratio)&&ratio.le(4),'sharp native packet family');sharp++;
 }
 return {complete_Laurent_wave_and_Gram_identities:4,finite_difference_square_checks:squares,finite_paired_wave_checks:waves,
  sharp_packet_checks:sharp,sharp_lower_defect:'1',sharp_upper_defect:'4',physical_mass_gap_identified:false};
});

const ring=(op,L,r=4,c=4)=>{const M=o.zeros(r*L*L,c*L*L);for(const [k,A]of op){const [dx,dy]=xy(k);
 for(let x=0;x<L;x++)for(let y=0;y<L;y++){const u=((x-dx)%L+L)%L,v=((y-dy)%L+L)%L,a=x*L+y,b=u*L+v;
  for(let i=0;i<r;i++)for(let j=0;j<c;j++)M[a*r+i][b*c+j]=M[a*r+i][b*c+j].add(A[i][j]);
 }}return M;};
const cyclic=[];
check('native_gap_obstructs_flat_transport_similarity_and_local_frame_removal',()=>{
 const rows=[];let frames=0;
 for(let L=1;L<=4;L++){
  const n=4*L*L,id=o.identity(n),fm=ring(F,L),wm=ring(W,L),gm=ring(G,L),flat=ring(minus(compose(bareW,bareW),unit(I4)),L);
  ensure(o.rank(fm)===n&&o.rank(flat)<n,'retained versus flat stationary obstruction');eq(o.mul(o.dagger(fm),fm),gm,'finite gap Gram');
  rows.push({period:L,retained_increment_rank:n,flat_increment_rank:o.rank(flat)});cyclic.push({L,n,id,fm,wm,gm});
 }
 const L=3,n=4*L*L,id=o.identity(n),D=o.zeros(n);
 for(let z=0;z<L*L;z++){const a=z+1,b=z%3,den=a*a+b*b,turn=o.add(o.scale(I,new o.F(a*a-b*b,den)),o.scale(R,new o.F(2*a*b,den))),g=o.kron(I,z%2?o.mul(turn,H):turn);
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)D[4*z+i][4*z+j]=g[i][j];}
 const change=A=>o.mul(o.mul(D,A),o.dagger(D)),xx=change(ring(X,L)),yy=change(ring(Y,L));
 eq(o.mul(o.dagger(D),D),id,'native local frame pairing');frames++;
 eq(o.mul(o.mul(o.mul(xx,yy),o.dagger(xx)),o.dagger(yy)),o.scale(id,-1),'minus-identity loop survives every frame');frames++;
 const ff=change(ring(F,L));eq(o.mul(o.dagger(ff),ff),change(ring(G,L)),'gap response is covariant');frames++;
 const a=ring(X,3),b=ring(TY,3);eq(o.commutator(a,b),o.zeros(n),'flat plaquettes can retain winding memory');
 ensure(!o.equal(o.power(a,3),id),'odd cyclic winding is not an identity frame');
 return {count_quotients:rows,exact_local_frame_identities:frames,periodic_flatness_also_requires_trivial_winding:true,
  bounded_time_independent_similarity_to_flat_source_excluded:true,all_factorizations_or_time_dependent_frames_excluded:false};
});

check('uniform_native_source_has_an_exact_twelve_block_cycle',()=>{
 const w0=ring(W,1),b0=o.mul(w0,w0),ranks=[];
 eq(o.sub(o.add(o.mul(b0,b0),I4),b0),o.zeros(4),'uniform internal cycle polynomial');
 eq(o.power(w0,6),o.scale(I4,-1),'uniform six-block sign reversal');eq(o.power(w0,12),I4,'uniform twelve-block return');
 for(let n=1;n<=12;n++){const rank=o.rank(o.sub(o.power(w0,n),I4));ensure(rank===(n===12?0:4),'least full signed period');ranks.push(rank);}
 return {cycle_polynomial_identities:3,increment_ranks_for_blocks_1_through_12:ranks,least_full_signed_block_period:12,
  original_source_events_per_full_return:48,physical_frequency_or_particle_mass_assigned:false};
});

check('signed_retained_increment_has_an_exact_inverse_and_uniform_decoder_error',()=>{
 let inverses=0,identities=0,bounds=0;const rows=[];
 for(const {L,n,id,fm,gm}of cyclic){if(L>3)continue;
  const qx=ring(QX,L),qy=ring(QY,L),ax=o.scale(ring(pairAverageX,L),'1/3'),ay=o.scale(ring(pairAverageY,L),'1/3');
  const decoder=o.mul(o.mul(o.inverse(qy),o.inverse(qx)),o.dagger(fm));eq(o.mul(decoder,fm),id,'complete signed source decoder');inverses++;
  const seed=col(Array.from({length:n},(_,i)=>(i%7)-3)),reading=o.mul(fm,seed);let px=id,py=id,sx=o.zeros(n),sy=o.zeros(n);
  for(let N=0;N<=5;N++){
   sx=o.add(sx,px);sy=o.add(sy,py);px=o.mul(ax,px);py=o.mul(ay,py);
   const finite=o.scale(o.mul(o.mul(sy,sx),o.dagger(fm)),'4/9');
   eq(o.mul(finite,fm),o.mul(o.sub(id,px),o.sub(id,py)),'factorized finite decoder remainder');identities++;
   const r=f(1).div(f(3).pow(N+1)),bound=r.mul(2).add(r.pow(2)),error=o.energy(o.sub(o.mul(finite,reading),seed));
   ensure(error.le(bound.pow(2).mul(o.energy(seed))),'source-norm finite decoder error');
   ensure(error.le(bound.pow(2).mul(o.energy(reading))),'observable-only finite decoder error');bounds+=2;
  }
  rows.push({period:L,source_roles:n,increment_observer_rank:n,forgotten_roles:0,finite_decoder_depths:6});
 }
 return {count_quotients:rows,exact_inverse_identities:inverses,complete_finite_remainder_identities:identities,certified_error_bound_checks:bounds,
  signed_record_roles_retained:true,intensity_only_observer_claimed_complete:false};
});

let rootBracket;
check('native_inverse_kernel_has_a_derived_exponential_count_tail',()=>{
 // C0 is an algebraic witness for the quadratic root relation only.
 // The positive scalar branch is selected separately by native rational cuts.
 const rho=o.sub(o.scale(I,3),o.scale(C0,2)),c=o.scale(C0,'1/2');
 eq(o.add(o.sub(o.mul(rho,rho),o.scale(rho,6)),I),o.zeros(2),'derived rho polynomial');
 eq(o.mul(c,o.add(I,rho)),o.sub(I,rho),'unit total inverse weight relation');
 const grad2=minus(compose(tx,tx),unit(I)),q=plus(unit(I),scale(compose(dag(grad2),grad2),'1/4'));let identities=0;
 for(let N=1;N<=6;N++){
  const kernel=new Map();for(let k=-N;k<=N;k++)kernel.set(key(2*k,0),o.mul(c,o.power(rho,Math.abs(k))));
  const expected=unit(I),inner=o.scale(o.mul(c,o.power(rho,N+1)),'1/4'),outer=o.scale(o.mul(c,o.power(rho,N)),'-1/4');
  for(const sign of [-1,1]){addAt(expected,key(sign*2*N,0),inner);addAt(expected,key(sign*2*(N+1),0),outer);}
  eqFields(compose(q,kernel),expected,2,2,'exact truncated native inverse kernel residual');identities++;
 }
 let lo=f(1),hi=f('3/2');for(let n=0;n<64;n++){const mid=lo.add(hi).div(2);if(mid.pow(2).le(2))lo=mid;else hi=mid;}
 ensure(lo.pow(2).le(2)&&f(2).le(hi.pow(2)),'positive native root cut');
 const rlo=f(3).sub(hi.mul(2)),rhi=f(3).sub(lo.mul(2));ensure(f(0).le(rlo)&&!rlo.eq(0)&&rhi.le(1)&&!rhi.eq(1),'native decaying root branch');
 rootBracket={sqrt2_lower:lo.toString(),sqrt2_upper:hi.toString(),rho_lower:rlo.toString(),rho_upper:rhi.toString()};
 return {native_quadratic_and_normalization_identities:2,complete_truncated_kernel_identities:identities,positive_root_bisections:64,
  selected_root:'3 - 2 sqrt(2)',positive_root_interval:rootBracket,physical_coupling_constant_identified:false};
});

check('retained_interaction_has_finite_native_propagation_range',()=>{
 let psi=unit(o.kron(e0,e0)),count=0;for(let n=0;n<=8;n++){
  ensure([...psi.keys()].every(k=>xy(k).every(a=>Math.abs(a)<=n)),'retained source count cone');ensure(energy(psi).eq(1),'retained total matching norm');count+=2;
  if(n<8)psi=compose(W,psi);
 }
 return {finite_support_and_norm_checks:count,max_four_event_blocks:8,count_radius_per_block:1,new_random_noise_or_external_force_supplied:false};
});

const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value}),add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),rr=word('R'),kk=word('K'),hh=sc(-1,word('R','K')),cp=add(hh,kk),cm=add(hh,sc(-1,kk)),rhoWord=add(sc(3,ii),sc(-2,cp));
const spec=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation,system=w.presentation(spec),audit=system.audit();
ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','canonical presentation');
const tasks=[['source_H_square',mul(hh,hh),ii],['source_K_square',mul(kk,kk),ii],['source_mixed_law',add(mul(hh,kk),mul(kk,hh)),sc(0,ii)],
 ['source_census_square',mul(cp,cp),sc(2,ii)],['reverse_role_square',mul(cm,cm),sc(2,ii)],['native_turn_square',mul(rr,rr),sc(-1,ii)],
 ['native_record_closed_square',mul(hh,kk,hh,kk),sc(-1,ii)],['native_record_reversed_square',mul(kk,hh,kk,hh),sc(-1,ii)],
 ['turn_reverses_census',mul(rr,cp),sc(-1,cm)],['census_then_turn',mul(cp,rr),cm],
 ['turn_reverses_reflection',mul(rr,cm),cp],['reflection_then_turn',mul(cm,rr),sc(-1,cp)],
 ['inverse_tail_quadratic',add(mul(rhoWord,rhoWord),sc(-6,rhoWord),ii),sc(0,ii)],
 ['inverse_tail_normalization',mul(sc('1/2',cp),add(ii,rhoWord)),add(ii,sc(-1,rhoWord))]];
const replays=[];for(const [name,left,right]of tasks){const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native replay '+name);replays.push({name,result,replay:'REPLAY_MATCH'});}
const altered=JSON.parse(JSON.stringify(replays[0].result.certificate));altered.output.terms=[[[],['99','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native certificate accepted');
const negatives={
 reciprocal_edges_force_every_two_direction_loop_flat:zeroOp(plus(compose(compose(compose(X,Y),dag(X)),dag(Y)),unit(I4))),
 one_local_frame_can_flatten_both_retained_shifts:!zeroOp(minus(compose(X,Y),compose(Y,X))),
 nonflat_record_cannot_be_factored_into_staggered_source_and_multiplicity:checks.find(r=>r.name==='area_parity_normal_form_changes_an_actual_source_reading').record_multiplicity_can_be_factored_with_staggered_signs,
 native_record_leaves_every_source_response_unchanged:checks.find(r=>r.name==='area_parity_normal_form_changes_an_actual_source_reading').origin_response_after_two_blocks.retained!=='9/64',
 record_preparation_independence_requires_erasing_the_record:checks.find(r=>r.name==='area_parity_normal_form_changes_an_actual_source_reading').normal_ordered_endpoint_operator_blocks>0,
 retained_paired_increment_has_the_flat_uniform_kernel:o.rank(cyclic[0].fm)===4,
 derived_defect_has_no_positive_identity_term:o.equal(ring(G,1),I4),
 mixed_difference_square_can_be_omitted:!zeroOp(minus(G,plus(unit(I4),scale(plus(LX,LY),'1/4')))),
 uniform_native_source_is_stationary:o.rank(o.sub(ring(W,1),I4))===4,
 full_signed_uniform_cycle_is_six_blocks:o.equal(o.power(ring(W,1),6),o.scale(I4,-1)),
 periodic_flatness_needs_no_winding_check:o.isZero(o.commutator(ring(X,3),ring(TY,3)))&&!o.equal(o.power(ring(X,3),3),o.identity(36)),
 complete_signed_increment_still_forgets_the_uniform_record:o.rank(cyclic[1].fm)===cyclic[1].n,
 one_generic_finite_decoder_term_is_exact:!o.equal(o.scale(cyclic[2].gm,'4/9'),cyclic[2].id),
 growing_root_is_an_admissible_localized_inverse_tail:f(3).add(f(2).mul(f(rootBracket.sqrt2_lower))).le(1)===false,
 unresolved_source_intensities_recover_an_arbitrary_unknown_record:(()=>{const a=compose(B,unit(o.kron(e0,e0))),b=compose(B,unit(o.kron(e0,e1)));return !zeroOp(minus(a,b))&&[...new Set([...a.keys(),...b.keys()])].every(k=>o.energy(a.get(k)||o.zeros(4,1)).eq(o.energy(b.get(k)||o.zeros(4,1))));})(),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
process.stdout.write(JSON.stringify({schema:'extra-ideas.r32.native-retained-loop-gap.v1',status:'PASS_R32_NATIVE_RETAINED_LOOP_GAP',
 input_sha256:inputHash,r31_native_certificate_sha256:parentHash,source_foundation:input.source_foundation,
 exact_checks:checks,check_count:checks.length,symbolic_presentations:[{name:'canonical_native_source',specification:spec,audit}],symbolic_replays:replays,symbolic_replay_count:replays.length,
 rejected_false_alternatives:negatives,altered_native_certificate_rejected:alteredRejected,
 native_reciprocal_loop_obstruction_derived:true,native_retained_propagation_gap_and_cycle_derived:true,
 native_source_response_contrast_derived:true,complete_retained_increment_inverse_and_tail_derived:true,
 ordinary_complex_or_classical_field_premise:false,physical_vacuum_identified:false,electromagnetism_derived:false,
 physical_metric_c_alpha_derived:false,physical_electric_charge_or_particle_mass_identified:false,physical_spatial_dimension_or_gauge_group_selected:false,
 primitive_physical_force_derived:false,formal_proof_assistant_verified:false,
 scope:'Seven written native results: reciprocal-loop obstruction, literal source/record transport, area-parity source response, exact gapped paired wave, obstruction to flat continuation and uniform 48-event cycle, complete retained-increment decoding and localized inverse tail. A staggered factorization remains available. No physical force/EM/mass/c/alpha identification.'},null,2)+'\n');
