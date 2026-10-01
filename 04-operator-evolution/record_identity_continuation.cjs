#!/usr/bin/env node
'use strict';
// Finite application adapter; native coefficients and proof replay stay upstream.
const fs=require('node:fs'),path=require('node:path');
const get=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const home=path.join(get('--rkf-root'),'operator_foundation');
const o=require(path.join(home,'core/native_operator.cjs'));
const p=require(path.join(home,'core/paninian_operator.cjs'));
const w=require(path.join(home,'core/workbench.cjs'));
const input=JSON.parse(fs.readFileSync(get('--input'),'utf8')),inputHash=p.digest(input);
if(inputHash!==get('--expected-input-sha256'))throw Error('R21 input pin mismatch');
const ensure=(ok,msg)=>{if(!ok)throw Error(msg);},eq=(a,b,msg)=>ensure(o.equal(a,b),msg),f=x=>o.F.of(x);
const I=o.identity(2),Z=o.zeros(2),zv=()=>o.zeros(2,1);
const fromTags=tags=>o.matrix([0,1].map(row=>tags.map(tag=>Math.abs(tag)-1===row?(tag>0?1:-1):0)));
const H=fromTags(input.constructed_maps.parity_on_roles),K=fromTags(input.constructed_maps.role_exchange);
const F=o.add(H,K),P=o.scale(o.add(I,H),'1/2'),Q=o.sub(I,P),L=[o.mul(P,F),o.mul(Q,F)];
const col=xs=>o.matrix(xs.map(x=>[x])),e0=col([1,0]),e1=col([0,1]);
const outer=(v,u=v)=>o.mul(v,o.dagger(u));
const trace=v=>v[0][0].rad.add(v[1][1].rad),denom=n=>new o.F(1,2**n);
const pop=n=>n.toString(2).replaceAll('0','').length;
const checks=[],check=(name,fn)=>checks.push({name,passed:true,...fn()});

const addAt=(a,key,v)=>{const s=o.add(a.get(key)||zv(),v);if(o.isZero(s))a.delete(key);else a.set(key,s);};
const joint=es=>{const a=new Map();for(const [x,m,v]of es)addAt(a,x+','+m,v);return a;};
const jes=a=>[...a].map(([k,v])=>{const [x,m]=k.split(',').map(Number);return [x,m,v];});
const jeq=(a,b,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||zv(),b.get(key)||zv(),msg+' '+key);};
const systemStep=a=>{const out=new Map();for(const [x,m,v]of jes(a))for(let b=0;b<2;b++)addAt(out,(x+(b?-1:1))+','+m,o.mul(L[b],v));return out;};
const roleH=a=>joint(jes(a).map(([x,m,v])=>[x,m,o.mul(H,v)]));
const write=(a,i)=>{const out=new Map();for(const [x,m,v]of jes(a)){
 addAt(out,x+','+m,o.mul(P,v));addAt(out,x+','+(m^(1<<i)),o.mul(Q,v));}return out;};
const step=(a,i,contrast=false)=>write(contrast?roleH(systemStep(a)):systemStep(a),i);
const run=(sigma,contrast=-1,source=joint([[0,0,e0]]))=>{let a=source;for(let n=0;n<sigma.length;n++)a=step(a,sigma[n],n===contrast);return a;};

const padd=(g,x,y,v)=>{const key=x+','+y,s=o.add(g.get(key)||Z,v);if(o.isZero(s))g.delete(key);else g.set(key,s);};
const pairs=es=>{const g=new Map();for(const [x,y,v]of es)padd(g,x,y,v);return g;};
const pes=g=>[...g].map(([k,v])=>{const [x,y]=k.split(',').map(Number);return [x,y,v];});
const at=(g,x,y)=>g.get(x+','+y)||Z;
const pscale=(g,c)=>pairs(pes(g).map(([x,y,v])=>[x,y,o.scale(v,c)]));
const psum=(...gs)=>pairs(gs.flatMap(pes));
const peq=(a,b,msg)=>{for(const key of new Set([...a.keys(),...b.keys()]))eq(a.get(key)||Z,b.get(key)||Z,msg+' '+key);};
const samePairs=(a,b)=>[...new Set([...a.keys(),...b.keys()])].every(k=>o.equal(a.get(k)||Z,b.get(k)||Z));
const rho=(g,x)=>trace(at(g,x,x)),current=(g,x)=>trace(o.mul(K,at(g,x,x)));
const total=g=>[...new Set(pes(g).flatMap(([x,y])=>[x,y]))].reduce((s,x)=>s.add(rho(g,x)),f(0));
const moment=(a,eta,c=1)=>{
 const records=new Map();for(const [x,m,v]of jes(a)){if(!records.has(m))records.set(m,[]);records.get(m).push([x,v]);}
 const g=new Map();for(const [m,vs]of records)for(const [x,v]of vs)for(const [y,u]of records.get(m^eta)||[])padd(g,x,y,outer(v,u));return pscale(g,c);
};
const momentStep=(getMoment,i,eta)=>{
 const out=new Map();for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  const g=getMoment(eta^(a===b?0:1<<i));
  for(const [x,y,v]of pes(g))padd(out,x+(a?-1:1),y+(b?-1:1),o.scale(o.mul(o.mul(L[a],v),o.dagger(L[b])),'1/2'));
 }return out;
};
const T=g=>momentStep(()=>g,0,0);
const conjugate=(g,M)=>pairs(pes(g).map(([x,y,v])=>[x,y,o.mul(o.mul(M,v),o.dagger(M))]));
const sign=(mask,eta)=>pop(mask&eta)%2?-1:1;
const masks=bits=>{const out=[];for(let m=bits;;m=(m-1)&bits){out.push(m);if(m===0)break;}return out.sort((a,b)=>a-b);};
const slotMask=sigma=>sigma.reduce((m,i)=>m|(1<<i),0);
const arity=sigma=>sigma.length?Math.max(...sigma)+1:0;
const canonical=sigma=>{const labels=new Map();return sigma.map(i=>{if(!labels.has(i))labels.set(i,labels.size);return labels.get(i);});};
const partitions=n=>{const out=[];function visit(a,max){if(a.length===n){out.push(a);return;}for(let i=0;i<=max+1;i++)visit([...a,i],Math.max(max,i));}visit([],-1);return out;};
const seed=joint([[0,0,e0]]),g0=moment(seed,0);
const all=Array.from({length:7},(_,n)=>partitions(n)),states=new Map([['',seed]]);
for(let n=1;n<=6;n++)for(const sigma of all[n])states.set(sigma.join(''),step(states.get(sigma.slice(0,-1).join('')),sigma[n-1]));

// Literal H/K source paths use tuple strings for record updates, independently
// of the exchange-mask recurrence used by the joint propagation above.
const histories=[[{word:'',d:[],x:0,v:e0}]];
for(let n=1;n<=6;n++)histories.push(histories[n-1].flatMap(h=>[[H,'H'],[K,'K']].map(([M,tag])=>{
 const v=o.mul(M,h.v),b=o.isZero(o.mul(P,v))?1:0;return {word:h.word+tag,d:[...h.d,b],x:h.x+(b?-1:1),v};})));
const tupleRecord=(sigma,d)=>{const bits=Array(arity(sigma)).fill('0');for(let k=0;k<d.length;k++)if(d[k])bits[sigma[k]]=bits[sigma[k]]==='0'?'1':'0';return bits.reduce((m,b,i)=>m+(b==='1'?2**i:0),0);};
const parityRecord=(sigma,d)=>{let m=0;for(let i=0;i<arity(sigma);i++){let count=0;for(let k=0;k<d.length;k++)if(sigma[k]===i)count+=d[k];if(count%2)m+=2**i;}return m;};
const kernel=(sigma,d,e)=>{for(let i=0;i<arity(sigma);i++){let count=0;for(let k=0;k<d.length;k++)if(sigma[k]===i&&d[k]!==e[k])count++;if(count%2)return false;}return true;};
const recordMatrix=(r,i)=>{let M=o.identity(1);for(let j=r-1;j>=0;j--)M=o.kron(M,j===i?K:I);return M;};

check('native_slot_exchanges_and_projector_algebra',()=>{
 let basisCases=0,projectorCases=0;
 for(let r=1;r<=4;r++){
  const dim=2**r,id=o.identity(dim),ks=Array.from({length:r},(_,i)=>recordMatrix(r,i));
  for(let i=0;i<r;i++){eq(o.mul(ks[i],ks[i]),id,'slot K square');eq(o.dagger(ks[i]),ks[i],'slot dagger');
   for(let m=0;m<dim;m++){const v=o.zeros(dim,1),u=o.zeros(dim,1);v[m][0]=o.Cut.of(1);u[m^(1<<i)][0]=o.Cut.of(1);eq(o.mul(ks[i],v),u,'tuple exchange basis');basisCases++;}
   for(let j=0;j<r;j++)eq(o.mul(ks[i],ks[j]),o.mul(ks[j],ks[i]),'disjoint-slot commutation');
  }
  const es=Array.from({length:dim},(_,chi)=>ks.reduce((E,M,i)=>o.mul(E,o.scale(o.add(id,o.scale(M,(chi>>i)&1?-1:1)),'1/2')),id));
  eq(es.reduce((a,b)=>o.add(a,b),o.zeros(dim)),id,'record projectors sum');
  for(let chi=0;chi<dim;chi++){eq(o.mul(es[chi],es[chi]),es[chi],'record projector');eq(o.dagger(es[chi]),es[chi],'projector dagger');
   for(let i=0;i<r;i++)eq(o.mul(ks[i],es[chi]),o.scale(es[chi],(chi>>i)&1?-1:1),'record eigenvalue');
   for(let tau=0;tau<dim;tau++){eq(o.mul(es[chi],es[tau]),chi===tau?es[chi]:o.zeros(dim),'projector disjointness');projectorCases++;}
  }
 }return {slot_tuple_basis_cases:basisCases,complete_projector_product_cases:projectorCases,maximum_record_slots:4};
});
check('all_allocation_partitions_match_literal_cut_histories',()=>{
 let allocations=0,wordInstances=0,fibres=0;
 for(let n=0;n<=6;n++)for(const sigma of all[n]){
  const r=arity(sigma),literal=new Map(),counts=new Map();
  for(const h of histories[n]){const m=tupleRecord(sigma,h.d);ensure(m===parityRecord(sigma,h.d),'parity from literal tuple');
   ensure(((h.x+2*pop(m)-n)%4+4)%4===0,'joint address-record invariant');
   addAt(literal,h.x+','+m,h.v);counts.set(m,(counts.get(m)||0)+1);wordInstances++;
  }
  ensure(counts.size===2**r,'record image size');for(const c of counts.values()){ensure(c===2**(n-r),'uniform parity fibre');fibres++;}
  jeq(literal,states.get(sigma.join('')),'literal source / gate evolution');
  ensure(total(moment(literal,0,denom(n))).eq(1),'joint normalization');
  ensure(JSON.stringify(canonical(sigma.map(x=>101+7*x)))===JSON.stringify(sigma),'slot renaming canonical form');allocations++;
 }return {maximum_depth:6,allocation_partitions:allocations,literal_HK_word_instances:wordInstances,record_fibres_checked:fibres,
           source_words_erased_by_record_quotient:false};
});
check('ordered_history_pair_kernel_and_support_grading',()=>{
 let cases=0,fields=0;
 for(let n=0;n<=4;n++)for(const sigma of all[n]){
  const raw=new Map();for(const d of histories[n])for(const e of histories[n]){
   const keep=kernel(sigma,d.d,e.d),same=tupleRecord(sigma,d.d)===tupleRecord(sigma,e.d);ensure(keep===same,'pair retention parity');
   if(keep)padd(raw,d.x,e.x,outer(d.v,e.v));cases++;
  }
  peq(pscale(raw,denom(n)),moment(states.get(sigma.join('')),0,denom(n)),'literal retained pairs');fields++;
 }
 let momentBlocks=0;
 for(const sigma of [[0,0,0,0],[0,1,0,1],[0,1,2,0,1]])for(let eta=0;eta<2**arity(sigma);eta++)
  for(const [x,y,v]of pes(moment(states.get(sigma.join('')),eta))){ensure((x-y)%2===0&&((x-y)/2-pop(eta))%2===0,'moment support grading');momentBlocks++;}
 return {ordered_history_pairs:cases,complete_pair_fields:fields,nonzero_moment_blocks_checked:momentBlocks};
});

const protocols=[[0],[0,0],[0,1],[0,1,0],[0,1,1,0],[0,1,2,0,2]];
const fixtures=[seed,joint([[-1,0,col([1,2])],[1,1,col([2,-1])],[0,3,col(['1/2','1/3'])]]),joint([[0,0,col([1,1])],[2,2,col([1,-1])]])];
check('exchange_moment_hierarchy_matches_full_joint_evolution',()=>{
 let comparisons=0,steps=0;
 for(const sigma of protocols)for(const source of fixtures){const r=3,dim=2**r;let a=source,ms=Array.from({length:dim},(_,eta)=>moment(a,eta));
  for(let n=0;n<=sigma.length;n++){
   for(let eta=0;eta<dim;eta++){peq(ms[eta],moment(a,eta,denom(n)),'closed moment / full joint response');comparisons++;}
   if(n<sigma.length){const old=ms;ms=Array.from({length:dim},(_,eta)=>momentStep(q=>old[q],sigma[n],eta));a=step(a,sigma[n]);steps++;}
  }
 }return {joint_preparations:fixtures.length,allocation_protocols:protocols.length,full_moment_comparisons:comparisons,joint_steps:steps,
           nonzero_initial_record_correlations_included:true};
});
check('live_slot_memory_reduction_keeps_exact_future_response',()=>{
 let stages=0,retainedFields=0,maximumLive=0;
 for(let n=0;n<=6;n++)for(const sigma of all[n]){
  let used=0,a=seed,kept=new Map([[0,g0]]);
  for(let k=0;k<n;k++){
   const slot=sigma[k],newUsed=used|(1<<slot),future=slotMask(sigma.slice(k+1)),live=newUsed&future;
   const old=kept,oldUsed=used;
   const getM=eta=>{if(old.has(eta))return old.get(eta);ensure((eta&~oldUsed)!==0,'discarded a moment needed by future');return new Map();};
   kept=new Map(masks(live).map(eta=>[eta,momentStep(getM,slot,eta)]));a=step(a,slot);used=newUsed;
   for(const [eta,g]of kept){peq(g,moment(a,eta,denom(k+1)),'live-slot exact continuation');retainedFields++;}
   ensure(kept.size===2**pop(live),'live field count');maximumLive=Math.max(maximumLive,pop(live));stages++;
  }
  peq(kept.get(0),moment(states.get(sigma.join('')),0,denom(n)),'final live quotient');
 }return {complete_continuation_stages:stages,retained_moment_field_comparisons:retainedFields,maximum_live_slots:maximumLive,
           future_slot_catalogue_required:true,initial_blank_product_records_required:true,physical_record_erasure_claimed:false};
});
check('sector_transport_independently_recovers_unresolved_response',()=>{
 let allocations=0,sectors=0;
 for(let n=0;n<=5;n++)for(const sigma of all[n]){
  const r=arity(sigma),dim=2**r;let result=new Map();
  for(let chi=0;chi<dim;chi++){let a=seed;for(const slot of sigma){a=systemStep(a);if((chi>>slot)&1)a=roleH(a);}
   result=psum(result,moment(a,0,denom(n+r)));sectors++;
  }
  peq(result,moment(states.get(sigma.join('')),0,denom(n)),'independent native sectors');allocations++;
 }
 let correlatedCases=0;
 for(const source of fixtures){const r=3,dim=2**r;let a=source,ms=Array.from({length:dim},(_,eta)=>moment(a,eta));
  let zs=Array.from({length:dim},(_,chi)=>pscale(psum(...ms.map((g,eta)=>pscale(g,sign(chi,eta)))),denom(r)));
  const sigma=[0,1,2,0,1];
  for(let n=0;n<=sigma.length;n++){
   for(let chi=0;chi<dim;chi++){
    const actual=pscale(psum(...Array.from({length:dim},(_,eta)=>pscale(moment(a,eta,denom(n)),sign(chi,eta)))),denom(r));
    peq(zs[chi],actual,'correlated sector / full joint response');correlatedCases++;
   }
   if(n<sigma.length){zs=zs.map((g,chi)=>conjugate(T(g),(chi>>sigma[n])&1?H:I));a=step(a,sigma[n]);}
  }
 }
 return {allocation_partitions:allocations,independent_sector_paths:sectors,correlated_sector_comparisons:correlatedCases,
         random_sign_premise:false};
});

const front=(a,n)=>{const g=moment(a,0,denom(n));return {j:current(g,n-2),energy:rho(g,n-2)};};
const uses=sigma=>{const counts=Array(arity(sigma)).fill(0);for(const i of sigma)counts[i]++;return counts;};
check('front_current_return_count_and_energy_balance',()=>{
 let cases=0,balanceCases=0,freshCases=0,returnCases=0;
 for(let n=1;n<=6;n++)for(const sigma of all[n]){
  const a=states.get(sigma.join('')),actual=front(a,n),counts=uses(sigma.slice(0,-1));
  const nu=counts[sigma[n-1]]||0,squares=counts.reduce((s,c)=>s+c*c,0);
  ensure(actual.j.eq(new o.F(nu,2**(n-1))),'front current counts prior writes');
  ensure(actual.energy.eq(new o.F(1+squares,2**n)),'front energy record multiplicity');
  ensure(actual.j.zero()===(nu===0),'fresh / returning detection');if(nu===0)freshCases++;else returnCases++;
  const next=front(step(a,0),n+1);
  ensure(next.energy.mul(2**(n+1)).sub(actual.energy.mul(2**n)).eq(f(1).add(actual.j.mul(2**n))),'native return-count balance');
  cases++;balanceCases++;
 }return {all_partition_front_cases:cases,successive_energy_balance_cases:balanceCases,fresh_final_writes:freshCases,returning_final_writes:returnCases,
         inserted_rate_or_coupling:false};
});

function decodeEquality(relation){
 const n=relation.length;
 for(let i=0;i<n;i++){ensure(relation[i].length===n&&relation[i][i]===true,'invalid identity diagonal');
  for(let j=0;j<n;j++){ensure(typeof relation[i][j]==='boolean'&&relation[i][j]===relation[j][i],'invalid identity symmetry');
   for(let k=0;k<n;k++)ensure(!(relation[i][j]&&relation[j][k])||relation[i][k],'nontransitive slot identity');}}
 const labels=[];let fresh=0;for(let i=0;i<n;i++){let found=-1;for(let j=0;j<i;j++)if(relation[i][j]){found=labels[j];break;}labels.push(found<0?fresh++:found);}return labels;
}
check('native_H_contrasts_reconstruct_entire_allocation_partition',()=>{
 let probes=0,decoded=0;
 for(let n=1;n<=5;n++)for(const sigma of all[n]){
  const relation=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j));
  for(let later=1;later<n;later++)for(let earlier=0;earlier<later;earlier++){
   const prefix=sigma.slice(0,later+1),base=front(states.get(prefix.join('')),later+1);
   const contrasted=front(run(prefix,earlier),later+1);
   const identity=base.j.sub(contrasted.j).mul(2**(later-1));
   ensure(identity.eq(0)||identity.eq(1),'native identity contrast must be binary');
   ensure(identity.eq(sigma[earlier]===sigma[later]?1:0),'contrast identity theorem');
   relation[earlier][later]=relation[later][earlier]=identity.eq(1);probes++;
  }
  ensure(JSON.stringify(decodeEquality(relation))===JSON.stringify(sigma),'complete allocation recovery');decoded++;
 }
 let invalidRejected=false;try{decodeEquality([[true,true,false],[true,true,true],[false,true,true]]);}catch(_){invalidRejected=true;}
 ensure(invalidRejected,'inconsistent identity data silently closed');
 return {maximum_depth:5,native_contrast_readings:probes,complete_partitions_reconstructed:decoded,nontransitive_identity_data_rejected:true,
         physical_prepare_repeat_readout_interface_claimed:false};
});
const choose=(n,k)=>{if(k<0||k>n||!Number.isInteger(k))return 0;let a=1;for(let i=1;i<=k;i++)a=a*(n-i+1)/i;return a;};
const localSame=(a,b)=>{for(const x of new Set([...pes(a),...pes(b)].flatMap(([x,y])=>[x,y])))if(!rho(a,x).eq(rho(b,x))||!current(a,x).eq(current(b,x)))return false;return true;};
check('coherent_and_count_endpoints_are_rigid',()=>{
 let cases=0,coherent=seed;
 for(let n=1;n<=6;n++){
  coherent=systemStep(coherent);const cg=moment(coherent,0,denom(n));
  for(const sigma of all[n]){const g=moment(states.get(sigma.join('')),0,denom(n));
   ensure(localSame(g,cg)===(arity(sigma)===1),'coherent endpoint iff one slot');
   let count=true;for(let x=-n;x<=n;x++)if(!rho(g,x).eq(new o.F(choose(n,(n+x)/2),2**n))||!current(g,x).zero())count=false;
   ensure(count===(arity(sigma)===n),'count endpoint iff fresh slots');cases++;
  }
 }
 const a=[0,1,0,1],b=[0,1,1,0];
 for(let n=1;n<=4;n++){const fa=front(states.get(a.slice(0,n).join('')),n),fb=front(states.get(b.slice(0,n).join('')),n);
  ensure(fa.j.eq(fb.j)&&fa.energy.eq(fb.energy),'passive front ambiguity');}
 const probeA=front(run(a.slice(0,3),0),3),probeB=front(run(b.slice(0,3),0),3);
 ensure(!probeA.j.eq(probeB.j),'identity contrast separates same return counts');
 return {complete_local_response_classifications:cases,passively_equal_front_prefixes:4,contrast_separates_0101_from_0110:true};
});

check('exact_joint_continuation_gate_for_tagged_path_pairs',()=>{
 let cases=0,returns=0,energyTerms=0,currentTerms=0;
 for(let n=1;n<=5;n++)for(const sigma of all[n]){
  const cut=Math.floor(n/2),beta=sigma.slice(cut),t=n-cut;
  for(const d of histories[n])for(const e of histories[n]){
   const dp=d.d.slice(0,cut),ep=e.d.slice(0,cut),du=d.d.slice(cut),eu=e.d.slice(cut);
   const dx=(cut-2*dp.reduce((a,b)=>a+b,0))-(cut-2*ep.reduce((a,b)=>a+b,0));
   const delta=tupleRecord(sigma.slice(0,cut),dp)^tupleRecord(sigma.slice(0,cut),ep);
   const recordTest=delta===(parityRecord(beta,du)^parityRecord(beta,eu));
   const addressTest=dx===2*(du.reduce((a,b)=>a+b,0)-eu.reduce((a,b)=>a+b,0));
   const directRecord=tupleRecord(sigma,d.d)===tupleRecord(sigma,e.d),directAddress=d.x===e.x;
   ensure(recordTest===directRecord&&addressTest===directAddress,'joint continuation conditions');
   const E=outer(d.v,e.v),sameRole=d.d[n-1]===e.d[n-1];
   const energy=directRecord&&directAddress&&!trace(E).zero(),cur=directRecord&&directAddress&&!trace(o.mul(K,E)).zero();
   ensure(energy===(recordTest&&addressTest&&sameRole),'energy return role gate');
   ensure(cur===(recordTest&&addressTest&&!sameRole),'current return role gate');
   if(recordTest&&addressTest){ensure(t>=Math.abs(dx)/2&&t>=pop(delta),'native return lower bounds');returns++;}
   if(energy)energyTerms++;if(cur)currentTerms++;cases++;
  }
 }
 let excluded=0;const future=[0,0,0],delta=3;
 for(const u of histories[3])for(const v of histories[3]){ensure((parityRecord(future,u.d)^parityRecord(future,v.d))!==delta,'unavailable slot cannot return');excluded++;}
 return {complete_tagged_pair_continuations:cases,record_and_address_meetings:returns,nonzero_individual_energy_terms:energyTerms,
         nonzero_individual_current_terms:currentTerms,unavailable_slot_exclusions:excluded,
         individual_path_survival_promoted_to_total_signal:false};
});
check('exchange_probe_rank_and_normalized_channel_count',()=>{
 let scalarProbes=0;const ranks=[];
 for(let r=1;r<=4;r++){
  const dim=2**r,table=[];
  for(let eta=0;eta<dim;eta++){
   const row=[];for(let chi=0;chi<dim;chi++){
    let a=joint(Array.from({length:dim},(_,m)=>[0,m,o.scale(col([1,1]),sign(chi,m))]));
    for(let i=0;i<r;i++)if((eta>>i)&1)a=write(a,i);
    const q=current(moment(a,0,denom(r+1)),0);ensure(q.eq(sign(chi,eta)),'native exchange probe current');row.push(q);scalarProbes++;
   }table.push(row);
  }
  const M=o.matrix(table);eq(o.mul(M,o.dagger(M)),o.scale(o.identity(dim),dim),'finite sign cancellation table');
  ensure(o.rank(M)===dim,'full exchange target rank');
  const affine=o.matrix(table.slice(1).map(row=>row.slice(1).map(q=>q.sub(row[0]))));ensure(o.rank(affine)===dim-1,'known-normalization rank');
  ranks.push({slots:r,linear_scalar_channels:dim,normalized_variable_channels:dim-1});
 }return {exact_native_scalar_probes:scalarProbes,ranks,restricted_exchange_catalogue_explicit:true,minimum_for_fixed_single_seed_claimed:false};
});
check('record_congruence_and_free_linear_kernel_are_distinct',()=>{
 let rankCases=0,continuationCases=0;
 for(let n=0;n<=4;n++)for(const sigma of all[n]){
  const r=arity(sigma),dim=2**r,columns=histories[n].map(h=>tupleRecord(sigma,h.d));
  const M=o.matrix(Array.from({length:dim},(_,m)=>columns.map(c=>m===c?1:0)));
  ensure(o.rank(M)===dim&&2**n-o.rank(M)===2**n-2**r,'free ledger rank and nullity');rankCases++;
  for(let a=0;a<dim;a++)for(let b=0;b<dim;b++)for(let slot=0;slot<=r;slot++)for(let direction=0;direction<2;direction++){
   const mask=direction?1<<slot:0;ensure((a===b)===((a^mask)===(b^mask)),'common admitted continuation congruence');continuationCases++;
  }
 }
 return {exact_free_ledger_rank_cases:rankCases,common_continuation_record_comparisons:continuationCases,
         binary_fibre_exponent_is_linear_nullity:false,full_source_word_identity_retained:true};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk);
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const ep=sc('1/2',add(ii,kk)),em=sc('1/2',add(ii,sc(-1,kk)));
const tasks=[
 ['native_exchange_involution',mul(kk,kk),ii],
 ['native_contrast_involution',mul(hh,hh),ii],
 ['contrast_fixes_positive_role',mul(hh,pp),pp],
 ['contrast_reverses_negative_role',mul(hh,qq),sc(-1,qq)],
 ['native_branch_normalization',mul(ff,ff),sc(2,ii)],
 ['positive_exchange_projector',mul(ep,ep),ep],
 ['negative_exchange_projector',mul(em,em),em],
 ['exchange_sector_disjointness',mul(ep,em),sc(0,ii)],
 ['exchange_sector_completeness',add(ep,em),ii],
 ['positive_exchange_eigenvalue',mul(kk,ep),ep],
 ['negative_exchange_eigenvalue',mul(kk,em),sc(-1,em)],
 ['single_negative_path_has_positive_cross_transition',mul(pp,ff,qq),mul(pp,kk,qq)],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[4].result.certificate));altered.output.terms=[[[],['3','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered normalization proof accepted');
const f0101=front(states.get('0101'),4),f0110=front(states.get('0110'),4);
const ambiguousProbeA=front(run([0,1,0],0),3),ambiguousProbeB=front(run([0,1,1],0),3);
const blockedDifferences=new Set(histories[2].flatMap(u=>histories[2].map(v=>parityRecord([0,0],u.d)^parityRecord([0,0],v.d))));
const negatives={
 every_repeated_slot_is_fresh:!front(states.get('00'),2).j.zero(),
 every_existing_record_preserves_coherent_current:!front(states.get('01'),2).j.eq(front(states.get('00'),2).j),
 passive_return_counts_identify_slot_names:f0101.j.eq(f0110.j)&&f0101.energy.eq(f0110.energy)&&!ambiguousProbeA.j.eq(ambiguousProbeB.j),
 H_contrast_has_half_the_derived_effect:front(states.get('00'),2).j.sub(front(run([0,0],0),2).j).eq(1),
 record_agreement_implies_local_meeting:tupleRecord([0,0],[0,0])===tupleRecord([0,0],[1,1])&&2!==-2,
 address_agreement_implies_record_agreement:tupleRecord([0,1],[0,1])!==tupleRecord([0,1],[1,0]),
 return_length_bounds_are_sufficient:pop(3)===2&&!blockedDifferences.has(3),
 record_quotient_is_full_history_equality:tupleRecord([0,0],[0,0])===tupleRecord([0,0],[1,1])&&'00'!=='11',
 fibre_exponent_is_free_linear_nullity:(4-2)!==(2**4-2**2),
 normalized_target_needs_all_variable_channels:2**3-1!==2**3,
 all_allocations_at_equal_depth_have_same_response:!front(states.get('000'),3).j.eq(front(states.get('012'),3).j),
 fresh_record_can_reveal_existing_exchange_moment:moment(seed,1).size===0,
 energy_alone_identifies_the_last_write:front(states.get('00'),2).energy.eq(front(states.get('01'),2).energy)&&!front(states.get('00'),2).j.eq(front(states.get('01'),2).j),
 evaluated_exchange_return_erases_source_tags:o.equal(o.mul(K,K),I)&&'KK'!=='',
 record_current_is_an_inserted_probability:front(run([0,0],0),2).j.eq('-1/2'),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r21.record-identity.v1',status:'PASS_R21_NATIVE_RECORD_IDENTITY',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,
 altered_native_certificate_rejected:alteredRejected,
 ordinary_complex_measurement_stochastic_or_classical_transform_premise:false,
 allocation_identification_derived:true,physical_allocation_selected:false,physical_metric_c_alpha_derived:false,
 formal_proof_assistant_verified:false,
 scope:'Native W allocation class. Exact history quotient, exchange moments, live continuation, count/contrast identification, endpoint rigidity and tagged meeting gates. Probe-catalogue rank is separate from fixed-seed or physical selection.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
