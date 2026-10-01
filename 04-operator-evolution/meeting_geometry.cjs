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
if(inputHash!==get('--expected-input-sha256'))throw Error('R22 input pin mismatch');
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
const mod=(n,k)=>((n%k)+k)%k;
const bitWords=t=>Array.from({length:2**t},(_,j)=>Array.from({length:t},(_,k)=>(j>>k)&1));
const words=(alphabet,t)=>t===0?[[]]:words(alphabet,t-1).flatMap(a=>alphabet.map(i=>[...a,i]));
const count=a=>a.reduce((s,b)=>s+b,0);
const rec=(beta,dirs)=>{const tuple=[];for(let k=0;k<beta.length;k++)if(dirs[k])tuple[beta[k]]=1-(tuple[beta[k]]||0);return tuple.reduce((s,b,i)=>s+(b||0)*2**i,0);};
const capacity=(beta,d,delta)=>{
 const counts=new Map();for(const i of beta)counts.set(i,(counts.get(i)||0)+1);
 for(let i=0;i<=Math.max(arity(beta),delta.toString(2).length);i++)if(((delta>>i)&1)>(counts.get(i)||0))return {ok:false};
 const w=pop(delta),Q=[...counts].reduce((s,[i,c])=>s+c-mod(c-((delta>>i)&1),2),0);
 return {ok:mod(d-w,2)===0&&Math.abs(d)<=Q,Q,w,counts};
};
const admits=(beta,d,delta)=>capacity(beta,d,delta).ok;
const capacityWitness=(beta,d,delta)=>{
 const data=capacity(beta,d,delta);if(!data.ok)return null;
 let extra=(Math.max(Math.abs(d),data.w)-data.w)/2;
 const chosen=new Map();for(const [i,c]of data.counts){const low=(delta>>i)&1,take=Math.min(extra,Math.floor((c-low)/2));chosen.set(i,low+2*take);extra-=take;}
 ensure(extra===0,'witness capacity allocation');let plus=(Math.max(Math.abs(d),data.w)+d)/2;
 const u=[],v=[];for(const i of beta){if(chosen.get(i)>0){chosen.set(i,chosen.get(i)-1);u.push(plus>0?1:0);v.push(plus>0?0:1);if(plus>0)plus--;}else{u.push(0);v.push(0);}}
 return {u,v};
};
const hasEnergy=(beta,d,delta)=>beta.length>0&&admits(beta.slice(0,-1),d,delta);
const hasCurrent=(beta,d,delta)=>beta.length>0&&[-1,1].some(s=>admits(beta.slice(0,-1),d+s,delta^(1<<beta.at(-1))));
const distance=(d,delta,J)=>{
 if(d===0&&delta===0)return 0;
 if(!J||((delta&~J)!==0)||mod(d-pop(delta),2)!==0)return Infinity;
 return Math.max(Math.abs(d),pop(delta));
};
const exhaustive=(beta)=>{
 const meeting=new Set(),energy=new Set(),curr=new Set(),bs=bitWords(beta.length);let cases=0;
 for(const u of bs)for(const v of bs){const key=(count(u)-count(v))+','+(rec(beta,u)^rec(beta,v));meeting.add(key);
  if(beta.length)(u.at(-1)===v.at(-1)?energy:curr).add(key);cases++;}
 return {meeting,energy,curr,cases};
};

check('native_relative_moves_and_parity_from_role_maps',()=>{
 let transitions=0;
 eq(F,o.matrix([[1,1],[1,-1]]),'evaluated source branch table');eq(o.mul(F,F),o.scale(I,2),'source normalization');
 for(let i=0;i<3;i++)for(let d=-2;d<=2;d++)for(let delta=0;delta<8;delta++)
 for(let a=0;a<2;a++)for(let b=0;b<2;b++)for(let u=0;u<2;u++)for(let v=0;v<2;v++){
  const coefficient=F[u][a].mul(F[v][b]);ensure(!coefficient.rad.zero(),'branch coefficient nonzero');
  const nextD=(2*d+(1-2*u)-(1-2*v))/2,nextDelta=delta^((u^v)<<i);
  ensure(nextD===d+v-u,'relative displacement');ensure(mod(nextD-pop(nextDelta),2)===mod(d-pop(delta),2),'relative parity');
  ensure((nextDelta&~(1<<i))===(delta&~(1<<i)),'untouched record bits');transitions++;
 }
 return {native_role_pair_transitions:transitions,derived_F_entries:true,primitive_complex_chart_used:false};
});

const wordFixtures=Array.from({length:6},(_,t)=>words([0,1,2],t)).flat();
const enumerations=wordFixtures.map(beta=>exhaustive(beta));
check('fixed_word_capacity_matches_exhaustive_direction_pairs',()=>{
 let decisions=0,pathPairs=0;for(let z=0;z<wordFixtures.length;z++){
  const beta=wordFixtures[z],brute=enumerations[z];pathPairs+=brute.cases;
  for(let d=-beta.length-1;d<=beta.length+1;d++)for(let delta=0;delta<16;delta++){
   ensure(admits(beta,d,delta)===brute.meeting.has(d+','+delta),'capacity / full direction enumeration');decisions++;
  }
 }
 return {allocation_words:wordFixtures.length,maximum_length:5,available_slots:3,record_bits_tested:4,
         exhaustive_direction_pairs:pathPairs,meeting_decisions:decisions};
});
check('capacity_solver_constructs_every_admitted_witness',()=>{
 let witnesses=0,orderChecks=0;
 for(const beta of wordFixtures)for(let d=-beta.length-1;d<=beta.length+1;d++)for(let delta=0;delta<16;delta++){
  const witness=capacityWitness(beta,d,delta);ensure(Boolean(witness)===admits(beta,d,delta),'witness existence');
  if(witness){ensure(count(witness.u)-count(witness.v)===d,'constructed address meeting');
   ensure((rec(beta,witness.u)^rec(beta,witness.v))===delta,'constructed record meeting');witnesses++;}
  ensure(admits(beta,d,delta)===admits([...beta].reverse(),d,delta),'order-independent fixed horizon capacity');orderChecks++;
 }
 ensure(!admits([0,0,1,1],4,3),'parity compatible but capacity exhausted');
 return {constructed_witnesses:witnesses,reversed_word_comparisons:orderChecks,exponential_search_used_by_solver:false};
});
check('energy_current_last_step_and_first_return_gates',()=>{
 let decisions=0,prefixCases=0,initialRoleCases=0;
 for(let z=0;z<wordFixtures.length;z++){
  const beta=wordFixtures[z],brute=enumerations[z];if(!beta.length)continue;
  for(let d=-beta.length-1;d<=beta.length+1;d++)for(let delta=0;delta<16;delta++){
   ensure(hasEnergy(beta,d,delta)===brute.energy.has(d+','+delta),'energy final role gate');
   ensure(hasCurrent(beta,d,delta)===brute.curr.has(d+','+delta),'current final role gate');decisions+=2;
   let firstMeet=Infinity,firstEnergy=Infinity,firstCurrent=Infinity;
   for(let t=0;t<=beta.length;t++){
    const pre=beta.slice(0,t);if(admits(pre,d,delta))firstMeet=Math.min(firstMeet,t);
    if(hasEnergy(pre,d,delta))firstEnergy=Math.min(firstEnergy,t);
    if(hasCurrent(pre,d,delta))firstCurrent=Math.min(firstCurrent,t);
   }
   if(d||delta){ensure(firstCurrent===firstMeet,'first nontrivial meeting is current');
    ensure(firstEnergy===(firstMeet<beta.length?firstMeet+1:Infinity),'first energy one event later');}
   else{const seen=new Set();let repeat=Infinity;for(let k=0;k<beta.length;k++){if(seen.has(beta[k])){repeat=k+1;break;}seen.add(beta[k]);}
    ensure(firstEnergy===1&&firstCurrent===repeat,'initial meeting positive return law');}
   prefixCases++;
  }
 }
 for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  const pair=outer(a?e1:e0,b?e1:e0);ensure(trace(pair).eq(a===b?1:0),'zero-step energy');
  ensure(trace(o.mul(K,pair)).eq(a!==b?1:0),'zero-step current');initialRoleCases++;
 }
 return {exact_horizon_target_decisions:decisions,prefix_latency_cases:prefixCases,zero_event_role_cases:initialRoleCases};
});
check('catalogue_distance_matches_independent_graph_exploration',()=>{
 let graphStates=0,comparisons=0,witnesses=0;
 for(let r=1;r<=4;r++){
  const distances=new Map([['0,0',0]]),queue=[[0,0]];for(let at=0;at<queue.length;at++){
   const [d,delta]=queue[at],depth=distances.get(d+','+delta);if(depth===8)continue;
   for(let i=0;i<r;i++)for(const sign of [-1,1]){const next=[d+sign,delta^(1<<i)],key=next.join(',');
    if(!distances.has(key)){distances.set(key,depth+1);queue.push(next);}}
  }
  graphStates+=distances.size;
  for(let d=-9;d<=9;d++)for(let delta=0;delta<2**r;delta++){
   const got=distances.get(d+','+delta)??Infinity,want=distance(d,delta,2**r-1);
   ensure(got===(want<=8?want:Infinity),'graph shortest path / closed distance');comparisons++;
  }
 }
 for(let J=0;J<8;J++)for(let d=-5;d<=5;d++)for(let delta=0;delta<16;delta++){
  const D=distance(d,delta,J);if(!Number.isFinite(D))continue;
  const beta=[];for(let i=0;i<4;i++)if((delta>>i)&1)beta.push(i);
  if(D>beta.length){const available=Array.from({length:3},(_,i)=>i).find(i=>(J>>i)&1);while(beta.length<D)beta.push(available,available);}
  ensure(beta.length===D&&beta.every(i=>(J>>i)&1),'catalogue geodesic uses available slots');
  ensure(admits(beta,d,delta),'catalogue distance attainment');
  if(D){ensure(hasCurrent(beta,d,delta),'first geodesic current');ensure(!hasEnergy(beta,d,delta),'no earlier energy');
   const available=beta[0];ensure(hasEnergy([...beta,available],d,delta),'next-step energy');}
  witnesses++;
 }
 return {graph_record_widths:[1,2,3,4],graph_depth:8,visited_graph_states:graphStates,distance_comparisons:comparisons,
         attained_catalogue_witnesses:witnesses,empty_catalogue_included:true};
});
const binomial=(n,k)=>{if(k<0||k>n)return 0;let a=1;for(let j=1;j<=k;j++)a=a*(n-j+1)/j;return a;};
check('endpoint_metric_ball_volume_and_record_fibres',()=>{
 let triples=0,balls=0,fibres=0;
 for(let r=1;r<=3;r++){
  const points=[];for(let z=-3;z<=3;z++)for(let m=0;m<2**r;m++)if(mod(z-pop(m),2)===0)points.push([z,m]);
  const D=(p,q)=>distance(p[0]-q[0],p[1]^q[1],2**r-1);
  for(const a of points)for(const b of points){ensure(D(a,b)===D(b,a),'distance symmetry');ensure((D(a,b)===0)===(a[0]===b[0]&&a[1]===b[1]),'endpoint separation');
   for(const c of points){ensure(D(a,c)<=D(a,b)+D(b,c),'triangle inequality');triples++;}}
 }
 for(let r=1;r<=6;r++){
  for(let L=0;L<=10;L++){
   let brute=0;for(let z=-L;z<=L;z++)for(let m=0;m<2**r;m++)if(distance(z,m,2**r-1)<=L)brute++;
   let formula=0;for(let w=0;w<=Math.min(r,L);w++)formula+=binomial(r,w)*(L+(mod(L-w,2)===0?1:0));
   ensure(brute===formula,'finite meeting ball volume');if(L>=r)ensure(brute===2**r*L+2**(r-1),'eventual exact linear volume');balls++;
  }
  for(let parity=0;parity<2;parity++){
   const ms=Array.from({length:2**r},(_,m)=>m).filter(m=>pop(m)%2===parity);ensure(ms.length===2**(r-1),'fixed address record cardinality');
   let max=0;for(const m of ms)for(const n of ms)max=Math.max(max,distance(0,m^n,2**r-1));
   ensure(max===2*Math.floor(r/2),'record fibre diameter');fibres++;
  }
 }
 return {metric_triangle_cases:triples,exact_ball_counts:balls,record_parity_fibres:fibres,maximum_record_width:6};
});

const drop=(m,omit)=>{let b=0,j=0;for(let i=0;i<Math.max(m.toString(2).length,omit+1);i++)if(i!==omit){b|=((m>>i)&1)<<j;j++;}return b;};
const recover=(b,x,n,omit)=>{ensure((n-x)%2===0,'compressed address parity');const last=mod((n-x)/2,2)^(pop(b)%2);
 let m=last<<omit,j=0;for(let i=0;i<Math.max(b.toString(2).length+1,omit+1);i++)if(i!==omit){m|=((b>>j)&1)<<i;j++;}return m;};
const compress=(a,n,omit)=>joint(jes(a).map(([x,m,v])=>{const b=drop(m,omit);ensure(recover(b,x,n,omit)===m,'source in native parity subspace');return [x,b,v];}));
const expand=(a,n,omit)=>joint(jes(a).map(([x,b,v])=>[x,recover(b,x,n,omit),v]));
const compressedStep=(a,i,omit)=>i===omit?systemStep(a):step(a,i<omit?i:i-1);
check('one_bit_compression_intertwines_arbitrary_parity_fields',()=>{
 let stages=0,basisLabels=0,steps=0,readoutBlocks=0;
 for(let r=1;r<=4;r++)for(let omit=0;omit<r;omit++){
  const rows=[];for(let x=-4;x<=4;x+=2)for(let m=0;m<2**r;m++)if(mod(-x/2-pop(m),2)===0){rows.push([x,m,col([m+1,x-m-1])]);basisLabels+=2;}
  const source=joint(rows),protocols=[Array.from({length:r},(_,i)=>i),Array.from({length:2*r+2},(_,i)=>r-1-i%r),[omit,omit,...Array.from({length:r},(_,i)=>i)]];
  for(const beta of protocols){let full=source,small=compress(source,0,omit);
   for(let n=0;n<=beta.length;n++){
    jeq(expand(small,n,omit),full,'compressed full field reconstruction');
    const g=moment(full,0),bar=moment(small,0),filtered=pairs(pes(bar).filter(([x,y])=>mod(x-y,4)===0));
    peq(filtered,g,'off-address equality filter');
    for(const x of new Set(pes(bar).flatMap(([x,y])=>[x,y]))){ensure(rho(g,x).eq(rho(bar,x))&&current(g,x).eq(current(bar,x)),'local compressed targets');readoutBlocks++;}
    stages++;if(n<beta.length){full=step(full,beta[n]);small=compressedStep(small,beta[n],omit);steps++;}
   }
  }
 }
 return {maximum_record_width:4,each_choice_of_omitted_slot_checked:true,arbitrary_parity_basis_labels:basisLabels,
         full_field_reconstructions:stages,compressed_native_steps:steps,local_readout_comparisons:readoutBlocks,
         off_address_readout_requires_mod_four_filter:true};
});

const all=Array.from({length:7},(_,n)=>partitions(n)),states=new Map([['',seed]]);
for(let n=1;n<=6;n++)for(const sigma of all[n])states.set(sigma.join(','),step(states.get(sigma.slice(0,-1).join(',')),sigma[n-1]));
const histories=[[{word:'',dirs:[],x:0,v:e0}]];
for(let n=1;n<=6;n++)histories.push(histories[n-1].flatMap(h=>[[H,'H'],[K,'K']].map(([M,tag])=>{
 const v=o.mul(M,h.v),b=o.isZero(o.mul(P,v))?1:0;return {word:h.word+tag,dirs:[...h.dirs,b],x:h.x+1-2*b,v};})));
check('literal_cut_histories_keep_tags_and_compress_exactly',()=>{
 let allocations=0,wordInstances=0,compressions=0;
 for(let n=0;n<=6;n++)for(const sigma of all[n]){
  const literal=new Map();for(const h of histories[n]){const m=rec(sigma,h.dirs);addAt(literal,h.x+','+m,h.v);
   ensure(mod(h.x+2*pop(m)-n,4)===0,'literal cut invariant');wordInstances++;}
  const full=states.get(sigma.join(','));jeq(literal,full,'literal H/K paths / native gates');
  ensure(total(moment(full,0,denom(n))).eq(1),'single seed normalized energy');
  for(let omit=0;omit<arity(sigma);omit++){let small=seed;for(const i of sigma)small=compressedStep(small,i,omit);
   jeq(small,compress(full,n,omit),'single seed compressed continuation');compressions++;}
  allocations++;
 }
 const sigma=[0,0,1,1],a=[1,1,0,0],b=[1,0,1,0];
 ensure(count(a)===count(b)&&rec(sigma,a)===0&&rec(sigma,b)===3,'hidden record distance witness');
 ensure(distance(0,3,3)===2,'same address positive meeting distance');
 return {allocation_partitions:allocations,maximum_depth:6,literal_HK_word_instances:wordInstances,
         compressed_protocol_comparisons:compressions,full_history_tags_erased:false};
});
check('first_reuse_has_exact_current_and_energy_onset',()=>{
 let allocations=0,freshStages=0,currentWitnesses=0,energyWitnesses=0;
 for(let n=0;n<=6;n++)for(const sigma of all[n]){
  let tau=Infinity;const seen=new Set();for(let k=0;k<n;k++){if(seen.has(sigma[k])){tau=k+1;break;}seen.add(sigma[k]);}
  for(let t=0;t<=Math.min(n,tau);t++){
   const g=moment(states.get(sigma.slice(0,t).join(',')),0);
   for(let x=-t-2;x<=t+2;x++){
    const countEnergy=(t-x)%2===0?binomial(t,(t-x)/2):0;ensure(rho(g,x).eq(countEnergy),'count energy through first reuse');
    if(t<tau)ensure(current(g,x).zero(),'no current before first reuse');
   }freshStages++;
  }
  if(tau<=n){const g=moment(states.get(sigma.slice(0,tau).join(',')),0);ensure(current(g,tau-2).eq(2),'first current exact normalized numerator');currentWitnesses++;}
  if(tau<n){const g=moment(states.get(sigma.slice(0,tau+1).join(',')),0);ensure(rho(g,tau-1).sub(tau+1).eq(2),'first energy excess numerator');energyWitnesses++;}
  allocations++;
 }
 return {allocation_partitions:allocations,count_energy_prefixes:freshStages,first_current_witnesses:currentWitnesses,
         next_event_energy_witnesses:energyWitnesses,normalization_denominator:'2^event',original_seed_and_blank_records_required:true};
});

const scalarAdd=(out,key,v)=>{const s=(out.get(key)||f(0)).add(v);if(s.zero())out.delete(key);else out.set(key,s);};
const pairAdvance=(g,i)=>{const out=new Map();for(const [key,c]of g){const [x,a,m,y,b,n]=key.split(',').map(Number);
 for(let u=0;u<2;u++)for(let v=0;v<2;v++)scalarAdd(out,[x+1-2*u,u,m^(u<<i),y+1-2*v,v,n^(v<<i)].join(','),
  c.mul(F[u][a].rad).mul(F[v][b].rad).mul('1/2'));}return out;};
const pairFromJoint=(a,scale)=>{const out=new Map();for(const [x,m,v]of jes(a))for(const [y,n,u]of jes(a))for(let b=0;b<2;b++)for(let c=0;c<2;c++)
 scalarAdd(out,[x,b,m,y,c,n].join(','),v[b][0].rad.mul(u[c][0].rad).mul(scale));return out;};
const scalarMapEq=(a,b,msg)=>{for(const k of new Set([...a.keys(),...b.keys()]))ensure((a.get(k)||f(0)).eq(b.get(k)||f(0)),msg+' '+k);};
check('signed_pair_evolution_and_cancellation_boundary',()=>{
 let comparisons=0;
 for(const source of [seed,joint([[0,0,col([1,1])]]),joint([[0,0,col([1,2])],[2,1,col([2,-1])]])])
 for(const beta of [[0,0],[0,1,0],[1,0,1,0]]){
  let a=source,g=pairFromJoint(a,1);for(let t=0;t<beta.length;t++){a=step(a,beta[t]);g=pairAdvance(g,beta[t]);
   scalarMapEq(g,pairFromJoint(a,denom(t+1)),'full signed pair / joint propagation');comparisons++;}
 }
 let channel=new Map([['0,0,0,0,1,0',f(1)]]);channel=pairAdvance(pairAdvance(channel,0),0);
 let j=f(0),nonzeroTerms=0;for(const [key,c]of channel){const [x,a,m,y,b,n]=key.split(',').map(Number);if(x===y&&m===n&&a!==b){j=j.add(c);nonzeroTerms++;}}
 ensure(j.zero()&&nonzeroTerms===2,'surviving opposite signed current terms cancel');
 ensure(channel.get('0,1,1,0,0,1').eq('-1/4')&&channel.get('0,0,1,0,1,1').eq('1/4'),'explicit cancellation coefficients');
 const balanced=moment(run([0,0],-1,joint([[0,0,col([1,1])]])),0);
 for(const x of [-2,0,2])ensure(current(balanced,x).zero(),'positive balanced preparation current cancellation');
 ensure(hasCurrent([0,0],0,0),'structural current possible despite zero response');
 return {full_pair_field_comparisons:comparisons,cancelling_nonzero_terms:nonzeroTerms,
         exact_cancellation_coefficients:['-1/4','1/4'],positive_preparation_zero_current_checked:true};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk);
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const tasks=[
 ['record_exchange_involution',mul(kk,kk),ii],
 ['role_parity_involution',mul(hh,hh),ii],
 ['native_branch_normalization',mul(ff,ff),sc(2,ii)],
 ['equal_role_positive_trace_channel',mul(pp,pp),pp],
 ['equal_role_negative_trace_channel',mul(qq,qq),qq],
 ['different_role_energy_channel_vanishes',mul(pp,qq),sc(0,ii)],
 ['current_exchanges_role_channels',mul(kk,pp,kk),qq],
 ['current_exchanges_negative_channel',mul(kk,qq,kk),pp],
 ['negative_persistent_branch_sign',mul(qq,ff,qq),sc(-1,qq)],
 ['positive_persistent_branch_sign',mul(pp,ff,pp),pp],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[2].result.certificate));altered.output.terms=[[[],['3','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered normalization proof accepted');
const oneSlot=states.get('0,0'),small=compress(oneSlot,2,0),fullReadout=moment(oneSlot,0),smallReadout=moment(small,0);
const negatives={
 address_gap_alone_is_meeting_distance:distance(0,3,3)===2,
 lower_length_bounds_alone_guarantee_return:Math.max(4,pop(3))<=4&&!admits([0,0,1,1],4,3),
 unmatched_slot_can_return_without_a_write:!admits([0,0],0,3),
 parity_mismatch_can_be_repaired:distance(0,1,1)===Infinity,
 odd_address_gap_can_meet:mod((1+1-(-1)),2)===1,
 all_catalogues_have_same_distance:distance(0,3,3)===2&&distance(0,3,1)===Infinity,
 current_and_energy_first_return_coincide:hasCurrent([0],1,1)&&!hasEnergy([0],1,1)&&hasEnergy([0,0],1,1),
 positive_current_return_occurs_after_one_event:!hasCurrent([0],0,0)&&hasCurrent([0,0],0,0),
 word_order_never_changes_prefix_latency:hasCurrent([0,0],0,0)&&!hasCurrent([0,1],0,0),
 omitted_bit_can_be_recovered_without_address:recover(0,0,2,0)!==recover(0,2,2,0),
 compressed_full_readout_needs_no_equality_filter:!samePairs(fullReadout,smallReadout),
 all_record_bits_are_independent_at_fixed_address:2**3!==2**(3-1),
 first_reuse_already_changes_energy:rho(moment(states.get('0,0'),0),0).eq(2)&&!current(moment(states.get('0,0'),0),0).zero(),
 every_geometric_return_has_nonzero_aggregate_current:hasCurrent([0,0],0,0)&&current(moment(run([0,0],-1,joint([[0,0,col([1,1])]])),0),0).zero(),
 full_history_metric_separates_every_source_word:rec([0,0,0],[0,1,0])===rec([0,0,0],[1,0,0])&&count([0,1,0])===count([1,0,0]),
 empty_catalogue_supplies_positive_events:distance(1,1,0)===Infinity&&distance(0,0,0)===0,
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r22.meeting-geometry.v1',status:'PASS_R22_NATIVE_MEETING_GEOMETRY',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,
 altered_native_certificate_rejected:alteredRejected,
 native_meeting_distance_derived:true,one_record_bit_reconstructed:true,
 ordinary_complex_measurement_stochastic_or_classical_metric_premise:false,
 physical_allocation_selected:false,physical_metric_c_alpha_derived:false,formal_proof_assistant_verified:false,
 scope:'Exact fixed-word meeting capacity, constructed complete-catalogue event distance, endpoint metric, ball counts, address-carried parity compression, first-reuse onset and signed cancellation. Written arbitrary-depth proofs with scoped finite checks; no physical metric selection.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
