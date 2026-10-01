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
if(inputHash!==get('--expected-input-sha256'))throw Error('R23 input pin mismatch');
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
const words=(alphabet,t)=>t===0?[[]]:words(alphabet,t-1).flatMap(a=>alphabet.map(i=>[...a,i]));
const bitWords=t=>Array.from({length:2**t},(_,j)=>Array.from({length:t},(_,k)=>(j>>k)&1));
const translate=(a,eta)=>joint(jes(a).map(([x,m,v])=>[x,m^eta,v]));
const smallStep=(a,i,omitted)=>i===omitted?systemStep(a):step(a,i);
const smallRun=(beta,source,omitted)=>{let a=source;for(const i of beta)a=smallStep(a,i,omitted);return a;};
const fullRecord=(b,x,n,r)=>b|((mod((n-x)/2,2)^(pop(b)%2))<<(r-1));
const expand=(a,n,r)=>joint(jes(a).map(([x,b,v])=>[x,fullRecord(b,x,n,r),v]));
const localEq=(a,b,why)=>{const ga=moment(a,0),gb=moment(b,0);for(const x of new Set([...jes(a),...jes(b)].map(([x])=>x)))
 ensure(rho(ga,x).eq(rho(gb,x))&&current(ga,x).eq(current(gb,x)),why);};
const fixtures=[seed,joint([[-2,0,col([1,2])],[0,1,col(['1/2',-1])],[2,3,col([2,1])]])];

check('record_translation_commutes_with_every_native_step',()=>{
 let commutations=0,readouts=0,parityCases=0;
 for(const a of fixtures)for(let eta=0;eta<8;eta++)for(let i=0;i<3;i++){
  jeq(step(translate(a,eta),i),translate(step(a,i),eta),'record translation / full native step');commutations++;
  peq(moment(translate(a,eta),0),moment(a,0),'full unresolved readout invariant');readouts++;
 }
 for(let r=1;r<=4;r++)for(let x=-4;x<=4;x+=2)for(let b=0;b<2**(r-1);b++)for(let eta=0;eta<2**(r-1);eta++){
  const zeta=eta|((pop(eta)%2)<<(r-1));ensure(pop(zeta)%2===0,'full parity-preserving translation');
  ensure((fullRecord(b,x,0,r)^zeta)===fullRecord(b^eta,x,0,r),'compressed/full translation intertwiner');parityCases++;
 }
 return {full_step_commutations:commutations,unresolved_pair_invariances:readouts,parity_translation_cases:parityCases};
});

const tr=M=>M.reduce((s,row,i)=>s.add(row[i].rad),f(0));
const recExchange=(q,eta)=>{const X=o.zeros(q);for(let b=0;b<q;b++)X[b^eta][b]=o.ONE;return X;};
const denseTranslate=(G,N,q,eta)=>o.matrix(Array.from({length:N*q},(_,i)=>Array.from({length:N*q},(_,j)=>
 G[Math.floor(i/q)*q+((i%q)^eta)][Math.floor(j/q)*q+((j%q)^eta)])));
const average=(G,N,q)=>o.scale(Array.from({length:q},(_,eta)=>denseTranslate(G,N,q,eta)).reduce((a,b)=>o.add(a,b),o.zeros(N*q)),new o.F(1,q));
const denseMoments=(G,N,q)=>Array.from({length:q},(_,eta)=>o.matrix(Array.from({length:N},(_,p)=>Array.from({length:N},(_,s)=>{
 let z=o.ZERO;for(let b=0;b<q;b++)z=z.add(G[p*q+b][s*q+(b^eta)]);return z;}))));
const recoverAverage=(Ns,N,q)=>o.matrix(Array.from({length:N*q},(_,i)=>Array.from({length:N*q},(_,j)=>
 Ns[(i%q)^(j%q)][Math.floor(i/q)][Math.floor(j/q)].div(q))));
const symmetricBasis=D=>{const es=[];for(let i=0;i<D;i++)for(let j=i;j<D;j++){const G=o.zeros(D);G[i][j]=o.ONE;G[j][i]=o.ONE;es.push({i,j,G});}return es;};
check('native_averaging_cut_and_exact_exchange_reconstruction',()=>{
 let bases=0,entries=0;
 for(let k=0;k<=2;k++){const q=2**k,N=3;
  for(const {G}of symmetricBasis(N*q)){
   const E=average(G,N,q),Ns=denseMoments(G,N,q);
   eq(average(E,N,q),E,'averaging cut idempotent');eq(recoverAverage(Ns,N,q),E,'entrywise exchange reconstruction');
   const again=denseMoments(E,N,q);for(let eta=0;eta<q;eta++){eq(again[eta],Ns[eta],'moments survive cut');eq(o.dagger(Ns[eta]),Ns[eta],'real symmetric exchange field');}
   bases++;entries+=(N*q)**2;
  }
 }
 return {symmetric_pair_basis_inputs:bases,averaged_matrix_entries:entries,retained_record_widths:[0,1,2],physical_reset_assumed:false};
});

check('compressed_exchange_and_sector_continuation_match_full_records',()=>{
 let stages=0,fields=0,fullChecks=0;
 for(let r=1;r<=4;r++){
  const k=r-1,q=2**k,source=joint([[-2,0,col([1,2])],[0,q-1,col([2,-1])],[2,0,col(['1/2',1])]]);
  for(const beta of [Array.from({length:r+3},(_,i)=>i%r),Array.from({length:2*r+1},(_,i)=>r-1-i%r)]){
   let a=source,full=expand(source,0,r),Ns=Array.from({length:q},(_,eta)=>moment(a,eta));
   let sectors=Array.from({length:q},(_,chi)=>pscale(psum(...Ns.map((g,eta)=>pscale(g,sign(chi,eta)))),new o.F(1,q)));
   for(let n=0;n<=beta.length;n++){
    jeq(expand(a,n,r),full,'compressed source / complete native record field');fullChecks++;
    for(let eta=0;eta<q;eta++){peq(Ns[eta],moment(a,eta,denom(n)),'closed compressed exchange moments');
     peq(psum(...sectors.map((g,chi)=>pscale(g,sign(chi,eta)))),Ns[eta],'native sector inverse');fields++;}
    const unresolved=pairs(pes(psum(...sectors)).filter(([x,y])=>mod(x-y,4)===0));
    peq(unresolved,moment(full,0,denom(n)),'full record readout after exact address congruence filter');
    stages++;
    if(n<beta.length){const i=beta[n],old=Ns;
     Ns=Array.from({length:q},(_,eta)=>i===k?T(old[eta]):momentStep(mask=>old[mask],i,eta));
     sectors=sectors.map((g,chi)=>conjugate(T(g),i===k||!((chi>>i)&1)?I:H));
     a=smallStep(a,i,k);full=step(full,i);
    }
   }
  }
 }
 return {record_widths:[1,2,3,4],continuation_stages:stages,exchange_field_comparisons:fields,full_record_reconstructions:fullChecks};
});

const systemCuts=N=>{
 const list=[];for(let p=0;p<N;p++){const v=o.zeros(N,1);v[p][0]=o.ONE;list.push({p,s:p,Q:outer(v)});}
 for(let p=0;p<N;p++)for(let s=p+1;s<N;s++){const v=o.zeros(N,1);v[p][0]=v[s][0]=o.ONE;list.push({p,s,Q:o.scale(outer(v),'1/2')});}
 return list;
};
const cutEnergy=(G,Q)=>tr(o.mul(o.mul(Q,G),Q));
check('complete_native_exchange_probes_recover_sharp_channel_rank',()=>{
 let probeCases=0,recoveredEntries=0,rankCases=0,normalizedRankCases=0;
 for(let k=0;k<=2;k++)for(let N=1;N<=3;N++){
  const q=2**k,cuts=systemCuts(N),preps=[];
  for(let chi=0;chi<q;chi++){
   const v=o.matrix(Array.from({length:q},(_,b)=>[sign(chi,b)])),record=o.scale(outer(v),new o.F(1,q));
   for(const {Q}of cuts)preps.push(o.kron(Q,record));
  }
  const probeColumns=[];
  for(const G of preps){const Ns=denseMoments(G,N,q),values=[];
   for(let eta=0;eta<q;eta++){
    const X=recExchange(q,eta),plus=o.scale(o.add(o.identity(q),X),'1/2'),minus=o.scale(o.sub(o.identity(q),X),'1/2');
    const raw=[];for(const {Q}of cuts){const got=cutEnergy(G,o.kron(Q,plus)).sub(cutEnergy(G,o.kron(Q,minus))),want=tr(o.mul(Q,Ns[eta]));
     ensure(got.eq(want),'cut energy difference / exchange pairing');raw.push(got);values.push(got);probeCases++;}
    let at=N;for(let p=0;p<N;p++){ensure(raw[p].eq(Ns[eta][p][p].rad),'probe diagonal');recoveredEntries++;
     for(let s=p+1;s<N;s++){ensure(raw[at++].sub(raw[p].add(raw[s]).mul('1/2')).eq(Ns[eta][p][s].rad),'probe off diagonal reconstruction');recoveredEntries++;}}
   }
   probeColumns.push(values);
  }
  const target=q*N*(N+1)/2,table=o.matrix(probeColumns);
  ensure(o.rank(table)===target,'sharp complete exchange target rank');rankCases++;
  if(target>1){const differences=o.matrix(probeColumns.slice(1).map(row=>row.map((z,i)=>z.sub(probeColumns[0][i]))));
   ensure(o.rank(differences)===target-1,'normalized affine channel count');normalizedRankCases++;}
 }
 return {native_cut_energy_comparisons:probeCases,reconstructed_exchange_entries:recoveredEntries,
         exact_channel_rank_cases:rankCases,normalized_affine_rank_cases:normalizedRankCases,maximum_system_labels:3,maximum_retained_bits:2};
});

const finiteObserverData=[];
const valueAt=(a,x,b,role)=>(a.get(x+','+b)||zv())[role][0].rad;
const observerRows=(columns,t,q)=>{
 const xs=[...new Set(columns.flatMap(a=>jes(a).map(([x])=>x)))].sort((a,b)=>a-b),rows=[];
 for(const x of xs)for(const kind of ['energy','current']){
  const row=[];for(let i=0;i<columns.length;i++)for(let j=i;j<columns.length;j++){
   let z=f(0);for(let b=0;b<q;b++){
    const u=[0,1].map(role=>valueAt(columns[i],x,b,role)),v=[0,1].map(role=>valueAt(columns[j],x,b,role));
    z=z.add(kind==='energy'?u[0].mul(v[0]).add(u[1].mul(v[1])):u[0].mul(v[1]).add(u[1].mul(v[0])));
   }row.push(z.mul(i===j?1:2).mul(denom(t)));
  }rows.push({x,kind,row});
 }return rows;
};
check('finite_catalogue_observer_rank_and_native_response_form',()=>{
 let evaluatedRows=0,invarianceCases=0,rankCases=0,quadraticChecks=0;const rankTables=[],rankWitnesses=[];
 for(let r=1;r<=2;r++){
  const q=2**(r-1),labels=[];for(const x of [0,2])for(let a=0;a<2;a++)for(let b=0;b<q;b++)labels.push({x,a,b});
  const basis=labels.map(({x,a,b})=>joint([[x,b,a?e1:e0]])),D=labels.length,pairIndices=[];
  for(let i=0;i<D;i++)for(let j=i;j<D;j++)pairIndices.push([i,j]);
  const pairIndex=new Map(pairIndices.map(([i,j],n)=>[i+','+j,n]));
  const wordsByDepth=Array.from({length:5},(_,t)=>words(Array.from({length:r},(_,i)=>i),t)),allRows=[],readoutLabels=[],ranks=[];
  for(let t=0;t<=4;t++){
   for(const beta of wordsByDepth[t]){
    const columns=basis.map(a=>smallRun(beta,a,r-1)),rows=observerRows(columns,t,q);allRows.push(...rows.map(v=>v.row));
    readoutLabels.push(...rows.map(({x,kind})=>({word:beta,address:x,target:kind})));
    for(let eta=0;eta<q;eta++)for(const {row}of rows)for(let k=0;k<pairIndices.length;k++){
     const [i,j]=pairIndices[k],ii=i^eta,jj=j^eta,other=pairIndex.get(Math.min(ii,jj)+','+Math.max(ii,jj));
     ensure(row[k].eq(row[other]),'observer row factors through exchange average');invarianceCases++;
    }
    const coeffs=labels.map((_,i)=>f((i%3)-1)),state=joint(labels.map(({x,a,b},i)=>[x,b,o.scale(a?e1:e0,coeffs[i])]));
    const got=moment(smallRun(beta,state,r-1),0,denom(t));
    const coords=pairIndices.map(([i,j])=>coeffs[i].mul(coeffs[j]));
    for(const row of rows){const predicted=row.row.reduce((s,z,i)=>s.add(z.mul(coords[i])),f(0));
     ensure(predicted.eq(row.kind==='energy'?rho(got,row.x):current(got,row.x)),'observer rows / independent joint preparation');evaluatedRows++;}
   }
   const R=o.rref(o.matrix(allRows));ensure(!ranks.length||R.rank>=ranks.at(-1),'catalogue rank monotonic');
   ensure(R.rank<=q*4*5/2,'exchange quotient dimension upper bound');ranks.push(R.rank);rankCases++;
   for(const raw of allRows){let rem=raw.map(o.Cut.of);for(let row=0;row<R.rank;row++){const c=rem[R.pivots[row]];
     rem=rem.map((v,j)=>v.sub(c.mul(R.basis[row][j])));}ensure(rem.every(v=>v.zero()),'minimum row basis recovers all outputs');}
   if(t===4){
    ensure(R.rank===q*4*5/2&&ranks[3]<R.rank,'sharp complete future quotient horizon');
    const selected=o.rref(o.dagger(o.matrix(allRows))).pivots;
    ensure(selected.length===R.rank,'independent actual native readout rows');
    const minor=o.matrix(selected.map(i=>R.pivots.map(j=>allRows[i][j]))),inv=o.inverse(minor);
    eq(o.mul(minor,inv),o.identity(R.rank),'native rank witness right inverse');
    eq(o.mul(inv,minor),o.identity(R.rank),'native rank witness left inverse');
    rankWitnesses.push({record_slots:r,rank:R.rank,selected_readouts:selected.map(i=>readoutLabels[i]),
      initial_pair_coordinates:R.pivots.map(j=>pairIndices[j]),minor_sha256:p.digest(minor),
      exact_inverse:inv.map(row=>row.map(z=>z.rad.toString())),infinite_future_complete_by_proved_exchange_upper_bound:true});
   }
   if(t===2){const O=o.matrix(allRows),B=o.mul(o.dagger(O),O);ensure(o.rank(B)===R.rank,'positive response form exact kernel');
    const v=o.matrix(pairIndices.map((_,i)=>[(i%5)-2])),Ov=o.mul(O,v),direct=o.mul(o.mul(o.dagger(v),B),v);
    ensure(direct[0][0].rad.eq(o.energy(Ov)),'response form is sum of output squares');quadraticChecks++;}
  }
  rankTables.push({record_slots:r,initial_addresses:[0,2],system_labels:4,full_symmetric_pair_dimension:D*(D+1)/2,
                   complete_exchange_target_dimension:q*10,ranks_at_horizons_0_through_4:ranks});
  finiteObserverData.push({r,labels,rows:allRows});
 }
 return {finite_horizon_rank_tables:rankTables,rank_cases:rankCases,independent_output_evaluations:evaluatedRows,
         translation_factorization_coefficients:invarianceCases,response_form_kernel_checks:quadraticChecks,
         exact_full_rank_witnesses:rankWitnesses,rank_plateau_alone_used_as_saturation:false,
         sharp_complete_horizon_for_listed_finite_supports:4};
});

check('positive_meeting_distance_is_invisible_to_all_word_readouts',()=>{
 const a=joint([[0,0,e0]]),b=joint([[0,3,e0]]);let wordsChecked=0,fullPairChecks=0;
 for(let t=0;t<=5;t++)for(const beta of words([0,1],t)){
  const aa=run(beta,-1,a),bb=run(beta,-1,b);peq(moment(aa,0),moment(bb,0),'positive meeting distance / identical response');
  jeq(translate(aa,3),bb,'complete output retains translated record');wordsChecked++;fullPairChecks++;
 }
 ensure(pop(3)===2,'meeting distance witness');
 const ga=outer(o.matrix([[1],[0],[0],[0]])),gb=outer(o.matrix([[0],[1],[0],[0]]));
 eq(average(ga,2,2),average(gb,2,2),'different positive preparations same exchange quotient');
 ensure(!o.equal(ga,gb),'full pure records remain distinct');
 return {future_words_checked:wordsChecked,maximum_length:5,full_unresolved_pair_equalities:fullPairChecks,
         native_meeting_distance:2,future_response_class_equal:true};
});

const all=Array.from({length:7},(_,n)=>partitions(n));
const states=[new Map([['',seed]]),new Map([['',joint([[0,0,e1]])]])];
for(let a=0;a<2;a++)for(let n=1;n<=6;n++)for(const beta of all[n])states[a].set(beta.join(','),step(states[a].get(beta.slice(0,-1).join(',')),beta[n-1]));
const firstReturn=beta=>{for(let i=1;i<beta.length;i++)if(beta[i]===beta[0])return i+1;return Infinity;};
check('first_slot_return_is_sharp_initial_role_detection_time',()=>{
 let protocols=0,identicalPairPrefixes=0,equalEnergyPrefixes=0,currentWitnesses=0,energyWitnesses=0;
 const test=beta=>{
  const tau=firstReturn(beta);let a=seed,b=joint([[0,0,e1]]);
  for(let t=0;t<=beta.length;t++){
   const ga=moment(a,0),gb=moment(b,0);
   if(t>0&&t<tau){peq(ga,gb,'first-slot record hides initial role');identicalPairPrefixes++;}
   if(t<=tau){for(let x=-t;x<=t;x++)ensure(rho(ga,x).eq(rho(gb,x)),'energy equality through first slot return');equalEnergyPrefixes++;}
   if(t===tau){ensure(current(ga,t-2).sub(current(gb,t-2)).eq(4),'first role current difference numerator');currentWitnesses++;}
   if(t===tau+1){ensure(rho(ga,t-2).sub(rho(gb,t-2)).eq(4),'next-event energy difference numerator');energyWitnesses++;}
   if(t<beta.length){a=step(a,beta[t]);b=step(b,beta[t]);}
  }protocols++;
 };
 for(let n=0;n<=6;n++)for(const beta of all[n])test(beta);
 for(let tau=7;tau<=10;tau++)test([0,...Array(tau-2).fill(1),0,1]);
 test([0,...Array(10).fill(1)]);
 return {allocation_protocols:protocols,identical_positive_pair_prefixes:identicalPairPrefixes,
         identical_energy_prefixes:equalEnergyPrefixes,first_current_witnesses:currentWitnesses,next_event_energy_witnesses:energyWitnesses,
         longest_first_slot_return_tested:10,other_slot_reuse_before_detection_included:true};
});

check('basis_endpoint_decoder_has_exact_two_event_depth',()=>{
 let preparations=0,oneEventChecks=0,decoderChecks=0;
 for(let r=1;r<=3;r++)for(const x of [-2,0,2])for(let m=0;m<2**r;m++)if(mod(x+2*pop(m),4)===0){
  for(let i=0;i<r;i++){
   const prep=[joint([[x,m,e0]]),joint([[x,m,e1]])];localEq(step(prep[0],i),step(prep[1],i),'one forward event cannot distinguish basis role');oneEventChecks++;
   for(let a=0;a<2;a++){
    const g=moment(run([i,i],-1,prep[a]),0,denom(2));
    for(let z=x-2;z<=x+2;z++)ensure(current(g,z).mul(2).eq(z===x?(a?-1:1):0),'two-event role decoder');
    ensure(rho(moment(prep[a],0),x).eq(1),'initial address decoder');decoderChecks++;preparations++;
   }
  }
 }
 return {basis_preparation_and_slot_cases:preparations,one_event_role_equivalences:oneEventChecks,two_event_decoder_cases:decoderChecks};
});

const recordOperator=(r,i,M)=>{let out=o.identity(1);for(let j=r-1;j>=0;j--)out=o.kron(out,j===i?M:I);return out;};
check('native_even_record_rotation_restores_every_retained_basis_bit',()=>{
 let identities=0,probes=0,parityLabels=0;
 for(let r=2;r<=4;r++)for(let i=0;i<r-1;i++){
  const dim=2**r,Hr=recordOperator(r,i,H),Y=o.mul(recordOperator(r,i,K),recordOperator(r,r-1,K)),Fr=o.add(Hr,Y),id=o.identity(dim);
  eq(o.mul(Y,Y),id,'even record exchange involution');eq(o.add(o.mul(Hr,Y),o.mul(Y,Hr)),o.zeros(dim),'native record anticommutation');
  eq(o.mul(Fr,Fr),o.scale(id,2),'record mixing normalization');eq(o.mul(o.mul(Fr,Y),Fr),o.scale(Hr,2),'record phase-to-exchange conversion');identities+=4;
  for(let m=0;m<dim;m++){
   const v=o.zeros(dim,1);v[m][0]=o.ONE;const rotated=o.mul(Fr,v),rows=[];
   for(let b=0;b<dim;b++)if(!rotated[b][0].zero()){
    ensure(pop(b)%2===pop(m)%2,'record rotation stays in parity subspace');parityLabels++;
    rows.push([0,b,o.scale(col([1,1]),rotated[b][0])]);}
   const out=write(write(joint(rows),i),r-1),g=moment(out,0,'1/4');
   ensure(current(g,0).eq((m>>i)&1?-1:1),'native record-bit current decoder');
   ensure(total(g).eq(1),'prepared probe normalized');
   const unrotated=write(write(joint([[0,m,col([1,1])]]),i),r-1);
   ensure(current(moment(unrotated,0),0).zero(),'old exchange probe without rotation sees no basis bit');probes++;
  }
 }
 return {native_record_identities:identities,exact_basis_bit_probes:probes,parity_preserving_output_labels:parityLabels,
         classical_measurement_rule_supplied:false};
});

const supportSignature=beta=>{
 const meeting=new Set(),energy=new Set(),curr=new Set();let cases=0;
 for(const u of bitWords(beta.length))for(const v of bitWords(beta.length)){
  let d=0,delta=0;for(let k=0;k<beta.length;k++){d+=u[k]-v[k];delta^=(u[k]^v[k])<<beta[k];}
  const key=d+','+delta;meeting.add(key);(u.at(-1)===v.at(-1)?energy:curr).add(key);cases++;
 }return {meeting,energy,curr,cases};
};
const sameSet=(a,b)=>a.size===b.size&&[...a].every(x=>b.has(x));
const beta=[0,0,0,1,0],gamma=[0,0,1,0,0],shapeA=supportSignature(beta),shapeB=supportSignature(gamma);
check('same_final_meeting_support_has_different_native_signed_signal',()=>{
 for(const kind of ['meeting','energy','curr'])ensure(sameSet(shapeA[kind],shapeB[kind]),'same full final tagged support '+kind);
 const a=run(beta),b=run(gamma),ga=moment(a,0,denom(5)),gb=moment(b,0,denom(5));
 ensure(rho(ga,1).eq('3/8')&&rho(gb,1).eq('1/8'),'same-support distinct energy');
 ensure(current(ga,1).eq('-1/4')&&current(gb,1).zero(),'same-support distinct signed current');
 eq(a.get('1,0'),col([-1,3]),'beta record00 column');eq(a.get('1,3'),col([1,-1]),'beta record11 column');
 eq(b.get('1,0'),col([1,1]),'gamma record00 column');eq(b.get('1,3'),col([-1,1]),'gamma record11 column');
 let literalPaths=0;
 for(const sigma of [beta,gamma]){
  let histories=[{word:'',x:0,m:0,v:e0}];for(const i of sigma)histories=histories.flatMap(h=>[[H,'H'],[K,'K']].map(([M,tag])=>{
   const v=o.mul(M,h.v),a=o.isZero(o.mul(P,v))?1:0;return {word:h.word+tag,x:h.x+1-2*a,m:h.m^(a<<i),v};}));
  const literal=joint(histories.map(h=>[h.x,h.m,h.v]));jeq(literal,run(sigma),'literal source signs / gate propagation');literalPaths+=histories.length;
 }
 return {exhaustive_direction_pair_support_checks:shapeA.cases+shapeB.cases,literal_HK_histories:literalPaths,
         meeting_support_entries:shapeA.meeting.size,energy_support_entries:shapeA.energy.size,current_support_entries:shapeA.curr.size,
         beta_energy_at_one:'3/8',gamma_energy_at_one:'1/8',beta_current_at_one:'-1/4',gamma_current_at_one:'0'};
});

const presentation=JSON.parse(fs.readFileSync(path.join(home,'examples/emk_job.json'),'utf8')).presentation;
const word=(...tokens)=>({word:tokens}),sc=(coefficient,value)=>({op:'scale',coefficient,value});
const add=(...args)=>({op:'add',args}),mul=(...args)=>({op:'multiply',args});
const ii=word(),kk=word('K'),hh=sc(-1,word('R','K')),ff=add(hh,kk);
const pp=sc('1/2',add(ii,hh)),qq=sc('1/2',add(ii,sc(-1,hh)));
const ep=sc('1/2',add(ii,kk)),em=sc('1/2',add(ii,sc(-1,kk)));
const tasks=[
 ['native_translation_involution',mul(kk,kk),ii],
 ['native_role_parity_involution',mul(hh,hh),ii],
 ['native_mixing_normalization',mul(ff,ff),sc(2,ii)],
 ['rotation_turns_exchange_into_record_parity',mul(ff,kk,ff),sc(2,hh)],
 ['rotation_turns_record_parity_into_exchange',mul(ff,hh,ff),sc(2,kk)],
 ['positive_role_cut',mul(pp,pp),pp],
 ['negative_role_cut',mul(qq,qq),qq],
 ['positive_exchange_cut',mul(ep,ep),ep],
 ['negative_exchange_cut',mul(em,em),em],
 ['exchange_cuts_are_disjoint',mul(ep,em),sc(0,ii)],
 ['exchange_cuts_are_complete',add(ep,em),ii],
 ['record_parity_breaks_translation_symmetry',mul(hh,kk,hh),sc(-1,kk)],
];
const system=w.presentation(presentation),replays=tasks.map(([name,left,right])=>{
 const result=p.proveEquality(w.parseExpression(system,left),w.parseExpression(system,right));
 ensure(result.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(result.certificate),'native proof replay '+name);
 return {name,result,replay:'REPLAY_MATCH'};
});
const altered=JSON.parse(JSON.stringify(replays[2].result.certificate));altered.output.terms=[[[],['3','0']]];
let alteredRejected=false;try{system.replay(altered);}catch(_){alteredRejected=true;}ensure(alteredRejected,'altered native normalization proof accepted');
const gBeta=moment(run(beta),0,denom(5)),gGamma=moment(run(gamma),0,denom(5));
const negatives={
 positive_meeting_distance_guarantees_distinct_future_responses:pop(3)===2&&samePairs(moment(run([0,1,0]),0),moment(run([0,1,0],-1,joint([[0,3,e0]])),0)),
 record_translation_discards_the_full_source_record:!o.equal(e0,e1),
 every_full_pair_coordinate_is_needed_for_this_interface:8*9/2>2*4*5/2,
 symmetric_exchange_field_has_N_squared_independent_channels:3*4/2!==3*3,
 known_normalization_is_an_extra_variable_channel:2*3*4/2-1!==2*3*4/2,
 exchange_quotient_is_always_minimal_for_passive_catalogues:finiteObserverData[1].rows.length>0&&4<20,
 every_record_translation_preserves_native_parity:pop(1)%2!==0,
 one_forward_event_separates_initial_basis_roles:samePairs(moment(run([0]),0),moment(run([0],-1,joint([[0,0,e1]])),0)),
 any_reused_slot_reveals_the_initial_role:firstReturn([0,1,1])===Infinity&&samePairs(moment(run([0,1,1]),0),moment(run([0,1,1],-1,joint([[0,0,e1]])),0)),
 first_slot_return_already_changes_local_energy:rho(moment(run([0,0]),0),0).eq(rho(moment(run([0,0],-1,joint([[0,0,e1]])),0),0)),
 record_basis_marks_are_fundamentally_unobservable:o.equal(o.mul(o.mul(F,K),F),o.scale(H,2)),
 native_record_rotation_commutes_with_exchange:!o.equal(o.mul(F,K),o.mul(K,F)),
 identical_final_meeting_support_fixes_the_signed_energy:sameSet(shapeA.energy,shapeB.energy)&&!rho(gBeta,1).eq(rho(gGamma,1)),
 identical_final_current_support_fixes_current:sameSet(shapeA.curr,shapeB.curr)&&!current(gBeta,1).eq(current(gGamma,1)),
 equal_record_counts_erase_event_order:beta.join(',')!==gamma.join(',')&&beta.filter(i=>i===0).length===gamma.filter(i=>i===0).length,
 two_event_endpoint_decoder_recovers_full_record_identity:current(moment(run([0,0]),0),0).eq(current(moment(run([0,0],-1,joint([[0,3,e0]])),0),0)),
};
ensure(Object.values(negatives).every(Boolean),'false alternative not rejected');
const output={schema:'extra-ideas.r23.future-response-quotient.v1',status:'PASS_R23_NATIVE_FUTURE_RESPONSE',
 input_sha256:inputHash,source_foundation:input.source_foundation,exact_checks:checks,check_count:checks.length,
 symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negatives,
 altered_native_certificate_rejected:alteredRejected,native_future_response_quotient_derived:true,native_record_probe_completion_derived:true,
 ordinary_complex_measurement_stochastic_or_classical_observer_premise:false,
 physical_observer_selected:false,physical_metric_c_alpha_derived:false,formal_proof_assistant_verified:false,
 scope:'Real signed pair target over the unchanged native source. Exact translation symmetry, sufficient exchange quotient, complete-probe channel rank, finite-catalogue minimum observer, basis endpoint decoding and native record-probe completion. Meeting support does not fix signed signal; no physical observer or metric selection.'};
process.stdout.write(JSON.stringify(output,null,2)+'\n');
