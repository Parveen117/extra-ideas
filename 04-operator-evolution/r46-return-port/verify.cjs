'use strict';
// R46 application verifier. Arithmetic, matrix elimination and the native
// return solver come from the separately pinned canonical RKF checkout.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const assert=require('node:assert/strict');
const here=__dirname,root=path.resolve(here,'../..');
const argv=process.argv.slice(2),options={};
for(let j=0;j<argv.length;j++){
  const k=argv[j];
  if(['--write','--check'].includes(k)) options[k]=true;
  else if(['--rkf-root','--publications-root','--materials-root','--thermo-root','--expected-protocol-sha256'].includes(k)){
    assert(argv[j+1]&&!argv[j+1].startsWith('--'),'Missing option value');options[k]=argv[++j];
  } else throw Error('Unknown option: '+k);
}
assert(!(options['--write']&&options['--check']),'Choose write or read-only check');
const homes={local:root};
for(const k of ['rkf','publications','materials','thermo']){
  assert(options['--'+k+'-root'],'Supply --'+k+'-root');homes[k]=path.resolve(options['--'+k+'-root']);
}
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const blob=b=>crypto.createHash('sha1').update(Buffer.concat([Buffer.from('blob '+b.length+'\0'),b])).digest('hex');
const bytes=x=>Buffer.from(JSON.stringify(x,null,2)+'\n');
const load=p=>JSON.parse(fs.readFileSync(p));
const pins=load(path.join(here,'SOURCE_PINS.json'));
function verifyPins(spec){
  assert.equal(spec.schema,'extra-ideas.r46.source-pins.v1');
  for(const [key,group] of Object.entries(spec.sources)){
    assert(homes[key],'Unknown source home');
    for(const row of group.files){
      const p=path.resolve(homes[key],row.path);
      assert(p.startsWith(homes[key]+path.sep),'Source outside checkout');
      const b=fs.readFileSync(p);
      assert.equal(sha(b),row.sha256,'Source SHA-256 mismatch: '+row.path);
      assert.equal(blob(b),row.git_blob_sha,'Source Git blob mismatch: '+row.path);
    }
  }
}
verifyPins(pins);
const protoPath=path.join(here,'PROTOCOL.json'),protocol=load(protoPath);
if(options['--expected-protocol-sha256'])assert.equal(sha(fs.readFileSync(protoPath)),options['--expected-protocol-sha256'],'Protocol pin mismatch');
const o=require(path.join(homes.rkf,'operator_foundation/core/native_operator.cjs'));
const r2=require(path.join(homes.rkf,'research/recognition_return/r2/variable_return.cjs'));
const r1=require(path.join(homes.rkf,'research/recognition_return/return_solver.cjs'));
const {F,Cut}=o,Q=F.of,Z=()=>new F(0),O=()=>new F(1);
const positive=(x,label='value')=>{x=Q(x);assert(x.n>0n,label+' must be positive');return x;};
const nonnegative=x=>{x=Q(x);assert(x.n>=0n,'Nonnegative value required');return x;};
const eq=(a,b)=>assert(Q(a).eq(b),'Exact equality failed: '+a+' != '+b);
const min=(a,b)=>Q(a).le(b)?Q(a):Q(b),max=(a,b)=>Q(a).le(b)?Q(b):Q(a);
function validateProtocol(p){
  assert.equal(p.data_kind,'DECLARED_IDEAL_COMPONENT_FIXTURE_NOT_MEASUREMENTS');
  positive(p.b);const d=positive(p.delta);assert(d.le(1)&&!d.eq(1));
  assert(Number.isSafeInteger(p.sections)&&p.sections>=1&&p.sections<=64);
  positive(p.reference_resistance_ohm);positive(p.source_voltage_volt);
  const e=nonnegative(p.relative_component_tolerance);assert(e.le(1)&&!e.eq(1));
  assert.equal(p.termination,'open_after_last_shunt');
  assert.equal(p.experiment_performed,false);assert.equal(p.fundamental_constant_selected,false);
  assert.deepEqual(p.controls.map(x=>x.id),['reference','minus','plus','third']);
  for(const [row,z] of p.controls.map((x,i)=>[x,[O(),O().sub(d),O().add(d),O().add(d.div(2))][i]]))eq(row.z,z);
}
validateProtocol(protocol);
const ledgerPath=path.join(here,'DERIVATION_LEDGER.json'),ledger=load(ledgerPath);
const proofText=fs.readFileSync(path.join(root,'02-relational-response/NATIVE_RETURN_PORT_R46.md'),'utf8');
function verifyLedger(l){
  assert.equal(l.native_engine_modified,false);
  assert.equal(l.physical_interface.status,'ASSUMED_INTERFACE_NOT_NATIVE_DERIVED');
  assert.equal(l.physical_interface.experimental_validation,false);
  assert.equal(l.evidence_scope,'WRITTEN_APPLICATION_PROOFS_AND_EXACT_FINITE_CHECKS_NOT_FORMAL_VERIFICATION');
  const sections=[...proofText.matchAll(/### R46\.(\d+) — ([^\n]+)\n([\s\S]*?)(?=\n#{2,3} |$)/g)];
  assert.equal(sections.length,5);assert.equal(l.claims.length,5);
  sections.forEach((s,j)=>{
    const c=l.claims[j];assert.equal(c.id,'R46.'+(j+1));assert.equal(c.title,s[2]);
    assert(s[3].includes('**Proof.**'));assert.equal(c.written_section_sha256,sha(Buffer.from(s[0])));
    assert(c.dependencies.includes('RKF_R2'));assert.equal(c.status,'PROVED_ON_DECLARED_TARGET');
  });
}
verifyLedger(ledger);
const checks=[],negative={},replays=[],counts={nodal_solves:0,node_balance_checks:0,power_balance_checks:0,component_corners:0,tail_controls:0};
function check(name,fn){fn();checks.push({name,passed:true});}
function rejects(name,fn){assert.throws(fn);negative[name]=true;}
function fixture(b,z,M,Rstar=1){
  b=positive(b);z=positive(z);Rstar=positive(Rstar);
  assert(Number.isSafeInteger(M)&&M>=1&&M<=64);
  const A=b.add(1).mul(z),B=b.mul(z),Rs=Rstar.div(A),Rp=Rstar.mul(B);
  return {s:Array(M).fill(Rs.div(Rstar)),g:Array(M).fill(Rstar.div(Rp)),Rs,Rp,cells:Array.from({length:2*M},(_,j)=>j%2?B:A)};
}
function reduceTarget(s,g,tail=0){
  assert(s.length===g.length&&s.length>0);let y=nonnegative(tail);
  for(let j=s.length-1;j>=0;j--){const q=positive(g[j]).add(y);y=q.div(O().add(positive(s[j]).mul(q)));}return y;
}
function nodeSolve(s,g,tail=0,v0=1){
  const M=s.length;assert(M===g.length&&M>0);tail=nonnegative(tail);v0=positive(v0);
  const A=o.zeros(M),rhs=o.zeros(M,1);
  for(let j=0;j<M;j++){
    const c=O().div(positive(s[j]));positive(g[j]);
    A[j][j]=A[j][j].add(c).add(g[j]);
    if(j===0)rhs[0][0]=new Cut(c.mul(v0));
    else {A[j-1][j-1]=A[j-1][j-1].add(c);A[j][j-1]=A[j][j-1].sub(c);A[j-1][j]=A[j-1][j].sub(c);}
  }
  A[M-1][M-1]=A[M-1][M-1].add(tail);
  const volt=o.mul(o.inverse(A),rhs);assert(o.equal(o.mul(A,volt),rhs));
  const v=volt.map(row=>{assert(row[0].turn.zero());return row[0].rad;});
  let power=tail.mul(v[M-1].pow(2));
  const flows=s.map((sj,j)=>(j?v[j-1]:v0).sub(v[j]).div(sj));
  for(let j=0;j<M;j++){
    eq(flows[j],Q(g[j]).mul(v[j]).add(j+1<M?flows[j+1]:tail.mul(v[j])));
    assert(Z().le(v[j])&&v[j].le(v0));
    power=power.add(Q(s[j]).mul(flows[j].pow(2))).add(Q(g[j]).mul(v[j].pow(2)));
    counts.node_balance_checks++;
  }
  eq(v0.mul(flows[0]),power);counts.nodal_solves++;counts.power_balance_checks++;
  return {response:flows[0].div(v0),power,voltages:v};
}
function controlGate(Rs,Rp,Rstar,b,z){
  [Rs,Rp,Rstar,b,z]=[Rs,Rp,Rstar,b,z].map(x=>positive(x));
  eq(Rstar.div(Rs),b.add(1).mul(z));eq(Rp.div(Rstar),b.mul(z));
  return true;
}
function componentBounds(s,g,e){
  e=nonnegative(e);assert(e.le(1)&&!e.eq(1));
  // Both physical series and shunt resistors vary within +/-e.
  const hi=O().add(e),lo=O().sub(e);
  return {lower:reduceTarget(s.map(x=>Q(x).mul(hi)),g.map(x=>Q(x).div(hi))),
    upper:reduceTarget(s.map(x=>Q(x).mul(lo)),g.map(x=>Q(x).div(lo)))};
}
function decimalInterval(lo,hi,places=12){
  const d=10n**BigInt(places),a=Q(lo),b=Q(hi);assert(a.n>=0n&&a.le(b));
  return {lower:new F(a.n*d/a.d,d).toString(),upper:new F((b.n*d+b.d-1n)/b.d,d).toString(),decimal_places:places};
}
check('native_pair_identities_replayed_without_engine_fork',()=>{
  for(const [a,x] of [['12/7','5/12'],['5/7',1]]){
    const w=r2.pairWitness(a,a,{u:1,v:x}),P=r1.emk();
    for(const [name,cert] of Object.entries(w.proofs)){P.replay(cert);replays.push({cell:a,tail:String(x),identity:name,sha256:sha(bytes(cert)),replay:'REPLAY_MATCH'});}
    const bad=JSON.parse(JSON.stringify(w.proofs.response));assert(bad.steps.length);bad.steps[0].rule='NONEXISTENT_RULE';
    rejects('altered_native_witness_'+a,()=>P.replay(bad));
  }
});
check('independent_node_elimination_matches_native_apertures',()=>{
  for(const b of ['1/3','5/7',2])for(const z of ['1/4',1,'7/4'])for(const M of [1,2,5])for(const tail of [0,'2/3']){
    const f=fixture(b,z,M),n=nodeSolve(f.s,f.g,tail,'3/7');
    eq(n.response,r2.fromTail(f.cells,tail));eq(n.response,reduceTarget(f.s,f.g,tail));
  }
  const s=['2/3','7/5','1/4','9/8'],g=['3/5','1/7','4/3','5/9'],cells=s.flatMap((x,j)=>[O().div(x),O().div(g[j])]);
  for(const tail of [0,'5/4',3])eq(nodeSolve(s,g,tail).response,r2.fromTail(cells,tail));
});
check('open_short_and_finite_loads_have_correct_endpoint_parity',()=>{
  for(const z of ['1/4',1,'11/8','7/4'])for(const M of [1,3,8]){
    const f=fixture('5/7',z,M),e=r2.enclose(f.cells).interval;eq(reduceTarget(f.s,f.g),e.lower);
    // Last node clamped to zero is a true short, not a huge numerical load.
    const short=M===1?O().div(f.s[0]):reduceTarget(f.s.slice(0,-1),f.g.slice(0,-1),O().div(f.s[M-1]));
    eq(short,e.upper);
    for(const tail of ['1/100',1,100]){const x=reduceTarget(f.s,f.g,tail);assert(e.lower.le(x)&&x.le(e.upper));counts.tail_controls++;}
  }
});
check('physical_units_and_common_cell_control_are_exact',()=>{
  const b=Q(protocol.b);
  for(const Rstar of [1,1000,'7/3'])for(const row of protocol.controls){
    const f=fixture(b,row.z,2,Rstar);controlGate(f.Rs,f.Rp,Rstar,b,row.z);
    eq(f.Rs.mul(f.Rp),Q(Rstar).pow(2).mul(b).div(b.add(1)));
    const x=reduceTarget(f.s,f.g);eq(Q(protocol.source_voltage_volt).mul(x).div(Rstar).mul(Rstar).div(protocol.source_voltage_volt),x);
  }
});
check('component_interval_encloses_every_two_section_corner',()=>{
  const f=fixture('5/7','11/8',2),e=Q('1/10'),bounds=componentBounds(f.s,f.g,e),x=reduceTarget(f.s,f.g);
  eq(bounds.lower,x.div(O().add(e)));eq(bounds.upper,x.div(O().sub(e)));
  for(let mask=0;mask<16;mask++){
    const factors=Array.from({length:4},(_,j)=>O().add(e.mul((mask>>j)&1?1:-1)));
    const s=f.s.map((v,j)=>v.mul(factors[j])),g=f.g.map((v,j)=>v.div(factors[j+2]));
    const actual=nodeSolve(s,g).response;assert(bounds.lower.le(actual)&&actual.le(bounds.upper));counts.component_corners++;
  }
});
check('source_amplitude_changes_current_but_not_return_coefficient',()=>{
  const f=fixture('5/7','11/8',3),a=nodeSolve(f.s,f.g,0,1),b=nodeSolve(f.s,f.g,0,'7/3');
  eq(a.response,b.response);eq(b.power,a.power.mul(Q('7/3').pow(2)));
});
const prediction=[];
check('prior_two_probe_fixture_and_third_polynomial_are_consumed',()=>{
  eq(protocol.b,'5/7');eq(protocol.delta,'3/4');
  const u=Q('6/5'),v=Q('2/5'),A=O().sub(u.pow(2)),B=O().sub(v.pow(2)),b=Q(protocol.b);
  eq(A.mul(B).mul(2).mul(b.pow(2)).add(A.add(B).mul(2).sub(u.mul(B)).sub(v.mul(A)).mul(b)).add(Q(2).sub(u).sub(v)),0);
  eq(u.div(O().add(b.mul(A))).sub(1),protocol.delta);
  for(const row of protocol.controls){
    const f=fixture(b,row.z,protocol.sections,protocol.reference_resistance_ohm),e=r2.enclose(f.cells).interval;
    const x=nodeSolve(f.s,f.g).response;eq(x,e.lower);
    if(row.prior_completed_response){eq(r2.periodTwoResidual(f.cells[0],f.cells[1],row.prior_completed_response),0);assert(e.lower.le(row.prior_completed_response)&&Q(row.prior_completed_response).le(e.upper));}
    if(row.id==='reference')assert(e.lower.le(1)&&O().le(e.upper));
    if(row.id==='third'){
      const poly=x=>Q(55).mul(x.pow(2)).add(Q(56).mul(x)).sub(132);
      assert(poly(e.lower).n<0n&&poly(e.upper).n>0n);assert(e.width.le('3/100000000000'));
      assert(Q('11216063508/10000000000').le(x)&&e.upper.le('11216063510/10000000000'));
      assert.deepEqual(row.prior_polynomial_coefficients,['55','56','-132']);
    }
    const currentFactor=Q(protocol.source_voltage_volt).mul(1000000).div(protocol.reference_resistance_ohm);
    const cb=componentBounds(f.s,f.g,protocol.relative_component_tolerance);
    prediction.push({id:row.id,z:row.z,sections:protocol.sections,series_ohm:f.Rs,shunt_ohm:f.Rp,
      finite_open_response:x,completed_enclosure:{lower:e.lower,upper:e.upper,width:e.width},
      completed_decimal_enclosure:decimalInterval(e.lower,e.upper),
      ideal_current_microamp:decimalInterval(x.mul(currentFactor),e.upper.mul(currentFactor),10),
      component_current_microamp:decimalInterval(cb.lower.mul(currentFactor),cb.upper.mul(currentFactor),8)});
  }
});
check('finite_baseline_and_gain_ratio_error_are_retained',()=>{
  const base=prediction[0].finite_open_response,e0=O().sub(base);assert(e0.n>0n&&e0.le(1));
  for(const row of prediction){
    const {lower,upper,width}=row.completed_enclosure,ratio=row.finite_open_response.div(base);
    const bound=width.add(upper.mul(e0)).div(O().sub(e0));
    const far=max(ratio.sub(lower).abs(),ratio.sub(upper).abs());assert(far.le(bound));
    for(const gain of ['1/3',7])eq(Q(gain).mul(row.finite_open_response).div(Q(gain).mul(base)),ratio);
    row.finite_reference_ratio=ratio;row.completed_ratio_error_bound=bound;
    const e=Q(protocol.relative_component_tolerance),lo=ratio.mul(O().sub(e)).div(O().add(e)),hi=ratio.mul(O().add(e)).div(O().sub(e));
    row.component_ratio_enclosure=decimalInterval(lo,hi,10);
  }
});
check('port_is_finite_at_the_native_inverse_pole',()=>{
  const P=r1.emk(),L=r1.normal(P.word(['K']).times(P.word(['R']))),I=P.one();
  assert.equal(P.reduce(I.plus(L).times(I.minus(L))).normal.terms.size,0);
  assert(I.minus(L).terms.size>0);eq(r2.periodTwoResidual('12/7','5/7',1),0);
  const f=fixture('5/7',1,8);assert(nodeSolve(f.s,f.g).response.n>0n);
  assert.equal(protocol.fundamental_constant_selected,false);
});
check('false_physical_controls_and_readout_shortcuts_are_rejected',()=>{
  const b=Q('5/7'),z=Q('11/8'),f=fixture(b,z,3,1000),base=fixture(b,1,3,1000);
  rejects('uniform_resistor_scaling_is_not_native_common_probe',()=>controlGate(base.Rs.mul(z),base.Rp.mul(z),1000,b,z));
  rejects('voltage_only_probe_does_not_change_admittance',()=>controlGate(base.Rs,base.Rp,1000,b,z));
  rejects('wrong_common_product_units',()=>controlGate(f.Rs,f.Rp,1,b,z));
  rejects('finite_port_is_not_exact_completed_cut',()=>eq(reduceTarget(base.s,base.g),1));
  rejects('odd_and_even_tail_bounds_may_not_be_swapped',()=>{const e=r2.enclose(f.cells).interval;eq(reduceTarget(f.s,f.g),e.upper);});
  rejects('additive_offset_does_not_cancel_like_gain',()=>{const x=reduceTarget(f.s,f.g),x0=reduceTarget(base.s,base.g);eq(x.add(1).div(x0.add(1)),x.div(x0));});
  rejects('configuration_gain_drift_does_not_cancel',()=>{const x=reduceTarget(f.s,f.g),x0=reduceTarget(base.s,base.g);eq(x.mul(2).div(x0),x.div(x0));});
  rejects('negative_component',()=>fixture(-1,1,2));
  rejects('invalid_tolerance',()=>componentBounds(f.s,f.g,1));
});
check('source_and_claim_status_mutations_are_rejected',()=>{
  const p=JSON.parse(JSON.stringify(pins));p.sources.local.files[0].sha256='0'.repeat(64);rejects('wrong_source_pin',()=>verifyPins(p));
  const p2=JSON.parse(JSON.stringify(pins));p2.sources.local.files[0].path='../outside';rejects('source_path_escape',()=>verifyPins(p2));
  const l=JSON.parse(JSON.stringify(ledger));l.physical_interface.status='NATIVE_DERIVED';rejects('physical_assumption_promoted_to_native_result',()=>verifyLedger(l));
  const l2=JSON.parse(JSON.stringify(ledger));l2.claims[0].written_section_sha256='0'.repeat(64);rejects('altered_written_proof_binding',()=>verifyLedger(l2));
  const p3=JSON.parse(JSON.stringify(protocol));p3.experiment_performed=true;rejects('constructed_fixture_promoted_to_measurement',()=>validateProtocol(p3));
});
verifyPins(pins);
const result={schema:'extra-ideas.r46.return-port-certificate.v1',status:'PASS_R46_NATIVE_RETURN_PORT_ADAPTER',
  data_kind:protocol.data_kind,source_pins_sha256:sha(fs.readFileSync(path.join(here,'SOURCE_PINS.json'))),
  protocol_sha256:sha(fs.readFileSync(protoPath)),written_results:5,exact_checks:checks,check_count:checks.length,counts,
  symbolic_replays:replays,symbolic_replay_count:replays.length,rejected_false_alternatives:negative,
  predictions:prediction,native_engine_modified:false,physical_experiment_performed:false,
  physical_interface:'Ideal steady-DC ohmic interpretation; target equations proved separately',
  full_emk_hardware_realization:false,physical_constants_selected:false,formal_proof_assistant_verified:false};
const output=bytes(result),digest=sha(output);
const record={schema:'extra-ideas.r46.verification.v1',status:result.status,written_results:5,
  source_pins_sha256:result.source_pins_sha256,native_certificate_sha256:digest,
  exact_check_groups:checks.length,counters:counts,native_symbolic_replays:replays.length,
  negative_controls_rejected:Object.keys(negative).length,native_engine_modified:false,
  physical_experiment_performed:false,physical_constants_selected:false,formal_proof_assistant_verified:false,
  scope:'Five written application proofs, exact finite target/port comparisons and explicit ideal electrical interface; no measured device result.'};
const files=[[path.join(here,'CERTIFICATE.json'),output],[path.join(here,'EXPECTED.sha256'),Buffer.from(digest+'\n')],
  [path.join(root,'04-operator-evolution/R46_VERIFICATION.json'),bytes(record)]];
if(options['--write'])for(const [p,b] of files)fs.writeFileSync(p,b);
else for(const [p,b] of files)assert(fs.readFileSync(p).equals(b),'Frozen evidence mismatch: '+p);
console.log(JSON.stringify({status:result.status,groups:checks.length,counts,native_replays:replays.length,
  negative_controls:Object.keys(negative).length,certificate_sha256:digest,
  third:prediction.find(x=>x.id==='third').ideal_current_microamp},null,2));
