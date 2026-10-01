'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),vm=require('node:vm');
const args=process.argv.slice(2),get=flag=>{const k=args.indexOf(flag);return k<0?null:args[k+1];};
const root=get('--rkf-root');
if(!root)throw new Error('Usage: node verify_r14.cjs --rkf-root PATH [--output PATH]');
const pins=JSON.parse(fs.readFileSync(path.join(__dirname,'R14_SOURCE_PINS.json'),'utf8'));
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const external=pins.runtime_sources.map(p=>{
  const b=fs.readFileSync(path.join(root,p.path));
  if(sha(b)!==p.sha256)throw new Error('Pinned canonical source mismatch: '+p.path);
  return {...p,matched:true};
});
const native=require(path.resolve(root,'operator_foundation/core/native_operator.cjs'));
const build=require('./native_source_foundation.cjs'),tests=require('./test_native_source_foundation.cjs');
const result=tests(build(native),native);
const input=fs.readFileSync(path.join(__dirname,'native_source_foundation.cjs'),'utf8');
const mutations=[
  ['wrong memory feedback sign','t.b.pow(2).neg().mul(t.a.pow(k-1-j))','t.b.pow(2).mul(t.a.pow(k-1-j))'],
  ['wrong hidden force sign','const force=t.b.neg().mul','const force=t.b.mul'],
  ['forced balanced preparation','const mean=sum(rows.map((x,j)=>weights[j].mul(x.y)));','const mean=f(0);'],
  ['quartic event rule','q:b.pow(2).div(d)','q:b.pow(4).div(a.pow(4).add(b.pow(4)))'],
  ['wrong conjugation projector','cut(z).add(cut(z).dagger()).div(2)','cut(z).sub(cut(z).dagger()).div(2)'],
  ['dropped record normalization','content:x.amplitude.norm2().div(totalEnergy)','content:x.amplitude.norm2()'],
  ['lost phase orientation','const phaseNumerator=t.a.mul(db).sub(t.b.mul(da));','const phaseNumerator=t.a.mul(db).add(t.b.mul(da));'],
  ['wrong local log Hessian factor','v[j].pow(2).div(x).sub(w[j])','v[j].pow(2).div(x).mul(2).sub(w[j])']
];
const controls=mutations.map(([name,before,after])=>{
  if(!input.includes(before))throw new Error('Mutation source missing: '+name);
  const context={module:{exports:{} }};vm.runInNewContext(input.replace(before,after),context);
  let rejected=false,reason='';try{tests(context.module.exports(native),native);}catch(e){rejected=true;reason=String(e.message).slice(0,160);}
  if(!rejected)throw new Error('Undetected mathematical mutation: '+name);
  return {name,rejected,reason};
});
const localFiles=['native_source_foundation.cjs','test_native_source_foundation.cjs','verify_r14.cjs','R14_SOURCE_PINS.json',
  '../02-relational-response/NATIVE_SOURCE_FOUNDATION_R14.md','../02-relational-response/DERIVATION_SOURCE_LEDGER_R14.md'];
const local=localFiles.map(p=>({path:p,sha256:sha(fs.readFileSync(path.join(__dirname,p)))}));
const report={schema:'extra-ideas.r14.finite-native-audit.v1',status:'PASS_R14_FINITE_NATIVE_SOURCE_FOUNDATION',
  test_count:result.tests.length,tests:result.tests,case_counts:result.counts,mathematical_mutations:controls,
  canonical_runtime_sources:external,local_sources:local,
  evidence_scope:'Exact finite source, history, record-content, covariance and metric identities; written infinite and analytic theorems remain separate.',
  full_source_selection_proved:false,physical_event_actualization_proved:false,
  all_repositories_pure_native:false,new_lean_formalization:false};
const out=get('--output');if(out)fs.writeFileSync(path.resolve(out),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({status:report.status,tests:report.test_count,case_counts:report.case_counts,mutations_rejected:controls.length,output:out||null}));
