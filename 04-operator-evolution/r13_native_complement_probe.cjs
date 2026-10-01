'use strict';
// Calls the existing KIR presentation and N03 observer routine unchanged.
const fs = require('node:fs');
const path = require('node:path');
function argument(name) {
  const at = process.argv.indexOf(name);
  return at < 0 ? null : process.argv[at+1];
}
const root = argument('--rkf-root');
if (!root) throw new Error('Use the pinned --rkf-root');
const core = path.join(path.resolve(root), 'operator_foundation', 'core');
const w = require(path.join(core, 'workbench.cjs'));
const p = require(path.join(core, 'paninian_operator.cjs'));
const native = require(path.join(core, 'native_operator.cjs'));
const source = JSON.parse(fs.readFileSync(path.join(path.resolve(root),
  'operator_foundation', 'examples', 'emk_job.json'), 'utf8'));
const word = (...tokens) => ({word: tokens});
const add = (...args) => ({op:'add',args});
const scale = (coefficient,value) => ({op:'scale',coefficient,value});
const mul = (...args) => ({op:'multiply',args});
const one = word(), zero = scale(0,one);
const R = word('R'), J = scale(-1,word('R','K'));
const P = scale(['1/2','0'],add(one,J));
const Q = scale(['1/2','0'],add(one,scale(-1,J)));
const T = add(scale(['3/5','0'],one),scale(['4/5','0'],R));
const Td = add(scale(['3/5','0'],one),scale(['-4/5','0'],R));
const survival = mul(P,T,P);
// Reuse the separately replayed P^2=P and PTP=(3/5)P below before
// squaring the branch norm; this avoids expanding duplicate projector sums.
const exit2 = scale(['9/25','0'],mul(Q,T,P));
const exit2D = scale(['9/25','0'],mul(P,Td,Q));
const excursion = mul(P,T,Q,T,P);
const contract = {
  schema:'extra-ideas.r13-complement-memory.v1',
  presentation:source.presentation,
  tasks:[
    {name:'observer_involution_is_derived_from_native_RK',left:mul(J,J),right:one},
    {name:'visible_cut_is_idempotent',left:mul(P,P),right:P},
    {name:'native_rotation_preserves_the_positive_chart_norm',left:mul(Td,T),right:one},
    {name:'visible_single_step_is_the_cosine_channel',left:survival,right:scale(['3/5','0'],P)},
    {name:'complement_return_is_negative_visible_identity',left:mul(P,R,Q,R,P),right:scale(-1,P)},
    {name:'coherent_two_step_return_keeps_the_sine_correction',left:mul(P,T,T,P),right:scale(['-7/25','0'],P)},
    {name:'two_step_compression_defect_is_exact_complement_memory',
      left:add(mul(P,T,T,P),scale(-1,mul(survival,survival))),right:excursion},
    {name:'first_exit_weight_is_the_native_squared_branch_amplitude',
      left:mul(exit2D,exit2),right:scale(['1296/15625','0'],P),
      depends_on:['visible_cut_is_idempotent','visible_single_step_is_the_cosine_channel']},
    {name:'cut_transport_noncommutation_is_the_sine_channel',
      left:add(mul(P,T),scale(-1,mul(T,P))),right:scale(['-4/5','0'],word('K'))},
    {name:'cut_commutator_square_is_complement_strength',
      left:mul(add(mul(P,T),scale(-1,mul(T,P))),add(mul(P,T),scale(-1,mul(T,P)))),
      right:scale(['16/25','0'],one)},
  ].map(t=>({...t,expected:'EQUAL_IN_DECLARED_QUOTIENT'})).concat([
    {name:'the_complement_excursion_is_not_zero',left:excursion,right:zero,
      expected:'DISTINCT_IN_DECLARED_QUOTIENT'},
  ]),
  observer_cases:[
    {name:'rotation_future_reveals_the_hidden_component',seed:[[1,0]],
      actions:[[['3/5','-4/5'],['4/5','3/5']]],expected_rank:2},
    {name:'no_sine_coupling_has_no_future_hidden_response',seed:[[1,0]],
      actions:[[[1,0],[0,1]]],expected_rank:1},
    {name:'one_way_hidden_record_does_not_return_to_the_observer',seed:[[1,0]],
      actions:[[[1,0],[1,1]]],expected_rank:1},
  ],
};
const hash = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256:hash,input:contract})+'\n');
} else {
  const expected = argument('--expected-input-sha256');
  function validateInput(input) {
    if (!expected || p.digest(input)!==expected) throw new Error('R13 native input pin mismatch');
  }
  validateInput(contract);
  const s = w.presentation(contract.presentation), audit = s.audit();
  if (audit.status!=='CONFLUENT_BY_CHECKED_DIAMONDS') throw new Error('R13 source audit failed');
  const replayedNames = new Set();
  const results = contract.tasks.map(task=>{
    if ((task.depends_on || []).some(name=>!replayedNames.has(name))) throw new Error('Missing prior proof dependency');
    const result = p.proveEquality(w.parseExpression(s,task.left),w.parseExpression(s,task.right));
    if (result.status!==task.expected || !s.replay(result.certificate)) throw new Error(task.name);
    replayedNames.add(task.name);
    return {name:task.name,expected:task.expected,depends_on:task.depends_on || [],result,replay:'REPLAY_MATCH'};
  });
  const observers = contract.observer_cases.map(c=>{
    const result = native.rowClosure(c.seed,c.actions);
    const actions = c.actions.map(a=>native.descendedAction(a,result.observer));
    if (result.rank!==c.expected_rank || actions.some(a=>a===null)) throw new Error(c.name);
    return {name:c.name,result,descended_actions:actions,expected_rank:c.expected_rank};
  });
  const changedInput = JSON.parse(JSON.stringify(contract));
  changedInput.tasks[5].right = scale(['9/25','0'],P);
  let inputRejected = false;
  try {validateInput(changedInput);} catch (_) {inputRejected=true;}
  const badRewrite = JSON.parse(JSON.stringify(results[6].result.certificate));
  badRewrite.output.terms = [[[],['1','0']]];
  let rewriteRejected = false;
  try {s.replay(badRewrite);} catch (_) {rewriteRejected=true;}
  const changedObserver = JSON.parse(JSON.stringify(observers[0].result));
  changedObserver.rank = 1;
  const observerRejected = p.digest(changedObserver)!==p.digest(
    native.rowClosure(contract.observer_cases[0].seed,contract.observer_cases[0].actions));
  if (!inputRejected || !rewriteRejected || !observerRejected) throw new Error('Altered R13 evidence accepted');
  process.stdout.write(JSON.stringify({input:contract,input_sha256:hash,
    presentation_sha256:s.hash(),audit,results,observers,
    abstract_rewrites_assume_finite_carrier:false,observer_completions_are_finite:true,
    negative_controls:{altered_input_rejected:inputRejected,altered_rewrite_rejected:rewriteRejected,
      altered_future_observer_rejected:observerRejected},status:'PASS_NATIVE_COMPLEMENT_MEMORY_REPLAYS'},null,2)+'\n');
}
