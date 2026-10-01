'use strict';
// R12 caller of unchanged canonical RKF proof and N03 observer routines.
const path = require('node:path');
function argument(name) {
  const at = process.argv.indexOf(name);
  return at < 0 ? null : process.argv[at + 1];
}
const root = argument('--rkf-root');
if (!root) throw new Error('Use the pinned --rkf-root');
const home = path.join(path.resolve(root), 'operator_foundation', 'core');
const w = require(path.join(home, 'workbench.cjs'));
const p = require(path.join(home, 'paninian_operator.cjs'));
const native = require(path.join(home, 'native_operator.cjs'));
const word = (...tokens) => ({word: tokens});
const add = (...args) => ({op: 'add', args});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const mul = (...args) => ({op: 'multiply', args});
const one = word(), zero = scale(0, one);
const P = word('P'), N = word('N'), D = word('D'), Q = word('Q');
const plus = add(one, scale(['2/3', '0'], N));
const minus = add(one, scale(['-2/3', '0'], N));
const plusD = add(one, scale(['2/3', '0'], D));
const minusD = add(one, scale(['-2/3', '0'], D));
const gram = add(Q, scale(['4/9', '0'], mul(D, Q, N)));
const cut = add(one, scale(-2, P));
const id4 = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
const shear = (i,j,a=1) => id4.map((row,r) => row.map((x,c) =>
  r === i && c === j ? (a === 1 ? x+1 : a === -1 ? x-1 : a) : x));
const visible = [[1,0,0,0],[0,1,0,0]];
const contract = {
  schema: 'extra-ideas.r12-observer-metric.v1',
  presentation: {tokens: ['P','N','D','Q'], rules: [
    {id:'PP',lhs:['P','P'],rhs:[[['P'],1]],source:'Hidden projector P^2=P'},
    {id:'PN',lhs:['P','N'],rhs:[[['N'],1]],source:'Native one-way record P N=N'},
    {id:'NP',lhs:['N','P'],rhs:[],source:'Native one-way record N P=0'},
    {id:'NN',lhs:['N','N'],rhs:[],source:'Nilpotent native record N^2=0'},
    {id:'DP',lhs:['D','P'],rhs:[[['D'],1]],source:'D represents the admitted transpose of N'},
    {id:'PD',lhs:['P','D'],rhs:[],source:'Adjoint record has visible output'},
    {id:'DD',lhs:['D','D'],rhs:[],source:'Transpose of N is also square-zero'},
  ]},
  tasks: [
    {name:'opposite_sheet_transports_are_exact_inverses',left:mul(plus,minus),right:one},
    {name:'balanced_first_moment_is_identity',left:scale(['1/2','0'],add(plus,minus)),right:one},
    {name:'second_moment_retains_the_hidden_record',
      left:scale(['1/2','0'],add(mul(plusD,Q,plus),mul(minusD,Q,minus))),right:gram},
    {name:'square_zero_stabilizes_the_recursive_record_cost',left:mul(D,gram,N),right:mul(D,Q,N)},
    {name:'derived_metric_solves_the_tagged_event_selection_equation',left:gram,
      right:add(scale(['1/2','0'],Q),scale(['1/4','0'],mul(plusD,gram,plus)),
                scale(['1/4','0'],mul(minusD,gram,minus)))},
    {name:'native_record_curvature_is_cut_odd',left:mul(cut,N,cut),right:scale(-1,N)},
    {name:'balanced_linear_curvature_disappears_in_this_sector',
      left:scale(['1/2','0'],add(N,mul(cut,N,cut))),right:zero},
    {name:'every_nilpotent_sheet_history_adds_signed_record_amplitudes',
      left:mul(add(one,scale(['1/3','0'],N)),add(one,scale(['2/5','0'],N))),
      right:add(one,scale(['11/15','0'],N))},
  ].map(t => ({...t,expected:'EQUAL_IN_DECLARED_QUOTIENT'})).concat([
    {name:'balanced_mean_does_not_determine_the_tagged_metric',left:gram,right:Q,
      expected:'DISTINCT_IN_DECLARED_QUOTIENT'},
  ]),
  observer_cases: [
    {name:'future_depth',seed:[[1,1]],actions:[[[1,0],[0,'1/2']]],expected_rank:2},
    {name:'paired_sheet_with_one_invisible_mode',seed:[...visible,[0,0,1,0]],
      actions:[shear(2,0,'2/3'),shear(2,0,'-2/3')],expected_rank:3},
    {name:'lower_record_without_reverse_feedback',seed:visible,
      actions:[shear(2,0),shear(2,0,-1),shear(3,1),shear(3,1,-1)],expected_rank:2},
    {name:'both_reverse_markers_admitted_as_transport',seed:visible,
      actions:[shear(2,0),shear(2,0,-1),shear(3,1),shear(3,1,-1),
               shear(0,2),shear(0,2,-1),shear(1,3),shear(1,3,-1)],expected_rank:4},
  ],
};
const hash = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256:hash,input:contract})+'\n');
} else {
  const expected = argument('--expected-input-sha256');
  const validateInput = input => {
    if (!expected || p.digest(input) !== expected) throw new Error('R12 input pin mismatch');
  };
  validateInput(contract);
  const s = w.presentation(contract.presentation), audit = s.audit();
  if (audit.status !== 'CONFLUENT_BY_CHECKED_DIAMONDS') throw new Error('R12 presentation audit failed');
  const results = contract.tasks.map(task => {
    const result = p.proveEquality(w.parseExpression(s,task.left),w.parseExpression(s,task.right));
    if (result.status !== task.expected || !s.replay(result.certificate)) throw new Error(task.name);
    return {name:task.name,expected:task.expected,result,replay:'REPLAY_MATCH'};
  });
  const observers = contract.observer_cases.map(c => {
    const result = native.rowClosure(c.seed,c.actions);
    const actions = c.actions.map(a => native.descendedAction(a,result.observer));
    if (result.rank !== c.expected_rank || actions.some(a => a === null)) throw new Error('N03 completion failed');
    return {name:c.name,result,descended_actions:actions,expected_rank:c.expected_rank};
  });
  const alteredInput = JSON.parse(JSON.stringify(contract));
  alteredInput.tasks[2].right = Q;
  let badInput = false;
  try {validateInput(alteredInput);} catch (_) {badInput=true;}
  const badProof = JSON.parse(JSON.stringify(results[2].result.certificate));
  badProof.output.terms = [[[],['1','0']]];
  let badWitness = false;
  try {s.replay(badProof);} catch (_) {badWitness=true;}
  const alteredObserver = JSON.parse(JSON.stringify(observers[3].result));
  alteredObserver.rank = 2;
  const badObserver = p.digest(alteredObserver) !== p.digest(
    native.rowClosure(contract.observer_cases[3].seed,contract.observer_cases[3].actions));
  if (!badInput || !badWitness || !badObserver) throw new Error('An altered R12 contract or result was accepted');
  process.stdout.write(JSON.stringify({input:contract,input_sha256:hash,
    presentation_sha256:s.hash(),audit,results,observers,
    algebra_replays_assume_finite_carrier:false,observer_completions_are_finite:true,
    negative_controls:{altered_input_rejected:badInput,altered_rewrite_rejected:badWitness,
      altered_future_observer_rejected:badObserver},status:'PASS_NATIVE_OBSERVER_METRIC_REPLAYS'},null,2)+'\n');
}
