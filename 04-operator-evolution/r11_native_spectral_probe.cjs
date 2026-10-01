'use strict';
// Calls the unchanged RKF engine; no engine implementation is duplicated.
const path = require('node:path');
function argument(name) {
  const at = process.argv.indexOf(name);
  return at < 0 ? null : process.argv[at+1];
}
const root = argument('--rkf-root');
if (!root) throw new Error('--rkf-root must name the pinned canonical source');
const engine = path.join(path.resolve(root), 'operator_foundation');
const w = require(path.join(engine, 'core', 'workbench.cjs'));
const p = require(path.join(engine, 'core', 'paninian_operator.cjs'));
const word = (...tokens) => ({word: tokens});
const add = (...args) => ({op: 'add', args});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const mul = (...args) => ({op: 'multiply', args});
const comm = (left, right) => ({op: 'commutator', left, right});
const one = word(), zero = scale(0, one);
const P = word('P'), N = word('N'), M = word('M');
const J = add(one, scale(-2, P));
const X = add(one, scale(['1/3', '0'], P));
const Y = add(one, scale(['2/5', '0'], N));
const Xi = add(one, scale(['-1/4', '0'], P));
const Yi = add(one, scale(['-2/5', '0'], N));
const contract = {
  schema: 'extra-ideas.r11-spectral-curvature.v1',
  presentation: {tokens: ['P', 'N', 'M'], rules: [
    {id: 'PP', lhs: ['P', 'P'], rhs: [[['P'], 1]], source: 'Declared hidden projector P^2=P'},
    {id: 'PN', lhs: ['P', 'N'], rhs: [[['N'], 1]], source: 'Declared visible-to-hidden map P N=N'},
    {id: 'NP', lhs: ['N', 'P'], rhs: [], source: 'Declared visible-to-hidden map N P=0'},
    {id: 'NN', lhs: ['N', 'N'], rhs: [], source: 'Strict one-way carrier map N^2=0'},
    {id: 'PM', lhs: ['P', 'M'], rhs: [], source: 'Feedback marker has visible output P M=0'},
    {id: 'MP', lhs: ['M', 'P'], rhs: [[['M'], 1]], source: 'Feedback marker has hidden input M P=M'},
    {id: 'MM', lhs: ['M', 'M'], rhs: [], source: 'Strict reverse-sector marker M^2=0'},
  ]},
  tasks: [
    {name: 'native_curvature_is_the_hidden_recording_map', left: comm(P, N), right: N},
    {name: 'curvature_square_is_zero', left: mul(N, N), right: zero},
    {name: 'derived_sector_cut_is_an_involution', left: mul(J, J), right: one},
    {name: 'curvature_is_cut_odd', left: mul(J, N, J), right: scale(-1, N)},
    {name: 'exchange_balanced_readout_annihilates_this_curvature',
      left: scale(['1/2', '0'], add(N, mul(J, N, J))), right: zero},
    {name: 'feedback_marker_closes_the_curvature_path', left: mul(M, comm(P, N)), right: mul(M, N)},
    {name: 'exact_rational_four_move_return', left: mul(X, Y, Xi, Yi),
      right: add(one, scale(['2/15', '0'], N))},
    {name: 'same_endpoint_histories_remain_collapsed', left: mul(add(one, N), add(one, scale(-1, N))), right: one},
  ].map(task => ({...task, expected: 'EQUAL_IN_DECLARED_QUOTIENT'})).concat([
    {name: 'nonzero_curvature_is_allowed_by_the_declared_algebra', left: N, right: zero,
      expected: 'DISTINCT_IN_DECLARED_QUOTIENT'},
  ]),
};
const hash = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: hash, input: contract})+'\n');
} else {
  const expected = argument('--expected-input-sha256');
  function validateInput(input) {
    if (!expected || p.digest(input) !== expected) throw new Error('Externally pinned R11 input mismatch');
  }
  validateInput(contract);
  const s = w.presentation(contract.presentation), audit = s.audit();
  if (audit.status !== 'CONFLUENT_BY_CHECKED_DIAMONDS') throw new Error('Uncertified R11 presentation');
  const results = contract.tasks.map(task => {
    const result = p.proveEquality(w.parseExpression(s, task.left), w.parseExpression(s, task.right));
    if (result.status !== task.expected || !s.replay(result.certificate)) throw new Error('R11 proof mismatch: '+task.name);
    return {name: task.name, expected: task.expected, result, replay: 'REPLAY_MATCH'};
  });
  const alteredInput = JSON.parse(JSON.stringify(contract));
  alteredInput.tasks[4].right = N;
  let inputRejected = false;
  try { validateInput(alteredInput); } catch (_) { inputRejected = true; }
  const alteredWitness = JSON.parse(JSON.stringify(results[4].result.certificate));
  alteredWitness.output.terms = [[[], ['1', '0']]];
  let witnessRejected = false;
  try { s.replay(alteredWitness); } catch (_) { witnessRejected = true; }
  if (!inputRejected || !witnessRejected) throw new Error('An altered R11 input or witness was accepted');
  process.stdout.write(JSON.stringify({input: contract, input_sha256: hash,
    presentation_sha256: s.hash(), audit, results, finite_carrier_assumed: false,
    negative_controls: {altered_contract_rejected: inputRejected, altered_witness_rejected: witnessRejected},
    status: 'PASS_NATIVE_SPECTRAL_CURVATURE_REPLAYS'}, null, 2)+'\n');
}
