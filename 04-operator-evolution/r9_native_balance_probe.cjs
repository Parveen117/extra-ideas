'use strict';
// A research caller of the unchanged canonical RKF engine, not an engine fork.
const fs = require('node:fs');
const path = require('node:path');
function argument(name) {
  const at = process.argv.indexOf(name);
  return at < 0 ? null : process.argv[at + 1];
}
const root = argument('--rkf-root');
if (!root) throw new Error('--rkf-root must name the pinned RKF checkout');
const engine = path.join(path.resolve(root), 'operator_foundation');
const w = require(path.join(engine, 'core', 'workbench.cjs'));
const p = require(path.join(engine, 'core', 'paninian_operator.cjs'));
const example = JSON.parse(fs.readFileSync(path.join(engine, 'examples', 'emk_job.json'), 'utf8'));
const word = (...tokens) => ({word: tokens});
const add = (...args) => ({op: 'add', args});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const subtract = (a, b) => add(a, scale(-1, b));
const multiply = (...args) => ({op: 'multiply', args});
const comm = (left, right) => ({op: 'commutator', left, right});
const alpha = a => multiply(word('K'), a, word('K'));
const even = a => scale(['1/2', '0'], add(a, alpha(a)));
const odd = a => scale(['1/2', '0'], subtract(a, alpha(a)));
const a = word('A'), b = word('B'), zero = scale(0, word());
const tasks = [
  {name: 'balanced_idempotence', left: even(even(a)), right: even(a)},
  {name: 'branch_exchange_invariance', left: even(alpha(a)), right: even(a)},
  {name: 'odd_sector_annihilated', left: even(odd(a)), right: zero},
  {name: 'product_memory_defect', left: subtract(even(multiply(a, b)), multiply(even(a), even(b))),
    right: multiply(odd(a), odd(b))},
  {name: 'commutator_memory_defect', left: subtract(even(comm(a, b)), comm(even(a), even(b))),
    right: comm(odd(a), odd(b))},
  {name: 'odd_second_moment_retained', left: even(multiply(odd(a), odd(a))),
    right: multiply(odd(a), odd(a))},
].map(task => ({...task, expected: 'EQUAL_IN_DECLARED_QUOTIENT'}));
const f = comm(word('K'), word('R'));
const contract = {
  schema: 'extra-ideas.r9-cut-balance.v1',
  generic: {presentation: {tokens: ['A', 'B', 'K'], rules: [
    {id: 'KK', lhs: ['K', 'K'], rhs: [[[], 1]], source: 'R9 native involution contract: K squared equals I'},
  ]}, tasks},
  kir: {presentation: example.presentation, tasks: [
    {name: 'raw_KIR_curvature_nonzero', left: f, right: zero, expected: 'DISTINCT_IN_DECLARED_QUOTIENT'},
    {name: 'radial_tangential_curvature_cancels', left: even(f), right: zero, expected: 'EQUAL_IN_DECLARED_QUOTIENT'},
    {name: 'tangential_tangential_curvature_survives', left: even(comm(word('R'), word('R', 'K'))),
      right: scale(-2, word('K')), expected: 'EQUAL_IN_DECLARED_QUOTIENT'},
  ]},
};
const hash = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: hash, input: contract}) + '\n');
} else {
  const expected = argument('--expected-input-sha256');
  function validateInput(input) {
    if (!expected || p.digest(input) !== expected) throw new Error('Externally pinned R9 input mismatch');
  }
  validateInput(contract);
  const groups = {};
  for (const name of ['generic', 'kir']) {
    const input = contract[name], s = w.presentation(input.presentation);
    const results = input.tasks.map(task => {
      const result = p.proveEquality(w.parseExpression(s, task.left), w.parseExpression(s, task.right));
      if (result.status !== task.expected || !s.replay(result.certificate)) throw new Error('R9 proof mismatch: ' + task.name);
      return {name: task.name, expected: task.expected, result, replay: 'REPLAY_MATCH'};
    });
    groups[name] = {presentation_sha256: s.hash(), audit: s.audit(), results,
      finite_carrier_assumed: false};
  }
  const alteredInput = JSON.parse(JSON.stringify(contract));
  alteredInput.generic.tasks[0].right = word('A');
  let inputRejected = false;
  try { validateInput(alteredInput); } catch (_) { inputRejected = true; }
  const alteredWitness = JSON.parse(JSON.stringify(groups.generic.results[0].result.certificate));
  alteredWitness.output.terms = [[[], ['1', '0']]];
  let witnessRejected = false;
  try { w.presentation(contract.generic.presentation).replay(alteredWitness); } catch (_) { witnessRejected = true; }
  if (!inputRejected || !witnessRejected) throw new Error('An altered R9 input or witness was accepted');
  process.stdout.write(JSON.stringify({input: contract, input_sha256: hash, groups,
    negative_controls: {altered_contract_rejected: inputRejected, altered_witness_rejected: witnessRejected},
    status: 'PASS_NATIVE_CUT_BALANCE_REPLAYS'}, null, 2) + '\n');
}
