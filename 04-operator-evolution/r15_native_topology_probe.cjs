'use strict';
// Reuse the canonical RKF symbolic engine. No engine source is copied here.
const path = require('node:path');
const fs = require('node:fs');
const arg = name => { const at = process.argv.indexOf(name); return at < 0 ? null : process.argv[at + 1]; };
const root = arg('--rkf-root');
if (!root) throw new Error('--rkf-root required');
const home = path.join(path.resolve(root), 'operator_foundation');
const w = require(path.join(home, 'core/workbench.cjs'));
const p = require(path.join(home, 'core/paninian_operator.cjs'));
const example = JSON.parse(fs.readFileSync(path.join(home, 'examples/emk_job.json'), 'utf8'));
const word = (...tokens) => ({word: tokens});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const comm = (left, right) => ({op: 'commutator', left, right});
const equal = 'EQUAL_IN_DECLARED_QUOTIENT', distinct = 'DISTINCT_IN_DECLARED_QUOTIENT';
const contract = {
  schema: 'extra-ideas.r15.operator-source-audit.v1',
  free: {presentation: {tokens: ['c', 'k'], rules: []}, tasks: [
    {name: 'J_squared_not_identity_without_closure_law', left: word('c', 'k', 'c', 'k'), right: word(), expected: distinct},
    {name: 'order_of_cut_and_self_cut_not_automatically_equal', left: word('c', 'k'), right: word('k', 'c'), expected: distinct},
  ]},
  kir: {presentation: example.presentation, tasks: [
    {name: 'reflection_squared', left: word('K', 'K'), right: word(), expected: equal},
    {name: 'rotation_squared', left: word('R', 'R'), right: scale(-1, word()), expected: equal},
    {name: 'four_rotations', left: word('R', 'R', 'R', 'R'), right: word(), expected: equal},
    {name: 'mixed_order_retained', left: word('K', 'R'), right: scale(-1, word('R', 'K')), expected: equal},
    {name: 'reflection_is_not_rotation', left: word('K'), right: word('R'), expected: distinct},
    {name: 'reflection_is_not_projector', left: word('K', 'K'), right: word('K'), expected: distinct},
    {name: 'self_commutator_zero', left: comm(word('K'), word('K')), right: scale(0, word()), expected: equal},
    {name: 'mixed_curvature_nonzero', left: comm(word('R'), word('K')), right: scale(2, word('R', 'K')), expected: equal},
  ]},
};
const digest = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: digest, input: contract}) + '\n');
} else {
  if (arg('--expected-input-sha256') !== digest) throw new Error('R15 contract hash mismatch');
  const groups = {};
  for (const [name, input] of Object.entries({free: contract.free, kir: contract.kir})) {
    const s = w.presentation(input.presentation);
    groups[name] = {presentation_sha256: s.hash(), audit: s.audit(), results: input.tasks.map(task => {
      const result = p.proveEquality(w.parseExpression(s, task.left), w.parseExpression(s, task.right));
      if (result.status !== task.expected || !s.replay(result.certificate)) throw new Error(task.name);
      return {name: task.name, expected: task.expected, result, replay: 'REPLAY_MATCH'};
    })};
  }
  const model = w.buildModel(example.presentation);
  if (model.basis.status !== 'COMPLETE_FINITE_NORMAL_BASIS' || model.basis.dimension !== 4) throw new Error('Native KIR basis mismatch');
  const altered = JSON.parse(JSON.stringify(groups.kir.results[1].result.certificate));
  altered.output.terms = [[[], ['1', '0']]];
  let rejected = false;
  try { w.presentation(contract.kir.presentation).replay(altered); } catch (_) { rejected = true; }
  if (!rejected) throw new Error('Altered rotation-square witness accepted');
  process.stdout.write(JSON.stringify({status: 'PASS_SCOPED_NATIVE_REPLAYS', input: contract,
    input_sha256: digest, groups, regular_basis: {dimension: model.basis.dimension, words: model.basis.words},
    negative_controls: {altered_rotation_square_rejected: rejected},
    scope: 'Consequences of the declared native presentations, not a derivation of their relations from the upload A0-A2.'}, null, 2) + '\n');
}
