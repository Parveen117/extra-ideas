'use strict';
// Scoped research caller. The only engine home is RKF/operator_foundation.
const fs = require('node:fs');
const path = require('node:path');

function argument(name) {
  const index = process.argv.indexOf(name);
  return index < 0 ? null : process.argv[index + 1];
}
const root = argument('--rkf-root');
if (!root) throw new Error('--rkf-root must name the separately pinned RKF checkout');
const core = path.join(path.resolve(root), 'operator_foundation', 'core');
const w = require(path.join(core, 'workbench.cjs'));
const p = require(path.join(core, 'paninian_operator.cjs'));
const example = JSON.parse(fs.readFileSync(path.join(path.resolve(root), 'operator_foundation', 'examples', 'emk_job.json'), 'utf8'));
const model = w.buildModel(example.presentation);
if (model.basis.status !== 'COMPLETE_FINITE_NORMAL_BASIS' || model.basis.dimension !== 4) {
  throw new Error('The unchanged declared EMK presentation must derive a four-dimensional basis');
}
const rows = [
  {'': 0, R: 1, K: 0, RK: -1},
  {'': 0, R: -1, K: 0, RK: -1},
].map(coefficients => model.basis.words.map(word => {
  const label = word.join('');
  if (!(label in coefficients)) throw new Error('Unexpected native basis word');
  return coefficients[label];
}));
const word = (...tokens) => ({word: tokens});
const f = {op: 'commutator', left: word('K'), right: word('R')};
const job = {
  schema: 'rkf.operator-job.v1',
  presentation: example.presentation,
  tasks: [
    {kind: 'equality', left: f, right: {op: 'scale', coefficient: -2, value: word('R', 'K')}},
    {kind: 'equality', left: {op: 'multiply', args: [word('K'), f, word('K')]},
      right: {op: 'scale', coefficient: 2, value: word('R', 'K')}},
    {kind: 'future_observer', observer: rows},
  ],
};
const inputHash = p.digest(job);
// A second contract proves R8.1 in a free associative presentation with only
// P^2=P. It uses the unchanged rewrite engine and assumes no finite carrier.
const multiply = (...args) => ({op: 'multiply', args});
const subtract = (left, right) => ({op: 'add', args: [left, {op: 'scale', coefficient: -1, value: right}]});
const genericP = word('P'), genericA = word('A'), genericB = word('B');
const genericQ = subtract(word(), genericP);
const compressionContract = {
  presentation: {tokens: ['A', 'B', 'P'], rules: [
    {id: 'PP', lhs: ['P', 'P'], rhs: [[['P'], 1]], source: 'R8.1 declared contract: P squared equals P'},
  ]},
  left: {op: 'commutator', left: multiply(genericP, genericA, genericP), right: multiply(genericP, genericB, genericP)},
  right: subtract(
    multiply(genericP, {op: 'commutator', left: genericA, right: genericB}, genericP),
    subtract(multiply(genericP, genericA, genericQ, genericB, genericP),
      multiply(genericP, genericB, genericQ, genericA, genericP))),
};
const compressionHash = p.digest(compressionContract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: inputHash, input: job,
    compression_input_sha256: compressionHash, compression_input: compressionContract}) + '\n');
} else {
  const expected = argument('--expected-input-sha256');
  if (!expected || inputHash !== expected) throw new Error('Externally pinned R8 job contract mismatch');
  const expectedCompression = argument('--expected-compression-sha256');
  if (!expectedCompression || compressionHash !== expectedCompression) throw new Error('Externally pinned compression contract mismatch');
  const generic = w.presentation(compressionContract.presentation);
  const genericResult = p.proveEquality(w.parseExpression(generic, compressionContract.left),
    w.parseExpression(generic, compressionContract.right));
  if (genericResult.status !== 'EQUAL_IN_DECLARED_QUOTIENT' || !generic.replay(genericResult.certificate)) {
    throw new Error('General idempotent compression identity did not replay');
  }
  const compressionCertificate = {input: compressionContract, input_sha256: compressionHash,
    presentation_sha256: generic.hash(), audit: generic.audit(), result: genericResult,
    replay: 'REPLAY_MATCH', finite_carrier_assumed: false};
  const packet = w.runJob(job);
  const replay = w.replayJob(packet, expected);
  if (packet.status !== 'VERIFIED_FINITE_JOB' ||
      packet.results.slice(0, 2).some(x => x.result.status !== 'EQUAL_IN_DECLARED_QUOTIENT') ||
      packet.results[2].result.initial_rank !== 2 || packet.results[2].result.completed_rank !== 4 ||
      packet.results[2].result.extra_channels !== 2 || replay.status !== 'REPLAY_MATCH') {
    throw new Error('Native equality or minimal future observer certificate failed');
  }
  function mustReject(value) {
    try { w.replayJob(value, expected); } catch (_) { return true; }
    throw new Error('Altered native certificate was accepted');
  }
  const alteredInput = JSON.parse(JSON.stringify(packet));
  alteredInput.input.tasks[0].right.coefficient = -3;
  const alteredResult = JSON.parse(JSON.stringify(packet));
  alteredResult.results[2].result.completed_rank = 2;
  const alteredRewrite = JSON.parse(JSON.stringify(genericResult.certificate));
  alteredRewrite.output.terms = [[[], ['1', '0']]];
  let alteredRewriteRejected = false;
  try { generic.replay(alteredRewrite); } catch (_) { alteredRewriteRejected = true; }
  if (!alteredRewriteRejected) throw new Error('Altered general rewrite proof was accepted');
  const negativeControls = {
    altered_job_contract_rejected: mustReject(alteredInput),
    altered_observer_result_rejected: mustReject(alteredResult),
    altered_general_rewrite_rejected: alteredRewriteRejected,
  };
  process.stdout.write(JSON.stringify({compression_certificate: compressionCertificate,
    packet, replay, negative_controls: negativeControls}, null, 2) + '\n');
}
