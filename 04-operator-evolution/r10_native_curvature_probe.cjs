'use strict';
// Research caller: presentation, reduction and replay remain in RKF's engine.
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
const word = (...tokens) => ({word: tokens});
const add = (...args) => ({op: 'add', args});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const subtract = (a, b) => add(a, scale(-1, b));
const multiply = (...args) => ({op: 'multiply', args});
const comm = (left, right) => ({op: 'commutator', left, right});
const alpha = a => multiply(word('K'), a, word('K'));
const even = a => scale(['1/2', '0'], add(a, alpha(a)));
const odd = a => scale(['1/2', '0'], subtract(a, alpha(a)));
const a = word('A'), b = word('B'), c = word('C'), d = word('D');
const zero = scale(0, word());
const cyclic = (f) => add(f(a, b, c), f(b, c, a), f(c, a, b));
const expectedEquality = 'EQUAL_IN_DECLARED_QUOTIENT';
// Cross-carrier inputs use the engine's typed Expr JSON interface. The
// workbench composite AST parser is restricted to its one-object interface.
const typed = (...terms) => ({from: 'V', to: 'W', terms});
const plus = ['1', '0'], minus = ['-1', '0'];
const contract = {
  schema: 'extra-ideas.r10-native-curvature-descent.v1',
  graded: {
    input_format: 'workbench_expression_AST',
    presentation: {tokens: ['A', 'B', 'C', 'D', 'K'], rules: [
      {id: 'KK', lhs: ['K', 'K'], rhs: [[[], 1]],
        source: 'R9/R10 constant native cut involution contract: K squared equals I'},
    ]},
    tasks: [
      {name: 'full_constant_Bianchi_Jacobi',
        left: cyclic((x, y, z) => comm(x, comm(y, z))), right: zero},
      {name: 'visible_Bianchi_plus_hidden_source',
        left: add(cyclic((x, y, z) => comm(even(x), even(comm(y, z)))),
                  cyclic((x, y, z) => comm(odd(x), odd(comm(y, z))))), right: zero},
      {name: 'Riemann_distortion_commutator_expansion',
        left: comm(add(a, c), add(b, d)),
        right: add(comm(a, b), comm(a, d), comm(c, b), comm(c, d))},
    ].map(task => ({...task, expected: expectedEquality})),
  },
  descent: {
    input_format: 'native_typed_Expr_JSON',
    presentation: {objects: ['V', 'W'], tokens: [
      {name: 'B1', from: 'W', to: 'W', charge: '0'},
      {name: 'B2', from: 'W', to: 'W', charge: '0'},
      {name: 'C', from: 'V', to: 'W', charge: '0'},
      {name: 'A1', from: 'V', to: 'V', charge: '0'},
      {name: 'A2', from: 'V', to: 'V', charge: '0'},
    ], rules: [
      {id: 'intertwine_1', lhs: ['C', 'A1'], rhs: [[['B1', 'C'], 1]],
        source: 'Declared typed observer contract C A1 = B1 C'},
      {id: 'intertwine_2', lhs: ['C', 'A2'], rhs: [[['B2', 'C'], 1]],
        source: 'Declared typed observer contract C A2 = B2 C'},
    ]},
    tasks: [
      {name: 'typed_observer_curvature_descent',
        left: typed([['C', 'A1', 'A2'], plus], [['C', 'A2', 'A1'], minus]),
        right: typed([['B1', 'B2', 'C'], plus], [['B2', 'B1', 'C'], minus])},
      {name: 'typed_ordered_composition_descent',
        left: typed([['C', 'A2', 'A1'], plus]),
        right: typed([['B2', 'B1', 'C'], plus])},
    ].map(task => ({...task, expected: expectedEquality})),
  },
};
const hash = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: hash, input: contract}) + '\n');
} else {
  const expected = argument('--expected-input-sha256');
  function validateInput(input) {
    if (!expected || p.digest(input) !== expected) throw new Error('Externally pinned R10 input mismatch');
  }
  validateInput(contract);
  const groups = {};
  for (const name of ['graded', 'descent']) {
    const input = contract[name], s = w.presentation(input.presentation), audit = s.audit();
    if (audit.status !== 'CONFLUENT_BY_CHECKED_DIAMONDS') throw new Error('Uncertified R10 presentation');
    const parse = name === 'graded' ? ast => w.parseExpression(s, ast) : ast => s.fromJSON(ast);
    const results = input.tasks.map(task => {
      const result = p.proveEquality(parse(task.left), parse(task.right));
      if (result.status !== task.expected || !s.replay(result.certificate)) throw new Error('R10 proof mismatch: '+task.name);
      return {name: task.name, expected: task.expected, result, replay: 'REPLAY_MATCH'};
    });
    groups[name] = {presentation_sha256: s.hash(), audit, results,
      finite_carrier_assumed: false};
  }
  const alteredInput = JSON.parse(JSON.stringify(contract));
  alteredInput.graded.tasks[1].right = word('A');
  let inputRejected = false;
  try { validateInput(alteredInput); } catch (_) { inputRejected = true; }
  const alteredWitness = JSON.parse(JSON.stringify(groups.graded.results[1].result.certificate));
  alteredWitness.output.terms = [[[], plus]];
  let witnessRejected = false;
  try { w.presentation(contract.graded.presentation).replay(alteredWitness); } catch (_) { witnessRejected = true; }
  if (!inputRejected || !witnessRejected) throw new Error('An altered R10 contract or witness was accepted');
  process.stdout.write(JSON.stringify({input: contract, input_sha256: hash, groups,
    negative_controls: {altered_contract_rejected: inputRejected, altered_witness_rejected: witnessRejected},
    status: 'PASS_NATIVE_CURVATURE_DESCENT_REPLAYS'}, null, 2)+'\n');
}
