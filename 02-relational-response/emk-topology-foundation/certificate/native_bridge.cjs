'use strict';
// Adapter only. Load the unchanged, externally pinned canonical RKF engine.
const path = require('node:path');
const fs = require('node:fs');
const arg = name => { const i = process.argv.indexOf(name); return i < 0 ? null : process.argv[i + 1]; };
if (!arg('--rkf-root')) throw new Error('--rkf-root required');
const home = path.join(path.resolve(arg('--rkf-root')), 'operator_foundation');
const o = require(path.join(home, 'core/native_operator.cjs'));
const p = require(path.join(home, 'core/paninian_operator.cjs'));
const w = require(path.join(home, 'core/workbench.cjs'));
const example = JSON.parse(fs.readFileSync(path.join(home, 'examples/emk_job.json'), 'utf8'));
const word = (...tokens) => ({word: tokens});
const scale = (coefficient, value) => ({op: 'scale', coefficient, value});
const add = (...args) => ({op: 'add', args});
const multiply = (...args) => ({op: 'multiply', args});
const comm = (left, right) => ({op: 'commutator', left, right});
const eq = 'EQUAL_IN_DECLARED_QUOTIENT', ne = 'DISTINCT_IN_DECLARED_QUOTIENT';
const id = word(), K = word('K'), R = word('R');
const P = scale('1/2', add(id, K)), Q = scale('1/2', add(id, scale(-1, K)));
const contract = {
  schema: 'emk-constructive-cut.native-bridge.v1',
  signed_role_construction: {ordered_roles: ['cut', 'self_cut_trace'], parity: [1, -2], exchange: [2, 1]},
  free: {presentation: {tokens: ['c', 'k'], rules: []}, tasks: [
    {name: 'primitive_composite_square_is_not_empty_history', left: word('c','k','c','k'), right: id, expected: ne},
    {name: 'primitive_order_is_retained', left: word('c','k'), right: word('k','c'), expected: ne},
  ]},
  roles: {presentation: example.presentation, tasks: [
    {name: 'reflection_square', left: word('K','K'), right: id, expected: eq},
    {name: 'derived_rotation_square', left: word('R','R'), right: scale(-1,id), expected: eq},
    {name: 'derived_rotation_fourth_power', left: word('R','R','R','R'), right: id, expected: eq},
    {name: 'derived_anticommutation', left: word('K','R'), right: scale(-1,word('R','K')), expected: eq},
    {name: 'reflection_differs_from_rotation', left: K, right: R, expected: ne},
    {name: 'reflection_differs_from_projector', left: word('K','K'), right: K, expected: ne},
    {name: 'mixed_order_defect', left: comm(R,K), right: scale(2,word('R','K')), expected: eq},
    {name: 'self_order_defect_zero', left: comm(K,K), right: scale(0,id), expected: eq},
    {name: 'seam_aperture_idempotent', left: multiply(P,P), right: P, expected: eq},
    {name: 'complement_aperture_idempotent', left: multiply(Q,Q), right: Q, expected: eq},
    {name: 'apertures_orthogonal', left: multiply(P,Q), right: scale(0,id), expected: eq},
  ]},
};
const digest = p.digest(contract);
if (process.argv.includes('--emit-contract')) {
  process.stdout.write(JSON.stringify({input_sha256: digest, input: contract})+'\n');
  process.exit(0);
}
if (arg('--expected-input-sha256') !== digest) throw new Error('Native contract hash mismatch');

// Coordinate arrays are OUTPUT of signed images on the two free generators.
const fromTags = tags => o.matrix([0,1].map(row => tags.map(tag =>
  Math.abs(tag) - 1 === row ? (tag > 0 ? 1 : -1) : 0)));
const h = fromTags(contract.signed_role_construction.parity);
const k = fromTags(contract.signed_role_construction.exchange);
const r = o.mul(k,h), i = o.identity(2), rk = o.mul(r,k);
const requireEquality = (a,b,message) => { if (!o.equal(a,b)) throw new Error(message); };
requireEquality(o.mul(h,h),i,'Parity square');
requireEquality(o.mul(k,k),i,'Exchange square');
requireEquality(o.mul(r,r),o.scale(i,-1),'Derived R square');
requireEquality(o.power(r,4),i,'Derived R fourth power');
requireEquality(o.mul(k,r),o.scale(rk,-1),'Derived anticommutation');
const algebra = o.operatorAlgebra([h,k]);
if (o.rank([i,k,r,rk].map(o.flatten)) !== 4 || algebra.dimension !== 4) throw new Error('Full native algebra rank');
for (const entry of [h,k,r,...algebra.basis].flat(2)) {
  if (!entry.turn.zero()) throw new Error('Role derivation used a pre-supplied complex coefficient');
}
const seam = o.scale(o.add(i,k),'1/2'), complement = o.sub(i,seam);
// All bilinear coefficient cases, hence the cut-corner identity for every A,B.
let cornerBasisCases = 0;
for (const a of [i,k,r,rk]) for (const b of [i,k,r,rk]) {
  requireEquality(o.sub(o.mul(o.mul(seam,o.mul(a,b)),seam),
    o.mul(o.mul(o.mul(seam,a),seam),o.mul(o.mul(seam,b),seam))),
    o.mul(o.mul(o.mul(o.mul(seam,a),complement),b),seam),'Cut-corner basis identity');
  cornerBasisCases++;
}
const groups = {};
for (const [name,input] of Object.entries({free: contract.free, roles: contract.roles})) {
  const s = w.presentation(input.presentation);
  groups[name] = {presentation_sha256:s.hash(), audit:s.audit(), results:input.tasks.map(task => {
    const result = p.proveEquality(w.parseExpression(s,task.left),w.parseExpression(s,task.right));
    if(result.status !== task.expected || !s.replay(result.certificate)) throw new Error(task.name);
    return {name:task.name,expected:task.expected,result,replay:'REPLAY_MATCH'};
  })};
}
const model = w.buildModel(example.presentation);
if(model.basis.status !== 'COMPLETE_FINITE_NORMAL_BASIS' || model.basis.dimension !== 4) throw new Error('Basis mismatch');
const altered = JSON.parse(JSON.stringify(groups.roles.results[1].result.certificate));
altered.output.terms = [[[],['1','0']]];
let rejected = false;
try { w.presentation(contract.roles.presentation).replay(altered); } catch (_) { rejected = true; }
if(!rejected) throw new Error('Altered rotation-square certificate accepted');
process.stdout.write(JSON.stringify({status:'PASS_SCOPED_NATIVE_BRIDGE',input_sha256:digest,
  derived_role_arrays:{h,k,r},full_algebra_dimension:algebra.dimension,
  all_coordinate_turn_components_zero:true,complete_cut_corner_basis_cases:cornerBasisCases,
  regular_basis:{dimension:model.basis.dimension,words:model.basis.words},groups,
  negative_controls:{altered_rotation_square_rejected:rejected},
  scope:'Presentation relations checked on constructed signed-role maps before replay; this does not identify arbitrary raw chi/kappa with them.'},null,2)+'\n');
