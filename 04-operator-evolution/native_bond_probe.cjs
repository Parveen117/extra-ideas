'use strict';
/** R4 adapter: consume the pinned RKF solver and proof replayer unchanged.
 * Run through verify_r4.py, which checks the four upstream runtime hashes.
 */
const path = require('node:path');
const crypto = require('node:crypto');

function run(rkfRoot) {
  const root = path.resolve(rkfRoot);
  const {F} = require(path.join(root, 'operator_foundation/core/native_operator.cjs'));
  const r1 = require(path.join(root, 'research/recognition_return/return_solver.cjs'));
  const source = require(path.join(root, 'research/recognition_return/r2/variable_return.cjs'));
  const P = r1.emk(), I = P.one(), R = P.word(['R']), K = P.word(['K']);
  const L = r1.normal(K.times(R));
  if (P.audit().status !== 'CONFLUENT_BY_CHECKED_DIAMONDS') {
    throw new Error('The source presentation did not pass its existing audit');
  }

  function identities(value) {
    const x = F.of(value);
    const B = I.plus(L.scale(x)).scale(new F(1, 2));
    const J = K.times(B).minus(B.times(K));
    const defect = B.times(B).minus(B);
    const claims = {
      complementary_exchange: K.times(B).times(K).minus(I.minus(B)),
      raw_commutator: J.minus(R.scale(x)),
      bond_defect: defect.minus(I.scale(x.pow(2).sub(1).div(4))),
      coupled_defects: J.times(J).plus(I).plus(defect.scale(4)),
    };
    const proofs = {};
    for (const [name, expression] of Object.entries(claims)) {
      const witness = P.reduce(expression, {witness: true});
      if (witness.normal.terms.size || P.replay(witness.certificate) !== true) {
        throw new Error('Native bridge proof failed: ' + name);
      }
      proofs[name] = {
        replayed: true, normal_terms: 0, steps: witness.certificate.steps.length,
        certificate_sha256: crypto.createHash('sha256')
          .update(JSON.stringify(witness.certificate)).digest('hex'),
      };
    }
    return {x, proofs};
  }

  function prefix(pattern, count) {
    const cells = Array.from({length: count}, (_, j) => pattern[j % pattern.length]);
    const packet = source.enclose(cells);
    return {pattern, count, actual_zero_tail: source.fromTail(cells),
      interval: packet.interval, aperture_height: packet.apertureHeight};
  }

  const designs = [source.synthesize('1', '1'), source.synthesize('1', '2'),
    source.synthesize('6/5', '20')];
  for (const design of designs) {
    if (!source.periodTwoResidual(...design.period, design.desiredReturn).zero()) {
      throw new Error('The unchanged source rejected its periodic response');
    }
  }
  const prefixes = {
    cut_32_short: prefix(['3', '2'], 6),
    cut_32_long: prefix(['3', '2'], 12),
    noncut_30_20: prefix(['30', '20'], 12),
    accidental_finite_cut: prefix(['1'], 1),
    reopened_same_profile: prefix(['1'], 3),
  };
  // Reuse a native Schur witness, rather than implementing another eliminator.
  const pairWitness = source.pairWitness('3', '3', {u: '1', v: '2/3'});
  return {
    protocol: 'EXTRA_IDEAS_NATIVE_BOND_BRIDGE_R4',
    runtime: process.version,
    presentation_sha256: P.hash(),
    algebra_checks: ['0', '1', '6/5', '3/7'].map(identities),
    source_schur_proof_names: Object.keys(pairWitness.proofs),
    source_schur_response: pairWitness.result,
    designs, prefixes,
    full_infinite_inverse_claimed: false,
    physical_constant_selected: false,
  };
}

module.exports = {run};
if (require.main === module) {
  if (process.argv.length !== 3) throw new Error('Supply the pinned RKF checkout path');
  process.stdout.write(JSON.stringify(run(process.argv[2])) + '\n');
}
