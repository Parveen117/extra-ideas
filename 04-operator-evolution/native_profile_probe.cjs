'use strict';
/** R5 adapter. Source recursion, synthesis and enclosures remain upstream.
 * The Python entry point checks the four original runtime hashes first.
 */
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');

function run(root, requests) {
  root = path.resolve(root);
  const {F} = require(path.join(root, 'operator_foundation/core/native_operator.cjs'));
  const source = require(path.join(root, 'research/recognition_return/r2/variable_return.cjs'));
  if (!Array.isArray(requests)) throw new Error('Supply an array of requests');
  const results = requests.map(request => {
    const z = F.of(request.z ?? '1');
    if (z.n <= 0n) throw new Error('The common product multiplier must be positive');
    const scaled = cells => cells.map(value => F.of(value).mul(z));
    let result;
    switch (request.kind) {
      case 'finite': {
        const cells = scaled(request.cells);
        result = {value: source.fromTail(cells, request.tail ?? '0'),
          transfer: source.transfer(cells)};
        break;
      }
      case 'prefix_interval':
        result = source.prefixEnclosure(scaled(request.prefix), request.lower, request.upper);
        break;
      case 'periodic_prefix': {
        const count = request.tail_cells;
        if (!Number.isSafeInteger(count) || count < 1 || count > 8192 ||
            !Array.isArray(request.pattern) || !request.pattern.length) {
          throw new Error('Supply a nonempty period and a finite tail aperture budget');
        }
        const cells = (request.prefix ?? []).concat(Array.from({length: count},
          (_, j) => request.pattern[j % request.pattern.length]));
        result = source.enclose(scaled(cells));
        break;
      }
      case 'synthesize':
        result = source.synthesize(request.desired, request.second);
        break;
      case 'schur': {
        const value = F.of(request.cell).mul(z);
        const witness = source.pairWitness(value, value, {u: '1', v: request.tail});
        const proofs = Object.fromEntries(Object.entries(witness.proofs).map(([name, proof]) =>
          [name, {replayed_by_source: true, certificate_sha256: crypto.createHash('sha256')
            .update(JSON.stringify(proof)).digest('hex')}]
        ));
        result = {response: witness.result, proofs,
          presentation_sha256: witness.presentation_sha256};
        break;
      }
      default:
        throw new Error('Unknown query kind: ' + request.kind);
    }
    return {id: request.id, kind: request.kind, result};
  });
  return {protocol: 'EXTRA_IDEAS_PROFILE_RECOVERY_R5', runtime: process.version, results};
}

module.exports = {run};
if (require.main === module) {
  if (process.argv.length !== 3) throw new Error('Supply the pinned RKF checkout path');
  process.stdout.write(JSON.stringify(run(process.argv[2],
    JSON.parse(fs.readFileSync(0, 'utf8')))) + '\n');
}
