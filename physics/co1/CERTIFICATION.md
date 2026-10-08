# CO1 reproduction certificate

Research input: `92e321c0c3b5c026fea357f3f595d07da8f69ce2` on
`physics-ph1-ph2-response-holonomy`.

Run with Python 3.12, SymPy 1.14.0 and mpmath 1.3.0:

```bash
python -B physics/co1/verify.py --publications-root /path/to/pinned/Publications
```

The Publications source is pinned to
`f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56`.
`SOURCE_PINS.json` hashes the unchanged CO1 proof, program, tests and frozen
result, its TP1/CV1 runtime dependencies, and the cited cosmology source.
The verifier replays all five symbolic statement groups, compares the complete
result with `CO1_RESULT.json`, runs the four existing tests (including the
wrong-law rejection), and rechecks source identity after execution. It writes
no research evidence. The dedicated workflow covers this working branch.

This certifies reproducibility of the stated symbolic identities within CO1's
declared homogeneous model. The matter density remains a supplied result of
the cited action, parameters remain free, and the redshift factorization is
restricted to constant expansion rate with `0 <= beta < 1`. This is not
empirical validation or proof-assistant verification. Earlier physics stages
are not promoted to certified status by this check.
