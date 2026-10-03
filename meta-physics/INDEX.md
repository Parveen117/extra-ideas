# Meta physics: current developments

Research owner: Monty Dabas. Updated 3 October 2026. Python 3.12 only.

The [98-PDF audit](SOURCE_AUDIT.md) is the original reading and triage
snapshot. It is frozen with MP-1's certificate. This live index records
later scoped developments without rewriting that evidence. No complete
source PDF is certified or republished.

| Packet | Actual advance | Evidence | Idea sources |
| --- | --- | --- | --- |
| [MP-1](mp1/THEOREM.md) | Hidden frame sign; exact winding; universal error-tube and sharp sampling conditions for the selected axis observer | 7 written results, 23 exact tests, 6,561 finite perturbations; [certificate](mp1/VERIFICATION.json) | 017, 078, in scope only |
| [MP-2](mp2/THEOREM.md) | An evolving Hessian response preserves native curvature while reducing its response budget; ordered free heat keeps a uniform decay floor; moving work and discarded-record variance have exact balances | 7 written results, 22 exact tests, 624 matrix controls; [certificate](mp2/VERIFICATION.json) | 033, 053, 093, in scope only |

MP-2 combines several suggestions into one missing bridge. Publications
YM54 already gave the frozen response-to-protocol construction. MP-2
selects an explicit time-dependent response law and proves that the
existing curvature/rate connection survives along that law. It imports
the unchanged native engine; it does not redevelop operator algebra.

The selected law projects radial relaxation onto constant oriented shape
area. It is a proved admissible model, not a derivation of nature's unique
constitutive law. Its scalar record accounts for paired response loss;
the separate observed noise term comes from discarded turn records.
The additional signed work term is necessary when the response changes.
At fixed curvature the full observer's best strength ratio is 3:1;
the sign-blind observer prefers balance. These optimize different declared
objectives, not physical coupling constants.

| Source gate from the frozen audit | MP-2 treatment | Still open |
| --- | --- | --- |
| 033: curvature alone does not give a spectral floor | T2/T4 preserve nonzero native area and bound the response budget, giving a common-reference free contraction floor | Arbitrary gauge fields, physical selection, interacting moving source |
| 053: dissipation/flux preservation need explicit conditions | T2/T3 give a selected tangent flow and exact paired loss; T5 includes moving work and actual record variance | General operator-splitting algorithms and physical stochastic calibration |
| 093: sectoral sums and action stationarity do not prove curvature conservation | T1/T2 provide a consistently typed Hessian realization with an explicit conservation proof | Unified force/sector equation and dimensional calibration |
| 051: cross-time memory/passivity | No promotion; scalar account and instantaneous record variance are distinguished from a retarded kernel | An augmented memory realization with a passivity proof |

The main source-law, clock and physical-unit choices stay declared.
R1--R46 and MP-1 remain unchanged. Yang--Mills development has resumed in
[Publications YM55](https://github.com/Parveen117/Publications/blob/4efb9dcdcee5685969f33b572ab04b1e94b8703d/papers/yang-mills-certified-benchmark/YM55_ANISOTROPIC_JOINT_HISTORY_LIMIT.md),
which supplies a stationary anisotropic joint-history limit. MP-2
proves a nonautonomous free contraction estimate, not a stationary
interacting mass gap or a four-dimensional Clay solution.

Replay from the Extra Ideas root with the pinned Publications checkout:

```sh
python3.12 -B meta-physics/mp2/verify.py --publications-root ../Publications --check
```

The verifier checks 28 unchanged runtime modules and the source hashes,
replays MP-1, and checks the complete frozen R1--R46 register. Exact finite
tests accompany written proofs; they are not proof-assistant certification,
independent peer review or experimental validation.
