# EMK master tensor: corrected assessment and source checks

Review date: 1 October 2026. Research owner: Monty Dabas.

The EMK master tensor, RTC connection, recognition-time sector, seam topology, metric geometry and advanced certificates already exist. The earlier roadmap gave too little weight to this work. R7 adds a typed finite implementation on one admitted KIR carrier; it does not establish that the rest of EMK lacks a tensor foundation.

## Reading scope and reproducibility

This review follows full-body reading of the standalone EMK Recognition Tensor Calculus master manuscript, its tensor/time, algebraic-core and topology appendices, and the original crown manuscript. It also reads the complete executable EMK-T1, EMK-T2, EMK-G1, EMK-G2, EMK-G3, CID-1, QTH-1 and EMK-TOP-1 certificates, and the RKF cut-covariance, metric-graph boundary-pairing and finite Onsager-Hodge theorems. This is a review of the central mathematical spine. It does not claim that every file in every connected repository has been read in full.

Public sources are pinned to [Publications at e1dc4e3](https://github.com/Parveen117/Publications/tree/e1dc4e3773f063f14e56cf30222c8e8a504519cd) and [Recognition-Kernel-Framework at 3cc5a33](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424). The review record stores SHA256 and Git blob identities for 43 Publications files and four RKF context documents. Private source bodies are not reproduced here.

All **191 existing tests across eight certificate families pass**. Their eight saved result hashes match, and each family's existing regeneration test passes. The verifier runs the unchanged no-argument test functions through Python's standard unittest runner. These successful checks coexist with the three cross-function or interpretation findings below.

Run against source checkouts at the pinned commits:

```bash
python3.12 -B 04-operator-evolution/verify_emk_master_review.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

The optional `--output` argument writes a review record outside the upstream checkouts. [Source pins](../04-operator-evolution/EMK_MASTER_REVIEW_SOURCE_PINS.json) and [verification record](../04-operator-evolution/EMK_MASTER_REVIEW_VERIFICATION.json) accompany the [verifier](../04-operator-evolution/verify_emk_master_review.py). All recorded R1–R7 hashes remain intact. A review pass means these checks and counterexamples reproduce; it does not mean every source claim is proved.

## Existing results that must be credited

| Layer | Existing mathematical content | Scope |
| --- | --- | --- |
| [EMK-T1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/emkt1_master_tensor_and_time.py) | Master channel-curvature decomposition; individually flat channels can have nonzero mixed curvature; recognition-time memory and closure; presentation agreement when active sectors are preserved | Exact finite matrix witnesses. The master tuple, channel assignments, tolerances and commit policy are declared structure. |
| [EMK-T2](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/emkt2_time_ordered_transport.py) | Ordered transport, noncommuting path order, exact nilpotent exponentials, refinement for the declared constant-generator partitions, and equal-clock/different-history witnesses | Does not certify general nonnilpotent ordered exponentials or continuum convergence. |
| [EMK-G1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-recognition-geometry/certificates/emkg1_rotational_seam_metric.py) | State-dependent seam metric, determinant and degeneracy, Christoffel symbols, Riemann/Gaussian curvature, and the geodesic seam criterion | Declared rotational metric families and their admissible domains; metric curvature is distinguished from recognition curvature. |
| [EMK-G2](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-recognition-geometry/certificates/emkg2_global_quotient_holonomy.py) / [EMK-G3](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-recognition-geometry/certificates/emkg3_helical_sheet_memory.py) | Global quotient versus double, gluing, Levi-Civita versus recognition holonomy, a mapping-torus lift, flat connection with global monodromy, and visible return with sheet nonreturn | The exact G3 certificate uses a declared abelian translation fibre; the general nonabelian and cut-localized sheet lift is outside its certification. |
| [CID-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/cid1_curvature_information_duality.py) / [QTH-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/qth1_quantum_recognition_information.py) | Exact covariance/Fisher metric, recognized plus discarded information, metric-dependent orthogonal decomposition with metric-independent loop residue; quantum tensor with symmetric information and antisymmetric commutator parts | CID-1 uses a declared finite response algebra. QTH-1 is an exact one-qubit witness family, not a general quantum-estimation theorem. |
| [EMK-TOP-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/recognition-seam-topology/certificates/emktop1_vault_topology_appendix.py) | Finite cochain complex with d1 d0 = 0, rank-one H1, conserved charge on homologous cycles, exact/closed-nonexact/open residues, determinant-charge and sheet-memory witnesses | Explicit finite complex and declared active sectors; no general continuum limit is certified. |
| RKF [cut covariance](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/32_cut_covariance_event_realization_theorem.md), [metric graph](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/35_native_metric_graph_boundary_pairing.md), [Onsager-Hodge](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/thermodynamics/06_onsager_hodge_no_leakage.md) | Conditional Gram/cut-covariance identities, boundary pairing from a native graph metric, quotient decoder theorem, and finite no-leakage under an admitted Hodge grading | Hypotheses and remaining domain-adapter gates are stated in the sources. These structures are already present. |

The native master is a structured recognition object with several typed sectors. It should not be assessed solely as an ordinary two-index matrix. Its component list is a framework declaration; its algebraic consequences and realization choices require their own proofs. The source certificates make this distinction explicitly.

## Correct scope of the R7 metric result

R7 proves that a fixed symmetric bilinear form H on the **standard real two-mode carrier** cannot satisfy both

\[
R^T H+HR=0,\qquad K^T H+HK=0
\]

unless H is zero. These are infinitesimal conditions for the two full continuous flows. The finite KIR generators themselves preserve the Euclidean form. The alternating form is also preserved by the traceless flow sector.

This result does not exclude state-dependent base metrics such as EMK-G1's metric, information metrics such as CID-1/QTH-1, restricted flows, larger carriers, or a moving metric. It does not show that EMK has no metric tensor. Connecting these existing metrics to a particular R7 tensor carrier is a compatibility problem for that chosen realization.

Likewise, even-rank blindness to I versus -I is a property of that tensor representation. The master framework already retains additional holonomy and sheet sectors that can distinguish returns lost by a projection.

## Three concrete findings

### 1. EMK-G3 T2 uses the opposite action from its own deck map

The implemented deck action is

\[
D_n(u,v,f)=(u+nL,v,\rho^{-n}f).
\]

Its support coordinate therefore changes by minus n beta. But `compensated()` returns sigma minus p_sigma u, and both `certify_T2()` and `test_compensated_coordinates_are_exactly_invariant()` test the forward shift sigma plus beta instead of calling `deck()`.

With the actual source constants L = 5 and beta = 2/5, take u = sigma = 0. One actual deck step gives u' = 5 and sigma' = -2/5. The compensated support coordinate changes from 0 to **-4/5**. Thus T2's claim of invariance under the implemented deck action is false despite its existing tests passing. The default half-period phase advance can conceal the phase sign error after reduction modulo six; the support witness has no such ambiguity.

Keeping the inverse deck convention, the consistent compensated coordinates are

\[
\widehat\phi=(\phi+p_\phi u)\bmod K,
\qquad \widehat\sigma=\sigma+p_\sigma u.
\]

The verifier checks this repair on 189 exact actual-deck cases, including positive and negative powers. Alternatively, one could change the deck convention and reconcile all transition signs. The upstream source is unchanged in this review. This finding concerns the coordinate-invariance claim; it does not refute the existence of the mapping torus or the flat/global-monodromy separation.

### 2. EMK-G2 T5 needs a nondegenerate domain check

The holonomy block tests the rational warp A(v) = 1 + c v^2, including c = -2 on the rectangle with v in [-1,1]. There A(-1) = A(1) = -1 and A(0) = 1, so A has interior zeros. At those zeros the metric determinant A^2 vanishes and the asserted coframe is not a frame.

Both polynomial routines still return 16: the formal boundary/area identity is correct. But it does not certify Levi-Civita holonomy across the degenerate points. A geometric test must require that the loop and spanning region stay in a declared nondegenerate domain, here A > 0. For c = -2, the tested interval [0,1/2] is admissible; [-1,1] is not. The connection/curvature formulas on admissible regions remain available.

### 3. EMK-T1 T6's determinant factor does not measure noncommutativity

For invertible T, B1 and B2, the factor actually computed is

\[
C=\frac{\det(TB_1B_2)}{\det(TB_1)\det(TB_2)}
=\frac{1}{\det T}.
\]

This identity is correct. Calling it the failure of transport to commute with the block product is not: T = 2I commutes with every block and gives C = 1/4. Conversely, the determinant-one shear used in the verifier fails to commute with the source block product and gives C = 1.

The factor measures the determinant difference between applying left transport once to a product and applying it separately to both factors. Correct the interpretation, or add a commutator-specific diagnostic if noncommutativity is the intended claim. The determinant arithmetic and its inverse ledger cancellation are unaffected.

## Development direction

Use the existing EMK master as the organizing framework. First reconcile the three identified claims and add tests that compose the actual source maps. Then specify how the chosen carrier's pairing, metric and derivative connect to the existing RTC, geometry, information and topology sectors. Extend only the operations or analytic hypotheses that remain outside those certificates' stated scope.

Christoffel, Riemann, holonomy, time, topology and conditional Hodge work should be credited and reused. Their names alone do not select one universal compatible realization, and the certificates do not claim that every physical identification is already proved. The next task is reconciliation and a precise extension of the existing foundation.
