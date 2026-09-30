# R7 — Vault tensor-foundation source audit

Research owner: Monty Dabas. Audit date: 30 September 2026.

The purpose is to develop a consistent EMK tensor foundation before applying it to physics. Existing tensor manuscripts are inputs to assess, not axioms that make every later claim true. New definitions and proofs are in [EMK_TENSOR_CALCULUS_R7.md](EMK_TENSOR_CALCULUS_R7.md).

## 1. Coverage and reading depth

The audit retrieved **143 mathematical source files** at immutable Vault commits: 74 current files and 69 recovered/branch files. The latter include all **43 LaTeX files in the recovered `04_EMK_CORE`**, plus selected UGD and Recognition-Seam Calculus foundations. The remaining core spreadsheet was not used as a tensor foundation. The public lineage check additionally retrieved six files from Publications.

This does **not** mean every Vault file, branch or page was read in full. Central tensor documents were read in depth; large books were inspected through relevant body sections and section maps. Current revision logs and peripheral application manuscripts were retrieved/indexed without a full proof audit. The [source index](../04-operator-evolution/R7_SOURCE_PINS.json) distinguishes:

- `FULL_FOUNDATION_BODY`: the mathematical body of a central foundation source was read, including its relevant proofs and limits.
- `TARGETED_BODY_AND_OUTLINE`: relevant definitions/claims were read and the wider outline inspected; other chapters remain outside this audit.
- `OUTLINE_AND_SELECTED_PASSAGES`: scope, headings and selected opening/definition passages were inspected; no full validation claim.
- `RETRIEVED_INDEXED`: fetched and pinned for the corpus index; not presented as a complete reading.

Reading depth is deliberately conservative. Retrieval and hashing certify which bytes were available, not the truth of a manuscript's claims. Number-theory, device, field-equation and experimental proposals are outside the R7 validation scope.

| Source layer | Immutable Vault commit | Role |
| --- | --- | --- |
| Current `master` math corpus | `7e89daed803570ab93f68a2cb636d4f8091000cd` | Current geometry, RTC appendices and clock-free RSC |
| Recovered EMK/UGD/RSC foundations | `55294a55e9bad98c3c8e4a6b9ec86be3d6c5b478` | Original tensor, Jacobian, geometry, algebra and calculus drafts |
| EMK RTC completion branch | `55518c9f3f45ccfbdb7d885fff32abb4f5408686` | Standalone master tensor/connection source |

The original core is absent from current master, so examining only that branch would miss the requested tensor drafts. This audit explicitly includes the recovery branch. No private source bodies are republished here; the index contains provenance metadata and hashes.

## 2. Central sources and what they actually supply

| Source path or family | Existing content | R7 treatment |
| --- | --- | --- |
| `04_EMK_CORE/emk_tensor.tex` | E/M tensor proposal, metric, seam and mixed curvature | Supply typed fibers and cross-fiber pairing; examine multilinearity and scalar-commutator claims |
| `04_EMK_CORE/emk_tensor_jacobian.tex` | Conjugation/cut swap, sector splitting and rotation/Jacobian rules | Use a real-linear carrier; fix the sector/bridge interpretation |
| `04_EMK_CORE/emk_geometry.tex` | Seam coordinates, metric, Laplacian and topology proposals | Keep only results with declared geometry and volume assumptions |
| `04_EMK_CORE/EMK_Recognition_Tensor_Calculus_Crown_v1.tex` and completion master | Structured master packet, RTC channels, recognition/Smriti derivative | Type each component; preserve existing recognition-policy boundaries |
| `04_EMK_CORE/emk_maths...tex` and foundational patches | Primitive state/cut language and intended reference independence | Treat linear, norm and analytic realizations as admissions needing proofs |
| `04_EMK_CORE/Morphic algebra complete.tex` | Partial operator system, collapse quotient and overlay axioms | Do not equate overlay with a multilinear tensor product without a representation theorem |
| `04_EMK_CORE/Morphic calculus complete.tex` | Clock/measure-relative derivatives and optional smooth representation | Separate increments from rates and check denominator/domain hypotheses |
| `04_EMK_CORE/ONSAGER GEOMETRIC (4).tex` | Hessian lambda metric, Riemann and four-coordinate tensor proposals | Targeted metric/curvature sections; no physical field-equation certification |
| `04_EMK_CORE/info-cur duality.tex` | Proposed information-gradient tensor relation | Identify it as an extra equation, not a consequence of tensor notation |
| `05_UGD_ENGINE/UGD_Tensor_Calculus_Algebraic_Integration_Theorem.tex` | Tensor operations, connections, exterior calculus and integration dictionary | Existing proposal credited; replace schematic claims by precise hypotheses and finite checks |
| UGD K-connection, Lie/BCH, exterior and admissible-module sources | Phase/scale/sheet operations and curvature/memory modules | Declare the admitted scalar carrier and distinguish discrete k from a continuous relaxation |
| RSC Chapters 2–4, 10–11, 14, 20–21 | Product rules, projections, crown connection, KIR, curvature and cohomology | Keep product defects, cross terms and the square-zero prerequisite |
| Current RTC and geometry appendices | Updated KIR basis, channel curvature, warped metrics and helical memory | Prefer their explicit realization/closure boundaries over earlier universal assertions |
| Current clock-free RSC sections | Naturality defect, functorial composition, lawful ledgers, product coupling | Credit existing finite derivative and conditional Leibniz theory |

## 3. Concrete mathematical corrections

These are local checks under the expressions actually used in the drafts. They do not establish or refute every possible future EMK realization.

### A. Multilinearity and E/M contraction

A tensor is multilinear on specified vector fibers, not on an arbitrary manifold as a set. A coefficient depending on the test arguments may destroy multilinearity. An E-covector cannot evaluate an M-vector merely because their component arrays have equal length. R7 requires B:M→E and proves invariance of `omega_e^T B v_m`. It refuses an unbridged cross-fiber contraction.

### B. Sector intersection and diagonal support

For an involution J, its +1 and −1 eigenspaces intersect only at zero in characteristic zero: if Jv=v and Jv=−v, then 2v=0. A nonzero fixed seam therefore cannot simultaneously be identified with `{e=m}` when e and m are declared to belong to these opposite subspaces. A graph of an admitted alignment map or a separate seam equation is needed.

A smooth coefficient supported only on the diagonal of a genuine positive-dimensional product vanishes by continuity. “Parallel seam sector” must therefore mean a typed projector block, not a nonzero smooth function existing only on that diagonal.

### C. The proposed scalar mixed curvature

For scalar kappa,

\[
(\nabla_X\kappa)(\nabla_Y\kappa)
-(\nabla_Y\kappa)(\nabla_X\kappa)=0.
\]

Multiplying this by a metric does not produce nonzero curvature. A noncommuting operator-valued replacement needs its own type and law. Also the trace of a proposed mixed Ricci tensor `1/2(dkappa tensor dkappa − g ||dkappa||²)` in dimension n is `(1−n)||dkappa||²/2`, not generally the stated positive scalar. R7 does not import either formula as a physical curvature theorem.

### D. A cut-involution trace cannot supply arbitrary rotation

If J²=I and the derivative is compatible with the usual endomorphism trace, differentiating gives `(nabla J)J+J(nabla J)=0`. Cyclicity gives `tr(J nabla J)=0`. Thus `omega_J=tr(J nabla J)/2` is identically zero under those assumptions; it cannot be asserted as a general nonzero rotation connection without changing the construction.

### E. A flat product metric does not produce a polar volume term

For `g=dr²+du²` in those coordinates, the metric Laplacian is `partial_r²+partial_u²`. The term `r^(-1)partial_r` requires a different volume/metric, such as a polar realization. Higher-dimensional normal directions likewise require angular components; `dkappa tensor dkappa` alone is degenerate along directions tangent to its level sets.

### F. Cut and Join from one generator commute

With an invertible U, define Cut=U−U^{-1} and Join=U+U^{-1}. Their commutator is zero because U commutes with its inverse. The exponential specialization in the old explicit K-curvature/helical sources therefore cannot yield the stated nonzero commutator. At generator zero its claimed identity has left side zero and right side −2I. R7 includes an exact check. Nonzero ordered curvature must involve genuinely different noncommuting transports or varying coefficients.

### G. Bianchi applies to the full connection

For a genuine associative connection, `F=dA+A wedge A` obeys `d_A F=0`. A nonzero projected Bianchi residual may expose omitted channels, noncompatible projections or a different ledger comparison. It is not failure of the identity for the full connection. In the finite sector R7 proves and checks the corresponding tetrahedron cancellation even when individual triangle curvatures are nonzero.

### H. Curved transport is not automatically a cochain complex

The old recognition topology chapter declares `H=kernel(d)/image(d)`. This quotient first needs `image(d)` contained in `kernel(d)`. R7's explicit bundle-valued difference satisfies `D_1D_0s=F s`, so curved transport does not supply that requirement automatically. Flat coefficients, an invariant curvature-null module, or a separately proved square-zero recognition differential is needed before cohomology and characteristic/index claims.

### I. Projection and quotient algebra require compatibility

For a projection P, the class `{A:PA=0}` is a **right** ideal in the full endomorphism algebra: P(AB)=0 for any B. It is generally not a left ideal, since PCA may be nonzero. A two-sided quotient algebra needs an admissible operator family preserving the relevant kernel. R7 retains full carrier order and the existing leakage correction instead of trusting a projected zero commutator.

### J. A Hessian needs an affine connection, and a metric needs nondegeneracy

The ordinary coordinate Hessian of a scalar is not a tensor under general nonlinear coordinate changes. For Phi=x²/2 and x=y², transforming the tensor `g_xx=1` gives `g_yy=4y²`, whereas taking the ordinary second derivative of Phi(y)=y⁴/2 gives 6y². A Hessian metric needs a specified affine connection/affine chart structure, or a covariant Hessian definition.

Twice differentiable does not ensure a nondegenerate or positive metric. The examined lambda Hessian section itself gives a negative PP component under its stated positive thermodynamic quantities, so positive Riemannian signature is not automatic. The equation of state, independent state coordinates and signature must be declared before using its inverse and Riemann/Ricci contractions. The four thermodynamic variables cannot be assumed independent on every equilibrium state space.

### K. Information-gradient tensor relations are extra equations

Defining positive I, J, K does not force `d ln I tensor d ln J + d ln J tensor d ln K + d ln K tensor d ln I=0`. For example, positive scalar functions I=J=exp(x), K=1 give `dx tensor dx`, which is nonzero. This illustrates the missing implication; it is not a physical equation-of-state counterexample. An admitted constitutive constraint might enforce such a relation, but it must be proved separately. Contraction loses information and cannot be reversed to recover the full tensor equation.

### L. A stationary clock does not define a zero derivative

For a clock-relative quotient `Delta f/Delta g`, zero clock increment makes normalization undefined. It does not set every internal variation to zero. A Radon–Nikodym formulation requires measures and absolute continuity of the measure being differentiated; an arbitrary scalar function is not itself a measure without another construction. R7 keeps its finite increment independent of such analytic admissions.

### M. Flux, curvature and winding have different types

Raw matrix curvature, its trace, ordered loop holonomy and integer sheet winding are not interchangeable. An arbitrary matrix surface integral is not automatically an integer. A classical characteristic form needs an invariant polynomial/trace, a genuine connection and the appropriate degree and integrality hypotheses. R7 supplies a separate integer winding register and does not infer it from a raw curvature integral.

## 4. Public lineage and duplication boundary

Publications was checked at `e1dc4e3773f063f14e56cf30222c8e8a504519cd`. Its unchanged EMK2 native carrier, EMKT1 tensor/time certificate and EMKT2 ordered-transport certificate are consumed directly in the exact tests. Existing KIR multiplication, cut-swap representation, RTC cross-curvature and ordered-loop witnesses are credited to those sources. The public recognition-geometry guide was also retrieved for lineage. RKF was checked at `3cc5a33b05c16d59c90994ddda69dedc0d392424`; prior R4 lineage already credits its KIR, compression/leakage and iota constructions.

R7 therefore does not claim a newly discovered KIR algebra, a newly discovered general tensor product, the first RSC derivative, or the first curvature/holonomy calculus. What is added here is a checked **typed tensor implementation on that source carrier**, a full rank-parity theorem, the exact simultaneous fixed-metric obstruction in the standard two-mode flow sector, and a linked distinction between tensor return, carrier return and independent sheet closure. Their general mathematical ingredients are familiar representation/bundle theory; their consistent EMK specification and tests are the development delivered here.

## 5. What must follow before physics

Complete the native direction carrier and soldering rule, choose a metric policy consistent with the admitted transport, then build torsion/Ricci and the higher exterior/Hodge operators on that sector. Prove smooth realization and analytic existence where needed. Only then add an action, units, constitutive dynamics and falsifiable physical predictions. Existing speculative field equations are not certified merely by their use of tensor notation.
