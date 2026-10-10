# Target reference, spectral reach and the next directional reading

10 October 2026. Review requested by Monty Dabas after YC34 (`3bc32aa`).
Scope: the pasted Claude exchange about moving the TVSP reference to a
target, RMG2, and NSB2. This is a correction and usable research map, not
a new Yang–Mills gap stage. Existing theorem packets remain unchanged.

**Decision.** Use a target-relative residual and record which admitted
operations can change it. Keep the distinction between relabelling zero,
changing the observation frame and modifying an actual response/dynamics.
NSB2 supports the spectral-channel obstruction and the force-data remedy.
Several stronger statements in the pasted interpretation do not follow.

## Sources actually inspected

- Publications, branch `codex/nsb2-nonlinear-response-closure`, pinned commit
  `bec73061e92b56f4772d5a31063dde4da0bcaa65`:
  `papers/nonlinear-spectral-blindness/NSB2_NONLINEAR_RESPONSE_CLOSURE.md`.
  Git blob `eb71d326439675f3b21d0a0a7237bdc3bb35d8e2`.
  The theorem text, hypotheses, proofs and YM boundary were read. Its old
  certificate suites were not rerun.
- extra-ideas at `3bc32aa`: UP1, UP3, UP6, CA1, LS1, LS2, YC20, YC23,
  YC30, YC33 and YC34, with their current interpretation/context notes.
- The RMG2 manuscript was not separately retrieved. Its eigenspace claim
  is checked here directly from NSB2 and elementary matrix algebra; this
  review does not certify all of RMG2 or the pasted quotation of its T2.
- Primary geometric background: Burt Totaro, *The curvature of a Hessian
  metric*, section2, <https://www.math.ucla.edu/~totaro/papers/public_html/hessian.pdf>.
  The calibrated Hessian metric and its third-derivative curvature are
  distinct from a quantum Hamiltonian and its excitation spectrum.

## 1. Which statements survive the audit?

| Pasted claim | Verified reading or correction |
|---|---|
| UP1 moves the silent point to depth1/k | Correct for a degree-k homogeneous potential and the declared diagonal interpolation. The point lies inside the depth interval only when1/k is admissible. Degree determines that zero; changing a label is a different operation. |
| A zero centre means an isotropic response with Cp=Cv | False in general. A scalar potential reading can vanish with a nonzero, anisotropic, coupled Hessian. Cp=Cv means the off-diagonal response vanishes in CA1's specified chart, not that both diagonal entries are equal. |
| Every polynomial of a2x2 H stays in span{I,H} | Correct pointwise. For symmetric H this preserves its eigenspaces; eigenvalues can merge. The coefficients in the span can depend on position when H is a field. |
| A commuting target is a function of H | For a nonscalar symmetric2x2 H, yes pointwise: its commutant among symmetric matrices is span{I,H}. At repeated eigenvalues commutation is insufficient. Global scalar-function compatibility is another restriction. |
| Same principal angle means the lambda tower can reach the target | Eigenspace agreement is necessary for a nondegenerate spectral target, not sufficient for a specified one-parameter tower or allowed iteration protocol. Check its coefficient and positivity constraints. |
| The UP6 bracket is already the defect [H,H_target] | These are different types: UP6 is a commutator of thermodynamic vector fields; the latter is a matrix commutator. A representation/intertwiner is needed before equating them. |
| A Hessian-preserving square must be flat/separable | NSB2 T2 proves flatness under C4 regularity, positive Hessian, fixed nonzero constant lambda and identities on an open set. Separation is local on simple-spectrum regions. |
| Therefore a potential-preserving tower cannot do anything | False. Eigenvalues and constitutive response can change nontrivially even in the flat, separable class. UP3's conditional-derivative recursion is not the spectral-power operation of NSB2. |
| F(U) supplies a new directional response | Correct: its Hessian has the rank-one term F'' dU tensor dU. It can rotate eigenspaces when that term fails to commute with H. It need not rotate them. |
| Choosing F is equivalent to choosing any target | False without a compatibility test. At one state the available response is aH+b gg^T. Across a neighborhood a,b must come from the same one-variable F and obey its integrability relations. |
| One vector measurement suffices in two dimensions | Correct for NSB2's local contract: simple spectrum, known eigenvalue fields and first derivatives, calibrated nonzero v, and the complete vector Hv. A scalar quadratic reading alone is insufficient in general. |
| NSB2 rules out reaching the YM continuum by a lambda route | False. It rules out a particular unjustified Hessian extension and notes the limitation of one YM82 sufficient bound. It proves no general continuum impossibility. |

Equilibrium also needs its declared meaning. A thermodynamic equation-of-state
manifold can parametrize equilibrium states, but a generic response field,
the LS1 closure seam, and YC20's balanced source are not automatically the
same equilibrium condition. Neither blanket statement about the centre
settles that identification.

## 2. Exact counterexamples to the overstatements

**Silent centre, coupled response.** Let

\[
 U(S,V)=\tfrac12(S^2+2SV+4V^2),\qquad
 H=\begin{pmatrix}1&1\\1&4\end{pmatrix}>0.
\]

The UP1 centre reading U-(S U_S+V U_V)/2 is identically zero, but
det H=3 and CA1 gives Cp/Cv=4/3 wherever the capacities are defined.
Silence of that scalar reading is not zero response or absence of coupling.

**Same axes, unreachable one-step target.** For H=diag(1,2), target
H_star=diag(2,3), the first entry of H+lambda H^2 demands lambda=1,
whereas the second demands lambda=1/4. Allowing an independent affine
spectral map can reach it, but that enlarges the declared operations.

**Flat does not mean an inactive spectral update.** For

\[
 U=e^x+e^y,\qquad
 \widetilde U=U+\tfrac\lambda4(e^{2x}+e^{2y}),
\]

D^2 U_tilde=H+lambda H^2 exactly. For lambda>0 the new response is positive,
flat and separable, yet its eigenvalues change. This satisfies NSB2 T2;
the contradiction is with the pasted interpretation, not with its theorem.

**Repeated spectrum.** If H=I, every target commutes with H, but every
scalar polynomial p(H) is scalar. A nonscalar target is still inaccessible
by that operation.

## 3. A precise target map in a calibrated two-channel response

Work on positive symmetric2x2 response matrices in a declared Euclidean
affine calibration. A different fixed background metric requires first
using its orthonormal representation. The Frobenius pairing below is
chosen data, not an invariant of arbitrary untransported unit changes.

For nonscalar H define

\[
 m=\tfrac12\operatorname{tr}H,\quad X=H-mI,\qquad
 m_* =\tfrac12\operatorname{tr}H_*,\quad X_*=H_*-m_*I,
\]
\[
 \beta={\operatorname{tr}(X X_*)\over\operatorname{tr}(X^2)},
 \quad \alpha=m_*-\beta m,\quad
 \boxed{H_*=\alpha I+\beta H+E_\perp.}                         \tag{1}
\]

Here E_perp is Frobenius-orthogonal to both I and H. It is the part of the
target that no scalar polynomial of this fixed H can supply. Since affine
polynomials span the whole spectral algebra for this nonscalar2x2 H,

\[
 \inf_{p}\|H_*-p(H)\|_F=\|E_\perp\|_F
 = {\|[H,H_*]\|_F\over|h_1-h_2|}.                            \tag{2}
\]

Proof: in an eigenbasis X=diag(r,-r) and E_perp has entries (0,c;c,0).
Its squared norm is2c^2; the commutator's is8r^2c^2, and the eigenvalue
separation is2|r|. Orthogonal projection proves the infimum. The projection
also preserves positivity of H_star, being its diagonal part in this
basis. At H=mI use span{I} directly instead of dividing by the gap.

For the specific single step H+lambda H^2, even E_perp=0 is not enough:

\[
 \beta=1+\lambda\operatorname{tr}H,
 \qquad \alpha=-\lambda\det H.                              \tag{3}
\]

Allowed lambda values and positivity must also be imposed. With iteration,
all steps stay in the original spectral algebra, but their reachable
subset requires the actual iteration law. With a field H(x), (1)-(2) are
pointwise diagnostics; they do not prove a global polynomial or integrable
potential representation of H_star(x).

Rebasing now has a clear meaning: use H-H_star as a comparison residual.
That does not physically replace the Hamiltonian by H-H_star, nor does
shifting a scalar energy origin change its Hessian or its spectral gap.
For a fixed nonscalar H, [H,H-H_star]=-[H,H_star] measures whether this
operator-valued residual brings in a different eigenspace channel. It is
not automatically curvature of a connection.

**Zero specialization.** Construct a declared family (H_lambda,H_star,lambda)
and its residual before evaluation. If H_0 is nonscalar, (1) is regular
near zero and evaluating its projection agrees with constructing that
projection from the evaluated data. If H_0=mI, the spectral algebra drops
to span{I}; the projection from nonzero levels need not have that same
limit. For example H_lambda=I+lambda diag(1,-1) retains a directional
eigenspace off zero which the value H_0 alone forgets. Use the scalar
branch at zero, or explicitly retain the directional jet/labelled cut as
additional input, as LS1/LS2 do. Do not conceal the missing source in a
division by the vanishing eigenvalue separation. Target rebasing alone
does not repair this degeneracy.

## 4. Use force data, but check reach before selecting F

NSB2's potential-preserving operation is

\[
 \widetilde H=F'(U)H+F''(U)gg^T,
 \qquad g=dU,\qquad
 [H,\widetilde H]=F''[(Hg)g^T-g(Hg)^T].                       \tag{4}
\]

If F'' is nonzero and g is not an eigenvector of H, this gives an actual
noncommuting response direction. If g=0 or g is an eigenvector, it does
not. F'>0 and F''>=0 suffice to retain strict convexity, but are not
necessary in every local case. Energy origin and force calibration matter:
adding an affine term to U can leave H unchanged while changing g and
therefore the constructed response.

At a fixed state a proposed target must lie in span{H,gg^T}, with any
chosen admissibility restrictions on a=F' and b=F''. This span need not
contain every symmetric target. Over a neighborhood, a and b must be
functions of U alone and satisfy da=b dU; compatible higher derivatives
and regularity are required. Independent pointwise fits are insufficient.

If additional operations independently provide the I direction, then
{I,H,gg^T} spans the real symmetric2x2 matrices precisely when H is
nonscalar and gg^T does not commute with H. This is a **pointwise algebraic
span** statement, not a theorem of global controllability or a license to
insert an arbitrary source into Yang–Mills.

The proposed next-direction rule can be made explicit. On an actual smooth
admitted response family H(c), let A_j=partial_(c_j)H at the present state,
E=H_star-H, and L=||E||_F^2/2. Define

\[
 G_{ij}=\operatorname{tr}(A_i A_j),\quad
 b_i=\operatorname{tr}(A_iE),\quad
 \dot c=G^+b,\qquad
 {dL\over d\epsilon}\bigg|_0=-\|\operatorname{proj}_{\mathrm{span}\{A_j\}}E\|_F^2.
                                                               \tag{5}
\]

This is the orthogonal projection/least-squares descent rule, using the
Moore--Penrose inverse. It identifies a locally helpful **admitted change**;
it is not physical time or the consequence of merely renaming zero. For
one-sided constraints use the allowed tangent cone instead. If the
projection vanishes, this test supplies no first-order descent; nonlinear
or bracket-generated routes may still exist. No global convergence claim.

A next **measurement** is a separate choice. When NSB2 T4's eigenvalue-field
data are already available, a calibrated vector Hv repairs the missing
two-dimensional orientation information. Without those prior data, one
vector does not determine all response derivatives or curvature. The
binomial(n,3) remainder counts constitutive coordinates, not spacetime
dimensions of the Yang–Mills problem.

## 5. Where the existing YM line sits on this map

| Existing operation | Information it uses | What it can change / next condition |
|---|---|---|
| NSB2 T1 spectral powers; YC30 normalization tau(A+tau)^-1 | One complete operator, scalar function | Eigenvalues/normalization; preserves the input operator's spectral subspaces. YC30 uses an already proved gap. |
| NSB2 T3 potential composition | U, dU, D^2U and chosen F | Can add a noncommuting response channel; calibrate forces, convexity and target compatibility. This is not automatically a YM Hamiltonian update. |
| UP6 cross-corner operations | State-variable charts and conditional derivatives | Ordered derivative defect; a target-specific map is needed to identify it with a matrix or YM operator commutator. |
| YC20/28 exact hidden return | Cross source B and complete hidden inverse (D-z)^-1, with metric | Can mix retained eigenvectors. This already contains information beyond powers of the retained block. Keep all source and metric contributions. |
| YC31/33 moving vacuum frame | Full family derivative and spectral connection | Transports the vacuum split; unitary conjugation alone preserves the physical spectrum. It is not a new gap-generating mechanism. |
| YC32/33 source reduction and boxes | Actual inverse response, source isometry and locality | Controlled one-factor reduced approximation; source selection is extra structure beyond scalar functional calculus. |
| YC34 simultaneous sources | Full excitation supports and vacuum-annihilating local interaction | Fixed-count spatial errors; the outstanding estimate is number/energy control of inverse responses, then box comparison and physical refinement. |

The hidden-return row has an elementary exact witness. Take the positive
matrix and its first-two-coordinate cut

\[
 A=\begin{pmatrix}2&0&1\\0&3&1\\1&1&4\end{pmatrix},\qquad
 A_{PP}=\operatorname{diag}(2,3),\qquad
 A_{\mathrm{eff}}=A_{PP}-{1\over4}
                   \begin{pmatrix}1&1\\1&1\end{pmatrix}.
                                                               \tag{6}
\]

Its principal minors are2,6,19, so A>0.
[A_PP,A_eff]=(0,1/4;-1/4,0) is nonzero, and
A_eff=(P A^-1 P)^-1 on the retained range. Normalization of the full A
is spectral; reduction through a source-selected cut need not be a scalar
function of A_PP. Thus our YM route has **not** remained confined to the
two-dimensional spectral algebra invoked in the pasted conversation.
This finite witness verifies the type distinction; it is not a YM estimate.

For YM the target remains a positive mass gap along a constructed physical
continuum trajectory. It has not been identified with one known2x2 target
matrix or a chosen F. Immediately after YC34, the useful proposed data are
the complete vectors Y_k=A^-1 V_k and their factor-number/energy tails,
not only eigenvalues of a response matrix. A bound on
||(1-P_(<=j))U*Y_k|| would measure their unresolved support, but additional
weighted control is needed to combine that tail with errors growing in j.
NSB2 guides this source-resolved choice; it does not supply that bound.

## 6. What was verified and retained

The source text was compared with the pasted claims. New direct exact
controls check the silent-centre counterexample, the constrained spectral
step, the nontrivial Hessian-preserving update, the target projection and
commutator norm, force-channel rank, local descent and the Schur witness.
No frozen predecessor tests or evidence were regenerated. RMG2's separate
paper and a universal reachability/continuum theorem are not certified.

Retain the owner's useful rule as: **declare the target and pairing, keep
the target residual, and choose the next observation or operation from
verified accessible directions, preserving source and history.** The map
is now explicit. It avoids treating a shifted zero as either a new physical
force or a proof that the target can be reached.
