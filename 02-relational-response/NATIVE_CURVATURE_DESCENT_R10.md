# R10: native curvature conservation and the classical tangent sector

Research owner: **Monty Dabas**. Development and certificate date: **1 October 2026**.

R10 develops the next layer after [R8](CURVATURE_OBSERVATION_R8.md) and [R9](CURVATURE_BALANCE_R9.md): a balance-selected curvature readout obeys a conservation law with a hidden-sector source, and classical Riemann curvature arises in a precisely specified tangent quotient. The same quotient can retain different full native curvatures and finite sheet records.

The result is conditional. The metric and the tangent observer are declared structures whose compatibility is checked; R10 does not derive their physical selection. The algebraic identities below share established connection mathematics. The contribution here is their integration with the existing native cut, typed observation, EMK geometry and retained-memory contracts.

## 1. Types and the meaning of curvature

Let a declared carrier bundle have connection operators

\[
D_i=\partial_i+A_i,\qquad
F_{ij}=\partial_iA_j-\partial_jA_i+[A_i,A_j].
\]

This is the smooth coordinate adapter admitted in R8. Its coordinate directions and its carrier modes are different types. In the primitive finite transport layer, comparison uses typed ordered paths, return operators and an independent integer sheet register. Identifying a finite group move with an exponential of a local connection requires a further integration contract.

For commuting coordinate vector fields, the local two-form has no pair of independent directions in one base dimension. A loop or an independent sheet record can still contain global memory. Adding or removing carrier modes is a different operation from changing the number of base directions. Curvature caused by embedding a curve is yet another geometric object and needs an ambient geometry.

Noncommuting generators alone do not certify this smooth curvature: derivative terms can cancel their commutator, as the R8 pure-gauge witness shows. Arbitrary off-diagonal entries do not certify order dependence either.

Use a constant involutive cut J, with J squared equal to the identity, and define

\[
\alpha(X)=JXJ,\quad
X_e=\tfrac12(X+\alpha X),\quad
X_o=\tfrac12(X-\alpha X),\quad
\Phi(X)=X_e.
\]

For the original two-mode chart, J is the already derived K. Larger carriers below are explicitly admitted examples, not a new derivation of EMK dimensionality.

## 2. The sourced native second Bianchi law

**Theorem R10.1.** For a smooth connection, commuting coordinate partials and a constant cut,

\[
\boxed{
\sum_{\mathrm{cyc}(i,j,k)}
\left(\partial_iF^e_{jk}+[A^e_i,F^e_{jk}]\right)
=-\sum_{\mathrm{cyc}(i,j,k)}[A^o_i,F^o_{jk}].
}
\tag{1}
\]

The visible curvature here is the balanced readout of the full curvature, \(F^e=\Phi(F)\). It is not automatically the curvature recomputed from balanced generators.

**Proof.** Expanding F and using commuting partials and Jacobi gives the established full second Bianchi identity

\[
\sum_{\mathrm{cyc}}(\partial_iF_{jk}+[A_i,F_{jk}])=0.
\]

The constant cut lets the derivative commute with grading. Products respect parity, so

\[
\Phi([A_i,F_{jk}])
=[A^e_i,F^e_{jk}]+[A^o_i,F^o_{jk}].
\]

Apply Phi to the full identity and move the second term to the right. No positivity, metric or physical clock is required. This proof is valid in an associative native algebra admitting the indicated derivations and cut involution.

Equation (1) does not violate ordinary Bianchi. Indeed,

\[
F^e_{ij}=F_{ij}(A_e)+[A^o_i,A^o_j],
\]

and the connection A_e still satisfies its own ordinary Bianchi identity for F(A_e). The extra source accounts for the retained odd-pair curvature in the balanced readout. Removing that contribution changes the object being differentiated.

**Exact nonzero witness.** On an admitted three-mode carrier take

\[
J=\operatorname{diag}(1,1,-1),\qquad
A_1=E_{12},\quad A_2=E_{23},\quad A_3=E_{31},
\]

with constant coefficients. The visible cyclic current and hidden source are

\[
\mathcal J_e=\operatorname{diag}(1,-1,0),\qquad
\mathcal J_o=\operatorname{diag}(-1,1,0).
\]

Each is nonzero; their sum and the full Jacobi current are zero. Thus balance can remove odd readout while its order-sensitive couplings remain necessary for a visible conservation equation.

The two-mode K-even algebra is abelian, so some cyclic source witnesses there vanish identically. The three-mode example demonstrates a nontrivial law without implying three physical spatial dimensions.

The executable certificate includes both constant coefficients and one smooth connection two-jet with nonzero curvature derivatives. Second partials are supplied and their coordinate symmetry is checked. If J varies, derivatives of J contribute additional terms; (1) must then be extended.

## 3. When Riemann curvature is a special case

**Theorem R10.2: connection descent.** On a coordinate neighborhood U, let V carry D and let a smooth surjective bundle map \(C:V\to TU\) intertwine a tangent connection nabla:

\[
\nabla_X(Cv)=C(D_Xv)
\quad\text{for all sections }v\text{ and directions }X.
\tag{2}
\]

Then

\[
\boxed{R^\nabla(X,Y)\,C=C\,F^D(X,Y).}
\tag{3}
\]

**Proof.** Apply (2) twice, subtract the reversed composition and subtract the connection along the vector-field bracket. Each term intertwines. Surjectivity ensures the target is the full declared tangent space rather than only a proper observable subspace.

Equivalently, the kernel of C is preserved by D and the induced quotient connection is the declared tangent connection. This is stronger than a one-time equality of curvature values.

Now admit a smooth nondegenerate symmetric tangent metric g. If the descended connection is metric-compatible and torsion-free **as connection laws on U**, uniqueness of the Levi-Civita connection gives

\[
\nabla=\nabla^{LC}(g),\qquad
\boxed{C F_{ij}=R_{ij}(g)C.}
\tag{4}
\]

For positive g this is Riemannian geometry; an indefinite g gives the corresponding pseudo-Riemannian sector. These established classical results are credited, not presented as a new discovery of Christoffel symbols.

This supplies the requested special-case relation: **classical curvature is the exact tangent quotient of the native connection under (2) and the Levi-Civita contract.** It neither identifies every native carrier with the tangent bundle nor erases native sheet or hidden-sector information.

### The certificate checks the complete local jet

For a constant C, (2) has coordinate form

\[
C A_i=\Gamma_i C.
\]

At a single point this equality is insufficient. The local curvature certificate also checks

\[
C\,\partial_j A_i=(\partial_j\Gamma_i)C
\quad\text{for every }i,j.
\tag{5}
\]

Using (5) for the derivative terms and multiplying the value intertwinings proves (4) at that point.

The implementation constructs Gamma and its first derivatives from an admitted metric two-jet using

\[
\Gamma^a{}_{ib}
=\tfrac12 g^{ac}
(\partial_i g_{cb}+\partial_b g_{ci}-\partial_c g_{ib}),
\qquad
\partial_jg^{-1}=-g^{-1}(\partial_jg)g^{-1}.
\]

It checks that the native jet descends to this entire Levi-Civita jet. It certifies a local two-jet relation; it does not certify an unspecified global metric field or solve its PDEs.

For a variable observer the coordinate law is
\(\partial_i C+\Gamma_iC=CA_i\).
The present executable adapter restricts C to be constant. Constant native and target reframings are checked with \(C'=TCS^{-1}\).

### Two necessary refusal cases

With constant Euclidean g, take A_u and A_v zero at the test point, but \(\partial_v A_u=K\) and all other first derivatives zero. Torsion and nonmetricity vanish at that point, yet

\[
F_{uv}=-K,\qquad R_{uv}(g)=0.
\]

Pointwise torsion and metric values therefore cannot certify a Levi-Civita curvature jet. The derivative laws are needed.

Conversely, the existing R8 noncommuting pure-gauge jet has F equal to zero. It matches the curvature of a flat metric but fails connection intertwining with its Levi-Civita jet. Curvature equality alone does not identify that connection.

## 4. Riemann, distortion and hidden excursions

Choose a declared visible/hidden splitting, with visible tangent dimension n, and write

\[
A_i=
\begin{pmatrix}
\Gamma_i(g)+S_i&E_i\\
L_i&H_i
\end{pmatrix}.
\]

Here E maps hidden modes into visible modes; L maps visible modes into hidden modes. S is the distortion of the visible connection from the admitted Levi-Civita connection.

**Theorem R10.3.** The visible corner of the full native curvature is

\[
\boxed{
(F_{ij})_{\mathrm{vis}}=R_{ij}(g)+\mathcal D_{ij}(S)+\mathcal E_{ij},
}
\tag{6}
\]

where

\[
\begin{aligned}
\mathcal D_{ij}(S)
&=\partial_i S_j-\partial_j S_i
+[\Gamma_i,S_j]+[S_i,\Gamma_j]+[S_i,S_j],\\
\mathcal E_{ij}
&=E_iL_j-E_jL_i.
\end{aligned}
\]

**Proof.** The visible derivative corner differentiates Gamma plus S. Block multiplication gives \((A_iA_j)_{\mathrm{vis}}=(\Gamma_i+S_i)(\Gamma_j+S_j)+E_iL_j\). Subtract the reversed product and expand the visible commutator. This is the connection extension of the R8 discarded-sector identity; its algebraic mechanism is inherited.

For the constant projection \(C=(I_n\;0)\), exact tangent descent requires E equal to zero and the visible connection equal to Gamma throughout the neighborhood. L can remain nonzero. In the local jet certificate their corresponding derivatives must also satisfy these conditions.

Thus visible-to-hidden recording is compatible with an exact classical quotient. Hidden-to-visible feedback generally prevents that quotient. Cancellation in (6) at one point is weaker than descent.

The decomposition is relative to the declared splitting and metric. It does not define a metric-independent scalar or identify a unique physical observer.

### Flat classical geometry with retained native curvature

Take a flat two-mode tangent metric, S and E and L zero, and hidden operators \(H_u=K,\ H_v=R\). Then

\[
R_{uv}=0,\qquad
F_{uv}=
\begin{pmatrix}
0&0\\
0&[K,R]
\end{pmatrix}
=\operatorname{diag}(0_2,-2RK).
\]

The full curvature is nonzero while the exact Riemann quotient is flat. Changing the hidden connection changes full native curvature without changing the classical geometry.

Allowing \(L_u=K,\ L_v=R\) preserves this quotient and produces a nonzero visible-to-hidden curvature block as well. The observer can have an exact tangent geometry while additional native records evolve.

### Hidden feedback produces visible curvature

Keep the visible metric flat and S zero, but set
\(E_u=I,\ E_v=0,\ L_u=0,\ L_v=K\).
Then

\[
(F_{uv})_{\mathrm{vis}}=\mathcal E_{uv}=K,
\]

although R and distortion are zero. The Riemann quotient gate correctly fails. This is a certified order-dependent excursion, not a classification by arbitrary off-diagonal entries.

A separate exact test makes all three terms in (6) nonzero and checks their complete reconstruction.

## 5. Bridge to the existing EMK classical geometry

The existing EMK-G1 paper already declares

\[
g=\operatorname{diag}(W(v),1),\qquad
W(v)=1+\kappa v^2.
\]

R10 uses its positive domain W greater than zero. It does not consume the known EMK-G2 extension across degeneracy or the EMK-G3 compensator-sign issue.

Its metric two-jet yields the source Christoffel operators

\[
\Gamma_u=\begin{pmatrix}0&a\\b&0\end{pmatrix},
\quad
\Gamma_v=\begin{pmatrix}a&0\\0&0\end{pmatrix},
\quad
a=\frac{W'}{2W},\quad b=-\frac{W'}2.
\]

The resulting curvature is

\[
R_{uv}=\begin{pmatrix}0&\mathcal K\\-W\mathcal K&0\end{pmatrix},
\qquad
\mathcal K=-\frac{\kappa}{W^2}.
\]

Lowering the first index gives
\((gR_{uv})_{uv}=W\mathcal K\), exactly the existing source component \(R_{uvuv}\).
Both existing Gaussian-curvature formulas and the source lowered component are compared at 23 rational positive-domain points.

At \(\kappa=3,v=1/3\), the certificate obtains
\(W=4/3,\ \mathcal K=-27/16,\ R_{uvuv}=-9/4\).
An admitted four-mode extension with this visible Levi-Civita jet, hidden K/R operators and nonzero L descends exactly to EMK-G1.

The source already contains metric curvature. The new result is the explicit native-to-classical quotient contract, its refusal tests and its hidden-sector decomposition.

## 6. Finite return and sheet memory beyond the tangent readout

**Theorem R10.4: finite composition.** If typed finite moves satisfy \(CT_i=V_iC\), then

\[
C T_m\cdots T_1=V_m\cdots V_1 C.
\tag{7}
\]

The proof is repeated substitution in the original order. This permits an exact visible return with a nontrivial full return.

Lift the existing finite R/K order loop to a four-mode carrier, acting trivially on the first two modes. Its return is

\[
T_{\mathrm{loop}}=\operatorname{diag}(I_2,-I_2),
\qquad
C T_{\mathrm{loop}}=C.
\]

The visible carrier returns; the hidden carrier does not. Separately, an identity four-mode transport can carry integer sheet residue two. The existing R7 return audit retains that residue until an explicit sheet ledger accounts for it.

These are finite-path certificates. They are not asserted to be the holonomy of the local metric jet without a supplied global path integration.

## 7. What is certified and what remains to select

| Evidence | R10 coverage |
| --- | --- |
| Written general proofs | Constant-cut sourced second Bianchi; bundle descent; block decomposition; finite composition |
| Native symbolic replay | Jacobi, graded source cancellation, distortion commutator expansion, typed curvature descent and typed composition |
| Exact rational tests | 25 tests, including smooth derivatives, source metric comparisons, positive and negative quotient cases, hidden feedback and finite memory |
| Mutation controls | Six incorrect mathematical variants rejected through assertion failures |
| Native alteration controls | Changed input contract and changed rewrite output rejected |
| Earlier evidence | R1–R9 and master-review records and their recorded hashes preserved |
| Further validation | No new Lean formalization, external peer review or physical experiment claimed |

The single canonical operator engine remains in **RKF/operator_foundation**. This repository adds research callers and admitted geometric adapters; it copies no engine implementation.

Reproduce from this repository with Python 3.11 or 3.12, Node and separate checkouts at the commits in the [pins](../04-operator-evolution/R10_SOURCE_PINS.json):

    python3.12 -B 04-operator-evolution/verify_r10.py \
      --publications-root ../Publications \
      --rkf-root ../Recognition-Kernel-Framework

The verifier checks pinned source SHA256 and Git blob hashes before executing them. See the [implementation](../04-operator-evolution/native_curvature_descent.py), [focused tests](../04-operator-evolution/test_native_curvature_descent.py), [verification record](../04-operator-evolution/R10_VERIFICATION.json) and [native proof packet](../04-operator-evolution/R10_NATIVE_CERTIFICATE.json).

The next substantive selection problem is to construct C and g from a specified native response/information sector and test whether its kernel is preserved under the actual transport. If the kernel is not preserved, the excursion term in (6) gives an explicit correction to the visible classical geometry. No universal metric, physical dimension, rate or flatness condition is selected by merely declaring the adapter.

## Source and result lineage

- [R7 typed tensor calculus](EMK_TENSOR_CALCULUS_R7.md), [R8 observation](CURVATURE_OBSERVATION_R8.md), [R9 balance](CURVATURE_BALANCE_R9.md) and the [master review](EMK_MASTER_TENSOR_REVIEW.md) are retained unchanged.
- [EMK recognition geometry](https://github.com/Parveen117/Publications/tree/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-recognition-geometry) supplies the existing metric family and classical Christoffel/Riemann formulas. Its unchanged EMK-G1 runtime is imported directly.
- [Native algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) supplies typed composition, observer-kernel/ideal conditions and the canonical proof engine.
- [Existing cut-graded curvature theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/48_emk_algebra_cut_graded_curvature_theorem.md) supplies native parity and curvature lineage.
- Jacobi, ordinary Bianchi, uniqueness of Levi-Civita and connection distortion are established mathematics. R10's certified development connects them to these native observation and memory contracts; it does not claim their classical discovery or universal physical novelty.
