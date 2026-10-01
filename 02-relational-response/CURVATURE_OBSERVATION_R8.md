# R8: Curvature, observation and retained sheet memory

R8 makes the proposed connection precise on the existing EMK KIR algebra: ordered mismatch, projected observation, cut sign, dimension restriction and response flatness are related by explicit maps. The main result is an exact correction law between **compressing full curvature** and **recomputing curvature after compression**. Either full or reduced curvature can vanish while the other is nonzero.

The general claims below have written algebraic proofs. The executable certificate supplies exact rational examples, counterexamples and proof replays against the unchanged canonical engine. This is a scoped research extension of the existing master, not an alternative engine or a universal physical law.

## 1. Native foundation and comparison contract

The active operator engine remains [RKF/operator_foundation](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation). Its cut scalars, typed actions, lawful composition and residue policy precede coordinate representations. Its existing EMK presentation derives

\[
K^2=I,\qquad R^2=-I,\qquad KR=-RK.
\]

The native adapter in R8 uses that presentation unchanged; the engine derives a four-word basis and certifies equalities before any two-mode interpretation. For exact coordinate examples we reuse R7's admitted rational chart

\[
K=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
RK=\begin{pmatrix}-1&0\\0&1\end{pmatrix}.
\]

Let two actions be lawfully comparable endomorphisms of the **same** carrier. Their algebraic ordered mismatch is \(F=[A,B]=AB-BA\). For differentiable coordinate comparison directions, the admitted connection convention is \(D_i=\partial_i+A_i\), so

\[
F_{ij}=\partial_iA_j-\partial_jA_i+[A_i,A_j].
\]

Only constant coefficients reduce this to a commutator alone. A matrix product, a scalar response form, a finite loop holonomy and an integer sheet residue remain separately typed quantities. Connecting them requires the stated adapter; a determinant ratio or trace alone does not identify the full mismatch.

Three dimensions must be distinguished: the number of base comparison directions, the dimension of a state representation, and the dimension of an algebra coefficient space. Here the native algebra has four coefficient directions, its admitted state chart has two modes, and the base may have one, two or three comparison directions. None equals an integer sheet label by definition.

## 2. Theorem R8.1: exact observation correction

In any unital associative algebra let \(P^2=P\), \(Q=I-P\), and define

\[
C_P=[PAP,PBP],\qquad
E_P=PAQBP-PBQAP.
\]

Then

\[
\boxed{C_P=P[A,B]P-E_P.}
\]

**Proof.** Insert \(I=P+Q\) between the two factors of \(PABP\) and \(PBAP\):
\(PABP=PAPBP+PAQBP\), and \(PBAP=PBPAP+PBQAP\). Subtract and rearrange. This is the antisymmetric specialization of the already established native compression defect, N06. It needs idempotence and lawful composition, not a metric, positivity or orthogonal projection. □

The native packet also certifies this identity in the free associative presentation on `A,B,P` with the sole rule `PP -> P`. The unchanged canonical rewrite engine checks its overlap/confluence audit and replays reduction of the identity's difference to zero. This symbolic certificate assumes no finite carrier and is separate from the finite matrix examples and from a Lean formalization.

Consequently \(C_P=0\) exactly when \(PFP=E_P\). This is not an equivalence with \(F=0\). The antisymmetrized excursion term records compositions that leave the retained sector between two actions.

For a constant projector applied to a connection, the derivative terms compress linearly and the same result becomes

\[
F^{\rm reduced}_{ij}=PF_{ij}P-
\big(PA_iQA_jP-PA_jQA_iP\big).
\]

A varying projector requires its derivative terms and is outside this certificate. Compression is a declared intervention between actions; it is not automatically the quotient action of an exact observer. In particular, an invariant discarded subspace or a multiplicative observation contract can force the excursion terms to vanish.

**Full curved, reduced flat.** Take \(A=K\), \(B=R\), \(P=\operatorname{diag}(1,0)\). Then

\[
F=\operatorname{diag}(2,-2),\quad PFP=2P,\quad E_P=2P,\quad C_P=0.
\]

**Full flat, reduced curved.** Take a three-mode admitted carrier,

\[
A=\operatorname{diag}(1,0,0),\quad B=\operatorname{diag}(0,1,0),\quad
P=I-\tfrac13\mathbf1\mathbf1^T.
\]

Here \(P^2=P\), \(\operatorname{rank}P=2\), and \([A,B]=0\), but

\[
C_P=\frac19\begin{pmatrix}0&-1&1\\1&0&-1\\-1&1&0\end{pmatrix}\ne0,
\qquad E_P=-C_P.
\]

Thus reduced ordered mismatch can also be introduced by the observation intervention. The transpose in this example specifies one concrete projector; it is not an assumed native inner product.

If \(A,B,P\) are all conjugated by the same invertible constant \(S\), each of \(F,PFP,C_P,E_P\) conjugates by \(S\). This follows by cancelling adjacent \(S^{-1}S\) factors and proves reference covariance of the whole identity.

## 3. Theorem R8.2: exactly two static odd-curvature channels

The existing cut-graded theorem gives, for \(G_e=aI+bK\) and \(G_o=cR+dRK\),

\[
[G_e,G_o]=-2b(dR+cRK).
\]

These KIR relations and this curvature formula are prior results. R8 specifies an observation and recovery contract on their odd target space. Write \(F=xR+yRK\), and introduce dual vectors without a norm assumption:

\[
w_+=(1,1)^T,\quad w_-=(1,-1)^T,\quad
\ell_+=\tfrac12(1,1),\quad\ell_-=\tfrac12(1,-1).
\]

They satisfy \(\ell_\alpha w_\beta=\delta_{\alpha\beta}\) and \(P_\pm=w_\pm\ell_\pm=(I\pm K)/2\). Define the oriented cross-sector channels

\[
q_+=\ell_+Fw_-=x-y,\qquad q_-=\ell_-Fw_+=-x-y.
\]

Then

\[
\boxed{x=(q_+-q_-)/2,\qquad y=-(q_++q_-)/2.}
\]

**Proof and minimality.** The readout matrix is \(\begin{pmatrix}1&-1\\-1&-1\end{pmatrix}\), with determinant \(-2\). It is invertible over the declared characteristic-zero field. Every single fixed linear scalar channel on a two-dimensional target has a nontrivial kernel, so two are necessary as well as sufficient. For example \(F=R+RK\ne0\) has \(q_+=0\). This minimality is for fixed **linear** channels on this target; arbitrary nonlinear or extra-prior protocols are not included. □

For any odd \(F\), \(P_\pm F P_\pm=0\), because \(KFK=-F\). Its cross-sector blocks need not vanish. Also \(\operatorname{tr}[A,B]=0\) for finite matrices, by cyclicity, including the nonzero \([K,R]=\operatorname{diag}(2,-2)\). The full trace and same-parity diagonal blocks are blind to this curvature. The two channels can instead be written as weighted traces, \(q_+=\operatorname{tr}(w_-\ell_+F)\) and \(q_-=\operatorname{tr}(w_+\ell_-F)\), with the weights explicitly declared.

An **active** cut action \(F\mapsto KFK=-F\), with probes held fixed, reverses both channel signs. A **passive** constant reference change \(F\mapsto SFS^{-1}\), with \(w\mapsto Sw\) and \(\ell\mapsto\ell S^{-1}\), leaves each pairing unchanged. A coordinate orientation reversal can also change a two-form component's sign without changing whether the full tensor vanishes. These are different operations.

The channels recover \(F\), not its generating factors: \((b,c,d)=(1,3,2)\) and \((2,3/2,1)\) produce the same curvature.

### Static recovery versus future-complete observation

Two channels suffice for the current odd target. Under the broader contract of unrestricted algebra-valued states and arbitrary further lawful native left actions, these same fixed rows are not future complete. For example their current values vanish on both \(I\) and \(0\), but a subsequent left action by \(R\) distinguishes these algebra states. The unchanged engine's N03 row-closure certificate expands these two rows on the discovered four-word basis to rank four, requiring **two extra channels**. This is an application of the existing minimum observer theorem, not a new theorem claiming every EMK observer has dimension four. An initially known odd target can instead be reconstructed and propagated with the full action history; that is a different observation contract.

The proof packet is [R8_NATIVE_CERTIFICATE.json](../04-operator-evolution/R8_NATIVE_CERTIFICATE.json). Besides the general compression identity, it replays \([K,R]=-2RK\), \(K[K,R]K=2RK\), and the rank-2 to rank-4 future observer completion against externally pinned input hashes.

## 4. Theorem R8.3: dimension restriction by pullback

For \(d\) independent coordinate comparison directions, an operator-valued curvature two-form has \(d(d-1)/2\) labelled component slots. This is a count of antisymmetric slots, not a claim that all components are nonzero or independent under extra equations.

For a constant affine base map \(x=Jy\), where \(J\) has \(d\) rows and \(m\) columns, the pulled-back connection jet is

\[
A'_a=\sum_iJ_{ia}A_i,\qquad
\partial'_b A'_a=\sum_{i,k}J_{ia}J_{kb}\,\partial_kA_i.
\]

Its curvature is

\[
\boxed{F'_{ab}=\sum_{i<j}(J_{ia}J_{jb}-J_{ja}J_{ib})F_{ij}.}
\]

**Proof.** Substitute the displayed values and derivatives into the connection formula. Bilinearity gives the double sum over \(i,j\); antisymmetry combines the two orderings into each displayed minor. □

An invertible base coordinate change preserves full local flatness, because its inverse pulls back to the original tensor. Restriction to one comparison direction has no local two-form component even when the original connection is curved. That restriction loses transverse comparison information; it does not prove full flatness. With \((A_r,A_{t1},A_{t2})=(K,R,RK)\), the existing relations give

\[
F_{r,t1}=-2RK,\qquad F_{r,t2}=-2R,\qquad F_{t1,t2}=-2K.
\]

The dimension-dependent result concerns available comparisons. It does not identify this two-form with extrinsic bending of a line or circle. The executable `ConnectionJet` stores values and first derivatives at a point, validates dimensions, and checks this constant-map rule. It does not infer a global field from finitely many jets.

## 5. Theorem R8.4: Onsager flatness in the declared scalar response sector

For constant linear response \(J=LX\), take the scalar one-form \(\omega=\sum_iJ_i\,dX_i\) on the independent force coordinates. Equivalently the connection coefficients are the one-by-one matrices \(A_i=J_i\). They commute, so

\[
\boxed{F_{ij}=\partial_iJ_j-\partial_jJ_i=L_{ji}-L_{ij}.}
\]

**Proof.** Differentiate \(J_i=\sum_kL_{ik}X_k\) and use the zero scalar commutator. Thus all \(F_{ij}\) vanish if and only if \(L=L^T\). □

This is the ordinary Onsager-symmetric condition for this constant response sector; selecting it as a physical model requires its physical hypotheses. It does not cover every Onsager-Casimir parity convention or imply that all components of an EMK master tensor vanish.

Two exact controls delimit the claim:

- For state-dependent \(L(x,y)=\operatorname{diag}(y,0)\), the matrix is symmetric at every point but \(\omega=xy\,dx\) has \(F_{xy}=-x\). Derivative compatibility matters.
- Noncommuting connection coefficients can nevertheless give zero total curvature. Let \(g=(I+sE_{12})(I+tE_{21})\), with determinant one. Direct differentiation gives
  \[
  A_s=g^{-1}\partial_sg=\begin{pmatrix}t&1\\-t^2&-t\end{pmatrix},\quad
  A_t=g^{-1}\partial_tg=\begin{pmatrix}0&0\\1&0\end{pmatrix}.
  \]
  Here \([A_s,A_t]=\partial_tA_s=\begin{pmatrix}1&0\\-2t&-1\end{pmatrix}\ne0\), but \(F_{st}=-\partial_tA_s+[A_s,A_t]=0\) identically. Conversely commuting scalar coefficients with asymmetric constant \(L\) have nonzero curvature.

Noncommutativity therefore diagnoses the constant-coefficient algebraic sector, not total connection curvature without its derivative terms.

## 6. Global sheet memory remains an independent target

R8 reuses R7's explicitly declared direct-product integer sheet sector. A closed move with identity carrier and winding two returns every tested local vector readout while retaining sheet residue two. A ledger value of two, fixed as part of the closure policy, removes that residue; an undeclared modulo-two alias does not.

The scope of that witness is the R7 product model. It does not derive a universal sheet law from local curvature, nor replace the existing EMK deck action. The previously identified EMKG3 compensator sign inconsistency is not consumed by this certificate. No upstream repair is applied here.

## 7. What is added and what already existed

| Source result | R8 use or addition |
|---|---|
| RKF [T48](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/48_emk_algebra_cut_graded_curvature_theorem.md): KIR grading, explicit curvature and cut sign | Reused unchanged; adds the two-channel static recovery contract and counterexamples to blind readouts |
| Canonical [N03–N06](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md): observer closure and compression defect | Antisymmetrizes the defect, certifies both directions of apparent curvature, and specializes minimum future observer repair |
| RKF [noncommutative holonomy](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/thermodynamics/lambda_geometry/11_noncommutative_holonomy_and_irreducibility.md) | Retains order sensitivity; does not equate every finite loop with an infinitesimal commutator |
| RKF [Onsager compass](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/thermodynamics/05_onsager_compass_and_constitutive_no_go.md) | Reuses response integrability; makes the constant-coefficient orientation and variable-response failure explicit |
| R7 tensor/transport and integer ledger; existing EMK master review | Reuses their chart and sheet witness unchanged; supplies a linked observation/curvature certificate |

The algebraic projection identity, linear rank argument, connection formula and pullback calculus are established mathematics. The contribution is this explicit EMK observation adapter, its counterexample suite, native-engine replay and linked claim boundary. Private Vault bodies and upstream engine files are not copied into this repository.

## 8. Reproduction and evidence levels

Use Python 3.11 or 3.12 and Node, with a separate RKF checkout at commit `3cc5a33b05c16d59c90994ddda69dedc0d392424`. From the `extra-ideas` root:

```bash
python3.12 -B 04-operator-evolution/verify_r8.py --rkf-root ../Recognition-Kernel-Framework
```

The verifier checks every listed upstream SHA256 and Git blob before upstream execution. It then runs 23 exact tests, including 100 compression cases, 25 odd-target reconstructions, 27 unchanged source comparisons and 81 constant-response cases. Six deliberate mathematical mutations must produce assertion failures; three altered native contracts/results must be rejected. The general rewrite proof and complete finite native result are recomputed under externally pinned input hashes. It checks recorded R1–R7 and master-review hashes without rerunning their unchanged suites.

[R8_SOURCE_PINS.json](../04-operator-evolution/R8_SOURCE_PINS.json) binds the upstream sources and job contract. [R8_VERIFICATION.json](../04-operator-evolution/R8_VERIFICATION.json) binds the proof, implementation, tests, adapter, native packet and preserved evidence. Its finite PASS is distinct from the general written proofs, a Lean formalization, independent peer review, global analytic existence or physical validation. No new Lean theorem or physical experiment is claimed.
