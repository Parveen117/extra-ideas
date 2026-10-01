# R28: native cut curvature, hidden response and propagation

Research owner: Monty Dabas. Development: 1 October 2026.

The useful physics route is to connect the existing native curvature to an
actual response law. R28 identifies R18's visible–memory coupling with the
source H/K commutator, eliminates the hidden channel exactly, and propagates
a signed curvature response with the already derived native wave identity.
It also proves a chart-free cyclic curvature ledger directly by cancellation
of words. No Maxwell equation, vacuum constitutive law, smooth connection,
spacetime dimension, physical metric or external clock enters these proofs.

This is a mathematical bridge toward a **vacuum-as-cut proposal**. It does
not identify a cut with physical vacuum or its two roles with electric and
magnetic fields. In particular, zero recognized output, zero operator,
zero external injection and an uncut ground are different statements.

## Source and target

R16 C3–C8 constructs signed/refined coefficients, the matching pairing,
the positive square root, and the two role maps

\[
H^2=K^2=I,\quad HK=-KH,\quad H^\dagger=H,\quad K^\dagger=K,
\qquad R=KH,\quad R^2=-I.
\]

R18 constructs the address of a retained source word from its successive
role contrasts. Let \((Sf)(x)=f(x-1)\),

\[
P=(I+H)/2,\quad Q=(I-H)/2,\quad
D=SP+S^{-1}Q,\quad C=(H+K)/\sqrt2,\quad U=DC.
\tag{28.0}
\]

S is an address shift, not a complex phase. P and Q are source-role cuts.
R18's signed aggregation at a common address and its finite sum of matching
pairings remain explicit response targets. Its normalization is derived
within that target; the target is not declared to be nature's unique readout.
We retain original word tags as provenance. These finite role calculations
do not replace the full generalized UGD source by ordinary complex numbers.

## Written results

### R28.1 — A zero cut reading can hide a nonzero curvature square

Define the signed order residue and the cut-even readout by

\[
\mathcal F=HK-KH=-2R,\qquad
\Phi_H(X)=(X+HXH)/2=PXP+QXQ.
\]

Then

\[
\boxed{\Phi_H(\mathcal F)=0,\qquad
\mathcal F^\dagger\mathcal F=4I,\qquad
\Phi_H(\mathcal F^\dagger\mathcal F)=4I.}
\tag{28.1}
\]

The actual four-arrow loop is also nontrivial:

\[
HKH^{-1}K^{-1}=-I,\qquad HKH^{-1}K^{-1}-I=-2I.
\tag{28.2}
\]

Thus a cut can annihilate the signed order reading while retaining its
quadratic response and the complete loop displacement. These are different
targets. A zero first reading is not a flatness theorem.

**Proof.** Anticommutation gives \(\mathcal F=-2KH\) and
\(H\mathcal F H=-\mathcal F\). Dagger reverses the two letters, so
\(\mathcal F^\dagger=-\mathcal F\) and its square is \(-4I\).
The displayed identities follow by multiplication. Since H and K are their
own inverses, the loop is \(HKHK=-I\). More generally its displacement is
\((HK-KH)H^{-1}K^{-1}\), explaining why it is not the same fixed-input
response as \(\mathcal F\). No homology or smooth small-loop limit is used.

### R28.2 — The native propagation coupling is exactly cut curvature

The finite transition differential of the role cut obeys

\[
\boxed{\nabla_U P:=PU-UP
       =D\mathcal F/(2\sqrt2),\qquad
(\nabla_U P)^\dagger\nabla_U P=I/2.}
\tag{28.3}
\]

Write the four blocks of U as \(A=PUP, B=PUQ, L=QUP, M=QUQ\).
As maps between their indicated native roles,

\[
A=SP/\sqrt2,\quad B=SPKQ/\sqrt2,\quad
L=S^{-1}QKP/\sqrt2,\quad M=-S^{-1}Q/\sqrt2.
\tag{28.4}
\]

In particular \(B^\dagger B=Q/2\), \(L^\dagger L=P/2\).
The full pairing splits as \(\|\psi\|^2=\|P\psi\|^2+\|Q\psi\|^2\)
and U preserves it. The transfer of a pure hidden or pure visible input
has exactly half its initial response energy in each output role. For a
mixed input the two contributions interfere; this is not a classical
two-state probability rule.

**Proof.** D commutes with P and is an isometry. Hence
\([P,U]=D[P,C]=D[H,K]/(2\sqrt2)\). Apply (28.1) to its dagger square.
Multiplication by the disjoint projectors gives (28.4), using
\(PK=KQ\). Squaring these blocks gives the stated halves. Orthogonality
of P and Q proves the energy split; \(C^\dagger C=I\) and
\(D^\dagger D=I\) prove conservation. No additional coupling coefficient
has been fitted or installed.

### R28.3 — Eliminating the hidden cut derives an exact memory kernel

For \(\psi_{n+1}=U\psi_n\), put \(a_n=P\psi_n\), \(b_n=Q\psi_n\).
For every finite n,

\[
\boxed{a_{n+1}=Aa_n+
 \sum_{k=0}^{n-1}BM^kL\,a_{n-1-k}+BM^n b_0,}
\tag{28.5}
\]

where the source determines every kernel coefficient:

\[
BM^kL=(-1)^k2^{-(k+2)/2}S^{-k}P,\qquad
BM^n b_0=(-1)^n2^{-(n+1)/2}S^{1-n}K b_0.
\tag{28.6}
\]

The formal event-lag series is therefore

\[
\sum_{k\ge0}z^k BM^kL
 =\tfrac12(I+zS^{-1}/\sqrt2)^{-1}P.
\tag{28.7}
\]

Here z records an event lag. The inverse is a formal series with constant
coefficient I; no analytic convergence, logarithm, physical frequency or
external time is asserted. Equation (28.5) is finite and needs no series.

A concrete zero-reading witness is
\(\psi_0=\delta_0\mathcal F e_0/2=-\delta_0e_1\):

\[
a_0=0,\quad b_0=-\delta_0e_1,\quad
\psi_1=(-\delta_{1}e_0+\delta_{-1}e_1)/\sqrt2.
\tag{28.8}
\]

The visible response becomes nonzero without a subsequent injected term.
Its energy comes from the retained hidden response: it is not creation from
the zero vector. A rule depending only on the current visible value cannot
describe both this initial state and the actual zero state.

**Proof.** Iterate \(b_{n+1}=La_n+Mb_n\) to obtain
\(b_n=M^n b_0+\sum_{k=0}^{n-1}M^kL a_{n-1-k}\), then substitute in
\(a_{n+1}=Aa_n+Bb_n\). From (28.4), \(M^k=(-1)^k2^{-k/2}S^{-k}Q\)
on the hidden role, and \(PKQKP=P\), proving (28.6). Multiplying the
formal series by \(I+zS^{-1}/\sqrt2\) cancels all positive powers and
leaves \(P/2\). Finally \(\mathcal F e_0/2=-e_1\); one application of
U gives (28.8). The zero vector instead stays zero by linearity.

### R28.4 — A signed curvature response propagates by the native wave law

For any finite native input field v, compare the two source protocols with
the same continuation:

\[
\chi_n=U^n(HK-KH)v/2.
\tag{28.9}
\]

This is the signed difference of two ordered preparations, not an assumed
classical curvature tensor field. It has conserved pairing
\(\|\chi_n\|^2=\|v\|^2\). Its support moves at most one constructed
address per source event. Every component, including P and Q readouts,
obeys for n at least two

\[
\boxed{\chi_{n+2}-2\chi_n+\chi_{n-2}
 =\tfrac12(S^2-2I+S^{-2})\chi_n.}
\tag{28.10}
\]

The same identity is inherited from R18.5; its novelty here is the explicit
curvature preparation and the coupling/memory identification, not discovery
of the existing native wave equation. Its coefficient is a source count
normalization, not a measured c in physical units.

The positive response intensity \(\rho_n(x)=\|\chi_n(x)\|^2\) generally
does **not** obey (28.10). For (28.8), at x=0 one has
\(\rho_0=1,\rho_2=1/2,\rho_4=1/8\), while
\(\rho_2(-2)=\rho_2(2)=1/4\). The two sides of the proposed intensity
wave law at n=2 are \(1/8\) and \(-1/4\).

**Proof.** By (28.1), \(\mathcal F/2\) preserves pairing. U preserves it
as well. Its two blocks shift by precisely plus/minus one; induction proves
the finite support statement. For an independent algebraic derivation write

\[
W=\sqrt2U=
\begin{pmatrix}S&S\\S^{-1}&-S^{-1}\end{pmatrix}.
\]

Direct multiplication gives \(W^2-(S-S^{-1})W-2I=0\). Substituting U
and squaring this identity gives
\(U^4-\tfrac12(S^2+2I+S^{-2})U^2+I=0\).
Apply it to \(\chi_{n-2}\) and rearrange. This derivation uses finite
operator products, without a supplied differential wave equation. The
intensity values follow by four literal applications of W with the
normalization \(2^{-n}\) on squared responses. They reject the nonlinear
replacement of a signed response by its information/energy size.

### R28.5 — Native composition gives a chart-free cyclic curvature ledger

Let \(X_a\) be endomorphisms of one finite native role target; no coordinate
axes or derivatives are presumed. Put \(F_{ab}=[X_a,X_b]\). Associativity
alone gives

\[
\boxed{[X_a,F_{bc}]+[X_b,F_{ca}]+[X_c,F_{ab}]=0.}
\tag{28.11}
\]

For a native involution J define
\(X^{e/o}=(X\mathbin{\pm}JXJ)/2\). Its cut-even reading obeys

\[
\boxed{\sum_{\rm cyc}[X_a^e,F_{bc}^e]
       =-\sum_{\rm cyc}[X_a^o,F_{bc}^o].}
\tag{28.12}
\]

In particular \(F_{ab}^e=[X_a^e,X_b^e]+[X_a^o,X_b^o]\), which is
not in general the curvature computed from the visible arrows alone.
The hidden term in (28.12) can be nonzero.

An exact witness uses two native role copies ordered 00,01,10,11.
With \(E_{01}=PK\), \(E_{10}=QK\), construct

\[
J=H\otimes I,\quad X_1=P\otimes E_{01},\quad
X_2=E_{01}\otimes E_{10},\quad X_3=E_{10}\otimes P.
\]

The left side of (28.12) is \(P\otimes H\ne0\) and its hidden
partner is \(-P\otimes H\). These are finite native tensor operators,
not three physical directions. The X arrows need not be isometries.

**Proof.** Expand the three nested commutators into twelve ordered words.
Each word occurs once with each sign, proving (28.11). Conjugation by J
preserves products, so even-even and odd-odd products are even, whereas
mixed products are odd. Applying the even projection to (28.11) proves
(28.12) and the curvature decomposition. For the witness, on its first
three roles the arrows are \(E_{01},E_{12},E_{20}\). Their three nested
commutators are \(E_{00}-E_{11}\), \(E_{11}-E_{22}\) and
\(E_{22}-E_{00}\). Only the first is an even-arrow/even-curvature term;
the other two supply its negative. The construction with P,Q,K derives
these matrix units from source roles. This is the finite algebraic cyclic
identity, not a claim to have derived spacetime Maxwell dynamics.

### R28.6 — Curvature alone does not select a unique dynamics or a physical vacuum

The same H/K source, nonzero \(\mathcal F\), native pairing and role cut
admit at least three explicitly constructible continuations: identity at
each address, local R at each address, and the R18 address transport U.
All preserve the pairing. The first two never spread a one-address input;
the third sends the witness (28.8) to two different addresses. Thus source
algebra and conservation alone do not select the transport target.

Likewise \(\Phi_H(\mathcal F)=0\) is an operator readout, while
\(P\psi=0\) is a state readout. Neither means \(\psi=0\), nor specifies
absence of external physical sources. The actual zero vector has zero
response under every linear continuation above.

**Proof.** I and R are products of source maps, with \(I^\dagger I=I\)
and \(R^\dagger R=I\). Their action is address-local. R28.2 proves the
same conservation for U, and (28.8) distinguishes its support. R28.1 and
R28.3 supply the two different zero-reading witnesses. Linearity proves
the zero-state statement. These counterexamples locate a real selection
obligation without importing a replacement physical axiom.

## Certification, attribution and physical target

Six written arguments are bound to exact native operator checks, native
symbolic replay, source hashes and a declared proof dependency graph. The
graph gate rejects imported/admitted/open premises and promotion of a target
definition to physical selection. This is not automatic semantic checking
of arbitrary prose or proof-assistant certification. Finite checks support
the written proofs; physical identification needs separate evidence.

| Pinned source | Premise status in this continuation |
| --- | --- |
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native-derived proof premise:** signed counts, cut roles, H/K, matching pairing and positive square-root completion. Finite equality, counting and induction are the stated metamathematical infrastructure. |
| [R18.1–R18.5](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | **Native-derived on an explicit response target:** word address, signed aggregation, normalization, conservation and wave identity. Physical choice of that target is not supplied. |
| [RKF 24](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/24_clock_free_recognition_seam_cut_calculus.md) | **Attribution, not a proof premise:** prior finite transition differential and composition residues on declared event lifts/projectors/transports. R28 reconstructs its required blocks on the source role target. |
| [RKF 48 and 51](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/48_emk_algebra_cut_graded_curvature_theorem.md) | **Attribution, not a new certificate replay:** existing cut-graded commutator and [native odd-channel exchange](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/51_odd_channel_exchange_law_theorem.md). The latter's native scalar-turn grading is distinct from R28's role-cut grading. |
| [MP Gold 07](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/mp_gold/07_zero_as_cut_topology.md) | **Model-scoped comparison only:** zero-as-cut/join defect and retained memory. Its circle/torus homology ledger is declared; numerical critical-response evidence has a user-local lineage qualifier. Neither topology nor a physical vacuum identification is imported here. |
| [R10](https://github.com/Parveen117/extra-ideas/blob/6035e66bd40fe4016eeded41a4d704ec4ff4b342/02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) | **Admitted classical/smooth adapter, comparison only:** tangent descent, smooth connection, metric and torsion/compatibility conditions give the Riemann special case. R28.5 proves its finite operator cyclic counterpart without those inputs. |
| [U30–U32 response tensor and cut ledger](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md) | **Declared-carrier construction, comparison only:** prior oriented response/Gram/seam distinction. U32 separately supplies classical ensemble weights; those weights and the ordinary complex Gram chart are not imported into R28. |
| [Canonical native engine](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native equality/replay premise:** unchanged canonical engine, scoped checks only. No whole-engine or formal-proof-assistant PASS is asserted. |

The next physical task is specific: derive a native field/readout whose
constraints and coupled propagation identify electric and magnetic response,
and determine whether a source-free physical state implements this cut law.
That requires a justified directional/loop structure, a distinction between
field and source, and an operational reading of charge and propagation.
The two role components here do not by themselves establish transverse
polarizations or electromagnetic gauge structure. No Maxwell law is used
as an input, and no EM-wave, physical c, alpha, mass or vacuum identification
is certified by this packet. Native cut curvature, its hidden memory kernel
and its signed propagating response are the completed mathematical bridge.
