# R27: native source-law selection of memory and copy-parity response

Research owner: Monty Dabas. Development: 1 October 2026.

R27 closes a specific freedom left by R26: a memory that continues under
the original native source law cannot have an independently chosen
commutator or coupling angle. Its two involutions and the source's balanced
normalization force the H/K relation. Every finite nonzero record carrying
that law is the same native two-role action with additional inactive
multiplicity. R26's gap, return coefficient and exact channel ratio then
hold throughout this class.

The distinction between **carrying the source law on a record** and
**copying the original source letters onto several records** also becomes
decisive. The latter has an exact parity obstruction and a minimum native
compensator. Counting the retained copies yields two different response
classes, rather than an adjustable continuous noise parameter.

This is selection of a mathematical source-law target. Physical selection
of a particular target, number of retained copies, control rule, clock or
mass/charge interpretation is not supplied by renaming that target.

## Native inputs and the selection target

Use the signed/refined count field, native matching pairing and H/K role
operators of R16 C3–C8. R18 supplies the constructed address shift and
balanced source transport. R20 supplies finite tuple records. R24–R26
supply the already derived return and feedback results. Original word tags
remain retained as provenance, even when their evaluated operators agree.

The main proofs concern finite real signed native role targets. The scalar
repair comparison in R27.5 explicitly uses the already derived native iota;
it does not replace the generalized UGD carrier with a complex primitive.

Write the matching pairing as \(\langle\ ,\ \rangle\). On a nonzero finite
record target \(\mathcal M\), let A and B be pairing-preserving
endomorphisms carrying the source's individual self-return laws
\(A^2=B^2=I\). Their inverses equal their daggers, so \(A^\dagger=A\)
and \(B^\dagger=B\). The source's equal two-tag coefficient census gives

\[
C_{\mathcal M}=(A+B)/\sqrt2,\qquad
G=(AB+BA)/2.
\tag{27.0}
\]

Here G is calculated from the record arrows. It is not an added interaction
parameter. A **source-law record** means that this same census preserves
the native pairing on the record. This target is inherited from the proved
source operations and normalization. R27.1 derives its full mixed relation;
it is not placed into the candidate record as a separate axiom. A generic
direction-controlled pairing-preserving transport need not realize this
additional source-continuation target, as R27.7 demonstrates.

On the moving source keep
\(R_s=K_sH_s\), \(F_s=H_s+K_s\),
\(\Pi_a=(I+(-1)^aH_s)/2\), \(L_a=\Pi_aF_s\).
The record-controlled transport is R26's

\[
V_{A,B}=\bigl(S\Pi_0\otimes A+S^{-1}\Pi_1\otimes B\bigr)
 (F_s\otimes I)/\sqrt2.
\tag{27.1}
\]

This control uses the outgoing direction, not the original H/K letter tag.

## Written results

### R27.1 — Source normalization forces the mixed law and measures its failure

The calculated defect G is self-adjoint and commutes with A and B. It obeys

\[
C_{\mathcal M}^\dagger C_{\mathcal M}=I+G,\qquad
I\pm G=(A\pm B)^\dagger(A\pm B)/2,
\tag{27.2}
\]

\[
(AB-BA)^\dagger(AB-BA)=4(I-G^2).
\tag{27.3}
\]

Consequently the following are equivalent:

1. The original balanced source census preserves the record pairing.
2. \(G=0\).
3. \(AB=-BA\).
4. \(H\mapsto A,\ K\mapsto B\) extends to a unital coefficient-linear,
   multiplicative, dagger-preserving action of the native role algebra.

For the R26 two-event return with source e0 and nonzero record \(\mu\),

\[
\boxed{j_2(0)+\tfrac12
       =\frac{\|G\mu\|^2}{\|\mu\|^2}.}
\tag{27.4}
\]

Thus the return current detects the square of the actual record
normalization defect. Its value is -1/2 for every record preparation
exactly when the record carries the source law. No maximum-curvature or
minimum-energy principle is postulated.

**Proof.** Expand the square in (27.2), using the two involutions.
For example \(AG=(B+ABA)/2=GA\); the other commutation and dagger
identities follow by the same finite multiplication. Put Z=AB, so
\(Z^\dagger=BA=Z^{-1}\). Then
\[
(Z-Z^\dagger)^\dagger(Z-Z^\dagger)
=2I-Z^2-(Z^\dagger)^2=4I-(Z+Z^\dagger)^2,
\]
which proves (27.3). Positivity in (27.2) also gives
\(-\|\mu\|^2\le\langle\mu,G\mu\rangle\le\|\mu\|^2\).

Pairing preservation of the census is its Gram identity, hence equivalent
to G=0. The source's finite relations reduce every word to a signed member
of I,H,K,KH. When G=0, the same reduction holds for I,A,B,BA, so replacing
the letters defines the asserted algebra action; finite sums and dagger
are preserved as well. Conversely such an action preserves HK+KH=0.
This derives a native action without importing a representation theorem.
Insert (27.3) into R26's commutator-current identity to obtain (27.4).
Positivity makes its vanishing for every preparation equivalent to G=0.
Even checking the zero defect on a complete native role basis suffices,
since \(\|G e_j\|^2=0\) forces each column to vanish.

### R27.2 — Every finite source-law record has one forced two-role factor

Every nonzero source-law record admits a pairing-preserving factorization

\[
\mathcal T:\mathcal V\otimes\mathcal W\longrightarrow\mathcal M,
\qquad
A=\mathcal T(H\otimes I)\mathcal T^{-1},\quad
B=\mathcal T(K\otimes I)\mathcal T^{-1}.
\tag{27.5}
\]

Here \(\mathcal V\) is the original two-role carrier and
\(\mathcal W=(I+A)\mathcal M/2\) has the inherited matching pairing.
In particular the record's rank is even, the smallest nonzero rank is two,
and its native role-algebra action is faithful. No continuous choice of a
different mixed source law remains in this class.

Every arrow intertwining two such source-law records has, in these
factorizations, the form \(I_{\mathcal V}\otimes T\). It preserves the
pairing exactly when T does. Thus the residual freedom consists of a change
of record coordinates and a multiplicity target on which every source word
acts as identity.

**Proof.** The already available native cuts
\(E_\pm=(I\pm A)/2\) are self-adjoint complementary projections.
Their ranges are orthogonal: for u in the plus range and v in the minus
range, \(\langle u,v\rangle=\langle Au,Av\rangle=-\langle u,v\rangle\).
Anticommutation gives \(BE_+=E_-B\), so B bijectively and isometrically
exchanges the two ranges. Define, without choosing an orthonormal basis,

\[
\mathcal T(e_0\otimes u+e_1\otimes v)=u+Bv,
\qquad
\mathcal T^{-1}m=e_0\otimes E_+m+e_1\otimes BE_-m.
\tag{27.6}
\]

Both inverse identities follow directly. Orthogonality and preservation by
B prove its pairing identity. Applying A and B to the displayed formula
gives (27.5). The plus range is nonzero: otherwise A=-I and AB=-BA would
force B=0 on a nonzero target, contrary to its involution.

A basis of \(\mathcal W\), followed by its B image, gives a basis of
\(\mathcal M\); hence its rank is twice that of \(\mathcal W\). The
original native record realizes rank two. If a native role operator X acts
as zero, then \(X\otimes I_{\mathcal W}=0\); applying it to each source
role tensored with one nonzero w shows every coefficient of X is zero.
This proves faithfulness directly.

Finally write an intertwiner in four source-role blocks. Commutation with
H removes the two off-diagonal blocks; commutation with K identifies the
two diagonal blocks. The remaining common block is T. Its pairing identity
is precisely the pairing identity of the full intertwiner. No classical
classification or Schur-type lemma is a premise.

### R27.3 — The record normalization defect enters propagation by an exact identity

For all the involutive record candidates of (27.0), including those that
fail source normalization, put \(\bar G=I_s\otimes G\). Then

\[
V^2+V^{-2}=\tfrac12(S^2+S^{-2})I+\bar G,
\qquad
V^4-\bigl[\tfrac12(S^2+S^{-2})I+\bar G\bigr]V^2+I=0,
\tag{27.7}
\]

\[
\boxed{\|(V^2-I)\Psi\|^2
 =\|\Psi\|^2+\tfrac12\|(S^2-I)\Psi\|^2
  -\langle\Psi,\bar G\Psi\rangle.}
\tag{27.8}
\]

For scalar defect \(G=\gamma I\), the sharp normalized bounds are
\(1-\gamma\) and \(3-\gamma\). Source-law continuation forces
\(\gamma=0\) and recovers R26's bounds one and three without separately
postulating its anticommuting memory.

For a fixed record preparation define
\(g_\mu=\langle\mu,G\mu\rangle/\|\mu\|^2\). Constant even-address
packets with this record have limiting normalized defect \(1-g_\mu\),
and the native pairing inequality gives the calibration relation

\[
g_\mu^2\le j_2(0)+\tfrac12.
\tag{27.9}
\]

**Proof.** Direction-controlled shifts and the source coin preserve the
pairing, so \(V^\dagger=V^{-1}\). With
\(a=S^2,b=S^{-2},Z=AB\), the unnormalized two-step block is

\[
2V^2=
\begin{pmatrix}
a+Z&a-Z\\
Z^\dagger-b&Z^\dagger+b
\end{pmatrix}.
\]

Add its dagger. Both off-diagonal blocks cancel and both diagonal blocks
become \(a+b+2G\), proving the first identity. G commutes with the
record arrows and shifts; multiplication by \(V^2\) proves the second.
Expanding the defect norm and the shift-difference norm proves (27.8).

For scalar G, the shift-difference ratio is between zero and four. The
constant and alternating even-address packets of R26 give respectively
\(1-\gamma+1/M\) and \(3-\gamma-1/M\), proving sharpness. With a
fixed \(\mu\), the same first packet calculation replaces gamma by
\(g_\mu\). Finally
\(\langle\mu,G\mu\rangle^2\le\|\mu\|^2\|G\mu\|^2\), derived
from positivity of \(\|G\mu-t\mu\|^2\) at its minimizing native
coefficient, together with (27.4), proves (27.9). These are native
propagation and readout identities, without an assigned physical mass.

### R27.4 — R26's response is universal across finite source-law records

For every nonzero finite source-law record, the conjugacy (27.5) intertwines
its joint transport with
\(V_{H,K}\otimes I_{\mathcal W}\), after reordering the existing tuple
factors. The same conjugacy intertwines the retained-arrival flags and the
separate source-controlled return gate.

In particular the following R26 results apply unchanged: the two-event
propagation identity, the signed first-return grammar, its p and proved
tail, the parity-dependent feedback recurrence, and

\[
D_* =\frac{1-p}{1+p^2/4}
       (I-H_s/2-R_s/4-pK_s/4),\qquad
\frac{\kappa_*(e_1)}{\kappa_*(e_0)}=3.
\tag{27.10}
\]

Record rank, multiplicity, record coordinates and record preparation do
not supply an extra coupling for these targets. For an arbitrary joint
source/record preparation at the origin, local source readouts and the
specified feedback depend only on its source pairing matrix. Correlation
with the record does not alter this statement.

**Proof.** Substitute (27.5) into each term of (27.1). The conjugating
arrow acts only on record coordinates, so it commutes with the shift,
source cuts, arrival flags and the separate source gate. Finite words and
finite matching sums therefore intertwine term by term. R26's exact
return formulas then carry over, as do their norm-Cauchy completions and
the same rational error bounds; no new numerical fit or limit is needed.

There is also a direct readout proof. Anticommutation normal-orders every
direction word at fixed n,x into a sign times \(A^aB^b\), with
\(a=(n+x)/2,b=(n-x)/2\). Thus the endpoint map from the origin is
\(\phi_n(x)\otimes A^aB^b\), where \(\phi_n\) is R26's same signed
source stencil. For a source readout X its pulled-back matching operator is
\(\phi_n(x)^\dagger X\phi_n(x)\otimes I\). This proves the claim
for arbitrary joint input by direct contraction of record indices. R26's
feedback operator likewise equals \(D_N\otimes I\) at every finite N.
The scalar ratio in (27.10) refers to the two stated source preparations;
it does not make all preparations equivalent or identify particle charges.
Extra multiplicity marks and full source histories remain available to
other readout targets.

### R27.5 — Copying both source letters requires a minimum native compensator

Let A,B be a source-law record and consider a literal operator copy
\(\widehat H=H_s\otimes A\), \(\widehat K=K_s\otimes B\).
These two letters are pairing-preserving involutions but commute, so they
do not carry the original mixed source law on the joint target.

For a further independent finite record \(\mathcal F\), consider the
factorized repair

\[
H'=H_s\otimes A\otimes U,\qquad
K'=K_s\otimes B\otimes W,
\tag{27.11}
\]

with pairing-preserving arrows U,W. This lift carries the full native
H/K law exactly when
\(U^2=W^2=I\) and \(UW=-WU\). Its smallest nonzero compensator has
two native roles, realized by U=H,W=K. Every finite repair in this class
has the normal form of R27.2. With a minimal first record and compensator,
the joint target has eight role marks and
\((H'+K')/\sqrt2\) preserves pairing.

A one-dimensional native-iota phase can repair the **balanced norm** of
the uncorrected double copy, as in R26.9, but cannot repair these two
involutions and their full mixed word law simultaneously.

**Proof.** Two sign reversals give
\(\widehat H\widehat K=\widehat K\widehat H\). Expanding the squares
of (27.11) forces U and W to be involutions. Expanding its anticommutator
gives the invertible common factor
\(H_sK_s\otimes AB\) tensored with \(UW+WU\). Hence its vanishing
is exactly anticommutation on the compensator. The converse follows by
the same calculations. R27.2 now gives the even-rank minimum and
classification. Source normalization follows from R27.1, and the explicit
third native copy realizes all the required identities.

On a one-dimensional target over the commutative native cut-scalar sector,
nonzero scalar arrows u,w commute. Their anticommutation would imply
\(2uw=0\), impossible for involutions. In particular multiplying the
second doubled letter by iota makes its square -I instead of I. This
distinguishes norm repair from preservation of the source word law. The
minimum is proved for the factorized operator-copy target (27.11); it is
not an assertion about arbitrary encodings, unknown-state copying or a
new physical particle.

### R27.6 — Retained copy parity selects two native response classes

For r identical native record copies put
\(A_r=H^{\otimes r}\), \(B_r=K^{\otimes r}\), with
\(A_0=B_0=1\). Then

\[
A_rB_r=(-1)^rB_rA_r.
\tag{27.12}
\]

Odd r carries the source law, and every such record has the R26 response
(27.10). Even r is commuting memory: from one origin its unresolved local
source readouts and retained-arrival feedback are those of the unrecorded
source. Its first-return energy is R24's eta and its completed feedback is
R25's

\[
D_*^{\rm even}=
\frac{1-\eta}{1+\eta^4}
(I+\eta H_s-\eta^2R_s-\eta^3K_s).
\tag{27.13}
\]

Thus the retained-copy parity, not the number of extra replicas within
one parity, distinguishes these response targets. No source word or
record mark is deleted in making this readout statement.

The same sign count says a literal lift onto a **total** of k source
copies carries the original law exactly when k is odd. For odd k the
full action is one native two-role action with \(2^{k-1}\) multiplicity
marks. For positive even k, the balanced literal lift has a kernel of
rank \(2^{k-1}\) and Gram value two on its complementary sector.

**Proof.** Every tensor factor contributes the already derived native
minus sign on interchanging H and K. This proves (27.12), with no new
tensor axiom beyond R20's finite tuple pairing. The odd case follows from
R27.1 and R27.4. In the even case all paths at fixed n,x carry the common
record factor \(A_r^{(n+x)/2}B_r^{(n-x)/2}\). At a return of length 2m
this is \((A_rB_r)^m\), independent of the path and preserving pairing.
It cancels from first-return norms and from the separate controller cross
pairing. R24 and R25 therefore apply directly, including their completion
bounds. The empty record r=0 gives the same equations.

For clarity the odd-copy normal form has a completely explicit native
tag encoder. For k odd and binary role tuple \(z_1,\ldots,z_k\), put
\[
p=z_1\mathbin\oplus\cdots\mathbin\oplus z_k,\qquad
d_i=z_i\mathbin\oplus z_k\quad(1\le i<k).
\tag{27.14}
\]
Here xor is just parity of finite mark counts. Recover
\(z_k=p\oplus d_1\oplus\cdots\oplus d_{k-1}\), then
\(z_i=d_i\oplus z_k\). This is a bijective relabelling, hence preserves
matching pairing. The copied H gives the sign \((-1)^p\), and the
copied K flips p while fixing all d_i. Thus the source acts only on the
one logical role, leaving the \(2^{k-1}\) relative marks fixed.

For k positive and even, \(A_k,B_k\) commute and
\(Q=A_kB_k\) is a self-adjoint involution. The balanced lift has Gram
I+Q, so its minus cut is the kernel and its plus cut has Gram value two.
H acting on just the first tuple coordinate anticommutes with Q and
isometrically exchanges the two cuts. Their ranks are therefore equal,
each \(2^{k-1}\). This extends R26's double-copy obstruction. The
total-copy count k and the ancillary record count r differ by one in a
literal source-plus-record lift; they must not be conflated.

### R27.7 — Native counterfamilies locate the remaining physical-selection freedom

The replica classification does not follow from causality and preservation
of the **joint transport** alone. A family constructed solely from native
count coefficients is, for signed counts m,n not both zero,

\[
A=H,\quad B=\gamma H+\sigma K,\qquad
\gamma=\frac{m^2-n^2}{m^2+n^2},\quad
\sigma=\frac{2mn}{m^2+n^2}.
\tag{27.15}
\]

Every member gives a pairing-preserving involutive record pair and a
pairing-preserving local transport (27.1). Its calculated defect is
\(G=\gamma I\). Only the equal-magnitude count case m squared equals
n squared carries the original balanced source law. For positive counts
this is m=n, giving B=K. Thus source-law normalization removes a genuine
native freedom; it is not an automatic consequence of having an isometric
transport.

For the native pair (m,n)=(2,1),
\[
G=\tfrac35I,\quad j_2(0)=-\tfrac7{50},\quad
\inf\frac{\|(V^2-I)\Psi\|^2}{\|\Psi\|^2}=\tfrac25.
\tag{27.16}
\]

Moreover B and -B give exactly the same fixed-event, local source
readouts from a single origin, and the same retained-arrival feedback.
Their two-event defect bounds are generally different: in this example
the lower bound becomes 8/5 after the sign change. A propagation defect
relative to the identity therefore also needs its continuation/phase
convention to be specified before any physical mass interpretation.

**Proof.** The native count identity
\((m^2-n^2)^2+(2mn)^2=(m^2+n^2)^2\) gives
\(\gamma^2+\sigma^2=1\). H and K are self-adjoint involutions and
anticommute, so B is a pairing-preserving involution and
\((HB+BH)/2=\gamma I\). The transport is a source coin followed by
pairing-preserving controlled shifts, for every one of these choices.
Its record census has Gram \((1+\gamma)I\), proving the selection
claim. Equation (27.4) gives \(j_2=\gamma^2-1/2\), and R27.3 gives
the sharp lower bound \(1-\gamma\). Substituting the stated counts
proves (27.16).

On replacing B by -B, every path with n events ending at x acquires the
same sign \((-1)^{(n-x)/2}\), because this counts its negative steps.
The sign cancels from each local matching readout and from both controller
branches at a fixed arrival ledger and endpoint. It does not change the
underlying tagged histories. But G changes sign, so (27.8)'s comparison
with an unchanged identity continuation changes as stated. A relabelling
that also changes the continuation target is not a proof that these
identity-relative gaps are equal.

This supplies exact witnesses of the remaining freedom. Neither source
normalization nor these counterexamples assert that a physical interaction
must use direction control, literal-letter copying, a particular copy
count, or a particular arrival-memory target. Within the source-law
direction-record target, R27.2–R27.4 do remove the continuous record law,
rank and preparation freedoms from R26's specified responses.

## Certification and source boundaries

The [ledger](../04-operator-evolution/R27_DERIVATION_LEDGER.json) binds
these seven written proofs and their premise boundaries. The
[adapter](../04-operator-evolution/native_replica_selection.cjs) uses the
unchanged native arithmetic/replayer. Its two-involution presentation is
used only for the general defect identities; anticommutation is not silently
added to that presentation. Canonical H/K identities are replayed separately.
The [wrapper](../04-operator-evolution/verify_r27.py) replays R26 and the
unchanged earlier chain. All 215 prior non-navigation files are preserved.

Finite exact checks and dependency-metadata gates accompany the written
proofs; they are not a formal proof assistant or a physical experiment.
Native source-law selection is distinguished from physical interaction
selection in the [verification](../04-operator-evolution/R27_VERIFICATION.json).

The certificate records seven written results, ten exact check groups,
fourteen native identity replays, twelve rejected false alternatives and
nine invalid dependency-graph mutations. It includes reconstruction of
records of ranks two, four and six from their own cut ranges, exact
intertwiner-space calculations, independent literal H/K histories and
both copied-memory return classes. R26's and R25's numerical completion
certificates are read through exact hash pins rather than recomputed with
new fitted constants. The two-involution and canonical-source symbolic
presentations are separately audited and identified in each replay.

| Pinned source | Premise status and use |
| --- | --- |
| [R16 C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native derived:** count coefficients, role law, pairing, derived iota and completion; full word history remains distinct. Equality, finite counting and induction remain explicit metamathematical infrastructure. |
| [R18.1–R18.4](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | **Native derived on an explicit target:** address shift and balanced transport. No physical clock or unique metric selection is supplied. |
| [R20.1–R20.4](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native derived:** finite tuple pairing, controlled record and matching readout. Record allocation and readout targets remain constructions. |
| [R24](https://github.com/Parveen117/extra-ideas/blob/d660bbf026e1de818fcb3dd61d46346a4ad618ae/02-relational-response/NATIVE_RETURN_COUPLING_R24.md) and [R25](https://github.com/Parveen117/extra-ideas/blob/d4f1ae2481a72a04c465f1c17b37d154a6122eff/02-relational-response/NATIVE_RETURN_FEEDBACK_R25.md) | **Native derived:** first-return eta, coherent/stopped distinction and retained-arrival feedback for the unrecorded/flat class. No fitted coefficient is used. |
| [R26.1–R26.9](https://github.com/Parveen117/extra-ideas/blob/6773afc38caff84ff0b71ed35e3aa8c3a4adbd45/02-relational-response/NATIVE_MEMORY_RESIDUE_R26.md) | **Native derived:** commutator readout, gap, signed return p, rigorous tail, parity feedback and literal-copy obstruction. R27 proves their universality for finite source-law records; it does not reidentify the gap as physical mass. |
| [Canonical RKF algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native equality/replay source:** scoped use of the unchanged engine; no new whole-engine or proof-assistant claim. |
| [Publications AS-1–AS-5](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-alpha-selection/THEOREM.md) | **Comparison only, mixed premises:** supplied cell weights and an admitted positive quadratic physical adapter in AS-3. That imported adapter and electromagnetic normalization remain excluded from R27 proof paths. |

The next physical selection problem is now more precise: which native
continuation target and retained-copy/control structure is realized. The
record law within the specified source-law target is classified; a physical
choice of target, c, alpha, particle mass and charge remain unproved. Earlier
RH and Yang–Mills claim boundaries are unchanged.
