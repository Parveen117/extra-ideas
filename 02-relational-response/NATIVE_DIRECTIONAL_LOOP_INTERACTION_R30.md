# R30: native directions, closed protocol loops and field interaction

Research owner: Monty Dabas. Development: 1 October 2026.

R29 derived a paired field, its source-compatibility residue, a positive
conserved form and boundary interfaces. R30 constructs two address labels
from the same source history while retaining one shared H/K role. Their
order residue has an exact native factorization and a closed protocol-loop
readout. Intersecting interfaces create a localized curvature response.
That response then realizes R29's compatibility residue and closes a local,
reversible coupled field system. Native changes of local role frame preserve
these results; an exact finite-count calculation classifies periodic blind
sectors and the resulting positive interaction quotient.

The directions are address records, not assumed physical spatial axes. The
loop is a commutator of full response transports, not an assumed elementary
spacetime plaquette. The interaction is an explicitly composed response
target, not a claimed unique physical force law. No Maxwell equation,
external gauge potential, ordinary complex phase, Fourier transform or
classical spectral theorem enters these proofs.

## Native source and targets

Use the source H,K, R=KH, P=(I+H)/2, Q=(I-H)/2 and matching pairing of
R16–R29. The two source-role combinations satisfy

\[
C_0=H+K,\quad J_0=H-K,\quad C_0^2=J_0^2=2I,
\quad C_0J_0=-J_0C_0.
\tag{30.0}
\]

Let \(T_x,T_y\) shift their two integer address records by plus one.
Define \(\delta_i=T_i-I\), and reuse the R29 operator along either label:

\[
N_i=\tfrac12(P+T_i^{-1}Q)C_0,\qquad V_i=I+N_i\delta_i.
\tag{30.1}
\]

The native directions are constructed in R30.1. The source-role maps are
shared between them. Finite signed response fields use the sum of the same
native matching pairing. Original H/K word tags remain provenance even
when their response values are aggregated. As in R18, this aggregation is
an explicit mathematical target, not a selected physical superposition law.

## Written results

### R30.1 — One retained source history constructs two directional transports

Divide the internally counted source events into pairs. Assign successive
pairs alternately to labels x and y. During a pair, add each event's native
role contrast, plus/minus one, to its assigned record. The increment over
the pair is even; divide it by two to obtain an integer coarse address.
The full word, its signs and its current native role remain retained.

Signed aggregation and the inherited two-event normalization give V_x for
an x pair and V_y for a y pair. A four-event block therefore has

\[
W=V_yV_x,
\]

while reversing the address-label order gives \(W_{\rm rev}=V_xV_y\).
Each preserves the finite native pairing. Address shifts commute, but the
two response transports need not: both act on the same source role.

If instead each label is assigned a separate source-role copy, the two
transports commute. Thus the retained shared role is a load-bearing part
of this target. The construction does not derive the physical number of
spatial dimensions or select this grouping of source events over all
possible readouts.

**Proof.** R18 constructs each event's sign from H acting on the occupied
role. Two such signs sum to -2,0 or 2, which gives the coarse integer record.
For one assigned label the two unnormalized steps are precisely R18's
operator squared. Dividing by two gives the R29 matrix V_i. Reading the
labels in chronological order yields W or W_rev. The pairing proof is
the product of the two R29 isometries. T_x and T_y commute because updating
two distinct record entries commutes; this does not interchange their
shared-role coefficients. On separate role copies all coefficients act on
different tuple entries as well, so those complete transports commute.
Finite literal H/K histories independently verify this constructed census.

### R30.2 — Directional order and the actual closed loop have one quadratic residue

Define the oriented order residue, its scalar difference factor and the
closed protocol loop by

\[
\Omega=[V_x,V_y],\qquad
g=\delta_x\delta_y(T_y^{-1}-T_x^{-1}),\qquad
\mathcal L=V_xV_yV_x^\dagger V_y^\dagger.
\]

Then, exactly,

\[
\boxed{\Omega=\tfrac14gJ_0,\qquad
       \Omega^\dagger\Omega=\tfrac18g^\dagger g\,I,}
\tag{30.2}
\]

\[
\mathcal L-I=\Omega V_x^\dagger V_y^\dagger,\qquad
\boxed{(\mathcal L-I)^\dagger(\mathcal L-I)
       =\Omega^\dagger\Omega.}
\tag{30.3}
\]

Thus the two signed targets are different operators but have the same
quadratic reading on the same input. This equality is proved for this
source, not assumed for arbitrary curvature/holonomy readouts. Neither
target is energy lost by a non-isometric transport; each protocol itself
preserves the source pairing.

On the unwrapped address plane, Omega annihilates no nonzero finite-support
field. For a unit packet at (0,0), its response energy is exactly 3/4.
Noncompact stripe profiles and finite periodic quotients have different
blind sectors, treated below.

**Proof.** The scalar differences commute with all native coefficients.
Since V_i=I+N_i delta_i,
\(\Omega=\delta_x\delta_y[N_x,N_y]\). Multiplying the two R29 matrices
gives \([N_x,N_y]=(T_y^{-1}-T_x^{-1})J_0/4\). Equation (30.0) proves
the dagger square in (30.2). Multiplying out the loop gives its first
identity. Its squared defect is a conjugate of Omega-dagger Omega; that
operator is a scalar shift polynomial and commutes with both transports,
which proves (30.3).

The difference factor is the six-term native Laurent polynomial

\[
g=T_x-T_y-T_xT_y^{-1}+T_y^{-1}+T_x^{-1}T_y-T_x^{-1}.
\tag{30.4}
\]

Its largest lexicographic shift is (1,0), with coefficient one. In a
nonzero finite field choose its largest occupied address; the largest
address of the product has the unique nonzero coefficient contributed by
these two largest terms. It cannot cancel. Since J_0 is invertible, this
proves the finite-support injectivity of Omega without a spectral theorem.
On a single packet the six shifts are distinct; (30.2) gives 6/8=3/4.

### R30.3 — Curvature has a derived dual continuation and a loop-response decoder

The native role involution \(J=J_0/\sqrt2\) reverses each directional
response:

\[
JV_iJ=V_i^\dagger,\qquad
\boxed{\Omega V_i=V_i^\dagger\Omega.}
\tag{30.5}
\]

For the alternating block W=V_y V_x, put
\(\widehat W=V_y^\dagger V_x^\dagger=(V_xV_y)^\dagger\). Then

\[
\Omega W^m=\widehat W^m\Omega,\qquad
\boxed{\mathcal L-I=V_xV_y\Omega.}
\tag{30.6}
\]

Both the oriented curvature response and the actual loop-defect response
therefore have conserved quadratic size under repeated W. The latter is
decoded into the former by the known inverse continuation
\(V_y^\dagger V_x^\dagger\). For two noncommuting directional blocks,
widehat W is the inverse of the reversed order, not generally W-dagger.

**Proof.** Insert (30.1) and multiply the source role coefficients to get
\(J_0V_iJ_0/2=V_i^\dagger\). The scalar g commutes with V_i, so (30.2)
gives (30.5). Iterating it proves the first identity in (30.6). It also
gives \(\Omega V_x^\dagger V_y^\dagger=V_xV_y\Omega\), proving the
decoder identity from (30.3). The dual transports and the decoder are
native isometries, which proves norm conservation. These are algebraic
inverse continuations, not assertions of reversed physical time.

### R30.4 — Intersecting cut interfaces create a localized oriented loop response

Let s(k)=1 for k<0 and zero otherwise. Its native difference is the unit
packet \(\delta s=\delta_0\). For a native role value v, form the corner
profile \(\phi(x,y)=s(x)s(y)v\). Its curvature response is finite:

\[
\boxed{\Omega\phi=\tfrac14
 (\delta_{(0,-1)}-\delta_{(-1,0)})J_0v,\qquad
 \|\Omega\phi\|^2=\tfrac14\|v\|^2.}
\tag{30.7}
\]

Reversing one interface orientation reverses the signed residue while
preserving its squared size. Its total signed response is zero. The
actual closed-loop response is the finite native continuation
\(V_xV_y\Omega\phi\), with the same norm.

Every profile depending on x alone, y alone, or x+y alone is invisible to
this loop reading, including a straight interface carrying the R29 end
residue. Thus a straight boundary residue and two-direction loop response
are different targets. A nonzero loop reading here detects the meeting
and relative direction of interfaces, not a declared electric charge.

**Proof.** Apply delta_x delta_y to the corner to obtain delta_(0,0) v;
the remaining two reverse shifts in g give (30.7). The two packets are
disjoint, and J_0-dagger J_0=2I proves the norm. Replacing s(x) by 1-s(x)
adds a profile independent of x and changes the sign of the corner term.
For a profile depending on only one address, one of the first two factors
of g vanishes. On a profile of x+y, T_x and T_y act identically, so its
third factor vanishes. The loop assertion follows from (30.6). All response
sums are finite; no norm is assigned to the infinite quadrant background.

### R30.5 — The loop reading closes a reversible native field interaction

Construct a joint response target with finite fields (b,phi). Use the
oriented loop reading as R29's compatibility residue:

\[
\kappa=\Omega\phi,\qquad e=N_xb+\kappa.
\tag{30.8}
\]

The following update is closed and reversible:

\[
\boxed{\phi'=V_x\phi,\qquad
       b'=V_xb+\delta_x\Omega\phi.}
\tag{30.9}
\]

It realizes the exact paired-field equations

\[
b'-b=\delta_xe,\quad e'-e=-\tfrac12\delta_x^\dagger b',\quad
\kappa'=V_x^\dagger\kappa.
\tag{30.10}
\]

The pulled-back R29 form

\[
\boxed{\mathcal E_{\rm loop}(b,\phi)
 =\mathcal E_{R29}(N_xb+\Omega\phi,b)}
\tag{30.11}
\]

is conserved. It is positive definite on finite-support joint fields on
the unwrapped plane. This form includes the paired cross term; the simple
sum \(\|b\|^2+\|\phi\|^2\) is generally not conserved.

For every finite row, \(\sum_xb(x,y)\) is unchanged. The local b-current
is R29.5 with y held fixed, plus the exact source contribution

\[
2\langle(V_xb)_{x,y},(\delta_x\kappa)_{x,y}\rangle
       +\|(\delta_x\kappa)_{x,y}\|^2.
\tag{30.12}
\]

For example b=0 and phi=delta_(0,0)e0 produce a nonzero first b response
with squared size 2; the complete form (30.11) remains 3/4. This is a
transfer within a coupled response ledger, not energy created from the
zero joint state or a physical energy claim.

**Proof.** Equation (30.5) gives
\(\Omega\phi'=V_x^\dagger\Omega\phi\). Substitute
c=Omega phi in R29.3; its triangular update is exactly (30.9), and R29.2
then gives (30.10). Equivalently, the source-built observable map
(b,phi) -> (N_x b+Omega phi,b) intertwines this joint update with the
full R29 paired map. R29.4 therefore proves (30.11) and its positivity
except on b=0, Omega phi=0. R30.2 makes that kernel trivial for compact
fields. The identities hold at every passive y index, so finite summation
extends the R29 proof without a new geometric assumption.

The inverse joint matrix is

\[
\begin{pmatrix}V_x^\dagger&-\delta_x\Omega\\0&V_x^\dagger\end{pmatrix}.
\]

Multiplication cancels the off-diagonal terms by (30.5), proving local
reversibility. Summing delta_x e telescopes on each row, proving the
boundary statement. Expanding \(\|V_xb+\delta_x\kappa\|^2\) gives
(30.12). The stated packet values follow by finite source multiplication.

The coupling target in (30.8) fixes its signed readout calibration by the
actual order difference. Its closure and coefficients are derived. A
rescaled readout would be a different calibration; these results do not
select a unique physical interaction strength or identify electric charge.

### R30.6 — Native local role-frame changes preserve curvature and cannot create flat-link flux

Construct a pairing-preserving local frame G_z from native H/K words and
native count turns, for example

\[
G_z(a,b)=\frac{(a^2-b^2)I+2abR}{a^2+b^2},\qquad(a,b)\ne(0,0).
\tag{30.13}
\]

These are source operators, not primitive ordinary-complex phases. Let
G act pointwise on the address field. Change all role and transport
readouts together by \(A^G=GAG^\dagger\), and fields by phi^G=G phi.
Then

\[
\boxed{\Omega^G=G\Omega G^\dagger,\quad
       \mathcal L^G=G\mathcal L G^\dagger,\quad
       \|\Omega^G\phi^G\|=\|\Omega\phi\|.}
\tag{30.14}
\]

The induced pure shift link is
\(\ell_i(z)=G_zG_{z-e_i}^\dagger\). Products along two orders telescope
to the same endpoint frame. Hence the transformed pure shifts commute
and have identity loop. A nonzero native V-loop cannot be removed by
changing local role labels.

The position-dependent image G_z q of a uniform background obeys
\(\delta_i^G(Gq)=0\), although its untransformed ordinary address
difference can be nonzero. Keeping numerical H/K or plain shifts fixed
while changing the state frame can therefore manufacture a false field.
The paired interaction and its complete energy form are covariant under
the same simultaneous transformation.

**Proof.** R-dagger=-R and R-squared=-I give G_z-dagger G_z=I because
\((a^2-b^2)^2+(2ab)^2=(a^2+b^2)^2\). Products with H/K words preserve
this identity. In every composed operator word adjacent G-dagger G
cancel, proving (30.14), the transformed field equations and energy
invariance. The same cancellation gives
\(\ell_x(z)\ell_y(z-e_x)=G_zG_{z-e_x-e_y}^\dagger
=\ell_y(z)\ell_x(z-e_y)\). Finally
\(\delta_i^G Gq=G\delta_iq=0\). This is local frame covariance
of explicitly constructed native readouts. It does not select a physical
gauge group or turn arbitrary frame parameters into a physical gauge field.

### R30.7 — Three native averaging cuts exactly classify periodic curvature blindness

Choose the finite count quotient \((\mathbb Z/L\mathbb Z)^2\), L>=1.
This is an explicit boundary target. Put \(T_d=T_xT_y^{-1}\) and

\[
A_i=\frac1L\sum_{k=0}^{L-1}T_i^k\quad(i=x,y,d),\qquad
A_\circ=A_xA_y,\qquad
\boxed{\Pi=A_x+A_y+A_d-2A_\circ.}
\tag{30.15}
\]

Pi is the exact native pairing-orthogonal projector onto ker(Omega).
Its image consists of sums of x-only, y-only and (x+y)-only profiles.
Counting their constant overlap gives

\[
\boxed{\dim\ker\Omega=2(3L-2),\qquad
       \operatorname{rank}\Omega=2(L-1)(L-2).}
\tag{30.16}
\]

The factor two counts native source roles. In particular the L=2 quotient
is completely blind to this loop, despite the nonzero compact response
on the unwrapped plane. It is a boundary/readout aliasing effect.

There is an exact reconstruction of every curvature-visible component.
Define the finite weighted-count sums

\[
B_i=\frac1L\sum_{k=0}^{L-1}kT_i^k.
\]

Then

\[
\boxed{(I-\Pi)\phi
       =2T_xB_xB_yB_dJ_0\Omega\phi.}
\tag{30.17}
\]

No spectral decomposition or Fourier modes are needed.

**Proof.** Reindexing the finite cyclic sum gives A_i-squared=A_i and
A_i-dagger=A_i. Any pair of the three directional cycles visits all L^2
addresses once, so A_i A_j=A_circle for distinct i,j. Consequently the
three A_i-A_circle are mutually orthogonal cuts, and
Pi=A_circle+sum_i(A_i-A_circle) is a cut. Each A_i has one free native
role value on each of L transverse classes; their common sector is the
two-role constant. Thus its rank is 2+3(2L-2)=2(3L-2).

Finite coefficient subtraction gives
\((T_i-I)B_i=I-A_i\). Since
\(g=T_x^{-1}(T_x-I)(T_y-I)(T_d-I)\),

\[
T_xB_xB_yB_dg=(I-A_x)(I-A_y)(I-A_d)=I-\Pi.
\]

Each A_i is annihilated by g, so g Pi=0; the last identity proves the
converse kernel inclusion. J_0 is invertible and commutes with the scalar
cuts, hence ker(Omega)=im(Pi). Insert J_0 Omega=g/2 to get (30.17).
The remaining rank follows from the finite source coordinate count.
Invariance under T_d means a profile depends on x+y, giving the stated
three classes. This proof is entirely a native finite-count construction.

### R30.8 — The coupled response descends to an exact positive observable quotient

On the L-by-L count quotient, the joint observable

\[
\mathcal O(b,\phi)=(N_xb+\Omega\phi,b)
\]

has exactly the kernel

\[
\ker\mathcal O=\{(0,\phi):\phi\in\operatorname{im}\Pi\}.
\tag{30.18}
\]

This kernel is invariant under the joint interaction. Its closed observable
quotient therefore has

\[
\boxed{4L^2-2(3L-2)=4L^2-6L+4}
\tag{30.19}
\]

native real signed roles. The pulled-back form (30.11) has exactly this
rank and becomes positive definite on that quotient. Curvature-visible
phi is recovered by (30.17); no invisible source value has to be erased
or declared nonexistent. This blind sector is generally larger than
R29's uniform-background kernel because the observation target is different.

**Proof.** The second component of O first forces b=0; its first then
forces Omega phi=0. R30.7 proves (30.18). The scalar averaging cuts
commute with V_x,V_y, and Omega Pi=0. Thus blind phi continues within
the same sector and contributes no drive to b. The interaction consequently
descends through this kernel. R29's positive paired form has no additional
null directions in its own field variables, so its pullback has exactly
the kernel (30.18). Subtracting its dimension from the 4L^2 joint source
coordinates gives (30.19), positivity on the quotient and completeness
of the curvature reconstruction. This is a response-specific quotient,
not a classification of all physical observers or of the full UGD source.

## Certification and source lineage

Eight written proofs are bound to complete native Laurent identities,
literal H/K histories, finite coupled evolution, local current balances,
count-quotient reconstruction/rank checks, native word replay and negative
controls. A declared proof graph rejects imported/admitted/open premises
and promotion of target definitions to physical selection. This metadata
check does not inspect arbitrary mathematical prose semantically. Written
proof, exact computation, source integrity, formal verification and physical
validation remain distinct evidence classes. The frozen R29 and earlier
source chain are replayed without changing their certificates.

| Pinned source | Premise status and use |
| --- | --- |
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native-derived:** signed/refined coefficients, source roles, H/K, pairing and positive square-root completion. Finite equality, counting and induction remain declared proof infrastructure. |
| [R18.1–R18.5](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | **Native-derived on an explicit target:** role-contrast address and signed source transport. R30 constructs a two-label history readout, not physical spatial dimension. |
| [R28.1–R28.6](https://github.com/Parveen117/extra-ideas/blob/9ed6ae27d3ddb005ed7f4f3f27417036729698fd/02-relational-response/NATIVE_CUT_CURVATURE_PROPAGATION_R28.md) | **Native-derived:** signed curvature preparation, cut/loop readout distinction and propagation. No physical vacuum or unique dynamics is supplied. |
| [R29.1–R29.7](https://github.com/Parveen117/extra-ideas/blob/16602a29c469ac96001977432fcf5fc4b548970b/02-relational-response/NATIVE_PAIRED_FIELD_DYNAMICS_R29.md) | **Native-derived on stated readouts:** V,N, paired induction, compatibility evolution, positive form, current and boundary residue. R30 reuses those proofs with a passive second address and constructs c=Omega phi. |
| [R26.1](https://github.com/Parveen117/extra-ideas/blob/6773afc38caff84ff0b71ed35e3aa8c3a4adbd45/02-relational-response/NATIVE_MEMORY_RESIDUE_R26.md) | **Native prior attribution:** reciprocal-record frame cancellation. R30.6 extends the cancellation to two address labels and keeps local role frames distinct from physical gauge fields. |
| [Canonical RKF algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native arithmetic/equality replay:** unchanged active engine, scoped use only; no whole-engine or proof-assistant PASS is asserted. |
| [R28 comparison audit](https://github.com/Parveen117/extra-ideas/blob/9ed6ae27d3ddb005ed7f4f3f27417036729698fd/04-operator-evolution/R28_DERIVATION_LEDGER.json) | **Excluded comparison premises:** declared topology, declared-carrier information models and admitted classical smooth/tangent/metric adapters retain their original status. None is imported into these new proof paths. |

The result is a native loop-sensitive interaction model with covariant
readouts, conserved complete energy and boundary residue, and an exactly
classified finite observation kernel. Electromagnetic interpretation still
requires a physical field/source identification, justified spatial structure,
gauge-group selection and operational charge/coupling readout. No physical
vacuum, EM field, c, alpha, charge quantum or particle mass is certified by
renaming the constructed targets. Those physical questions remain open;
the displayed native identities and closed response constructions are the
completed mathematical development.
