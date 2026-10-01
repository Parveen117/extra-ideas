# R32: retained native loop, propagation gap and localized inverse response

Research owner: Monty Dabas. Development: 1 October 2026.

R31 showed that its triangular loop-driven target can be unmixed into flat
source continuations. R32 constructs a different retained-record interface:
the same native source-law record is carried reciprocally along both count
directions. Its record-link loop is -I. Individual reverse edges still
cancel, but one local frame cannot flatten both directions simultaneously.

The resulting source transport has an exact paired wave with a positive
identity term and mixed difference square. Its normalized paired-increment
defect lies sharply between one and four, independently of count volume.
This obstructs similarity to the flat continuation in a stated class.
Uniform signed states have a 48-source-event cycle. The full retained
increment is an invertible observer, with a finite decoder error bound and
a localized inverse whose native tail ratio is 3-2 sqrt(2).

The interface and its address/control convention are explicit constructions.
The record law and all displayed identities are derived within them. This
does not select a physical force, gauge group, mass, c or alpha. A staggered
source/record factorization remains available and is proved below; the
obstruction is to flattening the continuation, not to every factorization.
No ordinary complex field, classical metric, gauge potential, Fourier
transform, wave equation or physical clock enters as a premise.

## Native source and retained interface

Use the R16 source roles and matching pairing; R20 constructs the joint
ledger on source/record mark tuples. The tensor notation below abbreviates
those finite coefficient arrays. It is not an independent rule for nature
or a primitive Hilbert-space tensor product.

The moving source has H_s,K_s,R_s=K_sH_s, cuts P_s,Q_s and
C_s=H_s+K_s. One memory has H_m,K_m, the same native source-law action.
R27 derives its anticommutation from the original balanced census and proves
the smallest nonzero finite real signed record has two roles. Thus the joint
target here has four signed roles per retained address.

Keep R30's two coarse count labels and scalar shifts T_x,T_y. Assign H_m
to an x edge and K_m to a y edge, with the inverse on the reverse edge.
The resulting source/record shift operators are

\[
X=T_x\otimes H_m,\qquad Y=T_y\otimes K_m,
\quad XY=-YX,\quad X^2=T_x^2I,\quad Y^2=T_y^2I.
\tag{32.0}
\]

Source roles are implicit identity factors in X,Y. Conversely, source
cuts and C_s below act as identity on the record. A dagger is the inherited
matching adjoint. Both X and Y preserve the pairing. Choosing this
two-direction reciprocal interface is a target construction, not a theorem
that every possible native record interface or physical process uses it.

## Written results

### R32.1 — Reciprocal edges can retain a loop that no common native frame flattens

For the source-law record in (32.0),

\[
\boxed{XYX^\dagger Y^\dagger=-I,\qquad
[X,Y]^\dagger[X,Y]=4I.}
\tag{32.1}
\]

These are record-link readings. They are distinct from a loop of the full
source transports constructed in R32.2. Reversing a single edge has zero
residue; comparing two different orders around a cell does not.

More generally, a reciprocal record on the unwrapped count plane admits
an endpoint frame that makes every edge identity exactly when every closed
path arrow is identity. Identity elementary cell loops suffice there.
On a finite periodic count target, winding loops must also be identity.
The minus-identity loop in (32.1) prevents a common frame from flattening
both X and Y. Even an arbitrary simultaneous similarity cannot turn these
two anticommuting shifts into two commuting bare shifts.

**Proof.** H_m K_m=-K_m H_m and both squares equal identity by R27's
source-law result. Reordering the four letters gives the first equality.
Also [X,Y]=2XY; XY preserves pairing, proving the second equality.
This fixes the loop residue without a supplied phase or coupling angle.

Under pointwise native frames G_z, every path arrow O becomes
G_end O G_start^(-1). A closed arrow is only conjugated. In particular
-I remains -I. Conversely, if all closed paths are identity, define G_z
as the arrow of any finite path from the origin to z. Two choices have
the same arrow because their joined closed path is identity. Each edge
then equals G_end G_start^(-1), proving the flattening. On the unwrapped
integer count plane, a finite closed word can be reduced by inverse-step
cancellation and interchanges of adjacent perpendicular steps. Identity
cell loops justify each interchange. Periodic words additionally have
winding words, so cell identities alone do not prove the periodic claim.
Finally simultaneous similarity preserves a commutator, excluding a
similarity to commuting shifts directly. None of this assumes a classical
topological or gauge-theoretic theorem.

### R32.2 — Original source events construct the reciprocal retained transport

For A=X or Y construct

\[
\mathcal N_A=\tfrac12(P_s+A^{-1}Q_s)C_s,\qquad
\boxed{\mathcal V_A=I+\mathcal N_A(A-I).}
\tag{32.2}
\]

One directional block consists of two original source events. Assign the
record arrow H_m or K_m to one of the two fine edges making each coarse
edge, and identity to the other; reverse crossings use inverses. One
explicit realization marks exactly the fine edges whose lower endpoint
count is odd. Two fine events ending one coarse address away record the
axis arrow, while a two-event return records identity.

Thus (32.2) is the signed aggregation of the original H/K event histories
with reciprocal retained edge records. It preserves joint pairing and has
a local inverse. Alternating the x and y blocks gives

\[
\mathcal W=\mathcal V_Y\mathcal V_X,
\qquad\mathcal B=\mathcal W^2.
\tag{32.3}
\]

W-script counts four original events; B-script counts eight. Finite fields
stay within one coarse address step per axis per W-script block. All full
source words remain provenance even when their coefficient values agree.

**Proof.** For a single direction, reciprocal fine-edge arrows telescope
to their endpoint record, as in R26.1. A two-step crossing has H_m or K_m
in either direction, since these native arrows are involutions; a return
has identity. Insert those three endpoint arrows into the R29 two-event
stencil to get exactly (32.2). Either allocation of the one nonidentity
fine edge yields the same coarse endpoint operator. Its fine-edge placement
is a realization of the declared coarse interface, not an extra physical law.

The source coefficients commute with A, which is pairing-preserving.
Multiplication of (32.2) with its dagger therefore repeats the R29 native
single-shift identity and gives V_A-dagger V_A=V_A V_A-dagger=I.
Its inverse uses the same finite inverse-shift stencil. Products prove
pairing preservation and locality of W-script and B-script. Literal fine
source histories independently verify this construction. The edge record
is controlled by the address crossing; it is not literal copying of the
original H/K event tag, whose different law was separated in R27.

### R32.3 — Native area parity changes source response while retaining record information

On the unwrapped count plane put

\[
G_{x,y}=H_m^xK_m^y,\qquad
(\mathcal D\Phi)(x,y)=(I_s\otimes G_{x,y})\Phi(x,y).
\]

Then

\[
\boxed{\mathcal D^\dagger X\mathcal D=T_xI_m,
\qquad\mathcal D^\dagger Y\mathcal D=(-1)^xT_yI_m.}
\tag{32.4}
\]

The record can consequently be factored into multiplicity with a staggered
source transport. The staggered shifts still anticommute and retain the
minus-identity cell loop. Factorization has not flattened the source.
This frame is periodic on even count periods; for odd periods its boundary
crossings retain the additional winding arrows. The unwrapped formula
must not be silently imposed as a periodic odd-count frame.

For a one-origin product preparation v_s tensor mu, every endpoint has

\[
\Psi_n(x,y)=\phi_n(x,y)\otimes G_{x,y}\mu.
\tag{32.5}
\]

Local source matching responses are therefore independent of a normalized
mu in this input class. Nevertheless they differ from the flat source.
For v_s=e0, after two W-script blocks the origin squared response is

\[
\boxed{\rho_{\rm retained}(0,0)=\tfrac1{64},\qquad
       \rho_{\rm flat}(0,0)=\tfrac9{64}.}
\tag{32.6}
\]

**Proof.** H_m G_(x-1,y)=G_(x,y), whereas
K_m G_(x,y-1)=(-1)^x G_(x,y), proving (32.4), including inverse shifts.
Apply the same conjugation to the source polynomial (32.2); all record
coefficients become identity and all y crossings acquire the displayed
sign. The initially factored record at the origin stays a common factor
in this transformed evolution, proving (32.5). Its pairing-preserving
endpoint G cancels from local source matching readings but is still retained.
No record-state erasure or stochastic noise law has been used.

An a-by-b native rectangle has record loop (-1)^(ab). More generally,
interchanging two perpendicular path steps changes the normal-order sign
once. Counting such elementary interchanges modulo two gives the relative
sign of paths with the same endpoints. This is a parity of native counted
cells, not an imported physical area or magnetic flux. Expanding two
four-event blocks in the source stencil gives (32.6). Independence of mu
is not claimed for arbitrary multi-origin preparations with unrelated
record factors, nor does it make unsigned intensities a complete observer.

### R32.4 — The retained source derives an exact gapped paired wave

Use scalar two-address differences, distinct from the record shifts X,Y:

\[
D_i=T_i^2-I,\quad\Lambda_i=D_i^\dagger D_i,
\quad Q_i=I+\Lambda_i/4,\qquad\mathscr G=Q_xQ_y.
\tag{32.7}
\]

The full retained transport has the exact native identity

\[
\boxed{\mathcal W^4-(2I-\mathscr G)\mathcal W^2+I=0,
\qquad\mathcal B+\mathcal B^\dagger=2I-\mathscr G.}
\tag{32.8}
\]

At the eight-event block resolution this gives the derived wave

\[
\Psi_{n+1}-2\Psi_n+\Psi_{n-1}=-\mathscr G\Psi_n,
\qquad\Psi_{n+1}=\mathcal B\Psi_n.
\tag{32.9}
\]

Its increment F=B-script-I has the positive factorization

\[
\boxed{F^\dagger F=\mathscr G
 =I+\tfrac14(\Lambda_x+\Lambda_y)+\tfrac1{16}\Lambda_x\Lambda_y.}
\tag{32.10}
\]

Equivalently, for every finite source/record field,

\[
\boxed{\|F\Psi\|^2=\|\Psi\|^2
 +\tfrac14\|D_x\Psi\|^2+\tfrac14\|D_y\Psi\|^2
 +\tfrac1{16}\|D_xD_y\Psi\|^2.}
\tag{32.11}
\]

The normalized squared defect has sharp bounds one and four. Its positive
identity term comes from the retained source law; no mass term was supplied.
This is an increment relative to identity in the fixed native protocol,
not dissipation, a physical Hamiltonian mass gap or a Yang--Mills result.

**Proof.** Put A_s=(I_s+R_s)/2,
Z_i=P_s T_i+Q_s T_i^(-1), and E_i=Z_i A_s-dagger. Multiplication of
(32.2) gives V_X=A_s tensor I+E_x tensor H_m and the analogous expression
with E_y,K_m. Hence

\[
\mathcal W=A_s^2\otimes I+A_sE_x\otimes H_m
 +E_yA_s\otimes K_m+E_yE_x\otimes R_m.
\]

Multiply these four terms twice. Use R_s-squared=R_m-squared=-I,
H_m K_m=-R_m, K_m H_m=R_m, and commuting scalar count shifts.
In the native sixteen-role-operator basis, the coefficients of
W-script-squared plus W-script-dagger-squared are:

| Basis term | Coefficient |
| --- | --- |
| I_s tensor I_m | 1-(Lambda_x+Lambda_y)/4-Lambda_x Lambda_y/16 |
| Every other term in {I,H_s,K_s,R_s} tensor {I,H_m,K_m,R_m} | 0 |

These finite multiplications prove the second identity in (32.8) and,
after multiplication by W-script-squared, the first. The certificate
checks every Laurent coefficient, not a list of substituted phases.
Expanding (B-I)-dagger(B-I) proves (32.10). The scalar differences commute,
so pairing with the field gives (32.11).

The lower bound follows by positivity. Each T_i-squared preserves pairing,
so ||D_i v||<=2||v|| and ||D_xD_y v||<=4||v||; (32.11) gives the upper
bound four. For sharpness on the unwrapped plane, put a fixed vector at
the M-by-M addresses (2a,2b), 1<=a,b<=M. The constant-sign packet gives
ratio (1+1/(2M))^2. The packet with signs (-1)^(a+b) gives
(2-1/(2M))^2. They approach one and four respectively. Periodic uniform
fields attain the lower value; period-four sign patterns with T_i^2=-I
attain the upper value. Completion means the native Cauchy completion of
finite matching fields already constructed in R26.4, on which the same
identities extend by their explicit norm bounds.

### R32.5 — Flat continuation is obstructed and uniform signed states have a fixed count cycle

The retained W-script cannot be transformed into the unrecorded R31 W
with inactive memory by an invertible similarity on a finite cyclic target.
On the completed unwrapped matching target, no bounded, boundedly invertible,
time-independent similarity does so. This is a stronger transport statement
than the local record-loop obstruction, but it still permits the staggered
factorization in R32.3.

On the uniform four-role sector of a finite cyclic target, W-script_0 obeys

\[
\boxed{\mathcal W_0^4-\mathcal W_0^2+I=0,\qquad
\mathcal W_0^6=-I,\qquad\mathcal W_0^{12}=I.}
\tag{32.12}
\]

Every nonzero fully signed uniform state has least period twelve four-event
blocks: exactly 48 original source events. An intensity-only observation
need not have this period; it already forgets the global sign at six blocks.
No physical time unit, frequency, rest energy or particle identity has been
assigned to this internal source cycle.

**Proof.** Equation (32.11) makes B-script-I injective on every finite
cyclic target, whereas the bare W-squared-I annihilates uniform states.
Similarity would preserve this kernel dimension, a contradiction.

On the unwrapped completion, normalized constant M-by-M bare source packets
have ||(W_bare^2-I)psi_M||^2<=16/M. To see this, use V_i-I=N_i delta_i,
N_i-dagger N_i=I/2, and the two address boundary differences, each of squared
size 2/M. The triangle inequality bounds W_bare-I by 2/sqrt(M), and its
second block by twice that value. A bounded invertible similarity would
carry these approximate stationary packets to approximate stationary
retained packets, contradicting the lower bound one in (32.11).

For a uniform state T_x=T_y=I, so G-script=I. Equation (32.8) gives the
first identity in (32.12). Multiplication by W_0-squared+I gives
W_0^6+I=0, then squaring gives the twelve-block return. A least positive
period must divide twelve by finite integer division. Periods 1,2,3 and 6
contradict W_0^6=-I. Period four would give W_0^2 v=2v from the first
polynomial, contradicting pairing preservation for a nonzero v. The least
period is therefore twelve. The gap is stated in this fixed protocol frame;
arbitrary event-dependent changes of frame are not ruled out or promoted
to a physical symmetry by the theorem.

### R32.6 — The full retained increment is a complete observer with a uniform decoder bound

Read the complete signed source/record field

\[
d=F\Psi=(\mathcal B-I)\Psi.
\]

Unlike the flat increment, this reading has no uniform kernel. The exact
inverse on a finite cyclic target, and on the native norm completion, is

\[
\boxed{\Psi=Q_y^{-1}Q_x^{-1}F^\dagger d.}
\tag{32.13}
\]

On an L-by-L count target its rank is 4L^2. All four signed source/record
roles are part of the reading; this is not a claim that unresolved source
intensities or one scalar sensor recover an arbitrary record state.

There is an explicit native finite decoder. Put

\[
A_i=\tfrac12(T_i^2+T_i^{-2}),\quad
P_{i,N}=\tfrac23\sum_{j=0}^{N}(A_i/3)^j,\quad
K_N=P_{y,N}P_{x,N}F^\dagger,\quad r_N=3^{-(N+1)}.
\tag{32.14}
\]

Then, independently of the number of retained addresses,

\[
\boxed{\|K_Nd-\Psi\|
 \le(2r_N+r_N^2)\|\Psi\|
 \le(2r_N+r_N^2)\|d\|.}
\tag{32.15}
\]

**Proof.** The native average A_i has norm at most one by the triangle
inequality. Since Q_i=(3/2)(I-A_i/3), multiplying the finite sum gives
P_(i,N) Q_i=I-(A_i/3)^(N+1). Its geometric norm tail is at most
3^(-(N+1)); hence its native Cauchy completion is a two-sided inverse
of Q_i. These scalar shift inverses commute with F. Equation (32.10)
proves (32.13); the corresponding identity FF-dagger=G-script follows
as well from B-script pairing preservation, so the inverse is two-sided.
On a finite cyclic target, native elimination gives the same inverse.

For the finite decoder, write E_i=(A_i/3)^(N+1). Exact multiplication gives

\[
K_NF=(I-E_y)(I-E_x).
\]

Both E norms are at most r_N. The remainder E_x+E_y-E_yE_x therefore
has norm at most 2r_N+r_N^2. Finally (32.11) gives ||Psi||<=||d||,
proving (32.15). This is a finite-error statement, not exactness at a
generic finite depth. The rank and completeness claims follow from the
two-sided inverse, without a classical spectral theorem.

### R32.7 — The inverse response has a source-derived localized count kernel

Let the positive native square root of two be the R16 completed cut scalar,
and define

\[
\boxed{\rho=3-2\sqrt2,\qquad 0<\rho<1.}
\tag{32.16}
\]

Then the inverse factors and the complete scalar inverse kernel are

\[
\boxed{Q_i^{-1}=\frac1{\sqrt2}
 \sum_{a\in\mathbb Z}\rho^{|a|}T_i^{2a},\qquad
\mathscr G^{-1}=\frac12\sum_{a,b\in\mathbb Z}
 \rho^{|a|+|b|}T_x^{2a}T_y^{2b}.}
\tag{32.17}
\]

Both are native Cauchy limits with positive coefficients of total weight
one. The inverse has an exponentially decreasing tail in counted two-edge
steps. On a periodic target the same absolutely controlled sums wrap around
the cyclic shifts and give the periodic inverse. No infinite non-periodic
kernel is silently substituted for a finite periodic boundary.

For a one-axis truncation |a|<=N the omitted coefficient weight is

\[
\tau_N=\frac{\sqrt2\,\rho^{N+1}}{1-\rho}.
\]

For the rectangular two-axis truncation the omitted weight is
2 tau_N-tau_N-squared. Consequently the reconstructed source error using
that truncated G-inverse and F-dagger is at most
2(2 tau_N-tau_N-squared)||d||. This gives a second finite decoder with an
explicit spatial count tail. Rho is the exact tail coefficient of this native target, not an
identification with alpha or a universal physical interaction strength.

**Proof.** Q_i=(3/2)I-(T_i^2+T_i^(-2))/4. Seek a candidate coefficient
c rho^|a| and verify it directly. The noncentral coefficient cancels when
rho-squared-6 rho+1=0. The two native cut roots are 3 plus/minus 2 sqrt(2);
the smaller is in (0,1), while the larger cannot give a summable tail.
The central coefficient is one for c=1/sqrt(2). Finite multiplication
therefore leaves only boundary terms at a=plus/minus N and
a=plus/minus(N+1). Their coefficients tend to zero with the proved
geometric bound. Thus the native completed sum is an inverse, and uniqueness
from R32.6 identifies it with Q_i^(-1).

The finite scalar identity c(1+rho)=1-rho gives total coefficient weight
one. Multiplying the two absolutely controlled one-axis sums proves the
two-axis formula. Summing the omitted geometric tails gives tau_N and
1-(1-tau_N)^2. Each shift preserves pairing, so these coefficient tails
bound operator error. Finally ||F-dagger v||<=2||v|| by (32.11), or by
B-script pairing preservation, proves the displayed reconstruction bound.
The exact replayer checks the quadratic and normalization identities in
the source census algebra; native rational square cuts separately isolate
the positive scalar root. The census operator itself is not declared a
positive square root on every source role.

## Certification and source lineage

Seven written proofs are bound to nine exact check groups, literal fine
source histories, complete Laurent identities, normal-order continuation,
sharp finite packet families, finite cyclic ranks, exact inverse/remainder
identities, native root cuts and symbolic replay. A declared dependency
graph rejects imported/admitted/open premises and promotion of targets to
physical selection. This checks metadata, not arbitrary mathematical prose.
Written proofs, exact computations, source integrity, formal verification
and physical validation remain separate evidence classes. The frozen R31
and earlier source chain replay without changes to their certificates or
to the canonical engine.

| Pinned source | Premise status and use |
| --- | --- |
| [R16 C3–C8 and C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native-derived:** signed/refined coefficients, roles, H/K, matching pairing and positive cut-root completion. Finite equality, records, counting and induction remain declared proof infrastructure. |
| [R20.1–R20.2 and tuple construction](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native-derived on an explicit target:** finite source/record tuples, controlled arrows and matching readouts. No Hilbert, tensor-for-nature or Born-probability premise is imported. |
| [R26.1–R26.4](https://github.com/Parveen117/extra-ideas/blob/6773afc38caff84ff0b71ed35e3aa8c3a4adbd45/02-relational-response/NATIVE_MEMORY_RESIDUE_R26.md) | **Native-derived on a stated record interface:** reciprocal frame cancellation, commutator response, retained signs and one-direction propagation gap/completion. R32 derives its different two-direction polynomial and mixed difference square. |
| [R27.1–R27.2](https://github.com/Parveen117/extra-ideas/blob/6035e66bd40fe4016eeded41a4d704ec4ff4b342/02-relational-response/NATIVE_REPLICA_SELECTION_R27.md) | **Native-derived source-law selection:** the balanced record census forces anticommutation and a minimum two-role factor. Physical choice of edge interface is not thereby selected. |
| [R29.1](https://github.com/Parveen117/extra-ideas/blob/16602a29c469ac96001977432fcf5fc4b548970b/02-relational-response/NATIVE_PAIRED_FIELD_DYNAMICS_R29.md) | **Native-derived:** two-event V/N factorization and matching preservation. R32 rederives its covariant version from literal source events and reciprocal record edges. |
| [R30.1 and R30.6](https://github.com/Parveen117/extra-ideas/blob/7ac5b3e5baba4dafb51f4330bfe84d001a9ab1d4/02-relational-response/NATIVE_DIRECTIONAL_LOOP_INTERACTION_R30.md) | **Native-derived on explicit targets:** two history-address labels and simultaneous role-frame covariance. Record-link curvature is distinguished from a full transport loop. |
| [R31.1–R31.7](https://github.com/Parveen117/extra-ideas/blob/63d040f031b6b36c36a4fa0475f5ce26ec765caf/02-relational-response/NATIVE_PROPAGATION_GEOMETRY_R31.md) | **Native-derived on the flat target:** local unmixing, propagation geometry and increment reconstruction modulo uniform source. R32 supplies a different retained interface and a precise obstruction to flat similarity. |
| [Canonical RKF algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native arithmetic/equality replay:** unchanged active engine, scoped application only. No whole-engine or proof-assistant PASS is asserted. |
| [R31 comparison ledger](https://github.com/Parveen117/extra-ideas/blob/63d040f031b6b36c36a4fa0475f5ce26ec765caf/04-operator-evolution/R31_DERIVATION_LEDGER.json) | **Excluded comparison premises:** declared topology/carriers and admitted classical smooth/tangent/metric adapters keep their labels and do not enter these proof paths. |

This development supplies a native retained loop with a measurable source
response, a derived volume-independent increment gap and internal count
cycle, and a complete localized inverse. Its staggered factorization and
its explicit edge/readout domains are part of the result. Physical field
selection, irreducible force identification, metric/clock calibration,
mass and c/alpha remain open; those claims cannot be obtained by renaming
the derived count coefficients.
