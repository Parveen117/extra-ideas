# R31: native propagation geometry and the complete increment observer

Research owner: Monty Dabas. Development: 1 October 2026.

R30 constructed a native loop reading that drives a closed paired-field
target. R31 resolves that drive into two independent source continuations
by an explicit local inverse. This refines the interpretation of R30:
its driven field is real within the stated readout, but its off-diagonal
term does not establish an irreducible new propagation law.

The alternating source itself has a richer directional structure. Its
wave operator has a mixed difference term, an exact positive native
factorization, a quartic transverse sector and only a rank-one leading
quadratic count form. On periodic count targets, the full stationary
sector is much smaller than the loop-blind sector. A single signed
four-event increment recovers every nonuniform source component, with
an exact native inverse and a certified finite-iteration error bound.

All directions, paired fields, finite boundaries and observations below
are explicitly constructed native targets. No classical metric, wave law,
Fourier modes, spectral theorem or physical clock is a premise. The native
matching pairing, finite equality/counting and scalar completion retain
their established source status. No claim of priority over all mathematics
or of formal proof-assistant verification is made.

## Source and notation

Reuse R30's H,K, R=KH, C_0=H+K, J_0=H-K, native shifts T_x,T_y,
delta_i=T_i-I and directional transports

\[
N_i=\tfrac12(P+T_i^{-1}Q)C_0,\qquad V_i=I+N_i\delta_i,
\qquad W=V_yV_x,\quad\Omega=[V_x,V_y].
\tag{31.0}
\]

Here P=(I+H)/2 and Q=(I-H)/2. The local source identities include
R-dagger=-R, R-squared=-I, C_0-squared=J_0-squared=2I,
RC_0=-J_0 and RJ_0=C_0. A dagger denotes the adjoint for the native
matching pairing; the inverse-shift rule follows by reindexing a finite sum.

Write

\[
\Delta_i=\delta_i^\dagger\delta_i=2I-T_i-T_i^{-1},\qquad
A_i^{\rm sk}=T_i-T_i^{-1}.
\tag{31.1}
\]

These are finite source-address differences, not assumed differential
operators. Initially fields have finite support on the unwrapped count
plane. Periodic counts and one-record directional quotients are specified
where used. Source roles are real signed native values in this target;
they are not a replacement definition of the full UGD state.

## Written results

### R31.1 — One local inverse unmixes the two directional loop drives

Retain R30's x update and construct the corresponding oriented y update:

\[
\mathcal C_x(b,\phi)=(V_xb+\delta_x\Omega\phi,V_x\phi),\qquad
\mathcal C_y(b,\phi)=(V_yb-\delta_y\Omega\phi,V_y\phi).
\tag{31.2}
\]

The y sign follows from [V_y,V_x]=-Omega. Its paired readout is
e_y=N_y b-Omega phi, just as the x readout is e_x=N_x b+Omega phi.
This two-direction protocol is an explicit extension of the R30 target.
It introduces no fitted interaction rate.

Define the finite operator and reversible change of variables

\[
S=\delta_xV_y+\delta_yV_x,\qquad B=b+S\phi,
\qquad
\mathcal D=\begin{pmatrix}I&S\\0&I\end{pmatrix},\quad
\mathcal D^{-1}=\begin{pmatrix}I&-S\\0&I\end{pmatrix}.
\tag{31.3}
\]

Then the same local change works for both directions:

\[
\boxed{\mathcal D\mathcal C_i\mathcal D^{-1}
       =\operatorname{diag}(V_i,V_i),\qquad i=x,y.}
\tag{31.4}
\]

Consequently the alternating interaction block
\(\mathcal C=\mathcal C_y\mathcal C_x\) is conjugate to diag(W,W),
and its directional order residue is conjugate to diag(Omega,Omega).
No input information is removed. There are two freely prepared source
copies, with the same propagation polynomial; the displayed drive does
not by itself supply an additional irreducible transport mechanism.

**Proof.** Scalar shifts commute with both directional coefficients, so

\[
[V_x,S]=\delta_x\Omega,\qquad
[V_y,S]=-\delta_y\Omega.
\]

Substitute b=B-S phi into (31.2) and these commutators cancel its drive.
The new variables satisfy B'=V_i B and phi'=V_i phi. Direct triangular
multiplication proves the two-sided inverse in (31.3), (31.4), the block
statement and the commutator statement. R30's x update alone also admits
the smaller local change B_x=b+delta_x V_y phi. None of these equalities
identifies a change of observable coordinates with a physical gauge law.

### R31.2 — The full retained source has a positive conserved energy and bounded drive

For either direction and for any sequence of them,

\[
\boxed{\mathcal E_{\rm full}(b,\phi)=\|b+S\phi\|^2+\|\phi\|^2}
\tag{31.5}
\]

is conserved and positive definite before taking an observation quotient.
It retains even the phi values invisible to R30's loop readout. It is
distinct from R30's pulled-back paired-field energy: that form remains
valid for its original x evolution and observation, whereas (31.5) is
proved common to the new two-direction protocol.

For n alternating blocks, the exact solution is

\[
\boxed{\phi_n=W^n\phi_0,\qquad
b_n=W^nb_0+(W^nS-SW^n)\phi_0.}
\tag{31.6}
\]

In particular,

\[
\|b_n-W^nb_0\|\le8\|\phi_0\|\quad(n\ge0).
\tag{31.7}
\]

This is a uniform bound, not a claim of its sharp value. For R30's
x-only protocol the smaller S_x gives coefficient four instead of eight.
An apparent source term cannot accumulate without bound in this target.

**Proof.** R30 proves each V_i preserves native pairing. Apply (31.4)
to (B,phi); both squared norms in (31.5) are separately preserved.
The form vanishes only if phi=0 and then b=0. Continuing B and phi and
substituting b=B-S phi gives (31.6). Since T_i preserves pairing,
the native triangle inequality gives ||delta_i v||<=2||v||. Therefore
||S v||<=4||v||, and the two terms of its commutator with W^n give (31.7).
These inequalities are native matching-square consequences, not a supplied
classical field-energy rule or a physical energy/mass identification.

### R31.3 — Alternating source events derive a mixed wave operator and a finite count cone

The source yields the scalar finite operator

\[
\boxed{\mathscr L=\tfrac12(\Delta_x+\Delta_y)
       -\tfrac14 A_x^{\rm sk}A_y^{\rm sk},\qquad
       W+W^\dagger=2I-\mathscr L.}
\tag{31.8}
\]

Thus each continued source satisfies

\[
\boxed{\psi_{n+1}-2\psi_n+\psi_{n-1}
       =-\mathscr L\psi_n.}
\tag{31.9}
\]

The same polynomial holds for the entire coupled block C, its b and phi
components, and the dual continuation of the curvature reading. The mixed
term is source-derived and cannot be dropped while retaining this identity.

A compact source initially supported on K stays within
K+[-n,n]^2 after n four-event blocks. The joint (b,phi) update stays within
K+[-n-1,n+1]^2 for initially joint support K. The extra fixed collar comes
from S; it is not an ever-growing independent propagation rate. A source
unit packet attains address (n,n), so the source outer count bound is attained.

**Proof.** Put d_i=Delta_i/4 and a_i=A_i^sk/4. Source multiplication gives

\[
V_i=(I-d_i)I+a_i C_0+d_iR.
\]

The scalar coefficient of V_y V_x is
1-(Delta_x+Delta_y)/4+A_x^sk A_y^sk/8. Its remaining three source-role
terms are skew under the native dagger. Adding W-dagger proves (31.8).
Multiplication by W uses W-dagger W=I and proves
W-squared-(2I-L)W+I=0, hence (31.9). L is a scalar shift polynomial,
so it commutes with S and (31.4) proves the coupled identity. Reversing
the direction order leaves the scalar coefficient unchanged; taking its
dagger gives R30's dual wave identity as well.

Each V_i shifts its own address by at most one. Products give the source
support bound by finite induction. The stencil of S lies in [-1,1]^2;
the two terms in (31.6) give the joint bound. For an initial e0 packet,
the (1,1) coefficient of W is (P C_0)^2/4=P C_0/4, and
P C_0 e0=e0. Hence the (n,n) value is 4^(-n)e0 and cannot cancel.
The count cone has no assigned physical length or duration, so it is not
yet a value of c in physical units.

### R31.4 — Positive native squares select exactly the stationary uniform sector

Construct

\[
\mathscr D=\tfrac14(\Delta_x+\Delta_y),\qquad
\mathscr Q=\tfrac14\bigl[(\Delta_x-\Delta_y)R
 +(A_x^{\rm sk}+A_y^{\rm sk})J_0\bigr].
\tag{31.10}
\]

Then Q-script is skew under dagger, D-script commutes with it, and

\[
\boxed{R^\dagger(V_x-V_y^\dagger)=\mathscr D+\mathscr Q,
\qquad \mathscr L=\mathscr D^2+\mathscr Q^\dagger\mathscr Q.}
\tag{31.11}
\]

The propagation reading is therefore an exact sum of native squares:

\[
\boxed{\|(W-I)\psi\|^2
 =\|\mathscr D\psi\|^2+\|\mathscr Q\psi\|^2.}
\tag{31.12}
\]

On the finite-pairing domains stated here, its stationary kernel consists
exactly of fields invariant under both address shifts. On the periodic
L-by-L target it has two native roles:
one uniform H/K value. On the unwrapped finite-support target it is zero.
The joint interaction has four stationary uniform roles on the periodic
target, by its reversible unmixing. The uniform statement is about finite
periodic norm sums; it assigns no norm to a nonzero infinite background.

**Proof.** Insert the V_i decomposition from R31.3. RC_0=-J_0 and
R-squared=-I give the first identity in (31.11). Every scalar Delta is
self-dagger and every A^sk is skew, proving Q-dagger=-Q. The scalar D
commutes with Q, so the cross terms in (D+Q)-dagger(D+Q) cancel.
Furthermore V_x-V_y-dagger=V_y-dagger(W-I), and both R and V_y
preserve pairing. This proves the second identity and (31.12).

If (W-I)psi=0, (31.12) forces D psi=0. Taking its pairing with psi gives

\[
0=\langle\psi,\mathscr D\psi\rangle
 =\tfrac14(\|\delta_x\psi\|^2+\|\delta_y\psi\|^2).
\]

Both differences vanish, hence the source is uniform on the connected
count target. The converse follows from (31.0). A compact uniform field
on the unwrapped plane is zero. For the joint statement apply (31.4);
S annihilates constants, so its four constant roles are exactly uniform
(b,phi). This argument uses no spectral decomposition.

The finite-pairing domain is essential: the unbounded affine preparation
psi(x,y)=(x-y)v is also stationary under W, although neither individual
address difference is zero. Its constant first differences cancel between
V_x and V_y. It has no finite global matching norm and is not part of the
stationary-kernel classification above.

An additional native identity supplies the upper bound used below:

\[
\mathscr L-\tfrac14\mathscr L^2
 =\left(\frac{W-W^\dagger}{2}\right)^\dagger
  \left(\frac{W-W^\dagger}{2}\right)\succeq0,
\qquad 0\preceq\mathscr L\preceq4I.
\tag{31.13}
\]

Multiply out using W-dagger W=I and (31.8) for the first equality.
For the second, ||(W-I)psi||<=2||psi|| proves the upper bound directly.

### R31.5 — Opposite count directions have an exact quartic wave and a degenerate quadratic jet

Identify the two address shifts with one native shift T in three stated
readouts, and put Delta=2I-T-T^(-1). Exact substitution in (31.8) gives

| Count target | Exact wave operator L-script |
| --- | --- |
| T_x=T, T_y=I | Delta/2 |
| T_x=T, T_y=T | 2 Delta-Delta-squared/4 |
| T_x=T, T_y=T^(-1) | Delta-squared/4 |

The third case is a fourth-order finite difference, without taking any
continuum limit. These are quotient readouts of the count algebra; their
existence does not assert a nonzero compact stripe on the unwrapped plane.
The associated curvature reductions are

\[
\Omega\big|_{T_y=T_x}=0,\qquad
\boxed{\Omega\big|_{T_y=T_x^{-1}}
 =\tfrac14\Delta(T-T^{-1})J_0.}
\tag{31.14}
\]

The leading native count coefficient form is also exact as a finite jet.
Let p=T_x-I and q=T_y-I, and discard monomials of total degree at least
three. Finite multiplication yields

\[
\mathscr L\equiv-\tfrac12(p+q)^2\pmod{(p,q)^3},\qquad
W-I\equiv\tfrac12 C_0(p+q)\pmod{(p,q)^2}.
\tag{31.15}
\]

The matching-square coefficient on these first differences is consequently

\[
\boxed{G_2=\tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix},
\quad\operatorname{rank}G_2=1,\quad G_2(1,-1)^t=0.}
\tag{31.16}
\]

Two history-address labels therefore do not force a nondegenerate
two-direction quadratic metric. This source target has a transverse
quartic response beyond its rank-one first-difference reading.

**Proof.** For equal or inverse shifts substitute directly in (31.8).
The native identity (T-T^(-1))^2=Delta^2-4 Delta gives all three rows;
substitution in R30.2 gives (31.14). For the finite jet,
(1+p)(1-p+p^2)=1+p^3, so T_x^(-1)=1-p+p^2 modulo degree three,
and likewise for q. Apply this finite equality to (31.8) and (31.0)
to get (31.15); no analytic Taylor theorem is an input. Finally
C_0-dagger C_0=2I gives (31.16). The negative sign in the scalar shift
jet is consistent with p-dagger=-p modulo degree two; G_2 is the positive
matching-square form, not a supplied physical metric tensor.

### R31.6 — Native finite averages give an explicit propagation gap bound

On the periodic count target choose L>=2. Reuse the finite averages

\[
A_i=\frac1L\sum_{k=0}^{L-1}T_i^k,\quad
B_i=\frac1L\sum_{k=0}^{L-1}kT_i^k,\quad
\Pi_0=A_xA_y,\quad\psi_0=(I-\Pi_0)\psi.
\]

Pi_0 is the uniform two-role average. With
\(\mu_L=(L-1)^{-4}\), the native propagation reading obeys

\[
\boxed{\mu_L\|\psi_0\|^2
 \le\|(W-I)\psi_0\|^2\le4\|\psi_0\|^2.}
\tag{31.17}
\]

The lower bound is explicit and sufficient; it is not claimed sharp.
It concerns the selected finite count boundary and is not a particle mass
gap or a volume-independent physical constant.

**Proof.** R30's finite coefficient subtraction gives
delta_i B_i=I-A_i. Weighted shifts preserve pairing, so
||B_i v||<=((L-1)/2)||v||. The two components in

\[
\psi_0=(I-A_x)\psi+A_x(I-A_y)\psi
\]

are orthogonal native cuts. Their squared norms and the count bound give

\[
\|\psi_0\|^2\le\tfrac14(L-1)^2
 (\|\delta_x\psi\|^2+\|\delta_y\psi\|^2).
\]

Apply this to psi_0 and use (31.10):
<psi_0,D psi_0> >= ||psi_0||^2/(L-1)^2. Native Cauchy--Schwarz,
derived by the matching-square argument in the source, then gives
||D psi_0|| >= ||psi_0||/(L-1)^2. Equation (31.12) proves the lower
bound in (31.17), while (31.13) proves the upper one. All averages and
weighted sums are finite; no assumed trigonometric eigenbasis is used.

### R31.7 — A single signed block increment is a complete observer modulo uniform source

Observe the full signed address field

\[
d=F\psi,\qquad F=W-I.
\tag{31.18}
\]

This is the change over one internally counted four-event block, including
both native roles at each retained address. It is not a single scalar
sensor reading or an intensity-only measurement. R31.4 implies

\[
\boxed{\ker F=\operatorname{im}\Pi_0,\quad
\operatorname{rank}F=2L^2-2,\qquad
\psi_0=(\mathscr L+\Pi_0)^{-1}F^\dagger d.}
\tag{31.19}
\]

The inverse is a finite native elimination, available for every L>=1;
L=1 has psi_0=0. It recovers the whole nonuniform source, including the
6(L-1) nonuniform directions invisible to R30's loop observation. A
linear observer reconstructing this full quotient needs at least 2L^2-2
independent signed scalar channels; native row elimination of F attains
that count. This is a target-specific observation rank, not a universal
number of physical observers.

There is also an explicit finite decoder with an error bound. For L>=2 put

\[
\mathcal A=I-\mathscr L/4,\qquad
K_N=\tfrac14\sum_{j=0}^{N}\mathcal A^jF^\dagger,
\qquad q_L^2=1-\mu_L/4.
\tag{31.20}
\]

Then

\[
\boxed{\|K_Nd-\psi_0\|
 \le q_L^{N+1}\|\psi_0\|
 \le (L-1)^2q_L^{N+1}\|d\|.}
\tag{31.21}
\]

The bound is derived from the reading itself, without a supplied unknown
source norm. It does not make a generic finite truncation an exact inverse.

**Proof.** R31.4 proves the kernel. F-dagger F=L-script by (31.12),
and both F and L annihilate Pi_0. The native form of L+Pi_0 is strictly
positive: its zero norm would force a uniform source with zero uniform
average. Finite native elimination therefore has a unique inverse.
Multiplying gives (L+Pi_0)^(-1)F-dagger F=I-Pi_0, proving (31.19).
The source coordinate count gives the rank and the lower bound on linear
observer channels; selecting independent native rows attains it.

For (31.21), all scalar shift operators commute. The finite geometric
identity gives

\[
K_NF=I-\mathcal A^{N+1}
 =(I-\mathcal A^{N+1})(I-\Pi_0).
\]

Using L-squared<=4L from (31.13), for every mean-free v,

\[
\|\mathcal Av\|^2
 =\|v\|^2-\tfrac12\langle v,\mathscr Lv\rangle
  +\tfrac1{16}\|\mathscr Lv\|^2
 \le(1-\mu_L/4)\|v\|^2.
\]

Induction bounds the finite remainder, and (31.17) replaces the unknown
source norm by ||d||/sqrt(mu_L). This proves the claimed error estimate
without a classical spectral theorem or an imported inverse transform.

Finally R30 gives dim ker(Omega)=2(3L-2) and
Omega W=widehat(W) Omega. Hence every loop-blind source stays blind under
all future repetitions of that same loop observation. Only two of those
roles are stationary, by R31.4. Subtracting gives 6(L-1) dynamic blind
directions recovered by the new signed increment. Waiting longer with
the old observer is therefore different from adding the new readout.

## Certification and source lineage

Seven written proofs are bound to nine exact check groups, complete
native Laurent identities, finite directional jets, finite source
continuation, exact positive-square witnesses, count ranks/inverses,
decoder remainder identities and native word replay. The declared
dependency graph rejects imported/admitted/open premises and promotion
of target constructions to physical selection. It audits declarations,
not arbitrary proof semantics. Written proof, exact computational checks,
provenance, formal verification and physical validation remain distinct.
The frozen R30 and earlier source chain are replayed without changing
their certificates or the canonical engine.

| Pinned source | Premise status and use |
| --- | --- |
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native-derived:** signed/refined coefficients, roles, H/K, matching pairing, its square inequalities and scalar completion. Finite equality, counting and induction remain declared proof infrastructure. |
| [R18.1–R18.5](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | **Native-derived on an explicit target:** event-count address and normalized signed transport. No physical clock, metric or length unit is supplied. |
| [R29.1–R29.5](https://github.com/Parveen117/extra-ideas/blob/16602a29c469ac96001977432fcf5fc4b548970b/02-relational-response/NATIVE_PAIRED_FIELD_DYNAMICS_R29.md) | **Native-derived on stated readouts:** V/N, paired induction, compatibility residue and positive field form. R31 preserves this x-only result and derives a separate form for the bidirectional extension. |
| [R30.1–R30.8](https://github.com/Parveen117/extra-ideas/blob/7ac5b3e5baba4dafb51f4330bfe84d001a9ab1d4/02-relational-response/NATIVE_DIRECTIONAL_LOOP_INTERACTION_R30.md) | **Native-derived on explicit direction/loop/interaction/boundary targets:** shared roles, exact curvature, dual transport and cyclic blindness. R31 proves a local unmixing and a stronger increment observer; the old physical-force/EM boundary is not promoted. |
| [R23.1–R23.9](https://github.com/Parveen117/extra-ideas/blob/dc985839f1c3d1ccbc29128b7aa76a06cc7a72ed/02-relational-response/FUTURE_RESPONSE_QUOTIENT_R23.md) | **Native prior attribution only:** finite observer ranks and complete response quotients. R31.7 provides its own proof for the increment target; this comparison does not import a physical observer count. |
| [Canonical RKF algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native arithmetic/equality replay:** unchanged canonical engine, scoped exact use. No whole-engine or proof-assistant PASS is claimed. |
| [R30 premise ledger](https://github.com/Parveen117/extra-ideas/blob/7ac5b3e5baba4dafb51f4330bfe84d001a9ab1d4/04-operator-evolution/R30_DERIVATION_LEDGER.json) | **Excluded comparison premises:** declared topology/carrier models and admitted classical smooth/tangent/metric adapters keep those labels and are not proof inputs here. |

The completed advance is a derived propagation operator, its native
positive factorization and directional degeneracy, and a complete source
observer beyond loop blindness. The precise unmixing theorem identifies
what the constructed interaction can establish. Physical irreducible
coupling, metric/dimension selection, EM, c, alpha and particle mass remain
separate open targets; no count gap or readout coefficient is renamed one
of those constants.
