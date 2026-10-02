# R39: a native local reverser, echo clock and distance reading

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R38 supplies the sharp reliable-signal cone of the retained source. R39
constructs an actual **local reversal operation** from that source's H/K
factors. It obeys

\[
\boxed{\mathcal C^\dagger=\mathcal C,\quad
\mathcal C^2=I,\quad \mathcal C Z\mathcal C=Z^\dagger.}
\tag{39.A}
\]

Consequently a pulse continued for a native count, reversed, and
continued for the same count returns its full signed information. A
second reversal restores the preparation exactly. The reversal can be
confined to a finite set of native links without losing unitarity.
Its localization error is bounded by the information not captured there.

With R37 packets, this constructs a high-retention remote echo, a return
tick, an error-controlled midpoint timestamp and a distance reading.
Native record counting supplies the memory needed to distinguish repeated
ticks. This is a derived **controlled source protocol**: the two reversal
controls and their record are explicit. It is not claimed to be an
autonomously selected material oscillator, a physical mirror, or the
unique physical choice of rods and clocks.

## Unchanged source and native control target

Use the pinned R32/R35 source

\[
X=T_xH_m,\quad Y=T_yK_m,\quad
V_A=I+\tfrac12(P_s+A^{-1}Q_s)(H_s+K_s)(A-I),
\qquad B=(V_YV_X)^2.
\tag{39.0}
\]

Here `P_s=(I+H_s)/2`, `Q_s=I-P_s`. The R35 signed address map and
parity regrouping give the sixteen-role envelope `Z=-chi B chi`.
All matching, scalar roots, tuples, shifts and finite word histories are
the already earned native R16/R20 constructions.

One Z block contains eight original source events. Reversal controls
are separately recorded arrows; they are not silently counted as free
ordinary Z events or assigned a material duration. The control below has
the explicit finite local composition `U (onsite cut operation) U^dagger`.
The protocol's flight counter counts Z blocks, and its control ledger
records those additional operations.

## Written results

### R39.1 — The source factors derive an exact local reversal involution

Define the native normalized cut contrast

\[
D=(H_s-K_s)/\sqrt2,\qquad D^\dagger=D,\quad D^2=I.
\]

For every unitary record/address arrow A commuting with the source-role
cuts, define `U_A=A P_s+Q_s`. Direct native cut multiplication gives

\[
\boxed{DV_AD=V_A^\dagger,\qquad
       V_AD=U_A D U_A^\dagger.}
\tag{39.1}
\]

Therefore the coarse reverser `C_B=V_YD` is a self-dagger involution and
`C_B B C_B=B^dagger`. In the signed, regrouped source let

\[
\mathcal U=\chi U_Y\chi,\qquad
\mathcal C=\chi V_YD\chi
             =\mathcal U D_{16}\mathcal U^\dagger.
\tag{39.2}
\]

It satisfies (39.A). Its component displacements are only `(0,0)`,
`(0,1)` and `(0,-1)`; U has only `(0,0)` and `(0,1)`.
Thus this is a local source/record operation, not a reflection of the
whole address plane.

**Proof.** `H_s^2=K_s^2=I` and their anticommutation give D's square
and dagger. The already earned positive root fixes its normalization.
Writing the source in its own plus/minus cuts, (39.0) becomes

\[
V_A=\frac12
\begin{pmatrix}I+A&A-I\\I-A^{-1}&I+A^{-1}\end{pmatrix},
\qquad
V_AD=\frac1{\sqrt2}
\begin{pmatrix}I&-A\\-A^{-1}&-I\end{pmatrix}.
\tag{39.3}
\]

These blocks are derived cut arrows. The second matrix is self-dagger
and squares to identity. Its equality with
`diag(A,I) D diag(A^dagger,I)` proves (39.1); self-daggerness also gives
`DV_AD=V_A^dagger`.

Put A=V_Y and F=V_X for this one calculation. Since the same D reverses
both, `(AD)(AF)(AD)=F^dagger A^dagger`. Thus C_B reverses V_Y V_X and
its square B. Conjugating by the source's own chi and regrouping proves
(39.A), including the minus sign in Z. The explicit conditional link
`U_Y=Y P_s+Q_s` moves only along y. Its regrouped coefficients have the
two stated displacements. Multiplying U, D_16 and U-dagger gives the
stated three-coefficient support. The certificate checks the complete
Laurent identities, not only their first jets or selected phases.

### R39.2 — Echo cancellation and finite-window reversal are exact native identities

For any nonnegative a,b,

\[
\boxed{\mathcal C Z^b\mathcal C Z^a=Z^{a-b}.}
\tag{39.4}
\]

Equal counts return every signed source field exactly. Unequal counts
leave precisely the unmatched source continuation; they are not erased.

For a finite set of component addresses T let P_T retain every role at
those addresses. Define the native paired-link cut and local gate

\[
Q_T=\mathcal U P_T\mathcal U^\dagger,\qquad
\boxed{\mathcal C_T=\mathcal C Q_T+I-Q_T
 =\mathcal U[D_{16}P_T+I-P_T]\mathcal U^\dagger.}
\tag{39.5}
\]

Then Q_T is an orthogonal cut commuting with C, and C_T is a
self-dagger unitary involution. It is identity outside the finite link
patch reached from T by U. For every field u,

\[
\boxed{\|(\mathcal C_T-\mathcal C)u\|
                      \le2\|(I-Q_T)u\|.}
\tag{39.6}
\]

**Proof.** The exact conjugation in (39.A) gives
`C Z^b C=Z^(-b)`, proving (39.4). U preserves matching and P_T is a
native address cut, so Q_T is an orthogonal cut. P_T commutes with the
onsite D_16. Hence Q_T commutes with C and the two orthogonal summands
in (39.5) respectively apply C and identity. This proves the gate's
dagger, inverse and square. Only the paired channels in the image of
U P_T are changed; the gate has finite local support.

Subtracting the full reverser gives
`C_T-C=(I-C)(I-Q_T)`. The native triangle bound `||I-C||<=2` proves
(39.6). No source amplitudes are thrown away. Simply replacing C by
`P_T C P_T+I-P_T` generally loses norm at the patch boundary; the
source-derived link cut in (39.5) is the necessary distinction.

### R39.3 — A local high-retention trigger produces a controlled remote echo

Take a fixed R37 carrier with drift 0<g<1 and its corrected compact
packet u_0. Set `M=L^2`, `N=L^3`, `W=3M-2`,
`d=floor(gN)` and `epsilon=epsilon_(N,M)` from R37. Its source bound
holds uniformly for 0<=n<=N and its initial support is inside `[0,W]^2`.

Let

\[
T_0=[0,W]\times[-1,W+1],\qquad T_d=T_0+(d,0),
\]
\[
\Delta_d=\lceil(W+3)/g\rceil,\qquad
\Delta_0=\lceil(W+1)/g\rceil.
\tag{39.7}
\]

Use Q_0,Q_d and C_0,C_d for the cuts and gates of these two patches.
Assume the explicit, eventually satisfied budgets
`epsilon^2<1/32` and `N>Delta_d+Delta_0`.

The first source count tau at which

\[
\|Q_d Z^\tau u_0\|^2\ge(1-2\epsilon^2)\|u_0\|^2
\]

exists and satisfies `N-Delta_d<tau<=N`. Applying C_d there gives a
state within `3 epsilon ||u_0||` of the ideal reversed state. After
another tau source blocks and the origin gate,

\[
\boxed{\|\mathcal C_0 Z^\tau\mathcal C_d Z^\tau u_0-u_0\|
                                  \le3\epsilon\|u_0\|.}
\tag{39.8}
\]

**Proof.** U and U-dagger do not move the x address. Their y range is
at most one. Therefore `Q_0u_0=u_0`. At N the translated R37 comparison
packet has shift d and is fully captured by Q_d. Its actual missing norm
is at most epsilon times the conserved source norm. Hence the capture
fraction is at least `1-epsilon^2`, above the chosen trigger.

For n<=N-Delta_d the comparison packet's x support is disjoint from
the remote patch: `floor(gN)-floor(gn)>=g Delta_d-1>W`. Its Q_d
fraction is at most `epsilon^2`, below the trigger. Finite earliest-event
selection gives the claimed tau interval. At the trigger, the missing
relative norm is at most sqrt(2) epsilon. Equation (39.6) makes the gate
error at most `2sqrt(2) epsilon<3epsilon`.

Subsequent source continuation preserves this error. The ideal identity
is `Z^tau C Z^tau u_0=C u_0`. Since Q_0 commutes with C and captures
u_0, `C_0 C u_0=u_0`. The origin gate is an isometry, so the error
remains at most 3epsilon, proving (39.8). The trigger is a declared
native norm readout controlling an available reversible arrow; no
physical measurement-update law is assumed or inferred.

### R39.4 — Return detection derives a midpoint timestamp and a distance interval

After the remote gate, let s_* be the first nonnegative source count
with at least half the prepared norm square back in Q_0. Let
`r=tau+s_*` be the recorded round-trip flight count. Then

\[
\boxed{\tau-\Delta_0<s_*\le\tau,\qquad
       2\tau-\Delta_0<r\le2\tau,}
\tag{39.9}
\]
\[
\boxed{0\le\tau-r/2<\Delta_0/2,}
\tag{39.10}
\]
\[
\boxed{-1\le d-gr/2<g(\Delta_d+\Delta_0/2).}
\tag{39.11}
\]

Thus r/2 estimates the remote trigger's source count, and gr/2 reads
the link-patch displacement with a derived finite uncertainty. Both
errors are O(M); their relative errors vanish when M=L^2 and N=L^3.

**Proof.** At return count s the actual state differs by at most
3epsilon in relative norm from
`Z^s C Z^tau u_0=C Z^(tau-s)u_0`, for 0<=s<=tau.
Because Q_0 commutes with C, its ideal captured norm is exactly that of
the outbound packet at count t=tau-s. If t>=Delta_0, then
`floor(gt)>=W+1`, so the R37 comparison packet is disjoint from Q_0.
The actual returned capture norm is at most `epsilon+3epsilon=4epsilon`.
Its fraction is at most `16epsilon^2<1/2`.

At s=tau the ideal state is C u_0, fully in Q_0. The actual outside
norm is at most 3epsilon, while the total norm is conserved. Its capture
fraction is at least `1-9epsilon^2>1/2`. Earliest-event selection proves
(39.9), and rearrangement gives (39.10).

Since `N-Delta_d<tau<=N`, (39.9) gives
`0<=N-r/2<Delta_d+Delta_0/2`. Insert
`-1<=floor(gN)-gN<=0` to obtain (39.11). The midpoint is a consequence
of the derived reciprocal continuation, not an imported synchronization
postulate. This is a timestamp in the common source-count protocol;
independent material-clock synchronization is not asserted.

### R39.5 — A repeatable echo clock has a derived error budget and a minimum tick memory

Fix a calibrated tau from R39.3 and use the same finite patches. The
scheduled cycle

\[
\mathcal E_\tau=\mathcal C_0 Z^\tau\mathcal C_d Z^\tau
\]

is unitary, and for every integer k>=0,

\[
\boxed{\|\mathcal E_\tau^k u_0-u_0\|
                         \le3k\epsilon\|u_0\|.}
\tag{39.12}
\]

Choosing, for example, `k=floor(sqrt(L))` gives an increasing number
of ticks with vanishing total relative state error. A cycle contains
exactly 2tau free source blocks and two separately retained reversal
controls. The origin may record its first return tick earlier, as in
(39.9); the scheduled reset is applied at 2tau, not silently at that
earlier threshold crossing.

An observer required to distinguish tick counts 0,...,k of an exact
return cycle needs at least k+1 distinct retained record labels. A native
cyclic successor on q>=k+1 labels achieves this bound on that horizon.
For binary marked-role storage the least p with `2^p>=k+1` suffices.

**Proof.** R39.3 gives the one-cycle defect. Finite telescoping gives
`E_tau^k-I=sum_(j=0)^(k-1) E_tau^j(E_tau-I)`. Unitarity bounds the
sum on u_0 by k times the one-cycle defect. R37 gives epsilon=O(1/L)
at a fixed carrier; hence the displayed tick family has error O(L^(-1/2)).
This repeats a fixed scheduled protocol. It does not assume that an
adaptive first-trigger controller has a linear cycle operator.

For the memory statement use the exact full reversal cycle in (39.4).
Its source field is identical after every tick, including signed phase.
A readout of that field alone assigns the same answer at all ticks.
Faithfully distinguishing k+1 counts therefore needs k+1 distinct
additional history labels. Construct q native marked count labels and
the permutation j->j+1 modulo q. It preserves matching and has inverse
j->j-1. Tensoring this available record arrow with the cycle gives the
claimed counter without copying an unknown source state. The finite
binary tuple count gives the p bound. After q ticks a q-label counter
aliases; retaining further winding history or enlarging the record is
necessary. Neither finite source recurrence nor phase alone stores an
unbounded elapsed tick count.
This is a capacity result for distinguishable marked-label counters;
no minimum dimension for arbitrary coherent phase encodings is asserted.

### R39.6 — Echo flight calibrates native information length and retains the sharp speed

For each fixed carrier, the constructed echoes obey

\[
\boxed{\frac{2d}{r}\longrightarrow g.}
\tag{39.13}
\]

Along the controlled R37 family `g->c_Sigma=1/sqrt(2)`, choose packet
scales to make the errors in (39.8)–(39.11) vanish relatively. Then

\[
\boxed{\frac{2d}{r}\longrightarrow c_\Sigma,\qquad
       \frac{r/2}{\sqrt2\,d}\longrightarrow1.}
\tag{39.14}
\]

Consequently half the recorded round-trip flight count becomes a native
information-length reading. R36's response-dual length, sqrt(2)d, and
the operational echo length agree in this limit. At fixed nonzero
carrier the ratio instead tends to `1/(sqrt(2)g)`; that derived
dispersion correction must be retained.

**Proof.** Equation (39.11) has error O(M) and d is of order N, with
M/N=1/L. This proves (39.13). R37's native rational family gives
`2/3,12/17,408/577,... ->1/sqrt(2)`. For each carrier first enlarge
L to control its own constants, then take a joint sequence. Substitution
gives (39.14) using R36's derived metric normalization, without assuming
a physical light-travel distance law.

The locality of the control matters. Q_d is contained in a patch of
diameter O(M)=o(d) at distance d, and C_d changes only adjacent retained
links there. R38 bounds the outbound fixed-fraction arrival by c_Sigma.
The high-retention state after reversal is again localized at the remote
patch up to vanishing norm, so R38 bounds its return leg as well. Thus
such two-leg signals cannot have limiting `2d/r>c_Sigma`; (39.14)
attains the supremum. This statement requires both legs to retain a
fixed positive fraction and local controls of negligible spatial extent
relative to d.

To check R38's localization contract on either leg, finite source
support first makes its free flight count at least a fixed positive
multiple of d whenever a fixed positive fraction crosses the separation.
A starting tail of vanishing norm cannot supply that fraction. Thus an
aperture of size o(d) is also o(the leg count), as R38 requires.

The history ledger is `(2tau free blocks, two reversal controls)` for a
reset cycle; the free part contains 16tau original events. No duration
for the control arrows is invented. If a separately derived control
implementation has uniformly bounded event overhead per reversal,
adding that recorded overhead leaves (39.13)–(39.14) unchanged. This
last accounting statement is not a proof of a physical implementation.

### R39.7 — Curvature-complete observation reads the echo and exposes the retained-link correction

Let O,D_o be the exact R36 observation and decoding maps. Transport the
source, reverser and cuts together:

\[
Z_o=OZD_o,\qquad C_o=O\mathcal C D_o,\qquad Q_{T,o}=OQ_TD_o.
\]

The same local echo, trigger fractions, state errors and clock-memory
requirements hold in this complete observer. At the zero carrier let
`C_* = C(1,1)` and retain R36's native curvature turn `J=2A_xA_y`.
Then

\[
C_*A_iC_*=-A_i,\quad [C_*,J]=0,\qquad
T_*=C_*J,\quad T_*^2=-I,\quad [T_*,A_i]=0,
\]
\[
\boxed{C_*=-T_*J.}
\tag{39.15}
\]

The exact reverser is supplied by the full retained source links.
R36's bare J is still not an exact finite-source time reversal.

**Proof.** O and D_o are mutual isometries, so adjacent inverse factors
cancel in every source, echo and cut product. All matching fractions
and error inequalities transport unchanged, as do simultaneously
transported pointwise native frames. Dropping a cut readout loses this
complete reconstruction and is a different observation target.

Differentiate `C(z)Z(z)C(z)=Z(z)^dagger` at the zero carrier, where
Z=I and C^2=I. The derivatives of the two C factors cancel by the
derivative of C^2; therefore C_* reverses both A_i. It then commutes
with their product J. Since C_* is a self-dagger involution and J is
a skew-dagger turn, their commuting product T_* has square -I.
Both factors reverse A_i, so their product commutes with it. This proves
(39.15). The complete Laurent defect `JZJ^dagger-Z^dagger` remains
nonzero, as in R36.6: the present retained-link construction does not
erase that earlier counterexample.

These are native metrological protocols and available source-derived
control arrows. They select no material field, physical reflection law,
unique autonomous clock or SI duration/length. Physical c, h and alpha
are not obtained by naming the count readings after those constants.
Finite-control costs and missing record information remain explicit.

## Certification and pinned source reuse

The [native application](../04-operator-evolution/native_echo_clock.cjs)
checks the complete local reversal and echo identities, paired-link
window gates, norm retention and clipping errors, exact count budgets,
record capacity, curvature transport and source/front boundaries.
Positive sqrt(2) belongs to R16's completion. Rational coefficient
identities are checked directly; finite local-gate witnesses use the
native matrix realization of a root satisfying `t^2=1/2`. This checks
both algebraic conjugates and is labelled as such, not as an empirical
simulation of a material clock.

The [verification](../04-operator-evolution/R39_VERIFICATION.json) binds
seven written proofs and scoped exact checks, replays frozen R38 and
preserves all 299 earlier non-navigation files. The
[ledger](../04-operator-evolution/R39_DERIVATION_LEDGER.json) labels the
controlled protocol, source laws and physical-selection boundary.
Metadata auditing and finite PASS are not semantic proof-assistant
verification or external physical validation.

| Pinned source | Premise status and use |
|---|---|
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** H/K, matching, positive roots, source histories and their distinction from evaluated states. Finite counting/induction remain explicit infrastructure. |
| [R20.1 and native tuple construction](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native construction:** available reversible marked-record arrows, not a physical measurement or allocation law. |
| [R32 retained source](https://github.com/Parveen117/extra-ideas/blob/1105a43209b144bd34a81bfbb99acc7876948e95/02-relational-response/NATIVE_RETAINED_LOOP_GAP_R32.md) | **Native-derived within its interface:** actual factorized record continuation and edge links. |
| [R35.1–R35.7](https://github.com/Parveen117/extra-ideas/blob/9c13e1752263c756d356dcf5851ac873d2e0d102/02-relational-response/NATIVE_SIGNED_ENVELOPE_R35.md) | **Native-derived:** signed regrouping, principal operators, source-event conversion and surviving faster front. |
| [R36.1, R36.3–R36.7](https://github.com/Parveen117/extra-ideas/blob/011b8224c92d95f7f726eb56bab3e8d8e6f7280c/02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md) | **Native-derived:** complete observer, information length and bare-curvature reversal defect. No smooth/Riemann adapter enters these proofs. |
| [R37.1–R37.7](https://github.com/Parveen117/extra-ideas/blob/f06932dd63a4f0827a31583bc42dfaa979ddc663/02-relational-response/NATIVE_PACKET_FLIGHT_R37.md) | **Native-derived:** compact source packets, all-count error, approaching speeds and finite arrival windows. |
| [R38.4–R38.7](https://github.com/Parveen117/extra-ideas/blob/79a4c4c2e15aad1391f7633ab0ab0b65ce580111/02-relational-response/NATIVE_UNIVERSAL_SIGNAL_CONE_R38.md) | **Native-derived:** universal localized fixed-threshold cone for free source legs; it is not assumed unchanged through arbitrary added controls. |
| [Canonical engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged native exact arithmetic/replay:** application evidence, not a second engine, whole-engine PASS or formal-assistant proof. |

```bash
python3.12 -B 04-operator-evolution/verify_r39.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```
