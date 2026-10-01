# R21: record identity from native current, and exact continuation memory

Research owner: Monty Dabas. Development: 1 October 2026.

R20 constructed fresh and reused records but did not identify an unknown
allocation from its response. R21 derives that identification and the memory
needed to continue it. Near the source front, current counts the previous
uses of the slot being written. A native H contrast at one earlier event
then decides whether that event used the same slot. All slot identities are
reconstructed up to renaming, with exact integer-valued tests.

No classical stochastic process, environment, measurement law, Fourier
theorem, physical clock or coupling strength is a premise. We expand the
existing H/K histories and their derived matching pairing. Allocation words
below enumerate available compositions of R20's gate; they are not a newly
postulated physical scheduling law. Identification from a response and
prediction of which protocol nature selects remain different achievements.

## Source and constructed protocol class

Consume unchanged R16 C1, C3–C8 and C13; R18's cut-address transport; R19's
pair record; and R20's native tuple gate and matching-record readout. Write
e_0=e_k, e_1=e_c, H e_a=(-1)^a e_a, K e_a=e_(1-a),

\[
F=H+K,\quad\Pi_a=(I+(-1)^aH)/2,\quad L_a=\Pi_a F,
\quad s_a=(-1)^a,
\]

\[
(U\psi)(x)=2^{-1/2}\sum_{a=0}^1L_a\psi(x-s_a),\qquad
W_i=\Pi_0\mathbin{\mathrm{tensor}}I+
     \Pi_1\mathbin{\mathrm{tensor}}K_i.
\tag{21.0}
\]

K_i exchanges the two marks in record slot i and leaves other marks alone.
This action is constructed on finite tuples, so K_i²=I and K_i K_j=K_j K_i
follow on every tuple. We do not assume a global grading on the entire UGD
carrier from these facts. The native tuple pairing is matching coefficients,
as proved in R20. The notation tensor is its derived product abbreviation.

An allocation word sigma=(sigma_1,...,sigma_N) composes U and then W_(sigma_n)
at event n. Each gate has a specified record target; distinct targets are
initially distinct typed slots. Renaming all slots consistently changes no
response. Unless otherwise stated the preparation is the single system seed
psi_0(0)=e_0 and all record marks are e_0. An event is a word position, not
a physical duration. The total number of slots used is r. General joint
preparations are explicitly allowed in the moment and observer statements.

Binary parity below means an even or odd number of applications of the
already derived exchange K. XOR abbreviates composition of those toggles;
no independently supplied binary field is required.

## Written results

### R21.1 — Slot allocation has a canonical history action and an exact quotient

For a direction history d=(d_1,...,d_n), define

\[
r_{\sigma,i}(d)=\Big(\sum_{k:\sigma_k=i}d_k\Big)\bmod2,
\quad x(d)=n-2|d|,\quad |d|=\sum_k d_k.
\tag{21.1}
\]

Its record is exactly e_(r_sigma(d)). The source word itself remains a
separate upstream tag. If b distinct slots occur in this prefix, the map
d -> r_sigma(d) has 2^b attained records, each with exactly 2^(n-b) histories.
The free linear record map therefore has rank 2^b and kernel dimension
2^n-2^b. The exponent n-b counts independent parity constraints removed
from binary words; it is not that linear kernel dimension.

Two record words agree precisely when every slot receives an even number
of differing negative directions. Their equality persists under a common
admitted continuation of writes. Conversely, unequal records are already
different for the record-equality target, so this is its exact continuation
congruence. It does not reconstruct full history identity.

For these single-origin histories there is also a joint invariant

\[
\boxed{x(d)+2|r_\sigma(d)|\equiv n\pmod4.}
\tag{21.2}
\]

**Proof.** A negative direction applies exactly one K to its designated slot;
a positive direction applies I. Cancel pairs of K on each slot using K²=I.
Exchanges on distinct tuple positions commute, giving (21.1). In every used
slot choose one occurrence. Assign the other n-b direction marks freely;
the chosen occurrence is then uniquely determined by the desired slot
parity. This proves the image and fibre counts without a rank theorem over
an imported field. The free ledger map has one nonzero basis image per
attained record. Differences of every other history in a fibre with a chosen
representative give an independent kernel basis of size 2^n-2^b.
Appending a common write toggles both records equally, proving continuation
compatibility. Finally |r_sigma(d)| and |d| have the same parity, and
x=n-2|d| gives (21.2). Allocation words with the same equality partition of
event positions are related by a unique bijection between their used slot
names. Numbering slots by first appearance gives their canonical form.

### R21.2 — The allocation partition determines the exact pair-retention kernel

Let epsilon(d) be the product of native F transition signs: only a transition
1 to 1 contributes -1. At n>=1 the joint response is

\[
A_n=2^{-n/2}\sum_{d\in\{0,1\}^n}
 \epsilon(d)e_{x(d),d_n}\mathbin{\mathrm{tensor}}e_{r_\sigma(d)}.
\tag{21.3}
\]

The unresolved pair field is the sum of ordered path-pair coefficients with
the exact retention factor

\[
\boxed{\kappa_\sigma(d,e)=
\mathbf1_{\forall i:\ \sum_{k:\sigma_k=i}(d_k\mathbin{\mathrm{xor}}e_k)
                          \text{ is even}}.}
\tag{21.4}
\]

This includes every allocation partition, not just the two R20 endpoints.

**Proof.** Every step appends a direction with native amplitude +/-1/sqrt(2)
and toggles its slot as in R21.1. Induction gives (21.3), with R18's bijection
retaining the original H/K word tags. Matching record marks in R20's readout
gives delta_(r(d),r(e)), which is (21.4). Signed addition occurs only after
these labelled contributions have been formed. In particular a surviving
path pair can cancel with another signed contribution; survival alone does
not prove a nonzero total measured response. The source normalization and
reversible gate preserve total pairing energy throughout the joint evolution.

### R21.3 — Exchange moments give a closed continuation and a lawful forgetting rule

For any finite joint response A(x,m), with a two-role system column at each
address x and record tuple m, define the exchange pair fields

\[
M_\eta(x,y)=\sum_m A(x,m)A(y,m\mathbin{\mathrm{xor}}\eta)^\dagger,
\qquad\eta\in\{0,1\}^r.
\tag{21.5}
\]

M_0 is R20's complete system pair readout. Writing slot i after U gives

\[
\boxed{M'_\eta(x,y)=\frac12\sum_{a,b=0}^1
 L_a M_{\eta\mathbin{\mathrm{xor}}((a\mathbin{\mathrm{xor}}b)e_i)}
      (x-s_a,y-s_b)L_b^\dagger.}
\tag{21.6}
\]

For a specified remaining continuation, only masks supported on slots that
will be addressed again can affect future M_0. With the stipulated initial
blank-product record preparation, slots never used yet are still blank,
so every moment containing such a slot vanishes until its first
write. Consequently one need store only the 2^ell moment fields supported
on **live slots**: slots used already and used again in that continuation.
An old slot can leave this ledger after its last admitted use. This contracts
its information into retained fields; it does not claim physical destruction
of its record. Future admissibility must be known as part of the stated
continuation contract, not inferred from present silence. This live-slot
reduction assumes those initially independent blanks; (21.6) itself admits
arbitrary initial joint correlations.

For the single-origin blank preparation, a useful support grading is

\[
M_\eta(x,y)\ne0\ \Longrightarrow\
x-y\text{ is even and }(x-y)/2\equiv|\eta|\pmod2.
\tag{21.7}
\]

**Proof.** The new joint column at (x,m) is
2^(-1/2) sum_a L_a A(x-s_a,m xor a e_i). Substitute twice into (21.5),
then change the summation mark to m xor a e_i. This gives (21.6), with no
omitted initial correlation. The equation uses only the same mask or that
mask toggled at the next written slot. Starting from mask zero and working
backward through the remaining word never introduces a slot outside that
word. A never-used e_0 factor has B(e_0,K e_0)=0, giving the zero moments.
At each cut, keeping all subsets of the live set therefore supplies every
nonzero predecessor needed by (21.6); induction proves exact future readout.
The lower bound for an explicitly richer exchange-probe catalogue is R21.9;
no minimality over all observers of a single fixed seed is claimed here.
For (21.7), a contributing pair has r(e)=r(d) xor eta. Apply (21.2) to both
histories and subtract. A finite sum cannot create support outside that rule.

### R21.4 — Native cut projectors provide an independent exact evolution

For signs chi_i in {+1,-1}, construct on the record ledger

\[
E_\chi=2^{-r}\prod_{i=1}^r(I+\chi_iK_i),\qquad
Z_\chi=2^{-r}\sum_\eta\Big(\prod_i\chi_i^{\eta_i}\Big)M_\eta.
\tag{21.8}
\]

Then E_chi are self-dagger mutually disjoint projections summing to I,
K_i E_chi=chi_i E_chi, and M_0=sum_chi Z_chi. On writing slot i,

\[
Z'_\chi=D_{\chi_i}\,\mathcal T Z_\chi\,D_{\chi_i},\qquad
D_{+1}=I,\quad D_{-1}=H.
\tag{21.9}
\]

For blank records initially, Z_chi,0=G_0/2^r. Thus the entire unresolved
response is the finite sum of native signed transports selected by the
record sectors. This is a derived algebraic decomposition, not random sign
noise assigned to events.

**Proof.** K_i²=I gives (I+chi_i K_i)²=2(I+chi_i K_i). Opposite signs give
zero product on that slot. Distinct slots commute; multiplying these
identities proves projection and disjointness. Summing each sign gives 2I,
so the projections sum to I. The eigenvalue identity follows by multiplying
one factor by K_i. Expansion yields (21.8). Summing chi eliminates every
nonempty mask by paired sign cancellation, proving M_0=sum Z_chi.
In (21.6), changing mask by e_i multiplies its signed sum by chi_i.
The factor for a,b is therefore chi_i^(a+b), which is left and right
multiplication by Pi_0+chi_i Pi_1. This gives (21.9). Initially only M_0
survives on blank records, giving equal coefficients 2^-r. These coefficients
come from the finite projection product. No probability independence axiom
or imported transform theorem enters the proof.

### R21.5 — Front current counts record returns exactly

Let c_i(n)=#{k<=n:sigma_k=i}, and let
nu_n=#{k<n:sigma_k=sigma_n} be the previous-use count of the current slot.
Read current and energy at the derived address one inward step from the
positive front:

\[
\jmath_n=j_n(n-2),\qquad\mathcal E_n=\rho_n(n-2),\qquad n>=1.
\]

Their exact values are

\[
\boxed{\jmath_n=\frac{\nu_n}{2^{n-1}},\qquad
\mathcal E_n=\frac{1+\sum_i c_i(n-1)^2}{2^n}.}
\tag{21.10}
\]

In particular jmath_n=0 precisely for a fresh slot, and its positive value
counts returns. At successive depths the same quantities satisfy

\[
\boxed{2^{n+1}\mathcal E_{n+1}-2^n\mathcal E_n
       =1+2^n\jmath_n.}
\tag{21.11}
\]

**Proof.** A path reaches n-2 exactly when it has one negative direction.
Call its position k. No such path has a 1-to-1 transition, so every native
coefficient is +2^(-n/2). Its record is the single marked slot e_(sigma_k).
For k<n the final role is 0; for k=n it is 1. Thus the two record-valued
role columns at this address are sum_(k<n)e_(sigma_k) and e_(sigma_n),
times 2^(-n/2). Their matching product is nu_n/2^n and current doubles it.
The squared norms are sum_i c_i(n-1)^2 and 1. This proves (21.10), also at
n=1 with an empty first column. Increasing c_(sigma_n) by one changes the
sum of squares by 2nu_n+1; subtracting the two energy formulas proves
(21.11). There are no sign cancellations in this extremal path family.

### R21.6 — One native contrast reconstructs equality of two slot identities

In a repeat of the same prepared protocol, insert H on the system just after
U at event i<n. Keep all other writes unchanged, and denote the resulting
front current at n by jmath_n^[i]. Then

\[
\boxed{2^{n-2}(\jmath_n-\jmath_n^{[i]})
       =\mathbf1_{\sigma_i=\sigma_n}.}
\tag{21.12}
\]

Baseline readings and these contrasts determine the entire allocation
partition up to slot renaming. Two allocation words give the same complete
family of these readings exactly when they differ only by such a renaming;
then their entire constructed joint evolutions are related by the record
permutation, and all unresolved responses agree.

**Proof.** On the one-negative histories used in R21.5, H reverses precisely
the history whose negative event is i. Because i<n it changes the role-0
record column by -2 e_(sigma_i), and leaves the final role-1 column unchanged.
Hence jmath_n^[i]=(nu_n-2 delta_(sigma_i,sigma_n))/2^(n-1), proving
(21.12). Each pair of event positions is recovered by taking the later as n.
Equality of these decisions is exactly equality of their partition, which
R21.1 identifies up to renaming. Conversely a consistent slot permutation
intertwines each W and preserves the blank tuple and its pairing, proving
equality of every response. The inserted H is a derived pairing-preserving
native arrow and commutes with W at that event. The protocol uses repeated
preparation, specified event positions and the signed current readout. It
does not assert an instrument can copy an unknown state or measure an
unavailable channel. Data failing the exact 0/1 and equivalence-relation
tests cannot be silently repaired into a certified partition.

### R21.7 — Fresh and fully coherent responses are rigid endpoints

At any fixed depth n>=1 in this protocol class, the entire local pair of
targets (rho_n,j_n) equals the unrecorded single-seed response if and only
if all n events use one slot. It equals the source-count response with zero
current if and only if all n events use distinct slots.

The uncontrasted front readings alone do not identify every intermediate
allocation: 0101 and 0110 have the same front current and energy at each
prefix, while their slot equality partitions differ.

**Proof.** For the coherent unrecorded response, R21.5's one-negative paths
sum without any distinction, giving jmath_n=(n-1)/2^(n-1). Equality with
an allocated response forces nu_n=n-1, so every earlier slot equals the
current one. Sufficiency is R20.5. For the count response the front energy
is n/2^n and current is zero. The energy equality forces
sum_i c_i(n-1)^2=n-1=sum_i c_i(n-1). Since each c_i is a nonnegative integer,
every nonzero c_i must be one. Zero current says that event n uses none of
them. Thus all n slots are distinct; sufficiency is R20.3–R20.4. The n=1
case satisfies both endpoint descriptions. For the intermediate examples,
the previous-use counts are (0,0,1,1) and the prefix multisets of slot-use
counts agree. Equations (21.10) give the stated identical front readings.
But at i=1,n=3, (21.12) returns 1 for 0101 and 0 for 0110. Thus a reuse
count is not a recovered identity. No equality of all their other unprobed
responses is asserted by this example.

### R21.8 — Returning local information has an exact joint continuation gate

At a common cut take two tagged histories h,g with addresses x_h,x_g and
record difference delta=r(h) xor r(g). Let beta be a common admitted future
allocation of t writes, and let u,v be two future direction words. Their
records and addresses coincide after this continuation exactly when

\[
\boxed{\delta=r_\beta(u)\mathbin{\mathrm{xor}}r_\beta(v),\qquad
       x_h-x_g=2(|u|-|v|).}
\tag{21.13}
\]

To contribute to the local energy trace, their final roles must additionally
agree; to the local current, their final roles must differ. Each admitted
F transition has nonzero native coefficient, so these conditions are also
sufficient for a nonzero **individual path-pair** contribution to that target.
They do not prevent cancellation in the whole signed sum.

If a nonzero component of delta belongs to a slot absent from beta, no such
record return is possible in that continuation. Any return also needs

\[
t\ge\max\{|x_h-x_g|/2,|\delta|\}.
\tag{21.14}
\]

**Proof.** R21.1 makes the final record difference delta xor r_beta(u)
xor r_beta(v); setting it to zero gives the first equality. A direction
word advances its address by t-2 times its negative count. Equating the
two final addresses gives the second equality. The trace of e_a e_b^dagger
is delta_ab and the trace after K is delta_(a,1-b), giving the role tests.
All source transition coefficients are +/-1 before the common normalization,
so an individual allowed term is nonzero. An unaddressed slot cannot be
toggled, proving the exclusion. At most t more negative directions can
occur on one side than the other, giving the address bound; each event can
toggle at most one bit of the record difference, giving the other bound.
These are necessary bounds; they are not sufficient when the future slot
catalogue or final-role condition prevents a return. Equal present shadows
do not authorize an unprovided write or identify its record target.

### R21.9 — The complete exchange-probe target has a sharp native memory count

For ell retained slots suppose the admitted probe catalogue includes every
finite product of their gates W_i, acting on a freshly prepared system role
C e_0=(e_0+e_1)/sqrt(2), with no intervening U. Let D be a record pair array
in the native finite linear pair ledger. The measured current after the
product indexed by a mask eta is

\[
q_\eta(D)=\sum_m D_{m,m\mathbin{\mathrm{xor}}\eta}.
\tag{21.15}
\]

The complete family has exactly 2^ell independent native scalar linear
channels. It is a minimum for this specified exchange-probe target; on the
known normalized affine slice, 2^ell-1 variable channels plus the known
normalization suffice and are necessary. This is not a lower bound for a
single known seed trajectory or arbitrary nonlinear coding.

**Proof.** The probe's two system roles acquire record arrows I and K^eta.
Their cross pairing is q_eta/2, so the doubled current is (21.15). For each
sign tuple chi form the native record vector
r_chi=2^(-ell/2) product_i(e_0+chi_i e_1) on tuples, and its pair array D_chi.
It has unit norm and K_i r_chi=chi_i r_chi. Hence
q_eta(D_chi)=product_i chi_i^eta_i. The square table of these values has
row product sum_chi product_i chi_i^(eta_i+zeta_i)=2^ell delta_(eta,zeta):
if the masks differ, pairing the two signs at a differing slot cancels the
sum; if they agree, all terms are one. Therefore its rows are independent.
Any linear observer determining all q_eta must have rank at least 2^ell
on the span of the D_chi; the q_eta themselves attain that rank. On the
unit-normalized slice q_0=1. The remaining independent variations are the
differences D_chi-D_chi0, a space of dimension 2^ell-1, proving the affine
statement. This is also a concrete instance of the native observer-rank
principle in RKF N03, rederived here rather than imported as a classical
tomography theorem. Its separately addressed probe catalogue is explicit.

## Certification and source citations

[R21_DERIVATION_LEDGER.json](../04-operator-evolution/R21_DERIVATION_LEDGER.json)
binds each written proof and its native dependencies. The gate rejects
classical or admitted model inputs on proof paths, unpinned sources, changed
citations, cycles and physical promotion of protocol definitions. It audits
declared provenance and bound text, not arbitrary proof semantics.

The [native adapter](../04-operator-evolution/record_identity_continuation.cjs)
checks literal original H/K histories against joint gate propagation,
independent pair/moment and sector evolution, all allocation partitions in
the stated finite domains, current contrasts, live-slot reduction, tagged
return witnesses and observer rank. General n and ell statements have the
written proofs above. Finite enumeration is not their replacement.

The [verifier](../04-operator-evolution/verify_r21.py) reproduces R20 and the
frozen R19/R18/R17/R16 chain. All 173 pre-existing non-navigation files and
the canonical engine remain byte-identical. Source provenance, written
proof, finite verification, formal proof-assistant checking and physical
validation have separate statuses; no formalization or new experiment is
claimed.

| Pinned source | Exact use and premise status |
| --- | --- |
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native free-history and marked-role constructions, refinements, pairing and completion; their finite target scope is retained. Logical equality and finite induction remain the explicitly declared metamathematical infrastructure. |
| [R18](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Native direction histories, constructed address and derived normalized transport; no physical coordinate/clock selection. |
| [R19](https://github.com/Parveen117/extra-ideas/blob/f46316e4407a65bca1e07f43baf32a419048f67e/02-relational-response/CURRENT_MEMORY_RETENTION_R19.md) | Native ordered pair propagation and distinction between hidden current and local energy. |
| [R20](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | Native tuple gate, matching-record readout and fresh/reused endpoints; allocation and preparation remain explicit constructions. |
| [RKF N03](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Native minimum observer for a declared catalogue and target. R21.9 supplies its own elementary finite proof; no full upstream-engine PASS is inferred. |
| [Publications U13–U16](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/ADMISSIBLE_CONTINUATION.md) | Comparison only: finite carrier, cut and partial-arrow domains are supplied. Availability is not automatically an observation. Those model inputs do not select the present allocation. Its certificate is not rerun. |
| [Publications NI-1–NI-4](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-return-identification/THEOREM.md) | Comparison only: a supplied paired-cell family and symmetric probes identify a parameter; identification does not predict the parameter or alpha. Those couplings and probe assumptions are not premises here. |

Binary parity codes, finite sign-transform tables and controlled exchanges
have established mathematical counterparts. Their occurrence here does not
establish historical priority. Their formulas needed by this packet were
derived above from native tuple operations and finite cancellation.

## What closes, and what follows

Within the source-constructed W protocol class, fresh versus reused is now
an exact response distinction, and slot identity has a complete contrast
decoder. The joint continuation law says which old distinctions can return
and which records a specified future can lawfully forget. This is stronger
than naming a slot or fitting an external noise strength.

The source algebra admits all these distinguishable allocation partitions.
These results do not select one physical partition from no further data.
A measured current would identify its native reuse count under the declared
interface, not predict that measured current from H/K alone. Next construct
the return/meeting relation and its event-distance directly from the joint
address-record continuation gate, retaining the allowed-arrow contract.
That is the next geometric target before a physical metric, clock, c or alpha.
