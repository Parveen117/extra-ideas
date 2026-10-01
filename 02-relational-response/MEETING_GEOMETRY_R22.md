# R22 — Native meeting geometry and address-carried record parity

Research owner: **Monty Dabas**. Development: 1 October 2026.

R21 gave the exact joint address/record return gate. Here that gate yields a
constructive decision for every finite future allocation word, a sharp minimum
number of events under the explicitly constructed complete slot catalogue,
and an exact reduction of the record representation by one binary mark.
Current return and energy return have different last-step requirements.

The distance is derived from native continuation, not postulated as a Euclidean
or spacetime metric. The complete slot catalogue is a mathematical construction
from already derived arrows, not a proof that every such arrow is physically
available. A particular physical allocation, rod, clock, c and alpha are not
selected by this packet. No classical metric, stochastic law, transform theorem
or complex-number primitive is a premise. Finite words, equality, integer
counting and induction are the metamathematical infrastructure declared in R16.
All full H/K word tags remain upstream of the endpoint observations.

## Inherited source and constructed target

Use R16's derived role parity H, exchange K and matching pairing; R18's
F=H+K, native normalized step U and address x; R20's finite record tuples
m and reversible write W_i; and R21's allocation word and joint gate.
The role labels a=0,1 have displacements s_a=1-2a. Evaluating the source
role maps gives F's entries (1,1;1,-1), so every role-to-role transition
has a nonzero coefficient. At a common cut the endpoint of a tagged history
is (x,a,m). Neither this endpoint nor any distance below identifies full words.

For two endpoints define delta=m xor n. If x-y is even set d=(x-y)/2.
Here |delta| counts differing marks; it is not a supplied norm. A future word
beta=(beta_1,...,beta_t) specifies which already constructed W is used at
each event, after U. Both histories share beta, but may take different native
direction branches u,v. A **meeting** means equal final address and record.
Energy additionally needs equal final roles; current needs different roles.
All existential return statements concern individual tagged pair terms.

## Derived results

### R22.1 — Relative native motion has a preserved parity component

At a write of slot i, outgoing roles (u,v) give exactly these relative moves:

| Outgoing roles | Relative state after the event |
| --- | --- |
| (0,0) or (1,1) | (d, delta) |
| (1,0) | (d-1, delta xor e_i) |
| (0,1) | (d+1, delta xor e_i) |

The component d-|delta| modulo two is invariant. Every record bit outside
the future slot set is invariant. Odd address separation can never meet.
For the single-origin blank preparation at depth n, every retained endpoint
satisfies x+2|m| congruent to n modulo four.

**Proof.** Subtract the two source displacements and divide by two:
(s_u-s_v)/2=v-u. R20's two controlled writes change the record difference
by (u xor v)e_i. This gives the table directly. A nonidle move changes d
by one and changes record weight by one with either sign; an idle move
changes neither. No write touches a different slot. The address difference
changes by an even integer, proving the odd-separation exclusion. A direction
word with k negative branches has x=n-2k, while its record weight has parity
k by pairwise cancellation of repeated slot toggles. This proves the last
congruence, also inherited from R21.1.

### R22.2 — Every fixed future word has an exact constructive capacity test

Let c_i be the number of occurrences of slot i in beta and w=|delta|.
First require **coverage**: c_i>=delta_i for every record bit, including
bits absent from beta. Under that coverage define

\[
Q_\beta(\delta)=\sum_i\left[c_i-((c_i-\delta_i)\bmod2)\right].
\tag{22.1}
\]

For even x-y, a meeting at exactly t events is possible if and only if

\[
\boxed{\text{coverage},\qquad d\equiv w\pmod2,\qquad |d|\le Q_\beta(\delta).}
\tag{22.2}
\]

Write A_beta(d,delta) for this predicate. At the empty word it is true
exactly for d=0,delta=0. For a fixed horizon the answer depends on slot
counts, not their order; earliest prefix return can depend on order.

**Proof.** Let k_i count disagreements between u and v at occurrences of
slot i. Meeting requires 0<=k_i<=c_i and k_i congruent to delta_i modulo
two. Their sum k must be at least |d| and have d's parity, because its k
nonidle moves each contribute +1 or -1 to the negative-count difference
|u|-|v|=d. This proves necessity. Under coverage, the possible k_i run
from delta_i to q_i=c_i-((c_i-delta_i) mod 2) in steps of two. Their sums
are every integer w,w+2,...,Q: start at each minimum and add pairs to any
slot with remaining capacity, until the desired total is reached. Choose
k=max(w,|d|), which has the required parity and lies in that range.
Pick k_i actual occurrences at each slot. Set (u,v)=(1,0) at (k+d)/2 of
the selected occurrences and (0,1) at the other (k-d)/2. Set (u,v)=(0,0)
elsewhere. This explicitly gives the required address sum and every record
parity. All these transitions have nonzero native coefficients. The argument
also covers t=0. Once a prefix meeting is possible, appending common
directions preserves it; no additional matching-event assumption is needed.

### R22.3 — Last roles separate current return from energy return

For t>=1 write beta^- for the word without its last slot i. At exactly t,
an individual local energy term exists precisely when

\[
A_{\beta^-}(d,\delta).
\tag{22.3}
\]

An individual local current term exists precisely when

\[
A_{\beta^-}(d-1,\delta\mathbin{\mathrm{xor}}e_i)
\quad\text{or}\quad
A_{\beta^-}(d+1,\delta\mathbin{\mathrm{xor}}e_i).
\tag{22.4}
\]

If initially (d,delta) differs from (0,0), and the earliest meeting prefix
is k, the first current term occurs at k and the first energy term at k+1,
provided that next event is in the stated continuation. If initially meeting,
the first positive energy return is event one; the first positive current
return is precisely the first repeated slot in beta, or never if none repeats.
At zero events, the initial roles themselves decide energy versus current.

**Proof.** The energy trace of e_u e_v^dagger is one exactly for u=v;
the current trace after K is one exactly for u!=v. Apply the three moves
of R22.1 at the last event. This proves (22.3)-(22.4). A first meeting
of an initially nonmeeting pair cannot have an idle last relative move,
since it would already have met at the preceding cut. Hence it is current;
one further equal-direction step makes it energy. No earlier energy is
possible by (22.3). From an initial meeting, a current return must have an
even positive number of disagreements at each slot. Its last differing
slot must therefore have occurred earlier. Conversely take the first two
occurrences of any repeated slot, use one +1 and one -1 relative move,
and use idle moves elsewhere. They return to zero with different last roles
at the second occurrence. A common first step supplies energy immediately.

### R22.4 — The complete native slot catalogue has a sharp event-distance

Construct the continuation catalogue allowing any one of the slots in a
fixed nonempty finite set J at every future event. No physical availability
is inferred from making this construction. Minimize the meeting length over
its words and the native branch pairs. For endpoint labels p=(x,m),q=(y,n),
the resulting distance is

\[
\boxed{D_J(p,q)=\max\{|x-y|/2,|m\mathbin{\mathrm{xor}}n|\}}
\tag{22.5}
\]

when x-y is even, delta is supported in J, and (x-y)/2 has delta's weight
parity. Otherwise D_J is infinite. If J is empty, only the empty continuation
exists: equal endpoints have distance zero and all other pairs infinite.

For compatible distinct endpoints, first current and energy term lengths
are D_J and D_J+1. At equal endpoints their first positive lengths are two
and one, respectively, for nonempty J.

**Proof.** Each event changes the relative address by at most one and at
most one record bit. Thus the maximum in (22.5) is a lower bound. Parity
and frozen-bit obstructions are R22.1. Put k=max(|d|,|delta|). For compatible
endpoints k-|delta| is even. Write each differing slot once, then insert
(k-|delta|)/2 pairs of writes to any one slot of J. This gives a word of
length k with the required record parity. Choose the nonidle branch signs
so that (k+d)/2 have negative-count difference +1 and (k-d)/2 have -1.
Their sum is d, so R21's gate gives a meeting at k. This attains the lower
bound. The empty-catalogue case follows from its definition, and the target
lengths follow from R22.3; two writes of one available slot attain the
initial-meeting current length. Positive target-return lengths have nonzero
diagonals and are not themselves metrics.

### R22.5 — The meeting distance is a metric on endpoints, not on full histories

For nonempty J the finite-distance components are exactly fixed values of
x+2|m| modulo four and fixed record bits outside J. D_J is a metric within
each such component and an extended metric across them. On endpoints that
also retain role a, or on full history tags, it is only a pseudometric within
a component: distinct roles or words can have distance zero.

**Proof.** Equality of the stated invariants is equivalent to the parity,
even-address and frozen-bit conditions in R22.4: |m xor n| has the same
parity as |m|-|n|. Formula (22.5) is symmetric and nonnegative, and vanishes
exactly when x=y and m=n. For three endpoints p,q,r in a component, every
record bit differing between p and r differs in at least one of (p,q) or
(q,r). Counting gives w_pr<=w_pq+w_qr. For address differences, concatenating
signed unit increments gives |x_p-x_r|<=|x_p-x_q|+|x_q-x_r|. Each summand
is bounded by its corresponding maximum in (22.5), yielding
D_J(p,r)<=D_J(p,q)+D_J(q,r). If an endpoint pair is in different components,
no third endpoint can have finite distance to both, proving the extended
inequality. Pullback to words or role-labelled endpoints keeps these
properties except separation, since the observed endpoint can coincide.
These are relative moves of two synchronous continuations. A single source
trajectory still advances address by one per event; the factor two in a
comparison distance is not a new propagation speed or a value of c.

### R22.6 — Record width determines exact meeting-ball volume

With all r>=1 record slots available, take a component and any endpoint as
origin. For integer L>=0 the number of endpoints at distance at most L is

\[
\boxed{V_r(L)=\sum_{w=0}^{\min(r,L)}{r\choose w}
             \left[L+\mathbf1_{L-w\ \mathrm{even}}\right].}
\tag{22.6}
\]

For L>=r this reduces exactly to

\[
V_r(L)=2^rL+2^{r-1}.
\tag{22.7}
\]

At a fixed address the component has 2^(r-1) record labels, with maximum
mutual distance 2 floor(r/2). For example records 00 and 11 at one address
are two events apart although their address separation vanishes.

**Proof.** Translate relative labels so that the origin has address zero and
record zero. Write x=2z. Component membership says z congruent to |m| modulo
two, and the ball requires |z|<=L and |m|<=L. Choose the w marked slots
in binomial(r,w) ways; this coefficient is defined by counting subsets.
There are L+1 integers in [-L,L] of parity L and L of the other parity,
proving (22.6). Once L>=r every record tuple is included. Flipping one fixed
bit pairs even-weight and odd-weight tuples, so exactly 2^(r-1) match L's
parity. This proves (22.7). At a fixed address precisely one record parity
is allowed, giving the same count. Record differences there have even
weight, whose largest possible value is 2 floor(r/2); all such masks occur.
Thus fixed finite record width gives eventual linear ball growth. This
does not derive a physical dimension. With an infinite complete slot
catalogue even a radius-one ball is infinite, so (22.6) has finite-r scope.

### R22.7 — Address and event already carry one record bit

At fixed event n restrict to the native parity subspace
x+2|m| congruent to n modulo four. Choose one of r>=1 slots to omit, called
slot r. Keep address, system role and the other r-1 bits b. Recover its bit by

\[
\boxed{m_r=\left((n-x)/2\bmod2\right)
               \mathbin{\mathrm{xor}}\bigoplus_{i<r}b_i.}
\tag{22.8}
\]

This defines an exact pairing-preserving bijection of basis labels between
that subspace and (x,a,b) with x congruent to n modulo two. Actual reachable
support for a fixed allocation can be smaller. In this representation,
writing a retained slot is U followed by its retained W_i; writing the
omitted slot is just U. Reconstruction at each new event restores the
complete record field. It removes a redundant coordinate, not a physical
record or a source-history tag.

Local energy and current can be read by matching b alone. If bar G is the
full system pair array formed by matching b, the original full-record readout
at fixed n is instead

\[
\boxed{G_0(x,y)=\mathbf1_{x-y\equiv0\pmod4}\,\bar G(x,y).}
\tag{22.9}
\]

Among encodings retaining address and each admissible record identity, r-1
binary marks is sharp on this entire parity subspace. It is not a minimum
for a single known trajectory or for arbitrary encodings of full histories.

**Proof.** R22.1 equates record parity to (n-x)/2 modulo two, so the missing
bit is uniquely (22.8). Conversely that formula constructs one admissible
full record for every retained label. This is a permutation of an
orthonormal finite-support basis, proving preservation of the matching
pairing. After an outgoing role a, x'=x+1-2a and n'=n+1, whence
(n'-x')/2=(n-x)/2+a. The required full-record parity toggles exactly when
a=1. A retained-slot write toggles its b_i by a; the recovered omitted bit
then stays fixed. An omitted-slot write keeps b unchanged, and the recovered
omitted bit toggles by a. These are exactly the original W actions. Native
F coefficients are unaffected, proving stepwise conjugacy and all finite
continuations. At a common address, equal retained records have equal
recovered bits, so local readout is unchanged. At x,y with equal b their
recovered bits agree exactly when (x-y)/2 is even, proving (22.9). Finally
R22.6 gives 2^(r-1) distinct records at a fixed admissible address; fewer
binary marks cannot label them injectively. When r=1 no explicit record
coordinate is needed for local response, recovering R20.5, while the
off-address factor in (22.9) remains essential.

### R22.8 — First slot reuse fixes the first current and energy departures

For the original e_0 seed at address zero and initially blank records, let
tau be the first event whose slot occurred earlier. If it never occurs,
the response stays the fresh-record source-count response. Otherwise:

* Current is zero everywhere for n<tau, and
  j_tau(tau-2)=2^(1-tau) is nonzero.
* Energy equals source counts everywhere through n=tau. If event tau+1
  is in the protocol, its first departure has the exact witness
  rho_(tau+1)(tau-1)-rho_count,(tau+1)(tau-1)=2^(-tau).

**Proof.** All prefixes before tau have fresh slots, so R20.3-R20.4 give
the count response and zero current. At tau the written slot has exactly
one previous occurrence; R21.5 gives the stated current. The write does
not change local role energies. The derived transport balance is
rho_(n+1)(x)=[rho_n(x-1)+rho_n(x+1)+j_n(x-1)-j_n(x+1)]/2.
Since j_(tau-1)=0, it gives the count energy also at tau. Among the first
tau writes, one slot has count two and tau-2 slots count one. R21.5 thus
gives rho_(tau+1)(tau-1)=(1+4+tau-2)/2^(tau+1)=(tau+3)/2^(tau+1).
The source-count value there is (tau+1)/2^(tau+1), because there are
tau+1 one-negative direction words. Subtraction gives 2^(-tau), regardless
of the next slot. This is a noncancelling extremal family for the specified
seed; R22.3 alone would not establish a nonzero total response.

### R22.9 — The signed pair propagator limits what geometry certifies

On an ordered pair basis |x,a,m><y,b,n|, a write to slot i after U gives

\[
\frac12\sum_{u,v=0}^1 F_{ua}F_{vb}
 |x+s_u,u,m\mathbin{\mathrm{xor}}u e_i\rangle
 \langle y+s_v,v,n\mathbin{\mathrm{xor}}v e_i|.
\tag{22.10}
\]

Geometry determines exactly which individual tagged terms can reach each
local target. It does not ensure their signed sum is nonzero. In particular
for the initial pair |0,0,0><0,1,0| and future word (i,i), the two current
terms at address zero have coefficients -1/4 and +1/4 and cancel.

**Proof.** Apply the two derived normalized branch maps to the two basis
endpoints. Native matching dagger fixes their real coefficients, giving
(22.10). Linearity extends it to any finite pair ledger, retaining all
path-pair tags until final addition. In the displayed example the path pairs
(u,v)=(01,10) and (10,01) both end at address zero with record e_i and
opposite final roles. The first path on the right of the first pair begins
with 1 to 1 and contributes the single minus sign. The second pair has no
1 to 1 transition. Both normalized pair magnitudes are 1/4, so the local
current sum is zero. All other paths fail the final current test. This pair
basis element is a linear coefficient channel, not by itself a positive
preparation. A positive preparation can also have a vanishing return:
(e_0+e_1)/sqrt(2) at the origin with blank records, followed by two writes
of i, has first-step role support only 0 because F(e_0+e_1)=2e_0.
After the second step its two roles occupy different addresses, so its
local current vanishes everywhere despite the catalogue admitting current
return terms. Thus target support, signed signal and preparation selection
remain distinct. No physical metric or constant follows just from (22.5).

## Certification and pinned source citations

The [derivation ledger](../04-operator-evolution/R22_DERIVATION_LEDGER.json)
binds these nine written sections to their native dependencies and scoped
exact checks. Its gate checks declared provenance and hashes, not arbitrary
proof semantics. The [adapter](../04-operator-evolution/meeting_geometry.cjs)
uses the unchanged canonical native arithmetic and proof replayer. It compares
the capacity construction with exhaustive direction pairs, the distance with
independent finite graph exploration, and compressed evolution with complete
record evolution. The [verifier](../04-operator-evolution/verify_r22.py)
replays the frozen R21/R20/R19/R18/R17/R16 chain in temporary outputs.
All 180 prior non-navigation files remain byte-identical. Finite checks do
not replace the arbitrary-depth proofs above. No Lean/Coq formalization or
new experimental validation is claimed.

| Pinned source | Exact use and premise status |
| --- | --- |
| [R16 C1, C3-C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native role maps, pairing and full tagged histories; declared logical equality, finite words and induction remain metamathematical infrastructure. |
| [R18](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Derived direction bijection and normalized transport on its explicitly constructed address target. No physical coordinate or clock selection. |
| [R19](https://github.com/Parveen117/extra-ideas/blob/f46316e4407a65bca1e07f43baf32a419048f67e/02-relational-response/CURRENT_MEMORY_RETENTION_R19.md) | Native signed pair propagation and local energy/current balance. |
| [R20](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | Derived reversible tuple write and matching-record readout; blank preparation and allocation are explicit constructions, not physical laws. |
| [R21](https://github.com/Parveen117/extra-ideas/blob/2888aa2a9ea1c63a89b3d695256c8b9405dbd7cc/02-relational-response/RECORD_IDENTITY_CONTINUATION_R21.md) | Native record invariant, exact tagged meeting gate and noncancelling front identities; no imported metric premise. |
| [RKF native algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Native identities replayed in their declared quotient; scoped replay is not a full upstream-engine PASS. |
| [Publications U13-U16](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/ADMISSIBLE_CONTINUATION.md) | Comparison only: carrier, cut, partial-arrow domains and availability contract are supplied inputs. They are not imported as a selection of J here; its certificate is not rerun. |
| [Publications NI-1-NI-4](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-return-identification/THEOREM.md) | Comparison only: supplied paired-cell family, couplings and probes identify a parameter. They do not predict alpha and are not proof premises here; its certificate is not rerun. |

Binary parity counting and shortest-path metrics have established mathematical
counterparts; historical priority is not claimed. Every needed formula and
metric inequality here is derived from the native move grammar above.

## Next development

The joint return gate now has an exact decision procedure, a sharp catalogue
distance and a memory-coordinate reduction. The next source question is which
different continuation catalogues have identical future energy/current
responses after this reduction, and what native response separates the
remaining ones. That will test whether a proposed metric choice is observable
or still an unselected target definition. A physically common clock/rod and
the route to c or alpha must follow such a justified selection, not a fitted
identification of this comparison distance with spacetime.
