# R43: native exchange, relative calibration and interaction memory

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R42 derived a complete phase generator for a finite native process.
R43 couples two copies through an exchange made from the existing native
record controls. The exchange conserves the sum of their phase readouts
and fixes their relative multiplicative calibration when that readout
is nonconstant. It also links interaction phase cost to native
distinguishability and to the memory visible after ignoring one ledger.

For normalized product inputs a and b, put kappa=|B(a,b)|^2. The central
identities, on the stated exchange phase branch, are

\[
\boxed{\frac{\mathcal I(a,b)}{\vartheta_-}
       =\frac{1-\kappa}{2},\qquad
 \mathcal M(\tau)=2\sin_\Sigma^2(\vartheta_-\tau)
                  \left(\frac{\mathcal I(a,b)}{\vartheta_-}\right)^2.}
\tag{43.A}
\]

Here I is a generator expectation and M is the explicitly defined
reduced-pair deficit `1-Tr(rho_A^2)`. Neither is introduced as physical
energy, thermodynamic entropy or random noise. Native tuple matching,
the existing record controls and R42's phase interpolation derive both.

The common physical unit remains open. Relative calibration is fixed
by a declared exchange-conservation target; this does not uniquely
select a material interaction or the positive half-turn branch.

## Native inputs and the interaction target

Use R16's earned H/K roles, iota, scalar completion and positive matching,
R20's finite tuple records and controlled K, and R42's finite process
generator. Let M be a finite d-role ledger, d>=2. The two copies share
the same specified role identification and the same chosen phase
generator E, with `V=Exp_Sigma(-iota E)`.

Write

\[
 E_A=\mathcal E\otimes I,\quad E_B=I\otimes\mathcal E,
 \quad E_0=E_A+E_B.
\]

The tuple symbol denotes R20's constructed matching ledger, not a
postulate about physical subsystems. The phase parameter tau is the
R42 interpolation of the retained native tick. The positive native
half-turn `theta_-=4 integral_0^1 (1+s^2)^(-1) ds` was constructed in
R42, with `Exp_Sigma(-iota theta_-)=-1` and half-angle value -iota.
Use theta_- to mean the displayed vartheta_- in plain-text formulas.

Every expectation below is normalized by the state norm square.
The coefficient norm `||A||_F^2=sum_ij |A_ij|^2` is just a finite
matching sum over matrix-entry marks. It is used for exact defect
certificates, without adding a physical norm or trace axiom.

## Written results

### R43.1 — Native record controls construct full exchange and its half-tick arrow

On tuple marks define S e_(a,b)=e_(b,a). It preserves matching and
satisfies S-square=I=S-dagger S and S-dagger=S. For two source roles,
the already available R20 controls give exactly

\[
 W_{AB}=\Pi_0\otimes I+\Pi_1\otimes K,\qquad
 W_{BA}=I\otimes\Pi_0+K\otimes\Pi_1,
 \qquad\boxed{S=W_{AB}W_{BA}W_{AB}.}
\tag{43.1}
\]

Set `P_-=(I-S)/2` and `P_+=(I+S)/2`. R42 applied to this involution
has zero midpoint generator and reversal cut P_-. Thus its selected
phase generator and full interpolation are

\[
\boxed{E_{\rm int}=\vartheta_-P_-,\qquad
 W(\tau)=P_++z(\tau)P_-,\quad
 z(\tau)=\operatorname{Exp}_{\Sigma}(-\iota\vartheta_-\tau).}
\tag{43.2}
\]

In particular W(1)=S, W(2)=I, and the half tick is the exact native
cut-rational arrow

\[
\boxed{W(1/2)=\frac{1-\iota}{2}I+
                    \frac{1+\iota}{2}S,\qquad W(1/2)^2=S.}
\tag{43.3}
\]

**Proof.** Swapping tuple marks twice returns every mark and preserves
their equality, proving the matching and inverse identities. On the
two source labels, the three controls successively change
`(a,b)` to `(a,b+a)`, then `(b,b+a)`, then `(b,a)`, with addition
modulo two describing repeated K. This proves (43.1) on every mark.

Expanding S-square=I proves that P_+,P_- are complementary self-dagger
cuts. Their ranges are spanned by the equal marks and the symmetric
or antisymmetric combinations of each unequal pair. Direct counting
gives ranks d(d+1)/2 and d(d-1)/2, respectively. These ranks do not
assert physical particle statistics.

The cut factorial identity `Exp_Sigma(tP)=I-P+Exp_Sigma(t)P`, proved
in R42.7 by expanding powers, gives (43.2). The inherited half-turn
and half-angle identities give its three special times and (43.3).
Thus neither an interaction Hamiltonian nor a continuous coupling
constant was supplied to obtain this particular native arrow and
branch. Other interpolation branches allowed by R42 remain distinct
constructions.

### R43.2 — Additive process history derives a conserved exchange generator and retains branch memory

For two identical native process copies,

\[
\boxed{[S,E_0]=0,\qquad
 U(\tau)=\operatorname{Exp}_{\Sigma}(-\iota\tau E_0)W(\tau)
       =\operatorname{Exp}_{\Sigma}(-\iota\tau E_{\rm lift}),
 \quad E_{\rm lift}=E_0+E_{\rm int}.}
\tag{43.4}
\]

At one tick, U=(V tensor V)S. E_0 and E_int commute with the complete
evolution, so their individual expectations and all their joint
polynomial moments are conserved. This includes an exactly balanced
redistribution of E_A and E_B.

The additive branch retains more information than a fresh endpoint-only
choice of the R42 generator. For the source process V=H, put
`Q=(I-H)/2`, `Q_A=Q tensor I`, `Q_B=I tensor Q`, and
`Q_11=Q tensor Q`. Then E=theta_- Q and

\[
\boxed{E_{\rm lift}=\vartheta_-(Q_A+Q_B+P_-),\qquad
 E_{\rm lift}-E_{\rm end}=2\vartheta_-(Q_{11}+P_-),
 \quad E_{\rm end}=\vartheta_-(I-U)/2.}
\tag{43.5}
\]

Both exponentiate to the same U. Their fractional histories can differ.

**Proof.** On every product mark, S(A tensor B)=(B tensor A)S.
Applying this to the two copies of E proves [S,E_0]=0 and hence
[E_int,E_0]=0. The factors E tensor I and I tensor E also commute.
Grouping finite coefficients of the absolutely bounded native
factorial series gives the commuting product law. This proves (43.4)
and `Exp_Sigma(-iota E_0)=V tensor V`. Matching preservation and
commutation prove all stated moment invariants.

For (43.5), the four native subspaces are the 00 mark, the symmetric
01/10 line, the antisymmetric 01/10 line and the 11 mark. On these,
U=(H tensor H)S has values +1,-1,+1,+1, while E_lift/theta_- has
values 0,1,2,2. The involution's R42 generator has values 0,1,0,0.
Their difference is therefore the claimed multiple of the orthogonal
cut Q_11+P_-. The native full-turn identity makes it invisible at
integer ticks. This is retained composition history, not permission
to silently reset the branch after combining processes.

### R43.3 — Exchange balance determines relative calibration and quantifies its defect

For any self-dagger native readout A, pure exchange W(tau) conserves
`A tensor I+I tensor A`. Its separate responses satisfy

\[
\boxed{\frac{d\langle A\otimes I\rangle}{d\tau}
 =-\frac{d\langle I\otimes A\rangle}{d\tau}
 =\iota\langle[E_{\rm int},A\otimes I]\rangle.}
\tag{43.6}
\]

For A=E, this balance also holds under the complete U(tau), because
the free generator commutes with each local E. Arbitrary A need not
have a conserved sum under the additional free process.

Now assign positive native readout factors b_A,b_B to the same
non-scalar E on the two ledgers. Conservation of
`C=b_A E_A+b_B E_B` under the full exchange, for every preparation,
holds exactly when b_A=b_B. More quantitatively,

\[
\boxed{\|[S,C]\|_F^2
 =2d(b_A-b_B)^2\left[\operatorname{Tr}(\mathcal E^2)
              -\frac{(\operatorname{Tr}\mathcal E)^2}{d}\right].}
\tag{43.7}
\]

Offsets proportional to identity remain free. A scalar E cannot fix
relative scale by this target.

**Proof.** The swap identity used in R43.2 gives commutation with
every identical-readout sum. R42's native pairing derivative proves
(43.6). Since each local E commutes with E_0, its free contribution
vanishes, proving the complete-process statement.

The same tuple identity gives
`[S,C]=(b_A-b_B)(E_B-E_A)S`. Summing native coefficient squares and
using S's matching permutation gives
`||E_B-E_A||_F^2=2d Tr(E^2)-2(Tr E)^2`, proving (43.7).
The bracket is the sum of all off-diagonal squared scalar norms plus
the squared differences of real diagonal entries from their common
mean. It is nonnegative and vanishes precisely for a scalar matrix.

Conservation for every preparation is equivalent to S-dagger C S=C:
if a self-dagger difference had zero quadratic readout everywhere,
testing individual marks, their sums and their native iota contrasts
would set every entry to zero. Thus non-scalar E forces b_A=b_B,
and equal factors plainly suffice. This fixes a relative multiplier
within the specified exchange-conservation target, not an absolute
physical energy unit or arbitrary calibration functions.

### R43.4 — A native two-role exchange has an exact transfer current and requires pair memory

Use the source-H process of (43.5). Define the one-mark-on-each-side
sector and its native response arrows

\[
 P_o=(I-H\otimes H)/2,\quad
 Z=(H\otimes I-I\otimes H)/2,
 \quad K_o=(K\otimes K)P_o,\quad J=\iota K_oZ.
\]

On P_o, Z,K_o are the source H/K action; off that sector they vanish.
They obey `Z^2=K_o^2=J^2=P_o`, and
`P_-=(P_o-K_o)/2`. Write p=<P_o>, z=<Z>, j=<J>.
Under either the pure exchange or the complete source-H process,

\[
\boxed{\dot z=-\vartheta_-j,\qquad
       \dot j=\vartheta_-z,\qquad \dot p=0,\qquad
 \frac{d\langle Q_A\rangle}{d\tau}
       =\frac{\vartheta_-}{2}j
       =-\frac{d\langle Q_B\rangle}{d\tau}.}
\tag{43.8}
\]

For initial e_01, the exact transfer is

\[
\boxed{\langle Q_A\rangle=\sin_\Sigma^2(\vartheta_-\tau/2),
 \quad\langle Q_B\rangle=\cos_\Sigma^2(\vartheta_-\tau/2).}
\tag{43.9}
\]

Half a tick gives equal local counts; one tick exchanges them completely.
Local pair readouts alone do not predict the current.

**Proof.** Multiply the inherited H/K arrows on the four tuple marks.
This proves the sector identities and the expression for P_-.
The commutators are
`iota[P_-,Z]=-J`, `iota[P_-,J]=Z`, and
`iota[P_-,Q_A]=J/2`. E_0 is constant on P_o and commutes with these
readouts, so R42's response derivative gives (43.8) for both flows.
Equivalently, with phi=theta_- tau/2 and the native circular functions,
the exact exchange is a common phase times `cos(phi)I+iota sin(phi)S`.
Its action on e_01 proves (43.9) and the displayed half/full tick values.

For the predictive-memory witness take
`psi_+=(e_01+iota e_10)/sqrt(2)` and
`psi_-=(e_01-iota e_10)/sqrt(2)`. Both local pair readouts equal I/2,
and both have z=0 and p=1. Their j values are +1 and -1. Thus their
next Q_A slopes are opposite under the same fixed exchange. Conversely
e_01 and e_10 have the same j=0 and opposite z, hence opposite next
j slopes. The two real readings z,j close this oscillator; either
one alone is insufficient for all these preparations. This statement
does not claim that two readings reconstruct the whole joint state.

### R43.5 — Interaction phase cost measures native distinguishability and its local response metric

For any normalized joint state define the specified nonnegative phase
cost `I=<E_int>=theta_-<P_->`. For normalized product inputs a tensor b,

\[
\boxed{\mathcal I(a,b)=\frac{\vartheta_-}{2}
              \left(1-|B(a,b)|^2\right),\qquad
       0\le\mathcal I(a,b)\le\vartheta_-/2.}
\tag{43.10}
\]

Parallel inputs have zero cost and orthogonal inputs attain the product
ceiling. General joint states have the larger ceiling theta_-; the
native antisymmetric state `(e_01-e_10)/sqrt(2)` attains it. Thus a cost
above theta_-/2 proves that a normalized pure input is not a product.

For a short native process displacement `b(h)=Exp_Sigma(-iota h E)a`,
R42's response metric gives the interaction expansion

\[
\boxed{\mathcal I(a,b(h))=
       \frac{\vartheta_-}{2}\operatorname{Var}_a(\mathcal E)h^2
       +O(|h|^3).}
\tag{43.11}
\]

**Proof.** The tuple matching identity yields
`B(a tensor b,S(a tensor b))=B(a,b)B(b,a)=|B(a,b)|^2`.
Inserting P_-=(I-S)/2 proves (43.10). The native pairing inequality
derived in R41 bounds the overlap between zero and one. Direct
substitution proves the endpoint cases. Positivity and complementarity
of P_-,P_+ give the full-state ceiling, and the antisymmetric state
lies entirely in P_-. Finally apply R42.4's controlled overlap expansion
to (43.10), proving (43.11).

This cost is a conserved exchange-generator expectation. It differs
from R42's sum of prediction-residue squares and from R33's spatial
cut-elimination response. No equality between those quantities or
identification with a physical coupling constant is inferred.

### R43.6 — Unresolved matching derives the full reduced-pair law and its memory deficit

For normalized product inputs define A=aa-dagger, B=bb-dagger and
`kappa=Tr(AB)=|B(a,b)|^2`. Let c=cos_Sigma(phi), s=sin_Sigma(phi),
phi=theta_- tau/2. After W(tau), match out the second or first ledger
as in R20 to obtain rho_A and rho_B. Then

\[
\boxed{\rho_A=c^2A+s^2B+\iota cs[B,A],\qquad
       \rho_B=s^2A+c^2B-\iota cs[B,A].}
\tag{43.12}
\]

In particular rho_A+rho_B=A+B. A mixing formula that omits the
commutator generally gives the wrong native readout.

Define the specific reduced-pair memory deficit
`M=1-Tr(rho_A^2)`. The two ledgers have equal deficits, and

\[
\boxed{\mathcal M=2c^2s^2(1-\kappa)^2
       =2\sin_\Sigma^2(\vartheta_-\tau)
                  (\mathcal I/\vartheta_-)^2.}
\tag{43.13}
\]

It reaches `(1-kappa)^2/2` at the half tick and returns to zero at
the full swap. These product-input formulas are not asserted for
arbitrary already correlated preparations.

**Proof.** The irrelevant common phase cancels from the joint pair.
Expand `(cI+iota sS)(A tensor B)(cI-iota sS)`. Direct matching of
the second index gives A from A tensor B, B from its double swap,
BA from S(A tensor B), and AB from (A tensor B)S. This proves the
first formula. Matching the first index gives the opposite cross term
and proves the second. No partial-trace or statistical independence
law is imported; these are finite coefficient sums on native tuples.

For the deficit, A-square=A, B-square=B and Tr A=Tr B=1.
Finite cyclic reindexing gives Tr AB=kappa, Tr ABAB=kappa squared,
and `Tr([B,A]^2)=2(kappa^2-kappa)`. The cross terms between either
A or B and [B,A] have zero trace. Expanding the square in (43.12)
therefore gives
`Tr(rho_A^2)=1-2c^2s^2(1-kappa)^2`. The same calculation applies
to rho_B. Substitute (43.10) and the derived circular double-angle
identity to obtain (43.13).

For a pure full state, positive deficit certifies failure of product
factorization: a normalized product would have a rank-one local pair
whose square equals itself. The full state remains reversible under
W(-tau). The deficit is a quadratic observation of omitted pair
information, not a newly postulated destruction or stochastic law.

For the complete U(tau), both local pairs are conjugated by the same
free native process. Their deficits and the initial overlap are
unchanged. The displayed pairs themselves are the exchange-frame
expressions, not an omission of the free evolution.

### R43.7 — The actual autonomous clock preserves exchange balance while a common action scale remains free

Use U=(V tensor V)S as the finite process on the last edge of R41's
actual autonomous echo/controller word. Every complete tick then
applies U, and all lifted E_0, E_int and E_lift moments are exactly
conserved at every controller update. The inherited carrier-return
error does not become an exchange-balance error. R36's complete source
observation, extended by both process ledgers and the clock marks,
preserves the same identities and matching deficits.

The implementation (43.1) contributes three native record-control
factors. For the source-H pair, the full process word additionally
contains one H on each ledger. Grouping these five factors into one
controller edge does not erase their word count or assign them a
physical duration.

With a common native time calibration t=a tau and an exchange-consistent
common readout factor b, define `E_phys,total=b E_lift`. The equation is

\[
\boxed{\iota(ab)\frac{d\psi}{dt}=E_{\rm phys,total}\psi.}
\tag{43.14}
\]

R43.3 fixes b_A=b_B for the nonconstant identical resource readout,
but neither a nor the remaining common b is selected physically here.

**Proof.** All earlier controller edges act as identity on the process
ledgers; the last applies U. Each of the three readouts commutes with
U by R43.2. Direct tuple multiplication therefore proves commutation
with the entire fixed controller arrow. This is exact independently
of how closely the prepared source carrier returns to its initial
field. Complete observation and decoder factors cancel in consecutive
arrows and preserve native matching, as already proved in R36/R42.
The reduced process pair is also preserved by matching out a source
that has undergone that complete isometry.

The native derivative chain rule under t=a tau yields (43.14).
Every common positive b leaves exchange conservation intact, so the
relative-calibration result does not determine the dimensional product
ab. A positive half-turn branch and the exchange allocation are stated
native constructions. Their availability is not a theorem that nature
selects this process, physical energy, a material time unit, h-bar, c
or alpha. No unrestricted remote interaction is asserted: these are
internal co-carried process roles of the constructed clock target.

## Certification and premise-labelled reuse

[The native application](../04-operator-evolution/native_exchange_calibration.cjs)
uses the unchanged canonical engine. It checks the control-word exchange,
finite cut ranks, additive branch witnesses, all-state calibration defects,
closed current equations, distinguishability costs, reduced-pair and
memory identities, bounded factorial reconstructions and actual clock
intertwinings. Irrational phases are supported by R42's written completion
proof and rational enclosures; they are not replaced by exact decimal data.

[Verification](../04-operator-evolution/R43_VERIFICATION.json) binds seven
written results, nine exact groups and fourteen native word replays, with
fifteen boundary controls and nine rejected graph mutations. Frozen R42
is replayed and all 327 earlier non-navigation files are preserved.
The [ledger](../04-operator-evolution/R43_DERIVATION_LEDGER.json) labels the
exchange, identical-process identification, phase branch and readout
targets. Metadata auditing is not semantic or proof-assistant verification.

| Pinned source | Premise status and use |
|---|---|
| [R16 C1, C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** cut roles, H/K, iota, matching and completion; finite counting remains explicit infrastructure. |
| [R20.1–R20.2, R20.6, R20.8](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native-derived on explicit tuples:** controlled exchange, specified unresolved matching, returning pair memory and reversible decoding. No physical tensor, measurement or randomness axiom is a premise. |
| [R33](https://github.com/Parveen117/extra-ideas/blob/5ddce4544de7111a80c93a8ade7cf2399c7c90f3/02-relational-response/NATIVE_INTERACTION_MEMORY_R33.md) | **Comparison of observation targets only:** its spatial return coefficient and elimination cost are not premises or identifications for the R43 exchange cost. |
| [R36](https://github.com/Parveen117/extra-ideas/blob/011b8224c92d95f7f726eb56bab3e8d8e6f7280c/02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md) | **Native-derived:** complete curvature observation and its matching decoder, without a Riemann adapter premise. |
| [R41](https://github.com/Parveen117/extra-ideas/blob/f0d8f254416dd05e8a837871e190ec905f0b0d26/02-relational-response/NATIVE_RELATIONAL_CLOCK_R41.md) | **Native-derived with explicit process/controller target:** actual clock coupling, carrier error and native pairing bound. |
| [R42.1–R42.7](https://github.com/Parveen117/extra-ideas/blob/0354de8bc1c3474568b1bb8f363c757a1cefcbf6/02-relational-response/NATIVE_PHASE_GENERATOR_R42.md) | **Native-derived:** phase generator, count integration, factorial flow, half-turn, response metric and branch/calibration boundaries. R43 retains the selected branch and derives its exchange consequences. |
| [Canonical engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged exact arithmetic/replay:** finite scoped evidence, not a whole-engine, formal-assistant or physical-validation claim. |

```bash
python3.12 -B 04-operator-evolution/verify_r43.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```
