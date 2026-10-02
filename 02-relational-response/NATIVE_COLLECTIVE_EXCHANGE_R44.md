# R44: collective native exchange, calibration geometry and higher memory

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R43 constructed an exchange from native record controls and related its
phase cost to distinguishability and pair memory. R44 joins such exchanges
on a finite list of ledger pairs. It derives a collective flow, common
calibration on connected components, local current balance, and a strict
failure of closure at the level of all two-ledger observations.

The same finite difference form appears twice: it measures calibration
failure and drives propagation of a single retained cut. For an edge list
G, let L_G be the matrix whose action is derived below. Then

\[
\boxed{(L_Gu)_a=\sum_{b:\{a,b\}\in G}(u_a-u_b),\qquad
 \mathcal E_G\big|_{\text{one cut}}=\frac{\vartheta_-}{2}L_G.}
\tag{44.A}
\]

No classical graph Laplacian, interaction Hamiltonian, probability law,
physical spatial metric or field equation is a premise. The finite edge
list, identical-role identification and equal R43 exchange allocation are
explicit native target constructions. Their physical selection is not
claimed. The native half-turn and analytic completion are inherited from
R42, with their branch information retained.

## Native target and notation

Take N identified d-role ledgers, N>=2 and d>=2, with the R20 tuple
matching. G is a specified simple unordered list of distinct pairs
{a,b}; m is its length. Isolated ledgers are allowed. A path is a finite
sequence of pairs sharing endpoints; connected components and distance
are defined by those finite sequences, without a metric-space axiom.

Let S_ab exchange the two specified marks and leave all other marks
unchanged. Put

\[
 P_{ab}=(I-S_{ab})/2,\qquad
 G_\Sigma=\sum_{\{a,b\}\in G}P_{ab},\qquad
 \mathcal E_G=\vartheta_-G_\Sigma.
\tag{44.1}
\]

Each P_ab is the R43 reversal cut. The sum in (44.1) selects one equal
exchange allocation per listed pair; no new adjustable edge weight is
inserted. Repeated allocations would be a different recorded target.
The phase parameter tau is native count interpolation, not a material
duration. All expectations are normalized by native matching norm square.

Write A_a for a specified self-dagger d-role readout acting at ledger a,
and C(b)=sum_a b_a A_a for positive native multipliers b_a. The coefficient
norm is the finite sum `||M||_F^2=sum_ij |M_ij|^2`. Dagger, trace, scalar
completion, iota and factorial flow retain their earned source meanings.

## Written results

### R44.1 — Native exchange words construct a collective flow with a count-refinement bound

The completed collective arrow is

\[
\boxed{U_G(\tau)=\operatorname{Exp}_\Sigma(-\iota\tau\mathcal E_G).}
\tag{44.2}
\]

It preserves matching and exactly conserves E_G and every identical
readout sum sum_a A_a. Its generator is the sum of retained R43 branches,
not a fresh endpoint-only branch selection.

For any fixed order e_1,...,e_m of the edge list, set
`F(h)=W_(e_m)(h)...W_(e_1)(h)`, where W_e is R43's exchange flow. Then

\[
\boxed{\|F(\tau/n)^n-U_G(\tau)\|
 \le \frac{m^2\vartheta_-^2\tau^2}{n}
       \operatorname{Exp}_\Sigma(m\vartheta_-|\tau|/n).}
\tag{44.3}
\]

The norm here is the gain bound from native matching. The estimate is
zero when the edge list is empty, with both arrows equal to identity.

**Proof.** R43 constructs every swap from matching permutations and,
for two source roles, three native record controls. Lifting these arrows
by identity preserves their cut and matching identities. Their finite
sum is self-dagger. R42's bounded finite-matrix factorial argument therefore
constructs (44.2), its derivative, inverse and all generator moments.
The swap identity `S_ab A_a=A_b S_ab` proves that each edge commutes with
sum_a A_a; the factorial series preserves that commutation.

For (44.3), put x=m theta_- |h|. Expand both F(h) and U_G(h) in their
earned factorial series. The constant and first-order terms agree. The
sum of coefficient gain bounds for either expansion, from order two on,
is at most `Exp_Sigma(x)-1-x <= x^2 Exp_Sigma(x)/2`; the last inequality
follows term by term from `(j+2)!>=2 j!`. Thus the one-step difference
has gain at most x^2 Exp_Sigma(x). Telescope the difference of the nth
powers. All exact factors preserve matching, so no additional gain is
introduced. Multiplication by n proves (44.3) and its count limit.

For executable rational approximants, an edge may instead use the native
inverse expression `C_e(h)=(I+iota h theta P_e/2)^(-1)
(I-iota h theta P_e/2)`. It is an exact matching isometry. Its phase on
P_e is `2 integral_0^(h theta/2) (1+s^2)^(-1) ds`. Comparing this with
h theta, using the native count integral, gives phase error at most
`|h theta|^3/12`. Replacing all mn factors adds at most
`m theta^3 |tau|^3/(12 n^2)` to (44.3). Replacing theta_- by a rational
bracket midpoint with error epsilon adds `m epsilon |tau|` to the
completed collective target. R42 supplies the remaining factorial tail.
These are derived bounds, not an imported product-limit theorem.

A finite ordered word generally fails to conserve E_G, although it
conserves every identical-readout sum exactly. The collective limit and
the finite word are not identified. Full swaps retain their three-control
cost; fractional exchange arrows are counted as the derived factors
they are. No physical duration per refined factor is assigned.

### R44.2 — Collective conservation fixes component calibration through a derived difference form

Let A be non-scalar and self-dagger, and define its native centered size

\[
 v(A)=\operatorname{Tr}(A^2)-(\operatorname{Tr}A)^2/d>0.
\]

Then

\[
\boxed{\|[\mathcal E_G,C(b)]\|_F^2
 =\frac{\vartheta_-^2d^{N-1}}{2}\,v(A)
       \sum_{\{a,b\}\in G}(b_a-b_b)^2.}
\tag{44.4}
\]

Consequently C(b) is conserved under the collective flow for every
preparation exactly when b is constant on each connected component.
The same condition is equivalent to conservation under every listed
exchange separately. On a connected target there is one common factor;
disconnected components retain independent factors. Scalar readouts and
offsets proportional to identity do not constrain them.

**Proof.** R43 gives
`[S_ab,C(b)]=(b_a-b_b)(A_b-A_a)S_ab`. Its squared coefficient norm is
`2 d^(N-1) v(A) (b_a-b_b)^2`: the two-ledger calculation acquires exactly
d^(N-2) spectator copies.

Distinct edge defects are orthogonal in the coefficient matching. For
disjoint edges the trace factors and a single commutator trace is zero.
For edges sharing one endpoint, expanding their two readout differences
leaves four traces around the resulting three-cycle. Direct summation
of its matched indices gives Tr(A^2) for each term, with signs +,-,-,+,
and hence zero; untouched ledgers only multiply by their count. This
argument uses finite index matching, not a supplied trace-permutation
or diagonalization theorem. Since `[G_Sigma,C]=-sum_e[S_e,C]/2`, summing
the orthogonal squared defects proves (44.4).

R43's entrywise centered-square expansion proves v(A)>0 for non-scalar
A. The right side of (44.4) vanishes precisely when endpoint factors
agree on every edge, which is equivalent to agreement along every
finite path. Conservation for every preparation is equivalent to the
commutator vanishing: differentiate the native quadratic readout and
test basis marks, sums and iota contrasts, as in R43.3. Conversely a
vanishing commutator commutes with every factorial term. This proves
both directions without a possible cancellation between different edges.

Thus the relative calibration theorem now holds for the complete
collective generator, not merely an assumed collection of pairwise tests.
The common scale itself remains free.

### R44.3 — Local exchange derives a continuity law and an explicit current hierarchy

For each oriented listed edge define

\[
 F_{ab}(A)=\iota[P_{ab},A_a],\qquad F_{ba}(A)=-F_{ab}(A).
\]

These are self-dagger native readouts. Under (44.2),

\[
\boxed{\frac{d\langle A_a\rangle}{d\tau}
 =\vartheta_-\sum_{b:\{a,b\}\in G}\langle F_{ab}(A)\rangle.}
\tag{44.5}
\]

For any finite subset of ledgers, its total derivative is the sum of
currents on edges crossing its boundary. Internal edge contributions
cancel exactly. The next current derivative is

\[
\boxed{\frac{d\langle F_{ab}(A)\rangle}{d\tau}
 =\iota\vartheta_-\sum_{e:\ e\cap\{a,b\}\ne\varnothing}
      \langle[P_e,F_{ab}(A)]\rangle.}
\tag{44.6}
\]

Thus a neighbouring exchange can require three-ledger information.

**Proof.** Disjoint native tuple actions commute. The identical-readout
sum on the two endpoints commutes with P_ab, proving antisymmetry of
F. Daggering its commutator proves self-daggerness. Apply R42's native
readout derivative to E_G and A_a: every edge not incident to a commutes
with A_a, yielding (44.5). Sum over a subset and pair the two orientations
of each internal edge. This derives the boundary law without a primitive
divergence or continuum field equation. Applying the same derivative to
F gives (44.6); disjoint terms vanish by tuple matching.

For two source roles take Q=(I-H)/2. Then F_ab(Q) is exactly one half
of the R43 pair-current arrow, lifted to its endpoints. This connects
the collective continuity law to the already derived native transfer.

### R44.4 — Three-ledger ordering curvature is invisible to all pair observations in a native witness

Take three source-role ledgers with edges {1,3} and {2,3}, and retain
one additional record r. For eta in {1,-1,iota,-iota}, prepare the native
normalized vector

\[
\boxed{\Psi_\eta=\tfrac12\big(
 |010;0\rangle+\eta|100;0\rangle+
 |011;1\rangle-\eta|101;1\rangle\big).}
\tag{44.7}
\]

Every two-ledger pair readout among ledgers 1,2,3 is independent of eta,
after matching out the remaining ledger and record. Yet, with
F_23=F_23(Q),

\[
\boxed{\langle F_{23}\rangle=0,\qquad
 \frac{1}{\vartheta_-}\frac{d\langle F_{23}\rangle}{d\tau}
 =-\frac{\operatorname{Re}_\Sigma\eta}{4}.}
\tag{44.8}
\]

Define the three-ledger ordering-curvature readout

\[
 \Omega_{13,23}=\iota[P_{13},P_{23}].
 \qquad
 \boxed{\langle\Omega_{13,23}\rangle
       =\frac{\operatorname{Im}_\Sigma\eta}{4}.}
\tag{44.9}
\]

Thus eta=+/-1 gives identical pair data and opposite next-current
responses; eta=+/-iota gives identical pair data and opposite curvature.
This disproves closure by all pair observations for this declared target.

**Proof.** The four tuple marks in (44.7) are mutually matching-orthogonal,
so the factor 1/2 normalizes the vector without an assumed probability
mixture. For the 1,2 readout, the two record sectors contribute opposite
off-diagonal coefficients, which cancel; the result is
`(|01><01|+|10><10|)/2`. For the 1,3 and 2,3 readouts, the unretained
first or second label kills the cross terms; each resulting pair is
I/4. These three pair matrices are independent of eta.

Apply the two swaps to the four displayed marks and form their cuts.
The coefficient sums give `<iota[G_Sigma,F_23]>=-Re_Sigma(eta)/4`
and `<iota[P_13,P_23]>=Im_Sigma(eta)/4`, while the current itself is
zero. Equations (44.6), (44.8) and (44.9) follow. The executable witness
retains the full four-ledger preparation and matches the actual indices;
no independent mixed-state or random-selection axiom is introduced.

This is also genuine ordering curvature. Write X=S_13 S_23. Direct
tuple permutation gives X^3=I and

\[
 \Omega_{13,23}=\frac{\iota}{4}(X-X^{-1}),\qquad
 \Omega_{13,23}^2=\frac{2I-X-X^{-1}}{16}
 =\frac{3}{16}P_{\rm circ},\quad
 P_{\rm circ}=I-\frac{I+X+X^{-1}}{3}.
\tag{44.10}
\]

Since X^3=I, direct multiplication makes `(I+X+X^(-1))/3` a self-dagger
cut onto the X-fixed sector. Its complement P_circ is therefore a cut.
The coefficient 3/16 is fixed by the three-cycle count and the two
exchange-cut halves; it is not identified with a physical coupling.
The squared expression vanishes on the fixed sector and is nonzero on
the remaining cycle sector. Expanding the four native factorial factors in the small loop
`W_13(s) W_23(t) W_13(-s) W_23(-t)` gives
`I-theta_-^2 s t[P_13,P_23]` through mixed degree two; the higher terms
carry factorial bounds. The coefficient is fixed by source composition,
not an independently supplied curvature tensor. The full-swap loop is
X^2, generally nonidentity. A loop in operation order need not be a cycle
in the ledger-pair list: this example's list is a tree.

The witness does not prove that three-ledger data suffice for every
larger network. Repeated response can involve more retained ledgers.

### R44.5 — Single-cut propagation and calibration share the same derived quadratic geometry

Use N two-role ledgers and restrict to marks with exactly one Q role.
Let |a> denote that role at ledger a. The swaps preserve this sector,
and direct action gives

\[
\boxed{\mathcal E_G\psi=\frac{\vartheta_-}{2}L_G\psi,
 \quad (L_G\psi)_a=\sum_{b:\{a,b\}\in G}(\psi_a-\psi_b).}
\tag{44.11}
\]

In particular

\[
\boxed{B(\psi,\mathcal E_G\psi)
 =\frac{\vartheta_-}{2}
    \sum_{\{a,b\}\in G}|\psi_a-\psi_b|^2.}
\tag{44.12}
\]

The kernel consists exactly of coefficients constant on each component.
Its dimension is the component count. On a connected target the constant
single-cut mode is the unique zero-cost line.

**Proof.** If a single Q mark is outside an edge, that exchange fixes it.
On the edge endpoints, `P_ab|a>=(|a>-|b>)/2` and
`P_ab|b>=(|b>-|a>)/2`. Summing proves (44.11) and derives the degree
diagonal and negative edge entries of L_G. Matching the result with psi
and collecting each unordered pair once gives (44.12). Every summand
is nonnegative by the native scalar norm. Zero cost is equivalent to
equality of endpoint coefficients, hence constancy along finite paths.
Independent component-constant vectors give the stated dimension.

For a real multiplier vector b, the same direct expansion gives
`b^T L_G b=sum_edges(b_a-b_b)^2`. Substitution in (44.4) proves the
shared difference geometry. It measures both failure of relative
calibration and phase cost of a single retained cut. This is a derived
finite quadratic form on the declared edge list; it does not by itself
select physical space, a material metric or a physical energy unit.

### R44.6 — A cyclic native target has exact collective modes, a quadratic gap scale and controlled propagation tails

For the simple ring of N>=3 ledgers, let
`z_m=Exp_Sigma(2 iota theta_- m/N)` and `u_m(a)=z_m^a`.
Direct native matching yields N orthogonal modes with

\[
\boxed{\mathcal E_Gu_m=\varepsilon_m u_m,\qquad
 \varepsilon_m=\vartheta_-(1-\cos_\Sigma(2\vartheta_-m/N))
 =2\vartheta_-\sin_\Sigma^2(\vartheta_-m/N).}
\tag{44.13}
\]

The smallest positive phase gap is

\[
\boxed{\gamma_N=2\vartheta_-\sin_\Sigma^2(\vartheta_-/N),\qquad
 \frac{N^2\gamma_N}{2\vartheta_-^3}\longrightarrow1.}
\tag{44.14}
\]

For N>=4 the ratio lies in `[1-16/(3N^2),1]`. The ring therefore has
a derived quadratic small-phase scale; the finite gap is not a positive
universal mass gap and tends to zero as the declared count size grows.

More generally, let Delta be the maximum edge degree of a finite simple
target and let ell be the edge distance from a to b in one component.
The single-cut collective arrow obeys

\[
\boxed{|(U_G(\tau))_{ba}|
 \le\frac{(\vartheta_-\Delta|\tau|)^\ell}{\ell!}
           \operatorname{Exp}_\Sigma(\vartheta_-\Delta|\tau|).}
\tag{44.15}
\]

Between distinct components the entry is identically zero. For distinct
connected a,b, if p_ab counts shortest edge paths, its leading term is
`(iota theta_- tau/2)^ell p_ab/ell!`. Thus this analytic flow generally
has small nonzero tails at every reachable distance; it is not assigned
R38's exact finite-word front or a physical light speed.

**Proof.** Substitution in the degree-two difference formula (44.11)
gives `theta_-(1-(z_m+z_m^(-1))/2)u_m`, proving the eigenvalue formula
using the earned circular addition identities. The modes are distinct:
R42's integral `2 integral_0^s (1+x^2)^(-1) dx` for 0<=s<=1 increases
strictly from zero to theta_-/2, and its Cayley circle coordinate is
injective on that arc. The four native quarter-turn rotations cover
the circle without another return to 1 before 2 theta_-. Thus z_m/z_k
is not 1 for m!=k. Summing the finite geometric progression gives zero
for their matching and N for a mode's own norm. N nonzero orthogonal
vectors span the N-mark sector by native finite elimination. No Fourier
or spectral theorem has been imported.

The same quarter arcs and the factorial derivative show that sin is
nonnegative on [0,theta_-] and increases on [0,theta_-/2]; reflection
then selects m=1,N-1 as the least positive eigenvalue. For 0<=x<=1,
the alternating native factorial terms give
`x-x^3/6 <= sin_Sigma(x) <= x`. Squaring and dividing by x^2 gives
`1-x^2/3 <= (sin_Sigma(x)/x)^2 <= 1`. Use x=theta_-/N and the inherited
theta_-<4 to obtain the displayed bound and limit.

For (44.15), the absolute row and column sums of L_G/2 are bounded
by Delta. Applying native matching Cauchy to each row and then summing
columns gives gain at most Delta. A product of j matrix entries can
only stay at a vertex or move along one listed edge j times. Therefore
`(L_G^j)_(ba)=0` for j<ell. Bound the remaining factorial terms by
`sum_(j>=ell) x^j/j! <= x^ell Exp_Sigma(x)/ell!`, using
`(ell+r)!>=ell! r!`. At the first nonzero power, no diagonal stay is
possible and each shortest path contributes `(-1/2)^ell` to
`((L_G/2)^ell)_(ba)`. This proves both the leading coefficient and
the tail statement directly from native words and count completion.

The finite refined word in R44.1 has its own exact word-depth support.
Its depth mn grows in the fixed-tau completion. Treating that limit as
a fixed-duration primitive control would discard the cost ledger.

### R44.7 — The actual clock and complete observation retain collective invariants and the calibration boundary

Use U_G(1) as the finite process on the last edge of R41's actual
autonomous source/controller. Then every controller update preserves
E_G and every identical-readout sum exactly. If a finite ordered exchange
word is used instead, the same sums remain exact but E_G need not.
This distinction is retained in the executable clock witnesses.

R36's complete source observer, extended by all process and record marks,
preserves the full process pair, the three-ledger witness and its current
and curvature values. Omitting a process ledger remains a different,
potentially blind observation.

For positive common calibrations t=a tau and E_phys=b E_G,

\[
\boxed{\iota(ab)\frac{d\psi}{dt}=E_{\rm phys}\psi.}
\tag{44.16}
\]

On each connected component R44.2 fixes relative multipliers; it does
not fix a, the common b, a physical length per edge or their action product.

**Proof.** Earlier controller edges act trivially on the process, and
the last applies the specified matching isometry. Every readout commuting
with that process therefore commutes with the complete fixed controller
arrow, exactly as in R41/R43. The carrier-return error has no role in
that commutation. A finite exchange word commutes with all identical
readout sums factor by factor, while a two-edge three-ledger word gives
a direct nonzero commutator with G_Sigma. This proves why its conservation
claims differ from the completed collective arrow.

The complete observer and decoder preserve source matching. Matching
out that source before or after the complete isometry therefore gives
the same full process pair, including a retained purification record.
Subsequent specified pair readouts and three-ledger operators consequently
agree. This does not invert an incomplete pair observation: the four
preparations of R44.4 remain its explicit blind witness.

Finally the native derivative chain rule gives (44.16), and every common
positive b on a component preserves (44.4). No physical h-bar, c, alpha,
mass, dimension or material interaction is selected by these internal
identities. The collective phase modes and finite gap are mathematically
derived on the chosen native target; their physical identification is a
separate unresolved step.

## Certification and premise-labelled reuse

[The native application](../04-operator-evolution/native_collective_exchange.cjs)
uses the unchanged canonical engine. It binds nine exact groups covering
native exchange lifts, finite count-refinement bounds, collective
calibration defects, current balance, pair-blind three-ledger witnesses,
ordering curvature, single-cut difference geometry, cyclic modes and
propagation bounds, and actual clock/complete-observer transport.

[Verification](../04-operator-evolution/R44_VERIFICATION.json) binds seven
written results and fourteen native word replays, with fifteen boundary
controls and nine rejected graph mutations. Frozen R43 is replayed;
all 334 earlier non-navigation files are preserved. The
[ledger](../04-operator-evolution/R44_DERIVATION_LEDGER.json) records native
target definitions and premise status. Finite exact checks and a declared
dependency audit are not semantic or proof-assistant verification.

| Pinned source | Premise status and use |
|---|---|
| [R16 C1, C3-C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** cut roles, H/K, iota, matching and count completion; finite counting is explicit infrastructure. |
| [R20](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native-derived on declared tuples:** controlled K, record preparation and unresolved index matching. No physical subsystem or randomness axiom is supplied. |
| [R36](https://github.com/Parveen117/extra-ideas/blob/011b8224c92d95f7f726eb56bab3e8d8e6f7280c/02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md) | **Native-derived:** complete curvature observer and matching decoder. Its classical Riemann adapter is not a premise. |
| [R38](https://github.com/Parveen117/extra-ideas/blob/79a4c4c2e15aad1391f7633ab0ab0b65ce580111/02-relational-response/NATIVE_UNIVERSAL_SIGNAL_CONE_R38.md) | **Comparison of target/cost conventions only:** its fixed-word exact front is not used as a premise or assigned to the completed R44 collective flow. |
| [R41](https://github.com/Parveen117/extra-ideas/blob/f0d8f254416dd05e8a837871e190ec905f0b0d26/02-relational-response/NATIVE_RELATIONAL_CLOCK_R41.md) | **Native-derived with explicit controller/process target:** actual autonomous clock coupling and preserved readouts. |
| [R42](https://github.com/Parveen117/extra-ideas/blob/0354de8bc1c3474568b1bb8f363c757a1cefcbf6/02-relational-response/NATIVE_PHASE_GENERATOR_R42.md) | **Native-derived:** factorial flow, count integration, half-turn, response differentiation and phase-branch boundaries. |
| [R43](https://github.com/Parveen117/extra-ideas/blob/8aa356cda8707e740cb713a680dc0f6aa2f8428f/02-relational-response/NATIVE_EXCHANGE_CALIBRATION_R43.md) | **Native-derived with explicit exchange/readout target:** source-controlled swap, interaction branch, calibration defect and pair current. R44 derives collective consequences without importing a graph operator or many-body law. |
| [Canonical engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged exact arithmetic/replay:** scoped application evidence; no whole-engine, formal-assistant or physical-validation claim. |

```bash
python3.12 -B 04-operator-evolution/verify_r44.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```
