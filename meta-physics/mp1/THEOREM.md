# MP-1: hidden frame memory and a certified observation protocol

Research owner: Monty Dabas. Development: 3 October 2026.

This is a scoped development of UAL source ideas 017 and 078, informed by
009. It certifies an axis observer, a retained sign, integer history and
an explicit tolerance to perturbations. It does not certify the source
PDFs as whole papers. The [source catalogue](../SOURCE_AUDIT.md) records
the full reading pass and the remaining work.

The operational result is

\[
  Q(g)=gKg^{-1}=Q(-g),\qquad
  g(1)=(-I)^{\nu}g(0),\qquad
  2(\epsilon+\eta)^2<m^2\ \Longrightarrow\ \nu_{\rm true}=\nu_{\rm polygon}.
\]

Here a closed **observed-axis loop** has winding nu, m is the measured
polygon's minimum distance from zero, epsilon bounds coordinate errors
at the vertices, and eta bounds the true path's departure from its own
chords. The last two bounds are declared inputs, not inferred from samples.

## Carrier, inputs and attribution

Use the existing EMK real two-role carrier, with `K^2=I`, `R^2=-I`,
`KR=-RK` and `L=RK`. Code imports R7's actual matrices and R2's exact
arithmetic unchanged. R7's K exchanges the two coordinates; we do not
replace it with another convention. Rational cut counts provide the
executable sector. For continuous statements use the previously earned
R16 completion and R42 phase interpolation. There is no new Hilbert-space,
spacetime, physical clock, probability law or entropy axiom.

A unit turn is `g=aI+bR`, `a^2+b^2=1`. The **declared observer** records
the labelled involution `Q=gKg^-1`, equivalently the labelled projection
`(I+Q)/2`. It forgets the overall sign of g, but does not identify Q with
-Q. Erasing the labels of the two eigenspaces would be a different
observer and is outside this packet. A nonzero raw coefficient pair
`z=(x,y)` is observed through `(xK+yL)/sqrt(x^2+y^2)`.

Winding, double covers and homotopy invariance are standard mathematics.
The circle/lifting comparison is credited to Hatcher, *Algebraic Topology*,
Theorem 1.7 and section 1.3 ([author's chapter](https://pi.math.cornell.edu/~hatcher/AT/ATch1.pdf)).
The proofs below specialize these mechanisms to the existing native
carrier and give the executable observation/error contract. No general
topological novelty is claimed.

Earlier native results are also used, not rediscovered:

| Source | Existing result | MP-1 addition or distinction |
|---|---|---|
| [R7](../../02-relational-response/EMK_TENSOR_CALCULUS_R7.md), [R14](../../02-relational-response/NATIVE_SOURCE_FOUNDATION_R14.md), [R15](../../02-relational-response/EMK_TOPOLOGY_R15.md) | Typed carriers; phase endpoint can erase a separate integer ledger | Derive the two-sheet fibre for a specific quadratic observer |
| [R39.5](../../02-relational-response/NATIVE_ECHO_CLOCK_R39.md) | Finite marked records cannot retain arbitrarily many counts | Apply the existing capacity argument separately to parity and integer winding |
| [R42](../../02-relational-response/NATIVE_PHASE_GENERATOR_R42.md) | Native phase interpolation and explicit branch contract | Sharp quarter-turn bound for recovery through this axis observer |
| [RKF T76](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/76_seam_flow_meter_theorem.md) | Pencil sector counts, determinant parity and finite-family protection | Different carrier/observable: a universal tube bound for polygon winding; not a strengthening of T76's pencil perturbation claim |

Extra Ideas dependencies are pinned at `6f131bab2b3b3cb813aed4d5c7d03b24bd1ae3d8`;
the verifier hashes the consumed local files and the frozen R register.

## T1 — The axis observer has exactly two frame sheets

For every unit turn,

\[
 Q(g)=(a^2-b^2)K+2abL,\qquad Q(g)^2=I.
\]

On this turn family, `Q(g)=Q(h)` if and only if `h=g` or `h=-g`.
In full-turn phase units t, its observed coefficient phase is `2t mod 1`.

**Proof.** Expand `(aI+bR)K(aI-bR)` using anticommutation. For equal
observations, `h^-1 g=cI+dR` commutes with K. Its commutator is `2dRK`,
so d=0; the unit identity forces c=+1 or -1. Conversely those two signs
give the same conjugation. The double-phase law is the multiplication
law of the earned native turn, applied to its square. This calculation
does not assume a classical metric curvature definition. Q is an observer
contrast, not the native curvature tensor.

## T2 — Rational samples give an exact integer for their polygon

Let `p_0,...,p_n=p_0` be rational coefficient pairs, with no segment
meeting zero. Self-intersections and repeated nonzero vertices are allowed.
For `a=p_j`, `b=p_(j+1)` and `D=det(a,b)`, give the edge contribution

\[
 c(a,b)=\begin{cases}
  +1 & a_y\leq0<b_y\ \text{and}\ D>0,\\
  -1 & b_y\leq0<a_y\ \text{and}\ D<0,\\
  0 & \text{otherwise}.
 \end{cases}\qquad \nu=\sum_j c(p_j,p_{j+1}).
\]

Then nu is the signed full-turn count of the normalized polygon. It is
additive on concatenation, negated on reversal, and unchanged by segment
subdivision or a deformation through closed nonzero paths.

**Proof.** A segment avoiding zero has a continuous phase branch whose
change is strictly between minus and plus half a turn. Across the
positive horizontal ray, the phase branch changes sheet by one signed
turn. At a transverse crossing its horizontal coordinate is
`D/(b_y-a_y)`, giving exactly the displayed signs. The half-open y
inequalities count a vertex crossing once; a contact returning to the
same side cancels. A horizontal edge on the ray contributes no turn and
its adjacent edges follow the same convention. Sum the local phase
changes: all intermediate endpoints cancel, leaving an integer because
the path is closed. Reversal and concatenation follow by reversing or
joining those changes, including the vertex convention.

For completeness, the required local phase/lift construction uses small
arcs only. On a unit arc `(x,y)` excluding `(-1,0)`, its square root with
positive first coordinate is
`(sqrt((1+x)/2), y/(2sqrt((1+x)/2)))`. Relative multiplication moves such
a chart to any basepoint. Finite subdivision of a uniformly continuous
nonzero path supplies these charts; choose successive signs to agree at
shared endpoints. Two agreeing lifts differ by a continuous sign and
therefore agree throughout. The same construction on small rectangles
of a continuous deformation makes its endpoint turn count continuous
and integer valued, hence constant. Roots/continuity here are in the
earned scalar completion. The executable count itself needs neither
roots, atan, floating point nor rounding.

## T3 — An observed return can hide a frame reversal

For a closed normalized axis path of winding nu, fix either initial
frame g_0 in its fibre. Its unique continuous frame lift satisfies

\[
                 g_1=(-I)^{\nu}g_0.
\]

One observed turn reverses the frame; two restore it. Both cases return
the observer to its initial Q. There is no continuous single-valued
choice of frame over the entire axis circle.

**Proof.** By T1 the frame phase changes by half the unwrapped axis
phase change. Thus its change is `nu/2` full turns, giving the sign.
Alternatively the matching local roots in T2 exchange sign at each
odd accumulated axis turn. If a continuous single-valued choice existed,
apply it to one observed turn: it would be a closed frame lift, contradicting
the odd sign just derived. This is a conditional observation mechanism:
observing does not force the source to execute the loop.

## T4 — A universal, computable protection radius

For an edge from a to b set `d=b-a`. If d is nonzero define

\[
 t_* = \min(1,\max(0,-a\cdot d/(d\cdot d))),\qquad
 m_{ab}^2=\|a+t_*d\|^2.
\]

For a repeated vertex use `m_ab^2=||a||^2`. Set `m^2=min_edges m_ab^2`.
These are exact rationals. Any closed continuous curve z with
`sup_t ||z(t)-p(t)|| < m`, using the same path parameter, has the same
winding and relative frame sign as p. In particular it suffices that
the true vertex errors are coordinatewise at most epsilon, the true
interpolation error from its own chords is coordinatewise at most eta,
and `2(epsilon+eta)^2 < m^2`.

**Proof.** Minimize the edge's quadratic squared norm: its unconstrained
minimizer is `-a.d/(d.d)` and clipping gives the segment minimum. The
linear deformation `(1-s)p+s z` stays at distance at least
`m-s sup||z-p|| > 0` from zero, so T2 preserves its count. Interpolating
the true endpoint errors costs at most epsilon in either coordinate;
adding the chord-departure bound costs eta. Squaring the resulting
two-coordinate bound proves the test. This protects against **every**
perturbation satisfying the contract, not just sampled noise draws.

Strictness is necessary for a uniform Euclidean guarantee: translating
the whole polygon by minus a closest point, a translation of norm m,
makes it meet zero. The coordinatewise sufficient bound can be
conservative. Failure of the test means uncertified, not a detected
transition. Nonzero vertices alone are insufficient: an edge may pass
through zero. Nor do samples alone certify an unseen continuous path;
a hidden loop can leave every sample unchanged and violate the tube bound.

## T5 — Axis observation halves the unambiguous sampling range

Let the actual unwrapped frame phases at consecutive samples be t_j,
and know the initial frame phase. If every actual increment obeys
`|t_(j+1)-t_j| < 1/4`, the axis samples `u_j=2t_j mod 1` recover all
frame sample phases uniquely: unwrap each axis increment into `(-1/2,1/2)`
and divide by two. For the directly observed frame phase the corresponding
bound is `|Delta t|<1/2`.

**Proof.** Two possible frame increments yielding the same axis increment
differ by an integer multiple of 1/2. The open interval `(-1/4,1/4)`
contains at most one representative; the actual bound supplies it.
The direct-frame argument uses integer differences instead. At the
boundary the two increments +1/4 and -1/4 give the same axis endpoint
but opposite frames. Thus the uniform bound cannot include equality.
Without an increment bound, even a hidden half-turn changes the frame
sign while returning the axis sample. The algorithm cannot test the
unwrapped bound from wrapped samples. A bound on total phase variation
within each interval is a stronger sufficient experimental contract.

## T6 — Retain the record required by the question

For loops based at the same Q and a known g_0:

| Recovery target using the axis endpoint plus a marked record | Necessary and sufficient extra record |
|---|---|
| Final frame g_1 | Two labels, `nu mod 2` |
| Winding nu restricted to `-N,...,N` | `2N+1` labels, the integer itself |
| Unbounded integer winding | No fixed finite alphabet of distinguishable labels suffices |

**Proof.** T3 supplies the parity decoder and shows that both signs
occur at the same observed endpoint. Repeating/reversing the diamond
polygon realizes every listed integer at that same endpoint. Distinct
required answers must receive distinct record labels, by R39's existing
count-capacity argument. An integer ledger attains the count. This is
not a bound on arbitrary real/coherent encodings; it is explicitly a
bound on distinguishable marked records. If the full frame endpoint
is retained already, it supplies parity and changes the additional
memory question.

The hidden sign can matter to a declared later reference probe. Keep a
nonzero reference v in a second labelled port and compare the matching
readout `||v+g_1 g_0^-1 v||^2`. It is `4||v||^2` for even nu and zero
for odd nu. Direct expansion proves this without a Born-rule premise.
Erasing the reference or quotienting all future probes by sign removes
this operational distinction. No global phase is declared observable
without a reference contract.

## T7 — Topological memory is not local non-Abelian curvature

A continuous family of closed nonzero axis loops cannot change nu.
A claimed change therefore requires a zero, failure of closure or
continuity, a changed observer, or a violated reconstruction contract.
A large scalar acceleration/skew-energy threshold alone does not prove
any of those events. Smooth changes of speed can make acceleration
arbitrarily large on a fixed loop while its winding is unchanged.

Moreover, this frame memory is compatible with locally flat transport.
On a patch with a single-valued frame choice h(q), define transport
`T(q_b,q_a)=h(q_b)h(q_a)^-1`. Products telescope, so every closed path
inside the patch gives I, including every sufficiently small rectangle.
Patching these local lifts around the puncture gives `(-I)^nu` by T3.
This is global sheet monodromy despite zero local loop defect. It does
not replace the existing native ordered-commutator curvature with circle
bending, entropy density or winding. In particular, all g=aI+bR commute.
MP-1 cannot by itself establish non-Abelian gauge dynamics.

## What is repaired, and what stays open

The source formula `exp(2 pi i nu)` equals one for integer nu and cannot
retain its value. MP-1 keeps the integer as a record and obtains a
different, two-sheet **frame** sign `(-I)^nu` through the axis observer.
This is a reconstruction of the source intuition, not validation of its
entire proposed detector. Fixing V also makes the pulled-back
`dT wedge dV` zero; it does not produce the source's claimed nonzero
area element. A single-valued exact phase differential has zero closed
integral; a wrapped phase/lift must be explicitly distinguished.

The source curvature-memory threshold, entropy cost, physical realization,
probabilistic noise law, thermodynamic arrow, awareness claims, quantum
gravity and Yang–Mills are unproved here. A useful next gate is to supply
an actual native observation protocol with a derived between-sample
error bound. The present certification ends at the declared observer,
completed phase path and exact error budget.

## Certification scope and replay

Seven written results, an exact rational verifier and adversarial tests
support this packet. They are not a proof-assistant certification,
independent peer review or an experimental result. Finite test cases
check the implementation; the universal statements rest on the written
proofs with their hypotheses.

From the repository root, using **Python 3.12 only**:

```sh
python3.12 -B -m unittest discover -s meta-physics/mp1 -p 'test_*.py' -v
python3.12 -B meta-physics/mp1/verify.py --check
```
