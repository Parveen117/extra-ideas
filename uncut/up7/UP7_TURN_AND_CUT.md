# UP7 — why the diagram exists: a turn and one self-dagger cut

Seventh stage of the uncut line. Nothing thermodynamic is used. Exact rational matrices on the primitive
carrier (sympy): `up7_turn_and_cut.py`, `test_up7.py` (6 tests).

The owner's statement: the diagram is only a diagram — no mathematics enforces it and it is not known why
nature follows it. The Morphic documents are the ones written for the uncut. The mathematics should say why
this exists in pure mathematics, whatever the thermodynamic diagram.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework `theorum/morphic_algebra/main.tex` | EMK-0: cut κ, curvature ι (ι² = −I), recognition χ = ι∘κ; EMK-3: the 4-cycle T →κ V →ι S →κ P →ι T; its Maxwell corollary; EMK-Rep; EMK-12 question 3 |
| Recognition-Kernel-Framework theorum/48, theorum/55 | K² = I, R² = −I, RK = −KR; dagger R → −R, K → K |
| response-geometry RMG2 T1(a), RMG9 T1 | X = pK + qS + rR, X² = (p² + q² − r²)·1 |
| extra-ideas `uncut/up4`, `up5`, `up6` | w; the four relations are one equation; the eight scale operations |

## Set-up

The primitive carrier has two real dimensions and a turn ι with ι² = −1. A cut is an operation κ on it with
κ² = 1, κ ≠ ±1: it has two sides. Write K for one self-dagger cut and S = ιK.

## Statements

**T1 (all cuts).** Every cut is κ = uK + vS + wι with u² + v² − w² = 1, and

  κι + ικ = −2w,  κ† = κ − 2wι.

A cut anticommutes with the turn exactly when w = 0, exactly when it is its own dagger. No cut commutes with
the turn: the only operations with square 1 that do are ±1.

**T2 (the four-step walk).** With χ = ικ,

  χ² = 1 − 2wχ.

The walk κ, ι, κ, ι returns to its start exactly when w = 0. Otherwise it misses by −2w times the state two
steps along, and never returns: for w = 3/4, χ has eigenvalues −2 and ½, and each further four steps scale by
4 and ¼.

**T3 (w = 0: the diagram).** The turn and one self-dagger cut generate a group of exactly eight operations,
with exactly four mirrors. Place a state T; then V = κT, S = ιV, P = κS, ιP = T — four ends on two
perpendicular lines, T opposite S and V opposite P. The four mirrors are those two lines and the two
diagonals; the ends of the diagonals are the four corners F = T+V, U = S+V, H = S+P, G = T+P. The eight
operations are:

| operation | on the ends (T V S P →) | on the corners (F U H G →) |
|---|---|---|
| identity | T V S P | F U H G |
| mirror in the V–P line (= χ) | S V T P | U F G H |
| mirror in the T–S line | T P S V | G H U F |
| half turn | S P T V | H G F U |
| mirror in the F–H diagonal (= κ) | V T P S | F G H U |
| mirror in the U–G diagonal | (the other exchange of the two lines) | |
| the two quarter turns | cycle the four ends | cycle the four corners |

The mirror in the V–P line exchanges one conjugate pair and keeps the other: it carries U to F and H to G.
The mirror in the T–S line carries F to G and U to H. These are the two partial changes of corner, and the
half turn is both. The diagonal mirrors exchange the two lines — the mirror between the thermal and the
mechanical side. The quarter turn is the step "next reading".

**T4 (eight).** The eight pairs (end, neighbouring end) are permuted by the eight operations with no
repetition. They are the eight scale operations E[x | y] of UP6. The count eight, the four reciprocal pairs
and the four corners are the order, the mirrors and the diagonals of this group.

**T5 (the same on any carrier).** Any self-dagger cut in place of K gives the same group. On a doubled
carrier a cut that commutes with the turn exists, and its walk needs eight steps, not four: the
four-step closure is the anticommutation.

**T6 (the response's own cut).** A response element L = [[A, B₁], [B₂, C]] is mean + √D·κ_L with κ_L a cut, and
κ_L ι + ι κ_L = −(B₂ − B₁)/√D. So κ_L is self-dagger exactly when the mixed responses agree — the one equation
of UP5-T1 — and its w is the w of UP4. The diagram closes for a response exactly when its cut is its own dagger.

## Answer to the question

The diagram is not a fact about heat. It is what a turn (ι² = −1) and one cut that is its own dagger generate
on the primitive carrier: two perpendicular lines, two diagonals, eight operations. Anything that carries
such a turn and such a cut carries the diagram. What is not forced is w = 0. A cut with w ≠ 0 is still a cut;
its walk opens into a spiral, and the amount by which it opens is the same 2w that UP4 found as the cost of
non-commuting cuts. "Nature follows the diagram" is the statement w = 0, and that statement is open.

## Relation to the Morphic manuscript

EMK-3 states this four-cycle and its ledger carries it as conditional. T1–T3 give the condition: the cut is
an operation of square 1 that is its own dagger; then four steps give (ικ)² = +1. The manuscript's proof
passes through χ² ≡ −I and χ⁴, and the representation it fixes for computations takes κ with κ² = κ, for which
(ικ)² = ικ and the walk stops instead of returning. Question 3 of EMK-12 — the four potentials from the
four-cycle — is answered by T3: they are the ends of the two diagonal mirrors, and the changes of corner are
the two line mirrors.

## What is put in, what is not claimed

* Two real dimensions, one turn, one cut of square 1. That a cut has square 1 (two sides) is taken as its
  definition here; the manuscript's requirement tr κ = 1 belongs to the aperture (1 + κ)/2, not to κ.
* Why the carrier has a turn at all, and why there is one cut and not several, is not derived.
* T3 gives the diagram's operations. That thermodynamic readings transform under them as labelled is the
  identification of UP6, not a new derivation here.
* Nothing selects w.

## Open gates

1. w: a statement that forces w = 0, or the structure for w ≠ 0 (the spiral; the manuscript's helical orbit).
2. Three cuts (FR1): the group in place of the eight operations, and the diagram it draws.
3. The sequence of UP3 as the action of the quarter turn on readings.
