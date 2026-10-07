# SC1 — The scale factor made local: a rate of loss, not a field; and why scale is the unit

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic.

GB1 found three factors in a block — phase, unit block, scale — and
physics for two of them: the phase made local gives light, the unit
block made local carries gravity's field. This stage makes the third one
local and reports what comes out.

Sources read before building: **GB1-B1** (the three factors and their
action on a reading), **OB1** U2–U5 (potential, field, conserved reading),
**IN1-T4/T5** (the invariant; frames are unit blocks), **QC3-T1/T2**
(even part real, odd part circular), **CL1-T5 / R41.2** (history
residue), and the owner's principle recorded in IN1: an invariant is an
invariant whatever the frame.

## 1. What a local scale does

Let W = w₀ + w·C be a self-dagger potential and the law (D + W)ψ = 0.

**S1 (the count is not kept).**

```text
∂_t n + ∇·r = − 2 ( w₀ n + w·r ) .
```

A scale potential is a local rate of gain or loss of the count.

**S2 (the phase keeps it).** With ι·q·𝒜 in place of W, the right-hand
side is zero identically. The two potentials are the two parities of
QC3: the phase potential is the circular one and keeps everything; the
scale potential is the real one and loses.

**S3 (an exact scale potential is a relabelling).** D(fχ) = (Df)χ + f·Dχ
for a scalar f, so W = −(Df)/f only rescales the reading of a free
solution; and for W = Dσ the scale field (∇w₀ − ∂_t w ; ∇×w) is zero
identically.

**S4 (a scale field makes the invariant depend on the history).** With a
scale field that is not zero, two histories with the same ends carry
different scale factors. Witness: link scales 2, 3 along one way round a
square and 1, 2 along the other. On a record,

```text
count        differs by the factor  9   = (ratio)²
invariant    differs by the factor  81  = (ratio)⁴
```

Two unit-block histories never change the invariant. A single reading
stays null under any scale: only records can tell.

## 2. Verdict

```text
local scale with a field       the invariant of a record depends on its history        REFUSED
local scale, exact             a relabelling of the count                              nothing physical
a self-dagger potential W      a local rate of loss or gain                            a medium, not a fundamental field
global scale                   the size of "one" for a count                           the unit
```

The refusal is native. IN1 states that the unrecoverable part is the same
in every frame, and GB1 that frames are unit blocks; a scale field would
make two readings with one history each, and equal in every other way,
disagree on the invariant. That is exactly what "an invariant is an
invariant" forbids. (In nature it would mean that two atoms of one kind
could differ in mass by their past; they do not — general knowledge.)

## 3. The three factors, complete

```text
factor         made local                       what it is in physics                 what a reading sees
phase          allowed; field F = E + ιB        light; content q is the charge        nothing, except through F
unit block     allowed; field of frame changes  gravity (thesis); boosts and turns    n and r move, invariant kept
scale          refused as a field               the unit of count                     n, r and the invariant all rescale
```

So the block has room for exactly two fields of frame change — one
weighted by content, one not — and one number that is not a field at
all. That number is where every dimensional constant of the line has
been waiting (QC2, QC5, HB1, EM1). It is not fixed by a law because a
linear law cannot fix an amplitude; it is fixed by what is counted as
one.

## 4. Certificate

S1–S2 on a polynomial column with polynomial potentials, by the balance
identity. S3 by the product rule and by the vanishing of the scale field
for an exact potential, with a witness potential whose field is not
zero. S4 on a rational record: two scale histories, two unit-block
histories, and a single reading. Five tests; a phase potential that is
not self-dagger is caught.

## 5. Claim boundary

```text
COUNT BALANCE WITH A SCALE POTENTIAL; PHASE KEEPS THE COUNT                       PROVED
EXACT SCALE POTENTIAL = RELABELLING; NO FIELD                                     PROVED
SCALE FIELD ⇒ HISTORY-DEPENDENT INVARIANT (COUNT ², INVARIANT ⁴)                  PROVED
A FUNDAMENTAL SCALE FIELD                                                         REFUSED (IN1 + GB1; and by identical atoms)
THE BLOCK HOLDS EXACTLY TWO FIELDS OF FRAME CHANGE AND ONE UNIT                   PROVED as a statement about its three factors
A VALUE FOR THE UNIT; A RELATION AMONG CONSTANTS                                  NO — and the stage says why no law in the block can give one
```

A local scale as a candidate field, and its failure on history-dependent
lengths, are known (general knowledge). The stage's content is that the
framework refuses it by its own invariance statement, and that the
refused part is the unit.

## 6. Reproduce

```text
python sc1_local_scale.py
python -m unittest test_sc1_exact
```
