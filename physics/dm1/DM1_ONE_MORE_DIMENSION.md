# DM1 — One more dimension reads the curvature: the hypotenuse, the unseen leg, and mass as the third momentum

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic over the
Gaussian rationals. Continues IN1.

Owner's intuition for this stage: the arc of a circle is measured by
Pythagoras; the hypotenuse is the complete information. But the circle
is curvature, and to read it one more dimension is needed. In two
dimensions we see only the flat line — the cosine, not the sine. The
sine is what gives the curvature, and it needs a dimension. A line
cannot be seen in one dimension: it needs thickness.

Sources read before building: **IN1-T1/T2/T4/T5** (reading tensor;
single reading null; frame invariance; records), **FR1-T3/T4** (two
cuts real, three cut-complex; the conjugate sheet), **PR1-T2/T4**
(speed² + (curvature/2)² = 1; the sector law with mass on R),
**PR2-T1/T3** (two sheets), **MS1-T1** (speed² + memory = 1).

## 1. The hypotenuse and its legs

**T1.** For a single reading, n² = r₁² + r₂² + r₃². An observer with the
two real cuts sees

```text
cos²α = ( r₁² + r₂² ) / n²        the seen leg — a speed
sin²α =   r₃² / n²                the unseen leg
```

and the two complete the hypotenuse n. For a real reading the unseen leg
is zero; for a cut-complex one it is not.

This is PR1's law for a coin and MS1's for a state — speed² + (what is
not seen)² = 1 — with the second term now identified: it is the reading
of the third cut.

**T4 (the ladder).** With one cut the unseen part is r₂² + r₃²; with two
it is r₃²; with three it is zero.

```text
reading (n; r₁, r₂, r₃)         one cut      two cuts      three cuts
(10; 8, 6, 0)   real            9/25         0             0
(5; 0, 3, 4)                    16/25        16/25         0
(15; −5, −2, −14)               8/9          196/225       0
```

## 2. The unseen leg is what the lower observer calls mass

**T2.** Every real frame change of determinant one changes n, r₁ and r₂
and leaves r₃ unchanged, and n² − r₁² − r₂² = r₃². For the two-cut
observer the reading is a massive one, and its invariant rest mass is
the reading of the cut that observer does not have. A cut-complex frame
change does move r₃.

**T3 (mass is the third momentum).** As operators,

```text
c ( ξ₁C₁ + ξ₂C₂ ) + c (ι k₃) C₃   =   c ( ξ₁K + ξ₂RK ) − c k₃ R .
```

The three-cut law with no mass, on a reading of third wave number k₃, is
the two-cut law of PR1 with flip rate g = −c·k₃. What turns in the plane
is what moves along the third cut. Reversing k₃ reverses the turn and
nothing else: the two sheets of PR2 — mass and antimass — are the two
directions along the cut that is not seen.

## 3. One level further

**T5.** A record ρ = AA† has n² − r·r = 4|det A|². It is zero exactly
when A is a single reading. On the doubled carrier — the frame together
with its second sheet — the record is a single reading again, and what
no cut of the first frame holds is the pairing between the two sheets.

So the pattern repeats: what is curvature, memory or mass at one level is
a plain reading at the next.

```text
seen with               looks like                          is, one level up
two real cuts           a turn, a mass, a speed below 1     a light-like reading tilted out of the plane
three cuts              a record with rest mass             a single reading on the doubled carrier
```

## 4. Reading

- The hypotenuse is the uncut. Each observer measures the leg it can and
  infers the other; Pythagoras is the conservation law of IN1.
- The cosine is what a frame reads; the sine is what it cannot, and what
  it cannot read it meets as curvature — as a turn, as mass, as a clock.
- Adding the cut makes the reading flat: nothing is lost, and everything
  becomes computable. That is why the frame with three cuts is the one
  in which single readings are light-like.
- A reading needs one more cut than it has to be read completely: the
  line needs thickness.

## 5. Certificate

T1, T4 on five rational readings (two real, three cut-complex). T2 on
three readings through four real frame changes (two boosts, a turn, a
shear) and one cut-complex change. T3 as an identity of cut-complex
matrices at three points, with the square law and the sheet reversal.
T5 on four coefficient blocks, one of rank one. Six tests; a cut-complex
frame passed off as real, and a third term without ι, are rejected.

## 6. Claim boundary

```text
SEEN² + UNSEEN² = UNCUT² ; THE LADDER                                           PROVED
r₃ INVARIANT UNDER THE TWO-CUT FRAME; TWO-CUT REST MASS = THIRD READING         PROVED
THREE-CUT MASSLESS LAW ON A k₃ MODE = TWO-CUT LAW WITH g = −c k₃ ; SHEETS       PROVED
RECORD: n² − r·r = 4|det A|² ; SINGLE ON THE DOUBLED CARRIER                    PROVED
THAT PHYSICAL MASS IS MOMENTUM IN AN UNSEEN DIRECTION                           NOT CLAIMED — the algebra allows the reading; nothing selects k₃
A VALUE OF ANY MASS; WHY k₃ WOULD BE FIXED OR DISCRETE                          NO
MO1's ASSUMPTION (PURE FRAMES SHARE ONE FLAT SPACE)                             NOT ADDRESSED
```

Reading mass as momentum along a direction that is not observed is a
known idea (general knowledge). Here it is a consequence of the block:
the turn R of the real frame is, on the cut-complex carrier, a cut.

## 7. Reproduce

```text
python dm1_one_more_dimension.py
python -m unittest test_dm1_exact
```
