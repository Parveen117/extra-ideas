# LC1 — The part lost to a cut is the part that is counted

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic.

LN1 ended on a statement it could not derive: that what a cut loses is
what gets counted. This stage derives it.

Sources read before building: **LN1-N1** (lost = (H²)_ii − (H_ii)²),
**QC1-T3** (content q is silent iff qΘ ∈ 2πℤ), **QC1-T4** (on a circle
of hyperbolic radius d: Θ = −π(cosh d − 1); silent circles
cosh d_n = 1 + 2n/q; heights k + n, k = q/2), **RMG1-T2**
(f = ¼ sinh(ℓ/2) dℓ∧dφ, with d = ℓ/2), **WQ1-W3**, **GR2** (once-around
turn 2π(1/N − 1)), **LT1-T3** (the level below).

## Setting

Normalised element H = Exp(d·n) and its root M = Exp((d/2)·n), both
read in a cut across the axis n:

```text
seen(H) = cosh d ,        lost(M) = sinh²(d/2) =: L ,        seen(H) = 1 + 2L .
```

## Results

**C1 — the return is the lost part.** Around a circle,

```text
Θ = −π ( cosh d − 1 ) = −2π · L .
```

In whole turns, the return of the element is the part its root loses to
a cut.

**C2 — the connection is the lost part.** A = L dφ, and
dA = dL∧dφ = ¼ sinh(ℓ/2) dℓ∧dφ: RMG1's curvature is the change of the
lost part against the turning of the axis.

**C3 — the count is the lost part.**

```text
content q is silent   ⇔   q · L is a whole number ,        L_n = n / q .
```

**C4 — the height is the seen part.**

```text
k · cosh d_n = k ( 1 + 2 L_n ) = k + n .
```

The ladder of QC1 is k times what is seen; n is q times what is lost;
and the lowest value, k = q/2, is what stays seen when nothing is lost.
For q = 1 the levels (n + ½) of WQ1 are half the seen part, with n the
lost part.

**C5 — gravity.** The once-around turn of GR2 is

```text
2π ( 1/N − 1 ) = 4π · sinh²(η/2) = 4π · x/(1 − x) ,     x = r_s/4s :
```

4π times the lost part at the lower level of the tower (LT1-T3). To
first order it is π r_s/r, the term named in the thesis.

## Reading

```text
observation in a cut      splits the element:   seen = 1 + 2·lost
the return of a cycle     is the lost part, in turns
recordability             makes the lost part a whole number over q
the level                 is k times the seen part ; its floor k is the seen part with nothing lost
```

The owner's sentence is now one chain of certified statements:
observation (a cut) produces a lost part; the tower returns it (LN1);
a cycle turns by it (C1); and a recordable cycle has it whole (C3).
Quantisation is the whole-number condition on what the cut lost.

## What is not shown

- C1–C4 are QC1-T4 and RMG1-T2 read through LN1's identity. No new
  number follows; what is new is that the counted quantity has been
  identified with a quantity defined by a cut.
- "Lost" here is the part lost by the *root* of the element, in a cut
  *across* its axis — the largest any cut can lose. The lost part of H
  itself is not the return (tested).
- C5 restates GR2's turn. Whole values of that turn occur at
  r = 4/3, 9/8, 16/15 r_s; no physical meaning is claimed for those
  radii.
- Why nature's records must be silent in content 1 is QC1's premise,
  not derived.

## Claim boundary

```text
seen(H) = 1 + 2·lost(M)                                              PROVED
RETURN AROUND A CIRCLE = −2π × LOST PART                             PROVED (QC1-T4 with LN1)
CONNECTION = LOST PART × dφ ; ITS CHANGE IS RMG1's CURVATURE          PROVED
SILENT ⇔ q × LOST PART WHOLE ; HEIGHT = k × SEEN PART                 PROVED, agrees with QC1's certified ladder
GR2's TURN = 4π × LOST PART AT THE LOWER TOWER LEVEL                  PROVED
THE PREMISE OF RECORDABILITY                                         TAKEN from QC1
A NEW PREDICTION                                                     NONE
```

## Reproduce

```text
python lc1_lost_is_counted.py
python -m unittest test_lc1_exact
```
