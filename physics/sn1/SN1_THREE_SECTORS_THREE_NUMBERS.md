# SN1 — The spread of a count through a cut, in the three sectors: 1/4, 1/2, 1/3

Monty Dabas. 9 October 2026. Python 3.12. sympy.

FD1-F4: an exclusive record seen through a cut counts like independent coins with biases q, the readings of the
seen kernel. Its mean is Σq and its spread Σq(1 − q) — QD1's F = R − D, the share law. The pure number is
spread / mean. This stage takes the seen share in each of the three sectors of the algebra, with the sector's
angle spread evenly, and computes it.

Sources read before building: **FD1-F4**, **QD1-D3** (spread n − n²), **DG1-D3/D4** (share law; memory = tanh²),
**IN1** (seen cos², lost sin²), **GR1** (N² + m = 1), **RMG2, RMG9** (the three sectors: R² = −1, the line between,
K² = +1), **DO1** (the observer is not fixed).

## Results

In one notation, lost / seen = t², with t = tan φ (turn), t = x (shear), t = sinh η (boost).

```text
sector              seen             angle spread evenly on     spread / mean       next number
turn   R² = −1      cos² φ           the quarter circle         1/4                 0
shear  ε² = 0       1/(1 + x²)       0 … X ,  X large           1/2                 1/4
boost  K² = +1      1/cosh² η        0 … s ,  s large           1/3                 1/15
```

**S1 — turn.** The seen share follows the arcsine law, mean ½; spread / mean = ¼ exactly. An observer anywhere on
the circle of cuts, none preferred.

**S2 — shear.** spread / mean = ½[1 − X/((1 + X²)·arctan X)], rising to ½.

**S3 — boost.** spread / mean = tanh²(s)/3, rising to ⅓: 0.310 at s = 2, 0.3333 at s = 5. Here seen = 1 − tanh² and
lost = tanh² is DG1's split.

**S4 — one cut.** A single cut on the diagonal has ½; along the reading 0; seeing nothing, 1 (independent marks).

## Against measurement

The count of charges passing a conductor is such a count (general knowledge); spread / mean is its measured
noise ratio.

```text
turn     1/4     open cavities                 measured 1/4 "in agreement with theory" (arXiv:cond-mat/0009087)
boost    1/3     long disordered wires         measured; 1/3 reached for large reservoirs (arXiv:cond-mat/9808042)
shear    1/2     symmetric double barrier      measured 1/2 (recalled, not re-read at source)
```

## What is put in

- That the sector's angle is spread evenly up to its largest value. In the known theory of these conductors this
  is a result; here it is taken. For the turn it is the statement that no cut is preferred.
- Which conductor belongs to which sector is read from the known forms of their seen shares.

## What is not shown

- 1/4, 1/2, 1/3 (and 0, 1/4, 1/15) are the known numbers. The line's part: they are the share law F = R − D
  averaged over the three sectors of one algebra, one number per sector.
- A derivation of the even spread from the line (for the boost: from the adding of angles when cells are
  joined, DU1-U3) is not done.

## Claim boundary

```text
SPREAD / MEAN = 1/4 (TURN), → 1/2 (SHEAR), tanh²(s)/3 → 1/3 (BOOST) ; NEXT NUMBERS 0, 1/4, 1/15     PROVED
MEASURED 1/4, 1/3, 1/2                                                                             MATCH ; KNOWN NUMBERS
EVEN SPREAD OF THE ANGLE                                                                           PUT IN
```

## Reproduce

```text
python sn1_three_sectors_three_numbers.py
python -m unittest test_sn1
```
