# CC1 — The cost of carrying a content: ratios 9/4, 5/2, 4, 9/2, 6, 25/4, 7

Monty Dabas. 9 October 2026. Python 3.12. Exact integer linear algebra (numpy integers, sympy ranks).

MG1: the mass count itself is open, because the cell and the coupling are not fixed. In a ratio of two costs
both cancel. This stage takes the ratios the Yang–Mills line does give.

Sources read before building: **YM-22** (the time-face log of content ½ is exactly 3a/4), **YM-40** (content
ladder e^(−a·C_c); the excess of content 1 over ½ lands on 5/4), **YM95** (three colours: Laplacian of a loop
−(4/3)·loop), **YM96** (completeness form; first Casimir (N² − 1)/2N), **YM98** (Casimir numbers enumerated), **MG1**.

## Results

**C1 — three colours.** The Casimir form C = ½[Σ E_ij E_ji − ⅓(Σ E_ii)²] is built on the tensor powers
(first content)^p ⊗ (dual)^q for p + q ≤ 4 as integer matrices. On each it is diagonal; its top value is
(p² + q² + pq)/3 + p + q and the space of the top value has exactly the dimension of the content (p, q).
First content 4/3, adjoint 3 (YM95, YM96).

**C2 — two colours.** Contents ½, 1, 3/2, 2 have j(j + 1); the difference of the first two is YM-40's 5/4, their
ratio 8/3.

**C3 — the ratios.** Cost of a content over the cost of the first:

```text
content      8      6      15a    10     27     24      15s
ratio        9/4    5/2    4      9/2    6      25/4    7
```

**C4 — against a simulation.** The standard theory with three colours has been simulated with static sources in
these eight contents (arXiv:hep-lat/0006022, Tables X and XI: ratios of the potentials, extrapolated to the
continuum). Weighted means over its 19 separations below 1 fm:

```text
content        8        6        15a      10       27       24       15s
ratio          2.25     2.5      4        4.5      6        6.25     7
simulation     2.246    2.499    3.965    4.474    6.187    5.957    7.069
difference     −0.2 %   −0.0 %   −0.9 %   −0.6 %   +3.1 %   −4.7 %   +1.0 %
```

All inside the 5 % that the paper states; 122 of the 133 single points are within two of their own errors.
As printed, the columns 27 and 24 lie in the opposite order to their ratios; with the errors taken as
independent the 24 column is several errors low. The paper reports no significant violation. That pair is the
one place where the table does not sit on the ratios.

## What this is and is not

- The ratios are exact at the free end of the Yang–Mills line (cost per cell of a content = its Casimir number).
  The simulation is at the other end, the continuum. That the ratios survive between the two ends is what the
  simulation shows; the line has not proved it. For two colours YM-40 certifies the first step (5/4) along its
  own trajectory.
- The numbers of C3 are the known "Casimir scaling" (general knowledge). The line's part is C1: they are computed
  from its own completeness form, with exact ranks.
- A simulation is a computation of the standard theory, not a measurement in a laboratory.
- These are ratios of forces between held contents, not ratios of masses of free kinds of matter. Mass ratios
  remain open.

## Claim boundary

```text
CASIMIR NUMBERS ON TENSOR POWERS, p + q ≤ 4 ; TOP SPACE HAS THE DIMENSION OF THE CONTENT   PROVED (exact)
RATIOS 9/4, 5/2, 4, 9/2, 6, 25/4, 7                                                        PROVED ; KNOWN NUMBERS
SIMULATION MEANS WITHIN 5 %                                                                MATCH ; 27/24 order noted
SURVIVAL OF THE RATIOS TO THE CONTINUUM                                                    NOT DERIVED
```

## Reproduce

```text
python cc1_cost_of_a_content.py
python -m unittest test_cc1
```
