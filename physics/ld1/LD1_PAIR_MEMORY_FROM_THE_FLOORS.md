# LD1 — The memory between two neutral centres, from the floors of their modes: 1/r⁶, the number ¾, and a rule for unlike pairs

Monty Dabas. 9 October 2026. Python 3.12. Exact algebra (sympy) and numbers.

DL1 ended on what the line lacks to say anything about a fluid: how two atoms act on each other. This stage gets
the far part of it from what the line already has.

Sources read before building: **WQ1-W3** and **LC1-C4** (a mode — a pair of readings turning at rate ω — has the
floor ½κω: what stays seen when nothing is lost), **ZP1 / ZP2** (floors between walls; explicit tails), **DO1**
(the block [[1, k], [k, 1]]: every cut reads alike, lost part k²), **MC1** (least-cost field r^−(d−2)), **CV1-T3**
(the pattern (2, −1, −1)/r³), **DC1** (memory between two masses), **DL1**, RKF **theorum/28** §4. Tabulated
values: **arXiv:0902.3929**, Tables III and D, read at source.

## Results

**F1 — two modes read together.** Two equal modes with coupling block [[1, k], [k, 1]] have own rates ω√(1 ± k).
Their floors add to

```text
κω · √( (1 + √F) / 2 ) ,        F = 1 − k²   ( = R − D of the block on the diagonal ) ,
```

instead of κω. The lost floor is κω·(k²/8 + 5k⁴/128 + …), and for every coupling it lies between κω·k²/8 and
κω·(1 − 1/√2)·k². Reading two modes together always lowers the floor: two neutral centres attract.

**F2 — the coupling through the least-cost field.** The second derivatives of r^−(d−2) between the two centres
have rates (d − 1, −1, …, −1) × (d − 2)/r^d: no trace, because the field is least-cost. In three cuts (2, −1, −1)/r³,
CV1-T3's pattern. The sum of the squares is d(d − 1).

**F3 — the law.** With one mode per cut on each centre,

```text
lost floor = κω · ( k₀² / r^(2d) ) · d(d − 1)/8 ,        three cuts:   (3/4) · κω · k₀² / r⁶ .
```

The sixth power is twice the number of cuts; ¾ is 3·2/8.

**F4, F5 — unlike centres.** For rates ω_A, ω_B the floor sum is √(ω_A² + ω_B² + 2√(ω_A²ω_B² − c²)) and the lost
floor at second order is κc²/4ω_Aω_B(ω_A + ω_B). With α the static response of each centre, the constant of the
pair is C_AB = (3/2)·κ·α_Aα_B·ω_Aω_B/(ω_A + ω_B), and the rates and κ drop out of

```text
C_AB = 2·C_AA·C_BB / ( C_AA·α_B/α_A + C_BB·α_A/α_B ) .
```

## Numbers

The rule of F5 against the tabulated constants of the ten unlike pairs of noble atoms (atomic units; the table
states 1% uncertainty):

```text
pair      rule     tabulated   deviation        pair      rule     tabulated   deviation
He–Ne     3.043    3.03        +0.4 %           Ne–Kr     27.28    27.3        −0.1 %
He–Ar     9.522    9.55        −0.3 %           Ne–Xe     39.23    39.7        −1.2 %
He–Kr     13.35    13.42       −0.5 %           Ar–Kr     91.24    91.1        +0.2 %
He–Xe     19.32    19.6        −1.4 %           Ar–Xe     134.1    134.5       −0.3 %
Ne–Ar     19.55    19.5        +0.2 %           Kr–Xe     192.1    192         +0.05 %
```

Eight of ten within half a percent; the two pairs of xenon with a light atom are 1.2 – 1.4% low.

## What this gives and what it does not

```text
the far action between two atoms      derived: attraction, 1/r⁶, the number ¾, the rule for unlike pairs
the constant of one kind of atom      not derived: it needs the atom's response α and its rate ω
a fluid (DL1's arch, T_c, the diagonal)   not reached: it also needs the size of the atom, where the attraction stops.
                                          With a size σ the count form gives a = (2π/3)·C/σ³ per pair; the line has no σ
```

## What is put in

- One mode per cut on each centre, each a pair with floor ½κω; the coupling through the least-cost field.
- For the numbers: the tabulated like-pair constants and static responses.

## What is not shown

- The attraction of neutral atoms from coupled floors, the ¾ and the rule for unlike pairs are known (general
  knowledge; no source re-read for the derivation). The line's part: the floor is its "seen with nothing lost",
  the coupled floor is a function of F = R − D of the diagonal block alone, and the pattern and the sixth power
  come from the least-cost field and the number of cuts.
- A real atom is many modes. Read as one, its rate comes out 1.1 – 1.5 times its ionisation energy (ionisation
  energies recalled). The rule for unlike pairs does not need the rate, which is why it does so well.
- Far apart the finite speed of light changes the sixth power to a seventh; not treated.

## Claim boundary

```text
FLOOR OF TWO COUPLED MODES = κω √((1 + √F)/2) ; BOUNDS FOR EVERY COUPLING                    PROVED (exact)
LEAST-COST COUPLING: RATES (d−1, −1, …, −1), NO TRACE, SQUARES d(d−1)                        PROVED (d = 3, 4, 5)
THREE CUTS: (3/4) κω k₀²/r⁶                                                                  PROVED
UNLIKE PAIRS: C_AB = 2C_AAC_BB/(C_AAα_B/α_A + C_BBα_A/α_B)                                    PROVED
THE RULE ON TEN PAIRS OF NOBLE ATOMS                                                         WITHIN 1.5 % (eight within 0.5 %)
THE SIZE OF AN ATOM ; A FLUID's NUMBERS                                                      NOT REACHED
```

## Reproduce

```text
python ld1_pair_memory_from_the_floors.py
python -m unittest test_ld1
```
