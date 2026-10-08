# QP1 — Three turns of a free circle on the exact frame, and the rational point 1/3 against measured rates

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy) for the rates; mpmath for the numbers.

RC1 states the owner's closure rule in the line's formulas for a centre at rest. This stage takes it to the
turning centre of SW2 and to measured rates.

Sources read before building: **SW2** (the exact frame; memory Re(r_s/R)), **SW1** (swirl r_s·a/r³), **MO1**
(free readings take the largest count; M5), **CL2**, **RC1**, **QC1-T3** (silence when a whole number of turns).
Measured rates and the published formulas: **arXiv:1408.0884** (eqs. 1–3 and Table 1, read at source), which
quotes the three rates of GRO J1655−40, the two of XTE J1550−564 and H 1743−322, and the masses from the
companion stars.

## 1. Three rates from the line's form

A free circle of radius r in the plane of the equator of SW2's frame (largest count; units GM/c² = r_s/2 = 1,
a in the same units):

```text
round rate        Ω = 1 / ( r^(3/2) + a )                                             Q1
in-out rate       κ² = Ω² ( 1 − 6/r + 8a/r^(3/2) − 3a²/r² )                           Q3
up-down rate      ν² = Ω² ( 1 − 4a/r^(3/2) + 3a²/r² )                                 Q2
```

They are the formulas printed in the cited paper. **Q4:** at a = 0 the second is MO1-M5 and the third equals the
first. **Q5:** to first order in a the plane of the circle turns at Ω − ν = r_s·a/r³, SW1's swirl.

So a free circle has three turns:

```text
the circuit itself                         Ω
the turn left over by the orbit's shape    Ω − κ          (RC1: zero left over ⇔ κ/Ω rational)
the turn of the orbit's plane              Ω − ν          (the swirl)
```

## 2. Measured

Three sharp rates are seen together in the X-rays of GRO J1655−40: 441 ± 2, 298 ± 4 and 17.3 ± 0.1 Hz. Reading
them as Ω, Ω − κ and Ω − ν (the reading of the cited paper, put in here):

```text
mass 5.300 ± 0.062 suns ,   a = 0.286 ± 0.003 ,   r = 5.680 ± 0.035         this stage
mass 5.31 ± 0.07 ,          a = 0.285 ± 0.003 ,   r = 5.68 ± 0.04           the paper
mass from the companion's motion   5.4 ± 0.3
```

**The ratio.** κ/Ω = 1 − 298/441 = 0.324 ± 0.010. The rational point 1/3 is 0.95 of that uncertainty away. In
H 1743−322 the pair 240 ± 3 and 165 (+9, −5) gives 0.31 (+0.03, −0.05).

## 3. The closure rule as a test

Take the rule literally: the pair is seen where the orbit closes, κ/Ω = 1/3 (three circuits per in-out period).
Then two measured rates fix the mass, with nothing else put in:

```text
GRO J1655−40     441 and 17.3 Hz      ⇒   mass 5.256 ± 0.030 ,  a = 0.289 ,  r = 5.712
                                          lower rate predicted 294.0 ± 1.3 ; measured 298 ± 4
                                          mass from the companion          5.4 ± 0.3
XTE J1550−564    183 and 13.08 Hz     ⇒   mass 8.84 ± 0.29 ,  a = 0.339 ,  r = 5.52
                                          mass from the companion          9.10 ± 0.61
```

Both masses agree with the masses measured from the companion stars, inside their uncertainties.

The ladder of RC1 for a centre at rest, for reference:

```text
κ/Ω     r / r_s     upper : lower     round rate × mass
1/2     4           2 : 1             1428 Hz·suns
1/3     3.375       3 : 2             1843
1/4     3.2         4 : 3             1996
2/5     3.571       5 : 3             1693
```

## What is put in

- SW2's frame and "free readings take the largest count".
- The reading of the three measured rates as Ω, Ω − κ, Ω − ν. It is the cited paper's; the line does not derive
  what makes the X-rays vary at those rates.
- The closure rule itself (RC1): tested here, not derived.
- G·M of the Sun and c. Uncertainties by first-order propagation.

## What is not shown

- That pairs of such rates sit near 3 : 2 is known (general knowledge). The rates and the three-rate fit are
  the cited paper's. What this stage adds: the rates come out of the line's own frame; 3 : 2 is the first-ladder
  point 1/3; and with closure two rates give a mass, which matches on the two sources that have one. Whether that
  last fit has been published is not checked.
- Why 1/3 and not 1/2 or 1/4: not derived.
- In two other observations of GRO J1655−40 (451 with 18.3 Hz, 446 with 18.1 Hz) the lower rate was not detected;
  with the mass and a above, the slow rate puts those circles at κ/Ω = 0.309 and 0.314. This is in the direction of
  the rule (no pair away from the rational point) but three observations are not a test.
- Spins of the same sources from other methods are recalled to be higher for GRO J1655−40 (not re-read here). If
  those are right, the reading of the rates is wrong for it.

## Refutation

```text
R8   a pair seen together whose κ/Ω is away from every low rational by several of its uncertainties
R9   a source where closure + two rates gives a mass outside the mass from its companion
```

## Claim boundary

```text
Ω, κ, ν ON THE EXACT FRAME FROM LARGEST COUNT                                PROVED (symbolic); equal to the published formulas
PLANE TURNS AT THE SWIRL r_s a/r³ AT FIRST ORDER ; a = 0 IS MO1              PROVED (symbolic)
THREE-RATE FIT OF GRO J1655−40                                               REPRODUCED (5.300 against 5.31)
MEASURED PAIR AT THE RATIONAL POINT 1/3                                      0.324 ± 0.010 : consistent (0.95 σ)
CLOSURE + TWO RATES ⇒ MASS                                                   5.256 ± 0.030 against 5.4 ± 0.3 ; 8.84 ± 0.29 against 9.10 ± 0.61
THE CLOSURE RULE                                                             NOT DERIVED ; survives this test
THE READING OF THE THREE RATES                                               PUT IN (cited paper)
```

## Reproduce

```text
python qp1_three_turns_and_the_rational_point.py
python -m unittest test_qp1
```

## Later note (FK1, 9 October)

The three rates satisfy ½[(in-out)² + (up-down)²]/round² = seen/(1 − aΩ)². For GRO J1655−40 the left side is
0.5141 ± 0.0031 from the rates alone and seen = 0.4931 ± 0.0030. Nothing above is changed.
