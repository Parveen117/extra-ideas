# GB1 — Where gravity sits in the block, and where it does not

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic over the
Gaussian rationals.

OB1 ended with a guess: light is the law of the traceless part of the
block; is gravity the law of the scalar / volume part that light leaves
unused? This stage checks the guess. **It is wrong**, with exact
witnesses, and the check shows where gravity does sit.

Sources read before building: **OB1** U1, U6; **IN1** T4–T5; **MO1**
M1, M6; **GR1** T2; **GR2** G1; **RMG1** (positive responses as points of
the Klein disc); **FR1-T4** (the conjugation and ι).

## 1. The three factors of a block

**B1.** Every invertible block is (scale) × (phase turn) × (unit block),
with det = scale² · phase². On a reading ρ → MρM†:

```text
phase turn     Exp(ιχ)         does nothing at all: a reading is blind to it
scale          a real factor   multiplies n and r; changes n² − r·r        (×3 gives n × 9, invariant × 81)
unit block     det = 1         changes n and r; keeps n² − r·r
```

So the scalar / volume part of the block is not free for gravity. Its
volume direction is the phase turn — the direction along which the
content q of OB1 is counted, and which no reading can see. Its scale
direction changes the one thing every frame agrees on.

## 2. A scalar gravity bends light too little

**B2.**

```text
clock factor only, flat space      u″ + u = r_s E² / ( 2 L² (1 − r_s u)² )      deflection   r_s / b      Sun: 0.876″
a common factor on time and space  light-like directions unchanged              deflection   0            Sun: 0
the form of MO1                    u″ + u = (3/2) r_s u²                        deflection   2 r_s / b    Sun: 1.751″
```

The measured value is the last one (general knowledge). A gravity
carried by the scalar part alone — either as a clock factor or as a
common scale — is refused by light.

## 3. Where it does sit

**B3.** The field of MO1 at a point, written as a block, is

```text
H = ( 1 + β r̂·C ) / N ,        N² = 1 − β² ,   β² = r_s / r ,
```

a self-dagger unit block: uncut and cuts, no turn, determinant one. It
carries the resting unit reading to (1/N ; β r̂/N) — the density and
proper density of GR1, now with a direction. It is a *response* in the
sense of RMG1: a point β·r̂ of the velocity disc at each place.

**B4 (algebra against group).**

```text
light      acts through the algebra:   q ( Fρ + ρF† )     proportional to the content q; opposite on the other sheet
gravity    acts through the group:     g ρ g†             the same on every reading; invariant kept
```

## 4. The assignment

```text
part of the block            what it is                                  status
phase turn  Exp(ιχ)          the direction of the content q (charge)     readings blind to it; its field is light (OB1)
unit block  det = 1          frames: boosts and turns                    gravity's field is a field of these (B3); static mass → boosts
scale                        changes the invariant                       unused — this is the open gate "what fixes a scale"
```

The guess is replaced by this: light and gravity are both fields of
frame changes of readings. Light's is weighted by the content and lives
in the algebra; gravity's is unweighted and lives in the group. The part
left over is not gravity but scale — and that is exactly the part
nothing in the line has been able to fix (QC5, HB1, EM1).

## 5. Claim boundary

```text
SCALE × PHASE × UNIT BLOCK ; THEIR ACTION ON A READING                         PROVED
CLOCK-FACTOR-ONLY DEFLECTION r_s/b ; COMMON-FACTOR DEFLECTION 0                PROVED (first order / exactly)
"GRAVITY IS THE SCALAR / VOLUME PART OF THE BLOCK"                             REFUSED
MO1's FIELD IS A SELF-DAGGER UNIT BLOCK ; IT REPRODUCES GR1's DENSITIES        PROVED
LIGHT THROUGH THE ALGEBRA × q ; GRAVITY THROUGH THE GROUP, UNWEIGHTED          PROVED as actions on a reading
A BLOCK LAW FOR THE GRAVITY FIELD (what D does for light)                      NOT FOUND — still the thesis of GR1/MC1/MO1 with its assumptions
ROTATING SOURCES (the turn part of the unit block)                             NOT TOUCHED
WHAT FIXES THE SCALE                                                           OPEN
```

## 6. Reproduce

```text
python gb1_gravity_in_the_block.py
python -m unittest test_gb1_exact
```
