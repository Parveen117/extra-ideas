# Synthesis — one block, one conservation law, two fields and a unit

Monty Dabas. 7 October 2026. Second edition, after FR1–DC1.
Certified core: every stage folder here; `python run_all_tests.py` runs
161 exact tests in 28 files (Python 3.12, standard library only).

## In one paragraph

Everything in the physics line is a statement about one object — the
block certified in EMK-1, read on the framework's cut-complex carrier —
and one law of conservation: what a cut observes and what it loses keep
the uncut. A reading is a part of the block; a frame is a choice of
cuts; geometry belongs to the frame. The block has exactly three
factors. Two of them can be made to vary from place to place, and they
are the two long-range fields: light, weighted by a counted content, and
gravity, the same for every reading. The third cannot, and it is the
unit in which everything is counted.

## 1. The object

```text
block      ( a + ι b ) + ( r + ι s )·C          a uncut,  r three cuts,  s three turns,  ι = C₁C₂C₃        OB1, FR1
channels   det = Δ∥ + Δ⊥                                                                                  EMK-1 T2
split      E = scalar part,  O = traceless part,  E² − O² = I on the unit quadric                         SY1, QC3
sectors    circular (turns), dual (light-like), split (boosts) — by the sign of O²                        SY1
law        D = ∂_t + ΣC_i∂_i ,  D·D̄ = wave operator                                                       OB1, PR1
```

A frame over the rationals carries two anticommuting cuts; on the
cut-complex carrier three, and no fourth. The number of cuts is the
number of readings after which nothing of a single reading is lost.
[FR1, IN1, DM1]

## 2. The conservation law

```text
reading tensor       ρ = ½ ( n + r·C )
single reading       n² = r₁² + r₂² + r₃²                 lost to one cut = observed by the others
every frame          n² − r·r is the same                 a form with one plus and three minus signs
record of two        n² − r·r = 2p(1−p)(n₁n₂ − r₁·r₂)     what no cut of the frame holds
```

Recoverable memory is frame-dependent and sums to the uncut;
unrecoverable memory is the invariant. [IN1] What a frame with fewer cuts
cannot read it meets as mass: the two-cut rest mass is the third
reading, the massless three-cut law on a mode of third wave number k₃ is
the two-cut law with flip rate −c·k₃, and the two sheets — mass and
antimass — are the two directions along the unseen cut. [DM1, PR2]

## 3. The two parities

```text
one arrow U and its reverse:      E = (U + U⁻¹)/2          O = (U − U⁻¹)/2          E² − O² = I ,  memory = −O²
                                  second order, real       first order, circular
                                  records, heat, loss      phase, clock, mass
```

Deterministic refinement keeps the odd part; symmetric records keep the
even part and turn the odd part into memory. [QC3, GE2] The same split
gives: the share law 4p(1−p) as the variance of a two-way choice [PH3,
QC4]; curvature of shared flat readings = −variance [QC4]; mass² = clock
curvature and speed = 1 − curvature/2 [CL1, PR1]; the unit of a frame as
the even part of its tick [HB1]; the character χ_N ≤ N for combined
masses and ≥ N for turned histories [MS1, CL2].

## 4. The three factors of a block

```text
factor        made local                         physics                               on a reading
phase         field F = E + ιB ; D(F·C) = J      light; content q is the charge        nothing, except through F
unit block    field H = (1 + β r̂·C)/N            gravity (thesis)                      n, r move; invariant kept
scale         refused as a field                 the unit of count                     everything rescales
```

- **Light.** Potential, field, source and force are parts of the block
  under the one law D: the four field equations are its four grades; the
  source is a reading, hence no magnetic source and a conserved charge;
  the force is the frame change Fρ + ρF†. u² − S·S = ¼|F·F|². [OB1, EM1]
- **Gravity.** Source: unrecoverable memory. Field: recoverable memory
  m = r_s/r, the least total variance between neighbouring shells in a
  frame of three cuts. Clock factor N = √(1 − r_s/r); horizon = balanced
  share. Removable at a point by one frame change; with the locally pure
  frames sharing one flat space, the largest-count histories are the
  known orbits. [GR1, MC1, GR2, MO1, GB1, GRAVITY_THESIS]
- **Unit.** Frames that exchange must count in one unit; its value is
  not set by any law of the block. [R43.3, DC1, SC1]

Light acts through the algebra and is odd between the sheets; gravity
acts through the group and is even. [GB1, FR1]

## 5. What comes out with numbers

```text
return angle of real fluids on near-critical cycles        0.9 – 1.1 rad, unit-free                 PH1
the same angle as the rotation left by a strain cycle      CO₂: 12 elements 41.45°, limit 53.06°    PH2
static acceleration at the Earth's surface                 9.82 m/s²                                GR1
advance of Mercury's orbit                                 42.98″ per century                       MO1
deflection of light at the Sun                             1.751″  (scalar-only gravity: 0.876″)    MO1, GB1
count lost on a turned history, legs (3,5), (−3,5)         8 against 10                             CL2
memory between two 10⁻¹⁴ kg masses, 2.5 s                  0.012                                    DC1
```

## 6. Where the account can be wrong

```text
R1  a static clock near a mass differing from √(1 − r_s/r) at second order
R2  matter and antimatter falling, or curving clocks, differently
R3  a source of gravity that is not a record: a single light-like reading must not gravitate, a pair must
R5  the once-around turn not containing π r_s/r at first order
R6  an isolated mass losing coherence exponentially at a rate set by its own gravity
R7  two masses coupled only by gravity that never share memory
```

## 7. Every stage in its place

| Stage | Statement | Place |
|---|---|---|
| PH1–PH3 | return angle of a fluid; strain-cycle rotation; equal share of two flat readings | §3: share law, curvature = −variance |
| QC1, QC2 | integer content; p ≙ κ∂/∂w in the ensemble | §3: counting; the unit |
| QC3, QC4 | even/odd parities; return is memory; compass law | §3 |
| QC5 | moving share between flat gauge readings; logistic law | §3 on a field of blocks |
| PR1–PR4 | sector speed and mass; two sheets; accelerated frame; the source τ | §2, §3 |
| CL1, CL2 | mass² = clock curvature; turned histories | §3 |
| MS1 | speed² + memory = 1; combined masses | §2, §3 |
| SY1 | the unit block and its split | §1 |
| FR1, IN1, DM1 | cuts and dimension; the invariant; the unseen leg | §1, §2 |
| HB1 | the unit of a frame | §3, §4 |
| EM1, OB1 | the wave under a cut; one block, one law | §4: light |
| GR1, MC1, GR2, MO1 | gravity as memory; 1/r from least cost; thesis core; motion | §4: gravity |
| GB1, SC1 | the three factors; scale is the unit | §4 |
| DC1 | one unit for all frames; memory between masses | §4, §6 |

## 8. Limits

Nothing here changes a classical equation or predicts a number that
existing theory does not give; the account reproduces known physics from
one certified algebra and takes sides where the known theory leaves a
choice (§6). The value of every dimensional constant is outside it, by
SC1. The gravity thesis rests on stated assumptions (the dictionary, the
least-cost reading, one flat space for the locally pure frames) and is
static or test-particle only. Fluid numbers use reference equations of
state through CoolProp. No stage has been checked by a proof assistant
or by external review.

## 9. Open

```text
the two pure numbers        the size of one unit of content for light (≈ 1/137); of gravity for a given mass
the assumption of MO1       derive or refuse that the locally pure frames fit together
many ledgers                when memory between systems becomes effectively permanent (R44)
a block law for gravity     what D is for light
an experiment               the programmed chain of PH2; the pair memory of DC1
```

## 10. Reproduce

```text
python run_all_tests.py
```
