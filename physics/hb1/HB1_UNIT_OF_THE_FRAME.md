# HB1 — The unit of a frame: count against response, and the minimum as one minus what is lost

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic over the
Gaussian rationals.

Owner's intuition for this stage: everything comes down to the minimum.
After the cut, the cosine part — what is its least value? That depends
on the frame, and on what the frame first measured. Numbers are tools
for measuring a frame; the scaling comes from its first fundamental
measure — time, and in this framework a clock can be anything, with no
external clock. And Planck's constant will come out of the frame's
invariant, 1 or √2, less the lost information.

Sources read before building: **R41.3** ([N, P] = ιR, R = C − (Q/2)·wrap),
**R41.4** (σ_N²σ_P² − Cov² ≥ ¼⟨R⟩², coefficient sharp), **R38** (c_Σ = 1/√2;
eight source events per block), **CL1-T2/T3** (mass² = clock curvature;
speed = 1 − curvature/2; [N, O] = E), **QC3-T1** (E² − O² = I),
**DM1-T1** (the seen leg is the cosine).

## 1. The first measure of a frame is its tick

A frame's clock is a tick T on marks, with count N. Its even part
C = (T + T†)/2 and its response P = ι(T − T†)/2 are all a frame needs to
have a unit.

**T1.** [N, P] = ι·C on readings that touch no wrap (R41.3; with the wrap,
[N, P] = ι·R, certified for Q = 3, 4, 5).

**T2.** var(N)·var(P) − Cov(N,P)² ≥ ¼⟨C⟩² (R41.4). R41's sharp instance is
reproduced exactly (variances 1 and ¼, covariance 0, R = −1); binomial
readings come within 15/14 of the bound.

So the least product of the spread of the count and the spread of the
response is **half the mean even part of the tick**. Nothing outside the
frame enters.

## 2. The minimum is one minus what is lost

**T3.** For a reading a_r = f_r·ζ^r with a real profile f and tick phase
ζ = c + ιs,

```text
⟨C⟩ = c · w_f ,        ⟨P⟩ = s · w_f ,        w_f = Σ f_r f_{r+1} / Σ f_r²  ≤ 1 ,
```

so the minimum is

```text
½ · c · w_f  =  ½ · ( 1 − loss to the tick's turn ) · ( 1 − loss to the finite profile ) ,
```

with c = 1 − (clock curvature)/2 by CL1. The frame's invariant is 1; the
unit it actually has is 1 less what it loses — to its own tick and to the
finiteness of what it holds.

```text
tick (c, s)        profile          overlap w_f      unit c·w_f
(3/5, 4/5)         box of 4         3/4              9/20
(3/5, 4/5)         box of 12        11/12            11/20
(4/5, 3/5)         binomial 8       8/9              32/45
(1, 0)             box of 12        11/12            11/12
(0, 1)             any              —                0
```

## 3. The ladder of clocks

```text
marks Q      even part of the tick      unit of the frame       CL1: speed, clock curvature
   2             −1                        (reversed)             −1 ,  4
   4              0                         0                      0 ,  2
   8             1/√2                      1/√2                   1/√2 , 2 − √2
   ∞              1                         1                      1 ,  0
```

- With the quarter-turn clock the unit is zero: count and response
  commute, and there is no minimum at all.
- With the eight-mark clock — R28's coin, R38's speed — the unit is 1/√2.
- A light-like frame has the full unit 1.

In this model **the unit of a frame, its cone speed, and one minus half
its clock curvature are the same number**: the even part of its tick.

## 4. Reading

- The minimum is not a constant of nature in this account; it is the
  even part of the frame's own first measure. Change the clock and the
  minimum changes.
- "The frame's invariant less the lost information" is exact here:
  1·(1 − (1 − c))·(1 − (1 − w_f)).
- More mass means a smaller unit: a frame whose clock is more curved
  has less room between count and response. At the quarter turn there is
  none.

## 5. Certificate

T1 on three open readings and on cyclic marks Q = 3, 4, 5 with the wrap.
T2 on five rational readings and R41's sharp instance. T3 on five tick
phases and four profiles: ⟨C⟩, ⟨P⟩, the overlap, the bound. Five tests;
dropping the wrap and a non-conjugating pairing are rejected.

## 6. Claim boundary

```text
[N, P] = ι C ; MINIMUM = ½|⟨C⟩|                                          PROVED (R41.3, R41.4 restated and re-certified)
⟨C⟩ = c·w_f ; MINIMUM = ½ (1 − turn loss)(1 − profile loss)              PROVED
UNIT OF THE FRAME = CONE SPEED = 1 − CURVATURE/2                         PROVED in the tick model
PLANCK'S CONSTANT, IN ANY PHYSICAL UNIT                                  NOT DERIVED — the unit here is a pure number
WHICH CLOCK NATURE USES; WHY ONE UNIT WOULD BE SHARED BY ALL FRAMES      OPEN
THAT A HEAVIER PARTICLE HAS A SMALLER PHYSICAL MINIMUM                   NOT CLAIMED — it is a statement about ticks, not yet about matter
```

A commutator whose right side is a cosine, on a discrete set of marks,
is known (general knowledge). The stage's content is R41's result read
as the frame's unit, and its identity with the speed and the clock
curvature of CL1.

## 7. Reproduce

```text
python hb1_unit_of_the_frame.py
python -m unittest test_hb1_exact
```
