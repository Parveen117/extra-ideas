# TH1 — At the rung 1/3 the three turns of a free circle are one rate, twice it and three times it

Monty Dabas. 8 October 2026. Python 3.12. Exact arithmetic and sympy.

QP1 found the measured pair of X-ray rates at the rung 1/3 of RC1's ladder and left open why that rung. This
stage does not derive the choice. It shows what is special about that rung in the line's own equations, and
what else it predicts.

Sources read before building: **RC1** (rungs p/q), **QP1** (three turns; closure + two rates ⇒ mass), **MO1-M5**
(u″ + u = r_s/2L² + (3/2)·r_s·u²), **QC1-T3** (content q is silent iff q·Θ is a whole number of turns), **GW1**
(light has content one, the frame waves content two). Measured: **arXiv:1408.0884** Table 1 and
**arXiv:astro-ph/0202305** abstract, both read at source.

## Results

**H1 — whole numbers.** On the rung p/q the in-out rate, the left-over-turn rate and the circuit are
p : q − p : q times one rate Ω/q. The rung 1/2 is the only one on which two of them coincide (1 : 1 : 2), so it
would show as a line and its double. The lowest rung with three distinct turns is 1/3, and there they are

```text
in-out : left-over turn : circuit   =   1 : 2 : 3 .
```

**H2 — the term that leaves the turn over also doubles the in-out motion.** Near a circle, to second order in
the in-out amplitude A (ω² = 1 − 3x):

```text
ε(φ) = A cos ωφ + A²·( 3r_s/4ω² − (r_s/4ω²)·cos 2ωφ ) .
```

The last term of M5 — the one that makes the ratio √(1 − 3x) instead of 1 — gives the in-out motion a second
harmonic of relative size x·e/(4(1 − 3x)), e the eccentricity. At the rung 1/3 (x = 8/27) that is (2/3)·e.

**H3 — where the lines meet.** Read through the round angle, the lines of such a history are |n₁ + n₂ω|·Ω. The
second harmonic of the in-out motion falls on the left-over-turn line exactly when 2ω = 1 − ω, that is ω = 1/3;
there the third falls on the circuit. On this rung every line is a multiple of the in-out rate.

**H4 — content.** With Θ = 2π·ω the in-out phase left over per circuit, QC1-T3 gives: content 1 is silent on no
rung inside the ladder; content 2 on 1/2; content 3 on 1/3 and 2/3.

## Measured

Closure at 1/3 says that below the measured pair there is a third rate, the in-out rate, at half the lower one.

```text
XTE J1550−564     lower rate 183 ± 5 Hz     ⇒   in-out 91.5 ± 2.5 ,   circuit 274.5 ± 7.5
                  reported together: features near 92, 184 and 276 Hz                       astro-ph/0202305
GRO J1655−40      circuit 441 ± 2 Hz        ⇒   in-out 147.0 ± 0.7 ,  left-over 294.0 ± 1.3  (measured 298 ± 4)
                  the source paper describes its pair as multiples of one rate; a feature at the single rate
                  is not stated in its abstract
```

In XTE J1550−564 all three turns are there. In the line's reading the broad 92 Hz feature is the in-out rate.

## What is put in

- RC1, QP1 and their inputs; MO1-M5 for a centre at rest in H2 (the turning centre is not redone here).

## What is not shown

- That such pairs are multiples 2 and 3 of one rate is the observers' own description (cited abstract). What this
  stage adds is which three motions the 1, 2, 3 are, and that H2–H3 put a real second harmonic on the middle line.
- Why the rung 1/3 is taken is still not derived. H1 (lowest rung with three distinct turns) and H4 (content 3)
  are properties, not reasons. The line has a carrier for content 1 (light) and 2 (frame waves, GW1) and none
  yet for 3.
- No statement about how strongly each line shows in the X-rays.

## Refutation

```text
R10   a source with the pair at 2 : 3 and a third sharp rate that is not at 1 (the in-out rate) within its uncertainty
```

## Claim boundary

```text
ON THE RUNG p/q THE THREE TURNS ARE p : q−p : q ; 1/3 IS THE LOWEST WITH THREE DISTINCT          PROVED
SECOND HARMONIC OF THE IN-OUT MOTION FROM M5's LAST TERM ; SIZE x e/4(1−3x)                      PROVED (symbolic, second order, centre at rest)
IT MEETS THE LEFT-OVER-TURN LINE EXACTLY AT 1/3                                                  PROVED
CONTENT 3 SILENT ON 1/3 AND 2/3                                                                  PROVED (QC1-T3 applied)
XTE J1550−564: THREE FEATURES AT 1 : 2 : 3                                                       CONSISTENT (91.5 ± 2.5 against "near 92")
WHY THE RUNG 1/3                                                                                 NOT DERIVED
```

## Reproduce

```text
python th1_three_turns_as_harmonics.py
python -m unittest test_th1
```
