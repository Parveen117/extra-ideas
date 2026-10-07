# EM1 — The electromagnetic wave under a cut: observed, lost, and what stays

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic over the
Gaussian rationals.

Owner's question: can Planck's constant not be derived from
electromagnetic waves — after the cut, the wave turning into lost and
observable?

Sources read before building: **IN1** (observed² + unrecoverable =
uncut², frame-independent), **HB1** (the unit of a frame is the even part
of its tick; 1 for a light-like frame), **QC2** (E = count × rate needs a
unit for the count), **PR1-T4** (the frame change of a sector along its
ray), **CL1/CL2** (time is a count; counts do not depend on the frame),
**NSR1** (the electromagnetic response in the framework).

## 1. The wave as one cut-complex reading

Join the two fields with the framework's own unit: F = E + ιB. The cut
is conjugation: E is one sheet, ιB the other.

**T1.** With u = ½(E·E + B·B) the energy density and S = E × B the flow,

```text
u² − S·S = ¼ |F·F|² ,        F·F = (E·E − B·B) + 2ι E·B .
```

This is IN1 for light: uncut u, observed flow S, and a part no frame can
turn into flow.

**T2.** A single wave has F·F = 0: u = |S|, nothing lost. Two waves
running against each other have u² − S·S = 2u₁u₂(1 − cos) — the share law
of IN1-T5. For amplitudes 2 and 3: u = 13, flow 5, unrecoverable 12.

**T3.** A change of frame along x with rapidity η turns F about x through
the cut-complex angle ιη. That turn is exactly the usual change of E and
B; it keeps F·F and changes u and S.

**T4 (what the cut reads).** E·E and B·B separately depend on the frame;
E�E − B·B and E·B do not. A reading that is purely electric in one frame
has a magnetic sheet in another (E² 1 → 25/16, B² 0 → 9/16, difference 1
in both).

So after the cut the wave does divide into observed and lost, the
division depends on the frame, and what every frame agrees on is F·F.

## 2. Energy over frequency

**T5.** For a wave along the frame change with factor D: amplitude × D,
energy density × D², frequency × D, and the length of a train of a fixed
number of crests × 1/D. Hence

```text
( energy of the train ) / ( frequency )   is the same in every frame.
```

The number of crests is a count, and counts do not depend on the frame
(CL1). The wave therefore carries one frame-independent quantity with the
dimension of an action. This much comes from the wave alone.

## 3. What the wave gives of Planck's constant, and what it does not

```text
that energy/frequency of a wave is frame-independent                  FROM THE WAVE (T5)
that the wave splits into observed and lost with an invariant rest    FROM THE WAVE (T1–T4)
that the least step of a count against its response is ½ × 1 for a
  light-like frame                                                    FROM THE FRAME (HB1)
that energy/frequency comes in whole steps                            FROM COUNTING (QC1, QC2) — not from the wave
the size of one step in joules × seconds                              NOT DERIVABLE: it is the unit in which a count is
                                                                      converted to energy × time
```

On the last line: in the present system of units the value of h is fixed
by definition, like the value of c; it converts frequency to energy the
way c converts time to length. A number of that kind cannot come out of
any theory. What a theory can give is a pure number. For light meeting
matter that number is the coupling of the wave to a charged reading —
about 1/137 — and in this language it would be a share: how much of the
wave a charged frame's cut observes. R41 already lists it among the
identification targets. Nothing here derives it.

## 4. Reading

- The owner's picture holds for the wave: the cut makes an observed part
  and a lost part, the two keep the uncut, and the frame decides the
  division.
- The wave supplies the *place* of Planck's constant — an invariant
  action per train — and the frame supplies its *minimum*. The wave does
  not supply the fact that the action is counted.
- The light-like frame has the full unit (HB1). If every frame reads
  action against light, the unit is common to all of them; that would
  answer why one constant serves every frame. Stated, not derived.

## 5. Certificate

T1 on five rational field pairs; T2 with parallel and counter-running
circular waves; T3–T4 on four fields and three boosts, against the usual
component formulas; T5 at three Doppler factors. Five tests; a real turn
in place of the cut-complex one is rejected.

## 6. Claim boundary

```text
u² − S·S = ¼|F·F|² ; SINGLE WAVE NULL ; TWO-WAVE SHARE LAW                    PROVED
FRAME CHANGE = TURN THROUGH ιη ; INVARIANTS E² − B², E·B                      PROVED
TRAIN ENERGY / FREQUENCY FRAME-INDEPENDENT                                    PROVED (crest count taken as a count)
PLANCK'S CONSTANT DERIVED FROM THE WAVE                                       NO — its place and its minimum, not its counting or its value
THE COUPLING OF LIGHT TO CHARGE (≈ 1/137)                                     NOT DERIVED; named as the pure number to aim at
```

The identity of T1 and the complex form of the fields are classical
(general knowledge). What the framework adds is that they are IN1 — the
same conservation of the uncut — with the electric and magnetic fields as
the two sheets of one cut.

## 7. Reproduce

```text
python em1_wave_cut.py
python -m unittest test_em1_exact
```
