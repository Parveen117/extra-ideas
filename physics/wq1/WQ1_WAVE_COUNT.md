# WQ1 — A wave is a pair of readings; its area comes in whole units

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic; floats only in the illustration.

The owner's statement: every curvature is made of waves, and in
discrete form it gives hν — like the energy of a photon.

Sources read before building: **QC1-T1/T2** (the two readings are a
canonical pair; the return is their bracket), **QC1-T3** (content q is
silent iff qΘ ∈ 2πℤ), **QC1-T4** (ladder k + n, k = q/2), **QC2**
(p acts as κ ∂/∂w), **EM1** (energy of a wave train divided by its
frequency is the same in every frame; the wave gives the place of the
unit, not its counting), **DC1-K1** (one unit for frames that
exchange).

EM1 left the counting open. QC1 supplies it.

## Results

**W1 — a wave mode is one pair.** A mode is a pair of readings (w, p)
turning at rate ω, with energy E = (ω/2)(w² + p²). Its orbit encloses an
area A with

```text
A / 2π = E / ω .
```

**W2 — the area comes in whole units.** The return around the orbit is
the bracket of the pair (QC1-T2): Θ = A/κ. A mode that can be recorded
— content 1 silent (QC1-T3) — needs

```text
A = 2πκ · n ,        E = n κ ω = n h ν ,        h = 2πκ ,  ν = ω/2π .
```

A fractional area is not silent. Content 2, the covariance, is already
silent at half counts.

**W3 — the ladder.** With raising z and lowering κ ∂/∂z (bracket κ),

```text
H = κω ( z ∂/∂z + ½ ) ,        H zⁿ = κω ( n + ½ ) zⁿ .
```

Equal steps κω = hν; the lowest level is half a step, the k = q/2 of
QC1-T4.

**W4 — the count is the same in every frame.** A change of frame scales
the energy and the frequency of a wave by the same factor (EM1), so n
does not change.

**W5 — many waves.** For a superposition each mode's area is counted
separately and the energies add.

## Reading

The area enclosed by a pair of readings is the return of QC1 — the
flux of the curvature form through the loop. Stated for a wave: the
curvature flux of each mode is a whole number of units, and its energy
is that number times hν. That is the owner's sentence. The unit is the
common unit of DC1; its size is not fixed by anything here (SC1).

What this applies to: any wave whose two readings form a pair. For the
light field the pair is its two parts under a cut (EM1) — the photon.
For a fluid the pair is the two readings of PH3/QC1 — the quantum of
sound. For the gravity field the pair would be the two parts of a wave
of the frame's order defect; those waves have not been built in this
line, so that case is not claimed.

## What this is not

- It is the old quantum condition for an oscillator, reached by the
  framework's route: recordability (silence of content 1) is what makes
  the area whole. Nothing new is predicted for light or sound.
- "Every curvature" is proved here only for a mode whose readings are a
  flat pair. A static field such as r_s/r is not a wave and carries no
  count by this argument.
- The value of h is not derived.

## Claim boundary

```text
A/2π = E/ω ; RETURN = A/κ                                            PROVED
RECORDABLE (CONTENT 1 SILENT) ⇔ A = 2πκ n ⇔ E = n κ ω                PROVED
LADDER κω(n + ½), EQUAL STEPS                                        PROVED exact on polynomials
COUNT INDEPENDENT OF FRAME                                           PROVED from EM1's scaling
SUPERPOSITION: COUNTS SEPARATE, ENERGIES ADD                         PROVED
WAVES OF THE GRAVITY FIELD AND THEIR COUNT                           NOT BUILT
THE VALUE OF THE UNIT                                                NOT DERIVED (SC1)
```

## Numbers

```text
green light, 500 nm      3.97·10⁻¹⁹ J per count
sound, 1 kHz             6.6·10⁻³¹ J per count
```

## Reproduce

```text
python wq1_wave_count.py
python -m unittest test_wq1_exact
```
