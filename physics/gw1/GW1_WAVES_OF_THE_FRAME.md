# GW1 — Waves of the frame's order defect: a pair, content two, the same count

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy).

WQ1 left one row empty: the waves of the gravity field. This stage
builds them from the law of TP1 and asks whether WQ1's argument applies.

Sources read before building: **TP1** (the law ¼I₁ + ½I₂ − I₃ on the
order defect of the frame, selected by equivalence), **WQ1** (a wave
mode is a pair; recordable ⇔ area = 2πκn ⇔ E = nκω), **QC1-T3**
(content q silent iff qΘ ∈ 2πℤ; the covariance has content 2),
**R43.3 / DC1-K1** (ledgers that exchange share one unit), **EM1**,
**TL1**, **EG1** (energy passes between the sides).

## The frame

Flat, plus a small transverse change travelling along z:

```text
e₁ = (1 − a/2) ∂_x − (b/2) ∂_y ,     e₂ = −(b/2) ∂_x + (1 + a/2) ∂_y ,     e₀ = ∂_t ,   e₃ = ∂_z ,
```

with a(t, z), b(t, z) of first order.

## Results

**V1 — the law to second order.** Pointwise,

```text
Q₂ = ½ [ (∂_t a)² − (∂_z a)² + (∂_t b)² − (∂_z b)² ] .
```

A difference of two squares — a rate part and a gradient part — as
E² − B² is for light. The rate term is positive.

**V2 — two waves.** Stationarity gives ∂_t²a − ∂_z²a = 0 and the same
for b: two waves, uncoupled, at the cone speed.

**V3 — content two.** Turning the frame about the direction of travel
by θ turns (a, b) by 2θ. The pattern returns after half a turn and
changes sign after a quarter. Light has content one; this has the
content of the covariance in QC1-T3.

**V4 — a mode is a pair.** For a = A(t) cos kz, averaged over a period,

```text
energy = ¼ [ (dA/dt)² + k² A² ] = (k/2)( w² + p² ) ,     w = A √(k/2) ,  p = (dA/dt)/√(2k) ,
```

an oscillator of rate ω = k: a pair in the sense of WQ1.

**V5 — the count.** By WQ1 a recordable mode has E = nκω. These waves
exchange energy with the other sides (TL1, EG1); by R43.3, ledgers that
exchange and keep the sum share one unit. So the unit is the same κ:

```text
E = n h ν ,   with the h of light.
```

**Control.** A turn of the frame that differs from place to place is
not a wave: the law gives no equation for it (TP1's equivalence, at
second order).

## Numbers

```text
a wave of 100 Hz and strain 10⁻²¹       energy per count 6.6·10⁻³² J
                                         flux 1.6·10⁻³ W/m²   →   2.4·10²⁸ counts per m² per second
```

(the flux formula is the standard one, used for illustration only).

## Reading

The quantum of the gravity wave is not put in. It follows from three
certified pieces: the law of the frame (TP1), the recordability of a
pair (QC1 → WQ1), and the shared unit of ledgers that exchange (R43.3).
The same theorem gives the photon and this; the only difference the
framework sees is the content — one for light, two here.

## What this is not

- Not a new prediction. That gravity waves carry hν per count is what
  physics assumes; single counts are far beyond detection (10²⁸ per m²
  per second in a detected wave).
- Second order only, on a flat frame, one direction of travel. Waves on
  a curved field, and the interaction of counts, are not built.
- V5 uses R43.3's hypothesis: exchange that keeps the summed readout
  for every preparation. Nothing weaker was shown to suffice.
- The static field r_s/r carries no count by this argument.

## Claim boundary

```text
Q₂ = ½[(∂_t a)² − (∂_z a)² + (∂_t b)² − (∂_z b)²], POINTWISE              PROVED (symbolic)
TWO UNCOUPLED WAVES AT THE CONE SPEED                                    PROVED
CONTENT TWO                                                              PROVED
ONE MODE = ONE PAIR, ω = k                                               PROVED
E = n κ ω FOR A RECORDABLE MODE                                          FOLLOWS from WQ1
SAME UNIT AS LIGHT                                                       FOLLOWS from R43.3 under its hypothesis
LOCAL TURN OF THE FRAME IS NOT A WAVE                                    PROVED at second order
BEYOND SECOND ORDER; CURVED BACKGROUND; DETECTION OF A SINGLE COUNT      NOT BUILT / not within reach
```

## Reproduce

```text
pip install sympy
python gw1_waves_of_the_frame.py
python -m unittest test_gw1
```
