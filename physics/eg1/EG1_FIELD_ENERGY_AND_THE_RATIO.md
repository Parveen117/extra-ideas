# EG1 — The light field's energy makes gravity as heat does; the ratio law with a source

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic for
the block identities; sympy for the curvature part.

The owner's statement: what the thermal side did, the electromagnetic
side will do — each symmetry has its diagram, and gravity connects them
all. This stage tests it for the light field.

Sources read before building: **EM1**, **OB1** (field F = E + ιB;
reading-type products FρF†), **IN1**, **CV1** (contracted curvature of
the fall frame), **NC1** (native ratio ρ), **TL1** (units follow the
clock), **QC1** (content q is a count), **NT-1** (response element as a
Hessian), and the correction to R3 in the gravity thesis.

## Results

**G1 — the field read through the cuts.** For a radial field F = E n·C,

```text
½ F F†      = u                          u = E²/2
½ F C_a F†  = u ( 2 n_a n·C − C_a )      pull along the field, push across it
```

the pattern (u ; −u, +u, +u). Its invariant part is zero.

**G2 — it curves the frame.** The memory

```text
m = r_s / r − q₂ / r²
```

has contracted curvature (q₂/r⁴) × (1 ; −1, 1, 1): exactly the pattern
of G1 for a field falling as 1/r². A reading whose invariant part is
zero gravitates. The correction to R3 — the source is the whole
reading, not its invariant — is confirmed by the framework's own
curvature.

**G3 — the ratio law with a source.** For the fall-frame family the
energy component of the law is (r m)′/r² = k·u for any energy density
u, and the native ratio of NC1 obeys

```text
ρ + ½  =  (r m)′ / 2m  =  k u r² / 2m .
```

ρ = −½ is the empty case. The excess of the ratio over −½ **is** the
energy density, of whatever kind, in units of m/r². For the charged
memory, ρ + ½ = q₂ / 2(r r_s − q₂).

**G4 — potential follows the clock.** Content is a count (QC1), the
same at every place; so the potential, energy per content, obeys
V_A N_A = V_B N_B (TL1-L2).

**G5 — the light field's response element.** The Hessian of
u(D, B) = D²/2ε + B²/2μ is the master element with

```text
determinant = 1/(εμ) = (speed)² ,       eigenvalue ratio = μ/ε = (impedance)² .
```

## Numbers

```text
one electron charge     length √q₂ = 1.38·10⁻³⁶ m      (its r_s: 1.35·10⁻⁵⁷ m)
one coulomb             length √q₂ = 8.6·10⁻¹⁸ m
```

## Reading

Heat entered gravity as energy, and gravity entered heat through the
clock (TL1). The light field does the same: its energy and stresses
curve the frame (G2), and its potential follows the clock (G4). In the
native ratio the two are not distinguished at all — G3 has one u.
Gravity connects the sides through energy, as stated; no side needs a
law of its own to gravitate.

G3 also sharpens NC1: the open question "why ρ = −½" becomes "why does
the energy component of the law equal the energy density". Empty space
then has ρ = −½ for the same reason matter has more.

## What is not shown

- G2 and G3 are the known law of a charged mass, reached through the
  framework's curvature; nothing new is predicted.
- "k·u" on the right of G3 is the law of CV1/TP1 with a source; the
  constant k is not derived.
- G5 is a dictionary entry. How a gravity field moves the light field's
  response element is not built here.
- Other symmetries (the internal ones of the Yang–Mills line) are not
  treated.

## Claim boundary

```text
FIELD READINGS: ENERGY u, PATTERN (u; −u, u, u), INVARIANT PART ZERO            PROVED exact
CHARGED MEMORY ⇒ CONTRACTED CURVATURE (q₂/r⁴)(1; −1, 1, 1)                       PROVED (symbolic)
ρ + ½ = (r m)′/2m ; ENERGY COMPONENT (r m)′/r²                                   PROVED
POTENTIAL FOLLOWS THE CLOCK                                                      PROVED from TL1
RESPONSE ELEMENT OF THE LIGHT FIELD: SPEED², IMPEDANCE²                          PROVED (dictionary)
THE CONSTANT k ; WHY THE ENERGY COMPONENT EQUALS THE ENERGY DENSITY              NOT DERIVED
A PREDICTION BEYOND THE KNOWN LAW                                                NONE
```

## Reproduce

```text
pip install sympy
python eg1_field_energy_and_the_ratio.py
python -m unittest test_eg1
```
