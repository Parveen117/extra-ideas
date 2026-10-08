# MO1 — Motion in the memory field: the largest-count histories are the known orbits

Monty Dabas. 7 October 2026. Python 3.12. Exact rational arithmetic for
the theorems; floats only in the illustration. First item of "Next" in
the gravity thesis.

Question: a reading moves through the memory field of GR1. CL2 says the
unturned history has the largest count. Does the largest count give the
orbit?

Sources read before building: **CL2-T1/T2** (count of a leg = g·τ; the
straight history is the largest), **GR2-G1** (at each place one frame
change makes the static reading pure), **GR1-T2** (N² = 1 − m), **MC1-T4**
(m = r_s/r), **PR2-T1** (the form).

## 1. The one assumption

By GR2-G1, at radius r the locally pure frame is the one boosted inward
with β² = m = r_s/r. **Assumed here:** those frames share one flat space
and one time t. Then the count of a leg is g·dτ with

```text
dτ² = dt² − ( dr + β dt )² − r² dφ² .
```

Everything below follows from this form and from "free readings take the
largest count" (CL2, applied leg by leg in the locally pure frames).

## 2. Results

**M1 (two ways of writing the count).** With static time
dT = dt − β·dr/N²,

```text
dτ² = N² dT² − dr²/N² − r² dφ² ,        N² = 1 − r_s/r .
```

**M2 (the two special readings).** A reading at rest counts N·dt — GR1's
clock factor, recovered from the boost picture. A reading moving with the
frame counts dt: it is the pure reading.

**M3 (radial law).** With the conserved E and L of the form,

```text
ṙ² = E² − N² ( 1 + L²/r² ) .
```

**M4 (circles).**

```text
(dφ/dt)² = r_s / (2r³)                      the period law, exactly
count factor² = 1 − (3/2)·r_s/r             for a reading on a circle
light circles at r = (3/2)·r_s              where that factor vanishes
last stable circle at r = 3·r_s , L² = 3·r_s²
```

**M5 (orbit equation).** With u = 1/r:

```text
u″ + u = r_s/(2L²) + (3/2)·r_s·u² ,          near a circle:  (radial / angular rate)² = 1 − 3·r_s/r .
```

The last term is what the memory field adds to the closed ellipse; the
orbit turns by 2π(1/√(1 − 3r_s/r) − 1) ≈ 3π·r_s/r per revolution.

**M6 (light).** For a light-like history the first-order solution gives a
deflection 2·r_s/b. A clock factor alone, with no stretching of the
radial direction, would give half of that.

Illustration (standard constants):

```text
advance of Mercury's orbit           42.98″ per century
deflection of light at the Sun       1.751″
count factor on the Earth's orbit    1 − 1.5·10⁻⁸
```

## 3. Reading

- The form of §1 is, exactly, the one whose largest-count histories are
  the known orbits around a point mass. So the thesis, with the one
  assumption of §1, reproduces the known motion to all orders in r_s/r.
  It predicts nothing new here.
- All the content is in the assumption: that the locally pure frames fit
  together into one flat space with one time. GR2-G1 proves that such a
  frame exists at each place; it does not prove that they fit together.
- The factor of two for light (M6) shows the assumption is doing real
  work: it supplies the stretching of the radial direction that the
  clock factor alone does not.

## 4. Claim boundary

```text
RIVER AND STATIC FORMS OF THE COUNT AGREE; N RECOVERED; PURE READING COUNTS dt       PROVED
RADIAL LAW; CIRCLES; LAST STABLE CIRCLE; ORBIT EQUATION; NEAR-CIRCULAR ADVANCE       PROVED from the form
LIGHT DEFLECTION 2 r_s / b                                                           PROVED to first order
"THE LOCALLY PURE FRAMES SHARE ONE FLAT SPACE AND ONE TIME"                          ASSUMED — not derived
"FREE READINGS TAKE THE LARGEST COUNT" IN A VARYING FIELD                            ASSUMED — CL2 proves it in one frame
ANYTHING THE KNOWN THEORY DOES NOT ALREADY GIVE                                      NOTHING
WAVES; MORE THAN ONE SOURCE; ROTATING SOURCES                                        NOT TOUCHED
```

## 5. Reproduce

```text
python mo1_motion.py
python -m unittest test_mo1_exact
```

## Later note (OR1, 8 October)

"Free readings take the largest count in a varying field" is listed above as assumed. OR1 shows that a record
shared by neighbouring histories (GE2-T2) keeps the one whose count is stationary, in any field: the assumption is
the same rule as QC1's premise of recordability. The rule itself remains a premise. Nothing above is changed.
