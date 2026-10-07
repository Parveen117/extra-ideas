# TP1 — All loops closed, gravity in the order defect: equivalence selects one law

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy —
outside the standard-library suite).

CV1 left one line open: why should the *contraction* of the curvature
vanish and not the whole curvature, when information invariance in the
paper's sense (every loop closed) asks for the whole? This stage closes
it, conditionally, and finds the place of general relativity among the
laws the framework allows.

Sources read before building: **EMK-C1** (curvature as order defect; the
frame term c_μν^ρ), the information-invariance paper **Theorem A**
(invariance ⇔ every loop closed ⇔ exact), **GR2-G1** (memory is removed
at a point by one frame change), **GB1** (the unit block acts alike on
every reading and keeps the invariant), **OB1** (a local phase gives the
light field), **CV1**, **MC1**.

## 1. Two connections on one frame field

A field of frames has its own connection: the one in which the frame is
constant from place to place.

**T1.** For the fall frame of the thesis, that connection has zero
curvature — every loop closes — and its torsion is minus the order
defect of the frame, [e_a, e_b] = c_ab^k e_k (for the fall frame:
β′, β/r).

So information invariance holds exactly, and the field has not
disappeared: it sits in the order defect of the frame. CV1 used the
other connection (no torsion), in which the same field appears as
curvature. They are two readings of one frame field.

## 2. The laws available

A law that does not depend on coordinates and is quadratic in the order
defect has three terms:

```text
Q = a₁ I₁ + a₂ I₂ + a₃ I₃ ,
I₁ = c_abk c^abk ,   I₂ = c_ab^k c_k^b_a ,   I₃ = c_b c^b  (c_b = c_ab^a) .
```

## 3. Equivalence selects the coefficients

GR2-G1: at each place a frame change removes the memory. Made into a
requirement: **the law must not depend on which frame is used at each
place.** Take flat space and read it in a frame that is boosted and
turned differently at each place; nothing is there, so Q must be a pure
boundary term.

**T2.** This holds if and only if

```text
a₁ : a₂ : a₃  =  1 : 2 : −4 ,        Q = k ( ¼ I₁ + ½ I₂ − I₃ ) .
```

Exact: with (1, 2, −4) every variation vanishes identically; the
conditions have rank two, so the combination is unique up to scale.

**T3.** For the fall frame, with that combination,

```text
¼ I₁ + ½ I₂ − I₃  =  −R  −  2 · div( c-trace ) ,
```

pointwise and exactly, R the twice-contracted curvature of CV1. The
selected law is therefore the law of CV1 — the known vacuum law — up to
a boundary term.

**T4.** Among all (a₁, a₂, a₃), the fall frame keeps the profile
m = r_s/r in the reduced variation only on the plane 2a₁ + a₂ + a₃ = 0.
The selected line lies in it.

## 4. Reading

```text
information invariance (all loops closed)   +   a field of frames      →  the field is the order defect
quadratic, coordinate-free law                                         →  three coefficients
equivalence: any frame at each place                                    →  one combination
that combination                                                        =  the curvature law of CV1, i.e. general relativity
```

General relativity is the member of a three-coefficient family that the
framework's own equivalence statement selects. The other members are
laws in which the frame at each place is itself physical; on the plane
2a₁ + a₂ + a₃ = 0 they still keep the profile r_s/r here.

This answers CV1's open line: the whole curvature *is* zero in the
frame's own connection, as information invariance asks; "contraction
zero" is the same law seen in the torsion-free connection.

## 5. What this is not

- Not a new prediction. Equivalence, as stated in GR2-G1, picks exactly
  the known law. A different prediction would need the framework to
  deny that the frame may be chosen freely at each place.
- Not derived from first statements alone: "quadratic" and "any frame at
  each place" are requirements stated here. The second rests on GR2-G1;
  the first is the lowest order available and is not derived.
- The family and its special combination are known in the literature on
  frame-based gravity. What is the framework's own is the route:
  information invariance supplies the flat connection, EMK-C1 the order
  defect, GR2 the equivalence.
- T3 and T4 are verified on the fall-frame family only; T2 on boosts
  along one axis and turns in one plane, separately and combined.

## 6. Claim boundary

```text
THE FRAME'S OWN CONNECTION IS FLAT; FIELD = ORDER DEFECT                    PROVED (fall frame)
LOCAL FRAME EQUIVALENCE ⇔ (1 : 2 : −4), UNIQUE UP TO SCALE                   PROVED on the stated class of local frame changes
SELECTED LAW = −R + BOUNDARY TERM                                           PROVED pointwise for the fall frame
PROFILE r_s/r KEPT ON THE PLANE 2a₁ + a₂ + a₃ = 0                            PROVED in the reduced variation
GENERAL FRAMES; SOURCES; THE SCALE k                                        NOT BUILT
WHY QUADRATIC                                                               NOT DERIVED
A PREDICTION DIFFERENT FROM GENERAL RELATIVITY                              NONE — would require a physical frame at each place
```

## 7. Reproduce

```text
pip install sympy
python tp1_frame_defect_law.py
python -m unittest test_tp1        (about three minutes)
```
