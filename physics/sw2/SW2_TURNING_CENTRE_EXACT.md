# SW2 — The turning centre exactly: the centre at rest, displaced by ι·a along its axis

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), exact; CV1's formulas for the connection and
the curvature of a frame.

SW1 found the swirl a/r³ at first order and showed that flat slices are refused at second order, with an exact
remainder. This stage finds the frame that the law accepts at every order.

Sources read before building: **SW1** (first-order swirl; remainder 9a² sin²θ/2r⁶), **MA1** (flat slices and one
time for the centre at rest), **MO1** (the form; a reading that moves with the frame counts dt), **CV1** (frame
curvature; with a centre at rest it depends on β only through the memory m, and linearly), **MC1** (the memory is
spread at least cost; in three dimensions r_s/r), **GR1** (N² = 1 − m; static acceleration), **OB1-U1** (ι = C₁C₂C₃;
a block has an uncut part, three cuts, three turns ι·C and the volume ι), **FR1** (three cuts on the cut-complex
carrier), **EG1**. Checked that no stage of the line treats this (dock, synthesis, both folders).

## The idea

On the cut-complex carrier a place is a block x·C. OB1-U1 says a block also has a turn part ι·s·C. Give the centre
one: move it from the origin to ι·a along its axis. Its distance from a place (X, Y, Z) is then

```text
R² = X² + Y² + ( Z − ι a )²      ⇒      R = r − ι·a·cos θ          (F1)
```

in the layout of flat space X + ιY = √(r² + a²)·sin θ·e^{ιφ}, Z = r cos θ. The memory of the centre at rest, r_s/r,
read at this distance is one cut-complex function

```text
Φ = r_s / R ,      Re Φ = r_s·r / ρ² ,      Im Φ = r_s·a·cos θ / ρ² ,      ρ² = r² + a² cos² θ = R·R̄ .
```

**F2.** Φ is still spread at least cost: it has no flat-space Laplacian away from the ring R = 0.
**F4.** Its multipoles are r_s·(ι a)^l: every one is fixed by r_s and a.

## The frame

```text
e⁰ = dt ,    e¹ = (ρ/s)·dr + b·ν ,    e² = ρ dθ ,    e³ = s·sin θ dφ ,
ν = dt − a sin² θ dφ ,    s² = r² + a² ,    b² = P(r)/ρ²   with P free.
```

With b = 0 this is flat space (G0). With a = 0 it is MO1's frame with memory P/r². The frame falls with speed b,
not against dt as in MO1 but against the turned form ν.

## Results

**G1 — one time, free fall.** For every P the frame direction e₀ is in free fall and e⁰ = dt: a reading that moves
with the frame counts dt (MO1-M2), exactly. Of MO1's assumption, "one time" survives the turning; "flat slices"
does not (SW1-T3).

**G2 — the contracted curvature is linear in the memory.** As in CV1, it depends on b only through P = b²ρ², and
linearly:

```text
(θ, θ):   ( r P′ − P ) / ρ⁴            (radial, radial):   ( 2P − 2rP′ + ρ² P″ ) / 2ρ⁴
```

and the others likewise. **G4:** at a = 0 these are CV1's components.

**G3 — the law.** All components vanish exactly when r·P′ = P, i.e. P = r_s·r:

```text
memory  =  b²  =  Re ( r_s / R ) .
```

The memory of the turning centre is the real part of the memory of the centre at rest, read at the displaced
distance. The law holds at every order in a; the remainder of SW1 is gone.

**W — the wider family.** With the fall function free in both r and θ (P = P(r, θ); long computation, stored in
`sw2_wide_family.json`): radial + (θ, θ) components give P_rr = 0; a combination of the (time, radial) and
(radial, round) components gives P·P_θ = 0. So the law forces P_θ = 0, and then G3. Within "flat layout, fall
against ν" the frame is unique.

**H — a held reading.** A reading held at fixed (r, θ, φ) has clock factor N² = 1 − Re Φ. Its acceleration a and
the turn Ω between its frame and carried directions are, exactly,

```text
a  =  ∇ Re Φ / 2N² ,        Ω  =  ∇ Im Φ / 2N² ,        so        a + ι·Ω  =  ∇Φ / 2N² .
```

(All four frame components, lower index, signature of IN1; orientation (r, θ, φ). The sense of Ω is that of SW1-T5.)

GR1's static acceleration and SW1's turn are the real part and the ι-part of one gradient. In OB1's words: the
rate of frame change of a held frame is one traceless block, cuts (boost) + turns, like the block of light.

**L1 — first order.** To first order in a, after relabelling the round angle (φ → φ + a·∫β/r² dr), the count form
is SW1's with swirl constant **r_s·a**.

**G5 — light-like shells.** The shell r is light-like where r² + a² = r_s·r. There is none when a > r_s/2.

## What this settles

```text
SW1's open item (the second-order frame)          found at every order                      G3
MO1's assumption for a turning centre             one time: exact; flat slices: replaced by "flat layout, fall against ν"    G0, G1
the swirl                                         not a second field: the ι-part of the same memory                         H
freedom of the turning centre                     two numbers, r_s and a; every multipole r_s (ιa)^l                        F4, W
```

## What is put in

- The family: the flat oblate layout, the turned form ν, a fall along e¹ only. The law then fixes the fall function
  (G3, W). That the law accepts no frame outside this family is not shown.
- The law of CV1/TP1 (zero contracted curvature in empty space).
- That a = (turning content)/(mass × c): through SW1's swirl constant, r_s·a·c = 2GJ/c², put in there.

## What is not shown

- This frame is the known exact field of a turning centre in its free-fall form, and the displaced distance is
  the known complex-shift construction of that field (general knowledge; no source re-read for this stage). What
  the stage adds is inside the line: the law of CV1 reaches it, the displacement is along the carrier's own ι, and
  the pair (acceleration, turn) is one block.
- Why a turning centre is a displaced one is not derived from a source. The ring R = 0 and the inside are not
  treated.
- H is for empty space. With a charge the turn is no longer a pure gradient (GM1).

## Claim boundary

```text
R = r − ι a cos θ ; Φ = r_s/R AT LEAST COST ; MULTIPOLES r_s (ιa)^l         PROVED (symbolic)
ONE TIME AND FREE FALL FOR EVERY FALL FUNCTION                             PROVED (symbolic)
CONTRACTED CURVATURE LINEAR IN THE MEMORY ; LAW ⇔ MEMORY = Re(r_s/R)       PROVED (symbolic, exact in a)
UNIQUENESS WITH THE FALL FUNCTION FREE IN r AND θ                          PROVED from the stored long computation
HELD READING:  a + ι Ω = ∇Φ / 2N²                                          PROVED (symbolic, exact)
FIRST ORDER = SW1 WITH SWIRL CONSTANT r_s a                                PROVED (symbolic)
NO FRAME OUTSIDE THE FAMILY                                                NOT SHOWN
WHY TURNING IS A DISPLACEMENT ALONG ι                                      NOT DERIVED
NEW MEASURED NUMBER                                                        NONE (see GM1 for the ratio that follows)
```

## Reproduce

```text
python sw2_turning_centre_exact.py          about four minutes
python -m unittest test_sw2
python sw2_wide_family.py                   about twenty minutes; rewrites sw2_wide_family.json
```
