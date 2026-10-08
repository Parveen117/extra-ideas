# The gravity thesis

Monty Dabas. 7 October 2026. Certified core: `gr2/` (exact rational
arithmetic, six tests) on top of `gr1/`, `mc1/`, `fr1/`, `in1/`, `ms1/`.

## The thesis

> Gravity is the memory of lost information that cannot be recovered.
> An observer does not read change as such; change appears only through
> what is lost. The lost and the observed together keep the uncut
> conserved. Geometry belongs to the frame that cuts, not to what is cut.

## The statement

Terms. The **uncut** is the count n of a reading tensor. A **cut** is an
involution of the frame; its **reading** is r. **Recoverable memory** is
what one cut loses and a complementary cut holds. **Unrecoverable
memory** is what no cut of the frame holds.

**1. Conservation.** In every frame, (observed by the three cuts)² +
(unrecoverable) = (uncut)², and the unrecoverable part is the same in
every frame. *[IN1-T2, T4, T5 — proved]*

**2. The source is unrecoverable memory.** It is zero for a single
reading and equals 2p(1−p)(n₁n₂ − r₁·r₂) for a record of two. No frame
change removes it. Read as energy–momentum it is the rest mass squared.
*[IN1-T5, GR2-G2 — proved; the reading as rest mass is PR2-T2]*

**3. The field is recoverable memory.** A static reading near a source
has clock factor N equal to the imbalance of its two light-like readings,
and memory m = 1 − N² against the radial cut. The potential is
Φ = −½c²·m. *[GR1-T1, T2 — proved in the model]*

**4. Equivalence.** That memory is removed at any one place by a single
change of frame: the static reading becomes a pure reading. What cannot
be removed is how the required frame change varies from place to place;
carried once around the source, a frame returns turned.
*[GR2-G1, G4 — proved in the model]*

**5. The law.** The memory is spread so that the total variance between
neighbouring readings is least. A frame on the cut-complex carrier has
three cuts, and with three the least-cost memory is r_s/r, with a flux
that is the same through every surface and adds when sources add.
*[MC1-T1–T4, FR1-T3, GR2-G3 — proved on scale-uniform shells; that the
cost is quadratic between neighbours, and that memory is what is spread,
are assumed]*

**6. No memory, no change.** ∂N = −2g·N·√m: where nothing is lost,
nothing varies, whatever the flip rate. *[GR1-T1 — proved in the model]*

**7. Both sheets alike.** Mass and antimass offer the same source, see
the same clock, and lose the same count on the same histories.
*[PR4-T4, PR3-T5, CL2-T4 — proved in the model]*

## What follows, with numbers

```text
clock factor              N = √(1 − r_s/r)                                  GR1
static acceleration       c² r_s / (2 r² N)            Earth: 9.82 m/s²     GR1
horizon                   balanced share p = ½, never reached at finite flip GR1
potential                 Φ = −½ c² × memory           Earth surface: memory 1.39·10⁻⁹
once around the source    a direction turns by 2π(1/N − 1) ≈ π r_s/r         GR2-G4
```

The once-around turn for a circuit of the Earth at radius 7020 km is
3.97·10⁻⁹ rad = 0.82 milliarcsecond. To first order this is the size of
the known space-curvature part of the precession of an orbiting
gyroscope; the full measured effect is half as large again, and the
remaining part is not in this model (general knowledge; not checked
against the literature here). At second order the model's 2π(1/N − 1)
exceeds the cone value 2π(1 − N) by exactly 2π(1 − N)²/N.

## What is proved, assumed, and put in

```text
PROVED (exact, in the framework's algebra)
  conservation in every frame; single reading null; record invariant            IN1
  three cuts, no fourth; mass with three cuts only through the conjugate sheet   FR1
  least-cost profile r^{−(d−2)}; flux constant and additive                      MC1, GR2-G3
  ∂N = −2gN√m; N² + m = 1; horizon = balanced share                              GR1
  removal of the static memory by one frame change; survival of the source       GR2-G1, G2
  the once-around turn                                                           GR2-G4, RMG1-T3

ASSUMED
  the static observer is the static flow reading with a light-like inflow        (the dictionary)
  the cost is quadratic and between neighbours                                   (motivated by QC4)
  what is spread at least cost is the memory, not the clock factor               (MC1-T5: the alternative fails at second order)
  the frame's readings are cut-complex                                           (the framework's carrier)

PUT IN
  the constant relating flux to source: 2G/c² = 1.485·10⁻²⁷ m/kg
  spherical symmetry; a static situation
```

## What would refute it

A statement with no refutation condition is a definition.

```text
R1  a static clock near a mass measured to differ from √(1 − r_s/r) at second order in r_s/r
R2  matter and antimatter shown to fall, or to curve clocks, differently
R3  a source of gravity with no rest mass that this account cannot express as a record
    (two non-parallel light-like readings do have unrecoverable memory; a single one has none)
R4  a certified derivation, inside the framework, that the least-cost quantity is not the memory
R5  the once-around turn measured and found not to contain π r_s/r at first order
R6  an isolated mass in superposition losing coherence exponentially at a rate set by its own gravity   (DC1)
R7  two masses coupled only by gravity that never share memory                                          (DC1)
```

R1 and R2 agree with what is known at present.

**Correction to R3 (same day).** Statement 2 makes the source the
unrecoverable memory alone, so a single light-like reading would not
gravitate. That cannot stand together with MO1: light is deflected by a
mass (2 r_s / b), so its momentum changes, and unless momentum is not
kept the mass must be pulled by the light in return. A source that is
only the invariant is therefore inconsistent with the thesis's own
motion result. The source has to be the whole reading (n; r), of which
the invariant is the rest part. Statement 2 is to be restated
accordingly; until then R3 is withdrawn as a prediction and recorded as
an error of this document.

## What is not here

```text
the value of G, or any relation between G and c alone         refused (dimensional; ledger Theorem F)
motion of bodies; bending of light; orbits                    MO1, on one assumption
anything time-dependent; waves                                not built
field equations beyond one static function                    not claimed (ledger D3′)
the remaining third of the gyroscope precession               not in the model
a unit: ħ, k_B, the Planck scale                              open (QC2, QC5)
```

## Motion (added with MO1)

With one further assumption — the locally pure frames of statement 4
share one flat space and one time — the count of a moving reading is
dτ² = dt² − (dr + β·dt)² − r²dφ², β² = r_s/r, and its largest-count
histories are the known orbits: the period law, the last stable circle
at 3·r_s, the advance of a near-circular orbit by 3π·r_s/r per turn
(Mercury: 42.98″ per century), and the deflection of light 2·r_s/b
(1.751″ at the Sun). *[MO1 — proved from the form; the assumption is not
derived, and nothing new is predicted]*

## Place in the block (added with OB1, GB1)

The static field of this thesis, as a block, is H = (1 + β r̂·C)/N: a
self-dagger unit block at each place, acting on every reading alike by
g ρ g†. Light acts through the algebra, weighted by the content q.
Gravity as the scalar or volume part of the block is refused: light
would bend by half the measured amount, or not at all. *[GB1 — proved]*

## Next

1. Derive, or refuse, the assumption of MO1: do the locally pure frames
   fit together into one flat space with one time?
2. R3: the record of two light-like readings as a source, exactly.
3. The second-order term of the once-around turn, now with the full form
   of MO1 instead of the static frames of GR2-G4.

## Reproduce

```text
python gr2/gr2_thesis_core.py
python -m unittest discover -s gr2 -p 'test_*.py'
```

## Curvature (CV1)

Statement 4 removes the memory at a point; it does not remove the
curvature of the frame field, which is the order defect of the connected
frame directions (EMK-C1 form). For the frame of MO1 its contraction
vanishes exactly when (r m)′ = 0, i.e. m = r_s/r, and the time–time part
is minus half the variance operator of MC1: the least-variance law of
statement 5 is the vanishing of the time–time contracted curvature, and
the remaining components remove the constant that least variance alone
allows. What every frame agrees on is the uncontracted curvature,
r_s/r³ × (1, −½, −½). That the contraction, and not something else, is
the law is taken here, not derived.

## The law, in the frame's own connection (TP1)

In the connection in which the frame field is constant, every loop
closes (information invariance holds exactly) and the field is the order
defect of the frame. Among the three-coefficient quadratic laws in that
order defect, the requirement that any frame may be used at each place
(GR2-G1 made exact) selects one combination, 1 : 2 : −4, and that
combination equals the curvature law of CV1 up to a boundary term. So
the open line of CV1 — why the contraction and not the whole curvature —
is answered: the whole curvature is zero in the frame's own connection;
the contraction is the same law in the torsion-free one. General
relativity is the member of the family that equivalence selects; no
prediction beyond it follows.

## The native curvature of the field (NC1)

With the framework's own curvature (NT-3: F = −¼[X_i,X_j], X = G⁻¹δG) the
field has two invariant sizes and one pure ratio between them,
ρ = rψ′/sinh ψ. The field of statements 3–5 is exactly ρ = −½ at every
radius. Making the total native curvature-square stationary instead
gives the same field at first order and a different one at second order,
with Mercury's advance 57.3″ or 32.2″ per century in place of 42.98″:
that law is refused by observation. Why the ratio is −½ is not derived.

## Sources of any kind (TL1, EG1)

Heat enters as energy and is met by the clock factor: κ·N is constant in
equilibrium (TL1). The light field enters the same way: a radial field
reads (u; −u, u, u) through the cuts, its invariant part is zero, and
the memory r_s/r − q₂/r² has exactly that pattern as contracted
curvature (EG1). This confirms the correction to R3: the source is the
whole reading. In the native ratio of NC1 the kinds are not
distinguished: ρ + ½ = k u r²/2m.

## Waves and their count (WQ1, GW1)

To second order on a flat frame the law of TP1 is a difference of two
squares; it gives two uncoupled waves at the cone speed with content
two. Each mode is a pair of readings, so a recordable mode has
E = nκω (WQ1), and because these waves exchange energy with the other
sides the unit is the one they share (R43.3): the h of light. The count
is a consequence here, not an assumption; it is far below detection.

## MO1's assumption (MA1)

Item 1 of "Next" is answered for the static field of one centre. With the radial stretch of the fall frame
left free, the law of CV1/TP1 gives A′ = 0 in empty space, and frames at rest far away give A = 1: flat slices,
one time, and β² = r_s/r from the same computation. The assumption is the statement N·S = 1 and is not
independent of the law. *[MA1 — proved for the static, spherically symmetric case; the general case is not treated]*
