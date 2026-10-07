# DC1 — One unit for all frames, and what gravity does to a superposition

Monty Dabas. 7 October 2026. Python 3.12. Exact arithmetic over the
Gaussian rationals for the theorems; floats only in the illustration.

This stage reads two of the owner's own results with today's thesis:
**R43** in this repository (exchange between two ledgers, relative
calibration, pair memory) and obligation **D4** of the
information-invariance ledger in Publications (the gravitational
decoherence rate). Both turn out to answer questions the physics line
had left open.

Sources read before building: **R43.3** (exchange balance fixes relative
calibration, eq. 43.7), **R43.5/R43.6** (interaction phase cost and
distinguishability; reduced-pair memory deficit), **R43.7** (a common
action scale remains free), **R44.4** (ordering curvature invisible to
pair observations); the information-invariance paper v2 with Theorems
**H**, **N** and the D4 verdict; **IN1-T2/T5**, **GR1-T2**, **CL1-T2**,
**HB1**, **SC1**, and the gravity thesis.

## 1. Why one unit serves every frame (R43.3)

HB1 found that each frame has its own unit and asked why nature shows one.
The owner's R43.3 already answers it.

**K1.** For two ledgers with the same readout E and unit factors b_A, b_B,
and the exchange S,

```text
‖ [ S , b_A E⊗1 + b_B 1⊗E ] ‖² = 2d ( b_A − b_B )² [ Tr E² − (Tr E)²/d ] .
```

Exchange keeps the summed readout for every preparation exactly when
b_A = b_B. Two frames that exchange anything must count in the same unit.
The common value stays free (R43.7) — it is the scale factor of SC1.

So the unit is universal because frames interact, and arbitrary because
nothing but counting sets its size. Re-certified here for two and three
roles.

## 2. Memory between two ledgers

Two ledgers, each with two branches, shares p and q, and a turn Δ on one
branch pair (Δ = φ₁₁ − φ₁₂ − φ₂₁ + φ₂₂ in general).

**K2.**

```text
memory of A      M = 1 − Tr ρ_A²  =  ½ · [ 4p(1−p) ] · [ 4q(1−q) ] · sin²(Δ/2) ,
IN1's invariant of A's reading     n² − r·r  =  2M .
```

The memory is the product of the two share laws and half the record
defect of the turn. And the frame-independent "rest" part of A's reading
(IN1) is this memory: what A cannot recover is what B holds.

**K4.** If B has one branch only, or Δ = 0, the memory is zero at all
times. A ledger alone loses nothing (IN1-T2).

## 3. The turn that gravity supplies

**K3.** By GR1 and CL1, a reading of rest turn g₁ at distance d from a
mass m₂ counts g₁(1 − r_s₂/2d) per unit time. Between branches this is the
phase

```text
φ = G m₁ m₂ τ / ( ħ d ) ,
```

so the gravitational Δ of two split masses is G m₁m₂τ/ħ times
(1/d₁₁ − 1/d₁₂ − 1/d₂₁ + 1/d₂₂).

## 4. The stance this gives on gravitational decoherence

```text
one mass in two places, nothing else        no memory, at any time
two masses, each in two places              memory ½ sin²(Δ/2): grows as τ², returns to zero, never exceeds ½
```

This differs from a law in which a single mass loses coherence by itself
at the rate E_G/ħ. In the thesis the field of a mass is *recoverable*
memory (statement 3) and a single reading loses nothing (IN1-T2), so an
isolated superposition has no partner to hold what it would lose. The
number G m²/(ħR) is still the right scale — as the rate of the *turn*
between two masses, not of a decay.

Illustration with the example of the owner's paper (m = 10⁻¹⁷ kg,
R = 10⁻⁷ m; scale 6.33·10⁻⁴ s⁻¹, as in D4):

```text
time        pair memory ½ sin²(Δ/2)        deficit under an exponential law at the same rate
10 s        5.0·10⁻⁶                       6.3·10⁻³
100 s       5.0·10⁻⁴                       5.9·10⁻²
1000 s      4.8·10⁻²                       3.6·10⁻¹
```

Two masses of 10⁻¹⁴ kg, each split by 250 μm, centres 450 μm apart, for
2.5 s: Δ = −0.31 rad, memory 0.012.

**Relation to D4.** The ledger's verdict stands as written: the printed
rate was not a rate, its repair is the Diósi–Penrose expression, and no
departure had been exhibited. This stage exhibits one. It does not
rewrite D4; it adds what the ledger said was missing — "a regime of
departure".

## 5. What would refute it

```text
R6  an isolated mass in superposition losing coherence at a rate set by its own gravity, exponentially
R7  two masses coupled only by gravity that never share memory (no turn Δ between their branches)
```

Either result contradicts the thesis. At present neither has been seen;
the parameter-free form of the single-mass law is reported to be
excluded by a search for the radiation it would cause (general
knowledge, not checked here), which is in the direction of R6 not being
met.

## 6. Claim boundary

```text
EXCHANGE CONSERVATION ⇒ EQUAL UNITS ; COMMON VALUE FREE                          PROVED (R43.3 re-certified)
PAIR MEMORY = ½ [4p(1−p)][4q(1−q)] sin²(Δ/2) ; INVARIANT = 2M                    PROVED
ONE LEDGER ALONE ⇒ NO MEMORY                                                     PROVED
GRAVITATIONAL TURN φ = G m₁m₂τ/(ħd) FROM THE CLOCK FACTOR                         DERIVED in the weak field from GR1 + CL1
NO INTRINSIC SINGLE-MASS GRAVITATIONAL DECOHERENCE                               STANCE of the thesis — follows from statements 1–3, not a theorem about nature
MANY LEDGERS (an environment); WHEN MEMORY BECOMES EFFECTIVELY PERMANENT         NOT BUILT (R44 is the place)
ANYTHING BEYOND ORDINARY QUANTUM MECHANICS WITH THE NEWTONIAN POTENTIAL          NOTHING — the stance is that there is nothing beyond it here
```

## 7. Reproduce

```text
python dc1_memory_between_masses.py
python -m unittest test_dc1_exact
```
