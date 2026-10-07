# PR3 — The sector law in an accelerated frame: acceleration changes the unit of time, and nothing else

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues PR2.

PR2 left a question open (its T6): in a truly accelerated frame, does a
split term survive next to the mass turn — a critical acceleration at
which the turn stops — or is it absorbed? And does the pair observer's
invariance hold there? This stage decides both.

Sources read before building: PR1-T4 (sector law, own boost), PR2-T4
(pair record), QC5 (the scale variable s = log r), RMG9-T2 (sectors by
the sign of a square).

## 1. Setting

Inertial law, c = 1: (∂_t − A∂_x − gR)ψ = 0. Uniformly accelerated frame:
x = ρ cosh η, t = ρ sinh η. In null coordinates u = x + t, v = x − t this
is ρ² = uv, e^{2η} = u/v. All identities are certified on monomials
u^a v^b with rational exponents.

## 2. Results

**T1 (boost alone leaves a split term).** With ψ = Exp(−ηA/2)φ,

```text
(1/ρ) ∂_η φ = [ A ( ∂_ρ + 1/(2ρ) ) + g R ] φ .
```

The local boost produces exactly the split term A/(2ρ): half the boost
rate, as PR2-T6 guessed.

**T2 (it is the half-density, and it is absorbed).** With φ = ρ^{−1/2}χ —
equivalently ψ₁ = u^{−1/2}χ₁, ψ₂ = v^{−1/2}χ₂, each light-like component
divided by the root of its own null coordinate —

```text
(1/ρ) ∂_η χ = ( A ∂_ρ + g R ) χ .
```

This is the inertial law with ∂_t replaced by (1/ρ)∂_η. No split term.
The whole effect of uniform acceleration on the sector is a
position-dependent unit of time.

**T3 (what the frame observer reads).**

```text
rest turn per unit frame time     g ρ
coordinate speed                  ρ
ratio between two heights         ρ₁ / ρ₂          (the same for the turn and for the speed)
at ρ → 0                          no rest turn per frame time: every sector is light-like there
```

**T4 (the turn is never stopped).** With s = log ρ the generator is
G = A∂_s + gρR and

```text
G² = ∂_s² − g²ρ² − gρ·K .
```

On the two exchange combinations χ₁ ± χ₂ the operator −G² factorizes,
with W = gρ:

```text
(∂_s + W)(−∂_s + W)     and     (−∂_s + W)(∂_s + W) .
```

The minus branch has a well W² − W of depth exactly 1/4 at gρ = ½, inside
the distance 1/g from ρ = 0; the factorization shows the well is exactly
compensated. PR2-T6's reading — a critical acceleration at which the
rest turn stops — is **refused**.

**T5 (both sheets share the clock).** A·G_g·A = G_{−g}: the reversed turn
obeys the same frame law with the same factor ρ. Mass and antimass see
the same change of the unit of time. With PR2-T4, whose frame changes
here are a boost at each point times a positive scalar, the pair
observer's even part and memory remain scalars in the accelerated frame;
the balanced pair reads g × (its own time) and nothing frame-dependent.

**T6 (general clock factor).** N^{1/2}∂N^{1/2} = N∂ + N′/2: for any
position-dependent clock factor N the split term is N′/2 and is absorbed
by N^{−1/2}. Uniform acceleration is N = ρ.

## 3. Reading

- In this sector model an accelerated frame is not a new force and not a
  new term: it is one position-dependent scale on the whole generator.
  Speed and rest turn scale together, so their ratio — the sector's
  coin — is the same at every height.
- The curvature of the coin against the cut (PR1) is untouched by
  acceleration. What changes from place to place is how much frame time
  one turn takes.
- The place where the factor vanishes is where every sector becomes
  flat and light-like in frame units.

## 4. Certificate

T1, T2 on three two-component test fields with half-integer exponents,
including the checks that the split term is needed without the
half-density and absent with it. T4: the square law, both
factorizations, the well values. T5: the conjugation identity. T6 on
sixteen (n, k) pairs. Eight tests; a wrong weight and a wrong generator
are rejected.

## 5. Claim boundary

```text
SPLIT TERM A/(2ρ) FROM THE LOCAL BOOST; ABSORBED BY THE HALF-DENSITY     PROVED
ACCELERATED LAW = INERTIAL LAW WITH ∂_t → (1/ρ)∂_η                       PROVED
REST TURN AND SPEED SCALE BY THE SAME FACTOR ρ                           PROVED
CRITICAL ACCELERATION STOPPING THE TURN (PR2-T6 reading)                 REFUSED
BOTH SHEETS SHARE THE CLOCK; PAIR OBSERVER INVARIANCE IN THE FRAME       PROVED
GRAVITY: A CLOCK FACTOR WITH N″ ≠ 0, FIELD EQUATIONS, WHAT SOURCES N     NOT TOUCHED
A TEMPERATURE FOR THE ACCELERATED OBSERVER                               NOT CLAIMED — needs a unit (QC2), as the information-invariance ledger's Theorem L says
ONE SPACE DIMENSION ONLY                                                 YES
```

## 6. Reproduce

```text
python pr3_accelerated_frame.py
python -m unittest test_pr3_exact
```
