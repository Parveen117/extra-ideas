# PR4 — What a sector can offer the clock: flip rate × proper density, and proper density is memory

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues PR3.

PR3 showed that acceleration is a position-dependent unit of time N with
N″ = 0, and left the gate: what fixes N when N″ ≠ 0. The hypothesis put
to this stage: **stored record memory sources the curvature of the clock.**

Sources read before building: PR1-T4, PR2-T4, PR3-T2/T6, QC4-T1
(return = −variance), CFE-1 (loop residue as the obstruction).

## 1. The two currents and the residue

Sector law (c = 1): ∂_tψ = (A∂_x + gR)ψ. Put

```text
n = ψ·ψ   (density)        j = ψ·Aψ   (current)        σ = ψ·Kψ = 2ψ₁ψ₂        τ = g σ .
```

**T1.**

```text
∂_t n = ∂_x j                       conserved
∂_t j − ∂_x n = −2 τ                the dual current fails to close by exactly 2τ
```

At g = 0 both close: the light-like sector has two conserved currents and
no residue. The mass term breaks the second one, and the amount is τ.

## 2. The residue is a scalar, and it is memory

**T2.** Under a boost (n, j) turns as a vector, σ is unchanged, and
σ² = n² − j²: σ is the proper density.

**T3 (memory).** Let p = ψ₁²/n be the share of the first light-like
reading. Then

```text
j = (2p − 1) n ,          σ² = n² · 4 p (1 − p) .
```

The proper density is the density times the root of the two-way memory
4p(1−p) of the two light-like readings (the share law of PH3 and QC4).
A pure reading — all weight in one mover — has zero memory, σ = 0, τ = 0.

**T4 (uniqueness and sheet parity).** Every boost-invariant bilinear of ψ
is a multiple of σ. Under the sheet map (ψ → Aψ, g → −g), n and j are
even, σ is odd, and τ = gσ is even: mass and antimass offer the same
source.

**T5 (the frame does not change it).** In the accelerated frame of PR3,

```text
∂_η n = ρ ∂_ρ j ,          (1/ρ) ∂_η j − ∂_ρ n = −2 τ .
```

The residue in local units is the same τ.

So the only scalar, sheet-even, local quantity a sector has to offer is

```text
τ = (flip rate) × (density) × √(memory of the two light-like readings) .
```

## 3. The hypothesis

**What is established.** If the clock is curved by anything this sector
carries, and the rule is local, frame-independent and the same for both
sheets, the source is a function of τ — the stored memory times the flip
rate. In that sense the hypothesis has exactly one candidate, and it is
the memory.

**What is refused.** A derivation of a law N″ = κ·τ·N from the sector
law. Two exact obstacles:

- PR3-T6: the sector law holds for every N. It contains no equation for N.
- QC4-T1 gives nothing here: with one space dimension the frame
  connection has a single generator A, and the variance term [A, A]
  vanishes identically. "Curvature = −variance" needs two boosts that do
  not commute.

With the law taken as a postulate (the simplest one; known in the
literature for one space dimension — general knowledge, unchecked), its
consequences inside this model are immediate from T1–T5: a pure reading
leaves the clock flat (N″ = 0, at most uniform acceleration); mass and
antimass curve it the same way; the weights between the sheets do not
enter.

## 4. Certificate

T1 and T5 on three two-component test fields with half-integer
exponents, with ∂_tψ replaced by the law; the g = 0 closure; T2–T3 at
three rational states and boosts (a pure reading included); T4 by the
invariance equations for two boosts and by the sheet map on the test
fields. Five tests; a split coupling in place of the circular one is
rejected.

## 5. Claim boundary

```text
ONE CONSERVED CURRENT; DUAL RESIDUE = −2τ                                 PROVED
σ SCALAR, σ² = n² − j² = n²·4p(1−p): PROPER DENSITY = DENSITY × √MEMORY    PROVED
τ UNIQUE (UP TO FUNCTIONS OF IT), SHEET-EVEN, FRAME-INDEPENDENT            PROVED
STORED MEMORY SOURCES THE CLOCK CURVATURE (N″ = κ τ N)                     NOT DERIVED — one candidate identified, law is a postulate
GRAVITY IN MORE THAN ONE SPACE DIMENSION                                  NOT TOUCHED
```

## 6. Next gate

Two space dimensions, where the algebra has two non-commuting boosts
(K and S) and their commutator is the rotation R. There QC4-T1 is not
empty: the curvature of a shared frame is minus the variance of its
boosts — the same object as the return angle of PH1. That is where a
derivation, or a refusal with a witness, can be had.

## 7. Reproduce

```text
python pr4_source_is_memory.py
python -m unittest test_pr4_exact
```
