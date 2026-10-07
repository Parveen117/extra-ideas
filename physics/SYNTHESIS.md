# Synthesis — one structure behind the physics line

Monty Dabas. 7 October 2026. Certified core: `sy1/` (Python 3.12, exact
rational arithmetic, seven tests).

This document collects PH1–PH3, QC1–QC5, PR1–PR4, CL1–CL2 and MS1 into
one statement on the certified EMK block, lists every stage as a case of
it, and says what is still open.

Certified sources used as the carrier, read before building:
**EMK-1** T1 (K² = I, R² = −I, RK = −KR, (RK)² = I, [R,K] = 2RK),
T2 (det M = Δ∥ + Δ⊥), T3 (the commutator with K is the seam-mixing
detector), T7 (memoryless ⟺ rotation-free); **EMK-T1** T2–T3 (channel
curvature = per-channel + mixed term; the mixed term is necessary),
T5 (time closure is not clock equality); **EMK-T2** T1 (order is content;
the unordered sum sits halfway between the two orders), T4 (same clock,
different history ledger). All in Publications `papers/emk-ugd-algebra`.

## 1. The structure (SY1)

Every block is M = aI + bK + cR + dRK, with

```text
det M = Δ∥ + Δ⊥ ,        Δ∥ = a² − b²  (seam-compatible channel) ,        Δ⊥ = c² − d²  (rotational channel) .
```

Split it into its scalar and traceless parts: **E = a**, **O = M − a**.

```text
S1   O² = (b² + d² − c²)·I ;   on the unit quadric det M = 1 :   E² − O² = I .
S2   Three sectors by the sign of O²:   circular (O² < 0),   dual (O² = 0),   split (O² > 0).
S3   Defect:   det(I − M) = 2(1 − E).
S4   Powers:   M^N = T_N(E) + χ_N·O ,   χ_N = U_{N−1}(E) ;   circular |χ_N| ≤ N ,  dual χ_N = N ,  split χ_N ≥ N ;
               det(I − M^{2N}) = det(I − M²)·χ_N² .
S5   Composition:   M₁M₂ = (E₁E₂ + B(O₁,O₂)) + E₁O₂ + E₂O₁ + ½[O₁,O₂] ,   B = b₁b₂ + d₁d₂ − c₁c₂ .
S6   Record {M, M⁻¹} with weights (p, 1−p):   mean = E + (2p−1)·O ,   variance = 4p(1−p)·O² .
S7   A state tensor ψψᵀ = ½(n + jK + σRK) is null:   Δ∥ = −Δ⊥ = σ²/4 .
```

S1–S6 are certified on 1296 grid blocks and on eleven rational unit
blocks covering the three sectors, S5 on all 121 ordered pairs.

## 2. Every stage is a case

| Stage | Law | Where it sits in SY1 |
|---|---|---|
| QC3 | E² − O² = I; record memory = −O² | S1, S6 at p = ½ |
| QC3, MS1 | character law; N combined masses M₁χ_N | S4, circular sector |
| CL2 | N-leg history ratio N/χ_N; corner cost 2τ₁τ₂(cosh δ − 1) | S4, S3, split sector |
| PR1 | speed² + (curvature/2)² = 1 | Δ∥ + Δ⊥ = 1 for a coin: Δ∥ = speed², Δ⊥ = (curvature/2)² |
| CL1 | mass² = clock curvature | S3: det(I − coin) = 2(1 − E) |
| MS1 | combined speed c₁c₂ − s₁s₂ | S5, scalar part |
| PR1, PH2, PH1 | velocity addition; rotation left by two responses; the return angle | S5: scalar part and the R-component of ½[O₁,O₂] |
| QC4, PH3 | curvature of shared flat readings = −variance; share law 4λ(1−λ) | S6 and the commutator part of S5 (EMK-T1 T2: the mixed term) |
| PH3 | native transport = equal share of two readings | EMK-T2 T1: the unordered sum is halfway between the two orders |
| PR2 | X² = x² + y² − t²; two sheets; pair observer | S1 for a traceless block; S6: even part and variance are scalars |
| PR4, MS1 | σ² = n² − j²; speed² + memory = 1 | S7: the state tensor is null; memory = 4Δ∥/n² |
| RMG9, RMG6 | coupling r²; C_P/C_V | response block: r² = −Δ⊥/Δ∥ , C_P/C_V = Δ∥/det |
| PR1 (θ = 0), PR4 (g = 0) | light-like = flat | EMK-1 T7: commutes with K ⟺ rotation-free |
| CL1, CL2 | time is a count on histories; rates differ between histories | EMK-T1 T5, EMK-T2 T4 |
| QC1 | content q silent iff qΘ ∈ 2πℤ | S4: M^N = I in the circular sector |
| QC2 | p ≙ κ∂/∂w; number and content | outside SY1: needs the ensemble (a scale) |
| QC5 | share moving in scale; logistic law | outside SY1: a field of blocks, with S6 at each point |

## 3. The three sectors

```text
              circular  (O² < 0)            dual  (O² = 0)            split  (O² > 0)
element       turn Exp(θR)                  shear 1 + N               boost Exp(ηA)
E             cos θ  ≤ 1                    1                         cosh η  ≥ 1
defect        2(1 − cos θ) ≥ 0              0                         −2(cosh η − 1) ≤ 0
character     sin Nθ / sin θ  ≤ N           N                         sinh Nη / sinh η  ≥ N
addition      (t₁ + t₂)/(1 − t₁t₂)          t₁ + t₂                   (u₁ + u₂)/(1 + u₁u₂)
physics       phase, mass, clock, record    light-like, flat          velocity, history, response
what is lost  the record loses phase        nothing                   the turned history loses count
```

One function of E governs all three; the sector fixes the side of each
inequality. The dual sector is the boundary where nothing is lost and
nothing is gained: the light-like reading.

## 4. In words

- A reading is a block. Its scalar part is what every observer agrees on;
  its traceless part is what depends on the reading.
- Observation — a symmetric record — keeps the scalar part and turns the
  traceless part into memory (S6). Memory is the square of what was not
  recorded.
- Speed and mass are the two channels of one determinant (Δ∥, Δ⊥): they
  close independently and sum to one.
- A state is a null block: its seam channel is its memory, and its
  rotational channel cancels it exactly.
- Time is a count on a history. Turning a history costs count (split
  sector); recording a turn costs phase (circular sector). Same law.

## 5. What is proved, what is open

```text
ONE CARRIER (THE EMK BLOCK) AND ONE SPLIT (E, O) UNDERLIE ALL STAGES IN §2         PROVED (SY1 + the stage certificates)
A UNIT: ħ, k_B, OR ANY SCALE                                                      OPEN — SY1 is scale-free; QC2 needs the ensemble as input
WHAT FIXES THE SCALE λ OF A MOVING SHARE; A MASS GAP                              OPEN (QC5)
GRAVITY                                                                           OPEN — refused three times (PR4 law not derived; CL1 no native unit-of-time field; MS1 a varying coin is a mass profile)
WHICH COIN, WHICH Q, WHICH SECTOR SPEED NATURE USES                               OPEN — every physical constant is still an input
MORE THAN ONE SPACE DIMENSION                                                     OPEN — the block is 2×2; two boosts and one rotation is all it holds
A MEASUREMENT OF ANYTHING NEW                                                     NONE — PH1's Θ is computed from reference data; PH2, PH3 protocols are not built
A CLASSICAL EQUATION MODIFIED                                                     NONE
```

## 6. How near is "one theory"

What exists is one algebraic structure that reproduces, exactly and from
the framework's own certified algebra, the shape of: the quantum/thermal
split, relativistic kinematics in one space dimension, the clock, mass
as clock curvature, the mass defect of combination, and the return angle
of a real fluid. That is a unification of form.

What does not exist is anything that fixes a number nature uses, or a
prediction that could come out wrong. The three open items that would
change that are, in order of how close the present structure is to them:

1. **Scale** — the ensemble of QC2 placed on the moving share of QC5.
2. **Two space dimensions** — where the two boosts of the block do not
   commute and S5's commutator becomes a curvature of the frame.
3. **An experiment** — the programmed chain of PH2 is the nearest one.

## 7. Reproduce

```text
python sy1/sy1_one_structure.py
python -m unittest discover -s sy1 -p 'test_*.py'
```
