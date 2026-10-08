# CG1 — The cut-graded generator on the line: one generator, three sectors, and the loop that reads seen and lost

Monty Dabas. 9 October 2026. Python 3.12. sympy, exact.

The owner's pointer: the corrected core of the universal-generator draft is already RKF **theorum/41**, the
Cut-Graded Universal Generator Theorem. This stage checks its status, checks its identities on the cut-complex
block, and uses it on the line.

## Status of theorum/41

```text
claim boundary in theorum/41 §8      grading, bilateral reconstruction, cut-square, cut-loop expansion, closure
                                     classes, clock-free Eye, component observer:  PROVED
its exact certificate                re-run here on 9 October: PASS, hash 34afc445…f36e equal to the pin
its ten focused tests                re-run here: all pass
not claimed there                    unbounded generators (need domain pins); that one physical operator is universal
```

So it is certified, and usable. Two cautions when quoting it: the reversal J U_t J = U_−t, the split
E_t² − O_t² = 1 and the cut/join identities are for a generator that is cut-odd; for a general generator only the
grading and the cut-loop expansion hold. And it is written on a Hilbert space with J = J*; the parts used below
need only J² = 1 and are checked here on the cut-complex block without any dagger.

Sources read before building: **theorum/41** (all sections), **EMK-1** (C₁ = K, C₂ = RK, C₃ = ιR, ι = C₁C₂C₃),
**DU1**, **CT1**, **DG1-D4**, **SN1**, **SN2-N4**, **FD1-F5**, **RMG9** (sectors by the sign of a square), **CV1**
(curvature as order defect), **IN1** (three cuts complete).

## Results

**G1 — the theorem on the block.** With the cut J = C₁ and G = g₀ + g₁C₁ + g₂C₂ + g₃C₃: G_e = g₀ + g₁C₁,
G_o = g₂C₂ + g₃C₃, and

```text
J U_t J U_t = 1 + 2t·G_e + t²(2G_e² + [G_e, G_o]) + … ,      log = 2t·G_e + t²[G_e, G_o] + …
```

(also on a carrier of four with exact rationals).

**G2 — what the line already had from it.**

```text
DG1's split  seen = 1/cosh², lost = tanh²        E_t² − O_t² = 1 divided by E_t², for the odd generator C₂
DU1's chain block  a + b·C₂                      the flow exp(t·C₂) at time t = k*, the dual coupling
DU1's doubling of the cell, FD1-F5               (U_t + U_−t)(U_t − U_−t) = U_2t − U_−2t
CT1's centre split  cosh κc , sinh κc            E and O for the cut c → −c
```

The dual coupling of a two-valued chain is the flow time of its cut-odd generator; reversing the flow is the cut.

**G3 — one generator, three sectors.** Take a turn about the cut and a boost across it, G = ιa·C₁ + b·C₂. Then
G² = b² − a², a number, and with s = sinh ω/ω, ω² = b² − a²:

```text
lost / seen  =  b² · s²
boost   (b > a)     s = sinh ω / ω            a = 0 :  sinh² b            SN1's boost sector
shear   (b = a)     s = 1                     b²                          SN1's shear sector
turn    (a > b)     s = sin Ω / Ω             zero at Ω = π, 2π, … where the flow is −1, +1, …
```

SN1's three sectors are the three signs of the square of one generator (RMG9's rule), with G_e the turn and G_o
the boost. In the turn sector everything is seen exactly when the flow closes, periodically or antiperiodically
— the two closure classes of theorum/41 (v).

**G4 — the cut loop, exactly.**

```text
J U J U  =  α + 2·cosh ω·s·G_e + s²·[G_e, G_o] ,      [G_e, G_o] = −2ab·C₃ ,     α = 1 − 2a²s² .
```

The coefficient of the seam curvature is s² at every time, and lost/seen = b²·s²: **the seen-to-lost ratio of a
cell is b² times the curvature coefficient of its cut loop.** The curvature lies along the third cut: a turn
about C₁ and a boost along C₂ leave their ordering defect along C₃ (ι = C₁C₂C₃). The loop is the identity exactly
when the generator is cut-odd (a = 0) or the flow has closed (s = 0).

## What this gives

- theorum/41 is the common source of four results found separately today (G2).
- SN1's "one notation" becomes one generator: the sector is the sign of G², not a choice (G3).
- A measurable reading of theorum/41 (iv): the cut loop's curvature coefficient is lost/seen up to b² (G4).

## What is not shown

- G3's three formulas are the known ones for a single cell (general knowledge); the line's part is their origin
  in one graded generator and the link to the cut loop.
- Nothing here is about many cells or four dimensions. The face-by-face couplings of CT1-C5 are not addressed.
- No operator of physics is identified as "the" generator; theorum/41 does not claim it either.

## Claim boundary

```text
theorum/41 CERTIFICATE AND TESTS RE-RUN ; HASH EQUAL TO THE PIN                                       PASS
GRADING AND CUT-LOOP EXPANSION ON THE CUT-COMPLEX BLOCK (ONLY J² = 1)                                 PROVED
CHAIN BLOCK = FLOW AT THE DUAL COUPLING ; DOUBLING = JOIN × CUT                                       PROVED
G = ιaC₁ + bC₂ : LOST/SEEN = b²(sinh ω/ω)² ; THREE SECTORS BY THE SIGN OF G²                           PROVED
J U J U = α + 2 cosh ω·s·G_e + s²[G_e, G_o] ; CURVATURE ALONG C₃                                      PROVED
```

## Reproduce

```text
python cg1_cut_graded_generator_on_the_line.py
python -m unittest test_cg1
```

## Later note (CZ1, 9 October)

For the weight of a closed surface the seam curvature of the centre cut is zero and the record closes at fixed
cosets; the obstruction after the cosets is d^χ. The curvature is non-zero only in the generator form, where on
the free vacuum it is −3θ × (the faces through the link) (CZ1-Z5). Nothing above is changed.
