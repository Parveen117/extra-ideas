# DU1 — The diagonal turn exchanges the two descriptions of a record; the rate is the dual coupling

Monty Dabas. 9 October 2026. Python 3.12. sympy and exact enumeration; numpy for finite strips (evidence only).

MG1 ended at a wall: a small mass count needs a rate that falls as an exponential of the coupling, and the chain
of turns gives only a power. The owner's statement — on the diagonal quantum and classical are the same — is
taken here as an instruction: read the record on the diagonal.

Sources read before building: **EMK-1** (cuts C₁ = K, C₂ = RK), **DO1-P1** (a block read at 45°: mean and
half-difference), **DG1-D4** (memory = tanh², seen + lost = 1), **CL1/HB1** (1/√2), **YM-15, YM-21** (decay of a
chain is exactly geometric, r^(cells); rate −Log r), **MG1**, **QD1**.

The record: a chain of two-valued marks with coupling k. Its block is a + b·C₂, a = e^k, b = e^(−k).

## Results

**U1 — the diagonal turn.** H = (C₁ + C₂)/√2 is its own inverse and exchanges the two cuts: H C₁ H = C₂. Read on the
diagonal the block is diag(a + b, a − b): DO1's mean and half-difference.

C₁ reads which mark — the description by configurations. C₂ reads the flip — the description by contents. The
diagonal turn exchanges the two descriptions.

**U2 — the dual record, and the record that is its own dual.** The diagonal reading is again a block of the same
kind with coupling k*, e^(−2k*) = tanh k. Then (k*)* = k and sinh 2k · sinh 2k* = 1. The record equal to its dual:

```text
sinh 2k = 1 ,   k = ½ ln(1 + √2) = 0.4407 ,   b/a = tan(π/8) ,
tanh² 2k = ½ = 1/cosh² 2k :      seen = lost ,  tanh 2k = 1/√2 .
```

The block points along the fixed direction of H, half-way between the cut and the diagonal. In DG1's split
(memory = tanh²) the self-dual record is exactly the diagonal seen = lost, with the eight-mark value. There the
two descriptions are one record: the owner's statement, for this record.

**U3 — the chain.** On a ring of N marks, exactly, ⟨s₀sₙ⟩ = (tⁿ + t^(N−n))/(1 + t^N), t = tanh k. So

```text
rate per cell  =  −Log tanh k  =  2 k* .
```

What the cut description calls the rate, the diagonal description calls the coupling. Doubling the cell (summing
over the middle mark) squares the ratio (YM-21's law) and

```text
k*  →  2 k*        exactly
k   →  ½ ln cosh 2k  =  k − ½ ln 2 + …        a log: each doubling of the cell takes ½ ln 2 from the coupling.
```

**U4 — the exponential law.** 2e^(−2k) < rate < 2e^(−2k)/(1 − e^(−4k)). A mass count of 7.7 × 10⁻²⁰ per cell is
k = 22.4; the electron's is k = 26.1. Numbers of order 20, not 10²⁰. For the chain of turns the same doubling
halves the coupling (0.5002 at κ = 1000): a power, no log.

**U5 — the plane.** On grids of 2×2, 2×3, 3×3, 3×4 marks the even subgraphs and the cuts of the dual grid are
equal in number, size by size (exact enumeration). So the record of a plane at k and its dual at k* are one
function, with the same map as U2.

**U6 — known, with evidence.** For the infinite plane the rate is 2|k* − k| (known result; not derived here):
zero on the diagonal. Rings of width 4 … 10: at k = 0.3 the rate is 0.6338 at width 10 against 0.6334; on the
diagonal width × rate falls to 0.7887 at width 10, toward π/4.

## What this gives

```text
a native exponential law          two-valued chain: mass count per cell = 2k* ≈ 2e^(−2k)        PROVED
where the log comes from          doubling the cell: k → k − ½ ln 2                              PROVED
what a small mass count is        closeness to the diagonal: 2k* (chain), 2|k* − k| (plane, known)
why the turn chain gave a power   continuous turns: doubling halves the coupling                 PROVED (numerically)
```

## What is not shown

- This is the two-valued record, not the colour record, and not four dimensions. The numbers 22.4 and 26.1 say
  what size of number an exponential law needs; they are not a derivation of any mass.
- The chain and plane results are the known ones for two-valued marks (general knowledge); 1/k = 2.2692 at the
  self-dual point is the known value. The line's part: the duality is DO1's diagonal turn, the self-dual record is
  seen = lost with the eight-mark value, and the rate is the dual coupling.
- The colour block has a two-valued record inside it: its centre (YM-39). Whether the doubling of that record
  carries the log in four dimensions is not tested. The two-valued gauge record in four dimensions is known to be
  its own dual under the same map (recalled, not re-read at source).
- U6 is evidence on finite rings; no statement about the infinite plane is proved here.

## Claim boundary

```text
H = (C₁ + C₂)/√2 EXCHANGES THE CUTS ; DUAL COUPLING e^(−2k*) = tanh k ; INVOLUTION                   PROVED
SELF-DUAL ⇔ sinh 2k = 1 ⇔ SEEN = LOST (tanh² 2k = ½) ; BLOCK ALONG π/8                              PROVED
CHAIN: RATE = 2k* ; DOUBLING: k* → 2k*, k → k − ½ ln 2 + … ; RATE ≈ 2e^(−2k)                         PROVED
PLANE DUALITY ON GRIDS TO 3×4                                                                       PROVED (enumeration)
PLANE RATE 2|k* − k| ; π/4 ON THE DIAGONAL                                                          KNOWN ; finite rings consistent
THE PROTON'S MASS COUNT                                                                             NOT DERIVED
```

## Reproduce

```text
python du1_diagonal_turn_is_the_duality.py
python -m unittest test_du1
```

## Later note (FD1, 9 October)

The chain's rate is the log-ratio of the two zeros of det(1 − zB); the diagonal turn leaves that function
unchanged, and doubling the cell is det(1 − z²B²) = det(1 − zB)·det(1 + zB). Nothing above is changed.

## Later note (CT1, 9 October): "the turn chain gives a power" is the same law in another coupling

The table above says the chain of turns gives a power because its turns are continuous. CT1 finds the two-valued
record inside that chain, its centre, with tanh k_c = I₂(κ)/I₁(κ) and k_c = ½ ln(4κ/3) + …. In k_c the chain of turns
obeys the exponential law of U4; the power in κ is that law rewritten. The statements U1–U6 are unchanged.
