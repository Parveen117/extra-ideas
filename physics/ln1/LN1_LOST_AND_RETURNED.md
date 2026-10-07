# LN1 — What a linear reading loses, the first non-linear layer returns

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic.

The owner's statement: an observation keeps only part of the
information (spectral blindness); the λ-tower moves toward the whole
without a new observation; the relation between the two should be the
bridge between linear and non-linear.

Sources read before building: **EMK-1 T2** (det = Δ∥ + Δ⊥, the
observed and lost channels), **NT-1**, **RMG1 Corollary 2.1**
(spectral blindness), **RMG2** (the λ-tower; spectral operations never
turn the axes), RKF **Theorem 42 §7** (a higher layer repairs what a
lower one is blind to), **QC4** (curvature of shared readings =
−variance), **GB1**, **NC1**, **PH3**.

A reading in a declared cut sees the diagonal entries a_i of the
response element H. The couplings b_ij are lost to it.

## Results

**N1 — the bridge identity.**

```text
( H² )_ii  =  a_i²  +  Σ_j b_ij² .
```

The non-linear layer, read in the same cut, is the square of the linear
reading plus exactly what the linear reading lost:

```text
lost_i = ( H² )_ii − ( H_ii )²        — the variance of the reading.
```

**N2 — every level carries it.** For any λ ≠ 0,
lost_i = ( T_λ(H)_ii − a_i − λ a_i² ) / λ. Any level of the tower can
be the reference.

**N3 — two cuts.** The linear reading gives (a, c). One generation adds
b². Then EMK-1's two channels are both known, det = ac − b², and with
them the trace, the spectrum and the anisotropy ℓ. The sign of b is not
returned by any generation.

**N4 — three cuts.** The second layer returns all three b_ij². The
third layer returns the sign of b₁₂b₁₃b₂₃. Changing two signs stays
hidden.

**N5 — the gravity element.** Read in a cut across its axis,
G = Exp(ψ n·C) shows cosh ψ in both slots — blind to the direction —
and loses sinh²ψ. One generation gives

```text
( G² )_ii − 1 = [ cosh²ψ − 1 ] + sinh²ψ = 2 sinh²ψ :
```

the part that was seen and the part that was lost are exactly equal.
This is the same one-to-one split as in GB1, where gravity through the
clock alone bends light by half: 0.876″ + 0.876″ = 1.751″ at the Sun.

## Reading

```text
linear reading in a cut        keeps a_i, loses the couplings                (spectral blindness, in a cut)
first non-linear layer         returns the lost part as a variance           (N1)
further layers                 return signs of closed products               (N4)
never returned                 a single sign — which sheet                   (N3)
```

So the two things the owner named are one identity read in two
directions. Observation in a cut splits the element into seen and lost
(EMK-1). The tower, without a new cut, adds the lost part back to what
is seen. Moving up the tower is moving toward the uncut element.

On counting: the lost part is a variance, and by QC4 the curvature of
shared readings is minus the variance — the quantity whose flux WQ1
counts in whole units. That the part lost to a cut is the part that is
counted is suggested by these three results together; it is not derived
here.

## What is not shown

- N1 is elementary algebra. What is the framework's own is the reading:
  EMK-1's lost channel is the first tower layer's excess.
- Convergence "to full information" is proved for two and three cuts.
  For more cuts the number of layers needed was not determined.
- N5's link to the bending of light is a match of two equal splits, not
  a derivation of the bending from the tower.
- The link to counting (last paragraph of the reading) is not derived.

## Claim boundary

```text
(H²)_ii = a_i² + Σ b_ij² ; LOST = VARIANCE OF THE READING                 PROVED
EVERY LEVEL OF THE TOWER CARRIES THE LOST PART                           PROVED
TWO CUTS: ONE GENERATION RETURNS b², HENCE det AND THE SPECTRUM          PROVED
THREE CUTS: SQUARES AT LAYER 2, SIGN OF THE TRIPLE PRODUCT AT LAYER 3    PROVED
ONE SIGN NEVER RETURNED                                                  PROVED
GRAVITY ELEMENT: SEEN PART = LOST PART EXACTLY                           PROVED
LOST PART = THE PART THAT IS COUNTED                                     NOT DERIVED
GENERAL NUMBER OF CUTS                                                   NOT BUILT
```

## Reproduce

```text
python ln1_lost_and_returned.py
python -m unittest test_ln1_exact
```
