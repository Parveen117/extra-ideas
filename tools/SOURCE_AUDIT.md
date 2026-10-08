# Source audit: two drafts of "new mathematical objects"

Both uploaded LaTeX drafts were read completely. They are recorded as **Draft A** (a catalog of fifteen
objects) and **Draft B** ("uncut horizons", eight objects). Their byte hashes are in
[SOURCE_PINS.json](SOURCE_PINS.json); the originals are not republished. A hash identifies evidence, not truth.

## How the drafts are read

The owner's rule: in these drafts, as in the whole framework, **the seam is the diagonal** — the 45° line, the
local reference on which one knows what remained and what is lost. The drafts' cut J(z) = i·z̄ exchanges the
real and the imaginary component; its fixed line is that diagonal.

Read so, almost every object is a member of one law (RKF theorum/24): S = R + D, F = R − D = Z*𝔧Z, with one
special line F = 0 and one count, the directions past it. [DS1](ds1/DS1_ONE_DIAGONAL.md) proves the common
part and gives the match with the Riemann line and the physics line. The tables say, for each object, which
member it is and what was needed for its formula. The drafts refer to earlier parts that were not supplied;
each object is read from what the drafts themselves state.

An object has a tool of its own only when its corrected statement is not already certified elsewhere.

## Draft A

| # | Object | On the diagonal it is | The draft's formula | Where it stands |
| --- | --- | --- | --- | --- |
| 1 | Seam Jacobian | A product of readings x + ι, x = w²/h (seen x², lost 1). With the winding of RW1, x rises from 0 toward K and never passes it; for K = 1 the factor turns toward the diagonal vector 1 + ι. | The modulus identity is right. The value is not an invariant: for K = 1 it grows from 1 to 2ⁿ. The limit iⁿ assumes w²/h → 0; the draft's own equation gives w²/h → K. | Instance: DS1-D2, D6 |
| 2 | Winding metric | S = R + D with R = 1 and a memory of rank one along the winding. Its flip 1 − wwᵗ/h² has the three signs the draft lists: + inside, 0 on the diagonal \|w\|² = h², − past it. | The signs belong to the flip, not to the displayed form (which is positive). Whole-number windings have derivative zero; the displayed "geodesic equation" is not one. | Instance: DS1-D3 with rank one |
| 3 | Cut cohomology | The count of the whole = the count on the seam + the count relative to it: the law at the level of cell counts. | χ(M) − χ(seam) is right for the pair. A derivative taken only across the seam does not define it. | Certified: Publications EMK-TOP-1 (T2, T4) |
| 4 | Cyclic tensor | Mirror = transpose. The mirror-odd (alternating) part is d ln(I/K) ∧ d ln(J/K): free of the reference. The mirror-even part needs the reference. | "Not symmetric" and the trace are right. The displayed cyclic law is ill-typed. | **Tool: [CY1](cy1/CY1_THREE_READING_CYCLIC_FORM.md)** |
| 5 | Prime lattice | The pairs (p, u_p) of PT1; the diagonal vector 1 + ι is the prime 2. | The spacing estimates are off; the Fourier sum does not converge. | **Tool: [PM1](pm1/PM1_PRIME_TURN_SERIES.md)**, DS1-D2 |
| 6 | Cut character | The harmonics of the cell between the two diagonals (a quarter turn). The cell does not see a unit-free harmonic except the mean; there are four cells. | Vanishing at the non-zero multiples of 4 and the decay are right. The orthogonality relation is not: the sum of squares is 4. | Instance: DS1-D7, PM1-P1 |
| 7 | Seam determinant | FD1's one function det(1 + zX) on the modes kept; each factor 1 + ι·k·ω is a reading with seen 1, lost (k·ω)². | A finite product is fine. The product over all modes does not converge; a constant factor does not change that. The framework's way to the infinite is theorum/28. | Instance: FD1, DS1-D2 |
| 8 | λ-residue | The eigen-reading of the change operator 1 + ι·D on tails with poles. | Zero on its stated domain; the example and the residue theorem need weights. | **Tool: [JR1](jr1/JR1_JET_READING_OF_THE_CHANGE_OPERATOR.md)** |
| 9 | Cut Fourier transform | The transform on the diagonal share: n^(−s)·n^(−(1−s)) = n^(−1) for every s, with equal factors only on rad s = ½. | Inversion and the shift rule are right. The bare sum for ζ(½ + it) has no value; the memory terms are what reach the strip (LAM-3-F2). | Certified: RKF F00-E §8–9, LAM-3 |
| 10 | Soliton winding number | A count that changes only when a zero crosses — the same rule as DS1-D3's count, which changes only when an eigenvalue crosses 1. | The kink example W = w is right. Conservation holds while no zero crosses. | Certified: Publications RST-1 T6, EMK-1 T6, RKF theorum/76 |
| 11 | Seam tangent bundle | The fields tangent to the diagonal. | A known object (log tangent bundle) when n is read as the scaling across the seam. The connection term is declared. | Not kept |
| 12 | Eye operator | Its claim has the shape of the certified reduction: the sign of a difference form decides. | The certified object of that shape is theorum/02's F = S − VV* with the five-matrix B and k_Σ = N₊(B − 1); the draft's operator is not constructed, and a spectrum of zeros is not real. | Certified shape: RKF theorum/02, 28; DS1-D3. RH OPEN |
| 13 | Prime-modular form | A series with a seam law under t → 1/t, fixed point t = 1. | The draft's series has terms of modulus tending to 1. The series with the law is the one on PT1's prime turns. | **Tool: [PM1](pm1/PM1_PRIME_TURN_SERIES.md)** |
| 14 | Seam zeta function | The same series read through its product over primes. | The twist e^(−2πi/p) is not an exact turn and not multiplicative. The expansion in shifted prime sums is right for Re s > 1. | **Tool: [PM1](pm1/PM1_PRIME_TURN_SERIES.md)** |
| 15 | Cut Grassmannian | One cut takes at most one dimension from any reading; a memory of rank r at most r. This is why a matrix of size five suffices in theorum/02. | The condition holds for every k-plane, so it defines no new space; the dimension is k(n−k). | Instance: DS1-D3 (count ≤ rank) |

## Draft B

| # | Object | On the diagonal it is | The draft's formula | Where it stands |
| --- | --- | --- | --- | --- |
| 1 | Seam fractal dimension | "The number of paths is the inverse of the measure weight": all histories over returned histories of two strands is Exp(w²) up to a power, w = a − b. | N = Exp(w²) holds in that count up to a power (the draft's ε rescales the amplitudes). "D = ∞" is then only exponential growth. | Instance: DS1-D5 |
| 2 | Emotional potential | The winding equation dw/dh = K − w²/h is the moment law of the returned count; the sign is that of K; spread between half the mean and the mean. | The conservation law is not defined. | **Tool: [RW1](rw1/RW1_RETURNED_WINDING_COUNT.md)** |
| 3 | Second cut operator | A mirror of mirrors: J ⊗ J is again a mirror, and the flip of a nested reading is the product of the flips — zero as soon as one factor is on the diagonal. | (J⁽ⁿ⁾)² = 1 is right. "J⁽ⁿ⁾ ∘ J⁽ᵐ⁾" should be ⊗. | Certified: Publications EMK-2 T1; DS1-D1 |
| 4 | Seam of the seam | The point where the cut is still: seen × lost is largest and stationary exactly on the diagonal (4RD = S² − F²). It is the diagonal observer. | Read as a set of points of the plane the definition is empty; read on the share it is the diagonal. | Instance: DS1-D2, DO1-P1, QD1-D1 |
| 5 | Eye of the eye | The law of an increasing chain of cuts, as the P_n of theorum/28. | No further definition is given. | Not kept |
| 6 | Path inequality measure | The loop sum of a differential: zero for a flat one. | For df/f it is 4π² times a squared whole number (the winding), also for analytic f with a zero inside. | Certified: winding as in row A10 |
| 7 | Regret operator | The pairing of a strand with the exchanged strand, F = Z*𝔧Z; its weight Exp(−w²) is the squared overlap of the strands — small far from the diagonal, as the draft says. | A derivative is written as a reflection. | Instance: DS1-D1, D5; RW1 |
| 8 | Emotion tensor | — | "Ricci tensor of a potential" is not defined. | Not kept |
| — | Seam measure e^(−w²) (used throughout) | The squared mirror overlap of two strands with amplitudes a, b, w = a − b: 1 on the diagonal, below 1 off it. | Exact in that form. | Instance: DS1-D5 |

## Count

Five tools: [DS1](ds1/DS1_ONE_DIAGONAL.md) (the common structure), RW1, PM1, JR1, CY1.
Of the twenty-three objects: six have a tool of their own or are replaced by one (A4, A5, A8, A13, A14, B2);
eight are instances of DS1 or of a certified stage with their formula corrected (A1, A2, A6, A7, A15, B1, B4, B7);
six are certified elsewhere in a sharper form (A3, A9, A10, A12 as a shape, B3, B6);
three are not kept (A11, B5, B8).
