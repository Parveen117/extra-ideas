# Source audit: two drafts of "new mathematical objects"

Both uploaded LaTeX drafts were read completely before any tool was chosen. They are recorded as **Draft A**
(a catalog of fifteen objects) and **Draft B** ("uncut horizons", eight objects). Their byte hashes are in
[SOURCE_PINS.json](SOURCE_PINS.json); the originals are not republished. A hash identifies evidence, not truth.

The drafts refer to earlier parts (cut complex analysis, seam equation theory, cut number theory, …) that
were not supplied. Each object was therefore judged on what the drafts themselves state, against the
repositories as they stand.

Rule applied: an object is kept only if a corrected statement (1) is true and checked here and (2) is not
already certified elsewhere. Whether a stage of the line uses it yet is said in each tool under "Use": RW1 and
PM1 bear on existing stages; JR1 and CY1 have no user yet.

Where a draft's own statement is right, the row says so.

## Draft A

| # | Object | Assessment | Where it stands |
| --- | --- | --- | --- |
| 1 | Seam Jacobian | Determinant of a diagonal matrix. The product formula drops the factors f_k. "Topological invariant" does not hold: the value changes continuously with w²/h. Its modulus identity is the norm law of RKF F00-E. | Not kept |
| 2 | Winding metric | Whole-number windings have derivative zero, so the defining integral vanishes. The stated form δ + wwᵗ/h² is positive definite, not "indefinite (+,−,0)". The displayed "geodesic equation" is linear in the velocities; a geodesic equation is quadratic in them, and for a metric constant in θ it is θ″ = 0. | Not kept |
| 3 | Cut cohomology | A derivative taken only across the seam does not give the stated relative cohomology. The Euler-number formula χ(M) − χ(seam) is the right one for the pair (known). | Not kept; seam cohomology is certified in Publications EMK-TOP-1 (T2, T4) |
| 4 | Cyclic tensor | The draft's "not symmetric" and its trace are right. The displayed cyclic law is ill-typed. The exact split, and which part depends on the reference, were missing. | **Kept, completed: [CY1](cy1/CY1_THREE_READING_CYCLIC_FORM.md)** |
| 5 | Prime lattice | The spacing estimates are off (the gap of log p is about (log p)/p). The Fourier sum does not converge and is not ζ(½+iu). The native pair is (p, u_p) of PT1. | Replaced: [PM1](pm1/PM1_PRIME_TURN_SERIES.md) reads PT1's pairs |
| 6 | Cut character | The Fourier coefficient of a quarter-turn window. It vanishes at the non-zero multiples of 4 (at k_j = 0 the factor is 1) and tends to zero, as the draft says. The orthogonality relation is false: the left side is χ(k)²δ. | Not kept: a window sum, nothing further |
| 7 | Seam determinant | The product over all modes diverges; a constant factor e^(−w²) does not change that. The theorem uses an undefined symbol. | Not kept |
| 8 | λ-residue | Zero on its stated domain; example and residue theorem do not hold as written. On tails with poles it is an exact eigen-reading. | **Kept, corrected: [JR1](jr1/JR1_JET_READING_OF_THE_CHANGE_OPERATOR.md)** |
| 9 | Cut Fourier transform | The Fourier transform in the logarithm of the scale (Mellin), which is not new; the inversion formula and the shift rule ("diagonalization") are right. ζ(½+it) as the transform of Σ n^(−½)δ(ω−n) is a divergent sum. | Not kept; the log-scale chart is RKF F00-E §8–9 and GE1's audit |
| 10 | Soliton winding number | The winding number of a loop that avoids zero; the kink example W = w is right. It is conserved only while no zero crosses the loop. No equation is given against which the "kink" could be checked as a solution. | Not kept; certified more sharply in Publications RST-1 T6 (sf = Wind = ind), EMK-1 T6, RKF theorum/76 |
| 11 | Seam tangent bundle | With n read as a normal vector, the fibre as written is the whole tangent space. With n read as the generator of scaling across the seam, it is the known bundle of fields tangent to the seam (log tangent bundle), and both of the draft's bullets hold for it. The connection term is declared. | Not kept: a known object, not a new one |
| 12 | Eye operator | Neither operator is constructed. "Spectrum = zeta zeros" and an equivalence with RH are asserted; the zeros are not real, a self-adjoint spectrum is. | Not kept. RH is OPEN in every ledger; this object must not be cited |
| 13 | Prime-modular form | The series converges for no τ: its terms tend to modulus 1. | Replaced: [PM1](pm1/PM1_PRIME_TURN_SERIES.md) (P6, P5) |
| 14 | Seam zeta function | The twist is not an exact turn, is not multiplicative, has no conductor. No s ↔ 1−s law is known or derived for it; the underlying prime sum has branch points inside the strip (classical). The expansion in shifted prime sums is right for Re s > 1. | Replaced: [PM1](pm1/PM1_PRIME_TURN_SERIES.md) (P6, P2–P4) |
| 15 | Cut Grassmannian | Every k-plane meets a hyperplane in dimension ≥ k−1, so the set is the whole Grassmannian and its dimension is k(n−k), not k(n−k)−1. | Not kept |

## Draft B

| # | Object | Assessment | Where it stands |
| --- | --- | --- | --- |
| 1 | Seam fractal dimension | The path count N(ε) = exp(w²/ε) is asserted, not derived; "D = ∞" is arithmetic on the assertion. | Not kept |
| 2 | Emotional potential | Its equation dw/dh = K − w²/h has an exact regular solution and a meaning, which the draft does not give. On that branch the sign is that of K (for K > 0 never negative, as the draft says for "forward"). The conservation law is not defined. | **Kept, completed: [RW1](rw1/RW1_RETURNED_WINDING_COUNT.md)** |
| 3 | Second cut operator | J ⊗ … ⊗ J is again an involution ((J⁽ⁿ⁾)² = I is right). "J⁽ⁿ⁾ ∘ J⁽ᵐ⁾ = J⁽ⁿ⁺ᵐ⁾" composes maps on different spaces. | Not kept; the parity of grades is Publications EMK-2 T1 |
| 4 | Seam of the seam | J is real-linear, so ∇J is nowhere zero and the set as defined is empty; ±1 is not on the 45° line. The seam meets the unit circle at the eighth mark, which needs √2. | Not kept; PT1-P3 |
| 5 | Eye of the eye | No definition from which a statement could be checked. | Not kept |
| 6 | Path inequality measure | For two paths with common ends and a function without zeros on them, analytic or not, the difference of the two integrals is 2πi times a whole number (its winding around the loop the paths form). The measure then takes only the values 4π²n²; it is non-zero for an analytic function with a zero between the paths, so "zero iff analytic" fails. The entropy relation has no defined right side. | Not kept; the whole number is the winding of row A10 |
| 7 | Regret operator | A derivative is stated to act as a reflection. | Not kept; the returned strand is in RW1 |
| 8 | Emotion tensor | "Ricci tensor of a potential" is not defined. | Not kept |

The closing "uncut horizon equation" of Draft B is not a statement.

## Count

Four tools from six of the twenty-three objects: A4 → CY1; A8 → JR1; A5, A13, A14 → PM1; B2 → RW1.
Six objects point to results already certified (A3, A9, A10, B3, B4, B6). Eleven are not kept.
