# CT1 — The two-valued record inside the chain of turns: its coupling is the log of the turn coupling

Monty Dabas. 9 October 2026. Python 3.12. Exact rational series, mpmath, exact quaternion arithmetic.

DU1: a chain of two-valued marks has rate 2k* ≈ 2e^(−2k), and doubling the cell takes ½ ln 2 from k. MG1-M4: the
chain of turns has rate −Log(I₂(κ)/I₁(κ)) → 3/(2κ), a power. The turn block has a two-valued record inside it,
its centre U → −U. This stage finds the coupling of that record.

Sources read before building: **YM-3/4** (face expansion exp[(κ/2)Tr U] = Σ d_j f_j χ_j, f_j = 2I_(2j+1)(κ)/κ),
**YM-15, YM-21** (decay ratio r = f_½/f₀, exactly geometric), **YM-F1** (the turn block is the unit sphere of the
EMK block), **YM-37** (its measure is given by rational moments of the sphere), **YM-39** (the centre grades the
fabric: contents carry (−1)^(2j); mixed quantities vanish), **DU1**, **FD1-F5**, **MG1**.

## Results

**C1 — the two readings of the centre.** Write c = ½Tr U. Under the centre c → −c, so the face weight e^(κc)
splits into an even part cosh(κc) and an odd part sinh(κc). With the moments of the sphere,
⟨c^(2n)⟩ = Catalan(n)/4ⁿ,

```text
f₀ = ⟨cosh κc⟩ ,      f_½ = ⟨c·sinh κc⟩          (series equal term by term, 25 terms, exact)
```

The top reading of each centre sector is the mean of that part.

**C2 — the centre coupling.** Read as a two-valued block (even reading a + b, odd reading a − b), the centre
record has

```text
tanh k_c  =  f_½ / f₀  =  I₂(κ)/I₁(κ)  =  ⟨ c·tanh(κc) ⟩   under the weight cosh(κc) .
```

By construction the chain's rate is then 2k_c*, DU1's form. The content is the next line.

**C3 — the log.**

```text
weak end      k_c = ½ ln(4κ/3) + O(1/κ)            κ = ¾ e^(2k_c)
strong end    k_c ≈ κ/4
```

The turn coupling is the exponential of the centre coupling. So the chain of turns obeys the same exponential
law as DU1 — rate ≈ 2e^(−2k_c) — and MG1's power 3/(2κ) is that law written in κ. Doubling the cell takes ½ ln 2
from k_c, which is κ → κ/2. The proton's count sits at κ = 2 × 10¹⁹ or at k_c = 22.4: one fact, two couplings.

**C4 — the self-dual centre record.** sinh 2k_c = 1 at κ = 1.886, where the rate is ln(1 + √2) = 0.8814. Below
it the description by contents is the short one, above it the description by marks.

**C5 — more dimensions.** For a face, Tr(Π s_l U_l) = (Π s_l)·Tr(Π U_l) (exact, all sixteen sign choices on a
face of four links). So the centre marks form a two-valued gauge record whose coupling on each face is κ·c_p, set
by the cosets and different from face to face.

## What this gives and what it does not

```text
the centre record of the chain            coupling k_c = artanh(I₂/I₁) ;  k_c = ½ ln(4κ/3) + …           PROVED
DU1's "turn chain gives a power"          corrected: the same exponential law, in the centre coupling
which coupling is the natural number      NOT DECIDED by the line: κ = 10¹⁹ and k_c = 22 are the same chain
four dimensions                           the centre record has face-by-face couplings κ·c_p ; no single
                                          doubling map ; its dual is not built                           OPEN
```

## What is not shown

- C2's last form is by construction; C1 and C3 carry the content.
- The value 1.886 is a property of one face of the chain. The known change of regime of the colour record in
  four dimensions lies near κ = 2.2 (recalled, not re-read); it is a different quantity and no link is claimed.
- Nothing here derives a mass count.

## Claim boundary

```text
f₀ = ⟨cosh κc⟩ , f_½ = ⟨c sinh κc⟩ ; tanh k_c = I₂/I₁ = ⟨c tanh κc⟩                                   PROVED
k_c = ½ ln(4κ/3) + O(1/κ) ; DOUBLING: k_c → k_c − ½ ln 2 ⇔ κ → κ/2                                    PROVED (numerically, to κ = 10⁶)
SELF-DUAL CENTRE RECORD AT κ = 1.886                                                                  PROVED
THE CENTRE RECORD IN FOUR DIMENSIONS: ONE DOUBLING MAP                                                NOT BUILT
```

## Reproduce

```text
python ct1_centre_record_of_the_turn_chain.py
python -m unittest test_ct1
```

## Later note (CG1, 9 October)

RKF theorum/41 (cut-graded generator; certificate re-run, pass) is the common source: the chain block is its flow
at the dual coupling, the doubling is its join × cut identity, the centre split is its even and odd channels, and
the three sectors are the three signs of the square of one generator (CG1). Nothing above is changed.

## Later note (CZ1, 9 October)

On a closed surface the centre record closes at fixed cosets, and the weight of a twist is an exponential of the
turn coupling: for a cube 4.27·κ·e^(−(6 − 3√3)κ), against the power 3/(4κ) for one open face. A twist weight of
7.7 × 10⁻²⁰ is κ = 61.7. This is not a mass count. Nothing above is changed.

## Later note (RW1, 9 October)

RW1 (`tools/rw1`) reads I₂/I₁ as the mean returned count of two strands and proves, for every κ > 0,
κ/2 < sinh 2k_c < 2κ/3, with sinh 2k_c = κ/(1 + spread/mean). C3's two ends are the two sides of this
inequality; the weak end, checked numerically above, now has a written proof
(|k_c − ½ ln(4κ/3)| < 1/(2(κ−2)) for κ ≥ 3). C4's point lies in 3/2 < κ < 2. Nothing above is changed.
