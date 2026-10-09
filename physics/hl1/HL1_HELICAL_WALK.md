# HL1 — The walk is a helical orbit: UGD digits of the pair, the record's own clock, the three-sector ledger of the turn block, and the toron core on the block

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integer series (order 400), exact
rationals, directed rational enclosures (tools/pm1, AG1).

The owner's pointer: how the memory of every step adds up over the couplings is already in the repositories —
UGD number, EMK geometry, the helical (double helical) torus. And his question: which space is the "four" of
four dimensions. This stage reads those sources and writes TW1 in their terms.

Sources read before building: Publications **emk-ugd-algebra** (README in full: EMK-1 determinant identity,
EMK-2 T2 Cayley law and multiplicative seam channel, UGD-1 digit = phase exponent, scale index, seam charge);
Publications **emk-recognition-geometry** (README in full, the headers of the EMK-G3 and UGD-G1 certificates:
helical lift, n circuits displace by (nα, nβ, nq), visible return with lifted non-return, three closure rules,
lossy projection); RKF **F00-G** §5, §7 (Exp(A(y)) = (1+y)/(1−y); Log(xy) = Log x + Log y); RKF **theorum/43**
§1–6 (bilateral strands, jet parity; "a later helical adapter"); **YM99** (torons; next gate "silenced by the
compact wrap-around"); **TW1**, **AG1**, **CZ1**, **CT1**. Nothing is imported from them; each statement used
is re-derived here on the record.

## The reading

```text
UGD / EMK                                   on the record (AG1, TW1)
block x + yK, K² = 1 (EMK-1)                K₊ + K₋·K : source and flip of the turn block ; K = the centre mark
Cayley law, Log adds (EMK-2 T2, F00-G 7.1)  TW1-W1, the step
deck transformation of the lift (EMK-G3)    one diagonal step = multiplication by 1 + ι
digit (phase, scale index, …) (UGD-1)       every whole cut-complex number = ι^φ·(1+ι)ⁿ·β
sheet                                       the contents of one scale index ; the lost part is the sheet before 0
scale sector (closes by equality)           2ⁿ ; the clock M(S, L)
lossy projection (UGD-G1 T5)                a reading in powers of the heat time u
```

## Results

**H1 — TW1's step is the UGD multiplicative law.** Write the turn block as the block K₊ + K₋·K. Then

```text
4·(a′K₊′ + b′K₋′·K) = (a − b·K)·(K₊ + K₋·K) ,          K² = 1 ,
16·(S′K₊′² − F′K₋′²) = (S − F)·(K₊² − K₋²)              the seam determinant channel multiplies ,
4·(a′K₊′ ± b′K₋′) = (a ∓ b)·(K₊ ± K₋)                   each seam channel steps by itself .
```

In the Cayley coordinate y = (K-part)/(1-part) the product is y₁ ⊕ y₂ = (y₁ + y₂)/(1 + y₁y₂), which is
multiplication in the native chart x = (1 + y)/(1 − y) (EMK-2 T2; F00-G Theorems 5.1, 7.1). So

```text
x(ρ′r′) = x(ρ)/x(r) ,        Log x(ρ′r′) = Log x(ρ) − Log x(r) .
```

This is how the memory of a step adds: in the Log chart, by subtraction of the pair's. At u = 1/4:
x(r) = 1.1892, x(ρ) = 2.6811, x(ρ′r′) = 2.2546. (x(ρ) is centre-even over centre-odd, the inverse of TW1's ε.)

**H2 — the diagonal step is a deck transformation; digits and sheets.** Every whole cut-complex number is,
in exactly one way,

```text
α = ι^φ · (1 + ι)ⁿ · β ,        φ mod 4 ,   n ≥ 0 ,   β primary (odd norm, β ≡ 1 mod (1+ι)³) .
```

Under multiplication the phases add mod 4, the scale indices add, the odd parts multiply (1680 numbers; 400
pairs). One step of the walk is n → n + 1 with the norm doubled. Two steps are a doubling and a quarter turn,
(1 + ι)² = 2ι. After eight steps the phase is back and the scale is not: (1 + ι)⁸ = 16 — EMK-G3's visible
return with lifted non-return; no earlier power is a positive whole number.

The contents of one scale index are a sheet, and the sheets are the lost parts of AG1's walk:

```text
sheet n = L_(n+1) ,        S = 1 + L₁ + L₂ + L₃ + … ,        F = 1 − L₁ + L₂ + L₃ + … .
```

The source is the origin and every later lost part, each once. The cut reverses sheet 0 and nothing else. The
lost part itself, L₀, is the sheet before 0: the half lattice is (1 + ι)⁻¹ times the odd norms, which no whole
content reaches; L₀² = 4·L₁·S′. (Lattice counts and series, order 400.)

**H3 — the record's own clock.** AG1 has S′ = (S + F)/2 and L′ = (S − F)/2. One step of the mean takes
(x + y, x − y) to (x, √(x² − y²)), so

```text
M(S′, F′) = M(S, F) = 1 ,        M(S′, L′) = ½·M(S, L) .
```

The walk has one number that never changes and one that halves at every step. The second is the record's
reading of its own heat time: at u = 1/4, 1, 3 the mean M(S, L) is 4, 1, 1/3 — u·M(S, L) = 1 to 100 digits.
The halving is proved here; that the constant is exactly 1 is the mirror law (classical, pinned as in AG1). On
the seam F = L and the two means are equal.

**H4 — three sectors in every coefficient of the turn block.** Write N = 2ⁿ·m, m odd, σ(m) the sum of the
divisors of m.

```text
K₋/b = Σ_N 2ⁿ·σ(m)·q^N ,        K₊/a = Σ_N (−1)^(N−1)·2ⁿ·σ(m)·q^N ,        L_n² = 16·Σ_(m odd) σ(m)·q^(2ⁿm) .
```

Every coefficient is mark × scale × count. The count is the block of four: 16·σ(m) is twice the number of
whole quaternions of norm m (lattice count, m ≤ 35), and L₀² = S² − F² is twice the odd-norm part of that
lattice. Scale and count are free of each other and of the mark. The mark is not free: it is − exactly off
sheet 0 (as in H2). A reading that drops the scale sector (2ⁿ → 1) is a different series, first at q².
(Series to order 400. The count formula is Jacobi's; it is checked, not re-proved.)

So TW1-W2 is a ledger over sheets: sheet n enters with its scale 2ⁿ, its count L_n², and the mark.

**H5 — what a power of u cannot see.** TW1 found lost/seen Λ → 1 − 4u/π at the weak end. The residue of that
reading is certified:

```text
u        (1 − Λ)·π/(4u) − 1         (8 − 4π/u)·Exp(−π/u)       below
1/4      −1.4740·10⁻⁴               −1.4739·10⁻⁴               u⁵
1/8      −1.1253·10⁻⁹               −1.1253·10⁻⁹               u⁹
1/16     −2.8555·10⁻²⁰              −2.8555·10⁻²⁰              u¹⁵
1/64     −3.8103·10⁻⁸⁵              −3.8103·10⁻⁸⁵              u⁴⁰
```

It is not zero, it is below high powers of u, and it is the first sheet term to a thousandth. The expansion of
the distance from the diagonal in powers of u ends at its first term; everything after it is a sheet term. A
reading in powers of u is UGD-G1's lossy projection: it keeps the scale sector, drops the sheets, and reports
the record closed on 4u/π. TW1's flip is made of sheet terms only. This is YM99's "silenced by the compact
wrap-around", exact for one turn.

**H6 — the two strands on the seam.** The forward walk (u → 2u) and the backward walk (u → u/2) are exchanged
by the mirror u → 1/u; TW1 used the first for positivity and the second for the weak end. They meet on the seam
u = 1. There, by theorum/43's jet parity for a mirror-even reading, the odd jet vanishes:

```text
8π·K₊ = a ,        1 − Λ = 4/(π·S²)        at u = 1 (enclosures, 100 digits ; off the seam they separate) .
```

theorum/43 asks that a helical adapter be "defined from physical data and proved to intertwine the strands".
Here the data are the deck transformation of H2 and the mirror; the intertwining is the mirror law, which this
line holds as classical — on the seam it is certified as agreement.

**H7 — the toron core on the block.** For two unit blocks U = (u₀; **u**), V = (v₀; **v**):

```text
scalar part of U·V·U⁻¹·V⁻¹  =  1 − 2·|u × v|² .
```

The commutator sees only the three cut components; the scalar part of each block is removed. It is blind to
both centre marks, so for three turns the eight centre points carry the same weight, and it is 1 exactly on a
common axis (the valley of YM99/YM100). 200 exact rational pairs.

```text
one turn      at its two centre points the content with m marks reads m and (−1)^(m−1)·m : the next jet of the
              abelian record. That is TW1's K₊, K₋ — the core of one turn is solved at every coupling.
two turns     ⟨χ_m(commutator)⟩ = 1/m for m ≤ 12, from ⟨|u × v|²ᵏ⟩ = 3/8, 5/24, 35/256, … (rational sphere
              moments, two routes). This is TW1's torus record, the one whose step inequality is reversed;
              U = i, V = j give commutator −1: the twist costs nothing there.
three turns   ⟨c₁₂⟩ = 1/4 ,  ⟨c₁₂c₂₃⟩ = 1/8 ,  ⟨c₁₂c₂₃c₃₁⟩ = 5/72 (two routes).
              The closed triple is above the open pair times a face (1/32) and above three free faces (1/64).
```

## Which four

In this record the "four" is the block: one count and three cuts, the whole quaternions. The lost part
squared counts them (H4). The turn block's closed surface is the unit sphere of that block — three dimensions
— and TW1-W2 writes it as one turn times the sheets of the block of four. The commutator of two blocks keeps
the three cut components and drops the count (H7). So "four to three by removing a block" is exact here in two
places. That the four directions of the lattice in the Yang–Mills problem are this four is not shown.

## What it gives and what it does not

```text
how the memory of the steps adds                 in the Log chart of the block K₊ + K₋K (H1) ; over sheets with
                                                 scale 2ⁿ (H4) ; the clock M(S, L) halves at every step (H3)          PROVED
the helical structure                            deck transformation, digits, sheets, visible return at eight (H2)    PROVED
the weak-coupling expansion                      keeps the scale sector only ; the flip and the residue are sheets     PROVED for one turn
the core at the centre points, one turn          the next jet ; TW1                                                   PROVED
two turns                                        TW1's torus record                                                   PROVED (m ≤ 12)
three turns (the toron core)                     the exact weight on the block, eight equal centre points,
                                                 the first three means                                                PROVED ; the record itself NOT BUILT
the fabric in four dimensions                    NOT SHOWN
```

## What is not shown

- No mass gap and nothing about the Riemann hypothesis.
- The three-turn record is not solved: no step law, no sheet decomposition, no flux energies. What is certified
  is its weight and three numbers. The Hamiltonian of the nine constants — Casimir on each block plus
  2·Σ|u_i × u_j|², with the centre grading — is the next stage; on the block it is a matrix of rational sphere
  moments.
- Old terms: H1 is the addition law of the hyperbolic tangent; H2 is the factorization of the Gaussian
  integers at the prime 1 + ι; H3 is the halving of the complementary mean under Landen's step; H4 is Jacobi's
  count of sums of four squares; H7's identity is elementary. The line's part: that these are one ledger —
  block, deck step, sheets, scale, clock — and that TW1's results are its entries; H5 as a certified residue;
  the two- and three-turn means on the block.
- u·M(S, L) = 1, the first sheet term in H5 and the seam values in H6 are agreement with the mirror description
  (classical, pinned). The halving in H3 and the non-vanishing and size of the residue in H5 do not use it.
- UGD-G1's sectors are pairwise independent on its state space; on this record the mark is tied to the sheet
  index (H4). The seam charge of UGD-1 and the sheet cocycle of UGD-G1 T2 are not used.
- "Double helical" is read as theorum/43's two strands (H6). A geometric double helix is not constructed.

## Claim boundary

```text
BLOCK STEP ; DETERMINANT CHANNEL ; CAYLEY / LOG LAW                         PROVED (series ; exact rationals)
DIGITS ι^φ(1+ι)ⁿβ ; SHEETS = LOST PARTS ; S = 1 + ΣL_n ; (1+ι)⁸ = 16         PROVED (box of 1680 ; lattice counts ; series)
M(S′, L′) = ½M(S, L) ; M(S, F) = 1                                          PROVED ; enclosures at three couplings
u·M(S, L) = 1                                                               AGREEMENT to 100 digits ; equality classical
K₋/b, K₊/a AS mark × scale × count ; L₀² ON THE BLOCK OF FOUR                PROVED to order 400 (count formula classical)
RESIDUE OF THE POWER READING ≠ 0, BELOW u⁵ … u⁴⁰                            PROVED at four couplings ; first sheet term: agreement
SEAM: 8πK₊ = a ; 1 − Λ = 4/(πS²)                                            AGREEMENT to 100 digits
COMMUTATOR = 1 − 2|u × v|² ; CENTRE-BLIND ; ⟨χ_m⟩ = 1/m ; 1/4, 1/8, 5/72     PROVED (exact)
THREE-TURN RECORD ; FOUR DIMENSIONS ; A MASS GAP ; RH                       NOT BUILT ; NOT SHOWN ; OPEN
```

## Reproduce

```text
python hl1_helical_walk.py
python -m unittest test_hl1
```

## Later note (TV1, 9 October)

[TV1](../tv1/TV1_VALLEY_OF_TURNS.md) takes H7 to d turns: the commutator weight splits exactly into a quadratic
layer and a quartic core; on the valley the layer leaves the weight (Σ sin²α_i)^−(d−2), a sum over closed
paths; it is integrable for three turns and logarithmic for four. Nothing above is changed.
