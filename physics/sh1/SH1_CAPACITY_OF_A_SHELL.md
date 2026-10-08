# SH1 — How many readings a shell holds: 2n², from least-cost readings and independence

Monty Dabas. 9 October 2026. Python 3.12. Exact linear algebra over the rationals (sympy).

BR1 gave one charge around a centre. A real atom has several. This stage asks what the line can say about that
without a rule for the kinds of matter: how many different readings one shell holds.

Sources read before building: **MC1** (least cost: no Laplacian), **QC1-T3** (a reading of a given degree splits
into contents; in two cuts ±q), **RC1-R7** (first-order bound histories of the 1/r attraction are great circles of
the 3-sphere in four dimensions), **BR1** (shells), **IN1-T1/T5** (a reading has two components; the record of two
readings), **DM1-T5** (a record ρ = AA† has n² − r·r = 4|det A|²; zero exactly for a single reading), **MM1** (memory
as an elementary symmetric function), **FR1** (three cuts).

## Results

**S1 — least-cost readings of degree l.**

```text
two cuts            2 for every l ≥ 1        the contents ±l of QC1-T3
three cuts          2l + 1                   1, 3, 5, 7, 9, …
four variables      (l + 1)²                 1, 4, 9, 16, …
```

**S2 — a shell.** Shell n holds the degrees 0 … n − 1 in three cuts: 1 + 3 + … + (2n − 1) = n². That is exactly the
number of least-cost readings of degree n − 1 on the four-dimensional sphere of RC1-R7: the shell is one degree
there.

**S3 — two components.** A reading on the cut-complex carrier has two components (IN1), so a shell holds 2n²:

```text
2 ,  8 ,  18 ,  32 .
```

**S4 — why no more.** IN1-T5's memory of two readings is the squared wedge, n₁n₂ − r₁·r₂ = 2|ψ₁∧ψ₂|²: two readings
that are the same leave no record of being two (DM1-T5). For N readings the N-fold memory of the record, e_N of
its tensor (MM1's form one step up), is (Π p_s) × the Gram determinant: zero exactly when the readings are not
independent. A record of N distinct readings exists only if they are independent, so a shell holds at most as
many as its dimension.

**S5 — the table of elements.** Its rows have lengths 2, 8, 8, 18, 18, 32: every one is a capacity 2n², and the
noble atoms sit at the sums 2, 10, 18, 36, 54, 86.

## What is put in

- Polynomial readings and least cost (MC1) as the readings of a shell; the degree range 0 … n − 1 of a shell,
  through RC1-R7's sphere; two components (IN1).
- The reading of "distinct" as "leaves an N-fold memory".

## What is not shown

- 2n² and the row lengths are the known counting (general knowledge). The line's part: the count is of
  least-cost readings, the shell is one degree on the four-dimensional sphere of its own orbits, and the limit is
  the vanishing of e_N — the same form as its other memories.
- Why each length after the first comes twice, and the order in which shells fill, are not derived: they need
  how the charges act on each other.
- That charges must be distinct readings — that their records must have N-fold memory — is the reading used here,
  not a derived law. The kinds of matter remain open.

## Claim boundary

```text
LEAST-COST READINGS: 2 ; 2l + 1 ; (l + 1)²                                           PROVED (exact, degrees to 6)
SHELL n = n² = ONE DEGREE ON THE FOUR-DIMENSIONAL SPHERE ; ×2 COMPONENTS = 2n²       PROVED
MEMORY OF TWO READINGS = SQUARED WEDGE ; N-FOLD MEMORY = GRAM DETERMINANT            PROVED (exact)
ROW LENGTHS 2, 8, 8, 18, 18, 32 ARE CAPACITIES                                       KNOWN NUMBERS ; MATCH
ORDER OF FILLING ; WHY DISTINCT                                                      NOT DERIVED
```

## Reproduce

```text
python sh1_capacity_of_a_shell.py
python -m unittest test_sh1
```

## Later note (FD1, 9 October)

e₂ and e_N are the z² and z^N coefficients of det(1 + zX); FD1 reads the line's other identities from the same
function. Nothing above is changed.
