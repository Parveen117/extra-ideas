# DS1 — One diagonal: the same structure in the framework, the Riemann line and the physics line

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact cut-complex rational arithmetic.

The owner's statement. In the drafts, and in the whole framework, the seam is the diagonal — the 45° line.
It is a local reference: from one state to the next one knows what remained and what is lost, and no
absolute reference is needed. The Riemann line was worked on this line; quantum and classical agree on it;
wherever there is an imaginary component against a real one the results have the same structure: the
infinite made finite, positivity turned into a charge, the Weil criterion, the Hermitian five-by-five matrix.

This stage matches that statement against what is certified, and proves the common part once.

Sources read before building: RKF **theorum/24** Theorem 6.1 (S = R + D; F = Z*𝔧Z = R − D, 𝔧 = P − Q),
**theorum/28** (finite signed cut form F_n = R_n − D_n; dim ker F = dim ker(1 − B); "decided by a matrix of
size at most five"; outward certificate u + e < 1), **theorum/02** (F = S − VV*, rank V ≤ 5, B = V*S⁻¹V,
k_Σ = N₊(B − 1)), **theorum/34–39** (odd sector: S − L*L ≥ 0.1709·S), **theorum/40** (W = B + A_Γ − P_½; the
zero-sum formula pinned, not rederived), physics **ONE_LAW**, **DO1**, **QD1**, **DU1**, **CT1**, **PT1**,
tools **RW1**, **PM1**, **CY1**, Publications **LAM-1/2/3**.

## Results

**D1 — the flip is the mirror pairing.** For a cut P + Q = 1 put 𝔧 = P − Q; it is a mirror (𝔧² = 1). Then

```text
S = Z*Z = R + D ,          F = Z*𝔧Z = R − D          (theorum/24) .
```

The signed form of the law is the pairing of a record with its own mirror image. For the cut that exchanges
two labels, P and Q project on the two diagonals. For a mirror of mirrors (𝔧 ⊗ 𝔧′ on a product reading) the
flip is the product of the flips: it vanishes as soon as one factor is on the diagonal.

**D2 — the diagonal is the 45° line.** One reading z = a + bι read by the cut (a seen, b lost):

```text
F/S = rad( z / z† ) ,       (F/S)² + turn( z / z† )² = 1 ,       4·R·D = S² − F² .
```

Seen = lost exactly when z/z† is a quarter turn: z on the 45° lines. There the product seen × lost is
largest. z/z† is PT1's exact turn B₂(z): *the flip is the rad-part of that turn.* The diagonal direction is
1 + ι. Its norm is the prime 2; its square is 2ι, the Gauss sum of the character mod 4 (LAM-2); its fourth
power is −4, PM1's value at 2. Its unit point is not rational (PT1-P3).

**D3 — positivity is a count, and the count is finite.** Let R be positive and D = VV* of rank r. Then

```text
n₋(R − D)  =  n₊(B − 1)  =  number of zeros of det(1 − zB) in 0 < z < 1 ,        B = V*R⁻¹V   (r × r) ,
dim ker(R − D) = dim ker(1 − B) .
```

Three routes, equal on 36 exact instances with counts 0 … 5. The count is at most r — a cut whose memory
has rank r takes at most r directions from any reading — and it does not change when the space is enlarged
by directions in which nothing is lost. F ≥ 0 exactly when no direction has lost/seen above 1: nothing past
the diagonal. On the diagonal itself (lost = seen in one direction) there is a kernel and no count; a
hundredth past it, one count.

The third route is the framework's one function (FD1) on the small matrix: *the charge is the number of zeros
of det(1 − zB) inside the diagonal z = 1.*

**D4 — the mirror form.** Points x with a mirror x → x′, weights m_x = m_x′ > 0, and

```text
W(f) = Σ_x m_x · f(x) · f(x′)† .
```

A point on the mirror's fixed line gives a square. A pair off the line gives
m·(|f(x) + f(x′)|² − |f(x) − f(x′)|²)/2: one plus and one minus. So

```text
inertia(W) = ( fixed + pairs , pairs ) ,          W ≥ 0  ⇔  every point is on the line ,
```

and the count is the number of pairs off the line. Verified for the mirror s → 1 − s† (line rad = ½) and for
J(z) = ι·z† (the 45° line) with the same result. They are one mirror in two charts: z = (1 − ι)(s − ½) carries
the first into the second, exactly — the chart is multiplication by the diagonal vector's dagger. Moving one
point off the line adds exactly one to the count.

**D5 — the seam measure.** For two strands with amplitudes a, b (factorial series of F00-E):

```text
E(ab)² = Exp(−(a − b)²) · E(a²)·E(b²) ,
E(a²)·E(b²) − E(ab)² = Σ_(i<j) (aⁱbʲ − aʲbⁱ)² / (i! j!)        a sum of squares .
```

The squared mirror overlap of the two strands is Exp(−w²), w = a − b: 1 on the diagonal a = b, below 1
off it. This is the drafts' seam measure e^(−w²). The share of returned histories (equal counts on the two
strands) lies between sinh(2ab)/(2ab)·Exp(−a² − b²) and ½[Exp(−(a−b)²) + Exp(−(a+b)²)]: Exp(−w²) up to a
power. On the diagonal the returned histories are a power-law share of all; off it an exponentially small one.

**D6 — no crossing.** For RW1's returned count put x = w²/y. Then

```text
0 < x < 1 ,          y·dx/dy = w·(2V − w)/y > 0        (by RW1-R4: V > w/2) .
```

The reading x + ι turns from the cut toward the diagonal 1 + ι, always closer, never across. RW1's bound
"spread above half the mean" *is* this monotone approach. (ν = 0: x = 0.09, 0.49, 0.84, 0.95, 0.984, 0.995 at
y = 0.1 … 10⁴.)

**D7 — the cell between the diagonals.** On a clock of 4M marks, the window of M marks (a quarter turn: one
cell of the four units) has harmonics that vanish exactly at the multiples of 4 other than the mean. The cell
does not see a unit-free harmonic — the counterpart of PM1-P1 — and Σ_k |W(k)|² = 4M²: four cells.

## The match

One law, S = R + D and F = R − D = Z*𝔧Z. One special line, F = 0. One count, the directions past it.

```text
line                                   seen R              lost D                on the diagonal                 where
one reading z = a + bι                 a²                  b²                    z/z† a quarter turn             D2, PT1
block read at angle χ                  mean²               half-difference²      χ = 45°: F is the determinant   DO1-P1, EMK-1
a reading through a cut                cos²(θ/2)           sin²(θ/2)             fair coin; the reading rule and  QD1-D1, D2
                                                                                 the fixed-observer rule agree
count of a mode                        n                   n²                    n = 1, quantum kT·ln 2          QD1-D3, D4
centre                                 1 − r_s/r           r_s/r                 r = 2r_s                         DG1
chain of two-valued marks              1/cosh² 2k          tanh² 2k              sinh 2k = 1: its own dual        DU1-U2
centre record of the chain of turns    —                   —                     sinh 2k_c = 1 ⇔ κ = 1 + spread/mean  RW1-R6
returned count of two strands          x² = (w²/y)²        1                     x → 1, never across              D6, RW1
two strands, amplitudes a, b           shared: Exp(−w²)    not shared:           a = b                            D5
                                                           1 − Exp(−w²)
three readings (mirror = transpose)    symmetric part      alternating part      symmetric: needs a reference;   CY1
                                                                                 the alternating part does not
generator and cut                      G_e (2t·G_e)        G_o (returns as       —                                theorum/41, CG1
                                                           t²[G_e, G_o])
finite to infinite                     R_n                 D_n                   λ_max(B) = 1 ; inside: u + e < 1 theorum/28
Riemann line, the count                S                   VV*, rank ≤ 5         an eigenvalue of B at 1 ;        theorum/02, D3
                                                                                 k_Σ = N₊(B − 1)
Riemann line, odd sector               S (full)            L*L                   lost/seen ≤ 0.8291 :             theorum/34–39
                                                                                 inside by 0.1709
Riemann line, the form                 B + A_Γ             P_½ (primes, weights  W = 0                            theorum/40
                                                           Λ(n)/√n, shifts ±log n)
zeros and their mirror s → 1 − s†      zeros on rad = ½:   pairs off the line:   the line rad = ½                 D4 (finite) ;
                                       squares             one minus each                                         zero-sum pinned
theta seam (mirror t → 1/t)            —                   —                     t = 1                            LAM-1/2, PM1-P5
```

What the owner's sentences are, in certified terms:

```text
"no absolute reference, the 45° line is local"     the reference-free part of three readings is exactly the
                                                    mirror-odd part (CY1-Y4) ; the flip needs only the record
                                                    and its mirror (D1)
"what remained and what is lost"                   R and D of one cut ; on the diagonal one cut's (seen, lost)
                                                    holds the whole block (DO1-P1)
"positivity turned into a charge"                  n₋(R − D) = n₊(B − 1) (theorum/02 ; D3)
"the infinite made finite"                         the count is at most the rank of the memory and is read from
                                                    an r × r matrix (D3) ; r ≤ 5 on the Riemann line (theorum/02, 28)
"Weil criterion"                                   a mirror form is non-negative iff every point is on the
                                                    line ; the count is the pairs off it (D4, finite form)
"real against imaginary; same structure"           the mirror s → 1 − s† with line rad = ½ and the mirror
                                                    J(z) = ι·z† with the 45° line are one mirror: the chart is
                                                    z = (1 − ι)(s − ½) (D4)
"quantum and classical agree there"                QD1-D2
```

## Where each line stands, in this reading

```text
returned count (RW1)                 inside for every coupling ; approaches the diagonal, never crosses     PROVED (D6)
face ladder (theorum/75, RW1-R6)     two-sided at every coupling                                             PROVED
Riemann line, odd sector             inside by 0.1709                                                        PROVED in RKF (cited, not re-run here)
Riemann line, the five-matrix        "RH outward five-matrix certificate"                                    NOT YET BUILT (theorum/28)
```

The open wall of the Riemann line, in this reading, is one statement: the five-by-five matrix stays inside
the diagonal down to the end of its path — equivalently (D3), det(1 − zB) has no zero in 0 < z < 1 there.

## A way of working this suggests (a proposal, not a result)

RW1 proved "never across" without any absolute bound: it wrote (diagonal − ratio) × (a positive series) as a
series whose every coefficient is ≥ 0, on the lattice of counts. The same question can be put to the
five-matrix: det(1 − zB) is FD1's one function on a 5 × 5 matrix, 1 at z = 0, and the claim is that it stays
positive up to the diagonal z = 1. Expanding it on the prime-memory lattice of theorum/40 (weights
Λ(n)/√n — the diagonal share n^(−½)) and looking for a sign that holds term by term is the same kind of
proof. Whether it exists there is not known.

## What is not shown

- D3 and D4 are finite. Nothing is proved about the infinite zero set, about the identification of the
  native form with the classical one (theorum/40 pins that formula), or about the endpoint. RH is OPEN and
  is not moved by this stage.
- In the row "the form", W has the shape seen − lost; that the prime part is by itself a Gram form is not
  claimed. The framework's own split into two positive forms is theorum/02's S − VV*.
- D3 is the inertia law of a bordered form and D4 is the polarization identity: known algebra. The line's
  part: the three routes with the determinant, the reading as one diagonal, and the match.
- The odd-sector numbers are quoted from RKF theorum/34–39; their certificates were not re-run here.
- D5's share bounds are those of this two-strand count; the drafts' "winding" is not derived to be a − b.

## Claim boundary

```text
F = Z*𝔧Z ; F/S = rad(z/z†) ; DIAGONAL = QUARTER TURN OF z/z† ; DIRECTION 1 + ι, NORM 2                PROVED
n₋(R − D) = n₊(B − 1) = ZEROS OF det(1 − zB) IN (0, 1) ; COUNT ≤ RANK ; SIZE-FREE                       PROVED (written ; 36 exact instances)
MIRROR FORM: INERTIA (FIXED + PAIRS, PAIRS) ; SAME FOR BOTH MIRRORS                                     PROVED (finite)
E(ab)² = Exp(−(a−b)²)E(a²)E(b²) ; GRAM DEFECT A SUM OF SQUARES ; SHARE BOUNDS                           PROVED
x = w²/y: 0 < x < 1, INCREASING                                                                         PROVED
THE MATCH (TABLE)                                                                                       each row certified where named
RH ; THE FIVE-MATRIX CERTIFICATE ; N1 / N2 / N3                                                         OPEN
```

## Reproduce

```text
python ds1_one_diagonal.py
python -m unittest test_ds1
```
