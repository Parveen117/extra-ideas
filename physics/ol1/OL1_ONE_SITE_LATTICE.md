# OL1 — The one-site lattice: the count on compact links, at strong coupling and at every coupling past 10⁷

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals; parameters of the trial
readings are found in floats and then used as exact rationals.

GC1, and SG1 independently with the same numbers, gave the free core — flat turns — a certified gap. YM99 ended
on the other side of the same question: at weak coupling a rate uniform in the volume "must come from the
compactness of the link groups and the non-quadratic part", and YM98 §5 said why no local comparison sees it at
the order √θ. This stage takes the count of GC1 to compact links, on the smallest lattice there is: the
periodic lattice of one site. It has three links and three plaquettes, every link is a toron, and nothing
is cut off — it is the Wilson operator of the YM line at volume one. The gap is certified on a strong window
and at **every** coupling from 10⁷ to infinity, where it grows like θ^(1/3).

Sources read before building: **SG1** (in full; merged here), **GC1** (G1 the two sign certificates and the
zero-average reading; G2 the count), **HL1** (H7 the commutator of two blocks), **CM1**; Publications **YM98**
(§5 the obstruction at order √θ; the normalization of A_θ) and **YM99** (T1 torons, T3 the flat Hessian, the
verdict); RKF **theorum/28** (an outward certificate needs its own floor).

## The operator

A link is a unit block U = (cos ψ, sin ψ·m), m a unit direction, 0 ≤ ψ ≤ π: a point of the unit 3-sphere. Its
three cut components are u = sin ψ·m. For two links, with W = ½·trace of U_i U_j U_i⁻¹ U_j⁻¹,

```text
O1      1 − W_ij = 2·|u_i × u_j|²                                           (HL1-H7 ; 40 exact pairs)

H = − Σ_i Δ_i + 2θ·Σ_{i<j} |u_i × u_j|²        on three unit 3-spheres,
```

Δ the Laplacian of the unit sphere. It is the core of TC1 with each flat turn rolled up into a sphere. The
Wilson operator of the YM line is A = H/4 − 3θ_YM with θ = 4θ_YM, so gap(A) = gap(H)/4.

Two symmetries carry the count. The **row turning**: one conjugation of all links, which turns all u_i
together (the gauge condition of the one site). The **centre flips**: U_i → −U_i for one link, which changes
no plaquette. Readings unchanged by the row turning and by the three flips are the sector of the vacuum; the
gap below is the gap in this sector.

## Results

**O2 — the layer on the sphere.** For one link with another held fixed, the weight is μ·t with
t = x₁² + x₂² (the part of u across the other link). On functions of t the Laplacian of the sphere is
4[t(1 − t)f″ + (1 − 2t)f′], and the reading Exp(−κt) has local rate 4κ + 4κ²t² when μ = 8κ + 4κ². So

```text
− Δ + μ·(x₁² + x₂²)   ≥   2·( √(4 + μ) − 2 ) :
```

μ/2 for small μ and 2√μ − 4 for large μ, the flat layer's 2√μ less a constant. (The true level is 2√μ − 2 at
large μ; tests.)

**O3 — the comparison.** Give half of each link's Laplacian to its two pairs, as in CR1-R4:

```text
H   ≥   k₁ + k₂ + k₃ ,        k = −½·Δ + √(4 + 4θ·sin²ψ) − 2 ,
```

one operator for each link, a function of its class angle only. Its potential has two wells, at U = 1 and at
U = −1; the centre flip exchanges them.

**O4 — the count.** Split each link into its average over its own conjugations (a class reading) and the
rest. As in GC1-G2, a reading unchanged by the row turning has no part with exactly one link outside its
average. On class readings even under the flip, k is −½u″ − ½u + V·u on [0, π/2] with a free end at π/2
(u = sin ψ × the reading); let e₀, e₁ be floors of its first two levels, and p₀ a floor of k on readings with
zero average (there the sphere adds at least 1/sin²ψ). Then at most **one** reading of the sector lies under

```text
z = min( 2·e₀ + e₁ , e₀ + 2·p₀ ) .
```

**O5 — the strong side: 0 ≤ θ ≤ 15/4.** The reading Exp(−κ·sin²ψ) has local rate
3κ − 4κy − 2κ²y(1 − y) + V with y = sin²ψ, and sin ψ·Exp(−κ·sin²ψ), for zero average,
3/2 + 5κ − 6κy − 2κ²y(1 − y) + V; their least values over y (1000 cells, exact) are e₀ and p₀. e₁ ≥ 4, the
second level without the potential. From above, the reading 1 + b·cos 2ψ on every link has mean
12b²/(2 − 2b + b²) + 4θ·S², S = (12 − 16b + 7b²)/(8(2 − 2b + b²)). The floors only rise with θ and the mean of
a fixed reading too, so the line at θ_k against the reading at θ_k + 1/8 covers the step:

```text
θ          e₀        p₀        line z     lowest rate at θ + 1/8 ≤     gap on the step ≥
0          0         1.5       3          0.2798                       2.7202
0.75       0.4793    2.0351    4.5496     1.8961                       2.6536
1.5        0.8686    2.4612    5.7373     3.4045                       2.3329
2.25       1.1980    2.8341    6.3961     4.8074                       1.5888
3          1.4975    3.1582    6.9950     6.1105                       0.8845
3.625      1.7177    3.4091    7.4354     7.1257                       0.3098
```

The gap is positive on the whole window (30 steps; least margin 0.3098). In the units of the YM line the
window is θ_YM ≤ 15/16.

**O6 — the weak side: every θ ≥ 10⁷.** √(4 + 4θ·sin²ψ) − 2 ≥ 2√θ·sin ψ − 2, and on [0, π/2]
sin ψ ≥ σ·min(ψ, a) with σ = sin a / a. With ψ = ℓ·s the class operator is then at least

```text
unit × ( −u″ + min(s, s_a)·u ) − 5/2 ,        unit = (2θσ²)^(1/3) ,      s_a³ = 4√θ·a²·sin a :
```

the half line of GC1-G1, cut flat at s_a. Take a = 2/5; from θ = 10⁷ on, s_a ≥ 9, past the reach of both
sign certificates of GC1. The same solution v serves: it is continued past s_a by a constant, where the
potential s_a is above λ and the slope can only drop. So the first two levels are at least unit × 2.3381 − 5/2
and unit × 4.0879 − 5/2. For zero average the reading s²·Exp(−k·s^(3/2)) of GC1 is kept up to s = 5 and
continued by cosh(6/5·(b − s)) to the free end b: its rate there is at least 5 − 36/25, and the join bends
downward. That gives unit × 3.245 − 5/2. Hence

```text
z   ≥   8.7641·(2θσ²)^(1/3) − 15/2   ≥   10.8464·θ^(1/3) − 15/2 .
```

From above, the reading cos^(2n)ψ on every link has mean 36n²/(4n − 1) + (9θ/4)/(n + 1)² (exact, by the
integrals of the sphere); with n the whole part of (θ/2)^(1/3) it is at most (27/2)·(θ/2)^(1/3) + 3. So

```text
for every θ ≥ 10⁷ :        gap  ≥  0.1312·θ^(1/3) − 21/2  >  0 ,        gap ≥ 0.08·θ^(1/3) ,
as θ → ∞ :                 gap / θ^(1/3)  ≥  8.7641·2^(1/3) − (27/2)·2^(−1/3)  =  0.3271 .
```

```text
θ          n          line z ≥          lowest rate ≤       gap/θ^(1/3) ≥
10⁷        170        2329.2996         2301.7209           0.1280
10⁹        793        10838.9633        10708.2094          0.1307
10¹²       7937       108457.1330       107142.8228         0.1314
10¹⁸       793700     10846455.8016     10714950.3508       0.1315
```

In the units of the YM line: for every θ_YM ≥ 2.5·10⁶, gap(A) ≥ 0.052·θ_YM^(1/3) − 2.625.

## What it says

```text
strong side      gap of order 1               the Laplacians of the links
weak side        gap of order θ^(1/3)         the quartic core of the torons, on compact links
quadratic rates  of order √θ (YM99-T4)        zero on the torons
```

At one site there are only torons, and their rate at weak coupling is not zero: it is θ^(1/3) times the
numbers of the free core, proved here from below with the constant 0.327 in the limit (0.669 = 2 × 0.3342
would follow from a reading with the correlations of CR1's ladder; the product reading used here loses half).
This is the order YM98 §5 named as invisible to comparisons of single links at order √θ; the count sees it
because the comparison keeps the layer — each link's zero-point rate across the others — and the row
condition removes the one cheap reading.

## Inputs taken from outside

```text
as in GC1: the kinetic part of one link in its class angle and its sphere of directions ; a reading with
    zero average on that sphere has mean of |gradient|² at least 2 × mean of its square ; levels of a
    self-adjoint operator and the count under one condition
the levels of −u″ on an interval (for e₁ ≥ 4) ; sin ψ ≤ ψ ; sin ψ/ψ falling on [0, π/2]
a positive reading with a downward bend at a join still bounds the form from below (one integration by parts
    on each side of the join)
```

## What the certificate checks

```text
O1    1 − W = 2|u × v|² on 40 rational pairs of unit blocks
O2    the Laplacian on functions of t (t^m, m ≤ 8, against the harmonic route) ; the local rate ; 4κ = 2(√(4+μ) − 2)
O3    the potential of one link
O5    the two readings' means by the integrals of the sphere ; 30 steps of the strong window ; θ = 0
O6    sin a from below ; s_a³ ≥ 729 ; the two sign certificates of GC1 (run again) ; the join of the
      zero-average reading ; the two coefficients ; the algebra of the power reading (n ≤ 60) ; the gap at
      10⁷ and beyond ; four spot values with the exact power reading ; the limit
```

Tests take other routes in floats: the plaquette with 2 × 2 matrices; the Laplacian of the sphere by
differences; the layer's true level against its floor; both local rates by differences; the readings by
quadrature; the two levels of the comparison at θ = 10⁷ by a direct solve against the floors by the half line.

## What is not shown

- **Between 15/4 and 10⁷ nothing is certified.** In floats the same comparison with its true levels and the
  best product reading holds up to about θ = 7 and again from about 2·10³; in between it does not decide.
  Using the true level of the layer (O2 loses about 2 at large μ) narrows the undecided stretch to roughly
  12 … 300, and there the comparison may simply be too weak: the lowest rate needs a reading with
  correlations, or the comparison a term between links. A certified grid from 2·10³ to 10⁷ needs only the
  shooting of GC1 with this potential; not built.
- **Only the sector of the vacuum.** Readings odd under a centre flip are the flux sectors. Their lowest rates
  are expected to approach the vacuum's at weak coupling (the two wells of each link), so the gap of the whole
  operator is a different, smaller number; not treated.
- **One site.** There are no modes that are not constant; YM99's transverse rates and their coupling to the
  torons begin at two sites. No statement about a larger lattice, a volume limit or a continuum; no mass gap.
- SU(2) only (the links are unit blocks).
- In older terms this is the small-volume picture of the theory — torons, their quartic valley, rates of
  order g^(2/3) — which is known from the physics literature (general knowledge, not re-read). The line's
  part: the exact weight on compact links, the closed floor of the layer on the sphere, the count in the
  sector of the vacuum, and a certificate that holds at every coupling past 10⁷ with the same two sign
  certificates as the free core.

## Claim boundary

```text
1 − W = 2|u_i × u_j|² ; −Δ + μ(x₁² + x₂²) ≥ 2(√(4 + μ) − 2) ON THE 3-SPHERE              PROVED
H ≥ k₁ + k₂ + k₃ ; AT MOST ONE READING OF THE VACUUM'S SECTOR UNDER z                    PROVED, given the inputs
GAP > 0 FOR 0 ≤ θ ≤ 15/4 (θ_YM ≤ 15/16)                                                  PROVED, given the inputs
GAP ≥ 0.1312·θ^(1/3) − 21/2 > 0 FOR EVERY θ ≥ 10⁷ ; LIMIT OF GAP/θ^(1/3) AT LEAST 0.3271   PROVED, given the inputs
15/4 < θ < 10⁷ ; THE FLUX SECTORS                                                       NOT OBTAINED
TWO OR MORE SITES ; ANY VOLUME OR CONTINUUM LIMIT ; A MASS GAP                           OPEN
```

## Reproduce

```text
python ol1_one_site_lattice.py          (about five seconds)
python -m unittest test_ol1
```
