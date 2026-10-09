# TV1 — The valley of d turns: the quadratic layer closes for three turns and is logarithmic for four

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rationals, exact polynomial identities,
exact whole-number path counts.

The owner's question: which space is the four of four dimensions, and the method "four to three by removing a
block". HL1 found the commutator of two unit blocks, 1 − 2|u × v|²: the count part of each block is removed and
three cut components are left. This stage takes d such turns at once — the constant modes (torons) of a torus
of d directions, YM99's flat directions — and asks what the quadratic layer leaves on their valley, as a
function of d.

Sources read before building: **HL1** (H7), **TW1**, **AG1**; **YM99** (torons: 3(N² − 1) flat directions of
the quadratic layer; "the non-quadratic, compact part" as the open target); **YM100** only as recorded in the
project notes (nine flat directions, three on the valley, lift zero on the centre orbit — its files are not in
this clone); **MG1** (verdict); tools **RW1** (returned count of two strands); **CZ1**.

## The record

d unit blocks U_i = (cos α_i ; **u**_i), |**u**_i| = s_i = sin α_i. Every pair carries the commutator weight,
a function of X_ij = |**u**_i × **u**_j|². The valley is X = 0: all vector parts on one axis **n**, the angles
α_i free. The centre points are s_i = 0 for every i: 2^d of them. Constant link fields on a torus of any size
have this weight (every plaquette of a constant field is a commutator).

## Results

**V1 — the exact split.** About any axis **n**, write **u**_i = p_i·**n** + **t**_i with **t**_i across the
axis. Then

```text
X_ij = | p_i·t_j − p_j·t_i |²  +  | t_i × t_j |² .
```

The first term is quadratic in the transverse parts: the layer. The second is quartic: the core. With no part
along the axis only the core is left. Near the valley, **u**_i = s_i(**n** + δ_i):

```text
X_ij = s_i²·s_j²·( |δ_i − δ_j|² + |δ_i × δ_j|² )        exactly, not to second order .
```

At a centre point the layer is zero — this is YM100's "the quadratic layer is blind there", in one line.
(200 exact rational pairs and 100 triples near the valley, four axes.)

**V2 — the layer on the valley is a graph form, and its determinant is a count of trees.** The layer is
Σ_(i<j) s_i²s_j²·|δ_i − δ_j|²: the d turns are the points of a complete graph, the pair (i, j) weighs
x_i·x_j with x_i = s_i². One common shift of all δ_i (turning the axis) is free. For the rest,

```text
reduced determinant  =  ( Π x_i ) · ( Σ x_i )^(d−2)        (polynomial identity, d = 2 … 6) ,
```

which at x = 1 counts the trees on d points: 1, 3, 16, 125, 1296. If no s_i is zero the common shift is the
only free direction; a turn sitting at a centre point frees more.

Each block's own measure is s_i²·dα_i (the sphere of three dimensions). The layer has two transverse
components, so its Gaussian weight is the inverse of the determinant; the product Π s_i² cancels:

```text
weight left on the valley by the quadratic layer  =  ( sin²α₁ + … + sin²α_d )^−(d−2) .
```

It is singular only where every turn is at a centre point.

**V3 — three turns close, four give a logarithm.** Near a centre point the weight is |α|^(−2(d−2)) against a
volume |α|^(d−1)·d|α|: radial exponent 3 − d.

```text
d = 2      weight 1                                  flat
d = 3      1/(sin²α₁ + sin²α₂ + sin²α₃)              integrable : sin²α ≥ (2α/π)² gives a mean below 11/4
d = 4      1/(sin²α₁ + … + sin²α₄)²                  not integrable : sin²α ≤ α² gives ∫ d⁴α/|α|⁴, a logarithm
d ≥ 5      a power
```

The same in whole numbers. With y = (1/d)·Σ cos 2α_i the sum of squares is (d/2)(1 − y), and the mean of yⁿ
over the angles is c_n(d)/(2d)ⁿ, where c_n(d) counts the closed paths of n steps on the lattice of d
directions (6, 90, 1860 for three; 8, 168, 5120 for four; three routes). So

```text
three turns     mean weight = (2/3)·Σ_n c_n(3)/6ⁿ                  the returns of one walk in three directions
four turns      mean weight = (1/4)·Σ_n (n + 1)·c_n(4)/8ⁿ          the meetings of two walks in four directions
```

and c_n(d) = n!·[xⁿ](Σ_j x^(2j)/(j!)²)^d : d returned counts of two strands (RW1 at ν = 0) sharing the steps.

```text
steps 2m           16        32        64        128       256       512       1024
three turns        0.9058    0.9349    0.9566    0.9723    0.9835    0.9915    0.9972       each doubling adds less
four turns         0.5191    0.5863    0.6550    0.7244    0.7942    0.8643    0.9344       each doubling adds more
four, per doubling           0.0672    0.0687    0.0694    0.0698    0.0700    0.0701
```

For three turns the partial sums rise to a finite mean. For four they rise by a fixed amount at every doubling
— a logarithm, with no end. For five each doubling adds about √2 times the one before.

## What this says about four and three

The quadratic layer of the constant modes is closed in itself for three turns: the centre points, where it is
blind, carry no weight of their own. For four turns it is not: the weight piles up at the centre points
logarithmically, and the quartic core of V1 — the part with no component along the axis — can no longer be
left out. Four is the number of turns at which the non-quadratic core becomes necessary for the record itself.

```text
three directions     one walk leaves and its returns are finite       the layer closes
four directions      two walks keep meeting, logarithmically          the core must enter
```

This is a statement about the constant modes and their quadratic layer. It places the logarithm; it does not
compute what the core does with it.

## What is not shown

- No mass gap. The Hamiltonian of the constant modes (Casimir on each block plus the commutator weight, with
  the centre grading) is not built; neither is the record of four turns with its core.
- The weight of V2–V3 is that of the quadratic layer only. The claim "the full record of four turns carries a
  logarithm of the coupling at the weak end" would follow if the core cuts the divergence off at its own
  scale; that is not certified.
- Only two colours (unit blocks). Sites of a lattice other than through their constant modes are not treated.
- In older terms (general knowledge, not re-read): the mean for three turns is twice Watson's integral for the
  cubic lattice, 1.010924…, which is 2(18 + 12√2 − 10√3 − 7√6)·S(√6)² with AG1's source at heat time √6 (a
  floating cross-check in the tests, partial sum plus estimated rest, four places); the per-doubling amount for
  four turns tends to Log 2/π²; the finiteness of the return count in three directions is Pólya's; the
  convergence of flat-space integrals of this weight only from five directions on is known for two colours. The
  line's part: the exact split on the block, the tree determinant on the valley, the reading by closed paths
  and returned counts, and the place of the logarithm at four turns in the terms of YM99/YM100.

## Claim boundary

```text
X_ij = |p_i t_j − p_j t_i|² + |t_i × t_j|² ; LAYER ZERO AT CENTRE POINTS          PROVED (exact)
REDUCED DETERMINANT (Π x_i)(Σ x_i)^(d−2) ; VALLEY WEIGHT (Σ sin²α_i)^−(d−2)        PROVED (polynomial identity, d ≤ 6 ; measure by one line)
MEAN WEIGHT AS A SUM OVER CLOSED PATHS ; c_n BY THREE ROUTES                      PROVED (exact counts)
THREE TURNS FINITE ; FOUR TURNS LOGARITHMIC ; FIVE A POWER                        PROVED by the two inequalities ; partial sums to 1024 steps
VALUES 1.0109… AND Log 2/π²                                                       CLASSICAL ; floating cross-check only
WHAT THE CORE DOES AT FOUR TURNS ; A MASS GAP                                     NOT SHOWN ; OPEN
```

## Reproduce

```text
python tv1_valley_of_turns.py
python -m unittest test_tv1
```

## Later note (TC1, 9 October)

[TC1](../tc1/TC1_CORE_OF_TURNS.md) takes the core. The weight of all pairs is the second symmetric function of
the three squared singular values of the block matrix — the one law X = R − D, zero on the valley — for any
number of directions; the number of directions is the exponent (d − 4)/2 of the measure, zero at four. For four
turns every doubling of the largest number adds (Log 2)/16 (exact bounds), and the logarithm is cut at the
core's own scale |C| ~ κ^(−1/4). Nothing above is changed.
