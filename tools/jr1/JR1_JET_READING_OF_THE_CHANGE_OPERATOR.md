# JR1 — The jet reading of the change operator

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact cut-complex rational arithmetic.

Source object: the "λ-residue" of the catalog draft,

```text
Res_λ(f, z₀) = e^(−i(λ−1)z₀) · Res( f·e^(i(λ−1)z), z₀ ) ,
```

stated there for "λ-analytic" functions (C f = λ f, C = I + i·d/dz), with a residue theorem and an example.

Sources read before building: **RKF F00-E** (cut-complex field, ι² = −1, factorial polynomials E_N),
**Publications GE1-T2** (the change operator C = I + D and its kernel), **RKF theorum/46** (first visible jet).

## What is wrong in the draft, in one paragraph

C f = λ f has only the solutions A·Exp(−ι(λ−1)z). They have no pole, so on its stated domain the object is
identically zero. The draft's example (a factor e^(−z₀/m) at a simple pole) does not hold: the expression
is defined on the tail in t = z − z₀ and does not depend on z₀ at all, and at a simple pole it equals the
ordinary first jet for every λ. The draft's residue theorem omits a weight Exp(ι(λ−1)z_k) at each pole.
What remains, and is exact, is below: the expression is a reading on *tails with poles*, and there it is the
eigen-reading of the change operator.

## Setting

Tails at a point over the exact cut-complex rationals, t = z − z₀, D t^k = k·t^(k−1):

```text
f = Σ_(j=1..m) a_(−j) t^(−j) + (regular part) ,          c = ι(λ − 1) ,
reading     Rd_c(f) = Σ_(j≥1) a_(−j) · c^(j−1) / (j−1)! .
```

Every sum is finite. No contour and no limit is used.

## Results

**J1 — it is the draft's expression.** Rd_c(f) is the t^(−1) coefficient of f·E_N(c·t) for any N ≥ m−1.

**J2 — eigen-reading.** Rd_c(D f) = −c·Rd_c(f), term by term (the factorial weights are what make it so). Hence

```text
Rd_c( C f ) = λ · Rd_c(f) ,        Rd_c( p(C) f ) = p(λ) · Rd_c(f)        for every tail f .
```

The change operator acts on the reading as multiplication by λ. Checked also on the tail of
1/((z−z₀)²(z−z₁)) with the derivative computed by the quotient rule.

**J3 — one jet carries no λ.** At a simple pole Rd_c(f) = a_(−1) for every λ. The λ-dependence starts at the
second polar jet: for a_(−2)t^(−2) + a_(−1)t^(−1) the reading is a_(−1) + ι(λ−1)·a_(−2).

**J4 — recovery.** A pole of order m is recovered exactly from m readings at distinct λ, and not from fewer
(a tail of order 3 with two vanishing readings is exhibited). The reading at λ = 1 is the first jet alone.

**J5 — partner of the eigenfunction.** The solutions of C f = λ f are the multiples of Exp(−c·z). A tail with
a pole of order m is never one: (C − λ)f contains −ι·m·a_(−m)·t^(−m−1). So the eigenfunctions have no pole
and reading zero. The reading of a tail is the first jet of its product with Exp(+c·t).

**J6 — shift and resolvent.** Rd_c(f·E(b·t)) = Rd_(c+b)(f). If (C − μ)g = f then Rd_c(g) = Rd_c(f)/(λ − μ)
(a corollary of J2).

## Use

No stage uses this tool yet. What it offers:

- A pole seen through the change operator: one number per λ, on which C is multiplication by λ. For a
  resolvent tail this isolates the λ-channel (J6).
- Counting how many readings a pole needs (J4): m, the order. Compare theorum/46, where the first visible jet
  decides a quotient; here the first jet is the λ = 1 reading and each further jet needs one further λ. This
  is a comparison, not a result about theorum/46.

## What is not shown

- No statement about a sum over all poles of a function is made. The classical theorem of that kind would
  carry the weights Exp(c·z_k); it is not used and not certified here.
- No link to mock modular forms. The draft's "mock theta case" is its example, which fails (J3).
- The carrier is finite tails. Essential singularities and infinite principal parts are outside.
- This is elementary algebra of Laurent tails; that Rd_c(Df) = −c·Rd_c(f) is the statement "the first jet of
  a derivative is zero" applied to f·Exp(c·t). The line's part is the correction of the draft's domain and
  example, and J2/J4 in the framework's terms.

## Claim boundary

```text
Rd_c(C f) = λ·Rd_c(f) ON EVERY TAIL ;  p(C) ↦ p(λ)                                PROVED
SIMPLE POLE: NO λ ;  ORDER m: m READINGS, NOT FEWER                               PROVED
EIGENFUNCTIONS OF C HAVE NO POLE ; READING ZERO ON THE DRAFT'S STATED DOMAIN      PROVED
THE DRAFT'S EXAMPLE AND ITS RESIDUE THEOREM AS WRITTEN                            NOT KEPT
A SUM-OVER-POLES THEOREM ; MOCK MODULARITY                                        NOT CLAIMED
```

## Reproduce

```text
python jr1_jet_reading_of_the_change_operator.py
python -m unittest test_jr1
```
