# HM1 — what a reading remembers of a far closure is the product of the returns on the way

Stage of the physics line, continuing RB1, TS1 and CO1. Exact rationals, stdlib only:
`hm1_horizon_memory_product.py`, `test_hm1_exact.py` (6 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework `research/recognition_return/r2` | RD2 backward recursion x_j = 1/(r_j + x_{j+1}), r_j = 1/a_j; RI1 (R2.6); RI2 criterion Σ 1/a_j; the class a_j > 0 |
| this line | RB1 (dictionary x = fall speed; horizon end remembers), TS1 (mirror sheet), CO1 (expanding space: x = hX), LT1 (tower) |

## Statements

**T1 (product law).** For one chain closed in two ways, u and v, with returns x_j(u), x_j(v) along it,

  x₀(u) − x₀(v) = (−1)ⁿ (u − v) · Π_j x_j(u) · Π_j x_j(v).

*Proof.* One cell: 1/(r + u) − 1/(r + v) = −(u − v)·[1/(r + u)]·[1/(r + v)]. Compose. ∎
On a profile itself the reading moves by Π x_j² per unit change of the closure. This is R2.6 written with
the returns instead of the continuants: what is remembered = the product of the returns passed through.

**T2 (the class ends at the horizon).** The cell between x_j and x_{j+1} is a_j = x_j/(1 − x_j x_{j+1}):
positive exactly while x_j x_{j+1} < 1. R2's positive class covers the region up to x = 1 and no further.

**T3 (toward a horizon along the tower).** The product converges to a positive number, below x₀²:

| reading at | x₀ = 1/3 (r = 9 r_s) | 1/10 | 1/100 |
|---|---|---|---|
| remembered fraction Π x_j² | 0.03066 | 2.13·10⁻⁵ | 2.2·10⁻¹⁶ |

In an expanding space (CO1, x = hX) the same chain runs from the observer's shell out to the horizon
hX = 1. The dictionary hX = √(r_s/r) maps the centre to the flat end: the central observer is at the
forgetting end, and what a shell at hX remembers of the horizon's closure is below (hX)².
The expected reversal between the two kinds of horizon does not occur — it is one chain.

**T4 (the verdict depends on how the cells are spaced).** Memory survives iff Σ (1 − x_j) converges
(equivalently R2's Σ 1/a_j < ∞). The tower converges. The spacing x_j = 1 − 1/(j + 2) also reaches the
horizon through positive cells, and its remembered fraction after n cells is exactly 1/(n + 1)² → 0.
So "the horizon end remembers" (RB1-T3) is a statement about the tower's spacing, not about every
approach to a horizon.

**T5 (past the horizon).** Continuing the profile literally to x > 1 gives negative cells — outside R2's
class by T2 — and a memory factor above 1 (3.624 for the six levels from x₀ = 1/3). Reading the same region
through the tower's mirror x → 1/x (TS1-T1) gives positive cells and the exact reciprocal factor, 0.2759.
The mirror sheet of TS1 is the continuation that stays inside the class.

## What is put in, what is not claimed

* What a "cell" is physically — hence which spacing nature uses near a horizon — is not derived. The tower
  is the line's own generation (LT1); T4 shows the conclusion of RB1 rests on that choice.
* T5 says the mirror stays in R2's class; it does not prove the literal continuation is forbidden, only that
  no theorem of R2 covers it.
* x = fall speed remains RB1's dictionary. No measured number.

## Open gate

Derive the spacing: a framework statement that fixes how many cells lie between a shell and the horizon.
Geometric approach ⇒ the horizon is remembered; harmonic ⇒ it is not.
