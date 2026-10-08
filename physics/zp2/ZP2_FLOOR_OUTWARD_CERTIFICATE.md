# ZP2 — the two floor constants with an explicit tail, and the hypothesis they need

Stage of the physics line, continuing ZP1. Exact rationals, stdlib only:
`zp2_floor_outward_certificate.py`, `test_zp2_exact.py` (6 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework theorum/28 §4 | a finite-to-infinite claim = finite packet + declared correction + tail that vanishes |
| theorum/28 §9 | a finite top without an explicit error e_n is a shadow certificate only |
| this line | ZP1-Z3, Z4 (the two constants, given there to O(ε²)); LC1-C4 / WQ1-W3 (floor ½κω of a mode); PS1-T4; RB1 (edge memory) |

ZP1 gave −1/12 and −1/360 with a remainder written O(ε²) and no bound. By §9 that is a shadow
certificate. This stage gives the packet, the correction and the bound.

## Set-up

Each mode's floor is weighted by η(ω/Ω) with η(t) = (1 − t)^m on [0, 1] and y = ΩL/πc an integer
(the cut-off lies on a cell boundary of the mode ladder). All sums are finite.

One cut (line of length L): E = (πκc/2L)·F_m(y), F_m(y) = Σ_{n≤y} n(1 − n/y)^m.

Three cuts (plates d apart): E/area = (κc/2π)(π/d)³ y³·S_m(y), S_m(y) = ½g(0) + Σ_{1≤n≤y} g(n/y),
g(u) = ∫_u¹ t²(1 − t)^m dt.

## Statements

**T1 (one cut, exact).**

| m | F_m(y) |
|---|---|
| 0 (sharp) | y²/2 + y/2 |
| 1 (closes with a kink) | y²/6 − 1/6 |
| 2 | y²/12 − 1/12, no remainder |
| 3 | y²/20 − 1/12 + 1/(30y²) |
| m ≥ 2 (checked to 7) | y²/((m+1)(m+2)) − 1/12 + r_m(y), \|r_m(y)\| ≤ e_m/y² |

with e_2 = 0, e_3 = 1/30, e_4 = 1/20, e_5 = 3/28, e_6 = 31/168, e_7 = 25/72.
Each row is a polynomial identity in y, fixed by interpolation and checked on further points.

**T2 (the certificate with e = 0).** For m = 2 and every integer y,

  E(L) = κΩ²L/(24πc) − πκc/(24L)  exactly.

For a wall at a inside a fixed length the first term does not depend on a; the force is that of
−πκc/24L with no cut-off left in it, at every finite cut-off. (Test: the wall energy on a lattice is
reproduced with zero remainder.)

**T3 (three cuts, exact).** The coefficient of y⁻³ in S_m is

| m | 0 | 1 | 2 | 3 … 7 |
|---|---|---|---|---|
| coefficient | 0 (a y⁻¹ term −1/12 instead) | −1/120 | 0 | −1/360 |

For m ≥ 3: S_m(y) = y·I_m − 1/(360y³) + r_m(y), |r_m(y)| ≤ e_m/y⁵, e_3 = e_4 = 1/252,
e_5 = 61/5040, e_6 = 17/840, e_7 = 119/2640. Then E/area = (bulk ∝ d) − π²κc/(720d³) + remainder
bounded by (κc/2π)(π/d)³·e_m/y², and the force per area is −π²κc/(240d⁴) up to that bound.

**T4 (the hypothesis).** The constant is the same for every cut-off that closes smoothly enough —
zero slope at the edge for the line (m ≥ 2), zero curvature as well for the plates (m ≥ 3) — and for
any mixture of such cut-offs with η(0) = 1 (linearity; checked on one mixture inside its bound).
It is **not** the same for a cut-off that closes with a kink: the line gives −1/6 (twice), the plates
give −1/120 (three times) for m = 1 and no d⁻³ term at all for m = 2. A sharp cut-off gives no
constant for the line.

## Reading

ZP1 said the constant "does not depend on the cut-off at all". That is too strong. It does not depend
on the cut-off within the class that closes smoothly; an edge that closes with a kink keeps a term of
its own, of the same size as the constant, however far away the edge is. In theorum/28's words the tail
does not vanish there. It is the same kind of fact as RB1's: an end that is not closed properly is
remembered. The measured force agrees with −1/360 (general knowledge), so the walls' loss of
reflection with frequency is of the smooth kind.

## What is put in, what is not claimed

* Cut-offs are polynomial, (1 − t)^m and their mixtures, with the edge on a cell boundary. Other
  shapes and non-integer y are not treated; for those the bound of T1/T3 is not claimed.
* m ≥ 2 (line) and m ≥ 3 (plates) are checked for m up to 7, not proved for all m here.
* Mode frequencies between ideal walls are taken from ZP1. No new measured number.
