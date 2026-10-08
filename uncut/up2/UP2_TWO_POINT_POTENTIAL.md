# UP2 — the two-point potential: one function on both halves of the diagram

Second stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up2_two_point_potential.py`, `test_up2.py` (6 tests).

The owner's statement: from the T–V–S–P diagram, the λ diagram and the higher-order diagrams, the potential in
(S, T) comes out after the cut; before the cut it must be some generalised potential.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Recognition-Kernel-Framework theorum/42, Theorem 3.1 and §7 | typed tower L_{n+1} = ∇L_n; a higher layer repairs what a lower one is blind to |
| Publications research branch, `coherence-first-thermodynamics` CF-3 (7) | C_V = T/A, C_P = TC/det H, K_S = VC, K_T = V det H/A |
| same branch, `thermo-compass-foundations` results index | C_P/C_V = K_S/K_T |
| response-geometry RMG2 T1, RMG9 T5 | number sectors (ρ = 1 is the dual, nilpotent sector); the Schur ratio |
| extra-ideas `uncut/up1`, `physics/lt1` | degree k and silent depth 1/k; the single null reading |

## The object

For a potential Φ(x) with first layer y = ∂Φ, define on pairs of points

  G(x; x₀) = Φ(x) − Φ(x₀) − y(x₀)·(x − x₀).

Read in the variables (x, y₀) it is a function on **both halves of the diagram at once** — for the compass,
on (S, V) and (T, P) together — with no choice of which half is independent.

## Statements

**T1 (the cut surface is the rest set).** G = 0 and every derivative of G in x and in y vanishes exactly on
y = ∂Φ(x). The equations of state are not put in beside the potential: they are where the two-point
potential is at rest. The rest set has half the dimension of the diagram (one pair: a curve in the plane;
two pairs: a surface of dimension 2 in 4).

**T2 (before the cut the response is a single null reading).** On the rest set the response of G in the
doubled variables is, for one pair,

  [[h, −1], [−1, 1/h]],  h = ∂²Φ,

of determinant 0: in the compass variables ρ = 1, the boundary (dual sector) of RMG2-T1. For two pairs it has
rank 2 of 4. Reading the x-half alone gives h; reading the y-half alone gives 1/h. The cut is what turns a
null doubled response into a full one, and the two readings are inverse because the doubled response is null.

**T3 (the layers of G are the λ level and every higher level).**

  G(x₀ + ξ; x₀) = Σ_{n ≥ 2} L_n(x₀) ξⁿ / n! .

The value of Φ and its first layer are absent — G does not change when any affine function is added to Φ.
The quadratic layer is the response level (the λ diagram), the cubic and higher layers are the higher
diagrams. Choosing the base point x₀ is the cut; changing it obeys
G(x; x₀) = G(x; x₁) + G(x₁; x₀) + (y(x₁) − y(x₀))·(x − x₁).

**T4 (the two halves exchange).** For Φ = x^K/K and its conjugate Φ° = y^{K*}/K*, K* = K/(K − 1):
G_Φ(x; x₀) = G_{Φ°}(y₀; y). In (x, y): G = x^K/K + y^{K*}/K* − xy.

**T5 (weights and the two depths).** G(s·x, s^{K−1}·y) = s^K·G: the x-half has weight 1, the y-half weight
K − 1, each pair weight K. And 1/K + 1/K* = 1: the silent depth of UP1 measured from the uncut corner, and
the same depth measured from the far corner by the conjugate degree, add to the whole diagonal.

| degree K | the two halves |
|---|---|
| 1 (classical) | the y-half has weight 0 — the intensive readings do not scale; no conjugate degree; Φ's own response is already null |
| 2 | equal weights; G = ½(x − y)², symmetric under exchange of the halves |
| 1 < K < 2 (CF-3's window) | the y-half lighter than the x-half; K* > 2 |

**T6 (cutting one pair).** Setting one pair of G to rest replaces the response C by (AC − B²)/A, the Schur
complement. C_P/C_V = K_S/K_T (CF-3 (7)) is this one ratio read on the two pairs.

**T7 (sign).** G ≥ 0 near the rest set exactly when the response is positive; outside the disc G takes both
signs (checked on exact points).

## Reading

The generalised potential asked for is a function of two points — equivalently of both halves of the
diagram. The diagram's relations are its rest set (T1). The cut is the choice of base point, and what the cut
supplies is exactly what G does not contain: the value and the first layer (T3). Everything from the response
level upward is already in G before any cut. Classical thermodynamics is the degree at which one half of
the diagram carries no weight at all (T5).

## What is put in, what is not claimed

* G is built from a potential Φ in an admitted chart. It is a candidate for "the potential before the cut"
  in the sense that it needs no choice of independent half; it is not the uncut ground of the manuscript.
* Lineage: in content this is the classical two-point function of a convex potential (in thermodynamics, the
  work available relative to a reservoir). The stage's content is its place here: rest set = diagram,
  null doubled response at ρ = 1, layers = tower from the second on, weights by degree.
* T4 is checked for pure-degree potentials (one pair symbolically at exact points, two pairs for the radial
  family); a general conjugate is not constructed. T7 is local.
* Why Φ exists at all — why the two-point potential should split as Φ(x) − Φ(x₀) − … — is not shown.
  A two-point function that is not of this form would have a rest set that is not a diagram.

## Later note (UP4)

This stage uses derivatives in a chart, so it assumes that the two cut operations commute. UP4 removes that
assumption and gives the term it hides (w = ½[D_S, D_V]U); the statements here about a potential's own
equation, a rest set, or a bound on the ratio are the w = 0 case.

## Open gates

1. Characterise the two-point functions whose rest set is a diagram (half-dimensional, with face relations):
   is the split form forced? This is "why the diagram is forced" stated as a theorem to prove or refuse.
2. The degree: with weights (1, K − 1) on the halves, what selects K.
3. The doubled null response and the tower: T_λ on the boundary ρ = 1 (RMG2-T3's fixed set) read in the
   doubled variables.
