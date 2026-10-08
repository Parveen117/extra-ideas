# UP1 — the potential before the cut: its degree decides where the diagram is silent and where the tower acts

First stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up1_degree_of_the_potential.py`, `test_up1.py` (6 tests).

The owner's statement: the physics done so far stands on a special case — the Euler relation itself assumes
flatness. Work with the generalised λ and the uncut potential: something like a flow of the uncut potential at
the centre of the diagram, depending on the degree in the λ-tower; λ can be taken as zero anywhere.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Publications research branch, `coherence-first-thermodynamics` | CF-1 the centre is a boundary of the cut chart; CF-3 W = κr^ν + εs²/2 with 1 < ν < 2; CF-4 response density ~ r^(ν−2) |
| same branch, `uncut-cut-measurement` MANUSCRIPT §1 | the uncut ground is not a set of states; a description declares its distinctions |
| Publications `cut-first-equivalence` | the classical laws are the flat, memoryless sector |
| Recognition-Kernel-Framework theorum/42, Theorem 3.1 | typed tower L_{n+1} = ∇L_n |
| response-geometry RMG2 T2, T3 | tower T_λ(H) = H + λH²; κ = (1 + 2λm)/(1 + λm(1 + ρ²)); fixed sets centre and boundary |
| extra-ideas physics LT1 | one generation at infinite λ is squaring |

## Set-up

Φ(x₁ … x_n) is the potential before any cut is chosen; y_i = ∂Φ/∂x_i. Cutting pair i replaces the reading
of Φ by Φ − x_i y_i. The 2ⁿ ways of cutting are the corners of the diagram (n = 2: the four corners of the
compass). A depth t_i ∈ [0, 1] on each pair gives the point Φ_t = Φ − Σ t_i x_i y_i inside it.
Φ has degree k when Φ(sx) = s^k Φ(x); the generator of the dilation, E = Σ x_i ∂_i, then reads E Φ = kΦ
(differentiate in s at s = 1).

## Statements

**T1 (each layer of the tower drops the degree by one).** E(∂Φ) = (k − 1)∂Φ, E(∂²Φ) = (k − 2)∂²Φ, and so on:
layer n has degree k − n.

**T2 (the corners).** Every face of the diagram closes: Φ_A + Φ_{A∪{i,j}} = Φ_{A∪{i}} + Φ_{A∪{j}}
(for the compass, U + G = H + F). The mean of all corners is (1 − k/2)Φ. The far corner is (1 − k)Φ.

**T3 (the silent point moves with the degree).** Along the diagonal, at depth t, the reading is
(1 − tk)Φ. It vanishes at t = 1/k:

| degree k | where the diagram is silent |
|---|---|
| 1 | the far corner — the Euler relation Φ = Σ x_i y_i |
| 2 | the centre — a quadratic potential, constant response |
| k > 1 | depth 1/k on the diagonal |
| k < 1 | nowhere inside |

The Euler relation of classical thermodynamics is the case k = 1 of this table and nothing more general.

**T4 (the diagonal sorts a potential by degree).** For a sum of parts of different degrees, depth 1/k removes
exactly the part of degree k and keeps the others with weight 1 − k′/k. For CF-3's potential:

  centre: ½(T₀s − P₀v) + (1 − ν/2)κr^ν — the regular quadratic part εs²/2 is gone;
  far corner: (1 − ν)κr^ν − εs²/2 — the linear part is gone.

CF-3's window 1 < ν < 2 is exactly the range in which the formation energy is positive at the centre and
negative at the far corner: strictly between the two silent cases. (Benchmark ν = 3/2, κ = 2/3, r = t⁴:
t⁶/6 at the centre, −t⁶/3 at the far corner.)

**T5 (λ and place are one parameter).** The response element H = ∂²Φ has degree k − 2, and

  T_λ(H(sx)) = s^(k−2) · T_{λ s^(k−2)}(H(x)).

Moving by the factor s in the diagram is the same as changing λ to λ s^(k−2). Only the product λ·m is ever
read (RMG2's formulas contain λ only as λm), so any value of λ can be taken as the reference: it is a choice
of the unit of place. For k = 2 — and only then — λ and place separate.

**T6 (the two ends, for k < 2).** Far from the centre λm → 0 and the tower is the identity: that end is the
λ = 0 special case. At the centre λm → ∞ and one generation is ρ → 2ρ/(1 + ρ²), with no λ left in it:
at the centre every λ gives the same map (LT1's squaring). For k > 2 the ends exchange.

**T7 (which layer is at rest).** Layer n is unchanged by the dilation exactly when n = k. For k = 1 it is the
first layer — the intensive readings. For k = 2 the second — the moduli. When k is not a whole number
no layer is at rest: layers below k grow outward, layers above k grow toward the centre (CF-4's divergent
response density is layer 2 with ν < 2).

## Reading

Three things that were treated as separate choices are one number, the degree k of the potential before
the cut: which point of the diagram is silent (T3), which layer of the tower can serve as the flat
reference (T7), and how λ is tied to place (T5). Classical thermodynamics is k = 1; linear response is
k = 2; CF-3's construction sits between them, where no corner, no centre and no layer is at rest.

## What is put in, what is not claimed

* Φ here is a potential in an admitted chart of variables, before a choice of which member of each pair
  is read. It is not the uncut ground itself, which the manuscript does not identify with any function or state.
* The depth t between corners is the straight interpolation Φ − Σ t_i x_i y_i. It is a definition, not a derived
  partial cut.
* T1–T3 are checked on a family with free coefficients, free real exponents and a radial term, in two and
  three variables; the general statement is the one-line argument in the set-up.
* Why the variables come in pairs — why the diagram exists — is not shown. Given the pairs, its corners and
  face relations are forced; the stronger statement is open.
* Nothing fixes k. No physical system is assigned a degree.

## Open gates

1. What fixes the degree: a framework statement that selects k (CF-3 leaves ν free in 1 < ν < 2).
2. Why pairs: derive the pairing from the cut itself (K² = 1 gives two sides; reading and lost part as the pair).
3. Non-integer degree: the tower with no layer at rest — its fixed sets and whether RMG2-T3's attractors survive.
