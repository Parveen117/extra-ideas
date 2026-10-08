# CO1 — an expanding space is a fall frame; the line's law gives the equations of the coupled solution

Stage of the physics line. Symbolic algebra (sympy): `co1_expansion_as_fall.py`, `test_co1.py` (4 tests).

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Publications research branch, `ugd-kahler-propagation/COUPLED_COSMOLOGICAL_SOLUTION.md` | CS.7 n = n₀/v; CS.9 ρ = mn + (3κ/16)n², p = (3κ/16)n²; CS.10 lapse form; CS.11; CS.13; CS.15 |
| this line | TP1 (law Q = ¼I₁ + ½I₂ − I₃ on the frame's order defect); MO1 (fall speed, β² = r_s/r, and its assumption: flat slices, one time); GB1; EN1 |

The research branch reached its cosmology from a declared action (matter with torsion eliminated).
The line reached its law from equivalence. Neither had been checked against the other.

## Statements

**T1 (the comoving frame is a fall frame).** In proper distances X = a·x the comoving frame is

  e₀ = ∂_t + h X·∂_X,  e_i = ∂_{X_i},  h = ȧ/a :

flat slices, one time, and a pure fall of speed β = hX directed outward. It is the same frame field
in two coordinate systems; I₁, I₂, I₃ are 6h², 3h², 9h² in both, and the law is Q = −6h².
MO1's assumption, which the line could not prove for a general field, holds exactly here.

**T2 (the line's law gives CS.11 with CS.9).** Varying N a³[Q/2κ − Λ/κ − ρ(a)] with
ρ(a) = m n₀/a³ + (3κ/16)n₀²/a⁶:

  lapse: 3h² = Λ + κρ,  scale: 2ḣ + 3h² = Λ − κp,  p = (3κ/16)n₀²/a⁶ .

The pressure is not put in; it comes from the a⁻⁶ dependence of the torsion contact term.

**T3 (the fall law, shell by shell).** The lapse equation is β² = r_s(X)/X with
r_s(X) = (Λ + κρ)X³/3. For Λ = 0 and κ = 8πG that is r_s = 2GM(X), M the mass inside the shell:
MO1's law for a point mass holds for every shell of the expanding space with the enclosed source.

**T4 (the closed solution).** v = a³ = d(cosh ωt − 1) + b sinh ωt (CS.15) satisfies the lapse
equation of T2 identically. The right side of v̇² is a perfect square exactly when m = √(3Λ)/2,
and then v = b(e^{ωt} − 1) — the branch's benchmark a³ = eᵗ − 1 is this case.

**T5 (constant rate: the redshift factorises).** For constant h, light from proper distance X
arrives with 1 + z = 1/(1 − β), β = hX, and

  1/(1 − β) = √((1 + β)/(1 − β)) · 1/√(1 − β²) :

the Doppler factor of the fall (EN1's bound, square-rooted) times the inverse clock factor of the
gravity element (GB1). The horizon is β = 1.

## What is put in, what is not claimed

* ρ(a) is taken from CS.9; the matter action and the torsion elimination behind it are the research
  branch's and are not re-derived here. 1/2κ in front of Q fixes what κ means.
* T1's agreement between the two forms is coordinate invariance, not a test of equivalence.
* T5 is for constant rate only. For CS.15 the redshift is (v_o/v_e)^{1/3} and is not a function of β alone.
* No measured number: n₀, m, Λ, κ are free. Nothing is fitted to the observed expansion.
* RB1/TS1's boundary-memory verdict for the horizon at β = 1 of an expanding space is not worked out here.
