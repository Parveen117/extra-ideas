# UP5 — the seed document's hierarchy, checked

Fifth stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up5_seed_hierarchy_checked.py`, `test_up5.py` (4 tests).

The owner supplied his seed document "Thermodynamic Hierarchy and EMK Evolution" (layers
L₀ = (T,V,S,P) → L₁ = (λ_t, λ_v, λ_s, λ_p) → L₂ → …, UGD higher-order numbers, higher-order lambda
thermodynamics) as the statement of what "before the cut" means. This stage sets its definitions and claims
beside UP3 and UP4.

## Seed definitions used

  λ_p = −ST/C_p,  λ_v = −ST/C_v,  λ_s = V(∂P/∂V)_S,  z_s = −PV/λ_s,
  λ_t = −P(∂V/∂P)_T,  z_t = −PV/λ_t,
  λ_p^{(n+1)} = −λ_s^{(n)}λ_t^{(n)}/Q_p^{(n)} and its three companions, Q^{(0)} = (C_p, C_v, z_s, z_t).

## Statements

**T1 (the four Maxwell-type relations of a layer are one equation).** For any four readings (t, v, s, p), each
of the seed's four relations F, G, H, U equals ({t,s} − {p,v}) divided by one bracket. A layer has its
Maxwell pattern exactly when its two axes carry the same area form. Layer 0 does.

**T2 (what the seed's layer 1 satisfies, with no potential assumed).**

  λ_p/λ_v = z_t/λ_s,  λ_p λ_s = λ_v z_t,  λ_p λ_s λ_t = −PV·λ_v.

With λ_t as written, the layer does not close on its own four readings: λ_pλ_s = λ_vλ_t fails, and the
relation that holds needs P and V of the layer below. The four readings that do close are
(z_t, λ_v, λ_s, λ_p) — the mirror in which both mechanical capacities are P·∂V/∂P; that is the choice
UP3 made and the corner layer of RMG2-T5.

**T3 (the capacities of the next layer).** The seed defines λC_p, λC_v, λZ_s, λZ_t only through the
quantities they are meant to produce. UP3's form supplies the missing definition —
Q_p^{(n)} = t_n(∂s_n/∂t_n) at fixed p_n, and so on, in brackets — and with it the recursion is a rule.

**T4 (layer 1 does not inherit the equation).** For a potential that is not a pure power, the closed layer 1
has {t,s} ≠ ±{p,v} (exact witness in the tests). The seed's "lambda-Maxwell pattern" is therefore a condition
on a layer — UP3-T4's criterion — and in general layer 1 is a diagram without a potential.

**T5 (the general defect of two cut operations).** With a frame part (α, β) and a weight part F,

  [D_S, D_V] = α D_S + β D_V + F,  B₂ − B₁ = α T − β P + F·U.

The defect pairs with the value of the potential and with its first layer, and with nothing higher: exactly
the two things that the two-point potential of UP2 does not contain. This is the seed's connection
∇ = ∂ + A read on a potential: F is its curvature, (α, β) the order defect of the frame.

**T6 (the seed's core statement).** "The visible state returns but the responses do not." With commuting
cuts every reading returns with the state, so the statement is false there. With a defect it is true:
UP4-T5 is an exact instance. The statement is a statement about non-commuting cuts and cannot be reached
from inside the commuting calculus.

**T7 (an instance of closure).** For a potential of pure power, UP3-T5 gives ratio 1 and coinciding axis ends
from the third layer on: every distinction between the ends has gone after finitely many layers. This is one
exact case of the seed's closure of the tower; the general case is not shown.

## Not established by the seed document or here

* The lambda entropy S_Λ = log‖Hol‖, its growth law, the curvature flow and the field equation of the seed
  are stated without derivation, and S_Λ needs a norm that the document does not fix. They are open.
* A second note supplied with the seed argues that reciprocity fails for the entropy-weighted couplings
  because C_p ≠ C_v. Its closing section reaches the standard result instead: the antisymmetric part of a
  response is reversible and is absent from the dissipative coefficients at zero frequency. The first argument
  takes entropy as odd under reversal of motion and passes over several steps; it is not a proof. What is
  exact is T5 above and RMG9-T3: an antisymmetric part exists exactly when the cuts do not commute, it is
  invisible to every quadratic reading, and it is read by loops.

## What is put in

* The seed's definitions as written. The frame and weight of T5 are free; nothing selects them.
* Two pairs. F2–F4 of UP4 remain assumed.

## Later note (UP6)

The two forms of λ_t are the two members of one reciprocal pair of scale operations (UP6-T1); open gate 1
below is answered there: both are kept, and their non-commutation is computed.

## Open gates

1. Which of the two layer-1 quadruples the owner intends — (λ_t, …) as written, or the closed one with z_t.
2. S_Λ with a native mass in place of the unspecified norm, for the loop of UP4-T5.
3. How the defect (α, β, F) passes from layer n to layer n + 1.
