# UP4 — the cost of assuming that the two cuts commute

Fourth stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up4_cost_of_commuting_cuts.py`, `test_up4.py` (5 tests).

The owner's objection to UP1–UP3: the tools were applied after flatness had been assumed; the assumption
has a cost, and it was ignored — the uncut was being looked at with the tools of the cut.

## What was assumed

UP1–UP3 used partial derivatives in a chart. That contains, unstated:

| | Assumption | Where it entered |
|---|---|---|
| F1 | the two cut operations commute, ∂_S∂_V = ∂_V∂_S | UP2 (rest set of G), UP3-T3 (the potential's equation), UP3-T4 |
| F2 | every reading has a value that does not depend on the way it was reached | the rule S·T/C itself; UP1's Φ |
| F3 | the states have exactly as many free directions as there are pairs | every bracket; UP3-T2 |
| F4 | the readings are smooth at the centre | UP1's degree at r = 0 (CF-1 already says the centre is a boundary) |

This stage removes F1 and computes what it had hidden. F2–F4 are not removed here.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| response-geometry DB1 | derivation algebras [D_a, D_b] = f_abc D_c: identities from the bracket and Leibniz only |
| response-geometry RMG9 T1, T3 | L = m(1 + uK + vS) + wR; det L = det L_sym + w²; a closed protocol returns 2w·Area |
| Recognition-Kernel-Framework theorum/55 T3, theorum/57 | the loop of two flows leaves h²[D₁, D₂] |
| Publications `cut-first-equivalence` | the classical laws are the memoryless sector; the loop residue is the obstruction |

## Set-up

Two operations D_S, D_V with Leibniz and

  [D_S, D_V] = α D_S + β D_V

(α, β: the order defect of the pair of cuts). For a potential U: T = D_S U, P = −D_V U,
A = D_S T, C = −D_V P, and the two mixed responses B₁ = D_V T (T along the V operation),
B₂ = −D_S P (−P along the S operation). Write B̄ = (B₁ + B₂)/2, w = (B₂ − B₁)/2.

## Statements

**T1 (the mixed responses differ by the defect times the first layer).**

  B₂ − B₁ = [D_S, D_V] U = α T − β P = 2w.

The four Maxwell relations each carry this term. It is made of the readings themselves — the first layer,
the one UP2's two-point potential does not contain.

**T2 (the ratio is no longer below one).** Both axes still share one ratio χ = det/(AC), but

  1 − χ = (B̄² − w²)/(AC),  det = det(symmetric part) + w².

With commuting cuts w = 0 and 1 − χ is a square: χ ≤ 1. With the defect, χ > 1 as soon as w² > B̄².
The inequality C_P ≥ C_V is a consequence of F1, at level 1 already. The response element is RMG9's
L = m(1 + uK + vS) + wR, and the R-part that RMG1–8 set to zero is w = ½[D_S, D_V]U.

**T3 (what does not use F1).** The shared ratio of ends and the three-term identity (UP3-T1, T2).

**T4 (witness).** [D_S, D_V] = D_V (α = 0, β = 1), U = S² + V² − SV/2 + 2V at the origin:
A = C = 2, B₁ = −1/2, B₂ = 3/2, P = −2, symmetric part positive, **χ = 19/16**. With the same defect and
P = 0 at that point, B₁ = B₂: the cost vanishes with the first layer.

**T5 (the loop of the two cuts does not close).** Four steps of size ε — S, V, back S, back V — end
ε²[D_S, D_V] away from the start, and the potential differs by

  ε²(α T − β P) = 2w·ε².

RMG9-T3's "2w·Area around a closed protocol" is this: the potential read across the gap that the two
cuts leave open.

## Reading

F1 hid one term, and the term is not small in principle: it is the order defect of the two cuts times
the first-layer readings. Everything in UP1–UP3 that mentioned a potential's own equation, a bound on χ,
or a rest set, is the w = 0 case. An earlier reading in this project judged the statement "∂_S∂_VΨ ≠ ∂_V∂_SΨ"
unfounded because mixed derivatives commute; that judgement assumed F1. With non-commuting cuts the
statement is exact and its value is [D_S, D_V]Ψ.

## What is put in, what is not claimed

* D_S and D_V are still two operations acting on readings. That is itself a cut: this stage prices one
  assumption of the cut-side calculus, it does not describe the uncut.
* F2, F3, F4 remain. In particular the rule of UP3 still needs values of S and T (F2).
* The concrete model (D_S = ∂_a, D_V = f∂_b + g∂_a) is used to check identities that follow from the bracket
  and Leibniz; it is not the definition.
* Nothing fixes α, β. No physical system is assigned a defect.
* UP2's statement "the diagram is the rest set of G" and UP3-T3, T4, T6 are to be read with w = 0.

## Open gates

1. Remove F2: readings as responses to steps, without values — the rule of UP3 and the degree of UP1 in that form.
2. The sequence of UP3 with the defect: how w passes from level n to level n + 1.
3. What a defect (α, β) would be selected by — the degree (UP1) is the first candidate, since [E, D] = −D for
   the dilation E and any first-order operation D of weight one.
