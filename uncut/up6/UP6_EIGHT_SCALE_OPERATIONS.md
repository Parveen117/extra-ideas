# UP6 — the eight scale operations of the diagram, and where they do not commute

Sixth stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up6_eight_scale_operations.py`, `test_up6.py` (5 tests).

The owner's statement: the two forms of the mechanical reading, −P(∂V/∂P)_T and V(∂P/∂V)_T, are not a slip
to be settled. Their mismatch is the reason for the connection of the λ-tower, and their non-commutation may
be exactly what is lost at λ = 0.

## Sources used (read, unchanged)

The seed document; `uncut/up1` (degree), `up3` (sequence), `up4` (cost of commuting cuts), `up5`.

## Set-up

A scale operation is E[x | y] f = x·(∂f/∂x) at fixed y. Each corner of the diagram has two:

| corner | operations |
|---|---|
| U(S, V) | E[S\|V], E[V\|S] |
| H(S, P) | E[S\|P], E[P\|S] |
| F(T, V) | E[T\|V], E[V\|T] |
| G(T, P) | E[T\|P], E[P\|T] |

## Statements

**T1 (eight operations, eight readings, four reciprocal pairs).** Each operation applied to the partner of
its own variable gives one reading:

| pair | readings | product |
|---|---|---|
| E[S\|V] T and E[T\|V] S | −λ_v and C_v | ST |
| E[S\|P] T and E[T\|P] S | −λ_p and C_p | ST |
| E[V\|S] P and E[P\|S] V | λ_s and −z_s | PV |
| E[V\|T] P and E[P\|T] V | form (B) and minus form (A) | PV |

Forms (A) and (B) are the two members of one reciprocal pair. Neither is an error; (B) scales V, (A) scales P.
The four λ of the closed layer are the four operations that scale S or V.

**T2 (inside a corner, the two operations commute).** All four corners. A corner is a flat chart.

**T3 (between corners they do not).** With m = (∂ ln S/∂ ln V)_T and e_T = (∂ ln P/∂ ln V)_T:

  E[V|T] = E[V|S] + m·E[S|V],
  [E[S|V], E[V|T]] = E[S|V](m) · E[S|V],
  [E[V|S], E[V|T]] = E[V|S](m) · E[S|V],
  E[P|T] = (1/e_T)·E[V|T],  [E[V|T], E[P|T]] = −E[V|T](ln e_T) · E[P|T].

The defect between two corners is a scale-derivative of a pure ratio. This is UP4's bracket
[D_S, D_V] = αD_S + βD_V, with α, β now computed.

**T4 (the cost appears inside ordinary thermodynamics).** Take one cut from the U corner and one from the F
corner, D_S = E[S|V], D_V = E[V|T]. Then D_S U = ST and

  B₂ − B₁ = E[S|V](m) · ST.

No exotic assumption is needed for an asymmetric pair of mixed responses: it is enough that the two cuts come
from different corners. The readings of the λ layer do come from different corners (T1).

**T5 (when everything commutes).** All the defects vanish exactly when the equations of state are power laws
(T and P each a product of powers of S and V). Witnesses at (S, V) = (1/2, 1/2):

| potential | E[S\|V](m) | E[V\|S](m) | E[V\|T](ln e_T) |
|---|---|---|---|
| S³/V | 0 | 0 | 0 |
| S³ + 2√V | 0 | 0 | 0 |
| S²/2 + SV/3 + V²/2 + S²V/5 + V³/7 | 5/18 | −10/27 | 1290275/1753182 |

For a gas with constant capacities, forms (A) and (B) commute (e_T = −1) while the U and F corners do not
(E[S|V](m) = −R/S, cost −RT).

## Reading

The owner's statement holds in this form. The mismatch between (A) and (B) is the relation between two
corners; its measure is E(ln e_T); and it is one of a family — every pair of corners has such a defect. These
defects are the connection that UP4 left free. They vanish together in the power-law case, which is where
the degree is pure (UP1) and where the sequence of diagrams stops (UP3-T5): the case in which the tower
does nothing.

## What is put in, what is not claimed

* The identification of "λ = 0" with the power-law case is by its effect (no defect, tower stops), not by a
  parameter set to zero. A statement in which a parameter λ multiplies the defect is not shown.
* T5's "exactly when" is argued in the chart of ln S, ln V (the corners' charts are related by a constant
  matrix iff ln T, ln P are affine there); the tests give witnesses, not the general proof.
* This is still the calculus of one surface of states: F2–F4 of UP4 remain. What has changed since UP4 is
  that the defect no longer has to be put in by hand — the diagram's own corners supply it.
* Two pairs only.

## Open gates

1. The sequence of UP3 written with the eight operations: which corner each reading of level n + 1 comes from,
   and the defect of level n + 1 in terms of level n.
2. A parameter form: a one-parameter family of potentials through a power law, with the defects to first order
   in the parameter.
3. The loop of two cross-corner operations (UP4-T5) in this setting: what a cycle S-at-V, V-at-T returns.
