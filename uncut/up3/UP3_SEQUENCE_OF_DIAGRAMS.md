# UP3 — the sequence of diagrams: what holds at every level with no potential

Third stage of the uncut line. No measurement and no physical constant.
Symbolic algebra (sympy): `up3_sequence_of_diagrams.py`, `test_up3.py` (5 tests).

The owner's statement: in the next diagram P is replaced by the product of the ends of the other axis over a
capacity, S T/C_p, and V by S T/C_v; likewise T by P V over a new property. Write it with indices —
P₂ = S₁T₁/C_{P,1} — so that each diagram is built from the one before, and generalise the sequence.

## Notation (proposed)

Level n is four readings (T_n, S_n, P_n, V_n) on one surface of states. Level 1 is the T–V–S–P diagram.

  C_{X,n} = T_n (∂S_n/∂T_n) at fixed X_n   (X = P or V) — the capacities of level n
  W_{X,n} = P_n (∂V_n/∂P_n) at fixed X_n   (X = T or S) — the "new property" on the other axis

  P_{n+1} = S_n T_n / C_{P,n}   V_{n+1} = S_n T_n / C_{V,n}
  T_{n+1} = P_n V_n / W_{T,n}   S_{n+1} = P_n V_n / W_{S,n}

Equivalently P_{n+1} = S_n (∂T_n/∂S_n)_{P_n}, T_{n+1} = V_n (∂P_n/∂V_n)_{T_n}, and so on. No sign is added;
for a stable fluid T₂ and S₂ come out negative (they are minus the two bulk moduli).

With the bracket {f, g} of two readings on the surface, (∂T/∂S)_X = {T, X}/{S, X}. The rule then needs
**no potential**: it is defined for any four readings.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Publications research branch, `thermo-compass-foundations` results index | Jacobian brackets; C_P/C_V = K_S/K_T |
| Recognition-Kernel-Framework theorum/42 | the λ-Jacobian tower |
| response-geometry RMG2 T5 | corner layer λ_t = VΔ/a, λ_s = −Vc, λ_p = −SΔ/c, λ_v = −Sa; closure shown only for separable U |
| response-geometry RMG6 | memory weight r² = 1 − C_V/C_P |
| extra-ideas `uncut/up1`, `uncut/up2` | degree; two-point potential |

## Statements

Write χ_n = {T_n, P_n}{S_n, V_n} / ({S_n, P_n}{T_n, V_n}) — the ratio C_{V,n}/C_{P,n}.

**T1 (every level, no potential).** P_{n+1}/V_{n+1} = T_{n+1}/S_{n+1} = χ_n. So at every level from the
second on, P·S = V·T, and the two axes have the same ratio of ends. χ is unchanged when the two axes are
exchanged and inverted when the ends of one axis are exchanged. (At level 1 this is C_P/C_V = K_S/K_T; it
does not use the potential.)

**T2 (the three-term identity).** {T,S}{P,V} − {T,P}{S,V} + {T,V}{S,P} = 0 for any four readings, hence
1 − χ = −{T,S}{P,V} / ({S,P}{T,V}).

**T3 (what the potential adds).** A potential makes the two axes carry the same area form, {T,S} = {P,V}.
Then 1 − χ₁ = B²/(AC) — a square over the two diagonal responses — and χ₁ ≤ 1 wherever the response is
positive. Level 2 of a potential U(S, V) is exactly RMG2-T5's corner layer.

**T4 (the shape of every later level).** From the second level on a diagram is (χ·s, s, χ·v, v): two readings
s, v and the ratio of the level below. Its area forms are {T,S} = s{χ, s} and {P,V} = v{χ, v}. So:

* the two ends of an axis are different readings only where χ of the level below varies;
* the level has a potential exactly when {χ, s² − εv²} = 0 for ε = +1 or −1;
* its own ratio is 1 − χ′ = −pq/((1 − q)(1 + p)), p = s{χ,v}/(χ{s,v}), q = v{χ,s}/(χ{s,v}).

**T5 (a potential of pure power form stops the sequence).** For U = S³/V: χ₁ = 1/4, constant; χ₂ = 1;
χ₃ = 1; at level 3 the two ends of each axis coincide. In general a constant χ_n gives χ_{n+1} = 1 (p = q = 0):
the next diagram has no distinction left between the ends of its axes.

**T6 (without a power form the ratio is no longer bounded).** For
U = S²/2 + SV/3 + V²/2 + S²V/5 + V³/7, inside the region of positive response, χ₁ stays in (0, 1] while

  χ₂ = 11636647/11469655 > 1 at (S, V) = (1/2, 1/2),  χ₂ = −70384/167315 < 0 at (3/2, 1/2).

The bound χ ≤ 1 belongs to the level that has a potential. It is not a property of the sequence.

## Reading

The rule that builds the next diagram does not need a potential, and neither do T1 and T2. What the
potential — the cut — adds is one equation, {T,S} = {P,V}, and with it the sign of 1 − χ. Level 1 has it.
Level 2 in general does not: it is a diagram before a potential. The sequence is driven by the variation of
one number per level, χ_n; where the potential is a pure power (UP1's pure degree) χ is constant and the
sequence dissolves in two steps.

## What is put in, what is not claimed

* The rule is taken from the owner's statement and his earlier definitions (λ_p = −sT/C_p, λ_v = −sT/C_v,
  v(∂p/∂v)_T, v(∂p/∂v)_S). If the "new property" on the T–S axis is meant differently, T_{n+1}, S_{n+1} change
  and T1 must be re-checked.
* Lineage: the bracket identities are the classical algebra of four directions in a plane (χ is their
  cross-ratio; T2 is the three-term relation). The stage's content is the level rule in that form, the
  shape (χs, s, χv, v), the potential criterion and the two witnesses.
* T5 is shown for one power potential and by the formula of T4; T6 is two exact witnesses, not a
  classification of where χ₂ leaves (0, 1].
* Two pairs of readings only. Nothing here says what the higher levels measure.

## Open gates

1. The generating object of the whole sequence: is there one function (as UP2's G is for the layers) whose
   successive readings are the levels?
2. The fixed points of the rule: which diagrams reproduce themselves up to scale.
3. Three pairs: the rule and the invariant that replaces χ.
