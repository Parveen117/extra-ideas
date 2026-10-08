# CI1 — one measure for both sides: the information shared between the two readings

Stage of the uncut line. No measurement. Exact (sympy), the logarithm checked through its argument:
`ci1_shared_information.py`, `test_ci1.py` (6 tests).

The owner's statement: measure space-time and thermo against something else that is present in both — not
against time. An invariant; something like entropy and the arrow of time.

## Sources used (read, unchanged)

| Source | Statement used |
|---|---|
| Publications `emk-ugd-algebra` EMK-1 | det(a + bK + cS + dR) = (a² − b²) + (d² − c²): seen channel + lost channel |
| Publications research branch, CF-3 (7) | C_V = T/A, C_P = TC/det H |
| response-geometry RMG1, RMG2 T2–T3, RMG6, RMG10 T1, T4 | rapidity ℓ, tanh(ℓ/2) = ρ; the tower; memory weight r²; det L = m²(1 − δ); the light cone |
| extra-ideas physics `gb1`, `ln1`, `in1`, `gr1`, `ms1` | unit block Exp(ψn) and its clock factor; seen and lost; speed² + memory = 1 |
| extra-ideas `uncut/up3`, `up4`, `tt1` | χ = C_V/C_P; w; one cycle = a boost |

## Definition

For an element M = a + bK + cS + dι:

  I_cut = ½ Log(seen / det), seen = a² − b²  — what the cut along K does not hold by itself;
  I₀ = ½ Log(a² / det)  — the same with no cut chosen.

Both are pure numbers. I₀ does not change under any change of frame on the carrier (turns and boosts);
I_cut depends on the cut.

## Statements

**T1 (thermo side).** For a response H = [[A, B], [B, C]]:

  exp(2 I_cut) = AC/det = C_P/C_V = 1/(1 − r²),  r² = B²/(AC).

r² is also the squared correlation of the fluctuations of the two axes, so I_cut is what one axis carries
about the other.

**T2 (turning the cut).** As the cut turns, I_cut runs from 0 — cut along the principal axes — to I₀ — cut
on the diagonal. And exp(I₀) = (mean of the two principal responses)/(their geometric mean) = cosh(ℓ/2).

**T3 (space-time side).** For a unit block Exp(ψn):

  exp(I₀) = cosh ψ = 1/(clock factor),  exp(−2 I₀) = 1 − speed².

A cut along the block sees everything (I_cut = 0); a cut across it has I_cut = I₀. For the fall,
exp(−2 I₀) = 1 − r_s/r, and far away I₀ = r_s/2r + …

**T4 (the same number on both sides).**

| | thermo | space-time |
|---|---|---|
| exp(−2I) | C_V/C_P = 1 − r² | (clock factor)² = 1 − speed² |
| I = 0 | no coupling between the axes | at rest, far from any mass |
| I → ∞ | C_P/C_V → ∞ (critical point, spinodal) | clock factor → 0 (horizon) |

**T5 (an arrow, without time).**
* Composing two blocks of the same sense: I₀(ψ₁ + ψ₂) ≥ I₀(ψ₁) + I₀(ψ₂), with equality only if one is trivial;
  opposite senses lower it.
* One generation of the tower at λ = ∞, or one cycle of an open diagram (TT1), at least doubles it:
  I₀(2ψ) ≥ 2 I₀(ψ).
* Along the tower: λ > 0 raises I₀, λ < 0 lowers it, λ = 0 keeps it. (Witness H = [[3,1],[1,2]]:
  exp(2I₀) = 289/244, 5/4, 16/11 for λ = −1/20, 0, 1/5.)

**T6 (with a turn-part).** For L = m(1 + uK + vS) + mt·ι: exp(−2 I₀) = 1 − ρ² + t². The turn-part lowers I₀.
I₀ > 0 where the response has a cut, I₀ = 0 on the light cone ρ = |t|, I₀ < 0 where it is a turn.

## Numbers (illustration)

```text
gas with C_P/C_V = 5/3            I_cut = 0.2554
the same I in space-time          speed 0.6325 of the cone speed, or a clock at 2.5 r_s
surface of the Earth              I = 7.0e-10          surface of the Sun    I = 2.1e-6
```

## Reading

There is a quantity present on both sides and not measured against time: how much the two readings of one
element share — equivalently, how much a single cut leaves out. On the thermo side its exponential is the
ratio of the two capacities. On the space-time side it is the inverse clock factor. It is zero exactly where
the element is flat for that cut, infinite at the critical point and at the horizon, and it is ordered: same-
sense composition, a generation of the tower and a cycle of an open diagram do not lower it.

## What is put in, what is not claimed

* That this number is the thermodynamic entropy, or that its ordering is the arrow of time, is not shown.
  T5 is an ordering under three named operations; nothing says nature only performs those.
* Lineage: for Gaussian fluctuations I_cut is in content the mutual information of two correlated variables;
  the stage's content is the definition through the determinant channels and the identity of the two sides.
* The space-time column uses the unit-block reading of GB1 and RB1 (the fall-speed dictionary). The thermo
  column is CF-3 (7).
* Two cuts on one carrier. For three cuts the measure is not built.
* I is a pure number; κ·I has the unit of fluctuation. No constant is derived.

## Open gates

1. A conservation or balance: what is exchanged for I when it rises (TL1's exchange law is the candidate).
2. I for three cuts (the reading tensor of IN1).
3. Whether a process on the thermo plane (RV1's river) raises I₀ whenever its strain is not zero.
