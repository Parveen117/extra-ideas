# FD1 — The frame defect of a substance as a measurable pure number, and a first check on argon

Sources read: NC1, SC1 (the defect is the response to a scale), DW1-T2 (w₂ = ½(S∂_S ε⁻¹)·TS), CI1
(exp(2I_cut) = C_P/C_V), ME1 (candidate 3), UP6.

## Exact statements (sympy)

In measurable terms V·P_V/P_S = −C_P/(αT) (α = thermal expansion). Define the pure number

  **D = (∂/∂S)[C_P/(αT)] at fixed V = (T/C_V)·∂[C_P/(αT)]/∂T at fixed V.**

It needs no zero for the entropy.

- **FD1-T1 (ideal gas).** C_P/(αT) = C_P, D = 0, and the turn-part is w₂ = C_P·T/2 — half the enthalpy;
  w₂/U = γ/2 = ½·exp(2I_cut).
- **FD1-T2 (who is flat).** Radiation (U ∝ S^{4/3}V^{−1/3}) and a horizon (M ∝ S^{1/2}) are pure powers: ε is
  constant, no defect at all. The ideal gas has ε = −S/(γC_V), which runs with S, and obeys SC1 with the
  scale b₀ = C_V: its built-in scale is the count unit. Systems without a particle count are the flat ones.
- **FD1-T3.** Hard cores alone (P = NkT/(V − b)) give D = 0.
- **FD1-T4 (attraction is the source).** van der Waals: D = 2a(V − b)/(NkT·V²); in reduced variables
  D = (9/4)(1 − 1/(3v_r))/(T_r·v_r) — one function for every such substance.
- **FD1-T5 (low density, any substance).** D·V → −(T²/c_v)·[T·B‴ + (c_v + 3)·B″], B(T) the second virial coefficient.

## Check on data (`fd1_argon_check.py`, `data/argon_isochore_1mol_per_l.csv`)

Argon on the isochore 1 mol/l, 160–400 K, tabulated C_v, C_p and Joule–Thomson coefficient (αT = 1 + μ·C_p·ρ).
D from central differences, against FD1-T4 with argon's usual a = 0.1355 Pa·m⁶/mol², b = 3.201·10⁻⁵ m³/mol
(nothing fitted):

| T (K) | D from data | D van der Waals | ratio |
|---|---|---|---|
| 180 | 0.15113 | 0.17528 | 0.862 |
| 240 | 0.11952 | 0.13146 | 0.909 |
| 300 | 0.09558 | 0.10517 | 0.909 |
| 380 | 0.07395 | 0.08303 | 0.891 |

D is positive at all eleven temperatures, falls with T, and sits within 9–14 % of the one-scale form.
The ideal value is 0.

Put in: the table was retrieved from the NIST Chemistry WebBook through a page reader on 8 October 2026 and
should be re-downloaded directly before any publication; 20 K steps for the derivative; a, b from standard tables.
Not claimed: a new law of gases — D follows from the equation of state by ordinary thermodynamics. What the
line adds is the quantity itself (the frame defect, zero for the ideal gas and for hard cores, sourced by
attraction) and its reading as a turn-part.
Open gate: other substances and densities on the reduced form of T4 (the collapse is SC1's statement);
the departure from 1 of the ratio as a second scale.

## Collapse check across substances (`fd1_collapse_check.py`, five tables in `data/`)

SC1 says a family with one built-in scale has one D-curve in reduced variables. Five substances at one reduced
density (ρ/ρ_c ≈ 0.0746), D against T_r = T/T_c:

| T_r | argon | krypton | xenon | spread of the three | nitrogen | carbon dioxide | van der Waals curve |
|---|---|---|---|---|---|---|---|
| 1.4 | 0.1343 | 0.1268 | 0.1239 | 8.1 % | 0.1696 | 0.5688 | 0.1169 |
| 1.6 | 0.1190 | 0.1178 | 0.1157 | 2.8 % | 0.1446 | 0.5332 | 0.1023 |
| 1.8 | 0.1061 | 0.1079 | 0.1064 | 1.7 % | 0.1264 | 0.5011 | 0.0909 |
| 2.0 | 0.0952 | 0.0983 | 0.0972 | 3.2 % | 0.1123 | 0.4720 | 0.0818 |
| 2.2 | 0.0860 | 0.0897 | 0.0888 | 4.2 % | 0.1026 | 0.4451 | 0.0744 |

- The three monatomic gases fall on one curve: within 5 % above T_r = 1.5, 8 % at 1.4.
- Nitrogen lies 16–32 % above them; carbon dioxide four to five times above. For carbon dioxide most of the
  excess is a heat capacity that itself changes with T (internal motion): removing that part leaves 0.22, 0.14,
  0.11 at T_r = 1.33, 1.73, 2.12 — still above the monatomic curve.
- The van der Waals curve has the right shape and lies 13–17 % below the monatomic data.

Reading: one scale ⇒ one curve holds where the substance has one scale (atoms); each departure marks a further
scale (shape, internal motion). D separates them without a model.

Put in: critical constants from standard tables (argon 150.687 K, 13.4074 mol/l; krypton 209.48, 10.85; xenon
289.733, 8.4; nitrogen 126.192, 11.1839; carbon dioxide 304.128, 10.6249); linear interpolation to common T_r;
tables retrieved through a page reader — to be re-downloaded directly before publication.
Not claimed: that the collapse of simple fluids is new (corresponding states is long known); the new item
is D as the collapsing quantity and its meaning as frame defect.
