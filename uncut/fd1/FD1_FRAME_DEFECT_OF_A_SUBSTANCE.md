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

## Step and rounding check (`fd1_step_check.py`)

D is a derivative taken from a table, so two errors compete: the step (truncation) and the printed digits
(rounding, amplified as the step shrinks). Argon, 300 K, 1 mol/l:

| span of the difference | D | rounding bound |
|---|---|---|
| 40 K (used in the tables above) | 0.095584 | 0.000054 |
| 20 K | 0.095209 | 0.000108 |
| 10 K | 0.095122 | 0.000215 |
| Richardson (20, 10) | 0.095093 | — |

The 40 K value is 0.5 % above the extrapolated one; the rounding bound is 0.06 % there and reaches 0.2 % at 10 K.
Both are far below the 2–8 % spread and the 16–32 % departures reported above.
The tables are values of a smooth reference equation, not sampled measurements: there is no sampling-rate
(aliasing) question here. For raw calorimetric or acoustic records the derivative would need the usual
conditioning — band limit, anti-alias filter, a stated step — before D is formed.

## The improved number D′, and seven substances on one curve (`fd1_second_scale_check.py`)

The departures of nitrogen and carbon dioxide reported in the collapse check above were mostly an artefact of
the definition: D still contains the temperature dependence of C_V itself. Removing it at the level where an
ideal gas has it gives

  **D′ = (T/C_V)·∂[C_P/(αT) − C_V]/∂T at fixed V.**

- **FD1-T6 (exact).** D′ = 0 for an ideal gas with any heat capacity C_V(T) (internal motion drops out);
  D′ = D for van der Waals; at low density D′·V → [2T·B′ − (c_v − 1)·T²·B″]/c_v.

Data, same reduced density (ρ/ρ_c ≈ 0.0746), oxygen and carbon monoxide added:

| T_r | argon | krypton | xenon | oxygen | nitrogen | carbon monoxide | carbon dioxide |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.1807 | 0.1800 | 0.1808 | 0.1884 | 0.1921 | 0.1878 | 0.1795 |
| 1.7 | 0.1501 | 0.1491 | 0.1495 | 0.1549 | 0.1590 | 0.1570 | 0.1465 |
| 1.9 | 0.1283 | 0.1275 | 0.1276 | 0.1308 | 0.1356 | 0.1347 | 0.1244 |
| 2.1 | 0.1118 | 0.1111 | 0.1113 | 0.1126 | 0.1179 | 0.1177 | 0.1085 |

- The three monatomic gases agree within 0.5 %.
- All seven lie within 7 % of the monatomic curve; carbon dioxide within 3 %, although its heat capacity is
  three times larger and changes with temperature.
- What remains: oxygen +1 to +4 %, carbon monoxide +4 to +6 %, nitrogen about +6 %, carbon dioxide −1 to −3 %.

Corrected reading: with D′ the frame defect at this density is close to one function of (T_r, ρ_r) for atoms
and small molecules alike. The earlier statement that nitrogen and carbon dioxide "mark further scales" by
16–32 % and by a factor of four does not stand; the residual marks are a few percent.

Put in: as above (page-reader tables; critical constants and acentric factors from standard tables; linear
interpolation). One density only. Not claimed: universality at other densities or near the critical point;
that this collapse is absent from the literature (not checked).
