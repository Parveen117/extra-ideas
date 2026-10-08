# Checklist — today's uncut stages against what the repository already had (8 October 2026)

Made by reading `physics/SYNTHESIS.md`, `PHYSICS_DOCK.md`, `GRAVITY_THESIS.md`, `UNUSED_RESULTS_MAP.md`, the claim
boundary of every stage in `physics/`, and in full: FR1, IN1, OB1, HB1, CL1, TL1, SC1, DC1, QC1, QC2, WQ1, BC1, GW1,
LC1, EG1, DG1, MC1, PH1, PR3. Top-level `CURRENT_STATUS.md`, `ASSESSMENT.md`, `PHYSICS_INGREDIENTS.md` were skimmed
(they concern R1–R46, GE1–GE4, MP1–MP2), not checked line by line.

## 1. What the line already had

| Subject | Stage | Statement | Status there |
|---|---|---|---|
| three directions | FR1 | two cuts over the rationals, three on the cut-complex carrier, no fourth | PROVED |
| the turn | OB1 | ι = C₁C₂C₃; block = (a + ιb) + (r + ιs)·C | PROVED |
| four components | IN1 | ρ = ½(n + r·C); det = ¼(n² − r·r); the same in every frame | PROVED |
| time | CL1, GE2 | time is a count along a history; a field equation for a unit of time is refused | PROVED / REFUSED |
| the unit | HB1, DC1, SC1 | a frame's tick unit; frames that exchange count in one unit; a local scale field is refused (history-dependent invariant) | PROVED / REFUSED |
| heat and the clock | TL1, PR3 | κ·N constant in equilibrium; accelerated frame κ·ρ constant; a temperature needs the unit | PROVED |
| whole numbers | QC1, WQ1, LC1, GW1 | content q silent ⇔ qΘ ∈ 2πℤ; recordable mode E = nκω; lost part counted; frame waves the same | PROVED |
| the horizon's count | BC1 | boundary adds, never returns; a count of the horizon | NOT DERIVED; quantisation there REFUSED |
| gravity | GR1, MC1, MO1, GR2, CV1, TP1, NC1 | memory dictionary; 1/r from least cost; motion; curvature; the law 1 : 2 : −4 | PROVED on stated assumptions |
| charge and gravity | EG1 | memory r_s/r − q₂/r²; energy component (rm)′/r² | PROVED (constant not derived) |
| the diagonal | DG1 | observed = lost at r = 2r_s, N = 1/√2 | PROVED |
| real fluids | PH1, QC2 | return angle for seven fluids on a reduced cycle; ideal gas flat for any c_v in the entropy representation | evidence packet / PROVED |

## 2. Today's stages

**Restate something the line already had** (kept, with a pointer):

| Today | Already in | Remark |
|---|---|---|
| WD2-T1/T2, FC4-T1/T2/T3/T6 | FR1, OB1, IN1 | three cuts, no fourth, ι as their product, four components and the invariant |
| AR1-T3, RT1-T2 (T·N constant) | TL1-L4, PR3 | including the accelerated frame |
| LP1-T1; LP1-T3's joining rule | MC1; EG1-G2/G3 | 1/r for a partner; "centre value at r = far value − cost outside" is EG1's energy component integrated |
| FD1: ideal gas flat for any heat capacity, attraction the source, monatomic agreement | PH1 Table 2, QC2 §2 | D′ is a local form of the same curvature: at low density f·du·dv = [D′/(2√c_v)]·d ln T·d ln v (`fd1_relation_to_response_curvature.py`) |
| CK1 gate (whole numbers) | BC1, WQ1 | answered there; QH1 adds the size of the step |
| ET1-T3 (T = 0 where N² = v²) | DG1 | the same diagonal, for the rim reader |

**Conflict with the line; corrected by notes in the stages:**

| Today | Conflicts with | Correction |
|---|---|---|
| UN1-T3/T4, UN2-T2/T5: "the unit of a frame is its clock factor; count = area/(4·c_f)" | DC1-K1 (frames that exchange count in one unit), TL1-L2 (u·N constant between places) | withdrawn as statements about the count. HB1's tick unit is not the unit of the count. UN2-T1 (legs of the coin) and the quarter with one unit stand. |
| PL2; XU1-T4…T6; XU2-T3; XU3; LG1; the weight part of ME1's master equation | SC1 (physics): a scale that varies is a scale field; it makes the invariant depend on the history and is refused as fundamental | the path-dependent count of XU3 is SC1's refused case, obtained by letting the unit follow the state. It is not an arrow. XU1-T1/T2 (the best frame's clock factor), XU2-T1/T2/T4 (the expansion plane; both horizons) stand. |
| AL1, ST1, QT1, SM1: a count for the horizon | BC1: not derived | no contradiction — these stages define the count through T = unit·rate per turn, an input BC1 does not grant. Recorded as a definition, not a derivation. |

**New in this repository** (their external counterparts, where known, are named in each stage):
ON1, SS1, CP1, DW1 (commutator of scale readings = running; scale-response law; turn-part from the potential);
BH1, BH2 (shape of the turn-part); PL1, XU2-T2 (count as areas in planes); HX1 (fall as a turning of cuts);
WD1 (boundary in whole turns only for d = 2, 3); WD2-T3/T4 and NA1-T5 (one eighth of a turn per octant);
LP1-T4, SP1, CF1, ET1, ET2, HD1, HL1, TS2, CT1 (centre potentials and the heat-reading law);
MT1, EQ1 (readings of a moving body as channels); SL1 (one unbiased step gives cosh 2η); CR1 (chromium);
QH1.

## 3. Names

Four of today's stages had taken identifiers already used in `physics/` (NC1, SC1, LC1) and one followed them
(LC2). They are renamed: NC1 → ON1, SC1 → SS1, LC1 → LP1, LC2 → LP2. In `uncut/`, every reference now uses the new
names; "NC1", "SC1", "LC1" mean the physics-line stages only.

## 4. Not read

The other repositories (Recognition-Kernel-Framework, RH-Framework, Publications) were not re-read for this
checklist; `UNUSED_RESULTS_MAP.md` lists results there that no stage has used yet.
