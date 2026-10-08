# CR1 — Chromium at its Néel transition, in the line's quantities

Sources read: FD1 (zero-free turn-part w = C_P/(2α); defect numbers D, D′), DW1 (the sense of a cycle is the sign
of w), BH1/BH2 (a pole of w without change of sign), PL2/XU3/LG1 (a path-dependent count; gain = enclosed area).

Published inputs:
- T_N = 311.4 K (single crystal); latent heat 1.10 ± 0.10 J/mol on heating, 0.97 ± 0.10 J/mol on cooling;
  peak a few tenths of a degree wide; transition first order — Benediktsson, Åström and Rao (1975),
  "Calorimetric studies of the order of magnetic phase transitions in Cr and some Cr alloys at the Néel point".
- dT_N/dP = −5.1 K/kbar; chromium contracts on going from the antiferromagnetic to the paramagnetic phase;
  expansion coefficient larger above T_N — University of Toronto advanced laboratory handout on the Néel
  transition, citing Matsumoto and Mitsui (1969).
- Standard-table values, from general knowledge and to be verified: molar volume 7.23·10⁻⁶ m³/mol, volume
  expansion coefficient ≈ 1.47·10⁻⁵ K⁻¹ and C_P ≈ 23.3 J/(mol·K) near room temperature.

## Statements

- **CR1-T1 (exact).** For a first-order step smeared over a width δ, the anomalous parts of C_P and of α have a
  ratio that does not depend on δ, and inside the window the turn-part tends to
  w → V·T_N / (2·dT_N/dP). Its sign is the sign of dT_N/dP.
- **CR1-T2 (exact).** α changes sign inside the window exactly when δ is below δ* = |ΔV/V| / α_background.
- **CR1-T3 (chromium).** From the latent heat and dT_N/dP: ΔV/V = −2.34·10⁻⁵ on heating (a contraction, as observed).
  - critical width δ* = 1.6 K;
  - background w = +7.9·10⁵ J/mol; limiting value inside the window −2.2·10⁴ J/mol;
  - width 0.3 K or 1 K: α is negative inside, w changes sign twice (two simple poles, at the two zeros of α);
    width 3 K: no change of sign, only a dip of α.

## Prediction, as a test

A chromium sample whose transition is narrower than about 1.6 K has a window of negative expansion in which
the turn-part w is negative, bounded on both sides by a pole of w; by DW1 the sense of a small cycle is
reversed inside that window. A sample broadened beyond 1.6 K (strain, impurities, grain size) has no reversal.
The threshold is fixed by three measured numbers and nothing else.

This is unlike the horizons of BH1/BH2, whose pole is double and keeps the sign: there the heat reading
vanishes without changing sign; here the expansion does change sign.

## What is and is not new

T1 is the known relation between the anomalies at a first-order transition, written for w. The line's own
content is the reading (sense of the cycle = sign of w; two simple poles against one double pole) and the
threshold δ* as a statement about samples.
The heating/cooling difference of the latent heat, 0.13 J/mol per loop (4.2·10⁻⁴ J/(mol·K) of count), is the
kind of loop quantity PL2/LG1 describe; no formula of the line predicts its size yet. The paper attributes it
to relaxation, and it is within the stated ±0.10 uncertainties.

Put in: a uniform smearing of the step; the background values above. Not claimed: the measured α(T) curve of
chromium (not retrieved here); any value of D′ across the transition.
Open gate: retrieve measured α(T) and C_P(T) through T_N for samples of different widths; raw records need the
step, band limit and noise statement fixed in advance, since w has poles and D′ is a derivative.
