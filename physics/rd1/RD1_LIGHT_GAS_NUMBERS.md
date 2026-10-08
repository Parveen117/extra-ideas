# RD1 — The pure numbers of a gas of light, from whole counts

Monty Dabas. 9 October 2026. Python 3.12. Exact rational arithmetic, sympy and mpmath.

TD1: a gas of light has seen = lost. This stage counts it.

Sources read before building: **WQ1** (a recordable mode has E = n·κ·ω), **OR1** (whole counts), **QC2** (the
ensemble; the thermal unit), **ZP1 / ZP2** (modes between walls; sums with explicit tails), **EM1** (two readings
per direction), **TD1**, RKF **theorum/28** §4 (a sum over all terms needs an explicit tail). Measured:
**arXiv:astro-ph/9605054**, abstract, read at source.

## Results

**L1 — the mean count.** Whole counts n with weights xⁿ, x = exp(−κω/θ), have mean x/(1 − x). Exact, with the tail
of the sum written out.

**L2 — two sums.** ∫y³/(e^y − 1) = Σ 6/k⁴ and ∫y²/(e^y − 1) = Σ 2/k³, each bracketed for every K with tails below
2/K³ and 1/K²: π⁴/15 and 2ζ(3).

**L3 — the numbers** (three cuts, two readings per direction; θ the thermal unit):

```text
energy density × (κc)³ / θ⁴          π²/15 = 0.6580
flux number                          π²/60 = 0.1645
count density × (κc/θ)³              2ζ(3)/π² = 0.2436
energy per count / θ                 π⁴/30ζ(3) = 2.7012
entropy per count                    2π⁴/45ζ(3) = 3.6016
```

**L4 — where the spectrum peaks.** At κω/θ = 2.8214 per unit rate and 4.9651 per unit length.

**L5 — ratios from TD1.** Pressure/energy density = 1/d; entropy·θ/energy = (d + 1)/d; energy density ∝ θ^(d+1):
the fourth power in three cuts.

## Measured

The light that fills the sky: deviations from this shape below 50 parts per million of the peak; count potential
over θ below 9 × 10⁻⁵. The second is TD1's "no count kept": the block has no third variable.

## What is put in

- Whole counts (WQ1), the ensemble's weights (QC2), modes between walls (ZP1), two readings per direction (EM1).

## What is not shown

- These are the known numbers of thermal light (general knowledge). The line's part is that they follow from its
  whole counts with explicit tails, and that "no kept count" is TD1's seen = lost.
- The thermal unit θ in terms of κ is not fixed here (QC2: k_BT in the energy representation).

## Claim boundary

```text
MEAN COUNT x/(1 − x) ; SUMS π⁴/15 AND 2ζ(3) WITH EXPLICIT TAILS                PROVED (exact, every K)
THE FIVE NUMBERS AND THE TWO PEAKS                                             PROVED
1/d , (d + 1)/d , POWER d + 1                                                  PROVED
SHAPE OF THE SKY'S LIGHT                                                       MEASURED to 50 parts per million
```

## Reproduce

```text
python rd1_light_gas_numbers.py
python -m unittest test_rd1
```

## Later note (QD1, 9 October)

The spread of the count of a mode is n + n² = S of T24-6.1 with R = n, D = n²; the diagonal R = D is the mode with
quantum kT·ln 2 (39.4 GHz for the sky). Nothing above is changed.
