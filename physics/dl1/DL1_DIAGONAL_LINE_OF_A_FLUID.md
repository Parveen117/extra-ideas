# DL1 — The diagonal line of a fluid: where lost = seen on the diagram

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy); fluid numbers from reference equations of
state (CoolProp).

TD1: for the diagonal observer the lost part of the block of U(S, V) is m = (C_P − C_V)/C_P, and the diagonal is
m = ½. TD1 found its place on the critical isochore (2.68 – 2.81 T_c for the noble fluids) and could not derive it.
This stage follows the whole line.

Sources read before building: **TD1**, **DO1**, **QC2-T1** (the ensemble exp(−wᵀHw/2κ)), **NT-1**, **PH1** (argon
against krypton and xenon), **DG1**, **CL1/HB1** (1/√2).

## 1. The same line, five ways

```text
lost = seen                              m = ½
heats                                    C_P = 2·C_V
stiffness                                K_S = 2·K_T          (sound speed² = 2 × its isothermal value)
fluctuations (QC2's ensemble)            squared correlation of entropy and volume offsets = ½   (and of T and P offsets)
the block's own rates                    ratio (1 + 1/√2)/(1 − 1/√2) = 3 + 2√2
```

In the ensemble m is exactly the squared correlation between the two readings: on the diagonal half of the
variance of one is carried by the other; at the critical point all of it; in a gas of light all of it (TD1).

## 2. The count form

**E2.** For van der Waals with f readings the diagonal line is the spinodal multiplied by a ratio of counts,

```text
T_diag(ρ) = T_spinodal(ρ) × f/(f − 2) ,
```

so it peaks exactly at the critical density, at f/(f − 2) = 3 for single atoms.

## 3. Measured: the line is an arch

T/T_c of the diagonal at reduced density ρ/ρ_c:

```text
ρ/ρ_c       0.2     0.3     0.4     0.5     0.6     0.8     1.0     1.2     1.4     1.6     1.8     2.0     2.2     2.4
neon        1.438   1.840   2.148   2.384   2.562   2.774   2.813   2.726   2.554   2.322   2.034   1.686   1.281   0.857
argon       1.437   1.839   2.150   2.385   2.555   2.740   2.767   2.679   2.505   2.266   1.982   1.666   1.330   0.987
krypton     1.442   1.844   2.142   2.361   2.516   2.675   2.678   2.578   2.414   2.211   1.979   1.719   1.429   1.104
xenon       1.442   1.845   2.144   2.364   2.519   2.674   2.675   2.573   2.408   2.205   1.974   1.714   1.422   1.093
count form  1.176   1.640   2.028   2.344   2.592   2.904   3.000   2.916   2.688   2.352   1.944   1.500   1.056   0.648
```

```text
gas side (ρ ≤ 0.5 ρ_c)        the four fluids are one curve to 1.0 %
krypton and xenon             one curve along the whole line to 1.0 %
the peak                      2.69 T_c at 0.90 ρ_c (krypton, xenon) ;  2.77 at 0.94 (argon) ;  2.82 at 0.96 (neon)
the gas-side end              on the coexistence curve at T = 0.73 – 0.76 T_c, vapour density 0.07 – 0.08 ρ_c
the liquid at its triple point   m = 0.497 (neon), 0.508 (argon), 0.510 (krypton), 0.504 (xenon)
peak / Boyle temperature      0.987 (krypton, xenon), 1.023 (argon), 1.048 (neon)     count form: 8/9
```

The line is an arch over the critical point. It starts on the saturated vapour at about three quarters of T_c,
rises to about 2.7 T_c just below the critical density, and comes down through the dense fluid to the liquid at
its triple point, which lies on the diagonal to within 2%: there C_P = 2·C_V and half of the variance of entropy
is carried by volume.

## What this settles and what it does not

```text
TD1's open number             it is one point of an arch that the noble fluids share in reduced units
the count form                right shape, peak 10 % high, gas side 20 % low at 0.2 ρ_c
derived from the law          nothing here: the arch is set by how two atoms act on each other, which the line does not have
```

## What is not shown

- The numbers come from reference equations of state, not from raw measurements. The equations for krypton and
  xenon are recalled to share one functional form, so their 1% agreement may be partly built in; argon's and
  neon's are independent. Heat capacities of dense liquids in such equations are good to about a percent or two.
- The four fluids follow one law of corresponding states, so "four fluids agree" is close to one fact, not four.
  That the triple-point liquid sits at m = 0.50 ± 0.01, the vapour end near 0.75 T_c and the peak near the Boyle
  temperature are observations; no reason is given, and the literature was not searched for them.
- Other fluids are not on this arch (nitrogen 0.41 and methane 0.36 at their triple points).

## Claim boundary

```text
FIVE EQUIVALENT READINGS OF THE DIAGONAL ; m = SQUARED CORRELATION OF THE TWO OFFSETS        PROVED (symbolic)
VAN DER WAALS: DIAGONAL = SPINODAL × f/(f − 2), PEAK AT THE CRITICAL DENSITY                 PROVED
THE ARCH OF THE NOBLE FLUIDS ; ITS PEAK, ENDS AND SPREADS                                    MEASURED (reference equations)
TRIPLE-POINT LIQUID ON THE DIAGONAL TO 2 %                                                   OBSERVED ; NOT DERIVED
THE ARCH FROM THE LAW                                                                        NOT DERIVED — needs the action between two atoms
```

## Reproduce

```text
python dl1_diagonal_line_of_a_fluid.py        fluid numbers need CoolProp
python -m unittest test_dl1
```

## Later note (LD1, 9 October)

The far part of the action between two atoms is derived in LD1 from the floors of their modes (attraction, 1/r⁶,
the number ¾, a rule for unlike pairs). The arch still needs the size of the atom. Nothing above is changed.
