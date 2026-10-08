# TD1 — The thermodynamic block on the diagonal: C_V/C_P is F = R − D; a gas of matter is 2 : f; a gas of light has seen = lost

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy); fluid numbers from reference equations of
state (CoolProp).

DO1 put the observer on the diagonal. This stage does it on the diagram the framework started from: the response
block of U(S, V).

Sources read before building: Publications **NT-1** (H = [[a, b], [b, c]], Δ = ac − b²; C_V = T/a, C_P = Tc/Δ,
K_S = Vc, K_T = VΔ/a; Γ_cΓ_m = (C_V/C_P)(K_S/K_T) = 1), **DO1**, **NC1** (C_V/C_P = 1 − tanh²η: the clock factor
squared; memory = (C_P − C_V)/C_P), **RMG1** (eigenvalue ratio e^ℓ), **QC2-T1** (in the ensemble each reading
carries the unit κ: ⟨p_i w_i⟩ = κ), **OR1** (whole counts), **ZP1** (modes between walls), **PR1** (slow sector:
energy quadratic in momentum), **DG1**, **PH1** (argon against krypton and xenon), **SC1** (scale).

## Results

**T1 — the pure number of the block.** m = b²/(ac) does not change under any choice of units for S and V, and

```text
C_V / C_P  =  K_T / K_S  =  1 − m .
```

In units in which the two cuts read alike the block is √(ac)·[[1, t], [t, 1]] with t² = m: the diagonal observer.
For him seen R = 1, lost D = m, and F = R − D = 1 − m. So C_V/C_P is theorum/24's F on the diagonal, and
(C_P − C_V)/C_P is the lost part. The block's own rates are 1 ± t, ratio e^ℓ with tanh²(ℓ/2) = m.

**T2 — a gas of matter is a ratio of counts.** With whole counts fixed at fixed S, momentum ∝ 1/L and energy ∝
momentum², U = g(S)·V^(−2/d); with each reading carrying the same unit, T ∝ U. Then

```text
m = 2 / ( d + 2 ) ;        with f readings in all:   m = 2 / ( f + 2 ) ,   lost : seen = 2 : f .
```

In two cuts a gas of single points is exactly on the diagonal: m = ½.

**T3 — a gas of light has seen = lost.** With energy ∝ momentum and no count kept, U ∝ S^((d+1)/d)·V^(−1/d), and

```text
m = 1 ,     F = R − D = 0 ,     in every number of cuts .
```

It is so for every U of degree one in (S, V): no scale is left in the block. Pressure/energy density = 1/d and
entropy × temperature / energy = (d + 1)/d. The empty-space law of gravity (DO1) and a gas of light are the same
statement about their blocks: the diagonal observer loses what he sees.

**T4 — the diagonal of a real fluid.** For van der Waals with f readings, on the critical isochore
C_P/C_V = 1 + (2/f)/(1 − T_c/T): the diagonal m = ½ is at T/T_c = f/(f − 2), that is 3 for single atoms.

## Numbers

```text
memory at low density, 300 K           helium, neon 0.40000   argon 0.40001   krypton 0.40002   xenon 0.40004      count: 2/5
                                       nitrogen 0.2855   carbon monoxide 0.2853   oxygen 0.2830                    count: 2/7 = 0.2857
diagonal (m = ½) on the critical       krypton 2.678   xenon 2.675   argon 2.767   neon 2.813      (T/T_c)         van der Waals count form: 3
isochore                               nitrogen 1.739   oxygen 1.727                                                van der Waals count form: 5/3
memory at 1.001·T_c, critical density  argon 0.9979   krypton 0.9980   xenon 0.9979   carbon dioxide 0.9965        → 1 at the critical point
```

The low-density numbers are the counts. The place of the diagonal on the critical isochore is a pure number of
each fluid; krypton and xenon agree to 0.13%, argon is 3.4% above them (the same order as PH1's argon difference),
and the count form of T4 is 10% above. It is measured here, not derived.

## What is put in

- NT-1's dictionary. For T2: whole counts at fixed S, momentum ∝ 1/L, energy ∝ momentum², equal unit per reading
  (the sources above). For T3: energy ∝ momentum and no kept count.
- The reference equations of state for the fluid numbers.

## What is not shown

- T2 and T3 are the known heat-capacity ratios and the known gas of light (general knowledge). The line's part:
  the ratio is F of the diagonal observer, a gas of matter is the ratio 2 : f of counts, and a gas of light obeys
  the same seen = lost as the empty-space law.
- The diagonal's place for real fluids (2.68 – 2.81) is not derived; nor is the way m goes to 1 at the critical
  point.

## Claim boundary

```text
C_V/C_P = K_T/K_S = 1 − m ; m UNIT-FREE ; = F OF THE DIAGONAL OBSERVER                PROVED (symbolic)
MATTER GAS m = 2/(f + 2), LOST : SEEN = 2 : f                                         PROVED from the stated inputs
LIGHT GAS m = 1 IN EVERY NUMBER OF CUTS ; ANY DEGREE-ONE U                            PROVED
LOW-DENSITY GASES AT 2/5 AND 2/7                                                      MEASURED (reference equations)
DIAGONAL ON THE CRITICAL ISOCHORE: 2.68 – 2.81 FOR THE NOBLE FLUIDS                   MEASURED ; NOT DERIVED
```

## Reproduce

```text
python td1_thermodynamic_block_on_the_diagonal.py        fluid numbers need CoolProp
python -m unittest test_td1
```
