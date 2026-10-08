# CS1 — The memory of the stretch is the count of what falls: why the source is energy

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), with the functions of TP1.

SE1/MM1: on flat slices the energy component of the law is e₂(K), the unrecoverable memory of the stretch, and
with a source it equals k·u (EG1). EG1 and MM1 list as not derived why it is the *energy* of the source.

Sources read before building: **OR1** (a history with clock g enters through its count g·∫e⁰; the count is what
is stationary for motion), **CL2-T1** (count of a leg = g·τ), **CL1** (mass is the clock), **TP1**, **SE1**,
**MM1**, **EG1** (ratio law with a source; constant k not derived), **CO1-T2/T3** (the same variation for the
comoving frame of the research branch: varying N·a³[Q/2κ − …]; law value −6h²), **MC1-T3** (flux), **GR2** (source
= rest mass; fluxes add).

## Setting

Let the time step of the frame vary from place to place: e₀ = (1/N)(∂_t + u·∇), e_i = ∂_i. N = 1 is SE1's family.
A reading that falls with the frame counts N per unit t.

## Results

**C1.** With the time step in, TP1's law times the frame's volume is

```text
N·Q  =  −2·e₂(K) / N   +   a pure boundary term ,
```

for every fall velocity u.

**C2 — response to the time step.** At N = 1 the response of the law to a change of the time step is 2·e₂(K):
the stretch memory is what the law holds against the frame's own time.

**C3 — the source.** Let ν histories per volume fall with the frame, each with clock g. Their count per volume is
g·ν·N (OR1, CL2). Stationarity of (1/2k)·N·Q − g·ν·N in the time step gives

```text
e₂(K)  =  k · g · ν .
```

The unrecoverable memory of the stretch equals the count rate of what falls there. It is the energy because the
count is the one quantity of a history that is linear in the frame's time step, and g — the clock — is CL1's mass.

**C4.** The response to the fall velocity is SE1's (time, cut) law, curl curl u = 0; a source at rest in the
frame adds nothing to it.

**C5 — a centre made of histories.** With SE1's e₂ = (r·m)′/r² for radial fall, a ball of n₀ histories per volume
and radius a has

```text
inside     m = k·g·n₀·r²/3            fall speed β = h·r ,  h² = k·g·n₀/3 ,  Q = −6h²  (CO1's value)
outside    m = r_s / r ,              r_s = k·g·(number of histories) / 4π .
```

The flux r_s of MC1 is the total count rate of the source: clocks add, as GR2's fluxes do.

## What this settles

```text
EG1 / MM1: "why the energy"        the source enters through its count, and the stretch memory is the law's
                                   response to the time step that the count is measured in
MC1: the flux constant             r_s = k·g·N/4π : proportional to the number of clocks and to their rate
```

One constant k remains: the size of one unit of count in the frame law.

## What is put in

- The functional: frame law plus the count of the histories (OR1). CO1-T2 used the same variation for one frame;
  here it is done for every fall velocity.
- Histories that fall with the frame, and no pressure. Other sources are not treated.
- Varying the time step leaves SE1's family (MA1: one time is N·S = 1); the variation is taken at N = 1.

## What is not shown

- This is the known way a source enters (general knowledge). The line's part: the source term is the same count
  whose stationarity gives the motion (OR1), and the quantity it fixes is the memory e₂ of MM1.
- The constant k (G in units of count) is not derived. SC1 says no law in the block can give it.
- Only the energy component. Momentum and stress of a source are not built.

## Claim boundary

```text
N·Q = −2e₂(K)/N + BOUNDARY TERM, ANY FALL VELOCITY                              PROVED (symbolic)
RESPONSE TO THE TIME STEP AT N = 1 IS 2e₂(K)                                    PROVED
WITH THE COUNT OF FALLING HISTORIES: e₂(K) = k·g·ν                              PROVED
BALL OF HISTORIES: r_s = k·g·NUMBER/4π ; INSIDE h² = k·g·n₀/3 (CO1)             PROVED
THE CONSTANT k                                                                  NOT DERIVED (SC1)
SOURCES WITH PRESSURE OR MOTION ACROSS THE FRAME                                NOT BUILT
```

## Reproduce

```text
python cs1_count_is_the_source.py
python -m unittest test_cs1
```
