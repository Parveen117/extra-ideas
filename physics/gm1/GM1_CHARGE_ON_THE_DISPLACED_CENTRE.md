# GM1 — A charge on the displaced centre: the same ι·a gives the turning and the magnetic moment, and their ratio is 2

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), exact. Uses the contracted curvature computed
and checked in SW2.

SW2: a turning centre is the centre at rest displaced by ι·a, and its memory is Re(r_s/R). This stage puts a
charge on the same centre.

Sources read before building: **SW2**, **EG1** (a radial light field reads (u; −u, u, u), its invariant part is
zero, and the memory r_s/r − q₂/r² has exactly that contracted curvature), **EM1-T1** (u² − S·S = ¼|F·F|²; what
every frame agrees on is F·F), **OB1** (F = E + ιB; U3: the four field equations are the four grades of D(F·C);
U6: E boosts, B turns), **SW1** (swirl constant and the turning content).

## Results

**C1 — two constants and no more.** For every a the invariant part (trace) of the contracted curvature of SW2's
frame is −P″/ρ². A source with no invariant part — a light field, by EG1 — leaves P″ = 0:

```text
P = r_s·r − q₂ ,        memory = Re( r_s / R ) − q₂ / ( R·R̄ ) .
```

This is EG1's memory r_s/r − q₂/r² read at the displaced distance.

**C2, C3 — the pattern.** With that memory the contracted curvature, seen from the frame that goes round the axis
with the centre at speed v = a·sin θ/√(r² + a²), is

```text
u · ( 1 ; −1 , 1 , 1 ) ,        u = q₂ / |R|⁴ ,
```

EG1's pattern at every a. In the falling frame itself it is that pattern boosted by v.

**C4, C5 — the light field.** On flat space the field of the displaced charge is F = E + ιB = −∇(q/R). It has no
source away from the ring (both grades of OB1-U3 vanish). Its frame-independent energy (EM1-T1) is ½|F·F| =
q²/2|R|⁴: the u of the pattern, with EG1's constant, at every a. Far away

```text
q/R  =  q/d  +  ι · q·a · cos ϑ / d²  +  …
```

a charge q in the cut part and a dipole of moment q·a in the turn part: E of a charge, B of a magnet.

**C6 — the ratio.** The same displacement in the memory gives r_s/R = r_s/d + ι·r_s·a·cos ϑ/d² + …: mass constant
r_s and swirl constant r_s·a (SW2-L1). So

```text
moment / charge   =   a   =   swirl constant / r_s .
```

With the dictionary of the thesis and SW1 (r_s = 2GM/c², r_s·a·c = 2GJ/c², moment = q·a·c):

```text
moment  =  ( q / M ) · J  =  2 × ( q / 2M ) · J .
```

A charged centre whose turning is a displacement along ι has twice the moment of a turning cloud of charge with
the same mass, charge and turning content.

**C7 — shells.** A light-like shell needs a² + q₂ ≤ r_s²/4.

## Numbers

```text
the electron ( J = ħ/2 )
    a = ħ / 2mc                      1.93·10⁻¹³ m
    r_s                              1.35·10⁻⁵⁷ m          a is 2.9·10⁴⁴ times r_s/2 : no light-like shell
    √q₂                              1.38·10⁻³⁶ m
    ratio in this stage              2
    measured                         2.002 319 304 36      (recalled; not re-read at source for this stage)
    measured / 2 − 1                 0.001 159 65
```

The first figure of the shortfall is the known α/2π = 0.001 161 with α ≈ 1/137 — the pure number the synthesis
lists as open. Here it shows up as the first departure from the bare displaced centre.

## What is put in

- SW2's frame family and CV1's law. EG1's statement that a light field has no invariant part.
- The constant between q₂ and q² (EG1's), and the dictionary r_s ↔ M, swirl ↔ J (thesis, SW1).
- For the numbers: ħ, the electron's mass and charge, G, c.

## What is not shown

- The ratio 2 for this field is known (general knowledge; no source re-read for this stage), and the measured
  electron value was already explained by the wave equation of the electron. Nothing new is predicted.
- The light field is checked on flat space (C4, C5) and through its pattern in the law (C2, C3). The field
  equations on the falling frame itself are not written out.
- The electron is not claimed to *be* this centre: at its scale a exceeds r_s by 44 orders, and the ring R = 0 is
  not treated. The statement is about the far field.
- Composite bodies (proton 5.59, neutron −3.83 in the same units; recalled) are not single displaced centres.
- The shortfall 0.00116 is not derived.

## Claim boundary

```text
NO INVARIANT PART ⇒ P = r_s r − q₂ , FOR EVERY a                           PROVED (symbolic)
PATTERN u(1; −1, 1, 1), u = q₂/|R|⁴, IN THE FRAME GOING ROUND AT a sin θ/√(r²+a²)   PROVED (symbolic)
F = −∇(q/R): NO SOURCE AWAY FROM THE RING ; CHARGE q + ι·DIPOLE q a        PROVED (symbolic, flat space)
MOMENT/CHARGE = a = SWIRL/r_s                                              PROVED
RATIO 2                                                                    FOLLOWS, with the dictionary r_s ↔ M, swirl ↔ J put in
THE ELECTRON IS SUCH A CENTRE                                              NOT CLAIMED
THE SHORTFALL α/2π                                                         NOT DERIVED — open gate
```

## Reproduce

```text
python gm1_charge_on_the_displaced_centre.py        needs sw2/sw2_contracted_curvature.json (written by SW2)
python -m unittest test_gm1
```
