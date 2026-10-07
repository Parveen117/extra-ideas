# PH1 — Response holonomy of real fluids: the first numbers

Monty Dabas. 7 October 2026. Evidence packet (reference equations of state
through CoolProp 8.0; Python 3.12; not part of the stdlib CI).

RMG1 proved that a closed process carries a return angle

```text
Θ(γ) = −½ ∮ (cosh(ℓ/2) − 1) dφ ,        |Θ(γ)| ≤ ½ tanh(r*/4) L(γ) < ½ L(γ),
```

computed from the energy Hessian H = D²U along the cycle. This stage
evaluates Θ for real substances. Per unit mass, in (s, v):

```text
a = T / c_v ,      b = −(T / c_v)(∂P/∂T)_v ,      c = ρ² w²        (w = speed of sound),
```

so Θ needs only c_v, the pressure coefficient and the speed of sound along
the cycle — no derivatives of them.

## 1. Results

**Table 1. Three supercritical cycles** (ellipse in (T, ρ), 8000 points).

```text
fluid   cycle T [K], ρ [kg/m³]        Θ reference   Θ van der Waals   Θ ideal, c_v⁰(T)   Θ ideal, const c_v    L
CO₂     310–350, 200–700              −0.9260       −0.6606           −0.0071            0                    11.71
Ar      155–200, 250–800              −1.0728       −0.9781            0                 0                    12.39
N₂      130–170, 150–480              −0.9622       −0.7251           −0.00001           0                    11.08
```

**Table 2. One reduced cycle** T/T_c ∈ [1.03, 1.30], ρ/ρ_c ∈ [0.40, 1.50].

```text
fluid      Θ reference   Θ van der Waals   Θ ideal, c_v⁰(T)   |Θ| / (L/2)
Argon      −1.0810       −1.0760            0                 0.162
Krypton    −1.0282       −1.0760            0                 0.157
Xenon      −1.0254       −1.0760            0                 0.157
Methane    −0.9278       −0.7685           −0.0041            0.167
Nitrogen   −0.9771       −0.8334           −0.00001           0.165
CO₂        −0.9176       −0.7360           −0.0139            0.171
Water      −0.7648       −0.7335           −0.0095            0.162
```

## 2. What the numbers say

1. **Θ is of order one radian** for a near-critical cycle of a real fluid
   (53°–62° here). It is not a small correction.
2. **Θ is a unit-free number.** Rescaling the entropy and volume units by
   (37, 0.004) leaves Θ and L unchanged to all printed digits; r* is
   unit-dependent, so the unit-free inequality is |Θ| < ½L.
3. **Two independent routes agree** to 10⁻⁶: the loop formula of RMG1-T3 and
   direct transport dw = −½H⁻¹dH·w, whose return is a rotation in the
   H-orthonormal frame (orthogonality defect below 10⁻⁶).
4. **The calorically perfect gas is flat** (Θ = 0 to rounding), as RMG1
   states. An ideal gas with temperature-dependent c_v is **not** exactly
   flat: molecular vibration alone gives Θ = −0.007 to −0.014 for CO₂ and
   −0.004 for methane. Interaction gives the remaining ~99%.
5. **Θ separates equations of state.** Van der Waals with the fluid's own
   critical constants misses the reference value by 0.5% to 30%, depending
   on the fluid and the cycle (argon: 0.5% in Table 2, 9% in Table 1; CO₂:
   20–29%). The bound is far from saturated (|Θ|/(L/2) ≈ 0.16).
6. **Corresponding states.** Krypton and xenon agree to 0.3%. Argon sits
   5% away from them. Whether that is physics or the different quality of
   the three reference equations near T = 1.03 T_c is not decided here.

## 3. Claim boundary

```text
Θ FOR REAL FLUIDS FROM REFERENCE EQUATIONS OF STATE              COMPUTED (two routes, unit-blind)
PERFECT GAS FLAT; VIBRATIONAL c_v(T) GIVES SMALL NONZERO Θ        COMPUTED
|Θ| < ½ L                                                        HOLDS (it is a theorem; data cannot violate it)
Θ MEASURED IN A LABORATORY                                       NOT DONE
A PHYSICAL PROCESS WHOSE OUTCOME IS Θ                            NOT IDENTIFIED
ANY CLASSICAL EQUATION MODIFIED                                  NO
```

Θ here is a derived quantity of equilibrium thermodynamics: every input is
classical. What is new is the quantity itself and its numbers. The open
physics question is the fifth line: which experiment returns Θ directly.

## 4. Reproduce

```text
pip install CoolProp numpy
python ph1_real_fluid_holonomy.py 8000
python ph1_corresponding_states.py 8000
```
