# PH3 — The fluid performs the return when a probe is shared between two flat readings

Monty Dabas. 7 October 2026. Python 3.12. Continues PH2.

PH2 left one gate: a process in which the fluid itself, not a programmed
chain, produces the return angle Θ. This stage gives one.

## 1. Two readings of a probe

Take a reference sample and a probe sample of the same fluid and mass,
differing by small offsets w = (δs, δv), with intensive offsets
(δT, −δP) = H·w. Carry the reference around a cycle. The probe can follow
in two ways:

```text
E reading   extensive offsets frozen: (δs, δv) constant            transport  dw = 0
M reading   intensive offsets frozen: (δT, δP) constant            transport  dw = −H⁻¹dH · w
```

**T1 (both readings are flat).** Around every closed cycle the E reading
returns w exactly and the M reading returns w exactly: neither has a
return angle.

**T2 (the native transport is their equal share).** dw = −½H⁻¹dH·w, the
transport of NT-3 and RMG1, is the mean of the two. This is the
metric-compatible mean of the two dual flat connections of NSB1.

**T3 (share law).** If the probe spends the fraction λ of every step in
the M reading, the transport is dw = −λH⁻¹dH·w and its curvature is

```text
F_λ = −λ(1−λ) [X_s, X_v] = 4λ(1−λ) · F_½ ,        X_i = H⁻¹ ∂_i H .
```

Zero at λ = 0 and λ = 1, largest at equal sharing, symmetric under
λ ↔ 1−λ.

**T4 (split legs).** Reading entropy moves in M and volume moves in E has
curvature F₁; the opposite assignment has F₂; and F₁ + F₂ = 4F_½.

*Proof.* H⁻¹dH is the derivative of a function of the state, so
∂_sX_v − ∂_vX_s + [X_s, X_v] = 0, which is T1 for M; E is constant. For
A = λX: dA + A∧A = λ(∂_sX_v − ∂_vX_s) + λ²[X_s, X_v] = (λ² − λ)[X_s, X_v].
For T4, F₁ = −∂_vX_s and F₂ = ∂_sX_v, whose sum is −[X_s, X_v]. ∎

*Certificate.* T1, T3 (five shares) and T4 as exact rational identities at
three points of a quartic convex energy; the M reading multiplied around
a rational pentagon returns the identity matrix exactly; a field whose
mixed derivative is not a derivative is rejected.

## 2. The protocol

Twin cells. In an E step both cells are sealed and receive the same
entropy and volume increments. In an M step the probe is held at fixed
temperature and pressure offsets from the reference. Alternate the two
along the cycle.

## 3. Numbers: CO₂, the cycle of PH1

Alternating E and M on successive segments:

```text
segments      angle         distance from a pure rotation
      8      −1.1522        1.5
     32      −0.9529        0.41
    128      −0.9275        0.079
    512      −0.9261        0.020
   2048      −0.92604       0.0052
  16384      −0.926030      0.0006          (PH1: Θ = −0.926030)
```

Share law:

```text
share λ     4λ(1−λ)     small cycle (5% size)     full cycle
 0           0           0                         0
 0.10        0.36        0.3603                    0.423
 0.25        0.75        0.7502                    0.796
 0.50        1           1                         1
 0.75        0.75        0.7502                    0.796
 0.90        0.36        0.3603                    0.423
 1           0           0                         0
```

The law is exact for the curvature and holds to 10⁻³ on the small cycle.
On the full cycle the λ ↔ 1−λ symmetry survives to all printed digits but
the values exceed the parabola: for λ ≠ ½ the return is not a pure
rotation (distance 0.2).

A concrete probe at the start state (350 K, 450 kg/m³), equal sharing:
an entropy offset of 1 J/(kg·K) with no volume offset returns as
δs = −0.607 J/(kg·K), δv = −7.24·10⁻⁶ m³/kg. The second-order energy
½wᵀHw is the same before and after (0.187284 to 10⁻⁷).

## 4. Claim boundary

```text
E AND M READINGS FLAT; NATIVE TRANSPORT = EQUAL SHARE            PROVED
SHARE LAW 4λ(1−λ) FOR THE CURVATURE; SPLIT-LEG SUM RULE          PROVED
ALTERNATING PROTOCOL → Θ FOR CO₂                                 COMPUTED (convergence 1/N)
SHARE LAW ON FINITE CYCLES                                       FAILS as an equality for λ ≠ ½; symmetry holds
TWIN-CELL PROTOCOL BUILT OR ITS FEASIBILITY ASSESSED             NO
ANY CLASSICAL EQUATION MODIFIED                                  NO
```

The protocol is linear-response and quasi-static: small offsets, no
dissipation. E steps on legs with heat exchange need equal entropy
increments in both cells, which is a calorimetric control problem this
stage does not address.

## 5. Reproduce

```text
pip install CoolProp numpy
python ph3_two_readings_transport.py
python -m unittest test_ph3_exact
```
