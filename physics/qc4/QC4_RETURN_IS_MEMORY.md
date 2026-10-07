# QC4 — Return is memory: one identity behind the thermodynamic return, the record memory and the Yang–Mills action

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues QC3.

QC3 found, for one arrow, M = −O²: the memory a symmetric record loses is
the squared odd part. This stage shows that the same identity, stated for
readings instead of turns, is the curvature that PH1–PH3 measured, and that
stated for a loop arrow it is the Yang–Mills action.

Sources read before building: GE2-T2 (M = I − S†S = Σ p_s‖U_s f − S f‖²),
PH3 (two flat readings, share law), QC3-T1.

## 1. The identity

Let A_s be flat readings of a probe (each returns exactly on every cycle)
with weights p_s, Σp_s = 1, and Ā = Σ p_s A_s their mean.

**T1 (curvature is minus the variance of the readings).**

```text
F(Ā) = − Σ_{s<t} p_s p_t [ A_s − A_t , A_s − A_t ]        (wedge commutator of the differences)
```

*Proof.* Each dA_s = −A_s∧A_s, so F(Ā) = Ā∧Ā − Σp_s A_s∧A_s, which is
minus the weighted variance; expand with Σp_s = 1. ∎

- Two readings with weights (1−λ, λ): the variance is λ(1−λ) — PH3's share
  law. The 4λ(1−λ) parabola is the variance of a two-way choice.
- A single reading has no variance and no return.
- GE2's record memory is the same object for arrows: M = Σ p_s‖U_s − S‖²
  (T2 below), and for a symmetric pair it is −O² (QC3).

So the return angle Θ of a real fluid is the memory of not having recorded
which reading the probe was in.

## 2. The compass law

Four flat readings of the probe (ds, dv), one per corner of the T–V–S–P
compass: U freezes (δs, δv); G freezes (δT, δP); F freezes (δT, δv);
H freezes (δs, δP). Returns are given as fractions of the native one.

**T3.**

```text
diagonal U–G, equal share        1          (the native transport)
diagonal F–H, any share          0          (curvature present, but no rotation)
four corners, equal share        1/2
the four edges, summed           1
```

The return lives on the U–G diagonal. The mixed diagonal F–H carries
curvature of the other kind only (no rotation part), at every share.

*Proof.* Direct expansion of T1 with the four explicit readings; the
identities hold for arbitrary response fields (checked symbolically) and
are certified exactly below. ∎

## 3. The loop arrow: Wilson density and field strength

Let U be the arrow of a closed loop (a plaquette) and read it without its
orientation: the symmetric record {U, U†}.

**T4.**

```text
E = (U + U†)/2 ,   O = (U − U†)/2 ,
record defect   1 − E            — its normalized trace is the Wilson density 1 − W,
record memory   M = I − E² = −O² — the squared field strength,
M = (1 − E)(1 + E) ,            M / (2(1 − E)) → 1 as the loop arrow → 1 .
```

Both are blind to the orientation and unchanged by a gauge turn at the
base point. The Wilson action is the record defect of the loop arrow; the
field-strength-squared action is its record memory; they differ by the
factor (1 + E) and agree in the small-loop limit.

With QC3-T4 for the links (record defect rate = Casimir), the Wilson
Hamiltonian is a sum of two record defects of native arrows:

```text
A_θ  =  [ defect rate of the link arrows under symmetric turn records ]
      + θ · [ defect of the loop arrows under the orientation record ] .
```

## 4. Where Yang–Mills differs from the share law

In PH3 the two readings are dual readings of one pair, and the share
λ ∈ [0, 1] is symmetric. In A_θ the two terms are record defects at two
nested levels — the loop arrow is a product of link arrows — and θ is the
ratio of two ledgers, not a share between dual readings. This is why the
λ ↔ 1−λ symmetry has no Yang–Mills counterpart, and it names the native
variable to carry forward: the gap per unit of loop ledger. Not a theorem;
the next gate.

## 5. Certificate

T1 at three points of a quartic convex energy with four flat readings
(flatness of each checked from derivative tables), five weightings, both
sides computed independently; T3 exactly at those points (and at random
rational points during development). T2: two- and three-record circular
examples. T4: an SU(2) loop of four rational quaternion links at three
turn sizes (E scalar, O vector, M = (1−E)(1+E), gauge and orientation
blindness, ratio 0.924 → 0.9938 → 0.99990) and an SO(3) loop where E is
not scalar. Six tests; weights not summing to one, a curved reading and a
non-isometric link are rejected.

## 6. Claim boundary

```text
CURVATURE OF A MEAN OF FLAT READINGS = −VARIANCE                     PROVED
SHARE LAW = TWO-WAY VARIANCE; COMPASS LAW (1, 0, 1/2, edges 1)       PROVED
RECORD MEMORY = VARIANCE OF RECORDS                                  PROVED (GE2-T2 restated, exact instances)
WILSON DENSITY = LOOP RECORD DEFECT; FIELD STRENGTH² = LOOP MEMORY   PROVED
WILSON HAMILTONIAN = TWO RECORD DEFECTS AT TWO LEVELS                STATED from QC3-T4 and T4; free terms only
MASS GAP AT WEAK COUPLING, CONTINUUM                                 NOT TOUCHED
k_BT IN TERMS OF ħ                                                   NO
```

## 7. Reproduce

```text
python qc4_return_is_memory.py
python -m unittest test_qc4_exact
```
