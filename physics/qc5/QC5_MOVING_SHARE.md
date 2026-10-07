# QC5 — A share that moves in scale: the winding reading, the logistic law and the cost of changing reading

Monty Dabas. 7 October 2026. Python 3.12, exact rational quaternion
arithmetic only. Continues QC4.

QC4 ended with the Yang–Mills gate: θ is not a share between dual
readings. This stage finds where the share law **does** live in
Yang–Mills — between two flat readings of the gauge field itself — and
what it costs to move from one to the other.

Sources read before building: GE2-T2 (typed record lift, M = A†QA),
PH3 (share law), QC4-T1 (curvature = −variance), YM100 (flat readings
with winding; the centre orbit), the quaternion carrier of YM-F1 and GE2.

## 1. Two flat readings on the quaternion carrier

For a quaternion point x ≠ 0 put g = x/|x| and L = g⁻¹dg = Im(x̄ dx)/|x|².
The trivial reading A = 0 and the winding reading A = L are both flat.
Share them with a weight that depends on the scale r = |x|:

```text
A = p(r) · L .
```

**T1 (gradient term + variance term).**

```text
F = dp ∧ L − p(1−p) L ∧ L .
```

The second term is QC4-T1 — the variance of a two-way choice. The first
is new: it appears because the share moves.

**T2 (velocity and variance).** With s = log r, velocity p_s = r·dp/dr and
V = 2p(1−p), at every point

```text
action density            = 3 ( p_s² + V² ) / r⁴ ,
topological density       = ± 6 p_s V / r⁴ ,
F = ± its dual   ⇔   p_s = ± V .
```

A pure reading (p ≡ 0 or 1) has no field. The equal share p ≡ ½ is pure
variance, density 3/(4r⁴), no velocity.

## 2. The logistic law

**T3 (cost of changing reading).** The reduced action is

```text
3 ∫ ( p_s² + V² ) ds  =  3 ∫ ( p_s − V )² ds  +  2 [ 3p² − 2p³ ] ,
```

so moving the share from 0 to 1 costs at least **2**, with equality iff

```text
dp/ds = 2 p (1 − p)        ⇔        p = r² / (r² + λ²) .
```

The share obeys the logistic law in the scale: its velocity equals its
own variance. The scale λ is free — the cost is 2 for every λ. Half the
change costs half: 3p² − 2p³ = ½ at p = ½.

**T4 (it is a record lift).** The logistic share is the two-record lift
u = (x, λ): with N = |x|² + λ² and Q = 1 − uu†/N,

```text
F = du† Q du / N = λ² ( ē_μ e_ν − ē_ν e_μ ) / N² ,        action density 24 λ⁴ / N⁴ .
```

The field is the memory form of the complement cut, GE2's M = A†QA.

## 3. Reading

```text
pure reading (p = 0 or 1)        no field, no cost                — a classical vacuum
equal share (p = ½)              pure variance, half the charge   — PH3's native transport, as a gauge field
logistic share (0 → 1)           velocity = variance, cost 2      — the cheapest way to change reading
```

- The share law of PH3 does live in Yang–Mills: between the windings of
  the flat readings, not between the electric and magnetic terms.
- Changing reading has a fixed minimal cost, independent of scale. In
  the classical normalization the reduced 2 is 4π² in the native
  coefficient norm, 8π² in the trace norm.
- A single reading cannot see this sector at any order: every term above
  vanishes at p = 0 and p = 1. That is the content of YM99's refusal and
  of YM100's centre orbit, now with a cost attached.
- The configurations of T3 and of the equal share are known classically
  (the one-instanton solution and the half-charge configuration). What
  is new here is their place: they are the moving and the frozen forms of
  the same share law that gives the return angle of a fluid.

## 4. Certificate

At three rational quaternion points and four share laws (logistic, a
steeper share, equal share, winding only): the full field F_μν by exact
differentiation, its dual, the action and topological densities against
T2, the square completion, and self-duality exactly when velocity =
variance. T4 at six (point, scale) pairs: field = memory form = closed
form, density 24λ⁴/N⁴. Reduced action and half-share value as exact
rationals. Six tests; a corrupted field, a wrong share law and a broken
dagger are rejected.

## 5. Claim boundary

```text
F = dp∧L − p(1−p)L∧L ; DENSITIES ; SELF-DUAL ⇔ VELOCITY = VARIANCE    PROVED
COST OF CHANGING READING ≥ 2, LOGISTIC LAW, SCALE FREE                PROVED (radial two-reading sector)
LOGISTIC SHARE = MEMORY FORM OF A TWO-RECORD LIFT                     PROVED
THIS SECTOR INVISIBLE TO A SINGLE READING                             PROVED (all terms vanish at p = 0, 1)
MASS GAP FROM THIS SECTOR                                             NO — the sum over scales λ is not controlled
EVERY GAUGE FIELD IS A RECORD LIFT                                    NOT CLAIMED (named classical statement, unchecked)
```

The cost is scale-free, so nothing here selects a scale; a mass is a
scale. The open gate is exactly that: what fixes λ. On the thermal side
the unit came from the ensemble (QC2); here it has to come from the sum
over shares at all scales.

## 6. Reproduce

```text
python qc5_moving_share.py
python -m unittest test_qc5_exact
```
