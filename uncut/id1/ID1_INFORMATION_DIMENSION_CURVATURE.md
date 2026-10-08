# ID1 — information, dimension and curvature

Stage of the uncut line. No measurement. Exact (sympy): `id1_information_dimension_curvature.py`, `test_id1.py` (5 tests).

The owner's statement: gravity may be a coupling of information and curvature. On a straight line the
dimension is one and the information is unbounded. As the dimension grows the information is distributed among
the phases; in three dimensions it is distributed further.

## Sources used (read, unchanged)

extra-ideas `physics/mc1` T1 (least-cost spreading between shells: r^−(d−2), log r, r), `cv1` (what no frame change
removes: r_s/r³ × (1, −½, −½)), `mo1` (fall speed² = memory m), `in1` T2 (a single reading: n² = Σ r_i²), `gb1`;
`uncut/ci1` (exp(−2I) = 1 − m), `up1` T6.

## Statements

**T1 (one cut among D).** For a single reading, Σ (r_i/n)² = 1, so one cut holds on average 1/D of it:
all of it on a line, half with two cuts, a third with three.

**T2 (spreading over directions).** A point source in d space dimensions has memory m with the same flux through
every shell:

| d | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| m | r | log r | 1/r | 1/r² | 1/r³ |
| far value | unbounded | unbounded | 0 | 0 | 0 |

On a line nothing is diluted: the shared information grows without bound with distance. Three is the first
dimension in which it goes to zero far away — the first in which a flat far region exists at all.

**T3 (curvature is the sharing among directions).** The part of the curvature that no change of frame removes
is half the Hessian of m (CV1's r_s/r³ × (1, −½, −½) in three dimensions). Along and across the radius:

| d | radial : each transverse |
|---|---|
| 1 | none — a line has no direction to share with, and no curvature away from the source |
| 2 | 1 : −1 |
| 3 | 1 : −½ |
| 4 | 1 : −⅓ |
| d | 1 : −1/(d − 1) |

One radial direction is balanced by the d − 1 transverse ones together. The ½ of three dimensions is one
radial direction shared between two transverse ones.

**T4 (the coupling, written for I).** With exp(−2I) = 1 − m:

  ½ Hessian(m) = exp(−2I) · ( Hessian(I) − 2 ∇I ∇I ),

and the free law (m without source) is

  ΔI = 2 |∇I|².

The shared information is not free to spread by itself: its spreading is fed by its own gradient. Far away, to
first order, ΔI = 0 and the curvature is the Hessian of I. Memories of two sources add; their I do not.

## Reading

The owner's sentence in this form: the information a point source shares with its surroundings is spread over
the directions available. On a line there are none, so it is never thinned and there is no curvature. Each
added dimension adds a transverse direction to share with; the curvature is exactly the record of that
sharing — one part outward against d − 1 parts sideways. Gravity as "a coupling of information and curvature"
is T4: curvature is the second spread of I, corrected by I's own gradient.

## What is put in, what is not claimed

* The spreading law is MC1's (least cost between neighbouring shells) and the curvature is CV1's, both for the
  fall frame with flat slices. T4 is those two rewritten with CI1's I; it is an exact change of variable,
  not a new law.
* Lineage: T2 and T3 are the classical facts about a point source in d dimensions.
* Why space has three directions is not derived. T2 says only what three is the first to allow.
* The strength of the source (r_s) and every constant are inputs. No number is predicted.

## Open gates

1. T1 against T3: one cut holds 1/D of a reading, one transverse direction carries 1/(d − 1) of the radial
   curvature. Whether these are one statement.
2. The law ΔI = 2|∇I|² with a source: what stands on the right for a body, in units of I.
3. Three cuts and three directions: FR1 has "no fourth cut"; whether that is why d = 3.
