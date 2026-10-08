# AC1 — The largest attraction that still has a first circle; the last stable circle is on the diagonal

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy).

BR1: the fine-structure number is the speed of the first recordable circle of a charge. The same question for a
mass, and what the law says about the two numbers together.

Sources read before building: **BR1**, **MO1-M4** (circles of a mass; last stable circle 3r_s), **OR1-D** (area of
the round pair), **OA1-T3** (count factor² = 1 − 3x/2), **CD1-D5** (special shells in d cuts), **DG1** (N = 1/√2 at
2r_s), **DO1** (the diagonal: seen = lost), **CL1/HB1** (the eight-mark value), **SC1**.

## Results

**A1 — a charge.** A circle of area l·κ has speed α/l, so it exists only for α < l. At α → l its radius and its
count factor go to zero. The first circle needs α < 1.

**A2 — a mass.** Over all circles of a mass the area is least at the last stable circle, L = √3·g·r_s. With
α_g := g·r_s/2κ, a circle of area l·κ exists only for

```text
α_g ≤ l / (2√3) = 0.2887·l .
```

**A3 — the diagonal.** On a free circle of a mass the count factor² is 1 − 3x/2. At the last stable circle it is
exactly ½: the circling reading sees as much as it loses. A free circle is stable exactly while it sees more than
it loses. In d cuts the value is (d − 2)/2: it is ½ only in three.

The eight-mark value 1/√2 marks two places in the field of a centre at rest: a held reading at 2r_s (DG1) and a
free circle at 3r_s.

**A4 — the two numbers.** For an electron and a proton: α = 7.30 × 10⁻³, α_g = 3.2 × 10⁻⁴², ratio 2.3 × 10³⁹.
Both are far below their limits. The law gives each number an upper limit and no relation between them.

## What this gives and what it does not

```text
limits            α < 1 for a charge ;  α_g ≤ 1/(2√3) for a mass — pure numbers from the law
the diagonal      the edge of stability of a free circle is seen = lost, in three cuts
a relation between α and α_g      REFUSED: none follows; they are two independent pure numbers (SC1)
```

## What is not shown

- A1 and A2 are the known limits for a point centre (general knowledge). The line's part is A3.
- A turning centre (SW2) changes A2 and A3; not done here.

## Claim boundary

```text
CHARGE: CIRCLE OF AREA lκ EXISTS ⇔ α < l                                             PROVED
MASS: LEAST AREA √3 g r_s AT THE LAST STABLE CIRCLE ; α_g ≤ l/(2√3)                  PROVED
LAST STABLE CIRCLE: COUNT FACTOR² = ½ ; (d − 2)/2 IN d CUTS                          PROVED
A RELATION BETWEEN THE TWO NUMBERS                                                  REFUSED
```

## Reproduce

```text
python ac1_largest_attraction_with_a_first_circle.py
python -m unittest test_ac1
```

## Later note (9 October): "two independent pure numbers" said too much

A4 above says the two numbers are independent. What was shown is less: nothing built so far relates them.
`ac1_two_numbers_one_gap.py` says exactly what stands between them.

**G1.** The turn term r_s·a of a centre does not contain its mass. So the charge number α = q₂/(r_s·a) knows
neither the mass nor the constant of the mass field: it is one number for every kind of matter.

**G2.** The mass number of a pair is not one number. It is the product of two mass counts,
n = r_s/2ℓ with ℓ² the record area in the units of the field: α_g = n₁·n₂.

**G3.** For one centre with the turn of one reading, α_g = r_s/4a, and α/α_g = 4q₂/r_s²: the large number is the
square of two lengths of the same centre, its charge length over half its mass length.

**G4.** The fall of a turning charged centre reaches 1 somewhere only if α_g ≥ (α + √(1 + α²))/2, a mass count of
0.71 at the measured α. The known kinds of matter are far on the other side.

**G5.** Electron 4.19 × 10⁻²³, proton 7.69 × 10⁻²⁰; their product is the 3.2 × 10⁻⁴² of A4.

So one thing is missing, not a relation between two constants: **what fixes the mass count of a kind of matter**.
The line has a centre with any r_s; it has no rule that says which centres exist. That is the entry "the kinds of
matter" of ONE_LAW.md, seen from here. The limits A1–A3 are unchanged.

```text
CHARGE NUMBER IS MASS-FREE ; MASS NUMBER = PRODUCT OF TWO MASS COUNTS               PROVED
α/α_g = (CHARGE LENGTH / HALF MASS LENGTH)²                                          PROVED
"INDEPENDENT"                                                                        WITHDRAWN: not related by anything built so far
WHAT FIXES A MASS COUNT                                                              OPEN — the one missing rule
```

```text
python ac1_two_numbers_one_gap.py
python -m unittest test_ac1_two_numbers
```

## Later note (FK1, MG1, 9 October)

FK1 does A3 for a turning centre: the last stable circle is at seen = ½ − (√6/12)a + …, on the diagonal only at
rest; and on every circle of a centre at rest seen − lost = (in-out/round)². MG1 finds the place of a mass count,
(ℓ/cell) × gap, and records that no certified gap gives the proton's value. Nothing above is changed.
