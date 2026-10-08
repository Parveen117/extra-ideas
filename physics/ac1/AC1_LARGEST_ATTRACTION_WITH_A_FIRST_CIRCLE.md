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
