# MC1 — 1/r from minimum cost: the least total variance between neighbouring shells

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues GR1.

GR1 reproduced the point-mass clock factor from the memory dictionary but
took the law m = r_s/r as input. Owner's suggestion: 1/r may come from a
minimum — of cost, of energy, of surface — and that minimum comes from
curvature.

Sources read before building: **QC4-T1** (the curvature of a shared
reading is minus the variance between the readings), **QC5-T3** (a share
moving in the scale s = log r; a cost with a minimum), **GR1-T1/T2**,
**SY1-S6**.

## 1. The cost

Take shells at radii r_k = b^k — equal steps of the scale s — in d space
dimensions, and a reading f_k on each. Charge every pair of neighbouring
cells the squared difference of their readings per unit distance (QC4:
the variance between two readings is what their sharing costs). A shell
has r^{d−1} cells and thickness proportional to r, so

```text
C[f] = Σ_k w_k (f_{k+1} − f_k)² ,        w_k = r_k^{d−2} .
```

Only d = 2 is scale-free. In d = 3 a difference costs more the farther
out it sits.

## 2. Results

**T1 (the minimizer).** With the two ends fixed, the minimum has a
constant flux Q = w_k(f_{k+1} − f_k) on every shell and is exactly

```text
f_k = A + B · r_k^{−(d−2)}        (d ≠ 2) ,          f_k = A + B·k = A + B′ log r_k        (d = 2) .
```

Exactly, at every node, for every ratio b — not approximately.

```text
d          1        2         3        4        5
law        r        log r     1/r      1/r²     1/r³
```

**T2 (it is the minimum).** For any change h vanishing at the ends,
C[f + h] − C[f] = C[h] > 0.

**T3 (gradient × area).** (f_{k+1} − f_k)/(r_{k+1} − r_k) · r_k^{d−1} is
the same on every shell: the difference between neighbours, per unit
distance, times the area of the shell, is conserved. What leaves through
one surface arrives at the next.

**T4 (three dimensions, reading = memory).** The least-cost memory
between m(r₀) = r_s/r₀ and its value far out is m = r_s/r at every node,
with the same flux wherever the chain is cut (4, 8, 16, 32 shells).
With GR1 this gives the clock factor N² = 1 − r_s/r without the 1/r law
as an input.

**T5 (which reading).** Spreading N itself at least cost gives
N = 1 − r_s/(2r), which differs from GR1's N at second order:
N² differs by exactly (r_s/2r)². The exact point-mass factor needs the
memory, not the clock factor, to be the quantity spread at least cost.

The four-dimensional gauge share of QC5 has the same tail: its deficit
1 − p falls as r^{−2} = r^{−(d−2)}.

## 3. Reading

- 1/r is what the least total variance between neighbouring readings
  looks like in three space dimensions. The exponent is d − 2 and comes
  from counting cells on a shell against the thickness of the shell.
- The conserved flux is the statement "what is lost is not recovered and
  not destroyed": the same amount crosses every surface. In GR1's
  dictionary that one constant is r_s — the mass.
- Gradient = flux / area: the lost information is spread over the
  surface. A sphere is where it is spread evenly.
- In two dimensions nothing sets a scale and the reading runs as log r;
  in three it must fall, because the same difference costs more farther
  out.

## 4. What went in

```text
MINIMIZER = POWER LAW r^{−(d−2)} ; UNIQUE ; GRADIENT × AREA CONSTANT          PROVED on scale-uniform shells
d = 3 AND READING = MEMORY ⇒ m = r_s/r ⇒ GR1's CLOCK FACTOR                  PROVED
THE COST IS QUADRATIC AND BETWEEN NEIGHBOURS ONLY                             ASSUMED (motivated by QC4, not derived from it)
THE READING SPREAD AT LEAST COST IS THE MEMORY                                ASSUMED (T5 shows the alternative fails at second order)
THREE SPACE DIMENSIONS                                                        INPUT
SPHERICAL SYMMETRY                                                            INPUT
THE VALUE OF THE FLUX IN TERMS OF MASS (2G/c²)                                INPUT
ANYTHING NOT STATIC                                                           NOT TOUCHED
```

The least-cost characterization of the 1/r law is known classically
(general knowledge). What this stage adds is that the cost is the
framework's own object — variance between neighbouring readings — and
that the reading which must be spread is the memory.

## 5. Reproduce

```text
python mc1_minimum_cost.py
python -m unittest test_mc1_exact
```

## Later note (CD1, 8 October)

"Three space dimensions" is listed above as an input. CD1 shows that, with this stage's least-cost memory in d
dimensions, bound histories close at first order only for d = 3. Given the closure rule of RC1 the input has a
reason. Nothing above is changed.

## Later note (PN1, 9 October)

"The reading spread at least cost is the memory" is listed above as assumed, with T5 showing the alternative
differs at second order. PN1: the alternative gives five sixths of Mercury's advance (35.8″ against the measured
42.98″) and a second-order number 3/2 against the measured 1 − (4.1 ± 7.8)·10⁻⁵. Measurement decides for the
memory. Nothing above is changed.
