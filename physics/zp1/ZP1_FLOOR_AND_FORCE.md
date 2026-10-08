# ZP1 — An indeterminate sum resolved: the floor of every mode, and the force it leaves

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Every
coefficient exact; floats only in the illustration.

The owner's programme for what remains: take a problem; it is either
solved or shown indeterminate, and indeterminacy is what the Bindu–Lopa
operators are for. First case: the sum of the floors of all modes
between two walls — a sum that is infinite term by term.

Sources read before building: **LC1-C4** and **WQ1-W3** (every mode has
the floor ½κω: what stays seen when nothing is lost), **PS1-T4** (Σ n wⁿ
in the two Exp sectors; the constant −1/12 is the same in both, the pole
changes sign), **EM1**, **FR1**.

## Results

**Z1.** Exactly, Σ_{n≤N} n wⁿ = w(1 − (N+1)w^N + N w^{N+1})/(1−w)².

**Z2 — the two sectors (PS1-T4, re-certified by exact series).**

```text
scale sector     Σ n e^{−nx}            =   1/x²  − 1/12 + x²/240 − …
phase sector     −1 / 4 sin²(φ/2)       =  −1/φ²  − 1/12 − φ²/240 − …
```

The constant is the same; the pole and the quadratic term change sign.
What depends on how the sum is cut off is exactly what differs between
the sectors; what does not, is common to both.

**Z3 — one cut.** A line of length L, modes ω_n = nπc/L, each with
floor ½κω_n:

```text
E = π κ c L / 2ε²   −   π κ c / 24 L   + O(ε²) .
```

The indeterminate part is proportional to L. For a wall at position a
inside a fixed length it does not depend on a at all: it is silent in
the force. What remains is E = −πκc/24L, and the wall is pulled toward
the nearer end.

**Z4 — three cuts.** Two plates a distance d apart:

```text
E / area = 3 d κ c / π² ε⁴  −  κ c / 2π ε³  −  π² κ c / 720 d³  + O(ε²) .
```

The first term is proportional to d and silent in the force by the same
argument; the second does not depend on d. What remains:

```text
force per area = − π² κ c / 240 d⁴ .
```

## Numbers (κ = ħ)

```text
plates 100 nm apart     13.0 Pa
plates 1 μm apart       1.3 mPa
```

This is the attraction between uncharged conducting plates. It has been
measured, and agrees with this value (general knowledge; the
measurements were not re-examined here).

## Reading

- The floor ½ of LC1 — what stays seen when nothing is lost — is not
  bookkeeping. Summed over the modes between two walls it pulls them
  together with a force that has been measured.
- The indeterminate sum splits by itself into two kinds of term: those
  that grow with the size of the region, which no wall position can
  feel, and one that does not depend on the cut-off at all. PS1's
  sector-blind constant is the second kind. "Indeterminate" here meant
  "silent", and what was not silent was determined.
- The numbers 24, 240, 720 come from one series, 1/(e^y − 1).

## What is not shown

- This is the known force between ideal plates. The route is the
  framework's (floor from LC1, separation from PS1); nothing new is
  predicted.
- The modes between the plates (ω_n = nπc/d, two kinds of wave, ideal
  walls) are taken as given, not derived from OB1's law with walls.
- Real plates, temperature, and other shapes are not treated.
- That every indeterminate quantity in physics splits this way is not
  claimed; this is one case.

## Claim boundary

```text
FINITE SUM AND ITS LIMIT                                              PROVED exact
SERIES IN BOTH SECTORS; CONSTANT −1/12 SECTOR-BLIND                   PROVED exact
ONE CUT: E = πκcL/2ε² − πκc/24L ; POLE SILENT IN THE FORCE            PROVED
THREE CUTS: E/area FINITE PART −π²κc/720d³ ; FORCE −π²κc/240d⁴        PROVED
THE MODES BETWEEN IDEAL WALLS                                        TAKEN as given
AGREEMENT WITH MEASUREMENT                                           general knowledge, not re-examined
```

## Reproduce

```text
python zp1_floor_and_force.py
python -m unittest test_zp1_exact
```
