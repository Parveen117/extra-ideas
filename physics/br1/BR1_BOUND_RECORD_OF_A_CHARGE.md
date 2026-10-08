# BR1 — The bound record of a charge: the fine-structure number is the speed of the first circle, and it sets the size of an atom

Monty Dabas. 9 October 2026. Python 3.12. Symbolic algebra (sympy) and mpmath.

FS1 placed the number 1/137 on one centre: α = q₂/(r_s·a). LD1 needed the size of an atom. This stage puts a
charge on a recordable history around a charged centre and reads where α appears.

Sources read before building: **OB1-U6** (the field changes a reading (n; r) by q(E·r ; nE + r×B)), **PR2-T2** (the
boosted mass is (n; r) = γ·g·(1; v)), **OR1** and **WQ1-W2** (a record is kept when the area of a pair is whole:
2πκ·n), **RC1** (closure and the turn left over), **MO1-M5**, **GR1** (N² + m = 1), **MS1** (speed² + memory = 1),
**GM1** (displacement a = κ/2g for turning κ/2), **FS1**, **LD1**. Units c = 1; K is the strength of the attraction
between the two contents, α := K/κ.

## Results

**B1 — every bound history.** With E = γg − K/r and L conserved, the history obeys

```text
u″ + ( 1 − K²/L² )·u = E·K / L² ,      u = 1/r :      in-out rate / round rate = √( 1 − K²/L² )   exactly, at every amplitude.
```

**B2 — circles.** On a circle the speed is K/L and the conserved energy is the rest count times the count factor.
With the round area whole, L = l·κ:

```text
speed = α / l ,        memory = α² / l² ,        E = g·√( 1 − α²/l² ) ,        radius = l²·( κ / gα )·√( 1 − α²/l² ) .
```

The fine-structure number is the speed of the first recordable circle in units of the cone speed; its square is
that circle's memory (MS1, GR1). Nothing is approximate in these lines.

**B3 — three lengths.** One charge has three lengths, each α times the next:

```text
charge length  K/g  =  α·κ/g        turn length  κ/g  ( = 2a of GM1 )        record length  κ/gα  ( = 2a/α ) .
```

The size of the atom is the electron's own turn length divided by α. The binding of the first circle is
1 − √(1 − α²) = α²/2 + α⁴/8 of the rest count.

**B4 — the in-out pair.** Its area is J_r = 2π[E·K/√(g² − E²) − √(L² − K²)] (checked by quadrature), and its slope in
L is −2π × (round rate / in-out rate). With both areas whole,

```text
E = g / √( 1 + α² / ( n_r + √(l² − α²) )² ) .
```

**B5.** In the second shell the two records (n_r, l) = (1, 1) and (0, 2) differ by α²/16 of the binding scale.

**B6 — the turn left over.** A bound charge never closes: per circuit it leaves 2π(1/√(1 − α²/l²) − 1) = π·α²/l².
At the same K/L the field of a mass leaves six times as much (MO1-M5).

## Numbers

```text
speed of the first circle / cone speed        0.0072974                      = α
memory of the first circle                    5.3251 × 10⁻⁵                  = α²
binding / rest count                          2.6626 × 10⁻⁵                  13.6057 eV on 510 999 eV
lengths                                       2.818 × 10⁻¹⁵ ,  3.862 × 10⁻¹³ ,  5.2918 × 10⁻¹¹ m
turn left over, first circle                  1.673 × 10⁻⁴ radian per circuit
split of the second shell / binding scale     line 3.3282 × 10⁻⁶ ;  measured 3.3342 × 10⁻⁶   (+0.18 %)
rates between circles (counts only)           (1→2) : (2→3) = 27 : 5 ;   (2→4) : (2→3) = 27 : 20 ;   (2→5) : (2→4) = 28 : 25
measured wavelength ratios                    1.35001 against 27/20 ;  1.11999 against 28/25
```

## What this gives and what it does not

```text
where α sits in a record        speed of the first circle; √memory; the ratio of the three lengths; the left-over turn π α²
the size of an atom             the turn length over α — for one charge around one centre
ratios of the spectrum          counts only: 27/5, 27/20, 28/25, …
the value of α                  not derived: it is K/κ, the strength of the attraction in units of the unit of count
a real noble atom               not reached: several charges on one centre need the rule for the kinds of matter
```

## What is put in

- OB1-U6's force, PR2-T2's energy–momentum, the whole-area rule of OR1/WQ1 for both pairs.
- A fixed centre; the attraction strength K; for the numbers, the constants of FS1.

## What is not shown

- These are the known circles, levels and spectrum of one charge around a centre (general knowledge). The line's
  part: each quantity is one of its own — speed, memory, count factor, left-over turn — and α is identified in
  each.
- The 0.18% in the second-shell split is the sum of two known effects not in this stage: the left-over turn of
  the moment (FS1-A3) and the finite mass of the centre.
- The measured split, binding and wavelengths are recalled, not re-read at source.

## Claim boundary

```text
IN-OUT / ROUND = √(1 − K²/L²) FOR EVERY BOUND HISTORY                                 PROVED (symbolic)
CIRCLES: SPEED α/l, MEMORY α²/l², E = g √(1 − α²/l²)                                  PROVED (exact)
THREE LENGTHS IN RATIO α ; SIZE = 2a/α                                               PROVED
BOTH AREAS WHOLE ⇒ THE LEVEL FORMULA ; SECOND-SHELL SPLIT α²/16                       PROVED (quadrature and symbolic)
LEFT-OVER TURN π α²/l² ; SIX TIMES AS MUCH IN THE FIELD OF A MASS                     PROVED (first order)
MEASURED SPLIT WITHIN 0.2 % ; WAVELENGTH RATIOS = COUNTS TO 10⁻⁵                      NUMBERS (measured values recalled)
THE VALUE OF α ; ATOMS WITH SEVERAL CHARGES                                          NOT DERIVED
```

## Reproduce

```text
python br1_bound_record_of_a_charge.py
python -m unittest test_br1
```

## Later note (SH1, AC1, 9 October)

SH1 counts what a shell holds: 2n² least-cost readings, limited by independence. AC1 puts the same first-circle
question to a mass: α < 1 here, g·r_s/2κ ≤ 1/(2√3) there, and no relation between the two. Nothing above is changed.
