# DG1 — The diagonal: where observed equals lost

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic.

The owner's remark: ½ is the observer on the diagonal — sin 45° =
cos 45° — and a link should be there.

Sources read before building: **GR1** (N² + m = 1; the horizon is all
memory), **MO1** (circular histories), **NC1** (the ratio ρ), **NT-3**
and **QC5** (share law q(1−q)), **R38 / HB1** (the eight-mark clock,
1/√2), **MC1**.

With x = r_s/r: observed N² = 1 − x, lost m = x.

**D1.** Observed equals lost at r = 2 r_s. There N = 1/√2 — the value
of the eight-mark clock of R38.

**D2.** At that radius, and only there,

```text
the speed of a circular history  =  the speed of the fall frame        ( both c/√2 )
a circular history has exactly the energy of rest far away              ( no binding )
```

Outside it circular histories are bound; inside it they are not. The
diagonal is the boundary between held and free.

**D3.** The share law q(1−q) is largest at q = ½, where it is ¼ — the
factor in the native curvature F = −¼[X_i, X_j].

**D4.** ρ = ½ · d ln m / d ln r. The ½ in ρ = −½ is the step from the
amplitude (tanh η, the odd part) to its square (the memory); the −1 is
the fall-off of memory in a frame of three cuts (MC1).

## Reading

The same ½ appears in four places for one reason each time: memory is
the square of the odd part, and the two shares are equal on the
diagonal. In the gravity field the diagonal is a real place, 2 r_s,
with a real meaning. The radius is known in orbital mechanics as the
innermost bound circular orbit; the link to the equal share and to the
eight-mark clock is the framework's.

## Claim boundary

```text
EQUAL SHARE AT r = 2 r_s, N = 1/√2                                  PROVED
ORBIT SPEED = FALL SPEED, AND ZERO BINDING, ONLY THERE              PROVED from MO1's circular histories
SHARE LAW MAXIMUM ¼ AT ½                                            PROVED
ρ = ½ d ln m / d ln r                                               PROVED (NC1)
A REASON WHY THE FIELD IS ρ = −½ FROM THE DIAGONAL ALONE            NOT FOUND
```

## Reproduce

```text
python dg1_the_diagonal.py
python -m unittest test_dg1_exact
```

## Later note (DO1, 8 October)

"A reason why the field is ρ = −½ from the diagonal alone — not found": DO1 finds it. For the observer whose cuts
are equally inclined to the fall, the stretch (ρ, 1, 1) has seen = lost exactly for ρ = −½ (in d cuts −(d − 2)/2).
Nothing above is changed.

## Later note (AC1, 9 October)

The value N = 1/√2 of the held reading at 2r_s has a partner: the free circle at the last stable circle, 3r_s,
has count factor² = ½ exactly. Nothing above is changed.
