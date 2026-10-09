# GC1 — The gap of the core as a count: a floor for the second rate, the lowest rate from both sides, and the walk in the number of turns

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact integers and rationals; one trial vector
is found in floats and then used as an exact rational vector.

CR1 ended with "a lower bound for the second rate: not obtained". CM1 added readings that turn under the
directions (7.5787 and 10.9738 from above), a written argument that the rates are discrete, and the same open
line: no number for the gap from below. This stage gives that number for three and for four turns, pins the
lowest rate from both sides, and finds one law that holds for every number of turns.

Sources read before building: **CR1** (R3 the ladder; R4 the comparison by the layer), **CM1** (§2 turning
readings, §4 the mean of h² on the ladder, §5 the count with a hidden floor), **TC1** (C4 the layer, C6 the
scale), the transfer map (§4, counts under a line); tools **DS1** (the count of directions past a line),
**AG1-A7** (a gap as a count); RKF **theorum/28** (a certificate from below needs its own floor).

## The operator and the comparison

C is a free matrix with columns c₁ … c_d (the turns), X the sum over pairs of |c_i × c_j|², and
h = −½·Δ + X. Readings unchanged when all cuts are turned together (the row condition) are the physical ones.
CR1-R4, step (i): give half of every turn's kinetic part to its d − 1 pairs; each pair, as an oscillator across
the other turn, is at least its ground rate. What is left is

```text
h  ≥  k₁ + … + k_d ,        k = −¼·Δ + β·|c| ,      β² = (d − 1)/2 ,
```

one copy of k for each turn, with no term between turns. CR1 used only the lowest level of this comparison.
This stage counts with it.

## Results

**G1 — one turn.** With |c| = ℓ·s, ℓ³ = 1/(4β), the operator k on readings u/s is

```text
unit × ( −u″ + [n(n + 1)/s²]·u + s·u ) ,        unit = ((d − 1)/8)^(1/3) ,
```

n the turning number of the reading on the sphere. One step proves floors: if v is positive on an interval and
−v″ + W·v = λ·v there, then for u = v·w vanishing at the ends ∫(u′² + W·u²) = ∫v²·w′² + λ·∫u² ≥ λ·∫u².

For n = 0 let v be the solution of v″ = (s − λ)·v with v(0) = 0, v′(0) = 1: an exact power series with rational
coefficients (114 and 133 terms kept; what is left adds less than 10⁻²⁰). On every step of a grid of 1/256 the
equation itself bounds v and v′, and past s = λ a solution keeps its sign once v and v′ agree in sign. Certified:

```text
λ = 2.3381 :   v > 0 on the whole half line                                first level  ≥ 2.3381
λ = 4.0879 :   v > 0, one node (between 1.746 and 1.75), v < 0 after it      second level ≥ 4.0879
```

(with one node, the form is at least λ on all readings that vanish at the node — one condition.) At 2.3382 and
at 4.0880 the same certificates are refused.

For readings with zero average over every sphere the term n(n + 1)/s² is at least 2/s². The reading
s²·Exp(−k·s^(3/2)) has local rate (27k/4)·s^(−1/2) + (1 − 9k²/4)·s, whose least value by the three-term mean is
3·(A²B/4)^(1/3); at k² = 2/9 this is

```text
level  ≥  (9/4)·3^(1/3)  =  3.2451 .
```

**G2 — d turns: at most one reading under a line.** Split every turn into its average over its own sphere and
the rest. These splittings commute with the comparison. A reading unchanged by turning all cuts together has no
part in which exactly one turn is outside its sphere average: such a part is unchanged by turning every other
turn alone and by turning all together, hence by turning that one turn alone. So a row-unchanged reading is

```text
all turns averaged                    levels of d radial operators :  one reading at the bottom, the rest
                                      at least (d − 1)·e₀ + e₁
two or more turns outside             at least (d − 2)·e₀ + 2·p₀
```

and at most **one** reading of the comparison, hence of h, lies under

```text
z = unit × min( (d − 1)·e₀ + e₁ , (d − 2)·e₀ + 2·p₀ ) :        5.5210 (three turns) ,   8.0060 (four turns) .
```

Control: without the row condition a single turn may leave its sphere average, and the same count stops at
4.9900 and 7.3981 — under the lowest rate itself. The floor of the second rate exists because the readings are
unchanged by the row turning.

**G3 — the gap.** The lowest rate is under the mean of h in any reading; in the reading of G4 that mean is under
5.1868 and 8.0029. The second rate is at least z. So, in the whole row-unchanged space,

```text
second rate − lowest rate   ≥   0.3342 (three turns) ,     0.0031 (four turns) .
```

**G4 — the lowest rate from both sides.** With one reading under z and all others at or above it,
(h − E)(h − z) has no negative mean, E the lowest rate. For any reading ψ with mean η of h below z this says

```text
E  ≥  η − ( mean of h² − η² ) / ( z − η ) .
```

ψ: one exact rational vector on CR1's ladder of degree 10 (67 readings); h·ψ lies in the ladder of degree 12,
so both means are exact (CM1's source Gram is this mean of h² with the retained part taken out). The spread
mean of h² − η² is 7.9·10⁻⁵ and 2.7·10⁻⁵:

```text
                 lowest rate               before (CR1)
three turns      [ 5.1865 , 5.1868 ]       [ 4.14 , 5.1868 ]
four turns       [ 7.9942 , 8.0029 ]       [ 6.32 , 8.0030 ]
```

**G5 — the window of the gap.** A reading unchanged by the directions and one that changes sign when two
directions are exchanged are orthogonal for 1 and for h, so the second rate is under the larger of their two
means: 7.5787 and 10.9738 (CM1).

```text
                 second rate               gap
three turns      [ 5.5210 , 7.5787 ]       [ 0.3342 , 2.3922 ]
four turns       [ 8.0060 , 10.9738 ]      [ 0.0031 , 2.9796 ]
```

**G6 — the walk in the number of turns.** The core of n turns is the sum of its n sub-cores of n − 1 turns,
each in the form a·T + b·X with a = 1/(n − 1), b = 1/(n − 2) (a turn is in n − 1 of them, a pair in n − 2). A
dilation makes a·T + b·X equal to (a²b)^(1/3) times the core of n − 1 turns, and the lowest rate of a sum is at
least the sum of the lowest rates. Write E(d) for the lowest rate over all readings (it is the lowest rate over
row-unchanged readings: passing to |ψ| does not raise the form, and the kinetic part is convex in ψ², so
averaging ψ² over the turnings does not raise it either). Then

```text
E(n)  ≥  n·(n − 1)^(−2/3)·(n − 2)^(−1/3)·E(n − 1) ,     which is :     ρ_d = E(d) / ( d·(d − 1)^(1/3) )  never falls as d grows .
```

The Gaussian reading gives ρ_d³ ≤ 729/256 for every d. So ρ_d rises to a limit not above 1.4175, and with G4:

```text
ρ₃ in [ 1.3721 , 1.3723 ] ,     ρ₄ in [ 1.3857 , 1.3873 ] ,     every d ≥ 4 :   1.3857 ≤ ρ_d ≤ 1.4175 .
five turns :  lowest rate in [ 10.9985 , 11.25 ] .          two turns :  lowest rate ≤ 2.7445 .
```

The two certified windows stand in the order the walk demands (ρ₃ < ρ₄), and one step from three turns gives
7.9160 for four, under the 7.9942 of G4. The same form d·(d − 1)^(1/3) is in CR1's bound from below, in the
comparison and in the Gaussian value: it is the form of the walk.

## The numbers

```text
                 lowest rate              second rate             gap
three turns      5.1865 … 5.1868          5.5210 … 7.5787         0.3342 … 2.3922
four turns       7.9942 … 8.0029          8.0060 … 10.9738        0.0031 … 2.9796
```

In the units of TC1 these are multiplied by g^(2/3), per size of the torus. With the rate operator h (time kept
apart) the turns are the directions of space, so three turns is the case of three directions of space.

## Inputs taken from outside

```text
the kinetic part of one turn in its radius and its sphere ; on a sphere a reading with zero average has
    mean of |gradient|² at least 2 × mean of its square
a self-adjoint rate operator splits its readings by level (used for the radial part of one turn and for h on
    the row space) ; a form that is ≥ z on all readings obeying one condition has at most one level under z
|∇|ψ|| ≤ |∇ψ| ; |∇ρ|²/ρ is convex in ρ                                    (G6 only)
CM1's certified turning readings ; CR1's ladder
```

The discreteness of all rates (CM1 §3) is not needed for G2 – G4: one level under z is all that is used.

## What the certificate checks

```text
G1    the kept polynomials solve the equation up to their last two degrees ; tails under 10⁻²⁰ ; the two sign
      certificates and their two refusals ; p₀³ ≤ 2187/64 ; e₀ ≤ e₁, e₀ ≤ p₀
G2    unit³ ≤ (d − 1)/8 ; the line ; the control without the row condition
G3    η under CR1's bound and under the line ; the floor of the gap
G4    ⟨ψ, hψ⟩ two ways (exact) ; the spread ; the floor, above CR1-R4 and above d copies of one turn
G5    CM1's turning value above the line ; the window
G6    the Gaussian mean of X and its value (d = 3, 4, 5) ; the step identity (n = 3 … 12) ; the order of the
      two windows ; the step from three turns under the floor of G4
```

## What this changes in the line

```text
CR1    second rate from below                 was: not obtained              now: ≥ 5.5210, ≥ 8.0060, in the whole
                                                                             row-unchanged space
CR1    lowest rate                            was: [4.14, 5.1868], [6.32, 8.0030]     now: [5.1865, 5.1868], [7.9942, 8.0029]
CM1    number for the gap from below          was: open                      now: 0.3342 and 0.0031
CM1    "no actual hidden floor"               the comparison is one: for the cut of one reading (the product of
                                              the lowest one-turn readings) the floor of everything else is z
```

The earlier notes are not rewritten; a later note in each points here.

## What is not shown

- **No mass gap.** The gap of the constant modes is g^(2/3) times these numbers per size of the torus: at a
  fixed coupling it goes to zero as the torus grows. The fabric (the modes that are not constant) is untouched.
- **The floor is far from the ceiling.** For three turns 0.334 against 2.392; by the older values the gap is
  2.391, so the ceiling is the sharp end. The comparison has no term between turns, and its second level is
  5.5211 whatever is done with it. A larger floor needs a cut that holds the comparison's low readings
  (products of one-turn readings): the comparison then gives that cut a true floor — the first level left out —
  and CM1's count with the source Gram of that cut can be run with it. Not built.
- **Four turns is thin (0.003), and five or more turns get nothing** from this comparison: by the older values
  of the two levels its line is under 1.344·d·(d − 1)^(1/3), below the lowest rate (ρ_d ≥ 1.3857).
- Whether ρ_d reaches the Gaussian value 1.4174 as d grows.
- The inputs listed above are taken, not derived.
- In older terms: the two levels of G1 are minus the first two zeros of the Airy function, 2.338107 and
  4.087949 (the tests bracket them); the step of G1 is the local-rate bound for a positive reading, G4 is
  Temple's inequality, and G6 is a decomposition of the kind used for many-body ground rates. With the usual
  normalization (2^(1/3)) the lowest rate 4.1167 is inside the window of G4. All general knowledge, not
  re-read; I did not search whether the floor of G3 or the monotone ρ_d are in print. The line's part: the
  count under the row condition, the exact sign certificates, the two-sided lowest rate on CR1's ladder, and
  ρ_d with its windows.

## Claim boundary

```text
FIRST TWO LEVELS OF −u″ + s·u ≥ 2.3381, 4.0879 ; ZERO-AVERAGE LEVEL ≥ 3.2451           PROVED (exact series ; local rate)
AT MOST ONE ROW-UNCHANGED READING OF h UNDER 5.5210 (THREE TURNS), 8.0060 (FOUR)       PROVED, given the inputs
GAP ≥ 0.3342 (THREE), ≥ 0.0031 (FOUR)                                                  PROVED, given the inputs
LOWEST RATE IN [5.1865, 5.1868] (THREE), [7.9942, 8.0029] (FOUR)                       PROVED, given the inputs
GAP ≤ 2.3922 (THREE), ≤ 2.9796 (FOUR)                                                  PROVED (CM1's readings)
ρ_d NEVER FALLS ; ρ_d³ ≤ 729/256 ; 1.3857 ≤ ρ_d ≤ 1.4175 FOR d ≥ 4                     PROVED, given the inputs
A FLOOR NEAR THE CEILING ; FIVE OR MORE TURNS ; THE LIMIT OF ρ_d                       NOT OBTAINED
A MASS GAP ; THE FABRIC                                                               OPEN
```

## Reproduce

```text
python gc1_core_gap_count.py          (about three seconds)
python -m unittest test_gc1
```

## Later note (SG1 and OL1, 9 October)

[SG1](../sg1/SG1_GAUGE_SINGLET_GAP.md) was written at the same time, independently, from the same comparison
and without sight of this note. It has the same floors — 5.5210 and 8.0060 for the second rate, 0.3342 and
0.0030 for the gap — by a different solve (rational intervals with a cut at 12 and a count with its endpoint
term), and a certified reading of the comparison without the row condition under the lowest rate, which makes
the control of G2 a statement about the comparison itself. G4 (the lowest rate from both sides) and G6 (the
walk in the number of turns) are only here. [OL1](../ol1/OL1_ONE_SITE_LATTICE.md) carries G1 and G2 to compact
links: on the lattice of one site the same two sign certificates give a gap at every coupling past 10⁷.
Nothing above is changed.
