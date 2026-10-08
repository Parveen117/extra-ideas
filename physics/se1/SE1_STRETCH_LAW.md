# SE1 — The frame law on flat slices is a law for the stretch alone: the plane determinants sum to zero

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy), with the functions of CV1 and TP1.

OR1 asks whether the line's rules are one rule. For the field of gravity that needs the law of TP1 in a form that
can be compared with a record memory. This stage writes it for every fall velocity on flat slices.

Sources read before building: **TP1** (Q = ¼I₁ + ½I₂ − I₃, selected by equivalence of local frames; = −R + boundary),
**CV1** (connection and curvature from the order defects), **SW1-T4** (frame e₀ = ∂_t + u·∇ on flat slices: turn
½ curl u, stretch = symmetric part of ∂u; T3: remainder at second order in the swirl), **MA1**, **MC1** (exponent
d − 2), **NC1** (the ratio ρ; N3: ρ = −½ ⇔ m = r_s/r; "why ρ = −½" not derived), **EG1** (ratio law with a source:
ρ + ½ = k·u·r²/2m), **EMK-1** (determinant of a block in a plane: seen product minus lost square), **CD1**.

## Results

K is the stretch block, K_ij = ½(∂_i u_j + ∂_j u_i); λ_i its rates along its own axes; e₂(K) = Σ_{i<j} (K_ii K_jj − K_ij²)
the sum of its determinants in the planes of the cuts.

**S1 — the law sees the stretch only.** For every u,

```text
Q  =  Σ K_ij² − ( Σ K_ii )²  =  −2·e₂(K) .
```

The turn of the frame (curl u) does not enter. The same holds with TP1's formulas written for two and for four
cuts.

**S2 — the two parts of the law.** With CV1's curvature:

```text
energy component        G₀₀ = e₂(K)
(time, cut) components  −½ · curl curl u
```

So in empty space: e₂(K) = 0, and the turn field ½ curl u is itself free of curl (a gradient — SW2-H).

**S3 — a ratio law.** e₂ = 0 says

```text
Σ λ_i²  =  ( Σ λ_i )²        ⇔        variance of the rates  =  ( d − 1 ) × ( their mean )² .
```

The stretch rates across the cuts add like perpendicular legs: the square of the sum is the sum of the squares.
The spread between the cuts, relative to their common part, is the number of other cuts: 2 in three dimensions.
It is a ratio between a variance and a mean, not a least variance.

**S4 — radial fall.** For u = −β r̂ in d cuts the rates are (ρ, 1, …, 1) × (−β/r) with ρ = r·β′/β, and

```text
e₂ = (β/r)² · [ (d − 1)·ρ + (d − 1)(d − 2)/2 ]        ⇒        ρ = −(d − 2)/2 ,     m = β² ∝ r^−(d−2) .
```

- ρ here is NC1's ratio (r·ψ′/sinh ψ with ψ = 2η, tanh η = β): the same number.
- In three cuts ρ = −½: rates (−½, 1, 1), sum 3/2, squares 9/4. NC1's open item is answered from TP1's law:
  −½ is −(number of cuts − 2)/2, and 2ρ is MC1's exponent.
- With a source, G₀₀ = e₂ = k·u is EG1's ratio law ρ + ½ = k·u·r²/2m.
- In two cuts e₂ = (β/r)²·ρ: the law leaves β constant. There is no falling memory in two dimensions; MC1's log r
  is not allowed by the frame law.

**S5 — the swirl.** For SW1's u = −β r̂ + W(r) ẑ×r: curl curl u = 0 is SW1-T1's (r⁴W′)′ = 0; and with W = a/r³

```text
e₂ = − ( shear of the swirl )² = − 9 a² sin²θ / 4r⁶ .
```

SW1's refusal in one line: the swirl adds a lost square to one plane's determinant, and on flat slices nothing
can balance it.

## What this settles

```text
NC1:   "why ρ = −½"                  from TP1's law: e₂ of the stretch vanishes; ρ = −(d − 2)/2
EG1:   ratio law with a source       the energy component of the same statement
SW1:   second-order refusal          e₂ of a swirl is minus a square
the law's kind                       a ratio (variance : mean² = d − 1), not a least variance      → OR1
```

## What is put in

- TP1's law (its premises: quadratic, coordinate-free, equivalence of local frames) and CV1's formulas.
- One time and flat slices (MA1: exact for a centre at rest; SW2: not beyond first order for a turning centre).

## What is not shown

- S1–S2 are the known form of the field equations for this family of frames (general knowledge; no source
  re-read). The line's part: the reading as plane determinants and as a ratio between variance and mean, the
  identification of NC1's ρ, and S5.
- Flat slices only. On SW2's exact frame the stretch is not the whole story.
- Why the ratio is d − 1 rests on TP1's equivalence premise; no further reason is given.

## Claim boundary

```text
Q = Σ K_ij² − (tr K)² = −2 e₂(K) FOR ANY FALL VELOCITY ; TURN ABSENT                PROVED (symbolic; 2, 3, 4 cuts)
G₀₀ = e₂(K) ; (TIME, CUT) COMPONENTS = −½ curl curl u                               PROVED (symbolic, three cuts)
e₂ = 0 ⇔ VARIANCE = (d − 1) MEAN²                                                   PROVED
RADIAL: ρ = −(d − 2)/2 ; = NC1's ρ ; −½ IN THREE CUTS ; EG1's LAW WITH A SOURCE      PROVED
TWO CUTS: β CONSTANT                                                                PROVED
SWIRL: e₂ = −9a² sin²θ/4r⁶                                                          PROVED
"WHY ρ = −½" (NC1)                                                                  ANSWERED from TP1's law
A REASON FOR TP1's PREMISES                                                         NOT GIVEN
```

## Reproduce

```text
python se1_stretch_law.py
python -m unittest test_se1
```
