# CZ1 — The centre record of a closed surface: the weight of a twist is an exponential of the turn coupling

Monty Dabas. 9 October 2026. Python 3.12. Exact quaternion algebra, an exact gluing engine, exact enumeration over
centre marks; mpmath for the sums.

CT1 left the centre record in more dimensions open: on every face the two-valued coupling is κ·c_p, set by the
cosets. This stage takes the smallest closed surface, a cube, and asks whether the centre record closes there.

Sources read before building: **theorum/41** (cut-graded generator; even and odd channels), **CG1**, **CT1**,
**DU1** (dual coupling; self-dual = seen = lost), **YM-3/4** (exact pairing of contents: joining along a link
costs 1/d), **YM-F1** (the turn block is the unit sphere), **YM-37** (rational moments of the sphere), **YM-39**
(centre grading), **YM95/96** (Laplacian of a face), **T24-6.1** (F = R − D).

## Results

**Z1 — joining faces.** From the second moments of the sphere, for the first content:
⟨c(AU)·c(U⁻¹B)⟩ = c(AB)/4 and ⟨c(AUBU⁻¹)⟩ = c(A)·c(B). In characters: joining two faces along a link costs 1/d,
and so does a link met twice by one face (YM-3/4 for every content).

**Z2 — a closed surface carries d^χ.** Gluing all faces of a closed surface (cube, boxes of 10 and 16 faces, tori
of 4 and 9 faces — exact engine): the surface with content j on every face weighs

```text
d_j^χ · (f_j / f₀)^F ,        χ the Euler number:  2 for a sphere, 0 for a torus.
```

**Z3 — at fixed cosets the record closes exactly.** Summing the centre marks of the twelve links of a cube:

```text
Σ over marks  Π_p e^(k_p·s_p)  =  2¹² ( Π_p cosh k_p + Π_p sinh k_p ) ,        k_p = κ·c_p .
```

Even channel times even channel plus odd times odd (theorum/41's two channels). The ratio is Π tanh k_p, so the
dual couplings add over the faces: K* = Σ k_p*. DU1's law for a chain, now for a surface. The seam curvature of
the weight is zero (the faces' weights commute).

**Z4 — after the cosets.** For the cube Z = Σ_j d_j² f_j⁶. Its centre-odd over centre-even part is 4r⁶ at the
strong end, not r⁶: the record of faces does not simply compose — the factor is d^χ = 4.

The flip of the cube, F = (even − odd)/(even + odd), is the weight of a cube with one face reversed: a twist.

```text
one face       F = (1 − r)/(1 + r)  →  3/(4κ)                              a power
a cube         F  ≈  4.27 · κ · e^(−(6 − 3√3)·κ)                           an exponential of the turn coupling
F faces        exponent  F·(1 − cos(π/F)) :    2 ,  4 − 2√2 ,  6 − 3√3 ,  … ,  π²/2F
```

The exponent is the cost of turning every face by one F-th of a full turn (their product is then −1); the sums
give it to 3 × 10⁻⁴ for 2, 4, 6, 10, 16 faces. An open face can undo a twist smoothly, a closed surface cannot.

So at the weak end a closed surface of the centre record sits on the diagonal — even = odd, seen = lost — up to
an exponentially small flip. The cube is self-dual at κ = 4.447.

**Z5 — the other form.** With the generator ΣΔ + θΣW and the cut "reverse one link", theorum/41's seam curvature
[G_e, G_o] acts on the free vacuum as −3θ × (the faces through the link): each of a face's four links gives −¾.

## What this gives

```text
does the centre record close on a cube?       at fixed cosets: yes, dual couplings add
                                              after the cosets: up to d^χ and the higher contents
a native exponential of the turn coupling     the weight of a twist of a closed surface: e^(−(6 − 3√3)κ) for a cube
size of number                                a twist weight of 7.7 × 10⁻²⁰ is κ = 61.7 — not 10¹⁹
```

In CT1 the proton's count needed κ = 2 × 10¹⁹ on an open chain. On a closed surface the same smallness is a turn
coupling near 60. The difference is the closing of the surface.

## What is not shown

- A twist weight is not a mass. No gap of the colour record is derived; the sum over all closed surfaces in four
  dimensions is not done.
- Z = Σ d_j² f_j^F uses YM-3/4's pairing for every content; the first content's two rules are verified here.
- The exponent F(1 − cos(π/F)) is the classical value for equal turns about one axis; that no other arrangement
  costs less is not proved. The sums agree with it numerically; the prefactor 4.27·κ is read off, not derived.
- Reversed cubes of the colour record and their exponentially small weight are known (general knowledge, not
  re-read at source). The line's part: the flip F = R − D of the closed centre record, the factor d^χ, and the
  exponent as one F-th of a turn per face.

## Correction to an earlier statement

After CG1 it was said that the seam curvature [G_e, G_o] would be the obstruction to one doubling map. For the
weight it is not: that curvature is zero and the record closes at fixed cosets (Z3). The obstruction after the
cosets is d^χ and the higher contents (Z4). The curvature is non-zero only in the generator form (Z5).

## Claim boundary

```text
JOINING RULES FROM THE MOMENTS OF THE SPHERE ; CLOSED SURFACE CARRIES d^χ                            PROVED (exact engine, five surfaces)
CENTRE SUM OF A CUBE AT FIXED COSETS = 2¹²(Π cosh + Π sinh) ; DUAL COUPLINGS ADD                     PROVED (enumeration)
CUBE: ODD/EVEN = 4r⁶ AT THE STRONG END ; FACE TWIST → 3/(4κ)                                         PROVED (numerically)
TWIST WEIGHT OF F FACES ~ e^(−F(1 − cos(π/F))κ) ; CUBE 6 − 3√3                                       NUMERICAL (3 × 10⁻⁴) ; classical value
A MASS COUNT                                                                                         NOT DERIVED
```

## Reproduce

```text
python cz1_centre_record_of_a_closed_surface.py
python -m unittest test_cz1
```
