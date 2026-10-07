# NC1 — The framework's own curvature of the gravity field: a ratio law, and a law that nature refuses

Monty Dabas. 8 October 2026. Python 3.12. Symbolic algebra (sympy —
outside the standard-library suite).

CV1 and TP1 used the curvature of a connection built from the frame and
the invariant form — the geometric curvature. The framework has its own
definition, and it is a different object. This stage uses that one.

Sources read before building: Publications
`papers/native-thermodynamic-curvature` **NT-1** (response element
H = ((a+c)/2)1 + ((a−c)/2)K + bRK; the six responses;
Γ_cΓ_m = (C_V/C_P)(K_S/K_T) = 1), **NT-2** (curvature = transport-square
defect), **NT-3** (X_i = H⁻¹δ_iH; F^{qX} = q(q−1)[X_i,X_j]; q = 1 flat;
strain choice q = ½, F^H = −¼[X_i,X_j]); **QC5** (the same q(q−1) as a
share law); **GB1** (the gravity field is a self-dagger unit block);
**FR1** (three cuts); **MO1** (radial law); **CV1**, **TP1**.

## 1. The native curvature of the field

The field of the thesis is the self-dagger unit block exp(η(r) n·C),
tanh η = β, m = β². Its response element is G = exp(ψ n·C), ψ = 2η.

**N1.** NT-3 holds on this field: F^{qX} = q(q−1)[X_i,X_j]; q = 1 is
flat; q = ½ gives F = −¼[X_i,X_j].

**N2.** The native curvature has two frame-independent sizes:

```text
radial–angular     −(ψ′ sinh ψ)² / 4r²            (twice)
angular–angular    −(sinh ψ)⁴ / 4r⁴
```

and one pure number between them,

```text
ρ = r ψ′ / sinh ψ ,
```

the ratio of radial to angular native curvature. It has no unit and does
not depend on the frame — a ratio of the same kind as C_P/C_V.

In the response dictionary of NT-1 the same field read across a cut has
C_V/C_P = 1 − tanh²η = 1 − m: the clock factor squared is the ratio
C_V/C_P, and the memory is (C_P − C_V)/C_P.

## 2. The ratio law

**N3.**

```text
ρ = −½ at every radius      ⇔      tanh²η = r_s / r   exactly.
```

The field of CV1 and MO1 — with Mercury's 42.98″ and the Sun's 1.751″ —
is the field whose radial native curvature is half its angular native
curvature, at every scale. (ρ = −1 would give memory falling as 1/r².)

## 3. The variational law

The other natural law is the one of the Yang–Mills line: make the total
invariant curvature-square stationary.

**N4.** With w = cosh ψ it reads

```text
w″ = w ( w² − 1 ) / r² .
```

To first order, w − 1 = A/r + B r²: the 1/r field and the
constant-curvature term, the same two solutions as CV1-T4. Beyond first
order, w = 1 + A/r + ¾A²/r² + …, and the field of N3 does not solve it.

**N5.** With MO1's radial law, a memory m = r_s/r + κ r_s²/r² advances
an orbit by (3 + 2κ)π r_s/p per turn.

```text
law                                   κ        Mercury, per century
ratio law (N3)                        0        42.98″      observed 42.98″
variational law, ψ = 2η               +½       57.3″       refused by observation
variational law, ψ = η                −⅜       32.2″       refused by observation
```

## 4. Reading

- The native curvature is not the geometric one. It is never zero for a
  field, it carries a pure ratio, and the observed field is a statement
  about that ratio.
- For the first time in this line a law built only from framework
  objects gives a prediction that differs from general relativity — at
  second order in r_s/r — and the difference is large enough to test.
  Nature refuses it: Mercury's advance is 42.98″, not 57.3″ or 32.2″.
- So the native law of the gravity field is not the Yang–Mills-type
  stationarity of its native curvature. It is the ratio law, or
  something equivalent to it.

## 5. What is not shown

```text
WHY ρ = −½                         not derived. It restates m ∝ 1/r in native terms; ½ is also the strain share q of NT-3,
                                   and the exponent of 1/r in three cuts — neither link is proved here.
OTHER NATIVE LAWS                  only these two were tried.
THE MOTION                         N5 uses MO1's radial law unchanged; a different field law could come with a different
                                   law of motion.
TIME-DEPENDENT OR NON-ROUND FIELDS not built.
```

## 6. Claim boundary

```text
NT-3 ON THE FIELD; TWO INVARIANT SIZES; THE RATIO ρ                       PROVED (symbolic)
ρ = −½ ⇔ tanh²η = r_s/r                                                   PROVED
VARIATIONAL LAW w″ = w(w²−1)/r²; FIRST ORDER 1/r AND r²                    PROVED
SECOND ORDER ¾A²; κ = +½ OR −⅜; ADVANCE (3+2κ)π r_s/p                      PROVED given MO1's radial law
VARIATIONAL LAW AS THE LAW OF GRAVITY                                      REFUSED by Mercury's orbit
ρ = −½ FROM A PRIOR STATEMENT OF THE FRAMEWORK                             NOT DERIVED
```

## 7. Reproduce

```text
pip install sympy
python nc1_native_curvature_of_the_field.py
python -m unittest test_nc1
```
