# CV1 — The curvature of the frame field: what no frame change removes

Monty Dabas. 7 October 2026. Python 3.12. Symbolic algebra (sympy —
outside the standard-library suite); two independent routes agree.

The owner's question: the thesis is about curvature, and curvature is
what every frame agrees on — so put the curvature of the frame field in,
and see what law comes out.

Sources read before building: **EMK-C1** (curvature as order defect,
with the frame term −c_μν^ρ A_ρ), **DB1/YM94** (derivation algebras),
**RMG4/RMG5** (curvature of response processes; intrinsic; threshold in
dimension), **IN1** (the invariant form), **GR2-G1** (static memory is
removed at a point by one frame change), **MC1**, **MO1**, **GB1**. The
repository `information-curvature-duality-provisional` was read: its
engine is a detector of symmetry residues across eight channels, not a
law for a frame field; it is not used here.

## 1. The construction — framework form only

The thesis's frame field (MO1, GB1), c = 1:

```text
e₀ = ∂_t − β(r) ∂_r       the locally pure frame's time direction
e₁ = ∂_r ,  e₂ = (1/r) ∂_θ ,  e₃ = (1/(r sin θ)) ∂_φ
invariant form  η = diag(+,−,−,−)                                   (IN1)
```

1. **Order defect of the directions**: [e_a, e_b] = c_ab^c e_c.
2. **The connection** that keeps η and has no torsion — fixed by the
   c's alone.
3. **Curvature = order defect of the connected directions**, in the form
   of EMK-C1:
   R(e_c, e_d) = e_c ω_d − e_d ω_c + [ω_c, ω_d] − c_cd^e ω_e.

The order defect of this frame is the gradient of the fall:

```text
c₀₁¹ = β′ ,     c₀₂² = c₀₃³ = β / r ,     and the flat-space angular terms.
```

## 2. Results

With m = β² (the memory):

**T1.** The order-defect route and the coordinate route (from the
interval of MO1) give the same curvature. Every component depends on β
only through m, m′, m″.

**T2 — the law.** The contracted curvature is diagonal:

```text
( −(rm)″/2r ,  +(rm)″/2r ,  (rm)′/r² ,  (rm)′/r² ) .
```

It vanishes exactly when (r m)′ = 0, that is

```text
m = r_s / r .
```

**T3 — what remains.** For m = r_s/r the curvature itself is not zero:

```text
on the planes (time, radial), (time, θ), (time, φ):     ( r_s/r³ ,  −r_s/2r³ ,  −r_s/2r³ )
sum zero ;   curvature square  12 r_s² / r⁶ .
```

GR2-G1 removes the memory at a point by one frame change; this is what
that change cannot remove. It is the same in every frame: the invariant
content of the field. A uniform fall (β constant) is also curved — a
constant fall speed toward a centre is not a change of frame.

**T4.** (r m)′ = Λ r² gives contracted curvature = −Λ η: the
constant-curvature case, m = r_s/r + Λ r²/3.

**T5 — two native routes meet.** The time–time contracted curvature is
exactly

```text
R₀₀ = −½ · (1/r²)( r² m′ )′ ,
```

minus half the variance operator of MC1. So MC1's law — least total
variance between neighbouring shells — is the vanishing of the
time–time contracted curvature. Least variance alone allows
m = A/r + B; the remaining components, (r m)′/r², remove the constant.
The curvature law is the sharper one.

## 3. What this is and is not

```text
is       the thesis's field law restated as one statement about the frame field:
         the contracted order defect of the connected frame directions is zero
is       an exact link between MC1 (variance) and curvature
is       the answer to "what is invariant": the uncontracted curvature, r_s/r³ × (1, −½, −½)
is not   new: this is the known vacuum law for this family of frames, and step 2 is the standard
         construction of a connection from a frame and a form
is not   derived: "contracted curvature = 0" is taken as the law; nothing here obtains it from a prior
         statement of the framework
is not   a result without relativity: the frame, the invariant form and the connection are its content,
         written in the framework's terms
```

## 4. Claim boundary

```text
ORDER-DEFECT CURVATURE OF THE FALL FRAME; TWO ROUTES AGREE                      PROVED (symbolic)
CONTRACTED CURVATURE ZERO ⇔ m = r_s/r                                           PROVED for this family
UNCONTRACTED CURVATURE r_s/r³ (1, −½, −½); SQUARE 12 r_s²/r⁶                     PROVED
R₀₀ = −½ × MC1's OPERATOR; CURVATURE LAW REMOVES THE CONSTANT                    PROVED
CONSTANT-CURVATURE CASE                                                          PROVED
WHY THE CONTRACTION AND NOT THE WHOLE CURVATURE                                  ANSWERED in TP1 (conditional on equivalence); originally: — information invariance in the
                                                                                 paper's sense (all loops closed) would force zero
                                                                                 curvature and no field at all
SOURCES INSIDE MATTER; FRAMES THAT ARE NOT STATIC OR NOT ROUND                   NOT BUILT
```

## 5. Reproduce

```text
pip install sympy
python cv1_frame_curvature.py
python -m unittest test_cv1
```
