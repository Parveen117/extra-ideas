# OB1 — One block, one law: matter, light, potential, source and force

Monty Dabas. 7 October 2026. Python 3.12. Exact arithmetic: polynomials
in (t, x, y, z) over the Gaussian rationals.

With a free hand, this stage asks the largest question the line has
earned: are the things physics keeps apart — matter, light, the
potential, the charge and current that source light, and the force light
exerts — different objects, or different parts of the one block the
framework already certifies?

Sources read before building: **EMK-1** T1–T3 (the block and its
relations), **FR1-T3/T4** (three cuts on the cut-complex carrier; the
conjugation anticommutes with ι), **IN1** (reading tensors; the
invariant), **PR1-T4** (propagation by cuts; the square law),
**PR4-T1** (the conserved current of a reading), **EM1** (F = E + ιB),
**SY1-S1** (scalar and traceless parts; O²), **QC1** (integer content).

## 1. The algebra has exactly four kinds of part

**U1.** C_iC_j = δ_ij + ι·ε_ijk·C_k, and ι = C₁C₂C₃. The cut-complex unit
is the product of the three cuts. A general block is

```text
( a + ι b )  +  ( r + ι s )·C          a: uncut      r: three cuts      s: three turns (ι C)      b: the volume ι
```

One operator, D = ∂_t + ΣC_i∂_i, and its reverse D̄ = ∂_t − ΣC_i∂_i, with
D·D̄ = ∂_t² − ∇²: the wave operator is the product of the law and its
reverse (PR1-T4 in three directions).

## 2. Each physical object is one kind of block

```text
potential      𝒜 = φ − A·C                    self-dagger: uncut + cuts        (reading-type)
light          F·C ,  F = E + ιB              traceless: cuts + turns
source         J = ρ − j·C                    self-dagger: uncut + cuts        (reading-type)
matter         ψ                              a column of the block
its reading    ψψ† = ½( n + r·C )             self-dagger: uncut + cuts        (reading-type)
```

## 3. One law gives all the equations

**U2 (the field is D̄ of the potential).** D̄𝒜 = L + F·C with
L = ∂_tφ + ∇·A the scalar part and E = −∇φ − ∂_tA, B = ∇×A.

**U3 (the four field equations are the four grades of one).**

```text
D( F·C ) = ( ∇·E + ι ∇·B )  +  ( ∂_tE − ∇×B  +  ι ( ∂_tB + ∇×E ) )·C .
```

For a field that comes from a reading-type potential, the two ι-parts
vanish identically.

**U4 (the source is a reading, and it is conserved).** J = D(F·C) is then
self-dagger — a real density and a real current, no ι-part — and
∂_tρ + ∇·j = 0 identically. There is no magnetic source because a source
is a reading, and readings have no ι-part.

**U5 (matter supplies exactly such a source).** For a column ψ,
∂_t n + ∇·r = ψ†(Dψ) + (Dψ)†ψ with (n; r) the reading of ψ. When Dψ = 0
the reading is conserved, self-dagger, and null — the kind of block U4
requires.

**U6 (the force is a frame change).** For a reading ρ = ½(n + r·C),

```text
F ρ + ρ F†   has scalar part  E·r   and vector part  n E + r × B ,
```

and it leaves n² − r·r unchanged. The field moves the reading and never
touches its unrecoverable part: E boosts, B turns. With SY1, F is the
traceless part O of a block and F·F = O².

A wave E = (0, p, 0), B = (0, 0, p) with p any polynomial of x − t solves
D(F·C) = 0 and has F·F = 0; the opposite handedness does not solve it.

## 4. The loop

```text
matter ψ  ──reading──▶  ρ = ψψ†   (conserved, self-dagger)
                           │  source
                           ▼
        D( F·C ) = q·ρ                light F from the reading of matter
                           │  frame change
                           ▼
        dρ ∝ q ( Fρ + ρF† )           the reading of matter moved by light
```

Everything in the loop is a part of the block; the only operator is D.
The number q is the content of the reading under the ι-turn — an integer
by QC1 — and it changes sign on the conjugate sheet, because the
conjugation anticommutes with ι (FR1-T4). So the coupling of light is odd
between the sheets; the clock factor of the gravity thesis is even.

## 5. What is new and what is not

Every equation in §3 is classical electrodynamics, and writing it in a
block algebra of this kind is known (general knowledge). Nothing is
predicted.

What this stage establishes is about the framework: the block certified
in EMK-1, read on the cut-complex carrier with the turn taken as a cut,
is by itself large enough to hold matter, light, potential, source and
force, with one law — and three facts that are usually separate
postulates come out as properties of the parts:

```text
no magnetic source          a source is a reading; readings have no ι-part
conservation of charge      the source of a traceless field has a conserved count
the force law               the field is a frame change of the reading
```

## 6. Certificate

U1 as matrix identities. U2–U4 on polynomial potentials of degree three
in four variables: the scalar part, the three field components against
the usual formulas, the four grades of D(F·C), reality and conservation
of the source, and D·D̄ against the wave operator. U5 on a polynomial
column with cut-complex coefficients. U6 at four rational points. A
polynomial wave and its wrong-handed partner. Seven tests; a third cut
without ι and a potential that is not a reading are rejected.

## 7. Claim boundary

```text
FOUR KINDS OF PART; ι = C₁C₂C₃ ; D·D̄ = WAVE OPERATOR                      PROVED
FIELD = D̄(POTENTIAL); FOUR EQUATIONS = FOUR GRADES                        PROVED
SOURCE SELF-DAGGER AND CONSERVED; MATTER'S READING IS SUCH A SOURCE        PROVED
FORCE = FRAME CHANGE; UNRECOVERABLE PART UNTOUCHED                         PROVED
THE STRENGTH OF THE COUPLING (≈ 1/137)                                     NOT DERIVED — q is an integer; the size of one unit is open
GRAVITY IN THE SAME BLOCK                                                  NOT DONE — the thesis is even between sheets; no block law for it yet
ANYTHING NEW ABOUT ELECTRODYNAMICS                                         NOTHING
```

## 8. Reproduce

```text
python ob1_one_block.py
python -m unittest test_ob1_exact
```
