# SD1 — T24's S = R + D on the objects of the physics line

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic and sympy.

OR1, MM1, GF1 and CS1 reduced the premises of the line to "no memory in a shared record". This stage ties that
statement to the framework's own primitive law and checks it object by object.

Sources read before building: RKF **theorum/24** §2–§6 (event lift Z, cut pair P + Q = 1; S = Z*Z, R = Z*PZ,
D = Z*QZ, F = R − D; **Theorem 6.1**: S = R + D, F = R − D, R ≥ 0, D ≥ 0, all Gram forms), Publications **GE2-T2**
(eq. 8: M = A†QA = I − S†S), **IN1-T1/T2**, **LN1-N1**, **SE1**, **MM1**, **GF1**, **TP1**.

## Results

**(a) A record of turns.** Z f = (u_s f)_s on the weighted sum over records, P = JJ*. Then S = 1, R = det of the
record mean, and D = 1 − det of the mean: GE2's memory is T24's D. OR1's closure is D = 0.

**(b) A single reading.** Z = ψ, cut i: P_i = ½(1 + C_i). Then S = n, R_i = ½(n + r_i), D_i = ½(n − r_i) and
F_i = r_i: IN1's reading tensor is (S; F₁, F₂, F₃). For a single reading S² = ΣF_i², and 4·R_i·D_i is the sum of
the other cuts' F²: what one cut loses is what the others observe (IN1-T2).

**(c) The stretch of the fall frame.** Z = K, cut i: P_i = e_i e_iᵀ. Then S_i = (K²)_ii, R_i = K_ii²,
D_i = Σ_{j≠i} K_ij²: LN1-N1 is Theorem 6.1. SE1's energy law reads

```text
2·e₂(K) = ( tr Z )² − Σ_i S_i ,        empty space:   Σ_i S_i = ( tr Z )² ,   i.e.   Σ_i D_i = 2·Σ_{i<j} K_ii K_jj .
```

The total source over the cuts equals the square of the uncut reading.

**(d) Any frame.** TP1's law is Q = G − I₃ with G = ¼I₁ + ½I₂ a quadratic form blind to turn rates and I₃ the
square of the uncut (trace) reading of the order defect. On an order defect in one plane G = I₃.

## What this settles

```text
"why is the law quadratic"      S, R, D, F are Gram forms of the event lift (T24); the laws of the line compare
                                the total source with the square of an uncut reading
the record rule                 D = 0 : the object is fully recognized. It is T24's own notion, not a new premise
```

## What is not shown

- (a)–(c) are instances of Theorem 6.1 with the event lift and the cut named. (d) is a rewriting of GF1; the general
  frame law is not exhibited as a single R − D with one cut pair (its three parts enter with weights ½, −1, −¼).
- Nothing new is predicted here.

## Claim boundary

```text
GE2's MEMORY = T24's D FOR A RECORD OF TURNS                                    PROVED (exact instance)
IN1's READING TENSOR = (S; F_i) ; 4 R_i D_i = OTHER CUTS' F²                     PROVED
LN1-N1 = THEOREM 6.1 ON THE STRETCH ; ENERGY LAW: Σ S_i = (tr Z)²                PROVED
ANY FRAME: Q = (TURN-BLIND GRAM PART) − (TRACE READING)²                        PROVED
THE GENERAL FRAME LAW AS ONE R − D                                              NOT SHOWN
```

## Reproduce

```text
python sd1_source_recognized_memory.py
python -m unittest test_sd1
```
