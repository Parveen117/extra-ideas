# RW1 — The returned winding count: the winding equation is a count of two strands

Monty Dabas. 9 October 2026. Python 3.12, standard library only. Exact rational arithmetic.

Source object: the winding equation of the "uncut horizons" draft,

```text
dw/dh = K − w²/h ,
```

from which the draft builds its "emotional potential" (dw/dh)/√h and its "regret operator" (backward winding).
The draft gives no solution and no meaning for w. This stage supplies both, and uses the result on two
existing stages.

Sources read before building: **RKF F00-E** (native factorial exponential; tail lemma 2.1), **RKF theorum/75-T1**
(ladder law f_(c+½)/f_c ≤ κ/(2(2c+2)), proved termwise on the positive series), **physics CT1** (centre record,
tanh k_c = I₂/I₁; it quotes the face expansion f_j = 2I_(2j+1)(κ)/κ from YM-3/4 and the rate of the chain of
turns from MG1-M4), **SN1/SN2** (spread/mean of a count).

## Setting

Two factorial series, one for each strand:

```text
E(u)·E(v) = Σ uᵃ vᵇ / (a! b!)        a = forward count ,  b = returned count .
```

Sector ν = a − b ≥ 0 (net winding). On the sector, with y = u·v and n = b,

```text
weight  t_n = yⁿ / (n! (n+ν)!) ,     Z_ν(y) = Σ t_n ,     w = ⟨n⟩ ,     V = ⟨n²⟩ − ⟨n⟩² .
```

Every series is dominated by the factorial exponential (F00-E, lemma 2.1). All statements are identities
between series with rational coefficients, or follow from the sign of every coefficient.

## Results

**R1 — the equation is a moment law.** n(n+ν)·t_n = y·t_(n−1), term by term. So

```text
⟨ n (n + ν) ⟩ = y      exactly ,          and with  y·dw/dy = V :        y·dw/dy = y − ν·w − w² .
```

For ν = 0 and y = K·h this is the draft's equation. A power-series solution must start at w(0) = 0 or
w(0) = −ν. With w(0) = 0 it is unique — (m+ν)·a_m = [m = 1] − Σ a_i a_(m−i) fixes every coefficient — and it is
the mean returned count. For ν ≥ 1 the start w(0) = −ν admits no power series (the equation at order ν reads
0 = 1, −1, ¼, −1/36, … ≠ 0).

```text
w = y − y²/2 + y³/3 − 11y⁴/48 + …   (ν = 0) .
```

(For ν = 0 the linear form (y·Z′)′ = Z has the double exponent 0, so its other solutions carry a logarithm at
h = 0: general knowledge, not checked here.)

**R2 — only the product of the two rates enters.** The sector weights depend on u and v through y = u·v
alone. R1 read on both counts is ⟨a·b⟩ = y in every sector: fixing the net winding does not change the mean
product of the two counts.

**R3 — ladder.** Z_ν′ = Z_(ν+1) and Z_ν − y·Z_(ν+2) = (ν+1)·Z_(ν+1), term by term. Hence

```text
y = w_ν · ( ν + 1 + w_(ν+1) ) .
```

**R4 — spread law.** With S_N = C(2N+2ν, N+ν) / (N! (N+2ν)!), both of the following hold as series, and every
coefficient is ≥ 0:

```text
(V − w/2) · Z²  =  Σ  S_N · N(2ν+1) / (4(2N+2ν−1)) · y^N
(w − V)   · Z²  =  Σ  S_N · N(N−1)  / (2(2N+2ν−1)) · y^N
```

Therefore, for every y > 0 and every ν ≥ 0,

```text
½ · mean  <  spread  <  mean ,
0 < V − w/2 < (2ν+1)/8   (ν ≥ 1) ,        0 < V − w/2 < 1/4   (ν = 0) .
```

Spread/mean runs between 1 (few counts) and ½ (many). Neither constant can be improved: spread/mean is
above 0.99 at y = 1/100 and below 0.51 at y = 2500 (ν = 0, certified enclosures); in general the bounded
remainder gives spread/mean − ½ < (2ν+1)/(8w) → 0.

*Proof of the coefficients.* w·Z = y·Z′ and V = y − νw − w². The coefficient of y^N in y·Z² is
Σ_(i+j=N) i(i+ν)·c_i c_j (shift i → i−1), so the first series has coefficient Σ c_i c_j · i(i − j − ½), which
after exchanging i and j is ½·Σ c_i c_j [(i−j)² − N/2]. The numbers c_i c_(N−i) are proportional to
C(N,i)·C(N+2ν, i+ν): i marked among N+ν drawn from 2N+2ν with N marked. Its spread is
N(N+2ν)/(4(2N+2ν−1)) (two binomial sums). Inserting it gives both lines.

**R5 — enclosure of the ladder ratio without a root.**

```text
w² + (ν + ½)·w  <  y  <  w² + (ν + 1)·w .
```

In the usual variables x = 2√y and r = I_(ν+1)(x)/I_ν(x) = 2w/x this reads

```text
2ν + 1  <  x(1 − r²)/r  <  2ν + 2 ,        and exactly     x(1 − r²)/r = 2ν + 2·(spread/mean) .
```

**R6 — use.**

*theorum/75 (face ladder; ν = 2c+1, x = κ, r = f_(c+½)/f_c).*

```text
4c + 3  <  κ(1 − r²)/r  <  4c + 4          for every content c and every coupling κ .
```

T75-T1's law r ≤ κ/(2(2c+2)) follows from R3 (w_(ν+1) ≥ 0); it is one-sided and loses its content once
κ > 4c+4. The line above is two-sided at every coupling.

*CT1 (centre record; c = 0, r = tanh k_c).*

```text
sinh 2k_c = κ / (1 + spread/mean) ,
κ/2  <  sinh 2k_c  <  2κ/3                              for every κ > 0 ,
2κ/3 − κ/(3(κ−2))  <  sinh 2k_c                         for κ > 2 .
```

*The third line.* ν = 1: spread/mean − ½ < 3/(8w) by R4, and w > √(1+y) − 1 > κ/2 − 1 by R5, so
2κ/3 − sinh 2k_c = κ·(spread/mean − ½)/((3/2)(1 + spread/mean)) < κ/(6w) < κ/(3(κ−2)).

CT1-C3's two ends (k_c ≈ κ/4 at strong coupling; k_c = ½ ln(4κ/3) + O(1/κ) at weak coupling, there checked
numerically to κ = 10⁶) are the two sides of one inequality. With 2s < e^(2k_c) < 2s + 1/(2s) for s = sinh 2k_c
the third line gives |k_c − ½ ln(4κ/3)| < 1/(2(κ−2)) for κ ≥ 3: the O(1/κ) with a constant (not the best one;
this last form uses the logarithm and is cross-checked in floats only). CT1-C4's self-dual point
sinh 2k_c = 1 lies in 3/2 < κ < 2, since κ/2 < 1 < 2κ/3 there; the enclosures put it between 1.886 and 1.887.

*MG1-M4 (chain of turns).* 3r/κ < 1 − r² < 4r/κ: a power at every κ.

**R7 — the other sign.** For K < 0 the weights alternate; there is no count, and w leaves at the first zero
of Z₀(−y), y₀ ∈ (36/25, 29/20) (no zero before 36/25, a change of sign after). On the regular branch:

```text
K > 0 :  dw/dh = V/h > 0 for every h ;  w² stays below K·h by V ≥ w/2 ("equilibrium" is not reached)
K = 0 :  w = 0
K < 0 :  dw/dh ≤ K < 0 until w leaves .
```

The sign of dw/dh is the sign of K, not a direction of time.

## What of the draft is not kept

```text
the sign table of E                             on the regular branch the sign is that of K, and for K > 0 it is
                                                never zero or negative (R7) ; the branches with a logarithm at
                                                h = 0 can have either sign ; no direction of time enters
conservation law  d/dh ∫ E dμ = 0               measure and domain not defined ; not derived ; not kept
regret operator  e^(−w²)·d/d(−h) "= reflection"  a derivative is not a reflection ; the returned strand is b above
the emotion tensor and its curvature equation   "Ricci tensor of a potential" is not defined ; not kept
names (emotion, joy, regret)                    not used
```

## What is not shown

- The two-sided bounds on I_(ν+1)/I_ν are known (Amos 1974; Segura 2011 — recalled, not re-read at source).
  The line's part: their reading as ½ < spread/mean < 1 of a returned count, the two positive series of R4,
  and R6.
- Why the draft's equation should hold for a physical winding is not derived here: K, h and "winding" are
  the draft's names. R1 says which count obeys it.
- SN1's shear sector also has spread/mean → ½. The number is the same; no identification is made.
- The sums over N in the proof of R4 are written for every N and checked exactly to N = 60, ν ≤ 6.
- The refined centre bound is checked exactly at nine couplings between 2.1 and 3000; its derivation is the
  two lines above.

## Claim boundary

```text
⟨n(n+ν)⟩ = y ;  y·dw/dy = y − νw − w² ;  REGULAR BRANCH UNIQUE, = MEAN RETURNED COUNT                PROVED
LADDER  y = w_ν(ν + 1 + w_(ν+1))                                                                    PROVED
½·MEAN < SPREAD < MEAN ;  REMAINDER BELOW (2ν+1)/8 ;  CONSTANTS SHARP                               PROVED (series to N = 60 exact ; all N by the written count)
4c+3 < κ(1 − r²)/r < 4c+4  AT EVERY COUPLING ;  κ/2 < sinh 2k_c < 2κ/3                               PROVED
2κ/3 − κ/(3(κ−2)) < sinh 2k_c  (κ > 2) ;  |k_c − ½ln(4κ/3)| < 1/(2(κ−2))  (κ ≥ 3)                   PROVED (written ; the log form cross-checked in floats)
THE DRAFT'S CONSERVATION LAW, REGRET OPERATOR, EMOTION TENSOR                                       NOT KEPT
A PHYSICAL MEANING OF K, h                                                                          NOT DERIVED
```

## Reproduce

```text
python rw1_returned_winding_count.py
python -m unittest test_rw1
```
