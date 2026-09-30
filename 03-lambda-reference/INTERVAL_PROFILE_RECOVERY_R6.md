# R6: Recover depth couplings with coefficient error bars

Monty Dabas · 30 September 2026

R5 recovered a depth prefix from exact weak-response coefficients. R6 makes
that inverse procedure usable with **supplied error intervals**. It encloses
the couplings that the data determine, stops where positivity cannot be
resolved, and propagates both coupling uncertainty and hidden-tail freedom
into a future response interval.

The input contract remains R5's calibrated common product multiplier `z`
and coefficients `c_k` of `z^(2k+1)` at zero coupling. The positive paired
native source and its bounded-profile completion remain dependencies. The
coefficient intervals must already account for the relevant observation and
coefficient-extraction errors. R6 does not establish those laboratory bounds
or remove unknown calibration gains.

## R6-1. Sound interval stripping and a depth stopping rule

Suppose each coefficient lies in a finite rational interval `C_k`. Real
coefficients may be correlated; their joint feasible set is assumed contained
in the supplied Cartesian box. Replace scalar operations in the R5 inverse
by enclosing interval operations. If the current constant interval `H_0`
is strictly positive, the reciprocal coefficients are enclosed recursively:

\[
V_0=1/H_0,\qquad
V_k=-V_0\sum_{i=1}^{k}H_iV_{k-i}. \tag{1}
\]

Record `H_0` as the next product interval and continue with `V_1,V_2,...`.
The reciprocal constant is discarded, exactly as in R5. All additions,
products and reciprocal endpoints use exact rational arithmetic.

**Theorem.** Every positive profile consistent with the input coefficient
box has each recovered cell in its recorded interval. At depth `j`:

- If the attempted interval has positive lower endpoint, depth `j` is certified.
- If its upper endpoint is nonpositive, no strictly positive profile of the
  requested depth fits any coefficient vector in the box.
- If it contains zero and positive values, positivity is unresolved. Keep
  the already certified prefix and stop before dividing by that interval.

**Proof.** The reciprocal of a positive scalar in `[l,u]` lies in
`[1/u,1/l]`; interval addition and multiplication contain all scalar
results. Induction on `k` proves inclusion in (1), and induction on stripping
depth proves the asserted coupling inclusion using R5's exact inverse.
Every earlier constant was certified positive, so the next scalar reciprocal
and stripped constant exist. A nonpositive upper endpoint therefore rules
out the next strictly positive product for every input realization. An
interval crossing zero supplies no such decision. QED.

Completing all `m` stages certifies a positive `m`-cell prefix for every
coefficient vector in the input box. R5 then gives a bounded positive
continuation for each vector. The cell intervals are generally conservative:
they can contain product combinations that no single original coefficient
vector realizes. A partial result or an incompatible result does not prove
that a physical source terminates. Increased precision can resolve an
undecided sign, but interval dependence can also cause overestimation.

For the first two coefficients with `0<L_0<=c_0<=U_0` and
`L_1<=c_1<=U_1<0`, the exact box ranges are already explicit:

\[
a_0\in[L_0,U_0],\qquad
a_1=-c_1/c_0^2\in[-U_1/U_0^2,-L_1/L_0^2]. \tag{2}
\]

This is a useful calibration check on the general algorithm.

## R6-2. An uncertainty-aware source distinction

Give each of the first three coefficients an independent absolute error
bound `1/100`. Apply R6-1 to the two R5 fixtures:

| Central coefficients | Certified third-cell enclosure, widened outward | Result |
| --- | --- | --- |
| `(3,-18,216)` | `[2.90,3.10]` | Contains the original third product `3` |
| `(3,-18,180)` | `[1.90,2.10]` | Contains the changed third product `2` |

The exact rational algorithm produces intervals inside these displayed
outward bounds. They are disjoint, so this coefficient-error contract still
distinguishes the third cell without a point fit. The first two couplings
and the unknown continuation remain separately reported. These are designed
arithmetic data, not readings from a performed physical experiment.

There is also an exact stopping witness. Fix `c_0=3`, `c_1=-18` and take
`c_2 in [107,109]`. The first two products are `3,2`, while the attempted
third product lies in `[-1/36,1/36]`. The tool returns
`UNRESOLVED_POSITIVITY` at depth two. With `c_2 in [100,107]`, its entire
third-product interval is negative, so the positive-depth model is
incompatible. A midpoint guess would erase this distinction.

## R6-3. Predict with both uncertain couplings and the remaining tail

Let the recovered positive prefix be enclosed by independent boxes
`a_j in [l_j,u_j]`, `j=0,...,m-1`. Supply a bound `Q>0` on the **first
unrecovered** product and retain nonnegative deeper responses. At probe
`z>0`, the source recurrence gives a terminal response in `[0,zQ]`.
Propagate an interval `[v,w]` outward through one cell by

\[
\boxed{[v,w]\ \mapsto\left[
\frac{zl_j}{1+zl_j w},
\frac{zu_j}{1+zu_j v}\right].} \tag{3}
\]

**Theorem.** The final interval contains every boundary response in the
declared coupling-box and tail contract. Its endpoints are the exact extrema
over that independent Cartesian box and tail interval.

**Proof.** For `f(a,x)=za/(1+za x)`, direct differentiation gives
`f_a=z/(1+za x)^2>0` and `f_x=-(za)^2/(1+za x)^2<0`. Thus (3) is the
exact one-cell range. At the boundary, dependence on the product at depth
`j` has sign `(-1)^j`; dependence on the final tail has sign `(-1)^m`.
The lower corner uses `l_j` at even depths and `u_j` at odd depths, with
tail `0` for even `m` and `zQ` for odd `m`. Reverse all choices for the
upper corner. These simultaneous corners attain the recursive endpoints,
proving exactness for the box. QED.

Recovered coupling intervals can be correlated, so the prediction is an
enclosure of their actual feasible responses rather than a claimed sharp
range for the original observation set. This preserves both errors that
R5 separated: imperfect coefficient knowledge and an unobserved tail.
The two corners are evaluated independently by the unchanged native solver,
with both even and odd prefix lengths checked.

`Q` is a supplied prior. On a partial recovery it refers to the first cell
at the reported stopping depth, not to the originally requested final depth.
The CLI reports the prediction prefix length and refuses a prediction from
data already certified incompatible. Infinite periodicity, full source
selection and a universal value of `lambda_*` still require an additional law.

## Implementation, evidence and lineage

[interval_profile_recovery.py](../04-operator-evolution/interval_profile_recovery.py)
implements the rational intervals, stripping statuses and prediction bounds.
For example:

```bash
python3.12 -B 04-operator-evolution/interval_profile_recovery.py \
  '[["2.99","3.01"],["-18.01","-17.99"],["215.99","216.01"]]' \
  --probe 1/10 --tail-bound 3
python3.12 -B 04-operator-evolution/verify_r6.py --rkf-root ../rkf-r4
```

Use Python 3.11 or 3.12 and the pinned R4 runtime. R6 reuses R5's native
adapter and the original RKF solver; their bytes are unchanged. The focused
checks cover exact-point recovery, all coefficient-box corners in the fixture,
noise discrimination, undecided and incompatible depths, narrowing input
errors, eight native prediction-corner comparisons, and completed native
responses inside the predicted intervals. The verification record also
checks every earlier R1–R5 proof/code hash.

R5's continued-fraction lineage and [source pins](../04-operator-evolution/R5_SOURCE_PINS.json)
remain authoritative dependencies. Publications NI-4 already propagates
uncertain symmetric readings to the single parameter of a declared
two-periodic model. R6 instead propagates coefficient intervals through an
arbitrary depth-prefix inverse and reports where that recovery stops.
General interval inclusion is established mathematics; see Siegfried M.
Rump, [*Verification methods: Rigorous results using floating-point
arithmetic*](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf), Acta
Numerica (2010), for interval inclusion and dependency overestimation.
R6's rational implementation does not depend on INTLAB. External novelty
of this integration is not claimed.

The added capability is a certified recovery boundary from supplied uncertain
data, followed by a native response prediction that keeps its remaining
uncertainty. Establishing experimental coefficient bounds and gain calibration
is the next observation-side task.
