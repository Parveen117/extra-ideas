# R5: What a boundary response determines about its depth profile

Monty Dabas · 30 September 2026

R4 connected the seam bond to an existing native source but left the source
profile free. R5 supplies an inverse procedure for **arbitrary positive depth
profiles**, without assuming a two-cell period. It also gives an exact limit:
any finite weak-response jet, even together with an exactly closed bond, can
leave different deeper profiles unresolved.

This is a constructive application of classical continued-fraction inversion.
It is not a new general inversion method, a universal source-selection law,
or a derivation of a physical constant. The contribution in this workspace is
the explicit depth-resolution theorem, its native implementation and error
contract, and a closed-bond ambiguity witness checked against the original
source solver.

## 1. Source, probe and observation contract

Use the unchanged RKF paired-depth source at commit
`3cc5a33b05c16d59c90994ddda69dedc0d392424`. Its native algebra has
`R^2=-I`, `K^2=I`, `KR=-RK`, and `L=KR`. RD1–RD2 give the boundary
inverse corner `F=I+xL` for supplied positive paired cell products `a_j`.
Each cell product occurs on two consecutive edges; individual opening and
closing amplitudes remain unresolved by these products.

Apply a **known common multiplier** `z` to every cell product. For `z>0`,
the source recurrence is

\[
x_j(z)=\frac{z a_j}{1+z a_j x_{j+1}(z)}. \tag{1}
\]

For a finite aperture the final scalar tail is zero. For completed responses
in this packet assume `0<a_j<=Q` for some finite `Q`. The source's RI1–RI2
then supplies a unique response at every positive `z`, independent of the
nonnegative tail. This is an inverse-corner completion, not a full infinite
operator inverse. The location of the boundary aperture is fixed.

Write `t=z^2` and `x_j(z)=z h_j(t)`. Formally, (1) becomes

\[
h_j(t)=\frac{a_j}{1+t a_j h_{j+1}(t)},\qquad h_j(0)=a_j. \tag{2}
\]

The observation is a finite list

\[
h_0(t)=c_0+c_1t+\cdots+c_{m-1}t^{m-1}+O(t^m),
\quad x_0(z)=\sum_{k=0}^{m-1}c_k z^{2k+1}+O(z^{2m+1}). \tag{3}
\]

These are coefficients at **zero coupling**, not response samples or
derivatives at the closed cut `z=1`. If derivatives are supplied, their
conversion is `c_k=x_0^(2k+1)(0)/(2k+1)!`. An actual experiment would need
a warranted common-product probe, gain calibration and coefficient errors.
No such experiment is asserted here.

**Why formal coefficients represent the bounded source near zero.** Each
extra cell in (2) delays the influence of its tail by one power of `t`, so
every coefficient stabilizes after finitely many cells. Write
`h_j=sum_k (-1)^k p_(j,k) t^k`. Multiplying (2) through its denominator gives

\[
p_{j,0}=a_j,\qquad
p_{j,k}=a_j\sum_{r+s=k-1}p_{j,r}p_{j+1,s}\quad(k\ge1).
\]

Induction and the Catalan convolution give
`0<=p_(j,k)<=a_j Cat_k Q^(2k)`. Consequently these series converge
absolutely for `|z|<1/(2Q)`, uniformly in depth after division by `a_j`.
For real `t>=0` in that disk,
`|h_j/a_j-1|<=sum_(k>=1) Cat_k(Q^2 t)^k<1`, so all tails are positive.
They solve (2), hence lie in every native prefix enclosure. RI2 identifies
them with the unique completed response. Thus (3) is an actual Taylor
expansion in this bounded case. It need not converge at `z=1`; the source's
`(3,2)` example already shows that failure. Closure below uses exact source
synthesis, never evaluation of a truncated weak-coupling series at one.

## 2. R5-1: Recover exactly one cell per odd coefficient

**Theorem.** The `m` coefficients in (3) determine exactly the first `m`
positive cell products. The inverse procedure is

\[
\boxed{a_j=h_j(0),\qquad
h_{j+1}(t)=\frac{h_j(t)^{-1}-a_j^{-1}}{t}.} \tag{4}
\]

At each step compute the reciprocal only to the available order, discard
its constant term, and shift the remaining coefficients down by one power.
Exactly one available coefficient is consumed. A nonpositive recovered
constant rejects a strictly positive profile of the requested depth.

**Proof.** Equation (4) is a rearrangement of (2). A formal series with
nonzero constant has a unique reciprocal, whose first `n` coefficients
depend only on its first `n` input coefficients. Therefore `h_0` modulo
`t^m` determines `a_0` and `h_1` modulo `t^(m-1)`, and induction determines
`a_0,...,a_(m-1)` uniquely. Conversely any positive continuation beyond
that prefix first affects `h_0` at order `t^m`, by (2). A finite list passes
this positivity test if and only if it is the `m`-coefficient jet of some
positive profile: reconstruct its positive prefix, attach any bounded
positive tail, and reverse (4) modulo the available powers. This also
supplies a bounded completion. QED.

Examples:

| Supplied odd coefficients of `x_0(z)` | Recovered prefix | Conclusion |
| --- | --- | --- |
| `3, -18` | `3, 2` | Two cell products, with tail free |
| `3, -18, 216` | `3, 2, 3` | Three cell products; periodicity unproved |
| `3, -18, 180` | `3, 2, 2` | A distinguishable third cell |
| `1, -1, 1/2` | Third stripped constant is `-1/2` | Incompatible with a positive third cell, despite alternating coefficient signs |

A zero stripped constant is a boundary of the strictly positive contract.
It can be compatible with a shorter terminated model; finite data do not
establish that physical termination. The implementation refuses to silently
switch to that different model.

## 3. R5-2: The first visible order is sharp

Suppose two profiles agree at depths `0,...,m-1`, and first differ at depth
`m`. Subtracting (2) at a common cell gives the exact formal identity

\[
h_j-\widetilde h_j=
\frac{-t a_j^2(h_{j+1}-\widetilde h_{j+1})}
{(1+t a_j h_{j+1})(1+t a_j\widetilde h_{j+1})}.
\]

All denominator constants are one. Iterating to the boundary proves

\[
\boxed{x_0(z)-\widetilde x_0(z)=
(-1)^m\!\left(\prod_{j<m}a_j^2\right)
(a_m-\widetilde a_m)z^{2m+1}+O(z^{2m+3}).} \tag{5}
\]

The displayed coefficient is nonzero. Thus coefficients through `z^(2m-1)`
cannot see the next cell, and `z^(2m+1)` is the first possible order that
does. Equality of the first `m` odd coefficients is equivalent to equality
of the first `m` products. A full formal series determines all products,
but a finite series never proves an infinite periodic continuation.

This is an exact identifiability statement. It is not a stability guarantee
for high derivatives estimated from noisy readings. The product in (5),
coefficient cancellations, and probe calibration can all affect conditioning.

## 4. R5-3: Exact closed bonds do not remove finite-jet ambiguity

The limitation survives the additional observation `x_0(1)=1`. Start with
the existing closed source `(3,2,3,2,...)` and retain its first `m` cells.
The required response at the new tail boundary is

\[
\tau_m=\begin{cases}1,&m\text{ even},\\2/3,&m\text{ odd}.\end{cases}
\]

RKF SY1 already proves that the period `(tau+b tau^2,b)` has response
`tau` at `z=1` for every `b>0`. Choose `b=1`. Attach the period `(2,1)`
when `m` is even, or `(10/9,1)` when `m` is odd. All products are positive
and at most three. The native recurrence through the unchanged prefix
gives `x_0(1)=1` exactly for both profiles. Their first `m` odd coefficients
coincide, but their next cell products differ, so (5) distinguishes the next
coefficient.

This constructs such a pair for **every finite `m`**. It proves that a
finite zero-coupling jet plus exact completed cut closure does not select
the entire profile, even within uniformly bounded positive sources.

For the concrete choice `m=2`:

| Quantity | Original source | Changed source |
| --- | --- | --- |
| Full cell sequence | `(3,2,3,2,...)` | `(3,2,2,1,2,1,...)` |
| Completed `x_0(1)` | `1` | `1` |
| Coefficient of `z` | `3` | `3` |
| Coefficient of `z^3` | `-18` | `-18` |
| Coefficient of `z^5` | `216` | `180` |
| Completed `x_0(1/2)` | `(sqrt(7)-1)/2` | `6-3 sqrt(3)` |
| Outward rational enclosure at `z=1/2` | `[0.822875655, 0.822875656]` | `[0.803847577, 0.803847578]` |

The original profile at `z=1/2` satisfies `2x^2+2x-3=0`. For the changed
profile, its `(2,1)` tail has response `sqrt(3)-1`; applying the prefix
with scaled cells `(3/2,1)` gives `6-3 sqrt(3)`, satisfying
`x^2-12x+9=0`. The native solver produces disjoint exact rational intervals;
independent polynomial sign checks identify the indicated roots. The table
uses outward rounding of those intervals, not floating-point fits.

Both sources give the same R4 bond `B=(I+L)/2`. A further response probe
distinguishes them. This witness does **not** claim to reproduce NI's two
specific symmetric readings. NI's reconstruction remains valid within its
declared globally two-periodic family; R5 studies a larger profile class
and a different, calibrated observation protocol.

## 5. R5-4: Propagate the unresolved tail into a prediction interval

After recovering a prefix, a future response prediction must keep its
unobserved tail. Suppose the first unobserved cell is independently bounded
by `a_m<=Q`, and deeper responses are nonnegative. Equation (1) implies
`0<=x_m(z)<=Qz`. Let the source transfer matrix for the scaled prefix be

\[
M_m(z)=\prod_{j<m}\begin{pmatrix}0&1\\1&(za_j)^{-1}\end{pmatrix}
=\begin{pmatrix}A&B\\C&D\end{pmatrix},
\quad\Phi_m(u)=\frac{Au+B}{Cu+D}.
\]

For `m>=1`, `C,D>0` and `det M_m=(-1)^m`. The predicted response lies
between `Phi_m(0)` and `Phi_m(Qz)`, sorted to account for parity. Reusing
the source's RI1 determinant identity gives the exact width of this
transported tail interval:

\[
\boxed{W=\frac{Qz}{D(CQz+D)}
\ \le\ Q\!\left(\prod_{j<m}a_j^2\right)z^{2m+1}.} \tag{6}
\]

For the inequality, subtraction of two instances of (1) gives a change
at most `(za_j)^2` times the tail change: both positive denominators are
at least one. Multiply these bounds through the prefix, starting with
tail width `Qz`. No derivative or depth-limit exchange is used.

The exact interval is evaluated by the existing native `prefixEnclosure`,
not a copied forward solver. The last inequality states how weak probing
hides deeper variation. `Q` is an additional tail prior, not something
learned from the finite jet. Without it the existing native `[0,infinity]`
tail enclosure remains available. Neither interval includes measurement
error in the supplied coefficients; that needs its own observation contract.

## 6. Calibration and scope

Known probe and readout scales matter here. If the observed function is
`y(z)=eta x_0(g z)` with `eta,g>0`, it is itself the boundary response of
the transformed profile

\[
\boxed{a'_j=g\,\eta^{(-1)^j}a_j.} \tag{7}
\]

Indeed set `x'_j(z)=eta^((-1)^j) x_j(gz)`. Adjacent factors of `eta`
cancel in the denominator of (1), proving (7), including finite zero tails.
Unknown probe gain therefore changes the common product scale, while
unknown readout gain alternately rescales even and odd depths. The code
checks both the native response identity and the recovered coefficients.
Unlike NI, R5 does not claim to cancel those gains through a symmetric
two-probe protocol.

Recovered products still do not specify individual directional amplitudes:
the existing freedom `b_h -> q_h b_h`, `c_h -> c_h/q_h` preserves them.
No physical time, frequency, photon, charge, action normalization or numerical
`lambda_*` is supplied. The boundedness assumption justifies the completed
Taylor response used here. For unbounded profiles, formal stripping still
has a coefficientwise meaning, but convergence and boundary selection
require the separate source criteria and are not inferred from it.

## 7. Lineage and what is added

The immutable [R4 lineage assessment](../CROSS_REPO_LINEAGE.md) remains
unchanged. The following comparison extends it for R5. Source refs and hashes
are recorded in [R5_SOURCE_PINS.json](../04-operator-evolution/R5_SOURCE_PINS.json).

| Existing work | Relationship to R5 |
| --- | --- |
| RKF RD1–RD2, RI1–RI2 | Native recurrence, finite inverse corners, transfer identity and bounded-profile completion are dependencies |
| RKF SY1–SY2 | The already proved target synthesis constructs the equal-cut tails in R5-3; it is not reintroduced as new |
| RKF ID1 | Recovers two local unpaired products with a known tail; R5 uses a weak-response jet at one boundary to recover a positive paired prefix with unknown deeper tail |
| Publications AS-1–AS-2 | Already distinguishes two-cell closed profiles by response derivatives and generates their fixed-model jets |
| Publications NI-2–NI-4 | Already identifies the parameter and probe gain inside a globally two-periodic model; R5 instead permits arbitrary depth dependence and requires calibrated coefficients |
| Publications CR-1–CR-4 | Already derives an inverse pole and a third-probe prediction; R5 does not claim either as new |
| Classical Stieltjes fractions and power-series inversion | The mathematical inversion method in (4) is established; R5 specializes it to the native source and supplies the depth and closure statements above |

For the last row, expanding (2) gives the Stieltjes-type fraction with
leading factor `a_0` and successive numerator weights `t a_j a_(j+1)`.
NIST [DLMF 3.10(ii)](https://dlmf.nist.gov/3.10.ii) documents the
power-series/continued-fraction correspondence and coefficient recovery.
Alan D. Sokal, [*A simple algorithm for expanding a power series as a
continued fraction*](https://arxiv.org/abs/2206.15434v2), Expositiones
Mathematicae 41 (2023), 245–287, discusses the older algorithmic lineage.
The argument here is self-contained; external novelty of this integration
is not established. The repository comparison is limited to the pinned
source packets, not an exhaustive audit of every branch.

The useful advance is operational: a response list can now yield a precise
recoverable prefix, reject an impossible positive profile, and carry its
remaining tail into a prediction interval. The all-depth construction in
R5-3 prevents a finite successful fit and exact closure from being promoted
to a universal source law.

## 8. Implementation and reproducibility

The [exact recovery module](../04-operator-evolution/profile_recovery.py)
accepts integer or rational coefficients, and explicitly reports an
unresolved tail. From the repository root:

```bash
python3.12 -B 04-operator-evolution/profile_recovery.py '[3, -18, 216]'
python3.12 -B 04-operator-evolution/verify_r5.py --rkf-root ../rkf-r4
```

Use Python 3.11 or 3.12 and Node.js, with the unchanged RKF runtime pinned
by R4. The R5 verifier checks those bytes before execution. Its 14 focused
tests cover exact recovery, nine polynomial/native finite-inverse
comparisons, minimal visible order, positivity refusals, both calibration
gains, source-designed closed tails, disjoint held-out intervals, transferred
tail bounds, and four native Schur proof replays. Earlier proof/code hashes
are checked without rerunning unrelated suites. Written proofs, finite
executable checks and empirical validation remain distinct evidence.
