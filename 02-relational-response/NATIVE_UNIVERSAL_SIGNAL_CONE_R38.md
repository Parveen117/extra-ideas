# R38: the sharp native cone for reliably detected localized signals

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R37 constructed compact source packets with speeds approaching
`c_Sigma=1/sqrt(2)`. R38 proves the matching upper result for **every
initially localized preparation of the same source**, including
preparations that change with the observation scale. The upper theorem
requires no chosen carrier, band, profile or source role.

If the initial matching norm is concentrated within radius `o(n)` of
its preparation centre x_n, then for every fixed eta>0,

\[
\boxed{\frac{\|1_{\{|x-x_n|>(c_\Sigma+\eta)n\}}Z^n f_n\|^2}
                    {\|f_n\|^2}\longrightarrow0.}
\tag{38.A}
\]

Compact support gives an explicit exponential outer-tail bound. A
detector at increasing distance cannot receive a fixed positive fraction
of the prepared matching norm at asymptotic speed greater than c_Sigma.
Together with R37, this makes c_Sigma the **sharp operational speed
supremum for this source and this detector convention**.

Finite cyclic phase resolution is derived from native roots and count
cancellation. No Fourier, spectral or propagation-bound theorem is a
premise. Positive weights are estimation tools, not a new source law.
The faster nonzero support front remains, and supplies a counterexample
for vanishing detection thresholds. Physical light and a universal speed
for other source laws are not identified.

## Source and target

Use R35's unchanged sixteen-role Laurent source
`Z=sum_d Z_d S^d`, supported on `d in {-1,0,1}^2`, with
`(S^d f)(x)=f(x-d)`, native matching and unitary continuation. It obeys

\[
Z+Z^\dagger=2I-L_*,\qquad
L_*=(\lambda_x+\lambda_y)/2-\lambda_x\lambda_y/16.
\tag{38.0}
\]

All scalars, positive roots and role tuples descend from R16/R20. Use
their derived matching completion. The length `|x|=sqrt(x_1^2+x_2^2)`
is component count length; R36's dual information length is sqrt(2)
times this. A component step contains four fine address increments;
one Z block contains eight original source events. Let
`E(t)=Exp_Sigma(t)` be F00-E's earned factorial exponential.

For a matrix put `m(A)=sum_ij(|rad A_ij|+|turn A_ij|)`. The complete
source coefficient ledger gives

\[
\mu=\sum_d m(Z_d)=68,\quad
\sum_d m(Z_d)|d|_1=64,\quad
\sum_d m(Z_d)|d|_1^2=88.
\tag{38.1}
\]

These conservative bounds are not optimized costs or physical couplings.

## Written results

### R38.1 — Native dyadic turns derive exact finite phase resolution

For every power of two Q>=4, construct a primitive native unit turn
xi_Q. On a Q-by-Q address ledger define

\[
\widehat f(j,k)=\sum_{x,y=0}^{Q-1}\xi_Q^{jx+ky}f(x,y).
\]

Then the following are derived finite identities:

\[
f(x,y)=Q^{-2}\sum_{j,k}\xi_Q^{-jx-ky}\widehat f(j,k),
\qquad
\boxed{\sum_{x,y}\|f(x,y)\|^2
       =Q^{-2}\sum_{j,k}\|\widehat f(j,k)\|^2.}
\tag{38.2}
\]

For a finite Laurent matrix A, cyclic continuation gives
`widehat(Af)(j,k)=A(xi_Q^j,xi_Q^k) widehat f(j,k)`.
Thus a bound `||A(z_1,z_2)||<=K` on native unit turns implies
`||Af||<=K||f||` on unwrapped finite fields and their norm completion.

**Proof.** Start with xi_4=iota. For xi_Q=a+iota b in the positive
quarter sector, define

\[
\xi_{2Q}=\sqrt{(1+a)/2}+\iota\sqrt{(1-a)/2}.
\]

The positive native roots give unit norm and square xi_Q: the turn
part is positive and its square is `1-a^2=b^2`. Inductively
`xi_(2Q)^Q=-1`. If a smaller positive exponent d gave one, squaring
would make Q divide d by the parent's primitivity, leaving only d=Q,
a contradiction. Division of integer exponents proves the divisibility
step directly. This constructs all needed finite phase alphabets.

For d not divisible by Q, finite geometric cancellation gives
`(1-xi_Q^d) sum_(j=0)^(Q-1) xi_Q^(jd)=0`; the first factor is
nonzero, so the sum is zero. If Q divides d the sum is Q. Expand
the inverse formula and squared matching sums. Both address differences
must vanish modulo Q, proving (38.2). Changing the summed address by
a shift proves the evaluation rule. Apply the phasewise bound inside
the finite norm identity.

For finite unwrapped f and A, choose a dyadic Q exceeding a common
coordinate span containing their input and output supports. Reduction
modulo Q then introduces no collisions in either norm. The cyclic
bound is exactly the unwrapped bound. Finite truncation and native
Cauchy completion give the extension. No circle measure, integral
transform or predeclared completeness theorem has been used.

### R38.2 — Repeated source blocks have a uniform sharp leading displacement bound

For a native unit direction e define the count derivation
`D_e A(z)=sum_d (e dot d) A_d z^d`. Addition of displacement counts
gives its product rule. For every integer m>=1 and every unit-turn pair z,

\[
\boxed{\|\mathcal D_e Z^m(z)\|
             \le m c_\Sigma+176\sqrt m.}
\tag{38.3}
\]

R38.1 gives the same unwrapped operator bound. It includes the touching
point Z(1,1)=I.

**Proof.** Put `r=(1-rad z_1)/2`, `s=(1-rad z_2)/2`,
`ell=2(r+s)-rs`, `a=1-ell/2`, `h^2=ell(4-ell)/4`.
R35's identities apply to any unit turn: their proof uses its unit-square
identity and local multiplication by an Euler turn, not a global angle.

For ell>0 let `B=Z-aI` and `U=D_e Z`. Then
`B^dagger=-B`, `B^2=-h^2 I`. The R35 orthogonal branch cuts give

\[
U_{\rm d}=\tfrac12(U-BUB/h^2),\qquad
U_{\rm o}=\tfrac12(U+BUB/h^2).
\tag{38.4}
\]

The first preserves each branch and the second switches them. If ell_e
is the local Euler derivative in direction e, differentiating the source
quadratic or its branch equations gives

\[
U_{\rm d}=-\frac{\iota\ell_e}{2h^2}ZB,
\qquad U_{\rm d}^\dagger U_{\rm d}
                 =\frac{\ell_e^2}{4h^2}I\le\tfrac12 I.
\tag{38.5}
\]

The inequality is the directional consequence of R35.4's complete
positive speed polynomial: for its velocity v and unit e,
`|v|^2-(v dot e)^2=(v_1 e_2-v_2 e_1)^2>=0`.
Hence `||U_d||<=c_Sigma`.
The two off-diagonal blocks have orthogonal inputs and outputs and
each is a compression of U, so `||U_o||<=||U||<=64` by (38.1).
Interference with the diagonal part of the full derivative is retained.

Define `s_0=0`, `s_1=1`, `s_(m+1)=2a s_m-s_(m-1)`. Product
differentiation and finite branch geometric cancellation give

\[
\mathcal D_e Z^m=mZ^{m-1}U_{\rm d}+s_m U_{\rm o},\qquad
s_m=\frac{\zeta_+^m-\zeta_-^m}{\zeta_+-\zeta_-}.
\tag{38.6}
\]

Both roots have unit norm, so `|s_m|<=1/h`. For ell>=delta^2,
`2h=sqrt(ell(4-ell))>=delta` because ell<=3. Therefore

\[
\|\mathcal D_e Z^m\|\le mc_\Sigma+128/\delta
                            \le mc_\Sigma+176/\delta.
\tag{38.7}
\]

Near the touching point no inverse gap is used. If ell<=delta^2,
then `r+s<=ell` and `|z_i-1|<=2delta`. Finite unit-turn telescoping
gives `|z^d-1|<=2delta |d|_1`. R35.2 gives
`D_e Z(1,1)=e_1 A_x+e_2 A_y` with square I/2. Thus (38.1) gives
`||D_e Z(z)||<=c_Sigma+176delta`. Differentiating m unitary factors
bounds the m-block derivative by `m(c_Sigma+176delta)`. Choose
`delta=1/sqrt(m)` and combine the regions to prove (38.3).
At ell=0 the direct derivative is `m D_e Z(1,1)`; division by the
zero gap has not been used.

The block estimate cannot be replaced by a universal single-block
bound c_Sigma. At `z_1=(3+4 iota)/5`, `z_2=(5+12 iota)/13`,
`e=(3/5,4/5)`, and the first unit source role v, exact source arithmetic
gives `||D_e Z(z) v||^2=1629/3250>1/2`. The cancellation in (38.6),
rather than an assumed instantaneous speed cap, controls long flight.

### R38.3 — Native positive count weights give a controlled block estimate

For beta>0 and a finite field set

\[
(W_{\beta,e}f)(x)=E(\beta e\cdot x)f(x),\qquad
Z_{\beta,e}=\sum_d E(\beta e\cdot d)Z_dS^d.
\]

Source composition gives `W_(beta,e) Z=Z_(beta,e) W_(beta,e)`.
For `0<beta<=1/(2m)`,

\[
\boxed{\|Z_{\beta,e}^{\,m}\|
 \le E\!\left(\beta m
       [c_\Sigma+176/\sqrt m+4m68^m\beta]\right).}
\tag{38.8}
\]

The weight is used on finite fields; it is not declared bounded on the
whole matching completion.

**Proof.** The factorial series gives E(t)>0 for t>=0, and its derived
product law gives `E(-t)=1/E(t)>0`. Consequently E is increasing;
`1+t<=E(t)` for t>=0. For |t|<=1 its factorial tail gives

\[
|E(t)-1-t|\le t^2,
\tag{38.9}
\]

since `sum_(j>=2)1/j!<=1`, using the native factorial/geometric bound
already proved in R35. Also `E(-t)<=1/(1+t)` for t>=0, so it tends
to zero at increasing positive t.

Let A=Z^m. Its support lies in `[-m,m]^2`; coefficient gauge and matrix
multiplication give `sum_d m(A_d)<=68^m`. A unit e satisfies
`|e dot d|<=|d|_1<=2m`. Tilting every coefficient of A gives
`Z_(beta,e)^m` by the product law. Equation (38.9) therefore yields

\[
\|Z_{\beta,e}^{\,m}-Z^m-\beta\mathcal D_eZ^m\|
                      \le4m^2 68^m\beta^2.
\tag{38.10}
\]

Every matrix unit and shift has norm one, so coefficient gauges bound
this remainder directly. Insert (38.3), use `||Z^m||=1`, and bound
one plus the resulting nonnegative sum by its native exponential.
This proves (38.8) without an imported large-deviation or spectral-radius
theorem.

### R38.4 — Every compact preparation has an explicit half-space tail bound

For eta>0 choose a positive integer k and radial q,beta by

\[
k\ge704/\eta,\quad m=k^2,\quad
q=\min\{1/(2m),\eta/(16m)\},\quad \beta=q68^{-m}.
\tag{38.11}
\]

Set `b=c_Sigma+eta/2` and `C=68^m E(2m beta)`. For every n>=0,

\[
\|Z_{\beta,e}^{\,n}\|\le C E(\beta bn).
\tag{38.12}
\]

For nonzero finite f supported in `e dot (x-x_0)<=R`, and any radial d,

\[
\boxed{\frac{\|1_{\{e\cdot(x-x_0)\ge d\}}Z^nf\|^2}{\|f\|^2}
       \le C^2 E\!\left(-2\beta[d-R-bn]\right).}
\tag{38.13}
\]

Constants are independent of packet shape, direction, phases and roles.
They are conservative; a small-event optimal estimate is not claimed.

**Proof.** Equation (38.11) gives `beta<=1/(2m)`,
`176/sqrt(m)<=eta/4` and `4m68^m beta=4mq<=eta/4`.
Each full m-block in (38.8) costs at most `E(beta b m)`.
Write n=lm+r with 0<=r<m. Direct coefficient mass gives
`||Z_(beta,e)||<=68 E(2beta)`, so the r remaining blocks cost at most C.
Multiply and replace lm by n to prove (38.12).

Use the weight centred at x_0. Its norm on f is at most
`E(beta R)||f||`; its multiplier on the detected half-space is at
least E(beta d). The exact conjugation identity and (38.12) give
(38.13) on squaring. No unbounded-weight operation on an arbitrary
completed field is needed. Density will extend the bounded propagation
conclusions below.

The powers 68^m may be retained as exact integer expressions. The
certificate checks parameter budgets without expanding huge integers
or presenting symbolic bounds as simulations of those event counts.

### R38.5 — A finite direction ledger gives the universal radial cone

Fix eta>0 and take beta,C from R38.4. Choose an integer Q>=2 with

\[
Q^2\ge4(c_\Sigma+\eta)/\eta,\qquad D_Q=(2Q+1)^2-1.
\tag{38.14}
\]

For support in `|x-x_0|<=R`, whenever `R<=eta n/8`,

\[
\boxed{\frac{\|1_{\{|x-x_0|>(c_\Sigma+\eta)n\}}Z^nf\|^2}
                    {\|f\|^2}
             \le D_Q C^2 E(-\beta\eta n/4).}
\tag{38.15}
\]

More generally let nonzero completed preparations f_n have centres x_n
and radii R_n such that

\[
R_n/n\longrightarrow0,\qquad
q_n=\frac{\|1_{\{|x-x_n|>R_n\}}f_n\|}{\|f_n\|}
                                                 \longrightarrow0.
\tag{38.16}
\]

Then (38.A) follows, allowing arbitrary scale-dependent roles, phases
and profiles without a fixed carrier or spectral gap.

**Proof.** Use directions `e_p=p/|p|` for all nonzero integer p in
`[-Q,Q]^2`. For a unit vector w, round Qw coordinatewise to p. Then
p is nonzero, `|p/Q-w|<=1/(sqrt(2)Q)` and the reverse triangle
inequality gives `|e_p-w|<=sqrt(2)/Q`. Squaring gives
`e_p dot w>=1-1/Q^2=gamma`.

The exterior of radius `(c_Sigma+eta)n` is therefore covered by D_Q
half-spaces at projected distance `gamma(c_Sigma+eta)n`. Adding their
nonnegative site norms bounds the exterior norm square; overlap only
increases the sum. Equation (38.14) gives
`gamma(c_Sigma+eta)>=c_Sigma+3eta/4`. Subtract R and
`bn=(c_Sigma+eta/2)n`; at least eta n/8 remains. Apply (38.13) and
sum to prove (38.15).

For (38.16), split f_n into its finite ball truncation and remainder.
Unitarity and detector contraction bound the propagated remainder by
`q_n||f_n||`. The normalized exterior norm is at most

\[
\sqrt{D_Q}\,C E(-\beta\eta n/8)+q_n.
\tag{38.17}
\]

Both terms vanish. Squaring proves (38.A), also for completed fields.
Localization is essential: information already prepared far from x_n
cannot be counted as information transported there from x_n.

### R38.6 — Fixed-threshold arrival has the sharp operational speed supremum

Let preparations f_j be concentrated about x_j within radius R_j, with
relative outside norm q_j. Let an address detector D_j have distance
at least d_j from x_j, where `d_j->infinity`, `R_j/d_j->0` and
`q_j->0`. Fix 0<theta<1. If a detected event count n_j obeys

\[
\|1_{D_j}Z^{n_j}f_j\|^2\ge\theta\|f_j\|^2,
\]

then

\[
\boxed{\limsup_j d_j/n_j\le c_\Sigma.}
\tag{38.18}
\]

For every such fixed theta the bound is sharp as a supremum: R37
supplies reliable compact packet speeds tending to c_Sigma.

**Proof.** Fix eta>0 and the constants above, with gamma=1-1/Q^2.
The detector lies in the union of half-spaces at distance gamma d_j.
For all integers `0<=n<=d_j/(c_Sigma+eta)`,

\[
\gamma d_j-R_j-bn
 \ge\frac{\eta d_j}{4(c_\Sigma+\eta)}-R_j
 \ge\frac{\eta d_j}{8(c_\Sigma+\eta)}
\]

eventually. The second step uses `R_j/d_j->0`. The half-space bound,
finite covering and initial-tail split give, uniformly over all those n,

\[
\frac{\|1_{D_j}Z^nf_j\|}{\|f_j\|}
 \le\sqrt{D_Q}C E\!\left(-\frac{\beta\eta d_j}
                                   {8(c_\Sigma+\eta)}\right)+q_j
 \longrightarrow0.
\tag{38.19}
\]

Eventually this is below sqrt(theta), excluding detection throughout
that early interval. Since eta is arbitrary, (38.18) follows. No sum
over event times or measurement-induced update is assumed: the readout
is of unchanged source continuation.

For sharpness, R37's fixed-carrier packet has width O(L^2), detector
distance `floor(g L^3)`, first-threshold arrival `L^3+O(L^2)` and
detector fraction tending to one. Thus its speed tends to g. Its native
rational family `2/3,12/17,408/577,...` tends to `1/sqrt(2)`.
Choose each next packet scale to control its carrier-dependent error
and width/distance ratio. This gives an admissible joint sequence at
the prescribed theta approaching c_Sigma. A fixed nonzero carrier is
not asserted to attain the limit. The proof establishes the speed of
these native information readouts, without selecting physical energy.

### R38.7 — Complete curvature observation preserves the cone and exact support remains faster

R36's complete pointwise reading `O u=(Pu,-PJu)` preserves every address
detector fraction, so all bounds above hold for that observer. They also
persist under R35.6's source-law record multiplicity and simultaneous
native frame transport.

The sharp speed is `1/sqrt(2)` component counts per eight-event block,
one in R36 dual information-length units per block, and `1/(2sqrt(2))`
fine address counts per original event.

**Proof.** R36 gives `||1_D O u||^2=||1_D u||^2` pointwise and hence
for every address set. It preserves initial and final tails and the
fixed-threshold contract. Finite inactive record multiplicity repeats
the same blocks with unchanged operator norm. Pointwise isometric frames
commute with scalar address cuts and weights; simultaneous conjugation
preserves every estimate. A changed edge law, schedule or nonlocal
observer requires its own intertwiner and is not included silently.
The R35 event/address conversions and R36 length form give the stated
count speeds.

The vanishing-threshold limitation has an explicit source witness.
R35's origin preparation has a nonzero corner at component address
(n,n) with relative norm square `256^(-n)`. No earlier block reaches
that address, because each coordinate changes by at most one. A detector
there with threshold `theta_n=256^(-n)` reads its first signal at n,
at speed sqrt(2). Its threshold tends to zero, whereas (38.18) fixes
theta>0. The exact support front has therefore not been replaced by
the reliable-signal cone.

Assigning physical length L_0 per component step and duration T_0 per
block would convert the coefficient to `L_0/(sqrt(2) T_0)`. The present
equations do not select those material observables, identify Z with
light, or impose this cone on independently specified fields. Physical
c, h and alpha require additional derived identifications.

## Certification and source status

The [native application](../04-operator-evolution/native_universal_signal_cone.cjs)
checks finite cyclic resolution, source coefficients, complete block
derivative identities, diagonal speed squares, finite tilts, parameter
budgets, directional covering, front witnesses and curvature detection.
Infinite-scale conclusions use the written proofs; they are not inferred
from a finite phase grid or simulated large-event limit.

The [verification](../04-operator-evolution/R38_VERIFICATION.json) binds
seven written proofs and scoped exact checks, replays frozen R37, and
preserves all 292 earlier non-navigation files. The
[ledger](../04-operator-evolution/R38_DERIVATION_LEDGER.json) separates
native constructions, source consequences and excluded comparisons.
Metadata auditing and finite PASS are not semantic proof-assistant
verification, external review or physical validation.

| Pinned source | Premise status and use |
|---|---|
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** signed/refined counts, H/K, earned iota, roots, pairing and completion. Finite equality, counting and induction remain explicit infrastructure. |
| [R20](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native construction:** coefficient tuples and matching; no supplied physical state space or measurement law. |
| [R35.1–R35.7](https://github.com/Parveen117/extra-ideas/blob/9c13e1752263c756d356dcf5851ac873d2e0d102/02-relational-response/NATIVE_SIGNED_ENVELOPE_R35.md) | **Native-derived within its retained source:** Laurent source, actual bands, complete positive slope polynomial, record/frame transport and faster front. Finite phase resolution and arbitrary-packet bounds are derived anew here. |
| [R36.1, R36.3, R36.7](https://github.com/Parveen117/extra-ideas/blob/011b8224c92d95f7f726eb56bab3e8d8e6f7280c/02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md) | **Native-derived:** dual length and complete pointwise curvature observer. No Riemann adapter is a premise. |
| [R37.1–R37.7](https://github.com/Parveen117/extra-ideas/blob/f06932dd63a4f0827a31583bc42dfaa979ddc663/02-relational-response/NATIVE_PACKET_FLIGHT_R37.md) | **Native-derived:** compact signals, controlled arrival and rational speed family. Supplies sharpness, not the universal upper bound. |
| [F00-E](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md) | **Native analytic theorem on the earned field:** factorial exponential, product law and local turn derivative. R38.3 derives the positive weighting properties it needs. |
| [R10 Riemann bridge](https://github.com/Parveen117/extra-ideas/blob/6035e66bd40fe4016eeded41a4d704ec4ff4b342/02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) | **Comparison only; admitted smooth/tangent/metric adapter.** Excluded from these proof paths. |
| [Canonical engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged native exact arithmetic/replay:** application certificate, not a second engine, whole-engine PASS or formal-assistant verification. |

```bash
python3.12 -B 04-operator-evolution/verify_r38.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```
