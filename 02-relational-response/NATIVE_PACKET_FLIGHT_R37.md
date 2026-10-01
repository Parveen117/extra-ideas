# R37: compact native packets, reliable arrival and a derived flight meter

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R35 derived a sharp phase-slope coefficient and a modewise limit. R36
constructed a complete curvature observer and a leading information
current. R37 supplies actual **finite-support source packets** whose
continued matching information stays close to a rigidly translated copy.
The error is controlled for every finite event count by native finite
products, differences and count sums.

An explicit rational carrier gives packet drift 2/3 component steps per
eight-event block. A native family gives 2/3, 12/17, 408/577, ... and
approaches 1/sqrt(2). A detector retaining a fixed positive fraction of
the packet's matching norm has a shrinking relative arrival window.
This turns the earlier phase coefficient into a constructive packet-flight
and observation result.

The result is a family of reliable signals with explicit errors. It is
not a universal ballistic bound for all initial fields, a zero-threshold
causal-front theorem or a physical identification of light. R35's faster
nonzero support front remains. No Fourier transform, stationary-phase
theorem, continuum wave equation, stochastic law or external clock is
assumed in the proofs below.

## Source and finite target

Use R35's sixteen-role Laurent transport
\(Z=\sum_{d\in\{-1,0,1\}^2}Z_dS_x^{d_x}S_y^{d_y}\),
constructed from the retained R32 H/K source. Scalar component shifts act
by \((S_if)(x)=f(x-e_i)\); \(\delta_i=S_i-I\). Work on finite
fields on the unwrapped component-count plane, with the source-derived
matching norm. Its completion is used only where already earned.
One component step is two coarse or four fine address increments; one
Z block contains eight original source events.

R35 proves

\[
Z^\dagger Z=I,\quad Z+Z^\dagger=2I-L_*,\quad
L_*=\tfrac12(\delta_x^\dagger\delta_x+
\delta_y^\dagger\delta_y)
-\delta_x^\dagger\delta_x\delta_y^\dagger\delta_y/16.
\tag{37.0}
\]

All cut scalars, unit turns, positive roots, finite sums and matrix
coefficients below are native R16 constructions. The phase carrier is a
specified source probe; it is not a fitted physical parameter. The
comparison averaging operator used below is an algebraic estimate for
the unchanged signed source, not a postulated random event law.

## Written results

### R37.1 — Rational native carriers realize two actual bands and a derived drift

Choose a positive native rational a with \(a^2>2\), and define

\[
t=\frac{a^2-2}{2a},\quad w=\frac{a^2+2}{2a},\quad
g=\frac1w,\quad D=1+t^2,
\]
\[
\boxed{z=\frac{1-t^2+2\iota t}{D},\qquad
\zeta=\frac{1+\iota tw}{D}.}
\tag{37.1}
\]

Both phases have native unit norm, since \(w^2-t^2=2\). Evaluate
the source at (z,1), writing \(Z_0=Z(z,1)\). Its actual branch cuts are

\[
\Pi_+=\frac{Z_0-\bar\zeta I}{\zeta-\bar\zeta},\qquad
\Pi_-=I-\Pi_+.
\tag{37.2}
\]

They are self-dagger disjoint rank-eight cuts with
\(Z_0\Pi_+=\zeta\Pi_+\),
\(Z_0\Pi_-=\bar\zeta\Pi_-\).
For the normalized carrier transport

\[
W(S_x,S_y)=\zeta^{-1}Z(zS_x,S_y),\quad
W_0=W(1,1),\quad
U_i=\sum_d d_i\zeta^{-1}Z_dz^{d_x},
\]

one has the exact first-jet identities

\[
\boxed{\Pi_+U_x\Pi_+=g\Pi_+,\qquad
       \Pi_+U_y\Pi_+=0.}
\tag{37.3}
\]

Thus the drift (g,0) is derived from the source branch. In particular,

| a | Spatial phase z | Event phase zeta | Derived x drift |
|---|---|---|---|
| 2 | (3+4 iota)/5 | (4+3 iota)/5 | 2/3 |
| 3/2 | (143+24 iota)/145 | (144+17 iota)/145 | 12/17 |
| 17/12 | (166463+816 iota)/166465 | (166464+577 iota)/166465 | 408/577 |

**Proof.** The definitions give w-squared=t-squared+2 by subtraction
of two finite squares. The squared numerators of z and zeta both equal
D-squared, proving their unit norm. In (37.0), z gives
\(\lambda_x=4t^2/D\), \(\lambda_y=0\), so
\(\ell=2t^2/D\). The two roots of the source quadratic are exactly
zeta and its dagger, by R35.3. They are distinct because t>0. Its finite
quadratic factorization gives (37.2). Daggering that formula and using
\(Z_0^\dagger=(\zeta+\bar\zeta)I-Z_0\) proves self-daggerness;
the source trace from R35.3 gives the two ranks.

For completeness, perturb the phases locally by the earned Euler turns
\(z\operatorname{Exp}_\Sigma(\iota q_x)\) and
\(\operatorname{Exp}_\Sigma(\iota q_y)\). Differentiating the finite
branch equation and multiplying by Pi_+ on both sides cancels the
derivatives of the branch cut. It gives
\(\Pi_+\partial_i Z\Pi_+=(\partial_i\zeta)\Pi_+\).
The source scalar quadratic, or R35.4's derived local turn differential,
gives at the carrier

\[
\frac{\partial_x\zeta}{\iota\zeta}
 =\frac{2t/D}{2tw/D}=1/w,
\qquad \partial_y\zeta=0.
\]

Since \(\partial_i Z=\iota\zeta U_i\), this proves (37.3).
Only a local native phase perturbation is used; no global argument angle
or plane-wave expansion of arbitrary fields is assumed. All displayed
carrier and jet identities are also checked exactly in the native engine.

### R37.2 — A finite source corrector removes the entire first-difference error

Set \(\eta=\bar\zeta/\zeta\), and construct

\[
C_i=\frac{\Pi_-U_i\Pi_+}{1-\eta},\qquad
\mathcal T=\Pi_++C_x\delta_x+C_y\delta_y,
\qquad E_g=I+g\delta_x.
\tag{37.4}
\]

The corrected packet residual
\(\mathcal R=W\mathcal T-\mathcal T E_g\) has an exact finite
factorization

\[
\boxed{\mathcal R=\mathcal R_{xx}\delta_x^2+
        \mathcal R_{xy}\delta_x\delta_y+
        \mathcal R_{yy}\delta_y^2.}
\tag{37.5}
\]

This is a complete Laurent identity, not an assertion from a sampled
phase or a first-order truncation. Each factor is source-computable.
For a coefficient matrix M let
\(m(M)=\sum_{ij}(|\operatorname{rad}M_{ij}|+
|\operatorname{turn}M_{ij}|)\), and sum m over Laurent coefficients
to define operator mass. Put

\[
\kappa=m(C_x)+m(C_y),\quad C_T=1+2\kappa,\quad
C_{ij}=m(\mathcal R_{ij}),\quad C_\Sigma=C_{xx}+C_{xy}+C_{yy}.
\tag{37.6}
\]

Then \(\|\mathcal T\|\le C_T\), and

\[
\|\mathcal R f\|\le C_{xx}\|\delta_x^2f\|
 +C_{xy}\|\delta_x\delta_y f\|+C_{yy}\|\delta_y^2 f\|.
\tag{37.7}
\]

For a=2 the exact constants from the coefficient construction are

\[
\kappa=214/9,\quad C_T=437/9,\quad
C_{xx}=13759/480,\quad C_{xy}=392/15,\quad C_{yy}=7819/240.
\tag{37.8}
\]

These are conservative bounds; no optimal packet width is claimed.

**Proof.** On the two branch cuts, W_0 is I and eta I. Therefore
\((W_0-I)C_i=-\Pi_-U_i\Pi_+\). Equation (37.3) shows that the
constant and both first jets of R vanish. To prove the stronger full
identity, define for every integer d the finite Laurent polynomial r_d by

\[
s^d=1+d(s-1)+(s-1)^2r_d(s),
\]
\[
r_d(s)=\sum_{j=0}^{d-2}(d-1-j)s^j\quad(d\ge2),\qquad
r_{-p}(s)=\sum_{j=1}^{p}(p-j+1)s^{-j}\quad(p\ge1),
\tag{37.9}
\]

with r_0=r_1=0. Multiplication by (s-1)^2 cancels interior terms and
proves each formula. For every pair of integer exponents a,b,

\[
s^at^b-1-a\delta_x-b\delta_y
=\delta_x^2r_a(s)t^b+ab\delta_x\delta_y
 +\delta_y^2(1+a\delta_x)r_b(t).
\]

Apply this to each coefficient R_ab. The cancelled jets remove the
three leading sums; the remaining three explicit sums define
R_xx,R_xy,R_yy and prove (37.5). The certificate reconstructs every
Laurent coefficient from them. Each matrix unit and each scalar shift
has matching norm one; triangle inequalities therefore bound a finite
operator by its mass. Also delta_i has norm at most two and Pi_+ has
norm one. These facts prove (37.6)–(37.7) and the computed constants.
Without the C_i corrector the first jets generally remain nonzero.

### R37.3 — Compact native count profiles have exact norm and difference bounds

For an integer M>=1 let \(b_M=\sum_{j=0}^{M-1}\delta_j\) on one
count line. Form the three-box convolution h_M=b_M*b_M*b_M, and then
the two-direction scalar profile \(f_M(x,y)=h_M(x)h_M(y)\).
These are finite native count arrays. No Gaussian, smooth bump or
probability distribution is imported.

Their exact one-direction counts are

\[
\boxed{N_M=\|h_M\|^2=(11M^5+5M^3+4M)/20,}
\]
\[
\boxed{\|\delta h_M\|^2=M(M^2+1),\qquad
       \|\delta^2h_M\|^2=6M.}
\tag{37.10}
\]

Consequently

\[
\frac{\|\delta_i f_M\|}{\|f_M\|}\le\frac2M,
\qquad
\frac{\|\delta_i\delta_j f_M\|}{\|f_M\|}\le\frac4{M^2}
\quad(i,j\in\{x,y\}).
\tag{37.11}
\]

For any fixed unit native role v in im(Pi_+), let
\(F_M=f_Mv\) and \(\phi_0=\mathcal T F_M\). If M>2 kappa,

\[
\|\phi_0\|\ge(1-2\kappa/M)\|F_M\|>0.
\tag{37.12}
\]

The corrected packet is supported inside \([0,3M-2]^2\).

**Proof.** Finite convolution counts triples of integers between zero
and M-1. Symmetry turns its squared norm into the central coefficient
of \((1+s+\cdots+s^{M-1})^6\). Multiplying by (1-s)^6 gives
(1-s^M)^6. The coefficients of its formal inverse (1-s)^(-6) count
six nonnegative integers with a given sum: inserting five separators
gives the finite count \(\binom{n+5}{5}\). Equivalently the same
coefficient follows by six repeated finite partial sums. Thus

\[
N_M=\binom{3M+2}{5}-6\binom{2M+2}{5}+15\binom{M+2}{5}.
\]

Expanding the three degree-five products proves the first formula of
(37.10), including M=1,2 where the small binomial counts vanish.

Since delta b_M=delta_M-delta_0, delta h_M is the difference of two
shifted two-box triangles. Their norm is M(2M^2+1)/3 and their overlap
is M(M^2-1)/6. These finite sums follow by telescoping adjacent squares
and cubes; subtracting twice the overlap from twice the norm gives
M(M^2+1). The second difference is three disjoint boxes with coefficients
1,-2,1, so its squared norm is 6M. This proves all of (37.10).

The two-direction norms factor into products of these finite sums.
Use N_M>=11M^5/20 and M(M^2+1)<=2M^3 to obtain the first bound of
(37.11), its mixed counterpart, and the second-difference bound.
Applying (37.4) to F_M, with Pi_+v=v, gives (37.12) by the triangle
inequality and (37.11). The three-box support is [0,3M-3] on each
axis; the correction adds at most one positive step, giving the stated
compact support. Normalizing a nonzero branch role uses the already
derived positive native root, not a new state-space axiom.

### R37.4 — Finite products give an all-event bound on rigid packet flight

Let \((\mathcal A_z f)(x,y)=z^{-x}f(x,y)\). Prepare
\(u_0=\mathcal A_z\mathcal T F_M\) and continue the actual source
\(u_n=Z^n u_0\). For M>2 kappa, n>=0 and
\(m_n=\lfloor ng\rfloor\),

\[
\boxed{\|u_n-\zeta^n z^{-m_n}S_x^{m_n}u_0\|
       \le\epsilon_{n,M}\|u_0\|,}
\tag{37.13}
\]
\[
\boxed{\epsilon_{n,M}
 =\frac{[4C_\Sigma+2C_Tg(1-g)]n/M^2+2C_T/M}
             {1-2\kappa/M}.}
\tag{37.14}
\]

This is an exact finite-time error estimate for source packets, not an
assumed continuum or asymptotic wave-packet law.

**Proof.** Both A_z and multiplication by zeta are pairing isometries.
Source composition gives
\(Z\mathcal A_z=\zeta\mathcal A_z W\), so W is also an
isometry. Because 0<g<1, \(E_g=(1-g)I+gS_x\) has norm at most
one by the matching triangle inequality. Ordered finite telescoping gives

\[
W^n\mathcal T-\mathcal T E_g^n
 =\sum_{j=0}^{n-1}W^{n-1-j}\mathcal R E_g^j.
\]

The scalar E_g commutes with the differences. Equations (37.7) and
(37.11) therefore bound this error on F_M by
\(4nC_\Sigma\|F_M\|/M^2\).

For the remaining comparison, multiply E_g n times. Native finite
binomial counts give nonnegative weights w_r with sum one, mean ng and
second centred sum ng(1-g). These follow directly by updating the three
finite sums under another multiplication by (1-g)+gs; they do not
assume independent random events. Put alpha=ng-m_n in [0,1). Applying
(37.9) to the integers d=r-m_n gives

\[
\|E_g^n f-S_x^{m_n}f\|
\le\alpha\|\delta_x f\|
 +\tfrac12ng(1-g)\|\delta_x^2 f\|.
\tag{37.15}
\]

Indeed r_d has nonnegative coefficient sum d(d-1)/2 for every integer d.
Its weighted sum is
\(\{ng(1-g)+\alpha^2-\alpha\}/2\le ng(1-g)/2\).
This proves (37.15) entirely through finite counts and shift identities.
Apply T, use its norm bound, then use (37.11) and alpha<=1. Adding the
two errors and dividing by (37.12) gives (37.14). Finally scalar shifts
commute with T and
\(\mathcal A_zS_x^m=z^{-m}S_x^m\mathcal A_z\), proving (37.13).
The unit scalar in that formula retains the internal signed phase while
leaving packet support and matching size unchanged.

### R37.5 — Compact information signals have reliable finite arrival windows

Fix a carrier. Choose integers L with L^2>2 kappa and set

\[
M=L^2,\qquad N=L^3,\qquad
\varepsilon_L=\frac{A/L+2C_T/L^2}{1-2\kappa/L^2},
\quad A=4C_\Sigma+2C_Tg(1-g).
\tag{37.16}
\]

The packet has initial width at most 3L^2-1 and travels
\(m_N=\lfloor gL^3\rfloor\) component counts. If

\[
\mathcal D_N=[m_N,m_N+3M-2]\times[0,3M-2],
\]

then its address cut obeys

\[
\boxed{\frac{\|1_{\mathcal D_N}u_N\|^2}{\|u_0\|^2}
       \ge1-\varepsilon_L^2.}
\tag{37.17}
\]

The initial width divided by travel time tends to zero, while the escaped
matching fraction tends to zero. These are spatially localized moving
signals in a source-count sense, beyond a mere plane-phase relation.

There is an operational arrival reading. Let
\(\Delta_N=\lceil(3M+1)/g\rceil\). For N>Delta_N and any
threshold \(\varepsilon_L^2<\theta<1-\varepsilon_L^2\), the
first integer n with
\(\|1_{\mathcal D_N}u_n\|^2/\|u_0\|^2\ge\theta\) satisfies

\[
\boxed{N-\Delta_N<n\le N.}
\tag{37.18}
\]

In particular any fixed threshold 0<theta<1 works for sufficiently large
L. The relative window Delta_N/N tends to zero, and the displacement
m_N divided by the detected arrival block tends to g.

**Proof.** At time N the translated comparison packet in (37.13) is
supported inside D_N. Its complement is a native orthogonal address cut,
so (37.13) bounds the outside norm by epsilon_L times the conserved
source norm. Subtracting its square from total norm proves (37.17).
Equation (37.16) gives epsilon_L=O(1/L), while M/N=1/L, proving the
localization and scale claims.

For any n<=N-Delta_N, the rounded shifts obey
\(m_N-m_n\ge g\Delta_N-1\ge3M\). Thus the translated comparison
packet at that earlier n is disjoint from D_N. Its detected fraction
is at most epsilon_{n,M}^2<=epsilon_L^2 by (37.13)–(37.14).
At N it is at least 1-epsilon_L^2. Finite earliest-event selection
therefore proves (37.18). The explicit window has size O(M), so its
relative size tends to zero. This derives the arrival instrument without
an assumed measurement probability law: theta is a declared matching-norm
readout threshold. No zero-threshold or first-nonzero-amplitude arrival
claim follows from these inequalities.

The certificate evaluates these rational bounds, including half-norm
thresholds for a=2 and L=1024,4096,16384. It does not numerically evolve
those enormous packets or label a symbolic bound a hardware experiment.

### R37.6 — The reliable packet-speed family approaches the native coefficient with controlled limit order

Starting from a_0=2, construct native rational upper cuts

\[
a_{k+1}=\frac{a_k^2+2}{2a_k},\qquad
g_k=1/a_{k+1}.
\tag{37.19}
\]

They give the carriers of R37.1 and reliable packet flights of R37.5.
Every g_k is strictly below 1/sqrt(2), and

\[
\boxed{g_k\longrightarrow1/\sqrt2.}
\tag{37.20}
\]

Thus the sharp R35 phase coefficient is approached by actual compact
packet families observed at a fixed positive information threshold.
This establishes achievability, not a universal upper theorem for every
possible native packet or observer.

**Proof.** Direct subtraction gives

\[
a_{k+1}^2-2=\frac{(a_k^2-2)^2}{4a_k^2},\qquad
\tfrac12-g_k^2=\frac{a_{k+1}^2-2}{2a_{k+1}^2}.
\]

The first positive defect is 2. Each next defect is positive and at most
one quarter of the preceding one, since the cuts remain between the
native positive root and 2. Thus the defects tend to zero with a derived
geometric bound. Native completion supplies that root; the equations
force (37.20). The associated w equals a_{k+1}, exactly as in (37.1).

For each fixed carrier, all constants in (37.6) are finite and R37.5
applies by increasing L. As k increases the band gap closes and the
corrector constants grow; they cannot be discarded from the error bound.
One may choose each next integer L large enough to make epsilon_L and
Delta_N/N smaller than any prescribed native rational tolerance. This
constructs a joint sequence with controlled errors. Substituting the
zero-gap carrier directly in (37.2) or (37.4) is invalid, and no fixed
nonzero carrier is claimed to attain the limiting speed.

### R37.7 — The curvature observer reads the same flight and the event ledger calibrates its counts

Retain R36's complete pointwise reading
\(\mathcal O u=(Pu,-PJu)\), with its exact decoder. For every
finite address detector D and every source event block,

\[
\boxed{\|1_DPu_n\|^2+\|1_DPJu_n\|^2
       =\|1_Du_n\|^2.}
\tag{37.21}
\]

The packet error, matching fractions and reliable arrival windows are
therefore unchanged in the curvature-complete observer. Simultaneously
transported native role frames preserve the statement as in R36.7.
Dropping a readout or changing its norm is a different detector target.

The flight meter uses retained source-event counts: N blocks are 8N
original events; displacement m_N is 2m_N coarse or 4m_N fine address
increments. The asymptotic fine-count speed per original event is g/2.
In the R36 dual information-length convention the speed is sqrt(2)g per
block, approaching one along (37.19). These are explicit conversions of
the same source flight, not independent constants.

**Proof.** R36.3 proves the pointwise matching isometry and inverse.
Address cuts commute with those pointwise role maps, so summing its
identity over D proves (37.21). Applying the same isometry to the
difference in (37.13) preserves the error as well. Frame covariance is
the same adjacent-inverse cancellation already proved in R30/R36.
R35's original event/address construction gives the integer conversions;
R36 gives information length sqrt(2) times component Euclidean length.
Substitution yields the stated speeds.

The threshold and scaling boundaries remain essential. R35's source
has a nonzero corner at component speed sqrt(2), with exponentially
small matching size; its unchanged word proof and finite witnesses are
replayed here. R37 controls retained information near a translated pulse,
not exact support or all thresholds. Nor does a fixed positive threshold
identify this mathematical source as physical light. A material length
L_0 per component step and duration T_0 per block would convert g to
L_0 g/T_0; those physical choices are not fixed by the present equations.
The source event clock, detector and distance readout are now explicit
native constructions. Physical field identification and calibration,
independent-field universality, h and alpha remain separate tasks.

## Certification and source reuse

The [native application](../04-operator-evolution/native_packet_flight.cjs)
reconstructs Z from the unchanged source tags. It checks rational carrier
bands, the complete corrected Laurent remainder, exact compact-profile
counts, binomial moments, finite literal packet evolution, large-flight
rational error bounds, speed cuts, the curvature detector and surviving
fast-front witnesses. The large-flight result is a written all-event
proof with exact bound evaluations; direct finite evolutions are separately
reported and are not claimed to simulate the large examples.

The [verification](../04-operator-evolution/R37_VERIFICATION.json) binds
seven written proofs, nine exact groups, fourteen native word replays,
fifteen false alternatives and nine graph mutation controls. It replays
the frozen R36 chain and preserves all 285 earlier non-navigation files.
The [ledger](../04-operator-evolution/R37_DERIVATION_LEDGER.json) distinguishes
native target constructions, inherited proofs and comparison-only adapters.
Metadata and finite PASS do not constitute a semantic proof assistant,
external review or physical experiment.

| Pinned source | Premise status and use |
|---|---|
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** signed/refined counts, H/K, matching pairing, earned iota and completion. Finite equality, counting and induction remain explicit infrastructure. |
| [R20 tuples](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native construction:** retained finite coefficient tuples and product matching; no physical tensor or measurement axiom is supplied. |
| [R32 retained source](https://github.com/Parveen117/extra-ideas/blob/1105a43209b144bd34a81bfbb99acc7876948e95/02-relational-response/NATIVE_RETAINED_LOOP_GAP_R32.md) | **Native-derived within its interface:** actual source-event continuation and noncommuting record links. Physical selection of the interface remains open. |
| [R35.1–R35.7](https://github.com/Parveen117/extra-ideas/blob/9c13e1752263c756d356dcf5851ac873d2e0d102/02-relational-response/NATIVE_SIGNED_ENVELOPE_R35.md) | **Native-derived:** exact signed envelope, phase bands and slopes, source coefficient ledger, finite error and faster-front boundary. No imported Fourier theorem is used. |
| [R36.1, R36.3–R36.7](https://github.com/Parveen117/extra-ideas/blob/011b8224c92d95f7f726eb56bab3e8d8e6f7280c/02-relational-response/NATIVE_CURVATURE_OBSERVER_R36.md) | **Native-derived:** response-dual length, complete curvature observer, exact observed dynamics, information current and frame transport. R37 now controls compact source flights. |
| [F00-E](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md) | **Native analytic theorem on the earned field:** used only for a local carrier-phase derivative; no classical exponential or wave-packet theorem is imported. |
| [R10 Riemann bridge](https://github.com/Parveen117/extra-ideas/blob/6035e66bd40fe4016eeded41a4d704ec4ff4b342/02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) | **Comparison only; admitted smooth/tangent/metric adapter.** Its existing conditional quotient is preserved and does not enter these packet proofs. |
| [Canonical operator engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged native exact arithmetic and replay:** scoped application certificate, not a second engine, whole-engine PASS or formal-assistant verification. |

```bash
python3.12 -B 04-operator-evolution/verify_r37.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```
