# CM2 — compact low levels converge to the core, in every centre sector

9 October 2026. Independent continuation of PC1 at
`4475afd7e548dac9106f23a1f226ee34e3d51646`, using OL1 and CB1.
PC1's 15 checks and 5 tests were replayed before this development.

**The new result is a spectral matching theorem for the actual one-site
compact operator.** For each of its eight fixed centre-flip characters and
each fixed eigenvalue index (k=0,1,\ldots),

\[
\boxed{\lim_{\theta\to\infty}
\frac{E_{k,\sigma}(\theta)}{2\theta^{1/3}}=e_k.}                  \tag{1}
\]

Here (e_k) is the (k)-th eigenvalue, counting multiplicity, of the free
three-turn **physical gauge-singlet** core. Consequently the internal gap
in each fixed centre sector has the same limit:

\[
\boxed{\lim_{\theta\to\infty}
\frac{\Delta_\sigma(\theta)}{\theta^{1/3}}
=2(e_1-e_0)>1.0112.}                                          \tag{2}
\]

The lower coefficient uses PC1's full-space count, with more of its already
certified digits retained, and CB1's frozen correlated ground reading.
CB1 previously had a vacuum-sector liminf lower of 0.6685. Equation (1)
also proves existence of the fixed-sector spectral and gap limits.

There are explicit uniform finite-coupling bounds on all eight internal gaps:

| Coupling range | Bound on each fixed-sector gap, in H units | Minimum at the starting point |
|---|---|---:|
| (\theta\ge2\cdot10^9) | (\Delta_\sigma\ge0.73\theta^{1/3}-878) | (>41.7423) |
| (\theta\ge10^{10}) | (\Delta_\sigma\ge0.82\theta^{1/3}-1002) | (>764.6364) |
| (\theta\ge10^{12}) | (\Delta_\sigma\ge0.93\theta^{1/3}-1476) | (\ge7824) |

These conservative windows complement CB1's earlier vacuum-sector windows;
they do not replace its stronger estimates where those are larger.

**Sector labels matter.** The eight sector ground energies have the same
leading limit. Thus the gap of the operator with all centre sectors included
satisfies

\[
\frac{\Delta_{\rm all}(\theta)}{\theta^{1/3}}\longrightarrow0.    \tag{3}
\]

No splitting rate between those sector bottoms is calculated here. This is
a fixed one-site result, with no spatial-volume or continuum limit.

## 1. Carrier, source, cut and the input from PC1

Keep OL1/CB1's operator and normalized Haar pairing:

\[
H_\theta=-\sum_{i=1}^3\Delta_{S^3,i}
+2\theta\sum_{i<j}|u_i\times u_j|^2,
\qquad U_i=(a_i,u_i)\in S^3.                                  \tag{4}
\]

Take its Friedrichs realization on the simultaneous-conjugation invariant
subspace of (L^2((S^3)^3,d\mu)). The commuting centre flips are
(Z_i:U_i\mapsto-U_i). For (\sigma\in\{\pm1\}^3), the sector
\(\mathscr H_\sigma\) consists of functions satisfying
(F(Z_iU)=\sigma_iF(U)). Each is a closed invariant subspace. Its ordered
levels are (E_{k,\sigma}); its internal gap is
(\Delta_\sigma=E_{1,\sigma}-E_{0,\sigma}). These are distinct from the gap
across all sectors. The unique positive ground of the full compact operator
belongs to the all-even sector.

The matching carrier is the free core

\[
h=-\tfrac12\Delta_C+X(C),\quad X=e_2(CC^T),\quad C\in\mathbb R^{3\times3},
\qquad \mathscr H^{\rm core}_{\rm phys}=L^2(\mathbb R^9,dC)^{SO(3)}.
\]

The core has compact resolvent and a unique positive gauge-invariant ground
by the confinement and positivity argument in CM1. Its levels on the
physical carrier are (e_k). PC1 gives a stronger input on the entire free
carrier, before any gauge restriction:

\[
\#\{\text{levels of }h\text{ strictly below }\rho\}\le1,
\qquad\rho=5.692389.                                         \tag{5}
\]

This is a downward rounding of PC1's exact cluster expression
(5.692389744906\ldots), not a new improved cluster theorem. CM2 freezes
PC1's 64-coefficient two-turn reading, replays its full Gaussian action and
source variance, and verifies the product overlap and cluster line without
rerunning a floating search. The required radial sign certificates are also
replayed. In particular (e_1\ge\rho).

CB1's unchanged 67-coefficient correlated reading gives

\[
e_0\le\eta<5.186743366131,
\qquad 2(\rho-\eta)>1.011291267738>1.0112.                      \tag{6}
\]

The new retained cut is geometric: an inner neighbourhood of the central
configurations, plus its full exterior. The flat core controls the inner
count; the transverse layer controls the exterior. Both norms and the entire
localization error are retained. The explicit intertwiner and its metric
correction are derived next.

## 2. The centre sheets fold to one chart per fixed character

Put (t=\theta^{1/3}), (\epsilon=t^{-1/2}), and
(r(U)^2=\sum_i|u_i|^2). In (r<\epsilon A<1), none of the (a_i) vanish.
There are eight disjoint patches labelled by
(s_i=\operatorname{sign}(a_i)\). Use the **folded coordinates**

\[
c_i=\epsilon^{-1}s_i u_i,
\qquad U_i=s_i\bigl(\sqrt{1-\epsilon^2|c_i|^2},\epsilon c_i\bigr).\tag{7}
\]

A centre flip changes both (s_i) and (u_i), so (c_i) stays fixed.
The character is (\chi_\sigma(s)=\prod_{i:s_i=-1}\sigma_i).
An inner function in the fixed sector is determined by a single function
(f(C)), with value (\chi_\sigma(s)f(C)) in patch (s). Thus there is one
inner copy of the core in each sector, not eight copies within one sector.
Gauge rotations act by (C\mapsto GC), so this folding exactly preserves
the declared gauge condition. It imposes no column parity restriction on
(f).

Up to the common normalization (8\epsilon^9/(2\pi^2)^3), which is removed
by a scalar unitary factor, the local measure and inverse metric are

\[
J_\epsilon(C)=\prod_i(1-\epsilon^2|c_i|^2)^{-1/2},
\qquad B_i=I-\epsilon^2c_ic_i^T.                               \tag{8}
\]

The character phases have modulus one and are constant on each patch.
Functions used below vanish before any patch boundary. Hence they introduce
no uncounted equator derivative or boundary term.

## 3. Exact half-density flattening, including the scalar term

For an inner function define (\phi=J_\epsilon^{1/2}f). This is the exact
unitary map from the weighted local norm to flat Lebesgue norm. Let
(q=\epsilon^2=1/t). For each column,

\[
b_i=\tfrac12\nabla_i\log J_\epsilon
=\frac{q c_i}{2(1-q|c_i|^2)},\quad
B_ib_i=\frac q2c_i,\quad
\operatorname{div}(B_ib_i)=\frac{3q}{2}.
\]

For compactly supported functions, expansion and integration by parts give

\[
\frac{\langle f,H_\theta f\rangle}{2t}
=\int\left\{\frac12\sum_i\nabla_i\bar\phi\cdot B_i\nabla_i\phi
+[X(C)+Q_q(C)]|\phi|^2\right\}\,dC,                           \tag{9}
\]

with the common scalar normalization understood, and

\[
\boxed{Q_q(C)=\frac{9q}{4}
+\sum_i\frac{q^2|c_i|^2}{8(1-q|c_i|^2)}\ge\frac{9q}{4}.}       \tag{10}
\]

Specifically the scalar term per column is
(\frac12\operatorname{div}(B_ib_i)+\frac12b_i^TB_ib_i).
Its sign is positive. The radial reduction of the inverse metric in (8)
remains in the kinetic term; replacing that metric by the identity would
give an invalid lower bound.

On the ball (|C|<A), provided (A^2<t), one has
(B_i\ge(1-A^2/t)I). Since (X\ge0), the Dirichlet inner form therefore
satisfies

\[
\frac{H_{\rm in}}{2t}
\ge(1-A^2/t)h_{B_A,D}+\frac{9}{4t}I.                          \tag{11}
\]

Extending the flat Dirichlet form by zero and applying min–max gives
(e_k(h_{B_A,D})\ge e_k(h)), also on the gauge-invariant subspace.
Thus (11) controls every inner count, not only a chosen trial.

## 4. The entire exterior, including commuting valleys

OL1's full-space comparison gives

\[
H_\theta\ge\sum_i\left[-\tfrac12\Delta_i
+\sqrt{4+4\theta|u_i|^2}-2\right]
\ge2\sqrt\theta\,r(U)-6.                                    \tag{12}
\]

The last step discards only a nonnegative kinetic form and uses
(\sum_i|u_i|\ge r(U)). On functions supported in (r\ge\epsilon R),

\[
H_{\rm out}\ge(2Rt-6)I.                                     \tag{13}
\]

This remains valid when the plaquette potential itself vanishes on a
commuting configuration. It is a quantum transverse-layer floor, not a
pointwise positive lower for the quartic potential.

For (0<R<A), choose (\chi=\cos\alpha(r)), (\xi=\sin\alpha(r)), where
(\alpha=0) up to (\epsilon R), increases linearly to (\pi/2) at
(\epsilon A), and is constant thereafter. These Lipschitz multipliers
preserve the form domain, all centre characters and gauge invariance.
They satisfy (\chi^2+\xi^2=1), and

\[
|\nabla_{(S^3)^3}r|^2
=1-\frac{\sum_i|u_i|^4}{r^2}\le1
\quad(r>0).
\]

The IMS form identity consequently gives

\[
\langle F,H_\theta F\rangle
=\langle\chi F,H_\theta\chi F\rangle
+\langle\xi F,H_\theta\xi F\rangle
-\int(|\nabla\chi|^2+|\nabla\xi|^2)|F|^2d\mu,
\]

\[
|\nabla\chi|^2+|\nabla\xi|^2
\le\frac{\pi^2t}{4(A-R)^2}.                                 \tag{14}
\]

The two-part map (F\mapsto(\chi F,\xi F)) is an isometry into the direct
sum of inner and exterior form carriers. Min–max, (11), and (13) prove,
for each fixed (k\ge0) and each character,

\[
\boxed{\frac{E_{k,\sigma}(\theta)}{2t}
\ge\min\left\{(1-A^2/t)e_k+\frac9{4t},\ R-\frac3t\right\}
-\frac{\pi^2}{8(A-R)^2}.}                                    \tag{15}
\]

Equivalently, below the minimum threshold the exterior contributes no
direction, and the inner form contributes at most (k) directions. This
also explains why the count is made separately for each character.

IMS is a classical localization identity: see Simon, *Semiclassical analysis
of low lying eigenvalues I*, [Lemma 3.1, page 301](https://www.numdam.org/item/AIHPA_1983__38_3_295_0.pdf).
Its proof is the product-rule cancellation displayed above. That paper's
nondegenerate-minimum theorem is not an input for this quartic core with
commuting valleys; (12) supplies the exterior control needed here.

## 5. Fixed-index spectral convergence

For the lower limit in (1), fix (R>e_k) and (A>R). First send
(t\to\infty) in (15), obtaining

\[
\liminf E_{k,\sigma}/(2t)\ge e_k-\frac{\pi^2}{8(A-R)^2}.
\]

Then let (A-R\to\infty). Each prior limit has fixed (A), so the chart
condition is eventually satisfied. This proves the lower limit without an
assumption about decay of the unknown eigenfunctions.

For the upper limit, approximate the first (k+1) flat physical
eigenfunctions in the form norm by compactly supported smooth physical
functions in one fixed ball. Such approximations follow from the usual
form core, cutoff approximation and averaging over the compact gauge group.
Choose an independent ((k+1))-dimensional subspace with maximal flat
Rayleigh value at most (e_k+\delta).

Transport that entire subspace by (7) and the inverse half-density map,
with the chosen character, and extend it by zero. On the fixed support,
(B_i\to I) and (Q_q\to0) uniformly. Norms are exact, and the finite
matrix of forms in (9) converges to the flat one. The maximal Rayleigh value
therefore tends to at most (e_k+\delta). Min–max and then
(\delta\downarrow0) give the upper limit. Together the bounds prove (1).

This is eigenvalue convergence at each fixed index. Uniformity over an
unbounded spectral index, norm-resolvent convergence, or convergence of a
spatial field is not claimed.

Subtracting the (k=0) and (k=1) limits proves the equality in (2);
equation (6) gives its numerical lower. The proof works for all eight fixed
characters because only constant patch signs differ in their local maps.

On the direct sum of the sectors, each limiting core level occurs eight
times its core multiplicity. In particular the seven other sector bottoms
and the all-even vacuum share the same leading scaled value (e_0).
The full compact operator has a unique positive ground at finite coupling;
using any other sector bottom as an excited upper and taking the limit
proves (3). No tunnelling exponent or absolute splitting rate follows from
this leading-order statement.

## 6. Explicit weak windows with all errors retained

For a numerical excited floor use (e_1\ge\rho) in (15). Set (b=A-R).
Our rational parameters satisfy (R\ge\rho) and (A^2\rho\ge21/4), so

\[
\left(R-\frac3t\right)
-\left[(1-A^2/t)\rho+\frac9{4t}\right]
=(R-\rho)+\frac{A^2\rho-21/4}{t}\ge0.
\]

The resulting explicit bound is

\[
E_{1,\sigma}(\theta)
\ge2t\left(\rho-\frac{\pi^2}{8b^2}\right)-2\rho A^2+\frac92.   \tag{16}
\]

For the ground upper in each character, use CB1's cutoff correlated trial
on each patch, multiplied by (\chi_\sigma(s)) in folded coordinates.
It vanishes before each equator. Its norm and form are independent of the
character, so CB1's entire Haar, metric, cutoff and tail estimate transfers:

\[
E_{0,\sigma}(\theta)\le2t\,G(\theta),                         \tag{17}
\]

where (G) is precisely CB1 equation (10). For the finite windows, the trial
cutoff is 1 through radius 6 and falls to zero over width (1/4); its tail
uses moment order 28 and Young parameter (1/500000). The tail mass bound
is less than (1.799285\cdot10^{-12}). These choices refer to the upper
trial, independently of the IMS partition radii.

Combining (16) and (17) gives

\[
\Delta_\sigma\ge
2t\left[\rho-G(\theta)-\frac{\pi^2}{8b^2}\right]
-2\rho(R+b)^2+\frac92.                                      \tag{18}
\]

All roots are rounded outwards. For the IMS coefficient the verifier uses
(\pi<22/7), justified by

\[
\frac{22}{7}-\pi=\int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx>0.
\]

It checks the polynomial division and rational integral part exactly.
For fixed cutoff parameters, (G(\theta)) decreases as (\theta) grows;
the chart condition (A^2<t) also improves. Therefore each row below
certifies its entire interval to infinity:

| θ start | IMS R | IMS b | A | Available gap coefficient lower | Subtracted constant upper |
|---:|---:|---:|---:|---:|---:|
| (2\cdot10^9) | 5.7 | 3.1 | 8.8 | 0.73719 | 877.1373 |
| (10^{10}) | 5.7 | 3.7 | 9.4 | 0.82107 | 1001.4590 |
| (10^{12}) | 5.7 | 5.7 | 11.4 | 0.93319 | 1475.0658 |

Rounding the coefficients down and constants up gives the stated three
uniform gap bounds. The record retains the exact rational computations.
Restoring YM-line units uses the unchanged map
(\theta=4\theta_{\rm YM}), (\Delta(A)=\Delta(H)/4).

## 7. Proof status, continuity and replay

| Item | Status |
|---|---|
| PC1 all-space core count | Replayed from a frozen actual pair trial and exact radial signs |
| Compact pairing and source | Exact folded chart and half-density form identity, including (Q_q) |
| Exterior commuting configurations | Included through OL1's full transverse-layer form floor |
| Localization defect | Entire IMS cost retained; no omitted transition region |
| Every fixed low compact level, in each fixed centre sector | Limit (1), by a written two-sided min–max proof |
| Eight internal sector gaps | Common limit (2) and three explicit finite-coupling windows |
| Gap across different centre sectors | Leading normalized limit zero; absolute splitting rate open |
| Spatial volume, nonconstant modes, 4D continuum | Not addressed by this one-site theorem |

The analytic inputs are the Friedrichs form/min–max principle, compact
resolvent and form-core density, positivity on the compact carrier, the
standard coordinate formula for a sphere, and the IMS product identity.
They are named and used in the written proof. Finite checks are arithmetic
evidence and controls, not substitutes for those analytic interfaces or
proof-assistant verification.

All PC1, CB1 and earlier files remain unchanged. DR2 remains a separate
free-core dimension result at `188b2606918404ac9a69137f10dade9bcfd7460e`;
CM2 does not require its all-d theorem. Here the number of core columns
stays three. The physical rate still needs the torus/lattice energy unit,
nonconstant-mode matching and the desired spatial-volume/cutoff limit.

The meaningful next physical target is now beyond this fixed-site spectral
adapter: local interactions between cells, retaining the full source and
centre-sector structure while controlling a reserve independent of volume.

```sh
python physics/cm2/cm2_compact_core_matching.py --check
python -m unittest discover -s physics/cm2 -p 'test_*.py'
```

`CM2_RESULT.json` pins the proof, code, tests, frozen pair reading and all
imported numerical/analytic stage inputs by SHA-256. Verification uses Python
3.12 standard library and exact fractions; 20 exact/outward checks and 8 tests
pass. The one-time PC1 floating proposal
is recorded explicitly in `CM2_PAIR_TRIAL.json`; replay uses only its frozen
rational coefficients. Tests independently differentiate the original
spherical operator through the half-density map and check the sector count
and IMS controls.
