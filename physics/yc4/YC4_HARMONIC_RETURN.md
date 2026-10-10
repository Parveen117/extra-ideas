# YC4 — larger harmonic cuts certify compact centre splitting through theta = 8

9 October 2026. Continuation of YC3 at `5904d33` on the integrated
CM2/DR2 baseline. A read-only remote refresh found no newer development.
Status: written analytic proof plus exact rational replay; not formal
proof-assistant verification or an independent expert review.

## 1. Result and what changes

On the **actual one-site SU(2)^3 gauge-singlet carrier**, including all eight
centre characters, YC4 extends YC3's identified-first-excitation interval
from [0,2] to **[0,8]**. Throughout that interval:

\[
\boxed{\Delta_{\rm all}(\theta)\ge\frac15.}
\tag{1}
\]

The first excitation consists of precisely the three equivalent
single-centre-odd sector bottoms; each is simple. The all-even ground is
simple too. The following **pointwise** gap windows are separately certified:

| theta | Lower bound | Upper bound |
|---:|---:|---:|
| 2 | 2.1542 | 2.1655 |
| 4 | 1.5272 | 1.6747 |
| 6 | 0.9540 | 1.5538 |
| 8 | 0.3063 | 1.8460 |

Decimals in this table are exact terminating rationals, rounded outwards.
The uniform 1/5 in (1) is smaller than the point bound at 8 because interval
coverage uses opposite ends of each coupling cell. No monotonicity of the
gap itself is assumed. YC3's expansion and its cubic remainder are preserved;
YC4 does not extrapolate that power series past its proved interval.

CB1 already proved positive one-site all-sector gaps at every finite coupling.
The development here is a useful numerical margin, actual excited-sector
identification, and a complete enlarged-cut certification method, not first
positivity. This is still not a spatial-volume-uniform or continuum theorem.

## 2. Exact carrier and normalization

Write U_i = a_i I + i u_i dot sigma, with a_i^2+|u_i|^2=1. Use normalized
product Haar measure and the simultaneous-conjugation invariant subspace
of L^2(SU(2)^3). In the unchanged OL1/YC3 normalization,

\[
H_\theta=H_0+\theta V,\quad
H_0=-\sum_{i=1}^3\Delta_{S^3,i},\quad
V=2\sum_{i<j}|u_i\times u_j|^2,\qquad0\le V\le6.
\tag{2}
\]

Let k denote the number of centre-odd links, chosen to be the first k.
There are binomial(3,k) equivalent centre characters. E_k denotes the bottom
of one such sector; E_(j,k) denotes its j-th eigenvalue, indexed from zero.
This distinguishes a centre label k from an excitation index j.

We retain **every gauge singlet in every free harmonic shell** with energy
strictly below Lambda_k = 3k+24. Let P_k be their orthogonal projector and
Q_k = 1-P_k. Compared with YC3:

| k | YC3 retained rank | YC4 retained rank | YC3 hidden free floor | YC4 hidden free floor |
|---:|---:|---:|---:|---:|
| 0 | 1 | 13 | 8 | 24 |
| 1 | 1 | 21 | 11 | 27 |
| 2 | 2 | 32 | 14 | 30 |
| 3 | 5 | 23 | 21 | 33 |

The new floor is not guessed from the finite matrix. Section 3 establishes
that the entire free spectrum below it has been retained. Since V >= 0,

\[
D_\theta:=Q_k H_\theta Q_k\ge\Lambda_k Q_k.
\tag{3}
\]

These compressions are defined by their restricted closed quadratic forms.
All retained functions are smooth, so the finite-to-hidden off-diagonal map
is bounded and the block elimination used below is well defined.

## 3. Complete harmonic shells, not selected trial polynomials

A one-link degree-n harmonic has free energy n(n+2), centre parity (-1)^n,
and under simultaneous conjugation contains each ordinary SO(3) spin
ell=0,...,n exactly once. One way to see the last assertion is the spherical
separation into Gegenbauer radial factors and spin-ell spherical harmonics;
their dimensions sum to (n+1)^2. Equivalently, the SU(2) left/right
(n/2,n/2) representation restricts diagonally to spins 0 through n.
The standard branching formula is explicitly recorded in Peter Kramer,
[arXiv:1504.01096, Section 4, equations (20)–(24)](https://arxiv.org/pdf/1504.01096).
No cosmological or topological claims from that paper are used.

Consequently the number of singlets in a shell n=(n_1,n_2,n_3) is

\[
m(n)=\sum_{0\le\ell_i\le n_i}
{\bf1}\{|\ell_1-\ell_2|\le\ell_3\le\ell_1+\ell_2\}.
\tag{4}
\]

The spin triangle has no extra even-sum restriction: the epsilon tensor
allows the odd case. Excluding it would wrongly discard SO(3) singlets.
The certificate enumerates all parity-allowed n with sum_i n_i(n_i+2)
< Lambda_k and constructs precisely m(n) independent singlets for each.

Construction uses contractions a_i, u_i dot u_j and
tau=u_1 dot (u_2 cross u_3), with tau exponent 0 or 1. On a polynomial
homogeneous of link degree n_i, the harmonic projection on that link is

\[
\Pi_{n_i}=\prod_{\substack{0\le r<n_i\\r\equiv n_i\pmod2}}
\frac{h_i-r(r+2)}{n_i(n_i+2)-r(r+2)},\qquad h_i=-\Delta_{S^3,i}.
\tag{5}
\]

The lower harmonics are exactly the other possible components of a restricted
homogeneous polynomial, so (5) is a genuine finite spectral projection.
Haar Gram–Schmidt, exact eigen-equations, gauge invariance and rank (4)
verify that the generated vectors fill the entire shell. Thus completeness
does not rest on assuming that a chosen monomial list spans it. Every omitted
harmonic has energy >= Lambda_k; this proves (3) on the infinite complement.

The scalar-coordinate reflections A_i:(a_i,u_i)->(-a_i,u_i) commute with
H_theta. They split the retained space into blocks with rational orthogonal
bases. All blocks are counted; no rotational or reflection excitation is
removed from the physical target. Centre flips and A_i are different maps.

## 4. The complete lost-source Gram and the matrix return bound

Let f_i be the retained orthogonal basis, G_ij=<f_i,f_j>,
K_ij=<f_i,H_0 f_j>, M_ij=<f_i,V f_j> and W_ij=<V f_i,V f_j>.
The norms in G need not be one. The full source Gram is

\[
\boxed{B=W-MG^{-1}M,\qquad B_{ij}=\langle Q_kVf_i,Q_kVf_j\rangle.}
\tag{6}
\]

This is an exact integral over the whole compact carrier. There is no
Galerkin cutoff inside W; all parts of V f_i are integrated. In particular,
replacing W by M G^(-1) M would set the genuine lost source to zero.
Different A-parities are orthogonal and V preserves them, so (6) is computed
block by block without omitting a cross-block return.

For a real threshold z < Lambda_k, define

\[
T_k(\theta,z)=K+\theta M-zG,
\qquad
L_k(\theta,z)=T_k(\theta,z)-\frac{\theta^2}{\Lambda_k-z}B.
\tag{7}
\]

Eliminating the entire Q space in H_theta-z gives the exact Schur form

\[
S_k=T_k-\theta^2P_kVQ_k(D_\theta-z)^{-1}Q_kVP_k.
\]

Equation (3) bounds the resolvent by (Lambda_k-z)^(-1), hence in the retained
Haar metric L_k <= S_k <= T_k. Completion of the block square gives
N_-(H_theta-z)=N_-(S_k) because D_theta-z is positive. Therefore

\[
\boxed{N_-(T_k)\le N_-(H_\theta-z)\le N_-(L_k).}
\tag{8}
\]

Here N_- counts strict negative directions; singular finite forms are handled
exactly. More generally, if T_k has at least j+1 nonpositive directions,
min–max gives E_(j,k) <= z. If L_k has at most j negative directions,
E_(j,k) >= z. Thus (8) supplies **lower** spectral bounds and eigenvalue
counts, unlike Ritz upper values alone. Congruence with G^(-1/2) is implicit
only for interpretation; the rational verifier works directly in G.

The underlying Schur/min–max mechanism is standard. The application-specific
content is the complete compact harmonic cut, its exact full source, the
certified floor and the explicit coupling coverage. No abstract-method
priority claim is made.

## 5. Exact endpoint signs and proof for every coupling in [0,8]

`YC4_ENDPOINTS.json` records rational lower and upper ground bounds for
k=0,1,2,3 at every theta_n=n/8, n=0,...,64. It also records lower bounds on
E_(1,0) and E_(1,1), the second levels in the vacuum and single-odd sectors.
Every endpoint is independently checked by exact inertia in (7)–(8).
The proposal generator used floating eigenvalues only to choose rational
thresholds; no floating calculation supplies a certificate sign.

Since V >= 0, min–max makes each ordered eigenvalue nondecreasing in theta.
For a cell [a,b] between neighboring mesh points,

\[
E_{(j,k)}(a)\le E_{(j,k)}(\theta)\le E_{(j,k)}(b).
\tag{9}
\]

Denote recorded ground bounds by l_k and u_k and second-level lower bounds
by s_0 and s_1. The exact replay checks on every one of the 64 closed cells:

\[
\begin{aligned}
l_1(a)-u_0(b)&\ge1/5,\\
\min\{l_2(a),l_3(a)\}-u_1(b)&>0,\\
\min\{s_0(a),s_1(a)\}-u_1(b)&>0.
\end{aligned}
\tag{10}
\]

The first line separates the vacuum from the lowest single-odd level. The
second excludes the other centre bottoms. The third excludes all vacuum
excitations and all higher single-odd levels, and proves the simplicity of
the relevant eigenvalues. Link permutations supply the three equivalent
single-odd copies and no others. Equations (9)–(10) prove (1) and the
multiplicity statement everywhere, not just at the mesh points. There is no
assumption that a nonlinear root or a gap varies monotonically.

The source packet includes every cell margin and every source hash, so a
missing mesh interval or stale threshold is rejected rather than silently
interpolated. Point gap windows use l_1(theta)-u_0(theta) and the reverse
subtraction u_1(theta)-l_0(theta).

The minimum cell gap reserve is exactly 2123/10000, above the displayed
1/5. The minimum other-centre ordering margin is 3339/5000, and the
minimum second-level separation is 29271/10000. The same cells give stronger
uniform gap floors 493/250 on [0,2], 861/625 on [0,4], and 8331/10000
on [0,6]. These are lower bounds, not fitted energies or optimal constants.

## 6. Native-framework interpretation and limits

| Native role | YC4 realization |
|---|---|
| Retained cut | Complete free harmonic shells below a declared threshold, including all gauge-singlet branches |
| Seen metric | Exact product-Haar Gram G, not an assumed identity |
| Lost source | QVP and its exact Gram B=W-MG^(-1)M |
| Hidden floor | QH_theta Q >= Lambda_k on the entire complement |
| Return | Positive energy-dependent resolvent correction in (7)–(8) |
| Target | Actual quantum sector bottoms, their ordering and the all-sector gap |

The observer is refined by bringing whole previously hidden harmonic shells
into its reading. The rest remains represented by its full source Gram and
proved floor. This is the concrete improvement over YC3's smaller cut.
For k=0 and k=1, the new retained space contains the entire first residual
resolved in YC3; higher returns remain and are still bounded, not declared zero.
For k=2 and k=3, a residual component at increment 24 lies at the new boundary
and is correctly left in B. A refinement does not license dropping the remainder.

TVSP -> VTSP remains a typed chart exchange with the metric/source transformed
alongside it. No response ratio is equated to an energy merely because both
are dimensionless. YC3 already related YC2's static variance to one source norm;
YC4 now enlarges the actual quantum cut and carries its kinetic operator too.

The old-to-new map is inclusion of nested finite spectral subspaces of the
**same** Hilbert carrier. It is not a transfer from one-site volume to a larger
lattice, nor from column count in DR2 to spatial dimension or volume.
CM2's common leading weak-limit sector energy stays intact. YC4 proves no
absolute weak-limit splitting rate or tunnelling law. In YM-line units,
theta=4 theta_YM and Delta(A)=Delta(H)/4: the interval becomes [0,2] and
the uniform A-gap bound is 1/20. This is not a dimensional mass calibration.

Exploration at theta=10 finds that this particular coarse return bound no
longer separates the vacuum and single-odd bounds. That is a failure of this
certificate's margin, **not** evidence for gap closing or a level crossing.
CB1's all-coupling finite-site positivity remains available. The current stage
claims no larger interval, thermodynamic limit, continuum field or Clay result.

## 7. Replay and next bounded development

From the repository root, using Python 3.12 and only the standard library:

```sh
PYTHONDONTWRITEBYTECODE=1 python physics/yc4/yc4_harmonic_return.py --check
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s physics/yc4 -p 'test_*.py' -v
```

The optional `explore_yc4.py --write-proposals` requires NumPy/SciPy and
regenerates candidate endpoint data; it also checks their exact signs before
writing. It is not called by replay. Tests cross-check Haar pairing, harmonic
rank, generator action, full-source projection, congruence inertia and negative
controls where a dropped source or omitted free state would make a false bound.

The next concrete target is a source-resolved hidden solve, retaining selected
Q dynamics instead of replacing every hidden denominator by Lambda_k-z.
Alternatively a further complete harmonic shell can raise the floor again.
Both are testable ways to extend the intermediate-coupling margin; neither
alone supplies the large-volume/4D limit.
