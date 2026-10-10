# YC3 — actual compact centre-sector energies and their hidden return

Status: analytic proof with exact rational replay, not proof-assistant verification.
One compact site, three SU(2) links, gauge singlets, **all eight centre characters**.
This is a finite-coupling certificate; it is not a continuum Yang–Mills gap.

## 1. Reviewed baseline and new statement

Claude's CM2 merge `9e582ae1775d3907915b4ebc7f4e16a52a576154`
contains the substantive CM2 already present here at `98eb5c8`. DR2 merge
`a4c4e2bf0532a5916323e42ed55f3ebac0409634` was followed by context repair
`4aad5a0ec0fbdac727894ad0c315ed822cc0350e`. YC3 uses that repaired head,
merged locally at `8f30b5a`; CM2 and DR2 checks/tests were replayed.
Their frozen evidence is unchanged.

CM2 matches every fixed low compact level to the free physical core in each
centre character as coupling grows. DR2 supplies the free physical core gap
for every integer column count d >= 3; d is not lattice volume. Neither result
computes the absolute difference between the compact centre-sector bottoms.
YC3 addresses that different observable on an explicit finite interval.
CB1 already proves positivity of the all-gauge-sector one-site gap for every
finite nonnegative coupling. YC3 does not claim first positivity; it identifies
the first excited sectors and gives sharper quantitative finite-window data.

Write E_k(theta) for the bottom in one centre sector with k odd links.
Link permutations make sectors with the same k isospectral. In OL1/CB1 units,

\[
H_\theta=H_0+\theta V,\quad H_0=-\sum_{i=1}^3\Delta_{S^3,i},\qquad
V=2\sum_{i<j}|u_i\times u_j|^2,\quad 0\le V\le6.
\]

For every 0 <= theta <= 2, the first excitation in the **whole gauge-singlet
space** consists of the bottoms of the three single-centre-odd sectors. Its
multiplicity is exactly three, and its gap above the all-even vacuum obeys

\[
\boxed{
3-\frac\theta2-\frac{(43/48)\theta^2}{8-7\theta/4}
\ \le\ \Delta_{\rm all}(\theta)\ \le\
3-\frac\theta2+\frac{(19/16)\theta^2}{8-9\theta/4}.}
\tag{1}
\]

In particular, Delta_all >= 65/54 throughout this interval. At theta = 1,
the enclosure is [707/300, 249/92]; at theta = 2 it is [65/54, 47/14].
There is also an exact second-order expansion, with a non-asymptotic remainder:

\[
\boxed{\Delta_{\rm all}=3-\theta/2+(77/1920)\theta^2+R_\Delta,}
\tag{2}
\]
\[
|R_\Delta|\le\theta^3\left[
\frac{43/8}{8(8-7\theta/4)}+
\frac{57/8}{8(8-9\theta/4)}\right],\qquad0\le\theta\le2.
\tag{3}
\]

The coefficient in (2) is exact, not a fit. The remainder can be loose near
the endpoint; intersecting its interval with (1) preserves both certificates.

## 2. Carrier, measure, symmetry and the complete lowest multiplet

Use U_i = a_i I + i u_i dot sigma, with a_i^2 + |u_i|^2 = 1, and normalized
product Haar measure on SU(2)^3. The carrier is
L^2(SU(2)^3)^(Ad SU(2)). Simultaneous conjugation rotates all three vectors
u_i by the same SO(3) rotation and leaves a_i fixed. The three centre operators
C_i reverse all four coordinates of link i. Fix their eigenvalues, with the
first k equal to -1 and the others +1. There are binomial(3,k) equivalent choices.

Degree-n harmonics on unit S^3 have eigenvalue n(n+2) and antipodal parity
(-1)^n. Consequently the lowest unprojected tensor shell in this sector has
energy e_k = 3k: degree one on odd links and degree zero on even links.
A degree-one link decomposes under gauge rotations as scalar plus vector.
The complete gauge-invariant part of this shell therefore has the following
basis (unnormalized):

| k | Complete basis for P_k | Rank | Free floor on Q_k = 1 - P_k |
|---|---|---:|---:|
| 0 | 1 | 1 | 8 |
| 1 | a_1 | 1 | 11 |
| 2 | a_1 a_2, u_1 dot u_2 | 2 | 14 |
| 3 | a_1 a_2 a_3; a_1(u_2 dot u_3), a_2(u_1 dot u_3), a_3(u_1 dot u_2); u_1 dot (u_2 cross u_3) | 5 | 21 |

Completeness is elementary invariant tensor counting: zero vector factors give
one scalar; two vectors can contract once; three vectors can contract once
with the SO(3) epsilon tensor. One vector cannot give a scalar. The triple
product must **not** be discarded: gauge conjugation is SO(3), not O(3).
These linearly independent tensors exhaust ranks 1,1,2,5. The independent
tests also compute the common kernel of the three infinitesimal gauge rotations.

The next parity-allowed degree on an even link raises its energy by 8, and
on an odd link by 12. Higher degrees only increase it. Thus, on the entire
gauge-singlet complement, Q_k H_0 Q_k >= e_k + delta_k, where
delta = (8,8,8,12). This is a spectral completeness argument, not an inference
from a finite list of trial vectors. Taking only a_1 a_2 as P in the k=2
sector would leave u_1 dot u_2 at the **same** free energy in Q and invalidate
that complement floor.

For reference, the elliptic SU(2) normalization n(n+2) agrees with Section 3.2,
equation (3.2), of Baudoin and Bonnefont,
[arXiv:0802.3320](https://arxiv.org/pdf/0802.3320). Only their elliptic
Laplace–Beltrami normalization is used, not their subelliptic estimates.
Here it also follows directly by restricting harmonic homogeneous polynomials
in four variables. The compact spectral theorem and min–max principle are
standard inputs; all new numerical constants below are derived in this packet.

## 3. Exact retained and hidden data

For one link, normalized Haar moments are

\[
\langle q_0^{2m_0}\cdots q_3^{2m_3}\rangle
=\frac{\prod_{j=0}^3(2m_j-1)!!}
 {\prod_{r=0}^{m_0+\cdots+m_3-1}(4+2r)};
\]

an odd power gives zero. Independent links factor. This follows from the
rotationally invariant Gaussian divided by its radius, or the Dirichlet
(1/2,1/2,1/2,1/2) distribution of squared coordinates. No numerical integration
or floating-point diagonalization is used.

In each displayed basis, the Haar Gram matrix, PVP and the residual Gram
matrix are diagonal. For basis vector f_j put
a_j = <f_j,Vf_j>/<f_j,f_j> and r_j = (V-a_j)f_j.
Every r_j is orthogonal to the **whole** P_k. The exact data are:

| k | Diagonal slopes a_j | Squared normalized residuals ||r_j||^2/||f_j||^2 |
|---|---|---|
| 0 | 9/4 | 19/16 |
| 1 | 7/4 | 43/48 |
| 2 | 4/3, 20/9 | 47/72, 715/648 |
| 3 | 1, 5/3, 5/3, 5/3, 10/3 | 11/24, 175/216, 175/216, 175/216, 25/24 |

Let a_k be the smallest slope and b_k^2 = ||Q_k V P_k||^2, the largest
normalized residual square. Thus

\[
(a_k)=(9/4,7/4,4/3,1),\qquad
(b_k^2)=(19/16,43/48,715/648,25/24).
\tag{4}
\]

## 4. Full-space energy enclosure, not just Ritz values

For theta >= 0 and a_k theta < delta_k,

\[
\boxed{\max\left\{e_k,\ e_k+a_k\theta-
\frac{b_k^2\theta^2}{\delta_k-a_k\theta}\right\}
\le E_k(\theta)\le e_k+a_k\theta.}
\tag{5}
\]

Proof: the upper bound is the Rayleigh quotient of the smallest-slope trial.
V >= 0 implies both H_theta >= H_0 and
D_theta := Q_k H_theta Q_k >= e_k + delta_k. Since the upper trial lies
strictly below this floor, the ground eigenfunction has a nonzero P_k component
and the Schur equation is exact:

\[
\left[e_kP_k+\theta P_kVP_k
-\theta^2P_kVQ_k(D_\theta-E_k)^{-1}Q_kVP_k\right]p=E_kp.
\tag{6}
\]

The return operator is positive and bounded above by
theta^2 b_k^2/(delta_k-a_k theta). Taking its Rayleigh quotient proves
(5). The uncomputed space is infinite-dimensional, but its floor is a proven
bound for **all** of it. Nothing is inferred from subtracting two upper bounds.

## 5. Which degenerate branch is the actual sector bottom?

There are further exact symmetries A_i:(a_i,u_i) -> (-a_i,u_i). They are
sphere isometries commuting with H_theta, gauge conjugation and the centre
operators, since V is independent of a_i. They are **not** the centre flips.
The lowest multiplet vectors have distinct triples of A_i parities.
Within each such reflection subspace there is exactly one lowest-shell vector;
its full complement still has floor e_k + delta_k. Subspaces with no lowest
vector start at that floor. Therefore the rank-one version of (5) applies
separately to every branch, with its own slope and residual square.

To show the scalar-product branch f_0 = product_(i<=k) a_i stays lowest, compare
its upper bound e_k+a_0 theta with each other branch's lower bound. The
difference, divided by theta > 0, is at least

\[
(a_j-a_0)-\frac{\theta b_j^2}{\delta_k-\theta a_j}.
\tag{7}
\]

This decreases with theta while the denominator is positive. At theta=2,
the values are 103/384 for k=2, and 449/936 (three times), 373/192 for k=3.
All denominators are positive, so every comparison holds on (0,2]. The
scalar-product branch is simple below the next free eigenvalue in its
reflection subspace. At theta=0, the k=2 and k=3 bottoms remain respectively
twofold and fivefold degenerate; the scalar branch is the right-hand choice.

## 6. Resolve the whole lost source: the exact second-order numbers

For normalized scalar-product f_0, the residual r = QVf_0 has only finitely
many free harmonic energies. Multiplying by V raises the polynomial degree
of at most two links by two. Each even link can acquire free increment 8;
each odd link can acquire increment 12. Thus the allowed residual increments
s are 0, one such increment, or a sum for two different links.

The replay applies the exact polynomial spectral projectors
product_(t != s) (H_0-e_k-t)/(s-t) to r. It verifies reconstruction and the
eigen-equations in exact Haar squared norm, hence as L^2 identities. Its
nonzero squared weights w_s are:

| k | s : w_s | c_k = sum_s w_s/s |
|---|---|---:|
| 0 | 8 : 3/4; 16 : 7/16 | 31/256 |
| 1 | 8 : 25/72; 12 : 1/4; 16 : 7/48; 20 : 11/72 | 311/3840 |
| 2 | 8 : 1/9; 12 : 25/72; 20 : 11/72; 24 : 1/24 | 451/8640 |
| 3 | 12 : 1/3; 24 : 1/8 | 19/576 |

The weight at zero vanishes. Distinct free eigenvalues are orthogonal, so the
sum of weights equals the full scalar residual square in Section 3. There
is no residual tail omitted from c_k.

The rank-one Schur equation in the scalar reflection subspace gives exactly

\[
E_k=e_k+a_k\theta-\theta^2
\langle r,(D_\theta-E_k)^{-1}r\rangle.
\]

Set D_0=QH_0Q in this reflection subspace. Since
0 <= E_k-e_k <= a_k theta <= 6 theta, and 0 <= QVQ <= 6Q,

\[
\|\theta QVQ-(E_k-e_k)Q\|\le6\theta.
\]

The resolvent identity, with inverse norms at most 1/delta_k and
1/(delta_k-a_k theta), now proves

\[
\boxed{E_k=3k+a_k\theta-c_k\theta^2+R_k,\qquad
|R_k|\le\frac{6b_{k,\mathrm{scalar}}^2\theta^3}
{\delta_k(\delta_k-a_k\theta)},\quad0\le\theta\le2.}
\tag{8}
\]

Here b_(k,scalar)^2 is the **first** residual square in the row, not its
maximum. The branch selection in Section 5 is what licenses calling this
an expansion of E_k rather than just an unspecified eigenbranch.

## 7. Identify the first excitation and prove (1)–(3)

The all-even second level is >= 8 by H_theta >= H_0 and min–max. The
single-odd second level is >= 11. Its first level has upper bound
3+7 theta/4 <= 13/2 on [0,2], so it is simple and below the all-even excitation.

For k=2 and k=3, subtract that same upper bound from the *raw* lower expression
in (5). Each difference decreases on [0,2]: its linear coefficient is negative,
and theta^2/(delta-a theta) increases. At theta=2 their energies have lower
bounds 6773/864 and 127/12, both greater than 13/2. Thus they cannot supply
the first excitation anywhere in the interval. Subtracting the vacuum upper
from the k=1 lower gives the left side of (1), positive throughout [0,2].
Hence the vacuum is the unique lowest state and the first excitation consists
of precisely three equivalent single-odd bottoms. At theta=0 the same
statement follows directly from the free spectrum.

The lower side of (1) decreases on [0,2] and ends at 65/54. The upper side
follows by the opposite subtraction in (5). Finally subtract (8) for k=0
from k=1: c_0-c_1=77/1920 and |R_1-R_0| <= |R_1|+|R_0|, giving (2)–(3).

In the YM-line normalization of OL1/CB1,
theta = 4 theta_YM and A = H_theta/4 - 3 theta_YM. Thus this window is
0 <= theta_YM <= 1/2 and the A-gap is >= 65/216. The scalar shift cancels
in a gap. These are one-site dimensionless spectral units, not a physical
mass calibration.

## 8. Framework map and necessary correction to YC2's wording

The structural transfer is completely typed here:

| Framework role | Actual object in this theorem |
|---|---|
| Carrier and metric | Gauge-singlet L^2 on three unit S^3 links; product Haar inner product |
| Retained reading | Whole free lowest multiplet P_k, or one branch after exact reflection decomposition |
| Lost source | Q_k V P_k; scalar residual r = (V-a_k) f_0 |
| Hidden floor | e_k + delta_k for every state in the complement |
| Return | theta^2 P_k V Q_k (D_theta-E)^(-1) Q_k V P_k |
| Target | Actual compact-sector energy E_k and difference E_1-E_0 |

YC2's Haar calculation Var(Q)=19/64 is useful: V=2Q gives vacuum residual
square 4 Var(Q)=19/16. But the spectral return is not that static variance.
The kinetic operator weights it as
c_0=(3/4)/8+(7/16)/16=31/256. This is the missing dynamic information.
Dimensionless response ratios alone do not determine these denominators.
TVSP -> VTSP remains the coordinate exchange of a typed response chart; it
does not replace H_0 or turn a Gibbs Hessian into a quantum resolvent.

YC2 correctly provides the character Fourier transform of twisted heat traces
as a sector-resolving observer. Its stronger phrase that this observer is
*required* is too narrow: sector-restricted spectral projectors and the Schur
calculation here also resolve the same target. A centre-even configuration
potential has no character label by itself, but its **quantum expectations in
different sector states can differ**, as the exact slopes above show.
Similarly, centre symmetry alone does not prove CM2's equality of leading
weak-limit energies: CM2's localization/core matching does that work. Centre
symmetry is already present at theta=0 when E_k=3k are unequal. YC2's frozen
packet is left intact; this note records the qualification explicitly.

## 9. Scope, reproducibility and next step

This packet proves neither the large-theta absolute centre splitting nor its
tunnelling scale. It adds no volume-uniform estimate, no lattice refinement,
and no 4D continuum construction. It does not replace CM2's asymptotic matching
or DR2's all-column free-core result. No literature-priority claim is made for
the standard harmonic/perturbative machinery or these finite-site expansions.

Run from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python physics/yc3/yc3_sector_splitting.py --check
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s physics/yc3 -p 'test_*.py' -v
```

`YC3_RESULT.json` pins the relevant predecessor evidence, this note, the code
and the independent tests. All arithmetic determining constants is rational.
The replay checks identities and endpoint margins; the arguments above supply
infinite-dimensional completeness, interval monotonicity and spectral ordering.

The next bounded target is to enlarge the parity-resolved retained space and
certify a higher complement floor, then enclose centre-bottom splittings beyond
theta=2. For the CM2 weak limit, a separate inter-well/tunnelling remainder is
still needed; a finite-coupling power series is not a substitute.
