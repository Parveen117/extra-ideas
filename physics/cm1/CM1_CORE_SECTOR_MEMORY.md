# CM1 — turning sectors and the core's missing-source memory

Independent continuation of TC1 → CR1, 2026-10-09. The fetched physics branch
ends at CR1 `80002cb2eed9fcbd810cb00477326ed146065072`; the proposed
valley-observer remainder stage in the owner's conversation is not a pushed
packet at this checkpoint. CM1 complements that route with a direction sector
and the actual retained-to-hidden source. It does not modify TC1 or CR1.

**Result.** A gauge-invariant quadrupole trial gives first-excitation upper
bounds **7.5787** for three turns and **10.9738** for four turns. CR1's scalar
second-rate upper bounds were 8.0480 and 11.5660. These are improved upper
bounds, not an ordering proof between the exact scalar and tensor energies.
The scalar degree-10 cut has 67 readings and an actual missing-source rank of
**26**, inside a 35-dimensional next shell. Its source Gram must use the
Gaussian metric. A bare hidden-floor comparison becomes more demanding with
larger trial cuts in the recorded examples. No numerical excitation lower
bound is obtained.

There is also a written proof, using declared standard analytical inputs, that
this **free constant-mode bosonic core** has compact resolvent and a unique
positive ground state. Its fixed-core gap therefore exists qualitatively.
This gives no explicit gap lower number, volume-uniform estimate, or continuum
Yang–Mills mass gap. The finite certificate does not certify the analytical
compactness and positivity theorems.

## 1. Carrier, cut, target and metric

The declared realization is

\[
C=(c_1,\ldots,c_d)\in\mathbb R^{3\times d},\quad d=3,4,\qquad
X=\sum_{i<j}|c_i\times c_j|^2=e_2(CC^T),\qquad
h=-\tfrac12\Delta_C+X.
\]

Use the Friedrichs operator of the closure of
\(q[\psi]=\frac12\int|\nabla\psi|^2+\int X|\psi|^2\) on
\(C_c^\infty(\mathbb R^{3d})\). Physical trial readings here are invariant
under the row action of SO(3). Column O(d) rotates direction sectors. The
scalar cut retains polynomials in \(e_1,e_2,e_3\), times
\(g_w=\exp(-w e_1/2)\). The pairing is the actual Euclidean C-space pairing,
or its normalized Gaussian form \(\langle p,q\rangle_w=\mathbb E_w[pq]\).
The scalar Gram S is not the identity in monomial coordinates.

CM1's targets are (a) valid excitation trial upper bounds outside CR1's scalar
cut; (b) the full source \(QhP\) and its Gram for the declared scalar cut;
(c) the precise extra hidden control needed by a Schur lower comparison.
No singular-value change-of-variables formula is assumed. Gaussian moments
come from CR1's integration-by-parts recursion and direct entry Wick checks.
No claim is made that this polynomial family exhausts the gauge-invariant
Hilbert space: determinant/orientation and other direction sectors remain.

The classical Schrödinger realization is an explicit analysis interface. In
the framework lens the primitive obligations remain retained reading, full
source, hidden memory, common metric, and an absolute lower floor. The operator
and measure declared above do not become new native axioms.

## 2. A turning reading that the scalar ladder misses

Write \(M=C^TC\), \(p_r=\operatorname{tr}M^r\), and

\[
t_r=(M^r)_{11}-(M^r)_{22}.
\]

The t's are row-gauge invariant, but belong to a traceless symmetric tensor
under column O(d). For d=3 this is spin 2 (five components); for d=4 the tensor
has nine components. One component and its invariant polynomial multiples
suffice for an upper bound. Interchanging columns 1 and 2 changes their sign,
so their average over directions is zero.

Set \(L=\Delta_C/2\). Differentiating the entries of M gives

\[
Lt_1=0,\qquad Lt_2=(d+8)t_1,\qquad
Lt_3=(2d+15)t_2+3e_1t_1.
\]

With \(\Gamma(f,g)=\nabla f\cdot\nabla g\), the cross derivatives are

\[
\begin{aligned}
\Gamma(t_r,e_1)&=4r t_r,\\
\Gamma(t_r,e_2)&=4r(e_1t_r-t_{r+1}),\\
\Gamma(t_r,e_3)&=4r(e_2t_r-e_1t_{r+1}+t_{r+2}).
\end{aligned}
\]

These follow from
\(\partial_C e_k=2C(e_{k-1}I-e_{k-2}M+\cdots+(-1)^{k-1}M^{k-1})\)
and differentiating \(\operatorname{tr}(A M^r)\) with
\(A=\operatorname{diag}(1,-1,0,\ldots)\). Thus
\(L(t_r p)=t_rLp+pLt_r+\sum_k\Gamma(t_r,e_k)\partial_kp\), retaining
every cross derivative. Cayley–Hamilton reduces t3 to e1 t2 − e2 t1 at d=3;
at d=4 the rank-three identity reduces t4 to e1 t3 − e2 t2 + e3 t1.

For a homogeneous polynomial p of C-degree 2m,

\[
g_w^{-1}h(pg_w)=-Lp+2wm p+
\left(\tfrac{3dw}{2}-\tfrac{w^2e_1}{2}+e_2\right)p.
\]

The exact entry-coordinate Laplacian is checked independently against these
reduced formulas, including mixed e1, e2 and e3 derivatives. Checking selected
fixtures is evidence for the implementation; the displayed differentiation
identities supply the general polynomial argument.

### Exact pairing without choosing a singular-value measure

For fixed eigenvalues, rotation averaging gives

\[
\mathbb E_{O(d)}[t_rt_s]
=\frac{4}{(d-1)(d+2)}\left(p_{r+s}-\frac{p_rp_s}{d}\right)=W_{rs}.
\]

To see the coefficient, the averaged traceless tensors have covariance
\(b(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}
-2\delta_{ij}\delta_{kl}/d)\). Contracting all i,j determines
\(b=(p_{r+s}-p_rp_s/d)/((d-1)(d+2))\); contracting with A twice gives 4b.
The full Gaussian measure is O(d)-invariant, so
\(\langle t_rp,t_sq\rangle_w=\mathbb E_w[W_{rs}pq]\). Newton identities
express every p_r in e1,e2,e3. This is valid even when d=4 and M has a zero
eigenvalue. Direct entry Wick pairings check the formula separately.

### Exact inertia certificate and scope

Retain \(t_r e_1^a e_2^b e_3^c\) with \(r+a+2b+3c\le D\), using
r=1,2 for d=3 and r=1,2,3 for d=4. CR1's widths are w=8/5 and 9/5.
All S and H entries are rational. Exact positive LDL pivots verify S>0;
exact inertia of H − μS verifies a trial eigenvalue below μ. Nested cuts give:

| d | D | Readings | Lowest finite Ritz bracket | True tensor energy upper |
|---|---:|---:|---:|---:|
| 3 | 2 | 3 | (7.9554, 7.9555) | 7.9555 |
| 3 | 4 | 11 | (7.6272, 7.6273) | 7.6273 |
| 3 | 6 | 27 | (7.5842, 7.5843) | 7.5843 |
| 3 | 8 | 54 | (7.5786, 7.5787) | 7.5787 |
| 4 | 2 | 3 | (11.2989, 11.2990) | 11.2990 |
| 4 | 4 | 13 | (11.0060, 11.0061) | 11.0061 |
| 4 | 6 | 34 | (10.9766, 10.9767) | 10.9767 |
| 4 | 8 | 70 | (10.9737, 10.9738) | 10.9738 |

The left endpoint bounds only the **finite Ritz value**. It is never a lower
bound for the infinite operator. Nor do two sector upper bounds order the
exact sector energies. Section 3 shows why these readings are orthogonal to
the true ground state; hence they also bound the full gauge-invariant first
excitation from above. Another physical sector may be lower.

## 3. Fixed-core compactness and the ground state

This is a written argument, separate from finite arithmetic. Its standard
analytical inputs are closed nonnegative form/Friedrichs realization, local
Rellich compactness, local elliptic regularity and the strong maximum principle
(equivalently local Harnack positivity). They are declared imports, not
theorems proved by the Python checks.

Let \(T_i=-\Delta_{c_i}/2\) and \(0<\alpha<1\). Allocate each ordered pair
\(i\ne j\) the operator
\(aT_j+\frac12|c_i\times c_j|^2\), where
\(a=(1-\alpha)/(d-1)\). Their kinetic terms total (1−α)T and their potentials
total X. For fixed ci the two transverse oscillators have joint ground
energy \(\sqrt a\,|c_i|\); the longitudinal kinetic term is nonnegative.
Applying that fiber inequality and summing yields the form bound

\[
h\ \ge\ \alpha T+\sqrt{(1-\alpha)(d-1)}\sum_i|c_i|.
\]

At α=1/2, set \(\beta=\sqrt{(d-1)/2}>0\). A unit vector with q[ψ]≤E
satisfies

\[
\int|\nabla\psi|^2\le4E,\qquad
\int_{|C|>R}|\psi|^2\le\frac{E}{\beta R},
\]

because Σ|ci|≥|C|. Rellich compactness on each bounded ball and this uniform
tail estimate imply compact embedding of the form domain into L2, hence
compact resolvent. This argument retains kinetic control; potential growth
alone would fail on the rank-one valley.

The lowest eigenvalue is attained. Taking |ψ| does not increase the form,
so a real nonnegative minimizer exists. Local regularity and positivity on the
connected space make it strictly positive. If the ground eigenspace had an
orthogonal real vector, that vector would change sign; its absolute value is
also a ground minimizer but would vanish at a nodal point, contradicting
strict positivity. Thus the ground is simple.

Every row SO(3) and column O(d) rotation preserves h and takes the normalized
positive ground into itself. The ground is consequently gauge invariant and
a direction scalar. A turning trial changes sign under a column swap while
the ground does not, proving their exact orthogonality. Min–max now transfers
the turning upper bounds to the physical first excitation.

Compact resolvent and a simple ground imply E1>E0 for this fixed bosonic
core. This closes CR1's named discreteness/ground-realization assumption for
this declared carrier, with the stated analytical inputs. It does **not**
close the missing explicit lower estimate for E1−E0. It is classical
confinement reasoning, not a claimed new solution of the YM problem.

## 4. The actual source and its metric

Let Φ be the scalar Gaussian basis of weighted degree ≤D. Let P be the
orthogonal projection onto its span and Q=I−P. In this nonorthogonal basis,

\[
S=\Phi^*\Phi,\quad H=\Phi^*h\Phi,\quad
T= (h\Phi)^*(h\Phi),\quad
K=h\Phi-\Phi S^{-1}H,\quad
G=K^*K=T-HS^{-1}H.
\]

The columns K are the **actual** QhΦ source, not a surrogate matrix. The action
lies in the polynomial ladder through degree D+2. If K is represented there
by coefficients k, its Gram is \(k^*S_{D+2}k\). CM1 checks exact equality with
T−HS^-1H and exact retained-source orthogonality. Replacing this by k*k in
Euclidean coefficient coordinates fails the negative controls.

Source columns of degree ≤D−2 vanish exactly, since h maps them inside P.
Above P the highest term of h is multiplication by e2. On the homogeneous
degree-D boundary it produces an injective degree-D+2 image; after that
boundary vanishes, the degree-D−1 boundary has an injective degree-D+1 image.
Thus the above-cut quotient has rank n(D)−n(D−2). Orthogonal projection does
not change that quotient, and the remaining columns vanish, proving

\[
\operatorname{rank}(QhP)=n(D)-n(D-2),\quad
n(D)=\#\{(a,b,c)\ge0:a+2b+3c\le D\}.
\]

The invariant polynomials are independent: positive distinct nonzero squared
singular values give an open set of (e1,e2,e3), so no nonzero polynomial
relation holds. This justifies the injective multiplication step.

| D | Retained n(D) | Next-shell dimension n(D+2)−n(D) | Actual source rank |
|---:|---:|---:|---:|
| 2 | 4 | 7 | 3 |
| 4 | 11 | 12 | 7 |
| 6 | 23 | 18 | 12 |
| 10 | 67 | 35 | 26 |

Full Gram reconstruction is recorded through D=6 for both d. At D=10 the
**actual above-cut action coefficient matrix** is exactly rank-checked; the
67-column full Gram is not computed in this packet. The rank is for this cut,
not a claim that the entire hidden dynamics is finite dimensional.

## 5. Observer completion as an exact return, with an open floor

In orthogonal retained/hidden coordinates write
\(h=\begin{pmatrix}A&B\\B^*&\mathcal D\end{pmatrix}\).
For z below the hidden spectrum the retained equation is

\[
\left(A-z-B(\mathcal D-z)^{-1}B^*\right)x
=f-B(\mathcal D-z)^{-1}g.
\]

The exact missing contribution is energy dependent and generally nonlocal.
It dresses both the retained operator and the source. The owner’s observer
idea therefore has a precise contract: reproduce this return in the retained
metric, including its source, and certify the approximation remainder. A
rotated chart or a local transverse ground potential alone does not discharge
that contract; TC1's valley potential can serve as an approximation to test.

If the **actual** hidden operator satisfies \(\mathcal D\ge d_0I\) with
d0>z, then its Schur complement is bounded below by
\(A-z-BB^*/(d_0-z)\). Congruence into the Φ coordinates gives

\[
H-zS-\frac{G}{d_0-z}.
\]

An exact negative count ≤1 here, together with the real hidden floor, would
bound the second **scalar-sector** rate below by z. CM1 bisects a sufficient
conditional floor for this bare-Gram comparison; it does not produce that
physical floor:

| d | z | D | Sufficient conditional d0 (rounded upward) |
|---|---:|---:|---:|
| 3 | 6 | 2 | 8.6600 |
| 3 | 6 | 4 | 12.8590 |
| 3 | 6 | 6 | 20.0643 |
| 4 | 9 | 2 | 11.4127 |
| 4 | 9 | 4 | 14.8972 |
| 4 | 9 | 6 | 20.4827 |

The JSON stores rational endpoints and counts. The global CR1 ground lower
bounds 4.14 and 6.32 lie below z: they supply neither hidden invertibility
there nor these sufficient floors. A hypothetical scalar bound would still
leave other physical sectors to control. Increasing the ladder alone makes
this crude floor requirement worse in these examples, despite improving the
Ritz upper. That is a measured limit of the adapter, not a no-go theorem for
the observer route.

The next constructive move is the source-resolved hidden solve already
identified in the [transfer map](../RH_YM_TRANSFER_MAP.md) and Publications
YM75: for F=B* and T=𝒟−z≥sI, a trial hidden solve Y has R=F−TY and
\(F^*T^{-1}F=(F^*Y+Y^*F-Y^*TY)+R^*T^{-1}R\).
Its residual bound is \(R^*R/s\) in the **correct** metric. Both a true s and
the full residual Gram remain obligatory. T28's recognition/memory norms and
tails must be attached to an explicitly bounded/resolvent or form adapter;
an unbounded h cannot simply be inserted into a bounded norm theorem.

## 6. Normalization, classical credit and claim ledger

TC1's dilation still gives \(H_g\simeq g^{2/3}h\). Physical torus rates need
the factor 1/L and matching. A fixed pure-number core gap consequently scales
as g^(2/3)/L and gives no volume-uniform positive limit at fixed g. Four turns
means four columns here, not a constructed four-dimensional continuum theory.

The spatially constant SU(2) matrix Hamiltonian, its rotation sectors and its
g^(2/3) law are classical. H.-P. Pavel's primary paper,
[SU(2) Yang–Mills quantum mechanics of spatially constant fields](https://arxiv.org/pdf/hep-th/0701283),
hep-th/0701283v1 (2007), eqs. (3)–(5), (15) and Table 2, uses potential X/2
at unit coupling. Dilation converts its unit energies into ours by 2^(1/3).
Its first spin-2 value 6.014500 is consistent with CM1's 7.5787 upper, divided
by that factor (about 6.0152). This is a literature cross-check, not an
outward lower certificate or a claim of new priority for the sector.
Flat-valley discrete spectra also have classical precedent in Barry Simon,
1983, DOI [10.1016/0003-4916(83)90057-X](https://doi.org/10.1016/0003-4916(83)90057-X).
Only its publication metadata was consulted here; the compactness proof used
in this packet is the explicit pair-form/Rellich argument above.

| Statement | Status |
|---|---|
| Turning generator and Gaussian tensor pairing | Written identities; exact independent entry checks |
| Turning trial upper values, nonorthogonal S>0 and inertia | Exact rational finite certificate |
| Compact resolvent and unique ground on the declared free core | Written proof with standard analytical inputs |
| Qualitative fixed-core gap existence | Consequence of that written proof |
| Actual scalar source Gram through D6 and source rank through D10 | Exact polynomial/rational certificate |
| Sufficient floor demands | Exact finite conditional comparisons; no actual hidden floor |
| Numerical second-rate/gap lower, lowest-sector identification | OPEN |
| Compact-block logarithm, matching, volume-uniform/continuum gap | OPEN |

Frozen sources: TC1/CR1 at extra-ideas CR1 commit above; RKF engine/theorem
catalog at e8e3089745bf7dc0532d01c62a5f660f5658c6f0; Publications YM75 at
e75e453b2fdb082570aa24e25a10488ec5d87b41. The transfer map names exact paths
and distinguishes native statements, classical interfaces and open bridges.
CM1_RESULT.json pins the exact local TC1/CR1 files, current map, this note,
implementation and tests by SHA-256. It does not hash itself.

Reproduce with Python 3.12, stdlib only:

```sh
python3.12 -B physics/cm1/cm1_core_sector_memory.py --check
python3.12 -B -m unittest discover -s physics/cm1 -p 'test_*.py'
```

`--check` recomputes and rejects stale results or source pins; `--write`
regenerates intentionally. The next handoff is source-preserving valley
trial solves, a genuine hidden lower bound and metric-aware outward residuals,
followed by the remaining physical sectors and then volume matching.

## Later note (GC1, 9 October)

[GC1](../gc1/GC1_CORE_GAP_COUNT.md) supplies the two lines left open in the ledger above. The comparison of
CR1-R4, counted under the row condition, has at most one reading under 5.5210 (three turns) and 8.0060 (four):
an explicit lower number for the first excitation, and a gap of at least 0.3342 and 0.0031. It is an actual
floor: for the cut of one reading (the product of the lowest one-turn readings) everything else is at or above
that line, and for a cut holding more products of one-turn readings the floor is the first level left out. The
mean of h² of §4, taken at degree 10, pins the lowest rate to [5.1865, 5.1868] and [7.9942, 8.0029]; with the
turning readings of §2 the gap is in [0.3342, 2.3922] and [0.0031, 2.9796]. The floor is far from the ceiling;
a cut of products of one-turn readings with this source Gram is the next step. Nothing above is changed.
