# SG1 — a certified gap after keeping the physical gauge-singlet target

Independent development, 9 October 2026. Built on CR1 and CM1 after reading
SL1 `4df12880d1f8f2a8567474e9323d06742350aa4b`.

**The free constant-mode physical core now has an explicit gap lower bound:**

| Number of turns d | Ground rate E0 | First excited rate E1 | Certified E1−E0 lower |
|---|---|---|---|
| 3 | 4.4187 ≤ E0 ≤ 5.1868 | 5.5210 ≤ E1 ≤ 7.5787 | **0.3342** |
| 4 | 6.7442 ≤ E0 ≤ 8.0030 | 8.0060 ≤ E1 ≤ 10.9738 | **0.0030** |

These are pure numbers for h below. CR1 supplies the ground uppers; CM1
supplies the excitation uppers. SG1 supplies the new ground/excitation lowers.
The new excitation lower applies to **every direction sector inside the
SO(3) gauge-invariant carrier**, including sectors absent from either trial
ladder. It does not identify which sector realizes E1. The d=4 gap bound is
small and conservative; it is not a computed approximation to the gap.

The key improvement is the target. CR1's roughly 5.1 comparison excitation
was a possible mode of the unrestricted comparison space. A single nonzero
orbital spin contains no gauge singlet. Once that forbidden mode is removed,
the same pair confiner has its second physical rate above CR1's ground upper.
No guessed vacuum, fitted Airy value or unproved hidden floor is used.

This is a fixed, finite-number-of-modes bosonic core result. Its physical
scale remains g^(2/3)/L. The continuum four-dimensional Yang–Mills gap remains
open, along with volume-uniform matching. Four turns here means four columns
of C and must not be equated with the continuum spacetime dimension.

## 1. The precise carrier and theorem

Let d=3 or 4, C=(c1,…,cd) in R^(3×d), and

\[
h=-\frac12\Delta_C+\sum_{i<j}|c_i\times c_j|^2.
\]

Use the Friedrichs realization of the nonnegative quadratic form initially
defined on smooth compactly supported functions. Restrict to

\[
\mathcal H_{\mathrm{phys}}
=\{\psi\in L^2(\mathbb R^{3d}):\psi(RC)=\psi(C)
\text{ for every }R\in SO(3)\}.
\]

This is the constant-mode SU(2) adjoint gauge-singlet condition. Column
rotations, column permutations and parity are not additionally imposed on
the target. The scalar CR1 trial functions belong to this space. CM1's
written confinement/positivity argument gives a discrete spectrum and a
unique positive ground on the declared free carrier. The estimates below
also separate its first two physical variational levels directly.

**Theorem.** For the ordered physical eigenvalues E0≤E1≤… (counting
multiplicity), the bounds in the opening table hold. In particular

\[
E_1-E_0\ge\frac{1671}{5000}\quad(d=3),\qquad
E_1-E_0\ge\frac3{1000}\quad(d=4).
\]

Analytical inputs: the closed-form/min–max principle, angular decomposition
in spherical harmonics, the one-dimensional Sturm oscillation theorem, and
Dirichlet–Neumann form bracketing. Their roles are specified below. They are
standard analytical interfaces; the finite code certifies the arithmetic
and ODE enclosures used by this argument, not these general theorems.

## 2. Keep the pair bound, then restrict the correct target

Let Ti=−Δci/2. For each ordered pair i≠j allocate

\[
\frac{1}{2(d-1)}T_j+\frac12|c_i\times c_j|^2.
\]

The two transverse oscillator ground energies sum to
\(|c_i|/\sqrt{2(d-1)}\); longitudinal kinetic energy is nonnegative.
Summing the allocations uses half the total kinetic energy and all the
quartic potential. Therefore, as quadratic forms,

\[
h\ge K_d:=\sum_{i=1}^d\left[-\frac14\Delta_{c_i}
          +\sqrt{\frac{d-1}{2}}\,|c_i|\right].
\tag{1}
\]

This is CR1's pair comparison, also written explicitly in CM1. Both sides
commute with simultaneous row rotations, so (1) restricts to H_phys. The
comparison Kd need not preserve column O(d) rotations; min–max only needs
the form order on the same physical carrier.

Scaling each radius turns every factor into

\[
c_d\,A_\ell,\qquad
c_d=\left(\frac{d-1}{8}\right)^{1/3},\qquad
A_\ell=-\frac{d^2}{dr^2}+\frac{\ell(\ell+1)}{r^2}+r,
\tag{2}
\]

on L2(0,∞), with its regular/Friedrichs condition at zero. Here u(r)=r f(r)
is the usual three-dimensional radial unitary transform. No singular-value
Jacobian of C is assumed.

## 3. The angular selection that changes the count

Let a0 and a1 be the first two eigenvalues of A0. SG1 certifies

\[
a_0\ge2.3381,\qquad a_1\ge4.0879.
\tag{3}
\]

For every ℓ≥1, let b be a lower bound for the bottom of Aℓ. We will use
b=3.245. A product angular channel has representation
V_(ℓ1)⊗…⊗V_(ℓd) under simultaneous SO(3). If exactly one ℓi is nonzero,
averaging that spherical harmonic over rotations gives zero: there is no
invariant vector. Thus every physical comparison channel is either:

1. all ℓi=0; or
2. at least two ℓi are nonzero.

In the first case the product radial ground is one-dimensional, and every
orthogonal radial excitation costs at least (d−1)a0+a1. In the second case
the energy is at least (d−2)a0+2b; extra nonzero spins only increase this
bound because a0≤5/2<b (the simple radial trial u=r exp(−r) gives 5/2).
Equivalently, after inserting the certified lower constants, each additional
nonzero spin replaces 2.3381 by 3.245. Hence

\[
E_1(K_d|_{\mathrm{phys}})\ge c_d\min\left\{
(d-1)a_0+a_1,\ (d-2)a_0+2b\right\}.
\tag{4}
\]

The individual angular labels diagonalize **Kd**, not necessarily h.
Superpositions of physical angular channels obey the same lower estimate;
the min–max transfer E1(h_phys)≥E1(Kd_phys) does not require h to preserve
every individual ℓi. This is why (4) covers all physical direction sectors.
The one-dimensional subspace removed for this comparison is Kd's product
ground, not an approximation asserted to equal h's vacuum.

### A uniform angular floor from a positive reading

For ℓ≥1 use φ(r)=r^(ℓ+1) exp(−k r^(3/2)), k²=2/9. Completing the square
against a compactly supported test function gives the ground-transform form
bound Aℓ≥inf(Aℓφ/φ). Direct differentiation gives

\[
\frac{A_\ell\phi}{\phi}
=\left(1-\frac{9k^2}{4}\right)r
 +\frac{3k(4\ell+5)}{4\sqrt r}.
\]

Minimizing Ar+B/√r gives 3(AB²/4)^(1/3). At the chosen k,

\[
\inf\operatorname{spec}A_\ell
\ge\frac34(4\ell+5)^{2/3}
\ge\frac34\,81^{1/3}>3.245.
\tag{5}
\]

The last comparison is exact after cubing:
(649/200)^3 < 2187/64. The identity extends from test functions to the
Friedrichs form by density. No assumption about a guessed eigenfunction is
needed. This is the positive-local-rate method already used in CR1, applied
to the correct nonzero angular sector.

## 4. Outward radial counts, independent of Airy tables

For (3), split the half-line at R=12 with a Neumann cut. On the tail,
−d²/dr²+r≥12, which is above both test energies. On [0,12], solve

\[
y''=(r-E)y,\qquad y(0)=0,\qquad y'(0)=1.
\tag{6}
\]

The Sturm count for Dirichlet at zero and Neumann at R is

\[
N_N(E)=\#\{y\text{ zeros in }(0,R)\}
       +\mathbf1_{y(R)y'(R)<0}.
\tag{7}
\]

Endpoints are verified nonzero. Formula (7) follows from the lifted Prüfer
angle atan2(y,y′), or from Dirichlet–Neumann interlacing with its endpoint
sign. Counting zeros alone would omit the possible final half-turn; a test
at E=2.3, R=2.5 explicitly detects that error.

| Test E | Interior zeros | Endpoint signs (y,y′) | N_N(E) |
|---:|---:|---|---:|
| 2.3381 | 0 | (+,+) | 0 |
| 4.0879 | 1, in (13/8,7/4) | (−,−) | 1 |

Neumann bracketing lowers the operator, so the half-line has at most these
counts below E. This proves (3). A separate count at 2.3382 has one interior
zero. The finite Dirichlet problem therefore has a trial below 2.3382;
extending it by zero gives a0<2.3382 on the half-line. That upper is only
used for the ungauged negative control below.

### The tail of each ODE step is bounded

For a cell based at a, write y(a+t)=Σ yn t^n. The exact recurrence is

\[
y_{n+2}=\frac{(a-E)y_n+y_{n-1}}{(n+2)(n+1)},\qquad y_{-1}=0.
\]

The code propagates rational intervals for (y,y′) with cell width h=1/8,
Taylor order N=32 and outward rounding to 160 binary places. It bounds the
whole cell as well as the endpoint. Each of the 96 cells either excludes a
zero of y or proves a strict sign for y′; endpoint signs then count crossings.
The initial Dirichlet zero is excluded. Any unresolved cell rejects the
certificate.

Here is the analytic remainder bound, so the truncated series is not being
treated as exact. On the complex disk |t|≤1, the first-order system matrix
for U=(y,y′) has infinity norm at most M=max(1,|a−E|+1). Integration along
a radial segment and Grönwall give ||U||≤exp(M)||U(0)||. Since e<4, set
B=4^ceil(M)||U(0)||, using an upper bound from the initial intervals.
Cauchy's coefficient estimate then gives, uniformly on |t|≤h,

\[
\left\|U(t)-\sum_{n=0}^N U_nt^n\right\|_\infty
\le \frac{B h^{N+1}}{1-h}.
\tag{8}
\]

Every arithmetic operation and rounding endpoint is rational and outward.
The JSON records the endpoint enclosures, zero-containing cells, certification
method counts, maximum local remainder and all parameters. Changing the mesh
to 1/16, order to 28 and precision to 144 places reproduces the counts and
overlapping endpoint enclosures. An independent Runge–Kutta implementation
cross-checks signs and counts in tests; its floats are not proof inputs.

For classical identification only, the eigenvalues of A0 are the absolute
negative zeros of Ai. [NIST DLMF §9.9](https://dlmf.nist.gov/9.9) supplies that
standard zero notation and numerical tables. No tabulated Ai value or special
function library enters this certificate. The linear-potential spectrum and
singlet angular selection are established mathematical tools; SG1's contribution
here is their explicit physical-carrier comparison and reproducible bounds
for the TC1/CR1 core.

## 5. Exact arithmetic closes both gaps

Use the rational lower scales
c3>0.629960 and c4>0.721124, certified by cubing against (d−1)/8.
Substitution of (3) and (5) into the two alternatives of (4) gives:

| d | Radial-excitation comparison lower | Angular-excitation comparison lower | Recorded E1 lower |
|---|---:|---:|---:|
| 3 | 5.521032436 | 5.561349876 | 5.5210 |
| 4 | 8.0060628728 | 8.0522148088 | 8.0060 |

The implementation replays the actual 67-reading CR1 Gram/Hamiltonian
inertia at E0 uppers 5.1868 and 8.0030: each has a negative direction.
These are valid physical trials, so subtraction of the **excitation lower**
and **ground upper** gives the gap bounds, with the correct inequality
directions. The smaller of two upper Ritz values is never used as a lower.

The same radial ground lower improves E0 to at least 4.418728428 (d3) and
6.7442400976 (d4). The initial table rounds those down.

### Negative control: the cheap unphysical mode really is below the threshold

For A1, u=r² exp(−a r) has Rayleigh value a²+5/(2a); choose a=27/25. Combine
it with d−1 radial ground factors, using a0<2.3382 and rational upper scales.
The resulting one-spin comparison trial is below 5.139 for d3 and below
7.569 for d4—both below CR1's ground upper. It is annihilated by gauge
averaging. Thus the previously observed unrestricted obstruction is real,
and its removal has a precise physical-target justification.

The numerical comparison did not need a sharper quartic approximation to
cross the ground upper once the gauge target was kept. CM1's source-memory
route remains useful for sharper bounds and matching, but a Schur hidden
floor is not an unproved hypothesis of SG1.

## 6. SL1 and the owner's TVSP → VTSP clarification

SL1 proves its pair-record inversion from the lattice diagonal substitutions,
the two arithmetic–geometric means and the Gaussian weak-end area. In
normalized readings f=F/S and l=L/S, its result is the exact exchange

\[
(f,l)(1/u)=(l,f)(u),\qquad f^2+l^2=1.
\]

This is read from SL1 at the pinned commit; it is not used to manufacture the
radial constants in SG1. Its proof supplies a concrete seam-exchange example
for the framework lens. The weak-end squared lower squeeze in SL1 is needed
only after u/2^n≤1; for larger values one keeps max(0,1−√t)^2 as the lower
bound. That qualification preserves its limit argument without promoting
an invalid squaring step for every t>0.

The owner's clarification on this date is recorded as the intended reading:
**TVSP is the pre-observation diagram; after observation T and V exchange
positions, and VTSP is the intended chart for observed phenomena.**
This is research interpretation supplied by the owner. Its algebraic chart
map can be written explicitly as

\[
x_{\rm pre}=(T,V,S,P)^T,\quad
x_{\rm obs}=\Pi x_{\rm pre}=(V,T,S,P)^T,\qquad
\Pi=\begin{pmatrix}0&1&0&0\\1&0&0&0\\0&0&1&0\\0&0&0&1\end{pmatrix},
\quad\Pi^2=I.
\]

Coordinate labels and their units travel together. For a common linear
four-slot representation, use dimensionless coordinates or the typed bundles
specified in RKF T42. If both state and response frames are reordered by Π,
their Jacobian transforms as J_obs=Π J_pre Π^T; the metric and cut transform
by the same congruence. Source, pairing and gauge action must travel with the
reading. This gives a precise meaning to a complete observer chart. The
physical observation dynamics and a thermodynamic-to-core intertwiner remain
to be supplied; the permutation alone is not such an intertwiner.

In SG1 the operative lesson is to keep the physical recognition target while
counting. The exact singlet condition removes the comparison reading that
caused the apparent obstruction. The TVSP/VTSP interpretation is not an
extra assumption in the numerical gap proof.

## 7. Claim ledger and continuation

| Claim | Status |
|---|---|
| Pair confiner on the same physical gauge carrier | Written form inequality; inherited CR1/CM1 allocation |
| Single nonzero orbital spin is excluded | Exact rotation-average argument |
| All physical comparison channels covered by the two alternatives | Written angular decomposition and min–max proof |
| Radial spectral thresholds, all ODE step errors, correct endpoint count | Outward rational certificate |
| Positive fixed-core gaps 0.3342 / 0.0030 | PROVED on the declared physical carrier, with named analytical inputs |
| Other direction-sector lower control needed for this core gap | DISCHARGED by the uniform singlet comparison |
| Exact lowest excited sector or sharp gap values | OPEN |
| Source-preserving valley approximation and sharper quantitative return | OPEN, no longer required to obtain a positive fixed-core gap |
| g^(2/3)/L matching, volume uniformity and continuum Yang–Mills gap | OPEN |
| TVSP before / VTSP after observation as physical identification | Owner's research interpretation; chart swap specified |

The next physical obligation is no longer simply “find a positive pure-number
core gap.” It is to connect the finite-mode gap to the full interacting theory,
including nonconstant modes, constraints, renormalization/matching and a bound
that survives the intended volume and cutoff limits. A sharp core certificate
can be pursued in parallel as a mathematical objective, but its fixed-g
1/L scaling alone cannot supply that limit.

## Sources and replay

- TC1/CR1 baseline: extra-ideas `80002cb2eed9fcbd810cb00477326ed146065072`.
- CM1: `de14f48c19cb5f84baac3e3d13e7de158317c738`,
  [core sector/source note](../cm1/CM1_CORE_SECTOR_MEMORY.md).
- SL1, fully read note, implementation and tests:
  [4df1288 seam-law packet](https://github.com/Parveen117/extra-ideas/blob/4df12880d1f8f2a8567474e9323d06742350aa4b/physics/sl1/SL1_SEAM_LAW_FROM_THE_WALK.md).
- RKF T42 at `e8e3089745bf7dc0532d01c62a5f660f5658c6f0`, sections 1–2 for
  the typed thermodynamic chart; its physical identification remains separate.
- [SG1_RESULT.json](SG1_RESULT.json) contains exact source SHA-256 pins,
  rational outcomes and outward ODE certificates. Earlier stage packets and
  their stored claim boundaries are preserved as historical checkpoints.

```sh
python3.12 -B physics/sg1/sg1_singlet_gap.py --check
python3.12 -B -m unittest discover -s physics/sg1 -p 'test_*.py'
```

Use `--write` only to intentionally regenerate the result. `--check` recomputes
the certificate and fails on stale outcomes or source pins.

## Later note (GC1 and OL1, 9 October)

[GC1](../gc1/GC1_CORE_GAP_COUNT.md) was written at the same time on the physics branch, without sight of this
note, from the same comparison. The two derivations agree to the last digit kept: first excited rate at least
5.5210 and 8.0060, gap at least 0.3342 and 0.0030 (GC1 keeps 0.0031, from its own reading under 8.0029). The
routes differ where it matters for a check: here an interval solve with a Neumann cut at 12 and a count with
its endpoint term; there one exact power series on a grid of 1/256 and the sign kept past the turning point.
GC1 adds the lowest rate from both sides, [5.1865, 5.1868] and [7.9942, 8.0029], and the walk in the number of
turns. The qualification of SL1's weak-end step in section 6 is taken into SL1 as a later note.
[OL1](../ol1/OL1_ONE_SITE_LATTICE.md) carries the count to the compact links of the one-site lattice. Nothing
above is changed.

## Later note (DR1, 9 October)

[DR1](../dr1/DR1_CORE_BY_ROWS.md) reads the weight by rows — three cuts in d directions — and counts with the
three half turns of the gauge group. Three turns: the same 5.5210. Four turns: first excited rate at least
8.4454, gap at least 0.4425. Five to ten turns: first floors. Nothing above is changed.
