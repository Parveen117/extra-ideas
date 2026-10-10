# YC5 — resolve the lost source, retain its first interaction

**Status:** one-site, all-centre-sector spectral enclosure, with a written
full-space comparison proof and replayable exact rational certificates.
Not a volume-uniform or four-dimensional Yang–Mills mass-gap theorem.

YC5 continues [YC4](../yc4/YC4_HARMONIC_RETURN.md), local commit `1172772`,
which continues YC3 and the repaired CM2/DR2 baseline integrated at `8f30b5a`.
The predecessor packets are unchanged. This stage does **not** enlarge YC4's
retained cut. It reads the complete lost source at its own free harmonic
energies, then retains one actual hidden-potential interaction.

## 1. Result and normalization

On the actual compact gauge-invariant carrier

\[
\mathcal H=L^2(SU(2)^3)^{\operatorname{Ad}SU(2)},\qquad
H_\theta=-\sum_{i=1}^3\Delta_{S^3,i}+\theta V,
\qquad V=2\sum_{i<j}|u_i\times u_j|^2,
\]

where \(U_i=(a_i,u_i)\in S^3\), the pointwise bounds are \(0\le V\le6\).
The centre flips are the three independent \(U_i\mapsto-U_i\).
Write \(E_{k,j}\) for ordered eigenvalues, counted with multiplicity, in
one representative centre sector with \(k\) odd links; \(j=0\) is its bottom.
Link permutations identify the three \(k=1\) sectors and the three \(k=2\)
sectors. The operator has compact resolvent; these eigenvalues are well-defined.

**Theorem YC5.** For every \(0\le\theta\le12\), the full gauge-invariant
ground state lies in the all-centre-even sector and is simple. The first
excited level consists precisely of the three equivalent, individually simple,
single-centre-odd bottoms. In particular,

\[
\Delta_{\rm all}(\theta)=E_{1,0}(\theta)-E_{0,0}(\theta)
\ge\frac1{20}.
\]

The exact point enclosures are:

| \(\theta\) | Lower bound on \(\Delta_{\rm all}\) | Upper bound |
|---:|---:|---:|
| 2 | 2.1607 | 2.1613 |
| 4 | 1.6020 | 1.6179 |
| 6 | 1.2120 | 1.3201 |
| 8 | 0.8801 | 1.2516 |
| 10 | 0.5244 | 1.4239 |
| 12 | 0.0863 | 1.8512 |

All displayed decimals terminate and denote exact rationals, not estimated
digits. The wide last interval is an error enclosure, not evidence that the
physical gap grows. YC4's lower bound at \(\theta=8\) was 0.3063; the new
one is 0.8801 with the **same retained ranks**. YC4's stronger uniform bound
\(1/5\) remains available on its smaller interval \([0,8]\).

In the earlier YM-line normalization \(A_{\theta_{\rm YM}}),
\(H_\theta=4A_{\theta/4}\). Thus the new window is
\(0\le\theta_{\rm YM}\le3\), with gap at least \(1/80\) in \(A\) units.
These are dimensionless finite-site statements. CB1 already gave positivity
for every finite coupling: the new content is a sharper quantitative return
certificate and an extended identification window, not first-ever positivity.

## 2. The cut, source and metric stay explicit

In each centre sector use YC4's projection \(P_k\) onto **all** gauge-singlet
free harmonic shells below \(\Lambda_k=3k+24\), and put \(Q_k=1-P_k\).
YC4 proves completeness of these shells; the ranks and full hidden floors are

| Odd links \(k\) | Retained rank | \(Q_kH_0Q_k\ge\Lambda_k Q_k\) |
|---:|---:|---:|
| 0 | 13 | 24 |
| 1 | 21 | 27 |
| 2 | 32 | 30 |
| 3 | 23 | 33 |

Scalar reflections \((a_i,u_i)\mapsto(-a_i,u_i)\) further split each cut
into invariant blocks. They are not centre flips. Every reflection block is
included in the counts below. In one such block, let \(f_i\) be the exact
Haar-orthogonal, unnormalized retained basis and define

\[
G_{ij}=\langle f_i,f_j\rangle,quad
K_{ij}=\langle f_i,H_0f_j\rangle,quad
M_{ij}=\langle f_i,Vf_j\rangle,quad
T_\theta(z)=K+\theta M-zG.
\]

The complete lost source is \(r_i=QVf_i\), not a finite selection of its
components. Its Gram matrix is YC4's
\(B=W-MG^{-1}M\), where \(W_{ij}=\langle Vf_i,Vf_j\rangle\).
Write \(D_0=QH_0Q\), \(D_\theta=QH_\theta Q\). For \(z<\Lambda\),
the exact retained Schur form is

\[
S_\theta(z)=T_\theta(z)-\theta^2\Sigma_\theta(z),\qquad
(\Sigma_\theta(z))_{ij}=\langle r_i,(D_\theta-z)^{-1}r_j\rangle.
\]

The native map is therefore: **carrier** = actual compact gauge singlets;
**cut** = complete harmonic shells; **source** = \(QVf\); **metric** = Haar
\(G\); **target** = actual sector eigenvalues; **return** = the hidden
resolvent. The explicit block elimination/reconstruction map, not a diagram
resemblance or a substitution of thermodynamic labels, supplies the transfer.

## 3. Full harmonic source measure

Every \(r_i\) is a polynomial on \((S^3)^3\), hence has finite free-harmonic
support. Resolve it exactly as

\[
r_i=\sum_\lambda r_{\lambda i},\qquad
H_0r_{\lambda i}=\lambda r_{\lambda i},\qquad
(B_\lambda)_{ij}=\langle r_{\lambda i},r_{\lambda j}\rangle.
\]

The possible link degrees are the monomial degrees decreased by even integers;
their free energies are \(\sum_i n_i(n_i+2)\). The replay constructs spectral
projections by Lagrange polynomials in \(H_0\). It independently checks every
eigen-equation and the full reconstruction in exact Haar squared norm.
Only components of zero norm are removed. Consequently

\[
B_\lambda\succeq0,\qquad \lambda\ge\Lambda,\qquad
\sum_\lambda B_\lambda=B,\qquad
R(z):=\langle r,(D_0-z)^{-1}r\rangle
=\sum_\lambda\frac{B_\lambda}{\lambda-z}.
\]

Thus \(R\) is the Stieltjes transform of the **full free source measure**.
It is not the interacting hidden resolvent. Since
\(D_0\le D_\theta\le D_0+6\theta\), inverse order gives

\[
\sum_\lambda\frac{B_\lambda}{\lambda+6\theta-z}
\preceq\Sigma_\theta(z)\preceq R(z).
\tag{1}
\]

These comparisons hold without any commutation assumption. The resulting
lower Schur form is stronger than YC4's, since
\(R(z)\preceq B/(\Lambda-z)\). The upper Schur form in (1) improves plain
Ritz. No componentwise minimum or entrywise matrix ordering is used.

## 4. One actual hidden interaction, without false closure

Refine only the reflection blocks \(000\) in \(k=0\) and \(100\) in
the chosen \(k=1\) representative. Their free-spectral source spans have
dimensions 27 and 33. Exact Gram–Schmidt within each free energy yields
Haar-orthogonal vectors \(h_a\) with norms \(N_a\) and energies \(\lambda_a\).
Coordinates \(c_{ai}=\langle h_a,r_i\rangle/N_a\) represent every source
component. Compute the actual form \(M^h_{ab}=\langle h_a,Vh_b\rangle\).

This is a representation for computing

\[
C(z)_{ij}=\langle (D_0-z)^{-1}r_i,
 V(D_0-z)^{-1}r_j\rangle
=\sum_{a,b}\frac{c_{ai}M^h_{ab}c_{bj}}
{(\lambda_a-z)(\lambda_b-z)}.
\tag{2}
\]

It is **not** an assertion that \(V\) preserves this finite span. Higher
hidden interactions are still controlled on the full infinite complement.
All selected vectors are in \(Q\), so inserting \(QVQ\) instead of \(V\)
in (2) is exact. The replay checks the Haar metric, reconstruction of every
\(B_\lambda\), and \(0\preceq M^h\preceq6\operatorname{diag}(N_a)\).

Set

\[
A=\theta(D_0-z)^{-1/2}QVQ(D_0-z)^{-1/2},\qquad
0\le A\le uI,\quad u=\frac{6\theta}{\Lambda-z}.
\]

For \(0\le x\le u\) and any \(s\ge0\), the scalar resolvent satisfies

\[
\frac{1+2s-x}{(1+s)^2}\le\frac1{1+x}
\le1-\frac{x}{1+u}.
\tag{3}
\]

The differences are respectively
\((x-s)^2/[(1+x)(1+s)^2]\) and
\(x(u-x)/[(1+u)(1+x)]\). Functional calculus in the **single positive
operator \(A\)** proves (3) for operators. This does not diagonalize
\(D_0\) and \(QVQ\) simultaneously. Compressing with
\((D_0-z)^{-1/2}r\), choose \(s=2\theta/(\Lambda-z)=u/3\) and obtain

\[
\frac{(1+2s)R(z)-\theta C(z)}{(1+s)^2}
\preceq\Sigma_\theta(z)
\preceq R(z)-\frac{\theta C(z)}{1+u}.
\tag{4}
\]

The tangent lower bound can be negative in some directions; the inequality
is still valid. It need not dominate the lower return in (1). We do not
assume it does. With (4), the two certified Schur forms are

\[
L=T-\theta^2R+\frac{\theta^3}{1+u}C,
\qquad
U=T-\theta^2\frac{1+2s}{(1+s)^2}R
       +\frac{\theta^3}{(1+s)^2}C,
\qquad L\preceq S\preceq U.
\tag{5}
\]

The lower form in (5) dominates the resolved lower form and YC4's lower
form. In all other reflection blocks use (1). This is an analytic bound
on the full return, not perturbative truncation at order \(\theta^3\).

## 5. Full-space signs and coverage of a continuous interval

Because \(D_\theta-z>0\), exact Schur congruence gives

\[
N_-(U)\le N_-(H_\theta-z)=N_-(S)\le N_-(L).
\tag{6}
\]

In the middle equality, the right count is summed over all retained blocks
of the sector. Negative inertia of \(L\) equal to zero proves a ground
lower bound; at most one negative direction proves a second-level lower
bound. A nonpositive direction of \(U\) proves a ground upper bound:
reconstruct the hidden component as
\(-\theta(D_\theta-z)^{-1}r(v)\). Its full quadratic form equals the
Schur form and is nonpositive. This handles exact zero eigenvalues at
\(\theta=0\) as well. No finite-span invariance is needed.

`YC5_ENDPOINTS.json` records 105 rational endpoints: eighths from 0 to 11,
then sixteenths to 12. At each node, (6) proves four sector-bottom
enclosures \([l_k,u_k]\) and two second-level lower bounds \(b_0,b_1\).
No root-finding accuracy or monotonicity of the comparison pencils in
\(z\) is used in this proof: only the stored rational signs are replayed.

Since \(V\ge0\), **each ordered eigenvalue**, not necessarily a difference
of eigenvalues, is nondecreasing in \(\theta\). On a mesh cell \([a,b]\),

\[
E_{1,0}(\theta)-E_{0,0}(\theta)\ge l_1(a)-u_0(b),
\]

and separation from the other centre bottoms and the second levels follows
from \(\min_{k=2,3}l_k(a)-u_1(b)\) and
\(\min_{k=0,1}b_k(a)-u_1(b)\). Across all 104 cells the respective minima
are exactly

\[
\frac{103}{2000}=0.0515>\frac1{20},\qquad
\frac{2557}{10000}=0.2557>0,\qquad
\frac{979}{500}=1.958>0.
\]

These also establish simplicity within each relevant centre sector. Permuting
the three links then gives multiplicity exactly three for the full first
excited level. This proves Theorem YC5 without a gap-monotonicity assumption.

The classical mechanism is spectral Schur/Feshbach reduction, inverse order,
functional calculus and min–max. For standard Feshbach–Schur terminology and
reconstruction see Dusson–Sigal–Stamm,
[The Feshbach–Schur map and perturbation theory](https://arxiv.org/abs/2105.02058),
Theorem 1.2 and equations (1.10)–(1.16). That paper is background for the
method, not a source of the YC5 arithmetic constants.

## 6. Replay and limits

From the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python physics/yc5/yc5_resolved_return.py --check
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s physics/yc5 -p 'test_*.py'
```

The certificate uses only Python's standard library and exact `Fraction`
arithmetic. `explore_yc5.py` optionally uses NumPy/SciPy to propose endpoints;
it is not imported by the certificate. Generating proposals is slower than
replaying them. Source hashes pin this note, code, tests, endpoints and the
frozen predecessor inputs. The evidence is exact finite arithmetic plus the
written analytic proof, not formal proof-assistant verification.

What remains open here:

- No absolute sector splitting in the weak-coupling limit is proved.
- No increasing-volume, lattice-spacing or four-dimensional limit is taken.
- No dimensional mass scale is supplied by these dimensionless ratios.
- The present comparison loses separating margin in exploratory calculations
  near \(\theta=14\). This is not evidence of physical gap closure.

The next useful target is a second hidden-interaction enclosure, or a larger
**complete** retained harmonic cut, with its own full-complement control.
The framework's lost information is now returned with its free energy and
one actual interaction; it has not been declared recovered merely by rotating
the observer's axes.
