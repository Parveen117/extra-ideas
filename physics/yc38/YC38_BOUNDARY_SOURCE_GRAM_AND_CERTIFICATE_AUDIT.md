# YC38 — boundary source Gram and the audit of stopping rules

10 October 2026. Research owner: Monty Dabas.
Parent: `a2c89f8e541e099caaec95fb0a9fb269e7a23ed3` (YC37).

**Result.** Reuse YC31's local vacuum transport, YC32's collective-source
principle and YC37's count-complement floor. For bounded local boundary
ports on an isolated spatial block, their centered vacuum-source Gram has
a norm bound independent of block size and number of ports. Its complete
return through the original count complement therefore obeys

\[
 0\preceq {\cal R}_B(z)
 \preceq {G_B\over cb-z}
 \preceq {\min\{N_\partial,C_\Gamma\}\over cb-z}I,
 \qquad z<cb.                                             \tag{1}
\]

The constant is uniform in the inherited coupling windows, but not
numerically optimized. For regular three-dimensional blocks, at zero
energy and bounded scalar port coefficients, this improves the *same
selected source return* from a crude O(b^(1/3)) upper bound to
O(b^(-1/3)). This is a written bound on vacuum-to-boundary source vectors.
It does not bound the full retained-space interaction or every boundary
history. No larger lattice gap window or continuum gap is established.

The accompanying [audit](STOPPING_RULE_AUDIT.md) found no incorrect
rejection in the inspected YM93 and NSB2 certificate guards. It found
several restrictions of representation or sufficient estimates that must
not be read as impossibility theorems. YC37 already kept the metric and
history; the remaining opportunity is to keep correlations in its bounds.

Theory-first development: no new executable certificate, routine test
suite, numerical source solve or frozen-packet regeneration accompanies
this note. Proofs below are ordinary mathematical derivations, conditional
on the stated inherited theorems, not formal machine verification.

## 1. Full carrier and the old endpoint retained

Use YC37's isolated b-factor block, b>=8, obtained by setting the original
interblock face coefficients to zero inside the admitted YC27 family.
Retain all factor harmonics, boundary representations and multiplicities.
Internal Gauss invariance is imposed as in YC37; boundary singlets are
not imposed separately. Let h_B be its exact Hamiltonian above its actual
ground energy and psi_B its normalized ground vector.

Keep the original domain-splitting count cut and hidden operator:

\[
 P_B=1_{[0,\lfloor b/8\rfloor]}(N_B),\quad Q_B=I-P_B,
 \quad D_B=Q_Bh_BQ_B\big|_{Q_B\mathcal K_B}\ge cbI.          \tag{2}
\]

The operator and its domain are those of YC37(7)-(8), not a compression
by the rotated complement of the graph isometry J_B. At the inherited
caps, c=47/800 or 1337/22800, respectively. These are dimensionless
lattice electric units and stay bounded away from zero as b increases.

YC31 applies to this block uniformly: set all other mixed coefficients
to zero and use its full tensor-space floor 1/2. The resulting family is
a tensor sum of the isolated blocks and untouched factors. Its unique
ground and spectral flow factorize. Thus

\[
 \psi_B=U_B\Omega_B
\]

with Omega_B the product of complete reference-factor vacua, and the
same quasi-local bounds for U_B and U_B*. This is a consequence of the
admitted inhomogeneous parent family, not a claim that an arbitrary open
box inherits a bulk gap. Locality is first proved on the complete tensor
space; internal gauge restriction is made afterward.

For any bounded O supported on X, YC31(12), with exponent n=12, supplies

\[
 \|U_B^*OU_B-\mathbb E_{X(r)\cap B}(U_B^*OU_B)\|
 \le \|O\|\,|X|\,L_{12}(r+1)^{-8}.                       \tag{3}
\]

Here the expectation is in the product reference outside the indicated
set. It is a unital contraction. One may use the supremum over the old
coupling windows of the explicit YC31 constant

\[
 L_{12}={4\kappa(\exp(2C_{12}J_{12})-1)\over8C_{12}}.
\]

The constants C_12 and J_12 are exactly YC31(11)-(12), not covariance
matrices from this note. Their finite filter integrals and sums are
inherited analytic bounds. No practical value is asserted here.

## 2. Specify the port family and retain its source metric

Let O_p, p=1,...,N_partial, be bounded block operators, with norm at most
one, supported on X_p of at most k_0 original factors and uniformly
bounded diameter. They preserve the internal Gauss carrier; they may
carry boundary charge. The local factors of YC37's Wilson boundary
expansion supply such ports. Keep every matrix-index label, with a
bounded number of labels per microscopic face. Normalization factors
are carried in the coefficients rather than erased.

Require the local density bound

\[
 \#\{q:d(X_p,X_q)\le r\}\le \nu(r+1)^3.                  \tag{4}
\]

The bounded face incidence, bounded support diameter, finite label
multiplicity and YC31's cubic factor-graph growth imply such a nu,
independent of block size. For example anchor each port in a factor of
its support. If there are at most mu ports per anchor, support diameter
at most d_0, and |X_p|<=k_0, then
nu=k_0 mu kappa(d_0+1)^3 is sufficient. Distances use the original
factor graph with all possible mixed faces, including faces whose
couplings have been turned off.

Construct the actual centered source vectors, their synthesis map and
the full complex Hermitian Gram:

\[
 \omega_p=\langle\psi_B,O_p\psi_B\rangle,\quad
 v_p=(O_p-\omega_p I)\psi_B,\quad
 \mathcal V_B a=\sum_p a_pv_p,\quad G_B=\mathcal V_B^*\mathcal V_B.                  \tag{5}
\]

Thus (G_B)_pq=omega(O_p* O_q)-conj(omega_p)omega_q and
||v_p||<=1. No orthogonality of interacting ports is assumed. This
source matrix retains relative orientations and phases lost by a sum
of individual norms. It can have a kernel; a positive lower Gram bound
or inversion of G_B is not assumed.

## 3. Summable covariance follows from the inherited transport

**Lemma.** Under (3)-(4), G_B has a size-independent upper bound.

For d=d(X_p,X_q)>=1 set r=floor((d-1)/2). The r-neighborhoods of the
two supports are disjoint. In reference coordinates put
A=U_B*O_pU_B, B=U_B*O_qU_B, and replace each by its product-state
localized approximation A_r,B_r. Their norms are at most one and
each replacement error is at most k_0 L_12(r+1)^(-8).

Product-state factorization gives Cov_Omega(A_r,B_r)=0, including
nonself-adjoint ports. Expanding the covariance difference shows

\[
 |\operatorname{Cov}_\Omega(A,B)|
 \le 2\|A-A_r\|+2\|B-B_r\|
 \le4k_0L_{12}(r+1)^{-8}.                               \tag{6}
\]

Cauchy-Schwarz applied to the centered vectors also bounds every
entry by one. Define

\[
 \gamma(0)=1,\qquad
 \gamma(d)=\min\{1,4k_0L_{12}
       (\lfloor(d-1)/2\rfloor+1)^{-8}\}\quad(d\ge1),
\]
\[
 C_\Gamma=\nu\sum_{d=0}^{\infty}(d+1)^3\gamma(d)<\infty.
                                                               \tag{7}
\]

The sum converges since its tail is O(d^(-5)). Bounding each distance
shell by the ball count (4) gives both row and column sums at most
C_Gamma. The matrix Schur test, positivity and trace G_B<=N_partial
then yield

\[
 \boxed{0\preceq G_B\preceq
       \min\{N_\partial,C_\Gamma\}I.}                    \tag{8}
\]

This is the specific reuse of YC31-32: pulled-back local approximants
factor in the product vacuum, and spatial overlap is summable. It is
not an assumption of product correlations in the actual vacuum.

## 4. The complete hidden inverse, energy response and metric

Project the synthesis map (5) and retain the entire inverse in (2):

\[
 S_B=Q_B\mathcal V_B,\qquad
 {\cal R}_B(z)=S_B^*(D_B-z)^{-1}S_B,\qquad z<cb.          \tag{9}
\]

The hidden-projected source Gram need not have entrywise local decay:
Q_B is a global count cut within the block. The argument does **not**
assume that it does. Instead use the positive-operator relations

\[
 S_B^*S_B=\mathcal V_B^*Q_B\mathcal V_B\preceq G_B,\quad
 0<(D_B-z)^{-1}\preceq(cb-z)^{-1}I.                       \tag{10}
\]

They prove (1). All matrix cross terms and all hidden harmonics remain.
This is the return for the specified centered source at the original
count split. It is not the complete self-energy of the metric-isometric
retained Hamiltonian, nor the exact ground-energy second derivative of
an arbitrary interblock perturbation. Those quantities have additional
source and retained-channel contributions.

Energy response is also controlled without a finite derivative cutoff.
For every integer k>=0 and real z<cb, functional calculus gives

\[
 {\cal R}_B^{(k)}(z)
   =k!S_B^*(D_B-z)^{-k-1}S_B
   \preceq {k!G_B\over(cb-z)^{k+1}}.                     \tag{11}
\]

In particular R_B'(z) is the Gram of the lifted hidden responses
(D_B-z)^(-1)S_B. At zero its upper bound has an extra factor 1/(cb)
compared with the energy-return bound. This retained metric information
must not be replaced by identity. The same resolvent is analytic for
|z|<cb; an absolute Taylor remainder is bounded by the geometric
resolvent expansion on any fixed |z|<=z_*<cb.

These are all energy derivatives of a two-port kernel. They are not
bounds for arbitrarily many boundary insertions.

## 5. A genuine improvement, with its crossover left explicit

Take scalar coefficients a_p with |a_p|<=g. Equation (1) gives

\[
 0\le a^*{\cal R}_B(0)a
 \le {g^2N_\partial\min\{N_\partial,C_\Gamma\}\over cb}.
                                                               \tag{12}
\]

The old triangle estimate on **these same centered source vectors** is
g^2 N_partial^2/(cb). Hence (12) is never a worse upper bound. The
improvement is strict at the level of bounds when N_partial>C_Gamma;
no numerical crossover size is certified.

For geometrically regular three-dimensional blocks with
N_partial<=C_partial b^(2/3), the size-independent branch gives

\[
 \boxed{a^*{\cal R}_B(0)a
       \le {g^2 C_\Gamma C_\partial\over c}\,b^{-1/3}.}    \tag{13}
\]

The corresponding crude estimate was O(g^2 b^(1/3)). Similarly the
selected lifted-response norm squared a*R_B'(0)a is O(g^2 b^(-4/3)).
At fixed |z|<=z_* and b>=2z_*/c, the analogous bounds change only by
the indicated factors of two in the denominator. Constants, g and the
energy window are held fixed in these statements.

This removes a growing bound for one selected vacuum-source piece of
the YC37 boundary problem. It does **not** replace YC37's full-operator
boundary warning by a general shrinking bound. Neighboring blocks in
an interacting join supply operator-valued, possibly entangled sources;
the scalar coefficient estimate (12) is not a bound for them.

## 6. Transport, source changes and zero specialization

A unitary frame change U must carry h, P, Q, D, psi and every O_p
together. Then V becomes UV, S becomes US, and G and R(z) are unchanged
on the labelled coefficient carrier. This is the invariant statement.
If the source labels are recombined by a matrix L, the exact rule is

\[
 (G,{\cal R})\longmapsto(L^*GL,L^*{\cal R}L).              \tag{14}
\]

The relative inequality R<=(cb-z)^(-1)G is preserved by this rule,
including noninvertible L. The simpler Euclidean coefficient bound
G<=C_Gamma I uses the original normalized local-port frame; arbitrary
nonunitary mixing must carry its metric or its norm cost. Spatial
locality is transported with the observable algebra, not inferred for
an arbitrary rotation in a fixed original locality chart.

At mixed coupling zero, U_B=I, psi_B=Omega_B and the count projection
commutes with H0,B. Equations (5),(9) recover the actual product-reference
source and complete reference hidden inverse. Disjoint source supports
then have zero covariance exactly. If the fixed port support cap obeys
k_0<=floor(b/8), each centered port acting on Omega_B excites at most
k_0 factors, so Q_B v_p=0 and R_B(z)=0 exactly. If this count condition
does not hold, retain the nonzero reference response rather than forcing
it to vanish. Boundary operator and coupling derivatives remain data.

Centering (5) is not a change to the physical Hamiltonian. In an actual
join write each O_p=omega_p I+delta O_p and retain its scalar, one-sided
and mixed pieces, including their source derivatives. Discarding the
mean terms would not be a legitimate application of this theorem.

The radius/size b is not a coordinate rotation and not a continuum
spacing. It physically changes the retained count cap. Under a positive
energy-unit change h_phys=e h, dimensionless ports give
R_phys(ez)=e^(-1)R(z), while physical coefficients a_phys=e a make
a_phys*R_phys a_phys=e a*R a. Relative and absolute energy claims thus
retain the calibration instead of normalizing it away.

## 7. Exact controls and the next unresolved estimate

Two elementary examples check why the source scope matters.

1. On b independent two-level factors, h=sum_i |1><1|_i, use the
   same count cut and O_p=X_p on distinct factors. Then G=I exactly.
   If floor(b/8)>=1, every centered source has count one and S=0.
   A norm bound of order N_partial^2/b misses this exact zero; the
   source rule preserves it. This is a mathematical control, not a
   simulation of a Yang-Mills block.
2. On the same product ground use O_p=Z_p. Every centered source
   vanishes and G=0, but sum_p Z_p is a nonzero operator of norm
   N_partial. Therefore even an exact vacuum-source bound cannot be
   promoted to a full retained-space operator bound. This counterexample
   blocks precisely the overclaim that (13) closes the entire iteration.

Both conclusions follow directly from X|0>=|1>, Z|0>=|0>; no numerical
test or code execution is being reported for these controls.

The next specific task is an operator-valued connected-history estimate
for YC37(23)-(24), controlling arbitrary retained inputs, neighboring
block entanglement, internal residuals and all word lengths in one
spatial norm. The present two-port vacuum Gram and all its energy
derivatives are a bounded part of that problem. YC36's constants still
depend on source count, and no all-word summation or repeated-block
contraction follows from (13).

**Established here:** a transported source-Gram rule, its uniform upper
bound from inherited locality, and decay of the selected original-cut
return at regular block size. **Still open:** practical constants,
computable vacuum ports, the all-input/all-history estimate, the next
scale's admissible interaction bounds, and physical continuum recovery.

## Sources and evidence boundary

Exact file digests and revisions are in [SOURCES.json](SOURCES.json).
The argument reuses YC31(12), YC32's source-overlap principle and
YC37(8), rather than restarting the Hamiltonian construction. Earlier
theorem packets and certificate evidence are unchanged.

The inherited locality framework is Nachtergaele--Sims--Young,
[arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2), including
unbounded onsite terms and quasi-local maps. Schur/Feshbach background:
Dusson--Sigal--Stamm,
[arXiv:2105.02058](https://arxiv.org/abs/2105.02058).
These references support the standard methods, not a priority claim or
an independent verification of this programme's Yang-Mills inputs.
