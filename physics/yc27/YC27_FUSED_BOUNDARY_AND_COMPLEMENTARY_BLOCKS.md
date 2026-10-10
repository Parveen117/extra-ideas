# YC27 — fuse the boundary charge, retain the complementary loops

10 October 2026. Research owner: Monty Dabas. Continues YC26 at `842e86f`.
Unit-S³ electric normalization, normalized Haar measure, Wp=Tr(Up)/2.
Predecessor proof and evidence packets remain unchanged.

**Result.** Several exterior links at a vertex carry a combined charge,
which can cancel. Their remaining kinetic and loop energies must still
be retained. For the existing rectangular tilings the exterior-link graph
has an exact, useful structure: its connected components are small cubes,
squares and individual links. Absorb their complete internal Hamiltonians
into the reference, including every charged and neutral excitation.

The new mixed-source estimate then gives actual, volume-uniform lattice
bounds on the following anisotropic profiles:

| Old reference block | Old internal caps | Complementary loop caps | Remaining mixed cap | Physical lattice gap lower bound |
|---|---|---:|---:|---:|
| 28 links, vertex shape 4×2×2 | cube theta<=1, x-join<=1 | **1/2** | **1/3200** | **592/175 > 3.3828** |
| 64 links, vertex shape 4×4×2 | cube theta<=1, x-join<=1, y-join<=1/2 | **1/2** | **1/5700** | **11188/3325 > 3.3648** |

Every complementary and mixed plaquette may have its own nonnegative
coupling below the stated cap. The old blocks use exactly YC26's admitted
internal profiles. The 28-link tori have Lx divisible by4 and at least8,
Ly,Lz even and at least4. The 64-link tori have Lx,Ly divisible by4 and
at least8, with Lz even and at least4.

Previously all exterior faces belonged to the weak perturbation. Here
the faces made entirely of old exterior links may be as strong as1/2.
The remaining mixed caps are smaller than YC26's old exterior cap. These
are different coupling profiles, not a uniformly enlarged window or an
improved isotropic theorem. The second row is a full spatial join using
YC26's actual64-link blocks. Repeated scale contraction and the4D
continuum gap are not established.

## 1. The full family and its two cuts

Let O denote the old blocks and D the connected components of their
complementary link graph, described in section2. They partition the links,
not the vertices. The complete Hilbert carrier is

\[
 \mathcal H=\bigotimes_{f\in O\cup D}L^2(SU(2)^{E_f},d\mu_f).
                                                               \tag{1}
\]

Keep each full factor space. Gauge invariance at a shared vertex acts
jointly on its old and complementary factors; it is not a restriction
to independent factor singlets. All electric domains are the inherited
finite-product Laplacian domains, and all potentials are bounded smooth
Wilson multipliers. The reference pairing is the full Haar pairing.

Construct the raw operator family first:

\[
 \widetilde H(u,\lambda)=\sum_{o\in O}H_o
 +\sum_{d\in D}\left[-\sum_{e\in E_d}\Delta_e
                         -\sum_{p\subset d}u_pW_p\right]
 -\lambda\sum_{p\in\mathcal P_{\rm mix}}\xi_pW_p,
 \quad 0\le u_p\le\tfrac12,\quad0\le\xi_p\le1.        \tag{2}
\]

The operator carrier includes composition, the gauge action, factor cuts
and their complete returns on these domains. The coordinates u and lambda
are an explicit YM adapter, not an identification of the universal
primitive operator space with scalar coupling space. The old couplings
are fixed in the table's windows.

Every factor has a unique positive normalized ground Omega_f, invariant
under its vertex gauge action. Put h_f=H_f-E_f and carry the scalar sum
of E_f separately. The shifted family becomes

\[
 H_\lambda(u)=H_0(u)+\lambda\Phi,\qquad
 H_0(u)=\sum_fh_f,\quad\Phi=-\sum_{p\in\mathcal P_{\rm mix}}\xi_pW_p,
 \qquad\Omega(u)=\bigotimes_f\Omega_f.                 \tag{3}
\]

At lambda=0 this is an exact product of correlated old and complementary
blocks. At u=0 each complementary factor is its full free-link product;
regrouping its tensor factors gives exactly the old reference with
individual Haar bridges. These evaluations commute on the raw family(2).
Scalar shifts must be restored when comparing unshifted operators.
Section7 states the corresponding source, state and metric compatibility.

## 2. What multi-link fusion leaves behind

At an old boundary vertex v, physical Gauss matching reads

\[
 \left(G_{o,v}+\sum_{e\ {m exterior\ at}\ v}J_{e,v}\right)\psi=0.
\]

Thus its block Casimir matches

\[
 C_{o,v}\psi=C_{{\rm ext},v}\psi,\qquad
 C_{{\rm ext},v}=-\left(\sum_eJ_{e,v}\right)^2.         \tag{4}
\]

It does not match the sum of the separate link Casimirs. For example,
two fundamental endpoint charges fuse to doubled charge q=0 or q=2,
with combined Casimir0 or8. Their electric energies still sum to6.
The q=0 possibility disproves summing YC26's one-bridge coefficient
independently over every exterior link at the same vertex.

Before regrouping, YC22's additive ground-form inequality gives the
valid replacement

\[
 \sum_oh_o+\sum_{e\ {\rm exterior}}(-\Delta_e)
 \succeq\sum_e(-\Delta_e)
       +\frac12\sum_{o,v}\frac{C_{{\rm ext},v}}{d_o(v)}.\tag{5}
\]

The neutral external loops in(5) are not an extensive unknown graph here.
For vertex side lengths a_i>=2, the exterior links in coordinate i join
the disjoint pairs {k a_i-1,k a_i}; other coordinate values are singletons.
Their full graph is the Cartesian product of these one-dimensional
components. Each connected component is therefore K2 raised to q factors,
with q=1,2,3, or an isolated vertex with q=0 and no link factor.

For our two shapes a_z=2, every vertex belongs to a nonempty component:

| Component | Vertices | Links | Internal plaquettes | Mixed-face incidence |
|---|---:|---:|---:|---:|
| Cube | 8 | 12 | 6 | 24 |
| Square | 4 | 4 | 1 | 12 |
| Link | 2 | 1 | 0 | 4 |

The 4×2×2 tiling has only cubes and squares. The 4×4×2 tiling also has
links. An old block and a complementary component share at most one
vertex: choosing either endpoint of a boundary pair chooses a different
old interval in that coordinate. At least two old intervals per periodic
coordinate ensure these are distinct across the periodic boundary too.

Make an incidence graph whose nodes are the old and complementary
factors and whose edges are their shared physical vertices. This graph
is **simple and bipartite**. At every shared vertex the exact degrees obey

\[
 \boxed{d_o(v)+d_d(v)=6.}                              \tag{6}
\]

A plaquette has either zero, one or two boundary-crossing coordinate
directions. It is respectively wholly old, mixed, or wholly complementary.
A mixed face contains one edge from each of exactly two old factors and
two distinct complementary factors. Its support size is4, even though
all four factors may now be internally correlated.

There are two mixed faces per boundary-surface internal edge of an old
block. Thus the old mixed incidences are 2(4+10+10)=48 for4×2×2, and
2(10+10+24)=88 for4×4×2. Complementary incidences follow because each
link has four lattice plaquettes: a cube link has two internal faces,
a square link has one, and a standalone link has none. This proves the
table for all admitted volumes, without extrapolating a finite drawing.

## 3. Neutral complementary gaps and a matched physical reference

On a free cube or square, a nonconstant gauge-invariant spin network
needs a cycle. A nontrivial edge cannot end alone at a gauge-invariant
vertex, and the shortest cycle has four edges. Each nontrivial edge
costs at least3, so the free physical gap is12. A fundamental loop attains
it. The gauge-invariant carrier of a single open link contains only constants.

For arbitrary complementary plaquette couplings at most u, the potential
oscillation is at most2 sum_p u_p. The min-max bound on both levels gives

\[
 \delta_{D,\rm cube}\ge12-12u\ge6,\qquad
 \delta_{D,\rm square}\ge12-2u\ge11\quad(u\le\tfrac12).
                                                               \tag{7}
\]

These gaps are above the true factor ground. Charged sectors are also
retained: YC22 gives their separate floors1,3/2,3 for a cube, square,
and link. YC26 gives the old neutral excited floors9/2 and4 in our two
rows. No charged factor is discarded when using these neutral bounds.

**F1 — extensive physical energy and its first level.** Let
N=sum_f(1-|Omega_f><Omega_f|) count excited factors in the new partition.
On the complete globally gauge-invariant carrier,

\[
 \boxed{H_0\succeq N,\qquad
       H_0|_{\Omega^\perp\cap\mathcal H_G}\succeq4I.} \tag{8}
\]

**Proof.** Decompose into complete vertex-charge isotypic carriers for
each factor. Their Casimirs commute with H_0 and the excitation cuts.
At a shared vertex the charges of the two factors must match; write
their doubled spin q_v and kappa_v=q_v(q_v+2). All multiplicities and
all neutral excitations inside each carrier remain.

YC22's additive Casimir bound and(6) give

\[
 H_{0,\rm charged}\succeq
 \frac12\sum_v\left(\frac1{d_o(v)}+\frac1{d_d(v)}\right)\kappa_v
 \succeq\frac13\sum_v\kappa_v,                       \tag{9}
\]

where d_o d_d<=9. The last coefficient is the gain from keeping the
two endpoint degrees together, rather than bounding each by an unrelated
maximum.

Let n_O,n_D count charged old and complementary factors and S=sum_v kappa_v.
Simultaneously applying the centre at every vertex of a pure-link factor
is the identity. A charged factor therefore has either at least two odd
q_v, costing at least6, or a nonzero even q_v, costing at least8. Since
each shared vertex is incident once on each side of the incidence graph,
S>=6n_O and S>=6n_D. Equation(9) is at least2 max(n_O,n_D), hence at
least n_O+n_D. Every neutral excited factor has its independent floor
at least4 by(7) and YC26. Adding them proves H_0>=N on every physical
support, then by form closure on the whole physical carrier.

For the first level, if an odd q_v occurs, odd-charge edges form a
nonempty even-degree subgraph of the simple bipartite incidence graph.
It contains a cycle of length at least4. Therefore S>=12 and(9) costs
at least4. If only even nonzero charges occur, choose one shared vertex.
Apply YC22's separate single-vertex form to its two different factors:

\[
 h_o+h_d\succeq\kappa_v\left(1/d_o+1/d_d\right)
                    \ge8\cdot\tfrac23=\tfrac{16}3>4.\tag{10}
\]

If no charged factor occurs, any excitation costs at least4. The neutral
ground product is the only zero vector. This proves the second assertion.
The use of a single-vertex bound in(10) does not sum that bound over
different vertices of the same factor. QED.

The statement includes all global centre characters of the physical
carrier. A mixed plaquette source has four excited factors and therefore
has the complete source floor4 as well. This is an all-harmonic bound,
not a restriction to fundamental intermediate states.

## 4. The mixed source with four correlated factors

For a mixed face p let r_p=Wp Omega. The single-edge marginal of each
gauge-invariant factor ground density is Haar. Since the four selected
edges belong to independent factors, the normalized trace calculation gives

\[
 \|r_p\|^2=\tfrac14,\qquad
 \langle\Omega,r_p\rangle=0,\qquad
 \langle r_p,H_0r_p\rangle=3.                         \tag{11}
\]

In fact r_p is excited in each of its four factors. Define
T_(f,e)=integral |grad_e Omega_f|² dmu_f. The full interacting residual
s_p=(H_0-12)r_p satisfies the exact identity

\[
 \boxed{\langle r_p,s_p\rangle=0,\qquad
        \|s_p\|^2=\sum_{(f,e)\ {m in}\ p}T_{f,e}.}  \tag{12}
\]

**Proof of the four-factor extension.** In one factor use the four
quaternion coordinates u_alpha of its selected oriented edge and put
xi_alpha=u_alpha Omega_f, eta_alpha=(h_f-3)xi_alpha. The endpoint
SU(2)×SU(2) representation on the coordinate indices is irreducible,
so each of the three Gram matrices below is scalar. The ground-form
identity and sum_alpha u_alpha²=1 imply

\[
 \langle\xi_\alpha,\xi_\beta\rangle=\delta_{\alpha\beta}/4,
 \quad\langle\xi_\alpha,\eta_\beta\rangle=0,
 \quad\langle\eta_\alpha,\eta_\beta\rangle=T_{f,e}\delta_{\alpha\beta}.
                                                               \tag{13}
\]

For the second matrix, the trace of <xi,h xi> is3. For the third, the
product rule gives eta_alpha=-2 grad u_alpha dot grad Omega_f, and the
spherical coordinate identity sums its norm squares to4 T_(f,e).
Expand Wp as a multilinear tensor in the four edges' quaternion
coordinates. Replacing one xi Gram by its eta Gram multiplies that
factor's contribution by4 T_(f,e), so its residual norm is T_(f,e)
because ||r_p||²=1/4. Cross terms with two different residual factors
vanish by the middle identity in(13). This proves(12).

No four-dimensional source span is declared invariant under h_f.
The eta vectors and the inverse can enter arbitrarily many copies of
the same endpoint representation. Equations(11)–(13) constrain their
moments without truncating that carrier.

For x>=d>0 the exact scalar identity

\[
 x^{-2}-12^{-2}+2(x-12)/12^3
 =(x-12)^2\left(\frac1{864x}+\frac1{144x^2}\right)
\]

and(11),(12) give, with the full source floor d=4,

\[
 \|H_0^{-1}r_p\|^2\le\frac1{576}
   +\|s_p\|^2\left(\frac1{864d}+\frac1{144d^2}\right).\tag{14}
\]

Thus the source mean12 is not substituted for the excitation operator.

## 5. Local forces and explicit source budgets

YC22's full-ground force bound is

\[
 T_{f,e}\le\frac{\min(d_f(v),d_f(w))^2}{64}
              \left(\sum_{p\subset f,\ p\ni e}|t_p|\right)^2.
                                                               \tag{15}
\]

For a complementary cube, its endpoint degrees are3 and its two
incident internal couplings sum to at most1. Thus T<=9/64. A square
has T<=1/64, and a link has T=0. The common9/64 cap suffices.

For old edges which actually occur in mixed faces, the local geometry
and the table's coupling caps give:

| Old shape | Edge parallel to x | Parallel to y | Parallel to z |
|---|---:|---:|---:|
| 4×2×2 | 1 | 9/4 | 9/4 |
| 4×4×2 | 625/256 | 9/4 | 9/4 |

Here is the local count behind the second row. An x edge has at most
three incident old faces. At an interior y value their caps sum to
1+1/2+1=5/2, and its endpoint degrees are at most5, giving625/256.
At a boundary y value the cap sum is at most2 and the degrees at most4.
A y edge inside a parent rectangle has cap sum at most3 and minimum
endpoint degree at most4; a y-joining edge has sum at most3/2 and
degrees at most5. Both bounds are at most9/4. A z edge in a mixed
face must lie on an x or y boundary, so its degrees are at most4 and
cap sum at most3. For the first row an x edge has cap sum at most2
and degrees at most4; y and z edges have sums at most3 and degrees
at most4. This supplies the table for every admitted tiling.

Each mixed face uses two old and two complementary factors. Consequently

\[
 v_{28}:=\|s_p\|^2\le\frac{153}{32},\qquad
 v_{64}:=\|s_p\|^2\le\frac{661}{128}.                 \tag{16}
\]

Equation(14) gives respectively

\[
 \|H_0^{-1}r_p\|^2\le\frac{383}{73728},\qquad
 \|H_0^{-1}r_p\|^2\le\frac{4841}{884736}.
\]

Both are strictly below (3/40)². With the support-size factor4, the
creator root norm therefore has first-source budget

\[
 s_0\le4I\frac3{40}|\lambda|,
 \qquad I=48\ \hbox{or}\ 88.                         \tag{17}
\]

The incidence bound also covers every complementary factor by section2.

## 6. Complete creator, physical resolvent and lattice windows

Use YC25's full invariant creator construction with this new partition.
For exact excited-support projections P_J, let

\[
 \widehat c_J=|c_J\rangle\langle\Omega_J|,\quad
 C=\sum_{J\ne\varnothing}\widehat c_J,\quad
 \|c\|_* =\max_f\sum_{J\ni f}|J|\|c_J\|,
\]
\[
 c_J=-\lambda H_J^{-1}P_Je^{-C}\Phi e^C\Omega.         \tag{18}
\]

The family acts on gauge-invariant creator targets. Its vacuum branch
also preserves global centre symmetry. Every projection, source word,
creator and inverse preserves the required gauge action. Equation(8)
therefore supplies ||H_J^-1||<=1/|J| on every nonempty physical support,
not merely the first source. The minimum nonlinear support is still1:
neutral single-factor excitations are allowed.

YC25's full root and excitation-quotient bounds now have mu=1, k=4.
Writing beta=I|lambda|, for a creator ball of radius r they are

\[
 q\le4\beta e^{8r}[2+8(1+2r)],\qquad
 \|\mathcal T(c)\|_*\le s_0+qr,
 \qquad b_G\le8\beta e^{8r}.                         \tag{19}
\]

The physical reference resolvent has first-level floor4 by(8), separate
from its extensive coefficient1. The relative Neumann argument on the
complete physical excitation quotient gives

\[
 \boxed{\Delta_G(H_\lambda)\ge4(1-b_G).}               \tag{20}
\]

This is the full physical quotient, including all global centre sectors.
As in YC25, the compact elliptic Hamiltonian's true positive ground is
gauge invariant. The constructed branch has no lower physical level by
the same resolvent argument, hence is this true ground. Independently
YC22 gives the nonphysical charged complement floor1/2 on the full torus.
If required, the complete unrestricted gap is at least1/2 on these windows.

Take r=1/32 and use exp(1/4)<=12491/9728<9/7. All endpoint budgets are
exact rational inequalities:

| Old block | Mixed cap | q upper | b_G upper | s0 upper | r-s0-qr lower | Physical gap lower |
|---|---:|---:|---:|---:|---:|---:|
| 28 links | 1/3200 | 81/100 | 27/175 | 9/2000 | 23/16000 | 592/175 |
| 64 links | 1/5700 | 396/475 | 528/3325 | 11/2375 | 43/76000 | 11188/3325 |

The ball reserves are positive and q<1. Monotonicity in |lambda| includes
every smaller mixed coupling and every allowed profile. Constants are
uniform in the admitted finite volumes and complementary couplings.
All support sizes, harmonic levels, charged exchanges and neutral
multiplicity mixing remain in(18). No extra harmonic cutoff enters(20).

## 7. Full-family and metric compatibility

For each fixed admitted real u, the uniform contraction constructs the
whole holomorphic creator family in the complex lambda disk. Its strict
endpoint reserves leave a small open neighborhood too. This is not a
finite series definition or a claim of an infinite-volume ground-state
construction. Finite volume has finitely many factor supports, with
complete infinite-dimensional carriers. The solutions belong to Dom H_J;
YC9/YC25's bounded-creator domain argument applies without change.

Let S_lambda=exp C_lambda. At the mixed-interface cut lambda=0,

\[
 C_0=0,\quad S_0=I,\quad\Psi_0=\Omega(u),\quad E_0=0,
 \quad A_0=H_0(u)|_Q,\quad Q=I-|\Omega\rangle\langle\Omega|.
                                                               \tag{21}
\]

Energy is measured relative to the carried reference scalar. The source
history survives as

\[
 c'_J(0)=-H_J^{-1}P_J\Phi\Omega,\qquad E'(0)=0.        \tag{22}
\]

Define the full excitation quotient by
A_lambda=Q[S_lambda^-1(H_lambda-E_lambda)S_lambda]Q on Q's carrier.
The vacuum is invariant under the dressed operator, so this compression
represents the quotient by that vacuum. For real lambda retain the full
similarity pairing and its excitation quotient metric,

\[
 M=S_\lambda^\dagger S_\lambda,\qquad
 G=M_{QQ}-M_{Q0}M_{00}^{-1}M_{0Q},\qquad A^\dagger G=GA.
                                                               \tag{23}
\]

At zero, M=I and G=I_Q. Equations(2),(18),(21)–(23) prove the typed
evaluation identity pi_0(Construct_lambda)=Construct_0 for the operator,
state, source return and metric in this representation. Internal old and
complementary loop interactions remain at this zero-reading.

At u=0 the original Hamiltonian, reference ground and complete physical
source inverse additionally recover the old free-bridge model under
tensor regrouping. Creator support coordinates change when individual
links are regrouped into factors; no equality of their nonzero-lambda
coordinate metrics is asserted without that coordinate map. At lambda=0
both readings have the original Haar pairing and identity metric. No
volume-uniform condition number of S or G is asserted.

## 8. What this resolves and what the next step needs

Fusion can remove a boundary charge while leaving a neutral loop. This
stage retains those loops in their exact finite connected components and
puts their full interactions into the reference. Their couplings no longer
need the small mixed-face cap. The gain is an exact new factorization and
two all-volume anisotropic profiles, rather than a declaration that every
neutral return has vanished.

The remaining source is specific: each mixed plaquette crosses two old
and two complementary correlated factors. Its first inverse is controlled
above; its later connected neutral returns remain in the complete creator
family but are still bounded by the general norm count. The present
incidence/extensive-energy budgets are48 and88, so these two estimates do
not constitute a contracting spatial iteration. Comparing factor counts
from different partitions would not establish one either.

The next useful development is a sharper connected return for these
mixed four-factor sources, retaining fused charges and the physical
metric, to improve the cost after another partition update. Alternating
two partitions alone supplies no RG law. Physical scale calibration,
running along a continuum trajectory and continuum observables remain
additional obligations.

## 9. Minimal verification and source route

```sh
python physics/yc27/yc27_complementary_join.py
```

The calculator enumerates only the two smallest admissible tori to check
the new partition, degree pairing, face supports, incidences and boundary
forces. It then reproduces the new exact source and endpoint budgets.
The all-volume geometry is proved in section2, the full support bound
in section3, and the all-order continuation in sections6–7. These finite
checks are not a lattice diagonalization or independent formal verification.
No new unit-test suite, frozen evidence regeneration or predecessor
regression run is used.

Inputs at `842e86f`: YC22 sections1–4 for the true-ground Casimir and
force forms; YC26 sections3–5 for the actual old block gap windows;
YC25 sections3–5, with its YC9/YC10/YC19 references, for the complete
creator/resolvent construction, domain control and quotient metric.
The complementary geometry, degree-paired support estimate and four-
correlated-factor source identity are derived explicitly here. This
written proof has not been independently peer reviewed.
