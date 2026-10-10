# YC8 — local vacuum dressing with volume-independent error density

9 October 2026. Continues YC7 at `baca59a`. The owner asks us to retain TVSP
as a source of structure and proceed toward the many-cell Yang–Mills
problem. The construction below accounts for the common vacuum energy by
a positive local reading, before estimating the excitation problem.

**Actual Yang–Mills result:** on every periodic rectangular three-dimensional
SU(2) lattice with each side length at least 3, the ground-energy density
has the uniform enclosure

\[
\boxed{
\theta-\frac{\theta^2}{48}-c_3\theta^3-c_4\theta^4
\ \le\ \frac{E_0(\theta)}{P}\ \le\
\theta-\frac{\theta^2}{48}+c_3\theta^3,
}\tag{1}
\]

for every θ≥0, where P is the number of plaquettes and

\[
c_3=\frac{205}{5616},\qquad
c_4=\frac{42025}{14017536}.
\]

The exact local-energy identity below has no omitted higher-order terms.
The bounds become coarse at large coupling. They may be intersected with
0≤E0/P≤θ. They do not assert existence of a thermodynamic or continuum
limit, but their constants are independent of all three side lengths.

A related **comparison Hamiltonian**, whose positive ground is explicitly
known, has a volume-independent gap ≥1/4 for 0≤θ≤1/2. The actual YM
operator differs from it, after a scalar shift, by an explicit local
residual starting at order θ³. That residual is not discarded; no
volume-uniform gap for the original operator is claimed.

## 1. Actual carrier and the local source

Let the periodic box be N_x×N_y×N_z with every N_i≥3. It has
n=N_x N_y N_z sites and E=P=3n links and plaquettes. Use the product of
unit S³ metrics and normalized SU(2) Haar measures on all links, with

\[
H_\theta=H_0+\theta V,\qquad
H_0=-\sum_e\Delta_e,\qquad
V=\sum_p(1-W_p)=P-S,\qquad S=\sum_pW_p,
\quad W_p=\tfrac12\operatorname{Tr}U_p.                     \tag{2}
\]

All link kinetic energies are kept. Local gauge invariance is imposed at
every site, and the vacuum sector has all three even global centre
characters. The positive full-space ground of (2) belongs to this sector.
All positive trial readings below are also gauge invariant and centre
even. Proving comparison curvature on the unreduced product avoids any
assumption about the regularity of a gauge quotient.

Unlike YC7's 2×1×1 geometry, each plaquette here contains four distinct
links, and distinct plaquettes share at most one link. Each link belongs
to four plaquettes. Thus

\[
H_0W_p=12W_p,\qquad
\langle W_p\rangle_\mu=0,\quad
\langle W_p^2\rangle_\mu=\tfrac14,\quad
\langle W_pW_q\rangle_\mu=0\ (p\ne q).                     \tag{3}
\]

Proof: a Wilson reading is linear in each of its four quaternion links,
whose degree-one electric energy is 3. A unique link in either of two
different plaquettes makes their product Haar mean zero. One Haar link
makes the complete holonomy Haar-distributed, giving the variance 1/4.
Consequently the constant reading's entire source has energy 12, square
P/4 and second-order return P/48. Equation (1) upgrades this coefficient
to a uniform finite-volume error statement, not just fixed-volume formal
perturbation theory.

The side-length condition matters: repeated-link commutator plaquettes
and the short winding channels in YC7 have different rates and coefficients.
YC8 does not overwrite that frozen packet or apply its ordinary-plaquette
formula to a collapsed torus.

## 2. Shared-edge algebra retains every local cross source

Write

\[
\Gamma(f,g)=\sum_e\langle\nabla_e f,\nabla_e g\rangle,
\qquad\Gamma(f)=\Gamma(f,f).
\]

The sign convention is H0=−Δ, hence
H0(fg)=fH0g+gH0f−2Γ(f,g).
For one plaquette,

\[
\Gamma(W_p)=4(1-W_p^2)=3-C_p,\quad
C_p=4W_p^2-1,\quad H_0C_p=32C_p.                           \tag{4}
\]

Let \(p\sim q\) mean an unordered pair sharing an edge e, with each pair
counted exactly once. In that edge's quaternion coordinate u∈S³ write

\[
W_p=a\cdot u,\quad W_q=b\cdot u,\quad |a|=|b|=1,
\quad A_{pq}=a\cdot b,\quad B_{pq}=W_pW_q.
\]

Signs from traversal of an inverse link are included in a and b. They are
not reset by the observer convention. The six other links are distinct,
and A is a normalized trace of the joined six-link word. Direct spherical
differentiation gives

\[
\Gamma_e(W_p,W_q)=A_{pq}-B_{pq},\qquad
H_0A_{pq}=18A_{pq},\quad
H_0B_{pq}=26B_{pq}-2A_{pq}.                                 \tag{5}
\]

Indeed the shared-edge quadratic has H_e B=8B−2A and the six other
fundamental links contribute 18B. Therefore

\[
\boxed{H_0^{-1}(A_{pq}-B_{pq})=
\frac{2}{39}A_{pq}-\frac1{26}B_{pq}.}                        \tag{6}
\]

Here the inverse is on the mean-zero finite source. Its two harmonic
components are (3/4)A at energy 18 and A/4−B at energy 26. Equation (6)
retains the orientation-dependent cross term; it is not a sum of separate
plaquette norms. Equations (4)–(5) give the full pointwise identity

\[
\Gamma(S)=3P-\sum_pC_p+2\sum_{p\sim q}(A_{pq}-B_{pq}).       \tag{7}
\]

## 3. A positive reading cancels the first two nonconstant orders

Define two real gauge-invariant functions, both with zero Haar mean,

\[
f_1=-\frac1{12}\sum_pW_p,
\qquad
\boxed{f_2=\frac1{4608}\sum_pC_p
-\frac1{1404}\sum_{p\sim q}A_{pq}
+\frac1{1872}\sum_{p\sim q}B_{pq}.}                          \tag{8}
\]

By (3)–(7),

\[
H_0f_1=V-P,\qquad
H_0f_2=\frac{P}{48}-\Gamma(f_1).                            \tag{9}
\]

Set Fθ=θf1+θ²f2 and ψθ=exp(−Fθ)>0, with any convenient positive
normalization. This is an explicitly constructed trial/dressed reading;
it is not identified with the actual YM ground. Direct differentiation
gives the **exact**, all-coupling local energy

\[
\boxed{\frac{H_\theta\psi_\theta}{\psi_\theta}
=P\left(\theta-\frac{\theta^2}{48}\right)+R_\theta,
\quad
R_\theta=-2\theta^3\Gamma(f_1,f_2)-\theta^4\Gamma(f_2).}
\tag{10}
\]

The exponential contains all powers; its local energy has exactly the
displayed cubic and quartic residual because Fθ is quadratic in θ.
There is no assumed convergence of an infinite perturbation series.

In the framework's language the first reading returns the single-plaquette
source; the second returns the complete adjacent-pair source generated by
the first. This is a concrete source/response construction. The degree
here is coupling degree in θ, not an identification of θ with UP1's λ or
its coordinate dilation.

## 4. Bounded local incidence removes the volume-square cost

For any one link e, precisely four plaquettes contain it. A plaquette has
12 distinct edge-sharing neighbours. Among pairs p∼q whose union contains
e, there are 6 sharing e itself and 36 with e in just one member, for a
total of 42. This is a local count, valid for every allowed box size:
4×12 counts pair incidences, and the 6 common-e pairs are counted twice,
giving 48−6=42 distinct pairs and 48−12=36 noncommon incidences.

On a unit S³ link,

\[
|\nabla_eW_p|\le1,\quad |\nabla_e C_p|\le4,\quad
|\nabla_e A_{pq}|\le1\text{ on its six links},
\]

and A is independent of the shared link. For B, the derivative norm is at
most 1 on a nonshared link and 2 on the shared link. Thus

\[
\boxed{|\nabla_e f_1|\le\frac13,\qquad
|\nabla_e f_2|\le
\frac{16}{4608}+\frac{36}{1404}+\frac{48}{1872}
=k:=\frac{205}{3744}.}                                     \tag{11}
\]

Consequently Rθ is a sum of local terms r_e with

\[
|r_e|\le\varepsilon(\theta):=\frac{2k}{3}\theta^3+k^2\theta^4,
\qquad |R_\theta|\le P\varepsilon(\theta).                 \tag{12}
\]

There is also an explicit local interaction norm, not just a density
statement. Let D(e) be the union of links in plaquettes containing e and
their edge-sharing neighbours. Then r_e depends only on D(e). At most
4+36=40 faces occur in D(e); adding them along shared edges to the initial
edge adds at most three links per face. Hence |D(e)|≤121. The relation
f∈D(e) is symmetric in e,f, so each link occurs in at most 121 such local
terms. Therefore

\[
\boxed{\sup_f\sum_{e:\,f\in D(e)}\|r_e\|_\infty
\le121\varepsilon(\theta),}                                \tag{13}
\]

independently of volume. The supports have bounded lattice range. The
finite geometry controls find smaller supports (51, 59 and 63 links on
the listed boxes); 121 is the conservative proved all-size bound.

**Proof of (1).** For a positive trial ψ and positive actual ground g,
self-adjointness gives
E0=∫gψ(Hψ/ψ)dμ / ∫gψdμ. It lies between the minimum and maximum local
energy. Apply (10)–(12); for the upper bound discard only the nonpositive
quartic term, retaining the two-sided cubic bound. This is the classical
Barta/local-energy argument, reproduced here. A primary reference is
Amaury Mouchet, *Upper and lower bounds for an eigenvalue associated with
a positive eigenvector*, theorem 1, https://arxiv.org/abs/math/0505541.

For illustration, outward decimal enclosures valid for **every** allowed
volume are:

| θ | Actual YM E0/P lower | Actual YM E0/P upper |
|---:|---:|---:|
| 1/10 | 0.0997548 | 0.0998282 |
| 1/4 | 0.2481158 | 0.2492683 |
| 1/2 | 0.4900414 | 0.4993546 |

The replay stores exact fractions. These are vacuum-energy densities in
the declared Hamiltonian units, not masses or excitation-gap estimates.

## 5. Information curvature gives a uniform gap for the explicit comparison

Define

\[
K_\theta=H_\theta-P(\theta-\theta^2/48)-R_\theta
=H_0+\frac{\Delta\psi_\theta}{\psi_\theta}.                 \tag{14}
\]

Kθ has ground ψθ at energy zero and is nonnegative. Its ground transform
is the weighted diffusion with probability measure
dνθ=ψθ²dμ/∫ψθ²dμ. This is the **comparison measure**, not the unknown
actual YM vacuum measure at nonzero θ. For this explicit operator,

\[
\langle\psi f,K_\theta\psi f\rangle
=\int|\nabla f|^2\psi^2d\mu.
\]

We now control its information curvature without any factor of P.
The following bounds use operator norms of Hessian blocks in the product
metric:

| Function | Same-link Hessian bound | Mixed-link block bound | Support |
|---|---:|---:|---:|
| W_p | 1 | 1 | 4 links |
| C_p=4W_p²−1 | 8 | 16 | 4 links |
| A_pq | 1 | 1 | 6 links |
| B_pq=W_pW_q | 4 | 4 | 7 links |

For W or A, linearity in one unit quaternion gives Hess_ee W=−W g_e;
in two distinct links, differentiating the unit-quaternion product has
norm at most 1. For C, differentiation gives
8 dW⊗dW+8W Hess W: the same-link eigenvalues have absolute value at
most 8, and a mixed block is at most 16. For B the product rule has at
most four contributions of norm at most 1. These prove the table.

The maximum absolute block-row sum bounds the full self-adjoint Hessian
operator norm. Using the same incidence counts as (11),

\[
\|\operatorname{Hess}f_1\|\le\frac{4\cdot4}{12}=\frac43,
\]
\[
\boxed{\|\operatorname{Hess}f_2\|\le
\frac{4(8+3\cdot16)}{4608}
+\frac{36\cdot6}{1404}
+\frac{42\cdot7\cdot4}{1872}
=\frac{1555}{1872}.}                                       \tag{15}
\]

Each unit S³ has Ricci tensor 2g. Therefore the product weighted curvature
is bounded by

\[
\operatorname{Ric}+2\operatorname{Hess}F_\theta
\succeq\rho(\theta)g,
\qquad
\rho(\theta)=2-\frac83\theta-\frac{1555}{936}\theta^2.
\tag{16}
\]

On 0≤θ≤1/2 this decreases to 941/3744>1/4. Weighted Bochner therefore
gives

\[
\boxed{\operatorname{gap}(K_\theta)\ge\rho(\theta)>\frac14
\qquad(0\le\theta\le1/2),}                                 \tag{17}
\]

on the whole product space, and hence also within its vacuum sector.
For completeness, if L=Δ−2∇F·∇, the integrated Bochner identity is
∫(Lf)²dν=∫[|Hess f|²+(Ric+2Hess F)(∇f,∇f)]dν.
For a nonconstant eigenfunction −Lf=λf, integration by parts yields
λ²∫f²≥ρλ∫f², hence λ≥ρ. Compactness supplies the spectral expansion and
the Poincare extension. This is the named Bakry–Émery/Bochner method; see
also the reversible-diffusion discussion in Laurent Veysseire,
https://arxiv.org/abs/1105.6080. Our displayed proof fixes the convention
without an implicit factor of 1/2 in the Laplacian.

## 6. The exact transfer obligation to original Yang–Mills

The decomposition is

\[
\boxed{H_\theta=K_\theta+P(\theta-\theta^2/48)I
+\sum_e r_e,\qquad
\|r_e\|_\infty\le\varepsilon(\theta),\quad
\sup_f\sum_{e:f\in D(e)}\|r_e\|_\infty\le121\varepsilon(\theta).}
\tag{18}
\]

This is an exact reduction with the original gauge symmetries, kinetic
metric and spatial interactions retained. The unaccounted interaction
starts at θ³ and its local norm has no volume factor. Its global operator
norm still has the bound Pε. Applying ordinary global min–max to (17)
would lose up to 2Pε from the gap, which does not give a useful uniform
result as P grows. A scalar shift removes the common bulk energy but
does not remove this nonconstant operator.

The next target is thus a **local stability or actual-vacuum return bound
for the explicit residual (18)**. Uniform comparison curvature by itself
is not such a theorem. No identification of νθ with the actual vacuum
measure, and no norm-resolvent or continuum limit, is assumed. The
physical spacing and coupling trajectory of YC6 remain additional data;
the small-θ window here is not the weak-coupling continuum trajectory.

This improves the many-cell obligation beyond an unspecified extensive
error: both the local remainder and a gapped comparison are explicit, with
all constants. It does not yet finish the original YM gap problem.

## 7. Certificate, sources and continuity

```text
python physics/yc8/yc8_local_vacuum_dressing.py --check
python -m unittest discover -s physics/yc8 -p 'test_*.py' -v
```

The packet checks seven exact local identities, twelve finite-geometry
controls on three boxes, and five exact constant/bound conditions. Eleven
independent tests include rational noncommuting SU(2) words with either
shared-edge orientation, direct geodesic differentiation, graph closure,
stencil incidence, and rejection of collapsed geometries. The analytic
all-size incidence and norm arguments above are part of the theorem;
finite checks alone are not claimed as its proof or as formal verification.

YC7 is kept unchanged. For the framework interpretation, see
`FRAMEWORK_TVSP_INTERPRETATION.md`: corner potentials, paired sources,
centre degree filters and response towers guide the construction, while
the geometric circle/square/corner-energy intuition remains a separate
target requiring its own pairing and energy map. The classical identities
used here are retained by name rather than presented as newly invented
general theorems.
