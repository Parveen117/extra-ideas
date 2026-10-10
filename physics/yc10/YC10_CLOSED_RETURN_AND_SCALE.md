# YC10 — closed returns enlarge the lattice window; concentration keeps its scale

9 October 2026. Continues YC9 at `9f5fb83`. Same SU(2) spatial lattice,
unit S³ link kinetic metric, complete link Hilbert spaces, Gauss laws and
global centre-even vacuum sector. Earlier packets remain unchanged.

**Actual YM result.** For Hθ=−Σ_e Δ_e+θΣ_p(1−Wp), uniformly in every
allowed finite rectangular periodic spatial volume,

| Side lengths | Coupling window | Vacuum-sector gap lower bound |
|---|---|---|
| Every N_i≥3 | 0≤θ≤1/64 | 12−(1248/5)θ ≥ 81/10 |
| Every N_i≥4 | 0≤θ≤3/160 | 12−(640/3)θ ≥ 8 |

These are respectively ten and twelve times YC9's coupling window.
The written proof improves its full-cluster return, using orthogonality
and the minimum support allowed by Gauss invariance. It is not a finite
harmonic approximation, a bound on a substituted comparison operator,
or a weak-coupling continuum theorem. The classical method remains
Kirkwood–Thomas / weak-interaction dressing, as credited in YC9.

## 1. The shrinking-circle hint: three different limiting objects

The owner proposes that increasing the native λ-tower compresses the
circle toward a point, concentrating energy at its corners, possibly
to infinity. Retain the concentration question, with its quantity and
scale specified. Three exact controls prevent assigning the wrong limit.

### 1.1 The native potential and its response can have different limits

UP1 has a homogeneous potential Φ(sx)=s^kΦ(x) and partial-cut corners
Φ_A=Φ−Σ_(i∈A) x_i∂_iΦ. All these corner readings scale as s^k.
For 1<k<2, they tend to zero as s→0, while nonzero Hessian components
scale as s^(k−2) and can diverge. For the concrete positive radial
Φ(x)=|x|^k, the radial/tangential response eigenvalues are
k(k−1)|x|^(k−2) and k|x|^(k−2).

Thus an infinite response density is compatible with a vanishing corner
potential. UP1's λ rescaling is λ→λs^(k−2); this is a typed response
parameter, not an identified bare lattice coupling.
LT1's separate gravity example approaches its nonzero horizon radius r_s,
and its normalized response becomes rank one. LT1 explicitly does not
define an energy density. Its inward tower is not already a YM continuum
map. UP3's derivative-defined diagram sequence is another distinct tower.

### 1.2 Classical Yang–Mills has an exact point-concentration test

For a smooth Lie-algebra-valued connection on Euclidean R^d with finite
curvature square integral, define, at fixed coupling,

\[
A_\varepsilon(x)=\varepsilon^{-1}
A((x-x_0)/\varepsilon),\qquad
F_\varepsilon(x)=\varepsilon^{-2}
F((x-x_0)/\varepsilon).                                   \tag{1}
\]

Both the derivative and bracket parts of F=dA+A∧A have this scaling.
With q=|F|²/(2g²), the d-dimensional curvature integral satisfies

\[
\int q_\varepsilon(x)d^d x
=\varepsilon^{d-4}\int q(y)d^d y.                         \tag{2}
\]

In four Euclidean dimensions the **action** is unchanged. If q(0)>0,
the value at the concentration point grows as ε^−4. More generally its
finite measures converge weakly to S δ_(x0), where S=∫q: for a bounded
continuous test f, substitute x=x0+εy and apply dominated convergence.
This gives a precise version of concentration to a point, without
requiring infinite total action or an instanton assumption.

For a purely spatial magnetic configuration in three dimensions, the
same scaling gives magnetic **energy** ε^−1 E. A shrinking classical
configuration can therefore have divergent energy. The Euclidean action,
spatial Hamiltonian energy and quantum vacuum excitation gap are different
objects. Neither (1) nor (2) determines E1−E0. No equation here selects
the corner of a TVSP chart as a spatial concentration point.

This scaling is classical; see the curvature functional and scaling
discussion in Fang-Hua Lin, *Lecture 3 — Yang–Mills theory*,
https://www.math.columbia.edu/~flin/Pisa03.pdf. Equations (1)–(2) above
are derived directly and fix the dimension and normalization being used.

### 1.3 A large common energy does not determine the excitation gap

For any scalar c(ε), replacing Hε by Hε+c(ε)I changes neither E1−E0
nor eigenvectors, even if c(ε) diverges. Conversely nonconstant
concentration cannot be removed just by calling it vacuum energy.
Its source and inverse must be bounded. This is precisely the distinction
between YC8's scalar term and its retained local residual.

The spectral choice in the rest of this stage is therefore to improve
the **closed source return**, while keeping the concentration lead as a
well-typed scale question. It is not used as an unproved gap premise.

## 2. Gauss invariance removes the smallest source supports

Use YC9's exact excitation decomposition P_I, creation operators ĉ_I,
and norm

\[
\|c\|_* = \max_e\sum_{I\ni e}|I|\,\|c_I\|.
\]

Run its fixed point only on the closed linear subspace in which each c_I
is gauge invariant and globally centre even. This subspace is preserved:
the seed is zero, every plaquette term and free inverse commutes with
the symmetries, and every P_I does too. A creator built from invariant
c_I and Ω_I commutes with these symmetries. No existence of the actual
ground is assumed in making this restriction.

In a spin-network component of any nonzero invariant c_I, every link of
I carries a nontrivial SU(2) representation. An endpoint with exactly one
such link cannot form a singlet. Hence the support graph has no leaf.
There are no self-loops or parallel links on the ordinary N_i≥3 tori.
It follows that

\[
\boxed{c_I=0\text{ for }|I|<s,\quad
s=3\text{ if all }N_i\ge3,\quad
s=4\text{ if all }N_i\ge4.}                               \tag{3}
\]

The second assertion uses absence of triangles when all sides are ≥4.
The first cannot generally be replaced by s=4: a length-three **adjoint**
winding on a side of length 3 is gauge invariant and centre even. Its
electric energy is 24. It is absent from the low free floor but not from
the full nonlinear return. We retain it. The vacuum free floor remains
12 in both geometries, by YC9's complete low-support argument.

For any four-link interaction support X, (3) improves YC9's creator sum:

\[
\boxed{\sum_{I:I\cap X\ne\varnothing}\|c_I\|
\le\sum_{e\in X}\sum_{I\ni e}\|c_I\|
\le\frac4s\|c\|_*.}                                     \tag{4}
\]

The s here is a minimum source support, not spatial dimension or tower
depth. This improvement holds for differences of two allowed collections
as well, which is essential for a contraction estimate.

## 3. Orthogonal readings cost four, not sixteen

YC9 proves that a fixed operator word from a local plaquette commutator
has at most 16 nonzero excitation components P_J: outside its original
four links every nonzero word has a fixed excitation pattern. The P_J
are mutually orthogonal. Therefore

\[
\boxed{\sum_J\|P_Jv\|\le
\sqrt{16}\left(\sum_J\|P_Jv\|^2\right)^{1/2}
\le4\|v\|.}                                               \tag{5}
\]

The bound holds when restricting to output components containing any
specified root link, and for words containing an extra excitation creator
û_J. The unbounded dimension of each excited link space changes neither
orthogonality nor the number of binary excitation patterns. The earlier
factor 16 was valid but unnecessarily large.

## 4. The refined all-order return and its explicit constants

Let β=max_e Σ_(p containing e)|t_p|, Φ=−Σ_p t_p Wp, and

\[
\mathcal T(c)_I=-H_I^{-1}P_I e^{-C}\Phi e^C\Omega,
\qquad H_I\ge3|I|,
\quad a_s=8/s.
\]

Repeat YC9's root count with (4) and (5). For ||c||_*≤r one obtains

\[
\boxed{\|\mathcal T(c)\|_*\le
\frac{4\beta}{3}(1+2r)e^{a_s r},}                         \tag{6}
\]
\[
\boxed{\operatorname{Lip}(\mathcal T)\le
\frac{4\beta}{3}e^{a_s r}
\big[2+a_s(1+2r)\big].}                                  \tag{7}
\]

Only the unmarked creator sum is reduced by s. The marked creator still
has Σ_(X meeting I)||Φ_X||≤β|I| and is bounded by the weighted norm r.
Thus (1+2r) in (6) is unchanged. The derivative of this factor supplies
the 2 in (7). Dividing this marked cost by s would be an invalid gain.

For the exact similarity-transformed excitation operator U, the same
argument gives

\[
\sum_I\|P_I Uu_J\|\le8\beta e^{a_s r}|J|\,\|u_J\|,
\qquad b:=\|QUQH_0^{-1}\|_{1\to1}
\le\frac{8\beta}{3}e^{a_s r}.                             \tag{8}
\]

The component-sum norm, full domain argument, block-triangular similarity
and Neumann resolvent proof are exactly those of YC9. Whenever (6)
maps the ball into itself, (7) is <1 and b<1, they give

\[
\boxed{\Delta_{\rm vac}\ge12(1-b),\qquad
\Delta_{\rm full}\ge3(1-b).}                              \tag{9}
\]

The full-product conclusion is valid even though the fixed point was
constructed in the invariant subspace: (8) permits an arbitrary u_J.
It proves that the constructed invariant eigenvector is the true unique
ground. This avoids assuming ground-state character to prove existence.

Choose r=1/4 in both cases and use the following outward rational bounds:

| Geometry | s | Bound for exp(a_s r) | β maximum | Map bound | Contraction bound | b bound | Vacuum gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| N_i≥3 | 3 | exp(2/3)<39/20 | 1/16 | 39/160<1/4 | 39/40 | 13/40 | 81/10 |
| N_i≥4 | 4 | exp(1/2)<5/3 | 3/40 | 1/4 | 5/6 | 1/3 | 8 |

For the first exponential, the series through n=5 plus the geometric
tail starting at n=6 gives

\[
e^{2/3}\le
\sum_{n=0}^5\frac{(2/3)^n}{n!}
+\frac{(2/3)^6}{6!}\frac1{1-(2/3)/7}
=\frac{404671}{207765}<\frac{39}{20}.
\]

The second is YC9's bound 33/20<5/3. These rational inequalities give
the endpoint estimates without floating-point root selection.
For arbitrary β within the corresponding interval, (8) gives

\[
\Delta_{\rm vac}\ge
\begin{cases}
12(1-26\beta/5),&N_i\ge3,\ \beta\le1/16,\\
12(1-40\beta/9),&N_i\ge4,\ \beta\le3/40.
\end{cases}                                                \tag{10}
\]

Uniform β=4θ yields the headline table. Nonuniform signed couplings and
the plane-join formula β=2max(|a|+|b|,|a|+|c|,|b|+|c|) remain valid.

The entire tower is controlled: iteration from zero has
||c−c^(m)||_*≤q^m β/[6(1−q)] using YC9's exact first-return norm.
The slower q=39/40 at the larger N_i≥3 endpoint is explicitly retained;
no finite number of iterates is declared equal to the vacuum.

## 5. What a genuine scale tower would have to preserve

Equations (6)–(9) are instances of a reusable product-reference criterion.
For single-factor gap g, interactions on at most k factors, orthogonal
projection count at most 2^k, and a symmetry-preserved creator space with
minimum support s, assume the excitation projections and constructed
creators preserve the sector in question. Put M=2^(k/2), a=2k/s.
Sufficient inequalities are

\[
\frac{M\beta}{g}(1+2r)e^{ar}\le r,\quad
\frac{M\beta}{g}e^{ar}[2+a(1+2r)]<1,\quad
b=\frac{2M\beta}{g}e^{ar}<1.                              \tag{11}
\]

Then a symmetry sector with free floor δ has gap ≥δ(1−b), by the
same proof. This allows correlated reference **blocks** only after their
unique reference ground, gap g, interaction norm β, support k, symmetry
selection s and invariant domains have actually been established. Gauss
constraints at block boundaries can change s; do not inherit s=4 without
proof. No such correlated block construction is claimed in this packet.

At a scale level n the quantities in (11) would become g_n, β_n, k_n,
s_n, r_n and δ_n. Merely shrinking coordinates multiplies both energies
by the same scale and leaves β_n/g_n unchanged. A useful renormalization
step must change or control this ratio by returning eliminated channels.
This is the scale-flow target suggested by the owner's tower language.
No equality of these parameters with UP1/LT1's λ is assumed.

For the declared YC6 physical convention
H_phys=(γg0²/a)H_(κ/g0⁴)+scalar, the present windows require θ=κ/g0⁴
small. A weak-coupling trajectory g0(a)→0 leaves them. The physical
gap lower bound would be (γg0²/a)δ_n(1−b_n), in the representation
actually obtained after each scale step. A finite continuum mass requires
both control of that expression and a continuum construction. Infinite
corner density supplies neither by itself.

Thus the concentration hint has a precise classical scale test, and the
new spectral result comes from retaining the closed, orthogonal return
structure. Next: a correlated reference block with its **actual** local
gap and interface-return norm, suitable for (11).

## 6. Verification and provenance

```text
python physics/yc10/yc10_closed_return_and_scale.py --check
python -m unittest discover -s physics/yc10 -p 'test_*.py' -v
```

The written proof establishes arbitrary finite volume and full local
Hilbert spaces. Exact checks cover both rational windows, support
selection on representative graph carriers, orthogonal projection
controls, native homogeneous scaling and the classical curvature
scaling identity. Tests include the important nonzero adjoint winding
on length-three sides, so the stronger s=4 cannot be used there.
Finite checks are not a machine formal proof or a construction of a
quantum continuum. YC9's algebra, source and domain arguments are
retained and pinned, not overwritten.
