# YC20 — retained return with the full metric and an energy-step certificate

9 October 2026. Continues YC19 at `c51ac1e`. Frozen packets are unchanged.

**Result.** A represented primitive cut gives complementary source channels;
balance is a condition on their readings, not disappearance of the cut.
YC19's exact excitation quotient can be reduced through a retained
cut without discarding its non-product metric. The retained metric is a Schur
complement. Orthogonalizing the cut in that metric changes the kinetic block
and the hidden source as well as the pairing. The resulting self-adjoint
energy pencil has an exact two-energy step, a positive returned-norm slope,
and a complete residual enclosure. Nested eliminations compose exactly.

This supplies a metric-aware spectral adapter for a prospective correlated
block update. It does **not** bound the new block metric uniformly in volume,
localize its square root, improve YC19's interaction budget, or establish an
iterated RG contraction. YC15's actual lattice gap window remains unchanged.
The new finite examples are controls, not computed spectra of a larger YM box.

## 0. Primitive order: cut, complementary readings, balance, then potential

Owner correction while this stage was being developed: what was called the
centre's "generating potential" is intended first as the **cut operation**.
Use the existing native cut algebra, rather than postulating an energy hill
or a scalar potential as the framework's first object. The order is cut ->
complementary channels -> their source, relations and balance -> response and,
where integrability allows it, a reconstructed potential. The earlier native
iota and carrier-dimension derivations are not being redone by this stage.

On a declared carrier, let the represented cut satisfy

\[
 K_{\rm cut}^2=I,\qquad K_{\rm cut}^{\dagger}M=MK_{\rm cut},\qquad M>0.
\]

The source v and positive pairing M must be specified or obtained from the
earlier construction; the bare involution alone does not select them. Then

\[
 P_\pm=\frac{I\pm K_{\rm cut}}2,\qquad P_+P_-=0,\qquad P_++P_-=I,
\]
\[
 w_\pm=\|P_\pm v\|_M^2,\qquad w_++w_-=\|v\|_M^2,\qquad
 F_{\rm cut}=w_+-w_-=v^\dagger MK_{\rm cut}v.\tag{0}
\]

The centre's balanced reading is Fcut=0. It is not Kcut=0: the nonzero
involution continues to make a distinction even when the two weights agree.
Both channels must actually be nontrivial for a nonzero balanced source.

An M-isometric exchange E with E Kcut E^-1=-Kcut makes the balance explicit:
an exchange-fixed source Ev=v has Fcut=0. Dynamics commuting with E preserves
that exchange-fixed class. Merely having opposite cut signs does not impose
equal source weights or prove thermodynamic/dynamical equilibrium. Such an
exchange requires compatible channel carriers and is not asserted for every
finite-retained/infinite-hidden YM cut.

All these readings are unchanged under v'=Z^-1 v, Kcut'=Z^-1 Kcut Z,
M'=Z* M Z. This is the cut-first interpretation of the metric transport below.
UP8's map from a *supplied* potential family to a cut remains a valid model
construction; it does not make that potential a primitive of the foundation.

## 1. Why the projective/winding drafts lead to this target

The owner supplied Transforma's projective eigengeometry, the Monti winding
sheet, and a merged lambda-space draft. Their useful questions are: what
survives removal of common scale, which frame is retained, and what information
is lost by identifying a complete turn with its endpoint? We use those as
questions, not as previously proved operator or physical identifications.

Three repairs are explicit.

* A projective pair is `[lambda_+:lambda_-]`, on its declared chart; one scalar
  modulo arbitrary nonzero scalar multiplication is not a coordinate on P^1.
  For a smooth nonzero scalar f, d(d log f)=0 locally. A nonzero winding around
  an excluded zero is possible, but is not a smooth curvature density there.
* An integer winding nu cannot be recovered from exp(2 pi i nu)=1. Retain the
  lift/integer separately. A smooth nonzero phase map on a fixed closed loop
  cannot change its integer degree under a homotopy through such maps.
* A response metric becoming singular does not by itself prove a mass, a
  curvature singularity, or a new spatial dimension. Pointwise normalization
  is not permission to discard derivatives of a moving normalization frame.

UP8 already gives a native nontrivial transport:
a=ell R dphi, F=R d ell wedge dphi. Its ell and phi come from an actual
cross-corner response and a declared half-share protocol. The projective
eigenvector control in section 8 uses a *different*, explicitly defined
connection; it is not silently identified with UP8 or Yang--Mills curvature.

## 2. Actual carrier and admissibility of a cut

Use YC19 on finite even spatial tori with sides >=4, real internal cube
couplings in [0,1] and real interface |zeta|<=rho=1/4320. It supplies

    A = H0 + Q U Q,     A* G = G A,
    G = M_QQ - M_Q0 M_00^-1 M_0Q > 0,     M = S* S.

A acts on the quotient by the exact ground line. Its spectrum is the physical
excitation spectrum, not the ground energy. At each fixed volume G is bounded
and boundedly invertible. No volume-independent bounds on G or G^-1 are used.

Set K=G A, as a Hermitian energy form (distinct from Kcut) with its transported domain, and choose
a retained subspace P and hidden complement Qh in this quotient. The symbols
Qh and G here are not the pre-ground-removal Q and full M above.

The following block statements are immediate for finite matrices. For the
actual unbounded carrier require: the split respects the closed form domain;
the bounded triangular frame below and its inverse respect that domain; the
retained columns and their metric correction admit the indicated source
pairings; and the hidden form has a proved lower bound. A finite source
matrix alone proves none of these hypotheses. One may work with form-dual
residuals; the displayed residual Gram assumes its columns are Hilbert vectors.

For a finite retained space, use admissible columns with their corrected
hidden components in the energy form domain. The hidden inverse always acts
on the *entire* complement, including all electric harmonics and spectators.

Write the blocks, in retained/hidden order, as

\[
 G=\begin{pmatrix}G_p&N\\N^\dagger&G_h\end{pmatrix},\qquad
 K=\begin{pmatrix}K_p&J\\J^\dagger&K_h\end{pmatrix}.
\]

Require G>0 and Kh>=d Gh in the form sense for some real d; work at real z<d.
These are the inputs of the adapter, not newly established uniform YM bounds.

## 3. The metric and the source both change under diagonal reset

Define the metric correction X, retained metric g, and triangular frame T:

\[
 X=G_h^{-1}N^\dagger,\qquad g=G_p-NG_h^{-1}N^\dagger,\qquad
 T=\begin{pmatrix}I&0\\-X&I\end{pmatrix}.
\]

Direct congruence gives

\[
 T^\dagger GT=\operatorname{diag}(g,G_h),\qquad
 T^\dagger KT=\begin{pmatrix}k&b\\b^\dagger&K_h\end{pmatrix},
\]
\[
 k=K_p-JX-X^\dagger J^\dagger+X^\dagger K_hX,\qquad
 b=J-X^\dagger K_h.\tag{1}
\]

The metric-orthogonal retained representative of u is (u,-Xu). It is the
unique hidden choice minimizing the G-norm for fixed retained u:

\[
 \|(u,v)\|_G^2=u^\dagger gu+\|v+Xu\|_{G_h}^2.\tag{2}
\]

There is an explicit compatible cut behind this formula. Start with the
diagonal cut J0=diag(I,-I) in the reset frame and transport it back:

\[
 K_{\rm cut}=TJ_0T^{-1}=\begin{pmatrix}I&0\\-2X&-I\end{pmatrix},\qquad
 P_+(u,v)=(u,-Xu),\quad P_-(u,v)=(0,v+Xu).
\]

Thus Kcut^2=I, Kcut*G=G Kcut, and (2) is exactly the source split (0).
The metric Schur complement is the retained channel's pairing. This is an
explicit realization of the primitive cut on YC19's represented carrier,
not a claim that Kcut alone fixes G, the source, the energy or a mass scale.

The transformed hidden source b* is not generally the old J*. Carrying g
but leaving the source unchanged would still change the problem.

The source has a direct cut-first characterization. Since Kcut acts as +I
and -I on its two channels,

\[
 \boxed{P_- A P_+=\tfrac12 P_-[A,K_{\rm cut}]P_+.}
\]

In the reset coordinates the lower off-diagonal of T^-1 A T is Gh^-1 b*.
Thus b records the dynamics crossing the cut, with its metric restored.
If A commutes with this cut then b=0 and the hidden return vanishes. This
operator/cut commutator is not by itself spacetime or connection curvature;
A still comes from the specified YM dynamics, not from the cut alone.

The exact retained pencil is

\[
 \boxed{F(z)=k-zg-b(K_h-zG_h)^{-1}b^\dagger.}\tag{3}
\]

It is also the direct Schur complement of K-zG before the triangular reset.
Thus scalar trial energy z couples to the original metric off-diagonal N;
simply eliminating K and then subtracting zGp is not equivalent.

If desired, whiten with diag(g^-1/2,Gh^-1/2) after T. This yields an ordinary
self-adjoint block operator with retained block g^-1/2 k g^-1/2 and source
g^-1/2 b Gh^-1/2. This is an exact fixed-volume change of coordinates, not a
local or uniformly bounded square-root construction.

## 4. A channelwise projective reserve

Let C=Gp^-1/2 N Gh^-1/2. Then

    g=Gp^1/2 (I-CC*) Gp^1/2.                             (4)

Consequently ||C||<=c<1 implies g>=(1-c^2)Gp. For one retained and one
hidden channel the factor is exactly

    chi = det(G)/(Gp Gh) = 1-|N|^2/(Gp Gh).

This is the same two-channel Schur algebra as CA1, now on the declared
excitation metric. It does not identify CA1's thermodynamic response with G.
For multiple channels the *smallest eigenvalue* of I-CC*, rather than merely
a determinant, controls the retained floor. Changes of invertible frames
separately in P and Qh preserve these normalized correlation singular values.

For a nested sequence of such eliminations, in consistently transported
coordinates, a proved c_j bound at each step yields

    g_final >= product_j(1-c_j^2) (G_initial)|_final.      (5)

For sum c_j^2<=q<1, every finite product is >=1-q. More generally a summable
sequence with every c_j<1 has a positive infinite product. This conditional
criterion does not establish a compatible infinite-volume operator limit.
No such scale-dependent c_j estimates for the YM block iteration are proved
here. In particular a positive reserve for each fixed step is insufficient
if their product tends to zero.

## 5. Energy step: the hidden norm returns to the retained reading

Put D_z=Kh-zGh, and in the original coordinates let

\[
 W_z=\begin{pmatrix}I\\-X-D_z^{-1}b^\dagger\end{pmatrix}.
\]

It solves the complete hidden equation `(K-zG)W_z=(F(z),0)^T`.
For x,y<d, multiplication on both sides gives the exact step identity

\[
 \boxed{F(y)-F(x)=-(y-x)W_y^\dagger GW_x.}\tag{6}
\]

On the diagonal this becomes

\[
 \boxed{-F'(z)=W_z^\dagger GW_z
       =g+bD_z^{-1}G_hD_z^{-1}b^\dagger\succeq g.}\tag{7}
\]

The hidden information is now a positive *norm return*, in addition to the
retained metric. This is the precise operator reading of retaining the missing
channel during an observer change. It is an energy-parameter identity, not a
spatial renormalization law or identification of z with observation lambda.

Furthermore

    F''(z)=-2 b D_z^-1 Gh D_z^-1 Gh D_z^-1 b* <=0.        (8)

The last sign follows by whitening Gh and using the spectral theorem for
Gh^-1/2 Kh Gh^-1/2. Thus F is decreasing and operator concave below the hidden
floor, even when the original A is nonnormal in the Euclidean pairing.

In the original full pencil, successive hidden eliminations are associative:
eliminating blocks C then B agrees with eliminating B+C together, whenever
the required inverses exist. The same holds for the metric alone. For the
energy-dependent lifts, W_total=W_first W_second, and

    -F_final'(z)=W_total*G W_total.                      (9)

After the first elimination the intermediate pencil is generally nonlinear
in z; do not replace its derivative metric by the original constant G block.
This distinction is essential if several updates are chained.

## 6. Complete residual bounds and a spectral decision

Given trial hidden columns Y for D_z^-1 b*, retain the entire residual

    R=b*-D_z Y,
    Sigma_trial=bY+Y*b*-Y*D_zY.

Completing the square gives the exact identity

    bD_z^-1b* = Sigma_trial + R*D_z^-1R.

Since D_z>=(d-z)Gh,

    0 <= R*D_z^-1R <= R*Gh^-1R/(d-z).                  (10)

Define E_R=R*Gh^-1R/(d-z). Then

    F_lower=k-zg-Sigma_trial-E_R <= F(z)
              <= k-zg-Sigma_trial=F_upper.              (11)

Every omitted harmonic enters R; no closure of the hidden space in the trial
span is assumed. Outward bounds on the finite residual Gram and retained
matrices can replace exact finite arithmetic. A non-Hermitian trial Y is fine;
Sigma_trial is Hermitian by construction.

With a finite retained cut and a positive hidden block, form congruence gives

    number of physical excitations below z = N_-(F(z)).  (12)

This uses the already removed ground line. Therefore F_lower>=0 proves no
excitation below z, and F_lower>0 excludes z as well. More generally

    N_-(F_upper)<=N_-(F(z))<=N_-(F_lower).

Once the two finite inertias agree, the exact count is known. This is the
generalized-metric Schur/Feshbach and inertia principle, not a newly invented
spectral theorem. The contribution here is the full YC19 metric/source adapter,
the energy-step norm ledger, and reproducible controls for its failure modes.

YC19's fourth-order error is in a component l1 norm relative to H0. It cannot
be inserted into (10) as a Hilbert-metric residual bound without a proved norm
adapter. Its existing l1 resolvent certificate remains valid separately.

## 7. What is concretely checked and what remains for YM

The certificate imports YC19's exact four-dimensional ground-similarity
control and reduces its three-dimensional excitation quotient through a
one-dimensional retained cut. It verifies the transported Hermitian pencil,
its compatible involution and positive cut weights, full Schur metric,
the source correction and an enclosure at z=1/2 using
trial columns solved at z=0. The hidden floor d=1 is checked exactly. A second
example uses complex off-diagonal entries, so conjugate transpose matters.
Direct generalized determinants and full matrix inertia independently check
the reduced spectral decisions. A three-block control checks both elimination
orders, lift composition and the derivative metric.

On the actual YM carrier equations (1)--(12) apply to admissible cuts with
proved full hidden floors and residual bounds. No new interacting block
metric or source columns have been numerically enclosed here. The current
actual uniform gap remains >9/40 for theta_c<=1 and |zeta|<=1/4320 in YC15's
physical coupling window (nonnegative interfaces in that theorem).

The next quantitative gate is to choose a spatial block cut, bound its full
metric correlation C and corrected source b, and show the new interaction
budget improves after paying those costs. This stage establishes the exact
map and a certificate it must pass; it does not claim that improvement.

## 8. Independent geometry/scale controls

**Eigenvalues are not the whole frame.** On C^2 choose the Pauli realization
of three anticommuting cuts and n=(sin theta cos phi,sin theta sin phi,cos theta).
For H=aI+epsilon n.sigma with epsilon>0, the eigenvalues are a+-epsilon and
P=(I+n.sigma)/2 is the positive-branch projector. In the convention
F_theta,phi=-i Tr(P[partial_theta P,partial_phi P]), direct algebra gives
F_theta,phi=sin(theta)/2 while the eigenvalue ratio is angle-independent.
This is a standard Berry-projector control, not a new YM curvature theorem.
The same curvature persists as epsilon tends to zero through positive values,
while the two-level gap 2epsilon tends to zero. At epsilon=0 the eigencut is
no longer singled out by H; no continuation through the degeneracy is claimed.

**A vanishing metric reserve need not close a physical gap.** Let
H=diag(1,4), S=[[1,t],[0,1]], A=S^-1 H S, G=S* S. Their physical gap above
zero is always 1 and their separation is 3. Yet the retained metric reserve
is 1/(1+t^2), which tends to zero. This is coordinate ill-conditioning. At
t=2 the ordinary Hermitian part of A has a negative eigenvalue, although its
physical spectrum remains {1,4}. Dropping the metric is demonstrably wrong.

**Capacity is not curvature.** The constant matrix [[1,r],[r,1]] is positive
for |r|<1 and flat as a spatial metric for each constant r. Its inverse norm
is 1/(1-|r|). Divergent inverse response alone does not force spatial curvature.

## 9. Sources and replay

Native sources: CA1, UP8, YC14 and YC19, with exact hashes in the result.
The uploaded idea sheets are registered by content hash in the code; their
unsupported physical identifications are not inputs to the proof.
Classical source for the spectral reduction: G. Dusson, I. M. Sigal and
B. Stamm, *The Feshbach--Schur map and perturbation theory* (2021),
https://arxiv.org/abs/2105.02058. The identities needed here are proved above.
Projector control: M. V. Berry, *Quantal phase factors accompanying adiabatic
changes* (1984), https://doi.org/10.1098/rspa.1984.0023.

Run with Python 3.12:

    python -B physics/yc20/yc20_metric_retained_return.py --check
    python -B -m unittest discover -s physics/yc20 -p 'test_*.py'

Exact finite checks complement the written full-carrier arguments. They are
not proof-assistant verification and are not a continuum mass-gap certificate.
