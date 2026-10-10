# UP9 — an actual SU(2) source selects the compass deformation and its flat limit

9 October 2026. Continues UP8 at `78157ee`, with YC1's unchanged face
source. The owner's direction is a general theory containing a flat
lambda-zero compass, with the nonzero-lambda cuts and seam retained.

**New result.** The deformation is now selected by a specified record:
the centred trace and trace-square of an SU(2) Haar face. Its complete
log-partition function, rather than a fitted polynomial, defines U_lambda.
The source has an exact heat identity and an independent-aggregation
step. Applying UP8 gives a flat zero limit and a nonzero native curvature
whose first visible order is **nine**, with full-integral rational
certificates at lambda=1/4 and 1/2. The order is source-dependent: UP8's
constructive potential had order three.

This closes the *selected statistical-source* gate left by UP8. It does
not select the unique deformation of all physical theories, identify a
source-dilation parameter with the YM coupling, or replace a spatial
block return by independent copies.

## 1. Source contract and the potential

Use YC1's c(U)=Tr(U)/2 on SU(2), normalized Haar measure, and

\[
 X=c,\qquad Y=c^2-\tfrac14,\qquad
 Z_c(k,h)=\langle e^{kX+hY}\rangle,\qquad W=\log Z_c.
\]

Both observables have mean zero. Their Haar covariance is diag(1/4,1/16).
Fix dimensionless source axes S,V and their Euclidean coefficient pairing.
Define

\[
 \boxed{U_\lambda(S,V)=\lambda^{-2}W(\lambda S,\lambda V)},
 \quad\lambda\ne0,\qquad
 \boxed{U_0=S^2/8+V^2/32}.                              \tag{1}
\]

The choice of observables, subtraction of their Haar means, source axes
and normalization in (1) are part of the experiment. Once fixed, all
coefficients are determined by the record. There is no claim that a
linear subtraction leaves every cross-corner construction unchanged.
At lambda=1 this is precisely YC1's Psi(S,V)-V/4.

| Required item | Declaration |
|---|---|
| statistical carrier | one SU(2) face, equivalently c in [-1,1] with density (2/pi)sqrt(1-c^2) |
| observables / source | X=c, Y=c^2-1/4; (k,h)=lambda(S,V) |
| state pairing | fixed dimensionless (S,V) axes; same real two-component cut carrier as UP8 |
| response metric | ordinary Hessian of (1), exactly the tilted covariance of (X,Y) |
| cross-corner cuts | D1=S partial_S, D2=V partial_V+m D1, m=-V U_SV/(S U_SS) |
| return metric / transport | B=(R kappa)^2; A=1/2 B^-1 dB, unchanged from UP8 |
| chart convention | T=U_S, P=-U_V, so dU=T dS-P dV; P is a signed source conjugate, not a fluid pressure |
| target | source-selected cut, native seam transport and its scale law |
| symmetry map | U -> -U sends X -> -X, Y -> Y, hence W(-k,h)=W(k,h) |

The ordinary Hessian is strictly positive for every real finite source:
its quadratic form is the variance of aX+bY, which is nonconstant on
the full Haar interval for every (a,b)!=0. Differentiation in (1)
cancels the two lambda factors against lambda^-2. Thus this Hessian
is the covariance without an extra scaling factor. In particular U_SS>0
and m is defined for S>0.

Since X,Y are bounded, Z_c is entire, positive on real sources, and
W is real analytic there. Its value and first derivative vanish at the
origin. Equation (1) therefore has a removable, analytic limit at
lambda=0, locally with every state derivative. The displayed lambda^-2
is not a physical singularity.

## 2. Exact source identities fix the flow and the diagram centre

**S1 — source heat identity.** The observable relation Y=X^2-1/4 gives

\[
 \partial_h Z_c=(\partial_k^2-\tfrac14)Z_c,
 \qquad W_h=W_{kk}+W_k^2-\tfrac14.
\]

The analytic initial source is fixed, not chosen later:

\[
 Z_c(k,0)=\frac{2I_1(k)}k,\quad Z_c(0,0)=1.
\]

This Bessel representation is classical (NIST DLMF 10.32.2 with nu=1);
the proof also follows immediately by integrating the Haar density.
The heat identity is in the **even source h**, not in physical time.
In the normalized potential it reads

\[
 \boxed{\lambda U_V=U_{SS}-\tfrac14+\lambda^2 U_S^2}.       \tag{2}
\]

Written this way it has a regular limit, U_0,SS=1/4. Dividing (2) by
lambda and then discarding the cancelling terms would create a false
singularity. The integral and analyticity select the full solution;
no claim about uniqueness for an unspecified PDE class is needed.

The Haar integration-by-parts identity additionally gives

\[
 k(Z_c-Z_{c,kk})+2h(Z_{c,k}-Z_{c,kkk})-3Z_{c,k}=0.       \tag{3}
\]

Indeed integrate the derivative of
(1-c^2)^(3/2) exp(kc+h(c^2-1/4)); the endpoint term vanishes.
Its zero-source recurrence is (n+3)<c^(n+1)>=n<c^(n-1)>.
This is a source moment/Ward identity, not a spectral gap equation.

**S2 — exact dilation and centre response.** Direct chain rule yields

\[
 \boxed{\lambda\partial_\lambda U=(S\partial_S+V\partial_V-2)U},
 \qquad
 \boxed{U_{\rm mid}=U-\tfrac12(SU_S+VU_V)
                   =-\tfrac\lambda2\partial_\lambda U}. \tag{4}
\]

Thus the scalar at the TVSP centre is exactly the source-dilation
response in this family. This composes UP1's degree filter with the
specified record, rather than equating a drawing with an energy gap.
U_mid can have either sign; it is not a positive matter-energy density.

**S3 — independent aggregation step.** For an integer n let

\[
 (X_n,Y_n)=n^{-1/2}\sum_{j=1}^n(X_j,Y_j),
\]

with independent copies of the Haar face. The log moment-generating
function of this record is exactly U_(1/sqrt(n)). More generally

\[
 \boxed{U_{\lambda/\sqrt b}(S,V)
        =b\,U_\lambda(S/\sqrt b,V/\sqrt b)},\quad b>0.    \tag{5}
\]

The probabilistic regrouping interpretation requires integer b and
appropriate independent blocks. The identity itself follows from (1)
for every positive b. In particular b=4 sends lambda to lambda/2.
Higher cumulants decrease with their usual central-limit powers;
lambda=0 is the Gaussian record with the stated covariance.

For fixed nonzero lambda, U_lambda is a normalized log-partition
potential. Do not assert that it is the log moment-generating function
of an independently divisible probability law for every real lambda:
the finite-copy interpretation above is proved at lambda=1/sqrt(n).

The cross-corner operations commute with the common source dilation in
(5). Therefore, at x=(S,V) and x'=x/sqrt(b),

\[
 L_{\lambda/\sqrt b}(x)=b L_\lambda(x'),\quad
 B_{\lambda/\sqrt b}(x)=B_\lambda(x'),\quad
 \boxed{F_{SV,\lambda/\sqrt b}(x)=b^{-1}F_{SV,\lambda}(x')}. \tag{6}
\]

Here F is the native curvature of the return metric, and the b^-1 is
the area Jacobian. Common energy scaling cancels from the normalized
cut. This is a genuine source step with its metric and two-form
transport, **not a spatial YM renormalization step**: neighbouring
plaquettes share links and are not independent Haar copies.

## 3. Exact first-visible orders of the selected compass

Haar moments are <c^(2n)>=Cat_n/4^n, with odd moments zero. The resulting
cumulant expansion, convergent near lambda=0, is

\[
\begin{aligned}
 U_\lambda={}&\frac{4S^2+V^2}{32}
 +\lambda\frac{V(12S^2+V^2)}{384}
 -\lambda^2\frac{S^2(2S^2-3V^2)}{768}\\
 &-\lambda^3\frac{V(60S^4+V^4)}{30720}
 +\lambda^4\frac{32S^6-192S^4V^2-24S^2V^4-V^6}{294912}
 +O(\lambda^5).                                         \tag{7}
\end{aligned}
\]

Every coefficient is generated from the exact moment/cumulant recursion
in the executable certificate. Source identities (2) and (4) independently
check the available coefficients.

Use exactly UP8's L and decomposition:

\[
 L=\begin{pmatrix}D_1^2U&D_2D_1U\\D_1D_2U&D_2^2U\end{pmatrix}
  =aI+pK+qRK+wR,\quad 2w=D_1(m)SU_S,
\]
\[
 \Delta=p^2+q^2-w^2,\quad \ell=w^2/\Delta,
 \quad\phi=\arg(p+iq),\quad \mathfrak F=R\,d\ell\wedge d\phi. \tag{8}
\]

Restrict to a positive patch S,V>0, **V!=2S**, and sufficiently small
lambda. At zero the ordinary Hessian is positive everywhere, but

\[
 L_0=\operatorname{diag}(S^2/2,V^2/8),\quad
 p_0=(4S^2-V^2)/16,\quad q_0=w_0=0.
\]

Thus Delta_0>0 on the declared patch, B_0=I and native curvature is
zero. At V=2S the traceless part vanishes and the normalized cut is
undefined; this is a cut degeneracy, not a singularity of the Haar
potential. All four ordinary Legendre corner charts remain regular
because the Hessian is positive. At zero their scale frames have
constant coefficients and commute.

**S4 — cancellation and the first cut defect.** Put d=4S^2-V^2. Exact
division of the source derivatives gives

\[
 m=-\lambda V/4+\lambda^3V^3/128
       -\lambda^4V^2d/1536+O(\lambda^5),
 \qquad
 \boxed{w=-\lambda^4S^4V^2/1536+O(\lambda^5)}.            \tag{9}
\]

The lower-order terms of m depend only on V, and D1 removes them.
The sixth cumulant in (7) is required for the lambda^4 term; dropping
it does not give the correct coefficient. Other frame commutators
appear earlier: [V partial_V,D2]=(V partial_V m)D1 already starts at
-lambda V D1/4. The order four assertion concerns the **selected
cross response and its compass**, not every pair of corner operations.

**S5 — seam, lost part and curvature.** Choose a smooth local argument
with phi_0=0 when d>0 or phi_0=pi when d<0. Then

\[
 \phi=\phi_0-\lambda\frac{S^2V}{d}+O(\lambda^2),\qquad
 \ell=\lambda^8\frac{S^8V^4}{9216d^2}+O(\lambda^9),
\]
\[
 \boxed{\mathfrak F_{SV}=
 -\lambda^9\frac{S^9V^4(2S^2+V^2)}{1152d^4}R
 +O(\lambda^{10}).}                                    \tag{10}
\]

Proof: apply (8) to (9), use q_1=-S^2 V/16 and p_0=d/16, and
differentiate the two leading coefficients of ell and phi. Their
Jacobian is precisely the rational coefficient in (10). Analyticity
and Delta_0>0 control the remainders, uniformly on each compact subset
of the declared patch. The leading coefficient is strictly negative.
Consequently every such compact subset has an epsilon>0 on which the
signed native curvature is nonzero for 0<lambda<epsilon. This is a
written local theorem, not a global sign claim for all lambda>0; no
numerical value of that uniform epsilon is asserted.

At S=V=1 the coefficient is exactly -1/31104. Under (5) it scales as
b^-9/2, in agreement with the two-form law (6). The centre reading
appears much earlier:

\[
 U_{\rm mid}=-\lambda V(12S^2+V^2)/768+O(\lambda^2).
\]

So a small or vanishing low-order curvature reading does not establish
that the record has no higher information. The centre, cut defect,
lost part and transport observe different orders of the same source.

## 4. Two finite full-integral curvature certificates

The power series above is **not** used as an uncontrolled approximation
at a finite lambda. The script encloses every required derivative of
the full integral and propagates it with outward rational arithmetic.

At S=V=1 put H=c+(c^2-1/4), so |H|<=7/4. For i+j<=4 compute

\[
 M_{ij}=\langle c^i(c^2-1/4)^j e^{\lambda H}\rangle.
\]

Expand only the exponential through degree N=24, integrate every term
by the exact Haar moments, and include the complete remainder

\[
 \left|M_{ij}-\sum_{n=0}^N\frac{\lambda^n}{n!}
       \langle c^i(c^2-1/4)^j H^n\rangle\right|
 \le (3/4)^j\frac{t^{N+1}}{(N+1)!(1-t)},\quad
 t=7|\lambda|/4<1.                                    \tag{11}
\]

This follows from the exponential Taylor remainder and e^t<=1/(1-t).
Normalize by the strictly positive interval for M00, form the fourth
state jet of log Z, and use it to compute L and its first derivatives.
Four derivatives of U suffice because m uses two, L uses at most three,
and F uses one further derivative. Every arithmetic operation rounds
outward onto a rational grid of spacing 2^-160. No floating-point
number is used for a sign or a bound.

| lambda, S=V=1 | certified enclosure of the coefficient F_SV/R |
|---|---|
| 1/4 | [-1.229129, -1.229128] times 10^-10 |
| 1/2 | [-6.081351, -6.081350] times 10^-8 |

The complete much narrower rational endpoints, positive cut
discriminants, positive ordinary Hessian determinants, signed w and
positive lost weights are retained in `UP9_RESULT.json`. These are two
point certificates with full tails, not an interval-wide or large-lambda
claim. A separate test feeds UP8's exact rational source through the
same jet transport and recovers its known full curvature witness.

## 5. Handoff and limits

This stage now provides a concrete answer to UP8's open source question:
Haar record -> exact generating integral -> cumulants and source flow
-> cross-corner response -> cut-return metric -> native curvature.
The statistical source is an actual SU(2) face. The resulting two-real-
component native connection remains so(2); it is not identified with
the full nonabelian Yang–Mills connection.

The specific normalization explains which lambda is being used:
**source amplitude**, or inverse square-root independent count on its
finite-copy subsequence. It is neither polynomial response squaring,
covariant derivative order, physical time, nor bare YM coupling theta.
Equation (4) is an exact source scaling identity, not a dynamically
derived beta function. Equation (2) changes a source, not a lattice
spacing. Analytic small-source powers do not decide the presence of
nonperturbative mass or logarithmic running in the interacting field.

For an interacting spatial block, independence in (5) is the next
obligation to replace. Its joint log-partition function has connected
mixed cumulants across shared links. Those terms must be retained in
the centre and in the source/inverse return. YC2 already shows nonzero
shared-link covariance on the compact carrier; YC7–YC10 supply actual
spatial source-return machinery. The next useful adapter is therefore
a **correlated block's complete source response and interface error**,
with the physical metric and scale fixed. The current spectral baseline
remains YC10; no new spectral gap bound is claimed here.

## Reproduce, credit and evidence

```sh
python up9_haar_source_flow.py --check
python -m unittest test_up9
```

Thirty-one exact/outward checks and twelve focused tests cover the
source coefficients, identities, finite full-integral signs, negative
controls and transported predecessor witness. Written analytic proofs
and rational computational certificates have separate roles; this is
not proof-assistant verification.

Classical ingredients: Haar/Bessel integral, moment-cumulant recursion,
exponential-family covariance, independent-sum/central-limit scaling,
and the logarithm of a heat identity. These are not claimed as newly
discovered mathematics. The new stage is their explicit source-to-native-
compass adapter, cancellations, first-visible orders, centre identity
and finite full-integral certificates in the declared framework.

- YC1 and CT1: unchanged SU(2) source and moments.
- UP1 / SS1: degree filter and typed scale response.
- UP6 / UP8 / LC1: unchanged cross-corner cut and native transport.
- NIST DLMF 10.32.2: <https://dlmf.nist.gov/10.32.E2>, checked 9 October 2026.

The JSON pins the exact source, note, implementation and tests by SHA-256.
