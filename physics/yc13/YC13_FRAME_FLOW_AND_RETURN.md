# YC13 — the observed frame, mixed lambda curvature and the complete cube step

9 October 2026. Continues YC12 at `89769b8`, with the native connection
of UP8 and the actual source families of UP9/YC11. Frozen packets remain
unchanged. This is a written analytic development with exact/outward
certificates, not proof-assistant verification.

**Results.** The same native transport extends over the lambda family.
An explicit globally admissible potential has zero native curvature on
every fixed-lambda state surface, but strictly nonzero mixed lambda/state
curvature for every lambda>0. Observation as the transported signed
TVSP→VTSP chart preserves the construction. On the actual cube, YC12's
parity gives an exact all-order positive hidden-return step. It certifies
the following full coupling intervals for the gauge-invariant cube:

| Coupling interval | Actual gauge-sector gap lower bound |
|---|---:|
| 0≤θ≤1 | 11 |
| 0≤θ≤3/2 | 42/5 |
| 0≤θ≤2 | 27/5 |

These are bounds on the ordinary twelve-link cube, not the collapsed
one-site model. The decreasing certified floors do not prove that its
actual gap decreases. YC12's charged-block gap and weak-interface lattice
window are unchanged. No continuum mass, physical time or spatial RG
identification follows from this stage.

## 1. Preserve the owner's movie convention

The owner's latest clarification is: **TVSP before observation; the screen
reading is VTSP; the film we observe is the lambda-zero reading.** Retain
this interpretation instead of exchanging the pre- and post-observation
labels in future discussions. The generalized lambda family describes
additional readings/frames to be investigated. It is not already a
derived physical clock, a measured movie-frame count or a quantum
measurement channel.

The four-slot layout swap is Π(T,V,S,P)=(V,T,S,P). For the signed arrow
convention clarified earlier, the normalized two-component carrier uses
J=[[0,−1],[−1,0]]. These two descriptions retain their different sign
actions. Units, conjugate pairings, cut and source travel with the chart;
the signed map does not negate a physical temperature or volume.

We keep two separate objects: a generating scalar potential U_lambda and
the corner-average reading U_mid=U−(S U_S+V U_V)/2. UP9/YC11 proved an
additional source-specific identity U_mid=−lambda ∂_lambda U/2. That is
not a universal identity for arbitrary trial families. Searching a general
potential family is allowed; its symmetry, positivity, flat limit and
source correspondence remain explicit tests.

## 2. The same transport connects neighboring family members

Take a jointly C⁵ potential U_lambda(S,V), with UP8's derivative
denominators nonzero and cut discriminant Δ=p²+q²−w²>0. On a patch with
a continuous seam angle φ=arg(p+iq), let

\[
\ell=\frac{w^2}{\Delta},\qquad
\mathfrak a=\ell R\,d\phi,\qquad R^2=-I,\qquad
d=dS\,\partial_S+dV\,\partial_V+d\lambda\,\partial_\lambda.
                                                               \tag{1}
\]

This extends the existing half-share connection over the family; it adds
no independent angular-velocity parameter. The matrix derivation from
B=(Rκ)² and A=B^−1dB/2 in UP8 is valid for every parameter derivative.
All orthonormal coefficients are multiples of the same R, hence

\[
\mathfrak F=R\,d\ell\wedge d\phi,\qquad
F_{\lambda i}/R=\ell_\lambda\phi_i-\ell_i\phi_\lambda,
\quad i\in\{S,V\}.                                      \tag{2}
\]

At a fixed state x=(S,V), the finite frame transition is

\[
\boxed{T_{b\leftarrow a}(x)=
\exp\!\left[-R\int_a^b\ell(\lambda,x)
                         \partial_\lambda\phi(\lambda,x)d\lambda\right].}
                                                               \tag{3}
\]

It composes exactly: T_(c←b)T_(b←a)=T_(c←a), and T_(a←b) is its inverse.
For a path changing both state and lambda, replace the integrand by
ell dφ. A closed contractible loop has returned angle
Θ=−∮ell dφ=−∫F/R. Open-path angles depend on endpoint frames. A single
one-dimensional sequence alone has no curvature two-form; the mixed
loop compares a lambda step with a state step.

The lambda-zero convention restricts the deformation coordinate. It does
not declare that all other state histories or actual quantum evolution
stop there. A physical movie clock would need its own evolution law.

## 3. Observation before or after the step gives the same transported result

For a fixed signed observation frame J, transform both response frames:
L'=JLJ^T. Then (p',q',w')=(−p,q,−w), Δ'=Δ, ell'=ell and locally
dφ'=−dφ. Also JRJ^T=−R. Thus B'=JBJ^T and

\[
\mathfrak a'=J\mathfrak aJ^T,\quad
\mathfrak F'=J\mathfrak FJ^T,\quad
\boxed{T'_{b\leftarrow a}J=JT_{b\leftarrow a}.}            \tag{4}
\]

The signed angle reverses; its magnitude and existence are unchanged.
A consistent fixed TVSP→VTSP chart therefore commutes with the frame
transition in the intertwining sense (4). The swap itself supplies no
extra physical curvature or energy.

If the screen frame O(lambda,x) moves too, z'=Oz gives instead

\[
\mathfrak a'=O\mathfrak aO^T-dO\,O^T,\qquad
\mathfrak F'=O\mathfrak FO^T,\qquad
T'_{b\leftarrow a}=O_bT_{b\leftarrow a}O_a^T.              \tag{5}
\]

The derivative term follows by differentiating z'=Oz in dz=−a z.
Dropping it would mistake motion of the observer's axes for a field.
This is the standard connection/frame covariance, specialized to the
native half-share protocol. A primary gauge-theory reference for the
connection transformation and field-strength covariance is David Tong,
*Gauge Theory*, section 2.1:
https://davidtong.org/pdfs/teaching/gauge-theory/gauge.pdf.
The two-real-component native connection is not identified with the full
SU(2) quantum gauge field by this similarity.

## 4. A family with flat state slices and curved frame-to-state transport

Choose fixed dimensionless S,V>0 and lambda≥0:

\[
\boxed{U_\lambda(S,V)=S^3/V+\lambda(S^2+V^2).}             \tag{6}
\]

At zero this is UP8's closed flat base. Its ordinary Hessian is strictly
positive globally: the Hessian of S³/V has upper-left entry 6S/V and
determinant 3S⁴/V⁴, while the added Hessian is 2lambda I. All source
denominators in UP8 are nonzero. Every member is homogeneous of degree
two, so U_mid=0 for the entire family even though ∂_lambda U=S²+V²>0.
This explicitly distinguishes the generating potential from its derived
centre average and from the special amplitude-flow identity.

Let z0=S³/V, r=V/S, x=lambda r and y=lambda r³. Computing the two corner
operations from U, not prescribing a response matrix, gives L=z0 Lbar.
For its traceless coefficients define

\[
p_0(x)=\frac{16x^4+188x^3+756x^2+1350x+945}{8(x+3)^3},
\quad p=p_0(x)-2y,
\]
\[
q=\frac{9(x+2)(2x+3)}{4(x+3)^2},\qquad
w=\frac{3x(2x+3)}{4(x+3)^2},\qquad
\Delta=p^2+q^2-w^2.                                      \tag{7}
\]

The common z0 cancels from ell and φ. The exact positivity certificate is

\[
q^2-w^2=\frac{9(2x+3)^3}{4(x+3)^3}>0.                    \tag{8}
\]

Consequently the native cut is admissible on the whole positive domain,
including lambda=0. Direct differentiation of (7) yields

\[
\boxed{f_{xy}:=F_{xy}/R
=\frac{81x(2x+3)^3}{8(x+3)^6\Delta^2}.}                   \tag{9}
\]

For any fixed lambda both x and y depend only on r=V/S, so the pullback
to the (S,V) state surface has **F_SV=0 identically**. At lambda>0 the
cut itself is still open (w>0, B≠I); zero connection curvature on that
surface is not the same assertion as zero cut defect.

Across the family,

\[
dx\wedge dy=2\lambda r^3d\lambda\wedge dr,
\qquad
\boxed{F_{\lambda r}/R=2\lambda r^3f_{xy}>0
                  \quad(\lambda,r>0).}                  \tag{10}
\]

Thus flatness of every state slice does not establish flatness of the
family's transport. At fixed positive r its leading term is

\[
F_{\lambda r}/R=\frac{3072}{1874161}\lambda^2r^4
                 +O(\lambda^3).                         \tag{11}
\]

At lambda=1/4,r=1/2 the exact value is
5184000000000000/959239754252475521293. A counterclockwise rectangle
[0,1/4]×[1/2,1] in (lambda,r), lifted for example by S=1,V=r, has a
strictly negative return angle. A 16×16 interval subdivision bounds
the entire surface integral, including every point between mesh nodes:

\[
\boxed{-2.061\cdot10^{-6}<\Theta<-1.103\cdot10^{-6}}.      \tag{12}
\]

The result file retains the tighter exact dyadic enclosure. The proof
of nonzero curvature is (8)–(10), not a plot or point sample. There is no
singularity at zero on this family: w and ell vanish there, and all mixed
curvature coefficients have regular limits.

**A failed simpler trial is informative.** Any quadratic potential
U=(aS²+2bSV+cV²)/2, with coefficients depending only on family parameters,
has D1D2U=0 for this corner protocol, hence q=−w. On an admissible patch
Δ=p², ell=(q/p)²=tan²φ. Therefore a=R d(tanφ−φ) locally and the entire
extended curvature vanishes, even when its coefficients change. Merely
varying a quadratic frame cannot produce (10). The nonquadratic S³/V
source in (6) is essential for this example.

## 5. What the actual Haar and cube sources permit

UP9 and YC11 have the special amplitude form

\[
U_\lambda(S,V)=\lambda^{-2}W(\lambda S,\lambda V).
\]

At lambda>0 put k=lambda S,h=lambda V. The corner response is lambda^−2
times the response built from W at (k,h); the normalized cut B, ell and
φ are therefore pullbacks from those two source coordinates. This gives
the exact extended-curvature identities

\[
\boxed{F_{\lambda S}=-\frac V\lambda F_{SV},\qquad
F_{\lambda V}=\frac S\lambda F_{SV}.}                     \tag{13}
\]

The vector lambda ∂lambda−S∂S−V∂V annihilates the normalized response and
contracts the full curvature to zero. This additional parameter is thus
a redundant dilation direction for these particular source families;
it does not independently supply a third physical direction or spatial
blocking law. Formula (6) is a shape deformation and need not satisfy
(13). The distinction is part of choosing the family, not a prohibition
on generalization.

YC11's F_SV begins at lambda^14; its two mixed components begin at
lambda^13. UP9 correspondingly gives orders nine and eight. The apparent
1/lambda in (13) is removable on their certified nondegenerate zero
patches. Reusing YC11's pinned full-integral certificates at S=2,V=1
gives F_lambda,S<0 and F_lambda,V>0 at lambda=1/4 and 1/2, with complete
outward errors. These are new components of the same known source
curvature, not an independent kinetic energy calculation.

## 6. Actual cube: all hidden interactions have an exact positive step

Use H_theta=−sum_(12 edges) Delta_e+theta(6−sum_p Wp), with
Wp=Tr(Up)/2, unit S³ link metrics and normalized product Haar measure.
The physical cube carrier obeys Gauss invariance at all eight vertices.

Keep YC12's gauge cut P: the vacuum and six fundamental face readings,
with complete hidden floor QH0Q≥18. Let r0=QVP, D0=QH0Q, T=QSQ,
S=sum Wp, ||T||≤6. At physical test energy z=6θ+ζ, define

\[
\Sigma_\theta(\zeta)=r_0^*(D_0-\theta T-\zeta)^{-1}r_0,
\quad R_0(\zeta)=\sum_{E=18,24,26,32}\frac{B_E}{E-\zeta}.
                                                               \tag{14}
\]

YC12's face-flipping involution commutes with D0, anticommutes with T,
and makes every source column even. With A0=D0−ζ, C0=A0^−1/2 T A0^−1/2
and v=A0^−1/2r0, the complete compressed Neumann expansion is

\[
\boxed{\Sigma_\theta(\zeta)=
v^*(I-\theta^2 C_0^2)^{-1}v,\qquad
u=6\theta/(18-\zeta)<1.}                                 \tag{15}
\]

Odd powers vanish by parity; every even power is positive. This identity
uses the complete hidden space, including excursions outside the four
source shells. It is not the inverse of a source-span truncation.

Let t=θ² and B0=C0². Its entire derivative tower satisfies

\[
\partial_t^n\Sigma=n!\,v^*B_0^n(I-tB_0)^{-n-1}v\succeq0.
                                                               \tag{16}
\]

By the scalar spectral inequality for 0≤B0≤[6/(18−ζ)]², if θ2≥θ1≥0
and u2<1 then

\[
\boxed{\Sigma_{\theta_1}\preceq\Sigma_{\theta_2}
\preceq\frac{1-u_1^2}{1-u_2^2}\Sigma_{\theta_1}.}          \tag{17}
\]

At ζ=0, θ1=1/2, θ2=1 the upper multiplier is exactly 35/32. For the
actual self-energy θ²Σ the multiplier is 35/8, including the source's
coupling factor. The comparison is at fixed **centred** ζ, not fixed
physical z; no monotonicity of an excitation gap is inferred from it.
This is an exact quantum coupling step. The native lambda flow is kept
separate; neither has been proved equal to a spatial coarse-graining map.

Any unitary change of retained and hidden observation frames transports
the complete expression (14) by conjugation, so its positive order and
the bounds survive the TVSP/VTSP interpretation where such an intertwiner
is supplied. Reversing an entire source changes neither its return nor
the spectral content. This reuses the covariance already established in
the signed-compass handoff and YC12.

## 7. A full-window gauge-gap certificate from that return

The uniform face vector maximizes each resolved source Gram. Hence

\[
\|R_0(\zeta)\|=
\frac{1/2}{18-\zeta}+\frac{1/2}{24-\zeta}
 +\frac{3/2}{26-\zeta}+\frac{1/4}{32-\zeta}=b(\zeta).
                                                               \tag{18}
\]

The face block of the seven-dimensional Schur operator is bounded below
by h(θ,ζ)I, where

\[
h(\theta,\zeta)=12-\zeta-
 \frac{\theta^2 b(\zeta)}{1-[6\theta/(18-\zeta)]^2}.       \tag{19}
\]

For fixed 0<ζ<12 this decreases with θ while the hidden reserve is
positive. If h(θmax,ζ)>0, the face block stays positive throughout the
interval. The retained vacuum diagonal is −ζ, its hidden-source column
is zero, and eliminating the positive face block leaves a negative
scalar. Thus there is exactly one eigenvalue below 6θ+ζ in the full
gauge carrier. The constant trial gives E0≤6θ, so E1−E0≥ζ uniformly
over the entire coupling interval, not merely its endpoints.

| θmax | ζ | u at endpoint | Exact positive face reserve h |
|---:|---:|---:|---:|
| 1 | 11 | 6/7 | 1663/10140 |
| 3/2 | 42/5 | 15/16 | 335196/1307735 |
| 2 | 27/5 | 20/21 | 2473726/12436735 |

The replay checks these fractions and the complete rank-seven inertias.
The analytic monotonicity of (19) certifies all intermediate couplings.
This closes the whole-window gauge-gap statement left separate from
YC12's point enclosures. It does not replace the all-boundary-charge
gap needed when joining cubes; using 27/5 as that block gap would be
invalid. The next lattice task remains improved boundary source/inverse
norms for stronger interfaces, with the observer and source transported
together.

## Replay

Run `python physics/yc13/yc13_frame_flow.py --check` and
`python -m unittest discover -s physics/yc13 -p 'test_*.py'` from the
repository root. Certificates store exact fractions and outward 160-bit
dyadic interval arithmetic, with SHA256 pins for all reused stages.
The two classical tools are connection covariance/holonomy and positive
operator functional calculus combined with the full Schur inertia law;
the explicit potential, mixed curvature and cube constants are derived
here. Curvature-energy normalization and a continuum mass remain open.
