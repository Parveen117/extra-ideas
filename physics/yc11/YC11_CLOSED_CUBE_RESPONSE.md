# YC11 — the closed cube's complete boundary return and deformed compass

9 October 2026. Continues UP9 at `b6e988a`. Reuses CZ1's exact character
gluing, YC8's ordinary-link distinction and UP8/UP9's source-to-transport
construction. Earlier packets are unchanged.

**New result.** On an actual twelve-link SU(2) cube, the six-face joint
record supplies a complete boundary return, an explicit bound on all its
omitted representation channels, and a source-selected native compass.
Every proper subset of face traces is independent under the base Haar
measure, yet the full six-face cumulant is 1/1024. It is the closure
that the lower-order observer misses. The selected cut first opens at
order lambda^4; native curvature first appears at **lambda^14**, after
an exact order-twelve cancellation. Complete-integral rational bounds
certify the curvature and centre-memory signs at two finite sources.

The carrier here is a finite Wilson **configuration** integral, with all
twelve links and all representation channels retained. A correlated
quantum block's gap and kinetic inverse are not constructed in this
stage. YC10 remains the spectral baseline.

## 1. The actual join: twelve links, six faces, one shared boundary

Take the oriented boundary of one ordinary unit cube in a spatial cubic
lattice: eight vertices, twelve independent SU(2) link variables, and
six elementary square plaquettes. Each link has one Haar integration;
its two incident faces do not get duplicate copies. Gauge transformations
act at the eight vertices. Write W_p=Tr(U_p)/2 in [-1,1].

The three faces meeting the lower corner form a disk L, and the three
meeting the opposite corner form a disk R. Their common boundary is a
six-link loop. Put

\[
 A=\sum_{p\in L}W_p,\qquad B=\sum_{p\in R}W_p,\qquad
 Z(k,h)=\left\langle e^{kA+hB}\right\rangle_{SU(2)^{12}}.
                                                               \tag{1}
\]

This is anisotropic Wilson configuration weighting on the cube. The
boundary is topologically a two-sphere embedded in the three-dimensional
lattice. Its exact evaluation does not integrate a general collection
of interacting cubes, nor a four-dimensional quantum field.

| Adapter item | Declaration |
|---|---|
| source carrier | actual cube with twelve Haar links and local vertex gauge action |
| selected observables | A,B, the sums of the opposite three-face disks |
| source axes | k=lambda S, h=lambda V, fixed dimensionless S,V and Euclidean coefficient pairing |
| ordinary response metric | Hessian of U_lambda=lambda^-2 log Z(lambda S,lambda V), the tilted covariance of (A,B) |
| retained interface | complete boundary holonomy G and every irreducible character |
| boundary pairing | L² of SU(2) class functions with normalized Haar measure |
| native transport | the same cross-corner L, positive return metric and half-share connection as UP8 |
| symmetry / gluing map | link gauge action plus character convolution over the common boundary |
| target | complete configuration source return, centre response and native curvature |

The ordinary response is strictly positive at every finite real source.
At Haar, Cov(A,B)=diag(3/4,3/4). Any affine relation between A and B
almost surely under a positive tilted density would also hold under
Haar, contradicting this covariance. The integrand is bounded on the
compact carrier, so its logarithm is real analytic on real sources.

## 2. Exact boundary memory, including every character

Let n=2j+1 be representation dimension, chi_n its character, and

\[
 f_n(t)=\frac{2I_n(t)}t,
 \quad f_1(0)=1,\quad f_{n>1}(0)=0,\quad
 u_n(t)=f_n(t)/f_1(t).
\]

The normalization is

\[
 e^{t\operatorname{Tr}G/2}
    =\sum_{n\ge1}n f_n(t)\chi_n(G),\qquad
 \int\chi_n\chi_m\,dG=\delta_{nm}.
\]

The SU(2) character convolution contributes 1/n for each joined edge.
Integrating the interior of each three-face disk gives

\[
 K_L(G;k)=\sum_{n\ge1}n f_n(k)^3\chi_n(G),
 \qquad K_R(G;h)=\sum_{n\ge1}n f_n(h)^3\chi_n(G).
\]

Consequently CZ1's closed-surface formula, now with two source values, is

\[
 \boxed{Z(k,h)=\int K_L(G;k)K_R(G^{-1};h)dG
       =\sum_{n\ge1}n^2 f_n(k)^3 f_n(h)^3}.                \tag{2}
\]

This formula is reused, not claimed as a new partition-function solution.
CZ1's exact gluing engine gives exponent -4 for six equal face characters;
the six n factors in the weights leave n². Classical character gluing is
also described by Witten, *Two Dimensional Gauge Theories Revisited*,
section 4.1, especially the character basis and boundary factorization.
The lattice normalization in (2) is derived here, not inferred from a
continuum heat-kernel convention.

Normalize the disk densities by their integrals:

\[
 \kappa_L=K_L/f_1(k)^3,\quad \kappa_R=K_R/f_1(h)^3,
 \quad v_L=\kappa_L-1,\quad v_R=\kappa_R-1.
\]

**C1 — complete boundary return.** With
Z_ind=f_1(k)^3 f_1(h)^3, the exact correction is

\[
 \boxed{\frac{Z}{Z_{\rm ind}}=1+\mathcal I(k,h)},\qquad
 \boxed{\mathcal I=\langle v_L,v_R\rangle
       =\sum_{n\ge2}n^2u_n(k)^3u_n(h)^3}.                 \tag{3}
\]

Every character contributes its own source product. For k,h>0 the return
is strictly positive. It is a configuration-kernel pairing, not the
Schur resolvent of the electric Hamiltonian. The kernel itself retains
the common boundary variable; the scalar I is its contracted reading.

**C2 — full-channel bounds.** Positive Bessel series show, for t>=0,

\[
 0\le u_n(t)\le\frac{(t/2)^{n-1}}{n!}.
\]

Proof: divide the series for I_n by that for I_1. After extracting
(t/2)^(n-1), each numerator/denominator coefficient ratio is
(r+1)!/(r+n)!<=1/n!. Thus the quotient is a weighted mean bounded by
1/n!, for all finite t. No large-order asymptotic is used.

Set z=kh/4 and a_n=n² z^(3(n-1))/(n!)^6. For a retained dimension D,

\[
 \boxed{0\le\sum_{n>D}n^2u_n(k)^3u_n(h)^3
 \le\frac{a_{D+1}}{1-z^3/[(D+1)^2(D+2)^4]}}             \tag{4}
\]

whenever the denominator is positive. Indeed
a_(n+1)/a_n=z³/[n²(n+1)^4], decreasing in n. In particular

\[
 \mathcal I\le\frac{(kh)^3/1024}{1-(kh/4)^3/324}.
                                                               \tag{5}
\]

This sufficient whole-return bound assumes (kh/4)^3<324. It does not
assert divergence when that condition fails. For each disk, |chi_n|<=n
similarly gives the useful uniform source bound

\[
 \boxed{\|\kappa_L-1\|_\infty
   \le\delta(k)=\frac{k^3}{16(1-k^3/96)}}\quad(k^3<96).   \tag{6}
\]

The right disk has delta(h). Hence also I<=delta(k)delta(h). The
log-potential correction satisfies I/(1+I)<=log(1+I)<=I. These bounds
include all representations, with no finite-span invariance assumption.

## 3. Which information is invisible before the cube closes?

**C3 — proper-subset independence and a six-face cumulant.** More
generally CZ1's gluing yields

\[
 Z(t_1,\ldots,t_6)=\sum_{n\ge1}n^2\prod_{p=1}^6 f_n(t_p).
                                                               \tag{7}
\]

If any t_p is zero, all n>1 terms vanish and the remaining integral
factorizes into the five or fewer single-face functions. Since these
bounded variables have analytic generating functions, **every proper
subset of the six face traces is jointly independent under Haar**.
Every nonempty proper face subset also has a boundary link; the exact
incidence test checks this for all 62 such subsets.

But f_2(t)=t/4+O(t³), so

\[
 \log Z-\sum_p\log f_1(t_p)
   =\frac{t_1t_2t_3t_4t_5t_6}{1024}+O(|t|^8).
\]

The joint six-face cumulant and moment are both 1/1024. In the grouped
channels, <A³B³>=9/256 while <A³><B³>=0. Lower mixed moments with either
power below three factorize. A pair-covariance observer can therefore
be completely blind to this closed record.

This qualifies the broad shared-link wording at UP9's handoff. Two
ordinary open faces can have independent Haar-distributed traces even
though they share a link. YC2's 1/64 covariance concerned repeated-link
commutator plaquettes on the one-site torus, a different carrier. Under
the complete cube tilt the traces need not remain independent.
Furthermore static Haar independence never removes the kinetic join:
YC8 already has nonzero pointwise shared-edge Gamma(W_p,W_q).

## 4. The source selects a flat limit and a correlated centre response

Define, using (1),

\[
 U_\lambda=\lambda^{-2}\log Z(\lambda S,\lambda V),\qquad
 U_0=3(S^2+V^2)/8.                                     \tag{8}
\]

Both source means vanish at Haar, so this is an analytic removable limit.
The exact UP9 source-dilation identity still holds, now with the complete
cube source:

\[
 \boxed{U_{\rm mid}=U-\tfrac12(SU_S+VU_V)
                    =-\tfrac\lambda2\partial_\lambda U}.
                                                               \tag{9}
\]

Let U_ind=lambda^-2 log Z_ind(lambda S,lambda V) and delta U=U-U_ind.
Equations (2)–(3) give the complete correlation correction

\[
 \boxed{\delta U=\lambda^{-2}\log(1+\mathcal I)},\qquad
 \boxed{\delta U_{\rm mid}
       =-\tfrac\lambda2\partial_\lambda\delta U}.         \tag{10}
\]

Its centre reading need not be positive. It is the signed degree-filter
response, not an assumed matter-energy density.

Write epsilon=lambda². Exact moments from (2) give

\[
\begin{aligned}
 U_\lambda={}&\frac38(S^2+V^2)
 -\frac{\lambda^2}{128}(S^4+V^4)\\
 &+\lambda^4\left[\frac{S^6+V^6}{3072}
                         +\frac{S^3V^3}{1024}\right]\\
 &-\lambda^6\left[\frac{S^8+V^8}{61440}
               +\frac{S^3V^3(S^2+V^2)}{8192}\right]
 +O(\lambda^8).                                        \tag{11}
\end{aligned}
\]

In particular delta U_mid=−lambda^4 S³V³/512+O(lambda^6).
The first closed-cube record contributes at order four after the
lambda^-2 normalization. Omitting I leaves a separable potential and
an identically closed selected compass, at every source value.

## 5. Source-to-compass transport: why order twelve cancels

Use UP8's unchanged definitions:

\[
 D_1=S\partial_S,\quad D_2=V\partial_V+mD_1,
 \quad m=-VU_{SV}/(SU_{SS}),
\]
\[
 L=\begin{pmatrix}D_1^2U&D_2D_1U\\D_1D_2U&D_2^2U\end{pmatrix}
   =aI+pK+qRK+wR,\quad\Delta=p^2+q^2-w^2,
\]
\[
 \ell=w^2/\Delta,\quad\phi=\arg(p+iq),\quad
 \mathfrak F=R\,d\ell\wedge d\phi.                      \tag{12}
\]

On the declared patch S>V>0 the zero limit has
p0=3(S²−V²)/4>0, q0=w0=0, B0=I. The line S=V is a traceless-cut
degeneracy, although the source Hessian is regular there. Put d=S²−V².
Exact division and differentiation of (11) yield

\[
\begin{aligned}
 w={}&-\lambda^4\frac{9S^3V^3}{2048}
   +\lambda^6\frac{3S^3V^3(7S^2+5V^2)}{16384}
   +O(\lambda^8),\\
 \ell={}&\lambda^8\frac{9S^6V^6}{262144d^2}
  -\lambda^{10}\frac{3S^6V^6(5S^2+3V^2)}{1048576d^2}
  +O(\lambda^{12}),\\
 \phi={}&-\lambda^4\frac{9S^3V^3}{512d}
    +\lambda^6\frac{S^3V^3(7S^2+9V^2)}{4096d}
    +O(\lambda^8).
\end{aligned}                                                \tag{13}
\]

**C4 — delayed curvature.** At leading order ell8=t² and phi4=3t,
where t=w4/p0. Hence d ell8 wedge d phi4=0. A nonzero cut defect
does not by itself give nonzero transport curvature at that order.
Retaining the next source shape gives

\[
 \boxed{\mathfrak F_{SV}
  =\lambda^{14}\frac{9S^{10}V^8(3S^2-V^2)}
                         {33554432(S^2-V^2)^4}R
     +O(\lambda^{16}).}                                 \tag{14}
\]

The coefficient is strictly positive on S>V>0; at S=2,V=1 it is exactly
11/294912. Analyticity and Delta0>0 prove positive curvature for all
sufficiently small positive lambda, uniformly on each compact subset
of that patch. This is a local theorem; no global sign or monotonicity
claim in lambda is made.

The return metric and transport remain the declared real two-component
objects of UP8. CZ1's commuting scalar centre-weight operators have
zero commutator; (14) is the curvature of a **different**, response-derived
matrix connection. They do not contradict each other and are not
silently identified with the full SU(2) gauge curvature.

## 6. Full-integral certificates at finite sources

The executable does not substitute (11) into a finite-source calculation.
It computes all needed derivatives of the complete twelve-link integral.
For each monomial, (2) has only finitely many contributing dimensions:

\[
 [k^a h^b]Z=\sum_{1\le n\le1+\min(a,b)/3}
        n^2[k^a]f_n(k)^3[h^b]f_n(h)^3.
\]

Multiplication by a!b! gives the exact rational moment <A^a B^b>.
At S=2,V=1, enclose <A^i B^j exp(lambda(2A+B))>, i+j<=4, by expanding
the exponential through degree N=40 and adding the entire tail

\[
 |\mathrm{tail}|\le
 3^{i+j+\lceil t\rceil}\frac{t^{N+1}}{(N+1)!},\qquad
 t=9|\lambda|,                                         \tag{15}
\]

using |A|,|B|<=3 and e^t<=3^ceil(t). The degree-four state jet of log Z
determines L and its first derivatives. All arithmetic rounds outward
to multiples of 2^-160, reusing UP9's tested jet transport.

An independent computation of I retains dimensions through six,
encloses each positive Bessel series with its complete tail, and uses
(4) for every higher representation. Its interval agrees with the
moment computation of Z/Z_ind−1. A bounded alternating logarithm series
and exact source-mean derivatives give delta U_mid from (10).

| lambda, S=2,V=1 | certified F_SV/R | certified delta U_mid |
|---|---|---|
| 1/4 | [1.205337,1.205338] times 10^-13 | [-0.0000575901,-0.0000575900] |
| 1/2 | [1.283772,1.283773] times 10^-9 | [-0.000778350,-0.000778349] |

These are outward rational bounds; the JSON retains much narrower
endpoints. The same packet certifies positive ordinary Hessians, positive
cut discriminants, negative w, positive mixed covariance and 0<chi<1.
For scale, at lambda=1/2 the full boundary return is about 0.0001049652,
below the all-channel bound 0.0001220711; these last decimals are only
readable summaries of the rational intervals.

## 7. What this supplies to the next Yang–Mills step

The adapter now retains a correlated block's complete configuration
boundary kernel, the shared representation label, exact source moments,
norm/tail budgets and native curvature. It demonstrates why low-order
pair response is insufficient: the first collective source can require
an entire closed surface.

Source rescaling still gives U_(lambda/sqrt(b))(x)=bU_lambda(x/sqrt(b))
and the two-form Jacobian b^-1. It is an independent-copy law only for
independent **whole cube records**, not for six independently substituted
faces. Spatial gluing of two cubes must retain their common face data;
the completely integrated scalar (2) cannot be reused as an independent
block weight without that additional construction.

Most importantly, Z is a configuration integral. Its covariance and
boundary pairing are not the actual quantum vacuum measure or the
resolvent of H0+theta V. The original electric operator couples shared
links even where the Haar pair covariance vanishes. To use YC10's
correlated-reference criterion one still needs an actual block ground,
gap g, the boundary electric form and its complete source/inverse budget.
This stage supplies the configuration-side source and a rigorous error
budget; it does not discharge those kinetic conditions or improve the
volume-uniform quantum gap. The physical scale and continuum remain open.

## Evidence and classical credit

```sh
python yc11_closed_cube_response.py --check
python -m unittest test_yc11
```

Forty exact/outward checks and fourteen tests cover the cube incidence,
character normalization, closed-source coefficients, full integral signs,
boundary tails and independent representations of the same return.
Written gluing/analytic proofs and finite arithmetic certificates are
separate; this is not proof-assistant verification.

CZ1 already supplies closed-surface character gluing. Standard ingredients
are Schur orthogonality / Peter–Weyl characters, the Bessel expansion,
cumulant generating functions and positive-series remainder bounds.
The new result is their complete joint-source and native-transport adapter,
with the order-twelve cancellation, order-fourteen curvature and finite
full-integral certificates. No new general character-gluing theorem is claimed.

External primary reference: E. Witten, *Two Dimensional Gauge Theories
Revisited* (1992), section 4.1, <https://arxiv.org/abs/hep-th/9204083>.
Native sources and the new note/code/tests are pinned by SHA-256 in
`YC11_RESULT.json`.
