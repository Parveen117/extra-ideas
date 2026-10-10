# YC29 — shared-block feedback is a neutral response Gram matrix

10 October 2026. Research owner: Monty Dabas. Continues YC28 at `ae1e7ee`.
Unit-S³ electric normalization, normalized Haar, Wp=Tr(Up)/2. Keep both
YC27 partitions, every charged and neutral multiplicity, and their domains.
No predecessor packet is changed.

**Result.** Resolve the shared-block part of the repeated-face return
before counting face pairs. Its complete neutral contribution to the
fourth mixed-coupling energy coefficient is

\[
 -x^*\Gamma x,
 \qquad x_p=\xi_p^2,\qquad
 \Gamma_{pq}=\langle\nu_p,H_N^{-1}\nu_q\rangle,
 \qquad \nu_p=(P-|\Omega\rangle\langle\Omega|)R_p(0)\Omega. \tag{1}
\]

This is a positive Gram matrix on face labels, with its sign reversed in
the energy. It includes same-face returns and repeated pairs on common
neutral excitation supports. Sharing a factor alone does not determine
its size. The full charged contribution and the folded energy term remain
separate; (1) is not the whole fourth coefficient or a mass-gap bound.

There is a useful exact first internal-coupling jet. Write each pure
internal plaquette coupling as eta*rho_r, with 0<=rho_r<=1. At eta=0,
every nu_p vanishes. Its derivative is

\[
 \boxed{\nu'_p(0)=\frac1{22464}
     \sum_{r\ {\rm pure}:\ r\cap p\ {\rm is\ an\ edge}}
                  \rho_r W_r\Omega_0.}               \tag{2}
\]

Thus a large common block is replaced, at this jet, by the internal
plaquettes actually touching the selected edges. Let M_rp be this
edge-incidence matrix. An additional exact centre cut makes Gamma even
in eta, and

\[
 \boxed{\Gamma(\eta)=\frac{\eta^2}{48\cdot22464^2}
                M^*\operatorname{diag}(\rho_r^2)M
                +O(\eta^4).}                         \tag{3}
\]

The leading matrix has a volume-independent face-label norm bound with
geometry constant **62 for the28-link tiling and60 for the64-link tiling**.
These replace the growing whole-block incidence count for this coefficient
only. The remainder in(3) is a local analytic expansion near the free
internal point, not a numerical bound at eta=1. Neither the complete
boundary budget nor its spatial iteration has been proved to contract.

An exact nested return and its lifted norm are constructed first below.
The selected neutral second-order state has norm kernel whose leading
coefficient is (3)'s coefficient divided by12. This is one part of the
full metric; it does not replace the other norm terms.

YC27's lattice gap windows are unchanged. The next estimate needed is
the finite-internal-coupling remainder of this shared response, together
with the retained charged return and full excitation metric. No4D
continuum theorem is claimed.

## 1. Full operator carrier and two independent sign cuts

On YC27's complete globally gauge-invariant carrier construct

\[
 H(\eta,\lambda)=H_0(\eta)-\lambda V,\quad
 H_0(\eta)=\sum_f\left(K_f-\eta\sum_{r\subset f}\rho_rW_r-E_f(\eta)\right),
 \quad V=\sum_{p\ {\rm mixed}}\xi_pW_p.                \tag{4}
\]

The tensor factors f are the old blocks and complementary cube/square/link
components. Kf is their full electric energy, Ef their actual ground
energy. Carry the scalar sum Ef separately when reading the raw operator.
The inherited product-Laplacian domain and Haar pairing are unchanged;
all Wilson potentials are bounded smooth multipliers. The scalar eta and
lambda are coordinates in this represented operator family, not a proposed
definition of the universal primitive lambda-space or of RG time.

Let Omega(eta) be the product of the normalized positive factor grounds,
P the projection onto all factor singlets, and Pi=|Omega><Omega|. Set
N=P-Pi. These projections retain every neutral multiplicity; Pi is the
additional ground cut, not a replacement for P. At each finite lattice
and finite real eta, the factor ground is simple and isolated, so HN,
the restriction of H0 to N, has an inverse. Its lower bound at general
internal coupling is not asserted uniform in block size. In YC27's
windows HN>=4. Sources and their finite-block returns are real analytic
near eta=0: isolated-ground resolvent contours and Neumann resolvent
expansions give this directly for the bounded perturbations in(4).

YC28's J preserves H0 and reverses V. There is a second cut. For a
positively oriented link (v,i), i=0,1,2, apply its centre sign

\[
 \epsilon(v,i)=\sum_{j<i}v_j\pmod2.                    \tag{5}
\]

Call the resulting unitary C. Every side length is even, so(5) is periodic.
The sum of its four signs around any elementary plaquette is1: C reverses
every W, preserves the electric operators, and commutes with the gauge
action and factor singlet projections. Set D=CJ. The centre operations
commute, so D is an involution, reversing every pure plaquette and leaving
every mixed plaquette fixed. Restricted to each factor it carries the
positive ground at eta to the positive ground at -eta. Hence

\[
 JH(\eta,\lambda)J=H(\eta,-\lambda),\qquad
 DH(\eta,\lambda)D=H(-\eta,\lambda),\qquad
 D\Omega(\eta)=\Omega(-\eta).                         \tag{6}
\]

Each Ef is even in eta. Negative eta is a legitimate bounded-potential
continuation used for this identity, not an assertion that the physical
coupling profile is negative. These cuts are link/factor operations, not
global gauge transformations that would act trivially on the carrier.

## 2. A complete nested return, including its norm

Write Pplus/minus=(I+/-J)/2, Hminus=Pminus H0 Pminus,
T=Pminus V Pplus. YC28 gives Hminus>=4 at finite nonnegative internal
couplings. For real z<4 put

\[
 B(z)=T^*(H_--z)^{-1}T,\quad
 Q_e=P_+-\Pi,\quad H_e=Q_eH_0Q_e,
\]
\[
 a(z)=\langle\Omega,B(z)\Omega\rangle,\quad
 b(z)=Q_eB(z)\Omega,\quad
 D_\lambda(z)=H_e-z-\lambda^2Q_eB(z)Q_e.              \tag{7}
\]

The symbol D_lambda in(7) is an operator pencil, distinct from the centre
involution D in(6). Whenever D_lambda(z) is invertible, two exact block
eliminations give the scalar ground-cut pencil

\[
 \boxed{F_0(\lambda,z)=-z-\lambda^2a(z)
          -\lambda^4\langle b(z),D_\lambda(z)^{-1}b(z)\rangle.} \tag{8}
\]

Such a neighborhood of (lambda,z)=(0,0) exists at each finite volume:
He has a positive gap there and B is bounded. No volume-independent
inverse norm or new mixed-coupling disk is inferred from this fact.
For real lambda,z, the lifted even vector and its odd companion are

\[
 f_+=\Omega+\lambda^2D_\lambda(z)^{-1}b(z),\qquad
 f_-=\lambda(H_--z)^{-1}Tf_+ .                        \tag{9}
\]

Their exact squared norm, the ground-cut metric, is

\[
 \boxed{-\partial_zF_0=
 1+\lambda^4\|D_\lambda(z)^{-1}b(z)\|^2
 +\lambda^2\|(H_--z)^{-1}T
          [\Omega+\lambda^2D_\lambda(z)^{-1}b(z)]\|^2.} \tag{10}
\]

This follows either by substituting(9) in the original pairing or by
differentiating the actual Schur identity. Energy dependence in both b
and D_lambda is retained. It is not legitimate to keep just a neutral
inverse-square term as the whole metric.

The even analytic ground energy above the carried reference is therefore

\[
 E(\lambda)=-a(0)\lambda^2+
 \left[a(0)a'(0)-\langle b(0),H_e^{-1}b(0)\rangle\right]\lambda^4
 +O(\lambda^6).                                      \tag{11}
\]

For clarity, expanding a(E)=a(0)+E a'(0)+... in F0=0 supplies the
positive folded term a(0)a'(0). It is essential for disconnected-cluster
cancellation. Omitting it would turn disconnected raw returns into
spurious energy interactions.

## 3. Gather the repeated neutral response before taking its norm

Recall Rp(z)=P Wp Gz Wp P from YC28. The charged inverse Gz acts on
I-P; on WpP its value is the odd inverse in(7). Different initial faces
have distinct corner-charge patterns, so

\[
 PB(z)P=\sum_p\xi_p^2R_p(z),\quad
 a(0)=\sum_p\xi_p^2\chi_p,\quad
 a'(0)=\sum_p\xi_p^2\langle\Omega,W_pG_0^2W_p\Omega\rangle. \tag{12}
\]

Decompose Qe=N+C_even, where C_even=Pplus-P is the complete remaining
even charged carrier. Both summands reduce H0. Then

\[
 Nb(0)=\sum_p\xi_p^2\nu_p,\qquad
 \langle b,H_e^{-1}b\rangle
 =\left\|H_N^{-1/2}\sum_p\xi_p^2\nu_p\right\|^2
   +\|H_{C_{\rm even}}^{-1/2}C_{\rm even}b\|^2.       \tag{13}
\]

This proves(1), exactly at each admitted finite internal profile. The
diagonal entry contributes to p^4 and the off-diagonal entries to p²q².
The charged term in(13) still contains repeated words and distinct tubes;
the folded term in(11) still contains its required cancellations.

There is a finer exact support statement. For a nonempty factor set I,
let EI excite precisely those neutral factors, using PG,f-Pi_f on I and
Pi_f elsewhere. Write nu_(p,I)=EI nu_p. If Xp is p's four-factor support,
then nu_(p,I)=0 unless I is a subset of Xp. Indeed Wp acts only on Xp,
and spectator reference grounds have zero energy in its exact resolvent.
The full resolvent on excited spectators is not replaced by I; this
argument is specifically for its source on Omega. Therefore

\[
 \boxed{\Gamma_{pq}=
 \sum_{\varnothing\ne I\subset X_p\cap X_q}
     \langle\nu_{p,I},H_I^{-1}\nu_{q,I}\rangle.}       \tag{14}
\]

H_I is the sum of the excited factor Hamiltonians on that support. If
Xp and Xq share exactly one factor, just that factor's complete neutral
response occurs in(14). If they share none, this particular contribution
vanishes. For overlapping factors it need not vanish when face corners
are disjoint. Formula(14) retains the multiplicities responsible for that
effect and shows exactly what needs bounding: a response Gram matrix,
not a count of all pairs that happen to meet a block.

The second-order neutral part of the intermediate-normalized actual
ground vector is HN^-1 sum xi_p² nu_p. Its squared norm is x* L x, where

\[
 L_{pq}=\langle\nu_p,H_N^{-2}\nu_q\rangle.            \tag{15}
\]

It enters the lifted norm at order lambda⁴. Equation(10) also contains
the odd source and its interference with higher returns, so(15) is not
the complete fourth-order coefficient of that norm.

## 4. Compute the first internal jet without a harmonic truncation

At eta=0 the reference ground is Omega0=1 and every Wp Omega0 has
electric energy12. Since P Wp²P=I/4,

\[
 R_p(0;0)\Omega_0=\Omega_0/48,\qquad \nu_p(0)=0.      \tag{16}
\]

Here the two arguments denote spectral z and internal eta. This is only
a statement about the neutral source on the vacuum. Rp remains an
operator on excited neutral states, and other fourth-order channels do
not vanish at the free point.

Perturb one pure plaquette r with potential -eta*rho_r Wr. The exact
first ground derivative is rho_r Wr Omega0/12, and E'_f(0)=0. Distinct
pure plaquette states are orthogonal: a link in their symmetric difference
has a single fundamental matrix element whose Haar mean is zero. Each
has norm squared1/4 and energy12.

A pure plaquette can share at most one edge with a mixed plaquette,
because the latter uses only one link from each factor. First consider
one shared edge, written as a unit quaternion with coordinates u_alpha.
On that edge Wr is a degree-one spherical harmonic. Multiplication by
u_alpha splits it into degree0 and degree2 pieces. The exact contractions
are

\[
 \sum_\alpha u_\alpha\Pi_0(u_\alpha W_r)=W_r/4,
 \qquad
 \sum_\alpha u_\alpha\Pi_2(u_\alpha W_r)=3W_r/4.       \tag{17}
\]

The first follows from integral u_alpha u_beta=delta_alpha,beta/4;
the second follows by subtraction from sum u_alpha²=1. There are no
other degrees in this product. The remaining six distinct links in
Wp Wr each contribute3, while the shared edge contributes0 or8.
Using the other three selected-link Haar contractions consequently gives

\[
 R_p(0;0)W_r\Omega_0=A_{pr}W_r\Omega_0,\qquad
 A_{pr}=\frac1{16}\left(\frac1{18}+\frac3{26}\right)
       =\frac5{468}                                  \tag{18}
\]

when r shares an edge. If r shares no edge, the electric energies simply
add to24, so A_pr=1/96, including when r and p belong to a common large
factor but use different links.

Differentiate the inverse defining Rp. The projection P is independent
of eta; the ground energy has zero first derivative. The inverse
derivative contributes A_pr/12, and differentiating the source ground
contributes another A_pr/12. Differentiating the subtracted ground line
contributes -1/576. The scalar susceptibility has zero first derivative:
all the remaining vectors are orthogonal to Omega0. Thus the coefficient
of rho_r Wr Omega0 in nu'_p is

\[
 \frac{2A_{pr}}{12}-\frac1{576}
 =\begin{cases}1/22464,&r\text{ shares an edge with }p,\\
                0,&r\text{ shares no edge with }p.
   \end{cases}                                      \tag{19}
\]

This proves(2), including the cancellation for an internal plaquette
elsewhere in the same factor or on a spectator. Only the exact product
of two fundamental harmonics was decomposed in(17); no omitted-harmonic
approximation is used for this derivative.

By(6), nu_p(-eta)=D nu_p(eta) and HN(-eta) is the corresponding conjugate
of HN(eta). Hence Gamma and L are both even in eta. Combining(19) with
the orthogonal energy12 states and their norm1/4 proves(3) and

\[
 \boxed{L(\eta)=\frac{\eta^2}{576\cdot22464^2}
             M^*\operatorname{diag}(\rho_r^2)M+O(\eta^4).} \tag{20}
\]

The O terms describe convergent Taylor expansions near zero at each
finite volume. Their coefficients depend only on bounded unions of
the specified factors by(14), but this stage does not supply a numerical
remainder at the existing internal profiles. In particular eta=1 cannot
be substituted into the leading term as an upper bound.

## 5. Sum the local geometry instead of the shared-block face count

Every entry of M is0 or1. An old pure face meets at most8 mixed faces;
a complementary cube face meets8, and a complementary square face12.
There are no internal faces in a link component. Every mixed face meets
at most8 pure faces. Consequently the elementary row/column estimate
already gives ||M||²<=12*8=96, independently of volume or the number of
faces elsewhere on an old block.

For these two tilings the exact local count improves96. Here is a
finite residue calculation with an all-volume interpretation. Let
b_i(n)=1 if n mod a_i is0 or a_i-1, and0 otherwise. For a pure face
r=(v,i,j), with k the remaining coordinate, its mixed incidence d_r is

\[
 d_r=\begin{cases}
 b_i(v_i)+b_i(v_i+1)+b_j(v_j)+b_j(v_j+1)+4b_k(v_k),
                                      &r\text{ old},\\
 12-4b_k(v_k),                         &r\text{ complementary}.
 \end{cases}                                             \tag{21}
\]

To prove the first line, each old edge sees one mixed plaquette for
each perpendicular coordinate on an old boundary. Sum over the four
edges of r. In the second line each complementary edge sees3-b_k
mixed faces; sum its four edges. No mixed face shares two edges with
r, so these are counts of distinct faces.

For a mixed p, enumerate the four plaquettes incident on each of its
four edges and keep the pure ones. Their set is N(p), of size at most8.
The row sum of M*M at p is exactly

\[
 s_p=\sum_{r\in N(p)}d_r.                              \tag{22}
\]

Equation(21) depends only on the residues of p's origin in one old
cell. The resulting exhaustive values are:

| Vertex cell | Values s_p: number of mixed faces per cell | Maximum |
|---|---|---:|
| 4×2×2 | 48:4, 56:8, 58:4, 62:8 | 62 |
| 4×4×2 | 22:4, 40:8, 46:8, 54:16, 56:4, 60:4 | 60 |

The checker evaluates(21) on the unwrapped integer lattice; it does not
extrapolate a finite torus spectrum. Repeating the cell gives every
admitted volume. The shortest periodic side is4, which does not identify
the distinct elementary faces around an edge in these counts. Repeated
second-neighbor paths to one face are intentionally counted with their
M*M multiplicity. Thus the same row sums hold on the shortest tori too.
The nonnegative symmetric matrix M*M has norm bounded by its largest
row sum. Multiplication by diag(rho²)<=I can only decrease the positive
Gram form, giving the all-volume bounds

\[
 \boxed{0\preceq\Gamma_2
 \preceq\frac{B_a}{48\cdot22464^2}I,\qquad
 0\preceq L_2\preceq\frac{B_a}{576\cdot22464^2}I,
 \quad B_{422}=62,\ B_{442}=60.}                     \tag{23}
\]

Here Gamma2 and L2 are the coefficients of eta², not values of the
full interacting matrices. In particular the leading selected energy
loss obeys

\[
 0\le x^*\Gamma_2x
 \le\frac{B_a}{48\cdot22464^2}\sum_p\xi_p^4.          \tag{24}
\]

For one distinct pair p,q, its coefficient of lambda⁴ eta² in the
selected neutral energy is explicitly

\[
 -\frac{\xi_p^2\xi_q^2}{24\cdot22464^2}
      \sum_{r\in N(p)\cap N(q)}\rho_r^2.              \tag{25}
\]

It is zero unless a single pure plaquette touches both faces. This is
stronger, at this order, than merely requiring a shared factor. It also
exhibits the retained internal interaction that creates the response.
Higher internal orders can propagate through more plaquettes and must
not be declared zero on the strength of(25).

## 6. What was recovered, checked, and left open

At lambda=0, the full family(4) is exactly the product of its correlated
factors; (8) becomes -z and (10) becomes1 on its ground line. The larger
even carrier recovers Hplus-z with identity graph metric, as in YC28.
At eta=0, those factor Hamiltonians become the full free-link operators,
nu_p=0 and Gamma=L=0. The nonzero source derivative(2), charged return,
and other mixed interactions remain recorded. The two scalar evaluations
commute on the raw family; restore the carried Ef when comparing raw
energies. These are typed zero specializations of the constructed
operators, states, and original pairing, not discarded histories.

```sh
python physics/yc29/yc29_shared_response.py
```

Only the new exact harmonic fractions and the new residue counts are
checked. The calculator verifies the primitive edge-sign curl for the
additional cut, enumerates neighboring faces on the integer lattice,
and compares those counts with(21). It checks the corresponding pure
degrees, at most8 neighbors per mixed face, the row-sum table, and the
coefficient/metric arithmetic. The nested return, support, parity and
analytic-jet proofs are written above. No generic tests, predecessor
suite, new spectral truncation, or frozen-certificate regeneration was run.

The surviving task is now precise: bound the full interacting Gram
matrix and its remainder at the chosen internal profiles, include the
charged part of(13) and the folded term of(11), and control the full
excitation operator with its transported metric. A small negative vacuum
energy contribution is not an excitation gap. Constants62 and60 are
geometry bounds for one leading response coefficient; their decrease
does not establish a contracting spatial map. YC27's lattice windows
and the open continuum obligations remain unchanged.

Sources at `ae1e7ee`: YC27 for the full correlated factor carrier and
admitted lattice windows; YC28 for J, the charged inverse, Rp, source
orthogonality and repeated-face selection; YC19 for the separate full
excitation/metric obligation. The new nested ground-cut formula,
internal sign cut, shared-neutral Gram resolution, exact jet and local
incidence bounds are derived here. No independent peer or formal
proof-assistant verification is claimed.
