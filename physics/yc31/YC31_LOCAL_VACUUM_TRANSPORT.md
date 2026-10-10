# YC31 — local vacuum transport and the centered source return

10 October 2026. Research owner: Monty Dabas. Continues YC30 at `f866faa`.
Unit-S³ electric normalization. The YC27 internal profiles, all physical
centre sectors and complete link harmonics remain. No frozen packet changes.

**Result.** Within the existing YC27 windows, the reference vacuum can be
transported to the actual vacuum by a unitary whose action on a bounded
local observable has spatial tails smaller than every inverse power of
distance, uniformly in lattice volume. This supplies an isometric,
quasi-local alternative to assuming that the square root of YC19's global
quotient metric is local. It does not set that older metric equal to I.

There is a second, separate result. For a bounded local gauge-invariant
source B, YC30's centered return has the representation

\[
 K_\tau Q B\Psi=\mathcal R_\tau(B)\Psi,\qquad
 \mathcal R_\tau(B)=\int_{\mathbb R} k_\tau(t)\,
              e^{itH}B e^{-itH}\,dt.                 \tag{1}
\]

The filtered observable on the right admits a local approximation with
an explicit tail bound. The complete physical norm and source Gram errors
are controlled. The approximation is on a specified local source class;
it is **not** a small operator-norm approximation of the entire global K.

Both results use an already established gap. The numerical lattice
windows are unchanged. This is progress on locality of the centered
construction, not a new gap proof, a contracting blocking map, or a
continuum construction.

## 1. Carrier, path, and the full-space gap

Use YC27's partition of links into old factors O and complementary factors
D. Each site of the present factor graph carries its full Hilbert space
\(\mathcal H_f=L^2(SU(2)^{E_f})\). The physical carrier is the globally
gauge-invariant subspace of their tensor product, not a tensor product of
factor singlets. Start with the represented operator family

\[
 H(s)=\sum_f h_f-s\ell\sum_{p\ {\rm mixed}}\xi_pW_p,
 \quad 0\le s\le1,\quad 0\le\xi_p\le1,\quad
 0\le\ell\le\ell_*.                                \tag{2}
\]

Here h_f is the complete correlated factor Hamiltonian above its ground,
with the inherited elliptic domain. Mixed Wilson multipliers are bounded
and self-adjoint, of norm at most one. Finite-volume H(s) has the same
domain as H(0). Keep its unique positive normalized ground Psi_s, actual
energy E_s and projection Pi_s. The source is H'(s), not E_s' alone.

The path s is an ordered coupling adapter in the full operator carrier.
It is not primitive lambda-space, physical time or a renormalization scale.
The auxiliary t in (1) has inverse-electric-energy units. Operators,
adjoints, source, complementary projections and the Haar pairing precede
either parameter evaluation. Path segments compose in their stated order.

| Old factor | Mixed cap ell_* | Physical gap m | Face incidence I |
|---|---:|---:|---:|
| 28 links | 1/3200 | 592/175 | 48 |
| 64 links | 1/5700 | 11188/3325 | 88 |

Old and complementary internal couplings and admitted torus sizes are
exactly those of YC27. Its physical gap is uniform along (2). To apply
locality on the tensor carrier one must also check the nonphysical part:
YC22's ground-form inequality gives H(s)-E_s >= 3/6 = 1/2 there, at
every finite real Wilson coupling. Therefore the **full tensor-space**
gap is at least 1/2. This is not inferred from the physical gap alone.
Choose a spectral filter threshold gamma=1/4, strictly below this floor.

## 2. A propagation bound that retains the unbounded electric operator

Join two factors in the graph when they occur in a common mixed face.
Each face has four factors and graph diameter one; a factor touches at
most I faces. Let d be this graph distance and X(r) its closed r-ball
about a support X. Use all possible mixed faces, even when their xi=0.

The interaction picture for sum h_f preserves every support and norm.
Consequently only the bounded mixed interactions enter the commutator
iteration. If A and B have disjoint supports X and Y, an interaction
path reaching Y needs at least d(X,Y) steps. The first face has at most
|X|I choices and each subsequent overlapping face at most 4I choices.
Iterating the norm commutator inequality gives the deliberately loose bound

\[
 \|[\tau_t^s(A),B]\|
 \le 2\|A\|\|B\||X|
       \sum_{n\ge d(X,Y)}\frac{(8I\ell |t|)^n}{n!}
 \le 2\|A\|\|B\||X|e^{-d(X,Y)+v|t|},
 \quad v=8eI\ell .                                 \tag{3}
\]

This remains valid for arbitrary bounded B on the complement of a ball;
there is no factor |Y| in the last bound. No electric-energy cutoff has
been inserted. Strong integrals, rather than an assumption of norm
continuity of the unbounded on-site dynamics, are used below.

Let E_Y be partial expectation outside Y in the normal product state
rho=product_f |Omega_f><Omega_f|. It is a norm-one positive map fixing
operators on Y. This is an auxiliary localization map, not a replacement
of the actual interacting vacuum by a product state. The
commutator-to-localization inequality, with its
factor two for a general normal product state, yields

\[
 \|\tau_t^s(A)-E_{X(r)}\tau_t^s(A)\|
 \le \|A\|\min\{2,\,4|X|e^{-r+v|t|}\}.             \tag{4}
\]

In particular for any weight w in L1, define L_w=integral |w| and
T_w(T)=integral_{|t|>T}|w(t)|. Splitting its integral at r/(2v) gives

\[
 \begin{split}
 \epsilon_w(r,k)&=\min\{2L_w,
             4kL_we^{-r/2}+2T_w(r/(2v))\},\quad r\ge1,\\
 \epsilon_w(0,k)&=2L_w,\\
 \left\|\int w(t)\tau_t^s(A)dt-
        E_{X(r)}\int w(t)\tau_t^s(A)dt\right\|
       &\le \|A\|\epsilon_w(r,|X|).                \tag{5}
 \end{split}
\]

At ell=0 the dynamics preserves X exactly and the actual localization
error is zero; handle this directly instead of dividing by v=0.
If w has finite absolute moments of all orders, (5) decays faster than
any power of r, for fixed admitted energy units and couplings.

The factor graph has polynomial, not exponential, ball growth despite
the crude degree bound 3I. Choose an anchor vertex in each factor. A graph
step moves anchors by at most seven original lattice units in each
coordinate, and at most two factors can have the same anchor. Hence a
safe common bound on either periodic tiling is

\[
 |B_d(f,r)|\le\kappa(r+1)^3,\qquad \kappa=6750.       \tag{6}
\]

Indeed 2(14r+1)^3 is already an upper bound by lifting a graph path to
the periodic cover. These constants are not optimized or physical speeds.
At the two coupling caps, v/e equals 3/25 and 176/1425 respectively.

## 3. Construct the exact vacuum transport, including its sign

For an explicit filter, let b_gamma(omega) be the real even smooth bump
exp[-1/(1-(2omega/gamma)^2)] on |omega|<gamma/2 and zero elsewhere.
Let g_gamma be its inverse Fourier transform and set
w_gamma(t)=|g_gamma(t)|²/integral |g_gamma|². This is an even nonnegative
Schwartz function, of integral one, with Fourier support in [-gamma,gamma].
Set, away from the immaterial value at zero,

\[
 W_\gamma(t)=\operatorname{sgn}(t)
                    \int_{|t|}^{\infty}w_\gamma(u)du .            \tag{7}
\]

W is real, odd, integrable and has every absolute time moment. Use the
explicit Fourier convention F_W(omega)=integral W(t)e^{i omega t}dt.
Since W'=delta_0-w in distributions,

\[
 F_W(\omega)=\frac{i(1-F_w(\omega))}{\omega}
             =\frac{i}{\omega}\quad(|\omega|\ge\gamma),
 \qquad F_W(0)=0.                                  \tag{8}
\]

Define the self-adjoint generator and its unitary flow by

\[
 D_s=\int W_\gamma(t)\tau_t^s(H'(s))dt,
 \qquad \partial_s U_s=iD_sU_s,\quad U_0=I.          \tag{9}
\]

The plus sign in (9) belongs to the Fourier convention (8). For an
excited eigenvector n and the ground 0, with omega=E_n-E_0,
(D_s)_{n0}=i(H')_{n0}/omega. Thus
\(i[D_s,\Pi_s]_{n0}=-(H')_{n0}/\omega=(\Pi_s')_{n0}\).
The adjoint block agrees and both diagonal blocks vanish. The same
identity follows from spectral integrals without choosing an eigenbasis:

\[
 \Pi_s'=i[D_s,\Pi_s],\qquad
 \Pi_s=U_s\Pi_0U_s^*,\qquad Q_s=U_sQ_0U_s^*.        \tag{10}
\]

The finite-volume generator is bounded and strongly continuous, which
suffices for its unitary propagator. Its full norm may be extensive;
locality is proved with an interaction norm, not that full norm.
Gauge and global centre symmetries commute with H, H' and D. Therefore
U preserves the physical carrier and every global centre sector.
The ground phase can be chosen so U_s Omega=Psi_s.

For completeness, locality of U follows quantitatively from (5).
Write D_s=sum_p d_p(s), where d_p is (9) with H' replaced by -ell xi_p Wp.
Resolve d_p into shells E_{Xp}d_p and
(E_{Xp(r)}-E_{Xp(r-1)})d_p, r>=1. Their norms are bounded by ell a_r,
where

\[
 a_0=L_W,\qquad
 a_r=\epsilon_W(r,4)+\epsilon_W(r-1,4)\quad(r\ge1).
\]

At most I kappa(r+1)^3 faces have a radius-r shell containing a given
factor. Its diameter is at most 2r+1. For any real n>4 the shell interaction
therefore has F_n-norm bounded uniformly in s and volume by the finite number

\[
 J_n=\ell I\kappa\sum_{r\ge0}(r+1)^3(2r+2)^n a_r,
 \qquad F_n(d)=(1+d)^{-n}.                          \tag{11}
\]

Here the F-norm is sup_{x,y} sum_{Z contains x,y} ||Phi(Z)||/F(d(x,y)).
Equation(6) gives S_n=sup_x sum_y F_n(d(x,y)) <=
kappa[1+1/(n-4)]. A convolution constant is
C_n=2^{n+1}kappa[1+1/(n-4)], by splitting an intermediate site according
to which endpoint is at distance at least d(x,y)/2. Applying the bounded
interaction commutator iteration to (9), then the same localization
inequality as (4), gives for either U_s A U_s* or U_s* A U_s:

\[
 \|U_s A U_s^*-E_{X(R)}(U_s A U_s^*)\|
 \le \|A\||X|\,
 \frac{4\kappa(e^{2C_nJ_n}-1)}{C_n(n-4)}
 (R+1)^{4-n},\quad 0\le s\le1.                    \tag{12}
\]

The reverse flow has the same interaction bound. Since n is arbitrary,
this proves the asserted faster-than-any-power locality. The constants
in (11)-(12) are explicit filter integrals/series, not optimized numerical
certificates. A strong gap does not make this deliberately coarse bound
a useful small-radius numerical estimate automatically.

The classical ingredients are Lieb–Robinson commutator bounds and
quasi-adiabatic spectral flow. For the unbounded-site setting and the
normal-state localizing map, see Nachtergaele–Sims–Young, *Quasi-Locality
Bounds … Part I*, [arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2),
sections3.2,4.1 and6. The phase-transport interpretation is also in
Bachmann–Michalakis–Nachtergaele–Sims,
[arXiv:1102.0842v2](https://arxiv.org/abs/1102.0842v2).
The carrier check, full-space floor, conservative constants and source
application here are for the declared YC27 adapter.

## 4. Cut and pairing transport; what a local approximation preserves

Restrict to the physical excitation carrier. In fixed reference coordinates,

\[
 \widehat A_s=Q_0U_s^*(H(s)-E_s)U_sQ_0,
 \qquad \widehat K_s=\tau(\widehat A_s+\tau)^{-1}.    \tag{13}
\]

The domain is the unitary pullback of the actual excitation domain;
preservation of the original electric domain by U is not assumed.
The pairing in these coordinates is exactly the original Hilbert pairing.
This differs from YC19's nonunitary creator coordinates, where G remains.

For any orthogonal reference excitation cut P0, its full source f0 and its complement,
use Ps=Us P0 Us*, fs=Us f0 and Qs-Ps. Functional calculus, Schur cuts and
their exact lifts for the bounded Fredholm family intertwine under U.
For an energy Schur map, also impose YC30's form/domain admissibility;
locality alone does not establish it for a proposed cut. On such cuts,
YC30's existing reserve and physical lift bounds carry over without an
extra metric condition number.
At tau=4 the bounds on ||K|| remain 175/323 and 3325/6122.

There are concrete local starting cuts. For example
P0=I-|Omega_f><Omega_f| on one complete factor, tensored with the identity
elsewhere, annihilates the product vacuum and commutes with global gauge
action. Its conjugate Ps is an exact projection in Qs and obeys (12).
For a local projection P0 on X, let Br=E_{X(r)} Ps and suppose
||Br-Ps||<=epsilon<1/2. Br is a positive contraction, but need not be a
projection. Its spectrum is within epsilon of {0,1}, so its local rounding

\[
 P_r=\mathbf1_{[1/2,1]}(B_r),\qquad
 \|P_r-P_s\|\le2\epsilon,\qquad
 \|P_r\Psi_s\|\le2\epsilon                         \tag{14}
\]

follows by adding ||Pr-Br||<=epsilon and ||Br-Ps||<=epsilon.
The product-state expectation respects gauge symmetry. However Pr need
not annihilate the actual vacuum exactly. Compressing by Qs repairs that
part but generally destroys the projection identity. An exact excitation
cut cannot silently be replaced by a locally rounded one in a Schur map.
Equation(14) states the remaining vacuum leakage, rather than discarding it.

## 5. Localize the actual centered resolvent on a declared source class

Now fix s, its physical gap floor m from the table, and tau>0. Let eta(r)
be a fixed C-infinity cutoff, zero for r<=m/2 and one for r>=m. One explicit
choice in the transition is h(x)/(h(x)+h(1-x)), x=2r/m-1,
h(x)=exp(-1/x) for x>0 and zero otherwise. Define

\[
 f_\tau(\omega)=\eta(|\omega|)\frac{\tau}{\tau+|\omega|},
 \qquad k_\tau(t)=\frac1{2\pi}\int_{\mathbb R}
                         f_\tau(\omega)e^{-it\omega}d\omega .   \tag{15}
\]

The last transform is initially an L2 transform, not an absolutely
convergent frequency integral at t=0. The even smooth f is in L2 and
f^(j) is in L1 for every j>=1. Plancherel on |t|<=1 and integration by
parts outside it establish

\[
 \begin{split}
 L_k&\le\|f\|_2/\sqrt\pi+\|f''\|_1/\pi<\infty,\\
 T_k(T)&\le\frac{\|f^{(j)}\|_1}
                {\pi(j-1)T^{j-1}},\qquad T\ge1,\ j\ge2 .       \tag{16}
 \end{split}
\]

Thus k has all absolute moments. Fourier inversion recovers f everywhere.
Since f(0)=0 and f(omega)=tau/(tau+omega) for omega>=m, functional calculus
on the **physical** spectrum gives f(H-E)=K_tau Q there. For a bounded
local gauge-invariant B, B Psi is physical and
e^{it(H-E)}B Psi= tau_t(B) Psi. This proves (1), including exact centering.
There is no assertion that (15) equals the same resolvent on the lower
nonphysical charged spectrum; those states are not sourced by B Psi.

Put R_r(B)=E_{X(r)} R(B) and center its source in the actual vacuum:

\[
 y=K_\tau Q B\Psi,\qquad
 y_r=Q R_r(B)\Psi
     =(R_r(B)-\langle\Psi,R_r(B)\Psi\rangle I)\Psi .
\]

This is a vector with a local observable representative; its expectation
still belongs to the actual state. Equations(5),(16) prove

\[
 \boxed{\|y-y_r\|\le e_B(r):=
            \|B\|\epsilon_{k_\tau}(r,|X|).}         \tag{17}
\]

No high harmonics were discarded. Subtracting the vacuum scalar in y_r
does not change its support and Q cannot increase the error.
On the exact vectors ||y||<=q||B||, q=tau/(m+tau). For a specified finite
source family B_i, retain the actual physical Gram matrix. With e_i from
(17), its selected return and inverse-square norm errors satisfy

\[
 \begin{split}
 |\langle Q B_i\Psi,y_j-y_{j,r}\rangle|
       &\le\|B_i\|e_j,\\
 |\langle y_i,y_j\rangle-\langle y_{i,r},y_{j,r}\rangle|
       &\le q\|B_i\|e_j+q\|B_j\|e_i+e_i e_j .       \tag{18}
 \end{split}
\]

The second line uses K_tau² and is a source norm kernel, not the entire
Schur derivative metric. For a linear combination the error bound is
sum_i |c_i|e_i. Consequently a number of sources growing with volume
can incur a growing cost. Equations(17)-(18) must not be inserted as an
unqualified full-input epsilon into YC30's global operator-norm theorem.

## 6. Zero recovery, verification, and the next gate

At ell=0, H(s)=H0, H'=0, D=0, U=I, Psi=Omega and E=0. Transported cuts,
sources and pairing recover their reference values. Equation(15) still
represents tau/(H0_exc+tau), and its local observable remains supported
on the original factors X. A zero mixed interaction does not erase the
internal dynamics or its nonzero response.

Retain the endpoint source jet as well. Differentiating at ell=0 gives

\[
 \left.\partial_\ell(U_1\Omega)\right|_0
       =H_{0,\mathrm{exc}}^{-1}\sum_p\xi_pW_p\Omega .             \tag{19}
\]

The source has zero vacuum component by YC28's grading. At the fully
free internal point each Wp Omega has energy12, so its coefficient is
1/12. This agrees with the ground equation and has the positive sign
appropriate to H=H0-ell V. Thus typed zero evaluation of the construction
and its retained source derivative agree with the reference construction.

```sh
python physics/yc31/yc31_local_transport.py
```

The direct check verifies the Fourier-convention sign through an exact
rotating two-level projection, its unitary pairing and sourced resolvent,
the zero ground-response sign, and the propagation/return fractions.
The unbounded-domain, all-volume locality and filter-tail statements have
written proofs above. No harmonic truncation, generic test suite or
predecessor certificate was run. The filter-series constants were not
numerically optimized.

**Next measurable task:** choose a spatially repeated retained source
space, control its Gram conditioning and vacuum leakage, and convert the
local errors to a volume-uniform bound in an appropriate interaction or
excitation norm. Then assess the actual reduced operator, not only a
finite list of its source vectors. Uniformly small errors for each local
reading alone do not establish that bound.

For compatible fixed-spacing infinite-lattice data the summable
interaction tails also permit the usual local-observable thermodynamic
limit of the flow. No norm limit of the extensive U or global K, or a
new continuum Hilbert space, follows from this observation. Refinement
still needs an energy conversion A_phys,n=c_n A_n and
tau_n=tau_phys/c_n, nontrivial convergence, and a physical gap along the
appropriate coupling trajectory. These have not been constructed here.
