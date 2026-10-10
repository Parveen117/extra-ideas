# YC36 — local coefficients on the complete count-restricted carrier

10 October 2026. Research owner: Monty Dabas. Continues YC35 at `c13e585`.
Keep the complete YC27 factors, admitted profiles, volumes and existing
physical gap floors. Frozen predecessor packets remain unchanged.

**Result.** YC35's conditional local-coefficient input can now be constructed.
Use the Hausdorff distance between complete excitation supports, rather
than just the nearest pair of excited factors. At each fixed reference
count j its balls contain a volume-independent number of supports.
Two-sided locality then gives a full-input operator-norm bound on P_j,
including all harmonics and all matched charged physical states.

Buffered unions of local boxes supply the coefficients without a box-gap
assumption or an interacting box ground. Symmetrization and an affine
correction preserve the exact old reserve delta=m/(m+tau). A finite
geometric expansion also makes the hidden elimination depend only on
local coefficient data. For integers 1<=k<=j and polynomial degree p>=0,
with tau>0, the resulting
reduced operator obeys

\[
 \boxed{\|S^{\mathrm{loc},p}_{k,j}-S_k\|
 \le \frac{\tau^2\mathcal B^2 k}{j+1}
       +\frac{\varepsilon_{\mathrm{loc}}}{\delta}
       +\frac{q^{p+3}}{\delta},\qquad
 \delta I\le S^{\mathrm{loc},p}_{k,j}\le I,\quad q=1-\delta.} \tag{1}
\]

Here B is YC35's inverse-count constant. All choices can be made independently
of total volume. The construction is not a numerical harmonic truncation:
its local coefficient operators remain infinite dimensional, its constants
are deliberately coarse, and no practical values of j,R,L,p are certified.
The existing lattice gap is an input, not a new enlarged window or a
continuum theorem.

## 1. Carrier, cuts and the bounded object to approximate

Begin with the full tensor product of complete YC27 factors and their local
ground vectors Omega_f. Put q_f=I-|Omega_f><Omega_f| and N=sum_f q_f.
The exact excitation-support decomposition, before taking gauge invariants,
is the orthogonal sum over all subsets S of factors. Taking global gauge
invariants gives H_S^G; its vectors may entangle charged factors with matched
vertex charges. Zero-dimensional physical summands may be omitted; the
vacuum support remains until the Q0 cut. No factor is first replaced by
its singlet subspace.

YC31's unitary U transports Omega=product_f Omega_f to the actual vacuum.
On the full tensor carrier define A_math=U*(H-E)U, and on the excited physical
carrier use Ahat=A_math|_(Q0 H_G)>=mI. Keep the actual unitary pullback domain
of H; the new cuts and functional calculus below are bounded and need no
additional energy-domain invariance assertion. The physical floors are
m=592/175 or11188/3325 in their respective profiles.

P_j=1_[1,j](N) restricts the complete excited physical carrier, and I_k is
the inclusion of P_k H_G. Write V_S for the inclusion of H_S^G. The bounded
family to approximate is

\[
 K=\tau(\widehat A+\tau)^{-1},\quad F=I-K,
 \quad 0\le K\le qI,\quad \delta I\le F\le I,
 \quad K_j=P_jKP_j,\quad F_j=I_j-K_j.                         \tag{2}
\]

These compressions are taken after the complete functional calculus. In
particular F_j is not identified with the normalization of P_j Ahat P_j.
The exact full-hidden-space return onto P_k is S_k as in YC35. Its j-cut
return is S_(k,j), the Schur complement of F_j onto P_k.

The primitive data are the full operator carrier, adjoint, composition,
gauge action, source and complementary cuts, ordered transport and Haar
pairing. The admitted scalar coupling path H(s)=H0-s ell sum_p xi_p Wp is
an adapter into those data. Count j, Hausdorff radius, filter time and
polynomial degree are distinct parameters, none automatically a physical
dimension, clock or renormalization scale. It suffices to write the endpoint
s=1; the same construction on an initial path segment gives every admitted s.

## 2. Complete support geometry has bounded degree at fixed count

For nonempty finite supports S,T in the factor graph define

\[
 d_H(S,T)=\max\{\max_{x\in S}d(x,T),\max_{y\in T}d(y,S)\}.
                                                               \tag{3}
\]

This is the usual Hausdorff metric. If d_H(S,T)<=r and |S|<=j, then
T is contained in S(r), the union of radius-r balls around S. YC31 gives
|S(r)|<=j kappa(r+1)^3 with kappa=6750. Hence, among supports of size1..j,
both a row and a column have at most

\[
 D_j(r)=\sum_{a=1}^j {j\kappa(r+1)^3\choose a}
          =O_j((r+1)^{3j})                                  \tag{4}
\]

neighbours at Hausdorff distance at most r. This bound does not count
internal harmonic dimensions. Arbitrarily far-separated factors within
one source support are allowed.

The nearest-pair distance min_(x in S,y in T)d(x,y) would fail here. For
S={0}, every T={0,x} has nearest-pair distance zero even when x is arbitrarily
far away; its row degree grows with the volume. Hausdorff distance records
that unmatched excitation. Both directed distances in (3) are needed.

## 3. Two-sided locality of the exact coefficient blocks

Use YC31's smooth even frequency function

\[
 f_\tau(\omega)=\chi(|\omega|)\frac{\tau}{\tau+|\omega|},
 \quad \chi=0\ (|\omega|\le m/2),\quad
 \chi=1\ (|\omega|\ge m),\quad
 f_\tau(\omega)=\int k_\tau(t)e^{it\omega}dt.                  \tag{5}
\]

Its time weight is integrable with all absolute moments, as proved in
YC31(15)-(16). Denote L_k=integral|k_tau| and T_k(t)=integral_(|u|>t)|k_tau(u)|.
On the physical carrier f_tau(A_math)=K Q0, so the vacuum is never inverted.
On lower nonphysical spectral values no resolvent identification is needed.

For phi in H_S^G, the bounded creator C_S(phi)=|phi><Omega_S| tensor I_(S^c)
has norm ||phi|| and retains its spectator identities. Put

\[
 \alpha(B)=UBU^*,\qquad
 \mathcal R(B)=\int k_\tau(t)\tau_t^H(B)dt,\qquad
 \mathcal D_S(\phi)=\alpha^{-1}\mathcal R\alpha(C_S(\phi)).
                                                               \tag{6}
\]

Since the actual vacuum is an eigenvector, K V_S phi=D_S(phi)Omega.
The latter vector is exactly orthogonal to Omega. Let E_Z be product-state
conditional expectation outside Z. Its action on Omega is orthogonal vector
projection onto excitations confined to Z.

Choose a YC31 flow exponent nu>4 and write
l_nu(r)=min{2,C_nu(r+1)^(4-nu)}, with C_nu the coefficient of YC31(12)
after removing ||B|| |support B|. It is not the convolution constant there.
For r>=3 put u=floor(r/3) and

\[
 h_j(r)=\min\{2L_k,\;
 2[L_k j\{1+\kappa(2u+1)^3\}l_\nu(u)
        +\epsilon_k(u,j\kappa(u+1)^3)]\}.                    \tag{7}
\]

Set h_j(0)=h_j(1)=h_j(2)=2L_k. Here
epsilon_w(r,n)=min{2L_w,4n L_w exp(-r/2)+2T_w(r/(2v))}, r>=1,
with epsilon_w(0,n)=2L_w and v=8e I_face ell. The face incidence I_face
is48 or88. All these symbols are YC31's filter/locality bounds.

Localize the inner transport in S(u), the physical filter in S(2u), and
the outer inverse transport in S(3u). Their errors are respectively
j L_k l_nu(u), epsilon_k(u,j kappa(u+1)^3), and
j kappa(2u+1)^3 L_k l_nu(u). Replacing the supported approximation by the
conditional expectation of the exact operator costs at most another
factor two. Thus

\[
 \|\mathcal D_S(\phi)-E_{S(r)}\mathcal D_S(\phi)\|
       \le h_j(r)\|\phi\|\quad (|S|\le j).                   \tag{8}
\]

With sufficiently high filter moments, h_j(r)=O_j(r^(7-nu)). If an excited
factor of T lies outside S(r-1), its support projection annihilates
E_(S(r-1))D_S(phi)Omega, hence ||K_(T,S)||<=h_j(r-1). If the other directed
distance in (3) attains r, apply this argument to K_(S,T) and use K*=K.
Therefore, for d_H(S,T)=r>=1,

\[
 \|K_{T,S}\|\le h_j(r-1).                                   \tag{9}
\]

Set K_j^(R)_(T,S)=K_(T,S) for d_H(S,T)<=R and zero otherwise. Its mask is
symmetric, so K_j^(R) is self-adjoint. Equations(4),(9) and the block-matrix
Schur test give the full retained-input norm bound

\[
 \boxed{\|K_j-K_j^{(R)}\|\le e_j(R),\qquad
 e_j(R)=\sum_{r>R}D_j(r)h_j(r-1).}                            \tag{10}
\]

Indeed both sums of block norms are bounded by that series; Cauchy--Schwarz
with these row/column sums proves the operator estimate on the direct sum,
even when each summand is infinite dimensional. This is not an entrywise
accuracy claim on a finite harmonic sample.

Take, for example, nu=3j+16. Then e_j(R)=O_j(R^-8). Constants depend strongly
on j and on the chosen filters, but not on total volume. There is no claim
of a bound uniform as j tends to infinity with R fixed.

## 4. Buffered unions of boxes supply the coefficients

For each source support S choose L>=max{4,R}, C=S(L) and B=S(2L). These
may be disconnected unions; their sizes obey
|C|<=j kappa(L+1)^3 and |B|<=j kappa(2L+1)^3. Include complete factors and
only complete mixed plaquettes:

\[
 H_B(s)=\sum_{f\in B}h_f-s\ell\sum_{X_p\subset B}\xi_pW_p,
 \quad J_C=-\ell\sum_{X_p\subset C}\xi_pW_p,
\]
\[
 D_{B,C}(s)=\int W_\gamma(t)\tau_t^{B,s}(J_C)dt,
 \quad U_{B,C}'=iD_{B,C}U_{B,C},\quad U_{B,C}(0)=I,
 \quad\gamma=1/4.                                           \tag{11}
\]

The finite-box electric domain is unchanged by its bounded potentials.
This flow exists without a box gap. It is not asserted to transport the
interacting box ground: its source J_C differs from H_B'. Gauge symmetry
and the inherited centre action remain because no plaquette was split.

Let alpha_BC be endpoint conjugation and R_B the same k_tau filter of
the box dynamics. Construct D_S^B=alpha_BC^-1 R_B alpha_BC(C_S). The only
state used for its coefficients is the known product Omega_B. The box may
contain all components of an entangled source without replacing that source
by a product of independent single-factor readings.

Here are explicit conservative bounds extending YC33 to these unions.
For d>=1 and any weight w define

\[
 b_w(d,n,M)=\min\{2L_w,\;
 2\ell I_{\rm face}M n M_1(w)e^{-d/2}+2T_w(d/(2v))\},        \tag{12}
\]

where M1(w)=integral |t w(t)|. Set b_w(0,n,M)=2L_w. It bounds the filtered
physical dynamics difference when the initial support has size n and
boundary distance at least d+1, and the box size is at most M. The proof
is bounded-boundary Duhamel, then a time split at d/(2v), exactly as YC33(8)-(10).

For d>=1 let r_d=floor((d-1)/3) and

\[
 c_\nu(d,n)=\min\{2L_W,\;
 2\epsilon_W(r_d,4)+2L_W n l_\nu(r_d)\},
\]
\[
 \begin{split}
 C_{j,L}(a,n)&=\ell I_{\rm face}j\kappa
                 \sum_{b\ge L}(b+2)^3 c_\nu(b-a,n),\\
 z_{j,L}&=\ell I_{\rm face}j\kappa(L+1)^3
                 b_W(L,4,j\kappa(2L+1)^3),\\
 a_{j,L}(a,n)&=\min\{2,C_{j,L}(a,n)+2z_{j,L}\}.
 \end{split}                                                 \tag{13}
\]

For an observable supported in S(a) with a<L and support size n, the true
and buffered transport, in either orientation, differ by at most
a_(j,L)(a,n) times its norm. To check the geometry, an omitted face has
minimum distance b>=L from S, and at most I_face j kappa(b+2)^3 faces have
that distance. Its distance from S(a) is at least b-a. Retained faces in C
have boundary clearance at least L in B. These are precisely the two
counts used in YC33's Duhamel argument; connectedness of S was not used.

Put a=floor(L/4), n_a=j kappa(a+1)^3, n_(2a)=j kappa(2a+1)^3 and
M_B=j kappa(2L+1)^3. The same two-intermediate composition proof gives

\[
 \boxed{\|\mathcal D_S(\phi)-\mathcal D_S^B(\phi)\|
                         \le\eta_{j,L}\|\phi\|,}             \tag{14}
\]
\[
 \begin{split}
 \eta_{j,L}=\min\{2L_k,\;&4j L_k l_\nu(a)
             +2\epsilon_k(a,n_a)+b_k(2L-a,n_a,M_B)\\
              &+L_k[a_{j,L}(0,j)+a_{j,L}(2a,n_{2a})]\}.
 \end{split}                                                 \tag{15}
\]

Specifically use E_(S(a)) alpha(C_S), then
E_(S(2a)) R E_(S(a)) alpha(C_S). Compare the outer transport on the
second supported intermediate, the physical filters on the first, and
the inner transports on C_S. The replacement errors give (15), including
the factor j in the first transport error. For fixed j, sufficiently high
moments give eta_(j,L)=O_j(L^(11-nu)); nu=3j+16 makes it tend to zero.

The auxiliary classical locality input is
[Nachtergaele--Sims--Young, arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2),
sections3.2,4,5.3 and6, through the explicit YC31/YC33 estimates. Equations
(3)-(15) specify their application to this complete support carrier.

## 5. Full-input assembly and an affine reserve-preserving correction

For d_H(S,T)<=R define the raw block by pairing D_S^B(phi)Omega_B with the
exact T-support subspace inside B, then imposing its inherited physical
gauge projection. Set all other blocks to zero. Vacuum centering may be
done first; it changes none of these nonempty-support blocks. Neither
the exact global U nor the interacting global ground is coefficient input.

Call this operator K_raw. Each retained block differs from K_j^(R) by at
most eta_(j,L), uniformly on its entire source Hilbert space. Every row
and column has at most D_j(R) retained blocks. Thus

\[
 \|K_{\rm raw}-K_j^{(R)}\|\le D_j(R)\eta_{j,L}.
\]

Different source boxes need not give adjoint blocks. Symmetrize first:

\[
 C_{j,R,L}=\tfrac12(K_{\rm raw}+K_{\rm raw}^*),\qquad
 \|C_{j,R,L}-K_j\|\le\epsilon,
 \quad\epsilon=e_j(R)+D_j(R)\eta_{j,L}.                      \tag{16}
\]

The row/column Schur test, rather than one-factor disjoint-ball orthogonality,
is the assembly estimate here. It remains valid for arbitrary entangled
vectors in P_j and keeps the full output within that retained carrier.

Since -epsilon I<=C<= (q+epsilon)I, define the affine correction

\[
 a_\epsilon=\frac q{q+2\epsilon},\qquad
 \widetilde K_j=a_\epsilon(C_{j,R,L}+\epsilon I),\qquad
 \widetilde F_j=I-\widetilde K_j.                             \tag{17}
\]

It uses only local coefficient blocks and a scalar, retaining Hausdorff
range R. No global spectral clipping is required. Moreover

\[
 \boxed{0\le\widetilde K_j\le qI,\quad
 \delta I\le\widetilde F_j\le I,\quad
 \|\widetilde F_j-F_j\|\le\varepsilon_{\rm loc}
       :=\frac{2q\epsilon}{q+2\epsilon}\le2\epsilon.}         \tag{18}
\]

For the error, write K_tilde-K_j=a_epsilon(C-K_j)
 +a_epsilon epsilon I-(1-a_epsilon)K_j. The last two terms together have
norm at most a_epsilon epsilon since 0<=K_j<=qI, and the first has the
same bound. This preserves the original reserve rather than replacing
it by delta minus an approximation error. The affine map is not asserted
to be a unitary change of the physical Hamiltonian or an RG transformation.

If a future numerical block solver has a proven complete-block operator
error at most zeta, D_j(R) zeta can be added to epsilon. Entrywise errors
on finitely many harmonics do not provide that hypothesis.

## 6. Local/count errors, physical sources and the hidden metric

Let S_tilde_(k,j) be the exact Schur reduction of F_tilde_j onto P_k and
W_tilde_(k,j) its minimizing lift, embedded into the full excited carrier.
YC35's variational comparison now applies with the **constructed** value
epsilon_loc from (18):

\[
 \|\widetilde S_{k,j}-S_k\|
 \le E_{k,j}+\varepsilon_{\rm loc}/\delta,
 \qquad E_{k,j}=\tau^2\mathcal B^2 k/(j+1).                  \tag{19}
\]

Both reduced operators retain delta I<=S<=I. If W_k is the full exact
minimizing lift, then

\[
 \|\widetilde W_{k,j}-W_k\|
 \le d_{k,j}:=\tau\mathcal B\sqrt{\frac{k}{\delta(j+1)}}
                       +\frac{\varepsilon_{\rm loc}}{\delta^{3/2}}.
                                                               \tag{20}
\]

The retained-source error on arbitrary b is <=d_(k,j)||b||, and the
physical Gram error is <=2 delta^-1/2 d_(k,j). Both exact minimizing lifts
have norm<=delta^-1/2. Physical vectors are identified by the same exact U,
so these remain their actual Haar-Hilbert errors, with zero vacuum leakage.
An arbitrary inhomogeneous full solution also needs its particular hidden
source response; no approximation of that entire hidden inverse is claimed.

## 7. A finite local elimination, with its residual retained

The Schur inverse in section6 can itself be approximated without a global
matrix inversion. Split P_j into P_k and P_j-P_k and write

\[
 \widetilde F_j=\begin{pmatrix}A&B\\ B^*&H_h\end{pmatrix},
 \quad \delta I\le H_h\le I,\quad
 C_h=I-H_h,\quad 0\le C_h\le qI,
 \quad R_p=\sum_{n=0}^p C_h^n.                               \tag{21}
\]

Take zero-dimensional hidden blocks in the usual empty-sum sense. Set

\[
 S^{\mathrm{loc},p}_{k,j}=A-BR_pB^*,\qquad
 W^{\mathrm{loc},p}_{k,j}=\binom{I}{-R_pB^*},                 \tag{22}
\]

and embed the lift by zero outside P_j. The geometric identity gives

\[
 H_h^{-1}-R_p=C_h^{p+1}H_h^{-1}\ge0,
\]
\[
 \boxed{\delta I\le\widetilde S_{k,j}
          \le S^{\mathrm{loc},p}_{k,j}\le A\le I,
 \quad 0\le S^{\mathrm{loc},p}_{k,j}-\widetilde S_{k,j}
              \le\frac{q^{p+3}}\delta I.}                   \tag{23}
\]

Here ||B||<=q because B is an off-diagonal block of -K_tilde. Combining
(19),(23) proves (1). The polynomial lift has additional norm error
l_p=q^(p+2)/delta. Its hidden stationarity residual is C_h^(p+1)B*, bounded
by q^(p+2). In particular W_p* F_tilde_j W_p is not silently equated with
S_p: the minimizing-lift identity applies to the exact lift only.
The full source/lift error relative to W_k is at most d_(k,j)+l_p, and
the Gram error is at most

\[
 (2\delta^{-1/2}+l_p)(d_{k,j}+l_p).                          \tag{24}
\]

Products in (22) follow at most p+2 coefficient edges. The Hausdorff triangle
inequality gives support range (p+2)R. A retained row based at S therefore
uses parent coefficients only within S((p+2)R+2L), including the boxes
needed by symmetrization. Its total factor count is at most
|S| kappa((p+2)R+2L+1)^3, regardless of the parent volume or separation
between factors of S. This is a union of bounded neighbourhoods, not a
claim that a widely separated source fits inside one bounded-diameter ball.
Local harmonics still have not been discretized.

For any target return error epsilon_target>0, first choose j>=k to make
E_(k,j) at most a third of it, then choose R and L to make
epsilon_loc/delta at most a third, then choose p for the last term.
Use nu=3j+16 after choosing j. Each choice exists uniformly in the admitted
volumes; no useful numerical complexity estimate follows from these coarse
constants. If j exceeds the number of factors, the actual count error is zero.

## 8. Typed zero specialization and the next genuine gate

At ell=0 both transports are identity and every physical/filter comparison
preserves the exact source support S. Consequently K_(T,S)=0 for T!=S,
the union-box coefficients are exact, and the actual bounds e_j,eta_(j,L)
may be set to zero. Use those exact zero values rather than the loose
positive-coupling envelopes with v in their denominators. Then epsilon=0,
the affine correction is identity, and on support S

\[
 K_S=\frac\tau{H_S+\tau},\qquad
 F_S=\frac{H_S}{H_S+\tau},\qquad H_S=\sum_{f\in S}h_f.       \tag{25}
\]

Energies add before filtering; products of separately normalized readings
would give different values. F is count-block diagonal, so B=0 in (21).
The hidden elimination and its polynomial are exact for every p, W=I_k,
and the source/pairing specialize to the original reference carrier. The
first source jet and ordered transport history are retained data, not erased
by this zero value. This proves the declared operator, cut, source, metric
and reduction compatibility along the full admitted family.

Under A_phys=c Ahat and tau_phys=c tau, K,F,delta,q and the coefficient
errors are dimensionless and unchanged; YC35's B scales as1/c. Thus (1)
is invariant under that unit conversion. It does not determine c or a
continuum coupling trajectory.

Written proofs establish the all-volume estimates above. The focused exact
script checks the nearest-pair/Hausdorff distinction, the affine reserve
repair and its sharp error, finite hidden elimination with its nonzero
residual, and the exact zero specialization. It does not certify numerical
Yang--Mills box spectra, harmonic tails or practical parameter choices.

**Next:** an actual spatial blocking map with complete coarse gauge factors,
followed by control of its generated interactions, sources and metrics under
iteration. The global count cut is not product-closed: two j-factor sources
can produce a 2j-factor source. Hausdorff locality at fixed j therefore does
not by itself supply a tensor-product coarse Hamiltonian or a contracting
spatial renormalization law. Continuum energy calibration remains separate.
