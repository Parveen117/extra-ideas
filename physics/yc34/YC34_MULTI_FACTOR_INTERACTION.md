# YC34 — simultaneous sources and a vacuum-annihilating interaction

10 October 2026. Research owner: Monty Dabas. Continues YC33 at `6407977`.
The complete YC27 factors, internal profiles, physical sectors and admitted
mixed-coupling windows remain. Frozen predecessor packets are unchanged.

**Result.** Construct the full graded excitation carrier first, including
entangled charge-matched states. In YC31's moving vacuum frame, the complete
Hamiltonian admits an exact expansion into self-adjoint local terms, each
annihilating the product reference vacuum. Truncating their spatial radius
at R has the bound

\[
 \boxed{\|[\widehat A_s-\widehat A_s^{(R)}]P_{\le k}\|
          \le s\,\mathcal E_k(R),\qquad
          \mathcal E_k(R)\longrightarrow0.}                 \tag{1}
\]

Here k counts excited **reference factors**, not physical particles. It is
fixed independently of volume; source supports may be arbitrarily far
apart, and the output is not projected back into that k-factor space.
Every factor retains all its harmonics. An explicit bound appears in (18).
For k=2, the flow exponent32 gives a conservative O(R^(-33/2)) tail.
No practical radius or numerical prefactor is certified.

This is a local interaction representation on the original factor graph,
not a spatial blocking step or closure into a new Wilson Hamiltonian.
The exact multi-source Fredholm reduction keeps its old reserve, but (1)
does not yet approximate that reduction's entire hidden inverse. The
existing lattice gap is an input; the physical continuum limit stays open.

## 1. The full source carrier and its products

On each complete factor put

\[
 p_f=|\Omega_f\rangle\langle\Omega_f|,\quad q_f=I-p_f,
 \quad \mathcal H_f^+=q_f\mathcal H_f,\quad N=\sum_f q_f.
                                                               \tag{2}
\]

Use the full factor space, not only its gauge singlets. For each subset I
of factors, let H_I^+ be the tensor product of these excited spaces on I,
with the vacuum on the other factors. The complete tensor carrier is the
orthogonal direct sum over I, including the empty support. Taking global
gauge invariants in each summand gives its physical carrier. Gauss matching
is imposed jointly, so an I-summand can include entangled charged factors.
The projections commute with gauge transformations because every Omega_f
is invariant. This retains all physical centre sectors as well.

For phi in the I-summand, define the bounded creator

\[
 C_I(\phi)=|\phi\rangle\langle\Omega_I|\otimes I_{I^c},
 \qquad \|C_I(\phi)\|=\|\phi\|.                              \tag{3}
\]

For disjoint I,J their product is C_(I union J)(phi tensor psi), with the
specified tensor ordering. If I and J overlap, their product is zero:
the vacuum bra on a shared factor annihilates its excited ket. Creators
therefore commute, with this overlap rule. Their adjoints and the full
local operator algebra remain available. At nonzero coupling all these
operations are conjugated by the same exact U_s; products, adjoints and
the Haar pairing are preserved. This is a graded source construction,
not an assertion that physical excitations are free bosons or fermions.

Let P_(<=k)=1_[0,k](N), initially on the full tensor space. It retains
vacuum and all states with at most k excited factors. Its excited physical
part will be denoted P_k^phys. It is an optional later cut of the full
carrier, not the primitive carrier. It is not invariant under the full
Hamiltonian, and products can exit it. The full source embedding U_s is
unitary; its restriction to any of these orthogonal source spaces is an
isometry with identity Gram pairing.

The lambda-space data are the admitted operator family, domains, product,
adjoint, source cuts, gauge action, pairing and ordered transport history.
The path s in H(s)=H0-s ell sum_p xi_p Wp remains a scoped coupling adapter,
not primitive lambda-space, physical time or an RG scale.

## 2. Differentiate in the moving frame before localizing

Use YC33's identity, including its sign convention U_s'=iD_s U_s:

\[
 B_s=H'(s)+i[H(s),D_s]=\int w_\gamma(t)\tau_t^s(H')dt,
 \qquad [B_s,\Pi_s]=0.                                       \tag{4}
\]

Here gamma=1/4, w_gamma is real, even and normalized, with Fourier support
inside the full tensor gap. Keep L_w=integral |w_gamma| in bounds (the
nonnegative filter of YC31 has L_w=1). Put

\[
 \widehat A_s=U_s^*(H(s)-E_s)U_s,
 \quad G_s=U_s^*B_sU_s-E_s'I,\quad
 \widehat A_s=H0+\int_0^sG_u\,du,
 \qquad G_s\Omega=0.                                        \tag{5}
\]

At finite volume these are identities on the common electric domain.
Equivalently, first differentiate on the form core as in YC33; the
right-hand derivative in (5) is bounded, and integration gives the bounded
difference from H0 and hence the same self-adjoint domain. In particular
this construction does not differentiate a guessed local ground state.

For each mixed face let j_p=-ell xi_p Wp and define

\[
 a_p(s)=U_s^*\left[\int w_\gamma(t)\tau_t^s(j_p)dt\right]U_s,
 \quad c_p(s)=\langle\Omega,a_p(s)\Omega\rangle.
                                                               \tag{6}
\]

Then sum_p c_p=E_s' and G_s=sum_p(a_p-c_p I). Each a_p is self-adjoint,
quasi-local and of norm <=ell L_w. The *sum* kills Omega; its individual
centered terms need not do so. Subtracting just their means is insufficient
for a volume-independent simultaneous-source bound.

## 3. Local shells, with a volume-independent size and incidence bound

Let Xp(r) be the graph r-neighbourhood of the four-factor mixed face Xp.
Write E_Z for partial expectation outside Z in the product reference.
It preserves the vacuum expectation, adjoints and inherited gauge symmetry.
Use YC31's constants kappa=6750 and I=48 or88. Uniformly in the torus,

\[
 |Xp(r)|\le n_r:=4\kappa(r+1)^3,\qquad
 \#\{p:f\in Xp(r)\}\le I\kappa(r+1)^3.                       \tag{7}
\]

The second count follows by selecting a face factor within distance r
of f and using the at-most-I incident faces there.

Choose a flow exponent nu>4 in YC31 and put
l_nu(r)=min{2,C_nu(r+1)^(4-nu)}, with C_nu its flow-localization constant.
Use epsilon_w(r,b) from YC31(5), for a b-factor source. For r>=2 let
u=floor(r/2) and set

\[
 t_\nu(r)=\min\{2L_w,\;
   2[\epsilon_w(u,4)+4\kappa(u+1)^3L_w l_\nu(u)]\}.
                                                               \tag{8}
\]

Set t_nu(0)=t_nu(1)=2L_w. Localize the physical filter in Xp(u), then
its inverse transport in Xp(2u); replacing the resulting supported
approximant by E_(Xp(r)) costs at most another factor two. Thus
||a_p-E_(Xp(r))a_p||<=ell t_nu(r). All filter moments are available, so
t_nu(r)=O(r^(7-nu)) for nu>7, taking a sufficiently high filter-tail bound.

Define the centered local shells

\[
 \begin{split}
 a_{p,0}&=E_{Xp}a_p-c_pI,\\
 a_{p,r}&=(E_{Xp(r)}-E_{Xp(r-1)})a_p\quad(r\ge1),\\
 b_0&=2L_w,\qquad b_r=t_\nu(r)+t_\nu(r-1)\quad(r\ge1).
 \end{split}                                                  \tag{9}
\]

They have support Xp(r), zero vacuum mean and norm <=ell b_r. At finite
volume the telescoping sum terminates once a ball is the whole connected
component. Therefore G_s=sum_(p,r) a_(p,r) exactly. The infinite positive
series used for estimates only extends this finite-volume sum.

## 4. A symmetric source subtraction that respects spectators

For a bounded operator a on a finite factor set X, decompose its source
vector into its exact excitation supports,

\[
 a\Omega_X=\sum_{J\subset X}v_J\otimes\Omega_{X\setminus J},
 \qquad v_J\in\bigotimes_{f\in J}\mathcal H_f^+.
\]

Define the creation lift, with identities on the spectator factors,

\[
 \Gamma_X(a\Omega_X)=v_\varnothing I+
           \sum_{\varnothing\ne J\subset X}C_J(v_J).          \tag{10}
\]

This is the source-support lift used in YC19, now on the full YC27 carrier.
No YC19 numerical constants or its l1-to-Hilbert-norm identification are
imported. The new estimate below is a Hilbert-space estimate.

For self-adjoint a with zero vacuum mean, put

\[
 \mathfrak n_X(a)=a-\Gamma_X(a\Omega_X)
                    -\Gamma_X(a\Omega_X)^* .                 \tag{11}
\]

This operator is self-adjoint, supported on X, and kills Omega_X from
both sides. The key cancellation is linear in the *complete* source.
Extending any local source to the global vacuum and then applying Gamma
is exactly its local Gamma with identities outside. Since G_s Omega=0,

\[
 \sum_{p,r}\Gamma_{Xp(r)}(a_{p,r}\Omega_{Xp(r)})=0,
 \qquad G_s=\sum_{p,r}\mathfrak n_{Xp(r)}(a_{p,r}).             \tag{12}
\]

Thus no physical interaction has been thrown away. Taking the adjoint
subtraction makes this a self-adjoint extension of the source-removal
mechanism, without a nonunitary metric change. Gauge invariance follows
because every excitation-support projection commutes with the gauge
action and each source component remains invariant under the joint action
on its support. The same argument preserves the inherited centre action.

A different tempting prescription, (I-p_X)a(I-p_X), does **not** have the
cancellation property (12) for different supports X. For example, on two
two-level factors take a_1=sigma_x on factor1 and a_12=-sigma_x tensor I.
Their sum is zero. Separate local-vacuum sandwiching leaves a nonzero
flip conditioned on factor2 being excited. The lift (10) instead cancels
exactly: it retains the spectator identity rather than a spectator-vacuum
projector. The direct check includes this counterexample.

## 5. The creation lift costs excitation count, not local dimension

Let n=|X| and

\[
 D(n,k)=\sum_{j=0}^{\min(k,n)}{n\choose j}.
\]

For any local vector v, with its creation lift defined as in (10),

\[
 \boxed{\|\Gamma_X(v)P_{\le k}\|\le\sqrt{D(n,k)}\,\|v\|.}    \tag{13}
\]

Proof: write an input as sum_J psi_J over its local excitation supports,
|J|<=k. The output with support K is
sum_(J subset K, |J|<=k) v_(K minus J) tensor psi_J. There are at most
D(n,k) summands. Cauchy--Schwarz, followed by summing K, gives
D(n,k) sum_(J,I disjoint J)||v_I||^2||psi_J||^2, at most
D(n,k)||v||^2||psi||^2. Tensoring the bound with the outside Hilbert space
preserves it. No basis or harmonic count occurs.

Gamma* only lowers excitation support, so its norm on P_(<=k) obeys the
same bound by taking the adjoint of the compression to that subspace.
Consequently, for the centered a in (11),

\[
 \|\mathfrak n_X(a)P_{\le k}\|
       \le[1+2\sqrt{D(n,k)}]\|a\|.                            \tag{14}
\]

There is no claim of this polynomial-in-n bound uniformly in k. The full
creation lift can have a much larger norm. At fixed k, D(n,k) grows at
most polynomially in n; this is sufficient to sum the spatial tails.

## 6. Sum overlapping interactions without an extensive factor

Here is the required collective estimate. Suppose self-adjoint local
operators z_X satisfy z_X Omega_X=0, |X|<=n, and
||z_X P_(<=k)||<=w_X. Put J_k=sup_f sum_(X containing f) w_X. Then

\[
 \left\|\sum_X z_X P_{\le k}\right\|
                \le J_k\sqrt{k(k+n)}.                       \tag{15}
\]

Proof: Q_X=I-p_X satisfies z_X=Q_X z_X Q_X and Q_X<=sum_(f in X)q_f.
For psi in P_(<=k), every output z_X psi lies in P_(<=k+n). It is enough
to test against chi in that latter subspace. Weighted Cauchy--Schwarz gives

\[
 \begin{split}
 |\langle\chi,\sum_Xz_X\psi\rangle|
 &\le\sum_Xw_X\|Q_X\chi\|\|Q_X\psi\|\\
 &\le J_k
   \langle\chi,N\chi\rangle^{1/2}
   \langle\psi,N\psi\rangle^{1/2}
 \le J_k\sqrt{k(k+n)}\|\chi\|\|\psi\|.
 \end{split}                                                  \tag{16}
\]

The restricted local norm bound applies because Q_X commutes with
P_(<=k). This argument does not assume disjoint sources, conserved factor
number, factor singlets or finite-dimensional local Hilbert spaces.

Set z_(p,r)=mathfrak n_(Xp(r))(a_(p,r)) and construct

\[
 \widehat A_s^{(R)}=H0+
             \int_0^s\sum_p\sum_{r=0}^Rz_{p,r}(u)\,du.       \tag{17}
\]

It is self-adjoint on Dom(H0) at every finite volume, is gauge invariant,
and has Omega as an exact eigenvector of energy zero. It has interaction
diameter at most2R+1 on the factor graph; its terms may be many-factor
operators. Each local term annihilates its local vacuum, but it need not
be positive. In particular Omega is not asserted to be the global ground
of this approximate operator on arbitrarily many excited factors.

These are supported local terms whose coefficients are still defined
through the exact global flow. YC33's buffered coefficients were for a
different, one-factor inverse-source construction. Replacing the present
coefficients by independently computed box data needs its own comparison
and a repair of the exact source cancellation (12); it is not included
in (18). A finite-volume identity with fixed-k uniform bounds also does
not by itself construct a thermodynamic limiting Hamiltonian.

Equations(7), (9), (14)-(16) prove (1), with the explicit budget

\[
 \boxed{\mathcal E_k(R)=\ell I\kappa
  \sum_{r>R}(r+1)^3 b_r
       [1+2\sqrt{D(n_r,k)}]\sqrt{k(k+n_r)},\qquad k\ge1.}      \tag{18}
\]

For k=0 the error is exactly zero. Bounds are uniform along 0<=s<=1
and independent of total volume. Since n_r=O(r^3) and b_r=O(r^(7-nu)),
the summand is O(r^((23+3k)/2-nu)). Hence, for
nu>(25+3k)/2,

\[
 \mathcal E_k(R)=O(R^{(25+3k)/2-\nu}).                        \tag{19}
\]

Choose nu=32 for example: k=1 gives O(R^-18), k=2 gives O(R^(-33/2)),
and k=3 gives O(R^-15). Any fixed k can be handled by choosing a higher
flow exponent and the required filter moments. Constants are deliberately
large and unoptimized; no useful numerical R is selected by these powers.

The locality input is YC31's spectral-flow and filtered-dynamics estimate.
The general unbounded-on-site setting is treated in sections3.2,5,6 of
[Nachtergaele--Sims--Young, arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2).
The source subtraction and bounds (13)-(19) are proved here rather than
attributed to an unverified stronger theorem in that paper.

## 7. What this controls, and the remaining inverse gate

The estimate is on inputs with at most k excited factors and keeps their
*entire* output. It therefore controls both the compressed interaction
and the leakage out of that source cut, to absolute error E_k(R).
It does not say that the leakage itself is small.

On P_k^phys define the ordinary Hamiltonian compressions A_k and A_k^(R)
of the two operators at s=1, using the restricted H0 domain. This domain
is well defined because H0 commutes with the factor-number projections.
Since the exact transported Hamiltonian has physical gap m,

\[
 A_k\ge mI,\qquad
 \|A_k-A_k^{(R)}\|\le\mathcal E_k(R),\qquad
 A_k^{(R)}\ge[m-\mathcal E_k(R)]I.                            \tag{20}
\]

The norm in (20) is of their bounded difference. If e=E_k(R)<m, their
bounded normalized returns obey the resolvent-identity estimate

\[
 \|\tau(A_k+\tau)^{-1}-\tau(A_k^{(R)}+\tau)^{-1}\|
       \le {\tau e\over(m+\tau)(m-e+\tau)}.                  \tag{21}
\]

These are compressions, **not exact elimination of the remaining factors**.
For that elimination retain the actual physical excitation operator
A=(H-E)|Q and the source isometry V_k=U|_(P_k^phys). YC30 gives instead

\[
 T_k=V_k^*A^{-1}V_k,\qquad S_k=(I+\tau T_k)^{-1},
 \qquad {m\over m+\tau}I\le S_k\le I.                       \tag{22}
\]

All other excitation supports remain in its hidden space Q_h=Q-V_k V_k*.
The exact lift is W_k=V_k+tau Q_h A^-1 V_k S_k. It retains the complete
input source through W_k* b and has physical pairing W_k*W_k<=
(m+tau)/m. The particular hidden solution is still needed for a general
inhomogeneous source. The isometric initial source Gram and this lifted
metric are different objects.

Equation(1) alone does not control the error after replacing A inside
(22): A^-1 V_k need not stay within any fixed factor-number cut. A next
step needs a quantitative factor-number/energy bound on these inverse
responses, or an alternative weighted inverse estimate. Inferring the
full Schur error from (21) would discard that hidden response.

There is a simple reason to retain the k dependence. The local error
epsilon sum_f q_f annihilates the vacuum and has norm epsilon on a
single excited factor, but norm epsilon min(k,|Lambda|) on P_(<=k).
Its unrestricted norm grows with volume. Local centering alone cannot
give a volume-independent full-operator error from such per-site errors.
This example does not rule out a sharper energy-weighted resolvent bound.

## 8. Typed zero specialization and verification

Evaluate the mixed-strength adapter at ell=0 after constructing the full
operator family. Then H=H0, E=0, U=I, all a_p and z_(p,r) vanish and
both (5) and (17) recover H0 exactly. Sources, source products, the full
graded carrier and Haar pairing remain; correlated internal factor
Hamiltonians are not set to zero. Set E_k(R)=0 there directly.

On an exact nonempty support J the zero operator is the **sum**
H_J=sum_(f in J) h_f, restricted to the globally matched physical space.
It preserves that support, so the zero-reading of (22) is

\[
 T_{k,0}=\bigoplus_{1\le|J|\le k}H_J^{-1},\qquad
 S_{k,0}=\bigoplus_{1\le|J|\le k}{H_J\over H_J+\tau}.         \tag{23}
\]

The direct sums here are over the physical support subspaces; empty ones
are omitted. Spectator energies add before taking the inverse. Even at
zero mixed coupling, this is not a product of the individual normalized
returns: for energies2,3 and tau4, the pair return is5/9, whereas the
product of the individual returns is1/7. Tensor source composition and
functional calculus are different typed operations.

Thus pi0(Construct_lambda)=Construct0(pi0(inputs)) holds for the carrier,
products, transport, interaction representation, source metric and exact
normalized reduction just constructed. A derivative source at zero is
separate data and is not erased by evaluating the interaction itself at
zero. No claim identifies this scoped adapter with a completed universal
lambda-space or with the physical continuum limit.

Verification: written full-family, domain, cancellation, collective and
zero-specialization proofs above. The small exact script checks the
spectator-sandwich counterexample, symmetric source cancellation with
overlapping supports, the creation-lift Hilbert bound on a finite model,
the excitation-count dependence, and additive rather than multiplicative
zero returns. Those examples check new algebra; they are not finite
harmonic evidence for Yang--Mills. No predecessor suite, new numerical
gap certificate or frozen packet was regenerated.

**Next obligation:** control factor-number growth in the exact hidden
inverse response and combine it with the spatial/box errors, while keeping
the full metric. New-scale interaction closure and physical energy/spacing
calibration remain subsequent obligations; YC27's gap windows are unchanged.
