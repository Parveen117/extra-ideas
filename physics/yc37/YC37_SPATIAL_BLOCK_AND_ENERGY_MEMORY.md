# YC37 — spatial gauge blocks, returned metric and energy memory

10 October 2026. Research owner: Monty Dabas. Continues YC36 at `aa3177f`.
**Theory-first development:** written constructions and derivations below;
no new executable certificate, test suite or numerical validation packet.
This follows the owner's request to prioritize theory and postpone additional
certification work. Frozen predecessor results are unchanged.

The new construction has three parts. First assemble complete spatial gauge
factors, retaining boundary charges. Then make a local count cut and restore
its returned metric before identifying a retained coarse Hamiltonian. Finally
retain boundary-insertion histories when joining these factors. A single
static effective Hamiltonian does not contain all the data for an exact join.

For a block of b original factors, the particular count cut below gives

\[
 D_B\ge c b I,\qquad I\le M_B\le(1+\beta^2)I,
 \qquad
 \|F_B(z)-(K_B-zM_B)\|
       \le\frac{\beta^2|z|^2}{cb-z_*},\quad |z|\le z_*<cb.   \tag{1}
\]

Thus at a fixed energy window the nonlinear remainder after retaining the
metric is O(1/b). This is a **single-block energy estimate**, with a retained
count cap growing as b/8. It is not a contraction of the interblock interaction
or a continuum limit. The full join retains explicit boundary histories.

## 1. Complete coarse gauge factors before the local cut

Partition YC27's complete factors into spatial patches B. Do not split an
original factor or discard one of its charged Hilbert spaces. The corresponding
microscopic edge sets partition the links; blocks can still share vertices.
Let V_int(B) contain vertices whose incident microscopic links all lie in B.
First impose Gauss invariance only at these internal vertices:

\[
 \mathcal K_B=
 \left(\bigotimes_{f\in B}\mathcal H_f\right)^{
                       \prod_{v\in V_{\rm int}(B)}SU(2)_v}.
                                                               \tag{2}
\]

Every remaining boundary vertex acts on the incident blocks jointly. With
G_boundary this product of boundary gauge groups, the exact physical carrier is

\[
 \mathcal H_G\simeq
           \left(\bigotimes_B\mathcal K_B\right)^{G_{\partial}}.
                                                               \tag{3}
\]

This follows by commuting Haar averages at distinct vertices and regrouping
the tensor factors. It preserves the original Haar pairing. In particular,
K_B is not restricted to boundary singlets: its decomposition under the
boundary group retains every boundary representation and its multiplicity
space. Physical gluing matches these representations at shared vertices.
This is a block factor with boundary gauge data, not an asserted single
L2(SU(2)) coarse link.

The represented operators, their domains, source cuts, gauge actions,
adjoints, pairing and ordered history are the underlying lambda-space data.
The mixed coupling ell is an adapter. Block size b, reference count and the
spectral variable z are different quantities; none is identified with a
primitive dimension, clock or a continuum running coupling.

## 2. Local block Hamiltonian and a concrete admissible count cut

Use the fixed YC27 profiles and 0<=xi_p<=1. On a block with b factors put

\[
 H_{0,B}=\sum_{f\in B}h_f,\quad
 V_B=-\ell\sum_{X_p\subset B}\xi_pW_p,\quad
 \bar h_B=H_{0,B}+V_B,\quad h_B=\bar h_B-E_B,
                                                               \tag{4}
\]

where E_B is the actual lowest block energy. Carry this scalar separately
when blocks are joined. All potentials are bounded; the inherited electric
domain is unchanged. Internal Gauss restriction commutes with these operators.

Two inherited facts give elementary bounds without a new block solve.
At mixed zero the full tensor reference has H0>=N/2, using YC31's full-space
floor, not only its physical-sector floor. Also the product reference has
zero expectation of every mixed plaquette: its four-factor source is excited
in all four factors, as in YC27 section4. Therefore E_B<=0.

Every mixed face meets four distinct original factors. If I_face is48 or88,
its incidence count implies

\[
 n_{\rm int}(B)\le I_{\rm face}b/4,\qquad
 \|V_B\|\le a b,\quad a=\ell I_{\rm face}/4.                 \tag{5}
\]

More generally a retained fraction sigma in(2a,1) gives a positive coefficient
c=sigma/2-a with r_B=floor(sigma b). Choose sigma=1/8 for a simple common
budget in both profiles; it is a retention choice, not a universal constant
derived from the framework. For b>=8, on K_B define

\[
 N_B=\sum_{f\in B}q_f,\qquad
 P_B=\mathbf1_{[0,r_B]}(N_B),\qquad Q_B=I-P_B.                \tag{6}
\]

This local cut includes the reference vacuum and all harmonics and matched
charges of its admitted supports. It is taken only after (2); omitted charge
and multiplicity channels remain in Q_B. Unlike a global count cut, the
tensor product of these retained factors allows arbitrarily many blocks to
be excited. Its dimension has not been fixed independently of b.

P_B commutes with H0,B and with the boundary gauge action. Its domains split
the reference electric domain. Hence the blocks in

\[
 h_B=\begin{pmatrix}A_B&T_B\\T_B^*&D_B\end{pmatrix}_{P_B,Q_B}
                                                               \tag{7}
\]

have self-adjoint diagonal parts on their inherited domains and a bounded
off-diagonal T_B=P_B V_B Q_B. From E_B<=0, (5), and r_B+1>b/8,

\[
 \boxed{D_B\ge cbI,\quad \|T_B\|\le ab,
               \qquad c=1/16-a>0.}                          \tag{8}
\]

At the two existing endpoint couplings one may use these rational values:

| Profile | a | c | beta=a/c |
|---|---:|---:|---:|
| 28-link, ell<=1/3200 | 3/800 | 47/800 | 3/47 |
| 64-link, ell<=1/5700 | 11/2850 | 1337/22800 | 88/1337 |

They are direct symbolic substitutions in (5),(8), not newly run numerical
certificates. For smaller ell use its own a and c, or these endpoint bounds.

There is also an inherited block gap. Turn every interblock xi_p to zero in
the original admitted family. The resulting full tensor Hamiltonian is the
sum of the independent block Hamiltonians. YC31's full-space gap remains
at least1/2 for that admitted inhomogeneous choice. Exciting just one block
therefore shows that each h_B has gap>=1/2 above its unique block ground.
This argument uses the already proved parent family; it is not a general
principle that arbitrary open boxes inherit a bulk gap.

## 3. Returned metric and the retained coarse Hamiltonian

Suppress B temporarily. Since D>=cb>0, define on the retained factor

\[
 K=A-TD^{-1}T^*,\quad Y=D^{-1}T^*,\quad
 W=\binom{I}{-Y},\quad M=W^*W=I+Y^*Y,
 \qquad J=WM^{-1/2}.                                        \tag{9}
\]

The exact quadratic-form decomposition is

\[
 \langle(u,v),h(u,v)\rangle
   =\langle u,Ku\rangle+\langle v+Yu,D(v+Yu)\rangle.         \tag{10}
\]

In particular K>=0 and W*hW=K. The returned metric and isometry obey

\[
 I\le M\le(1+\beta^2)I,\qquad J^*J=I,\qquad
 h^{\rm ret}=M^{-1/2}KM^{-1/2}=J^*hJ.                       \tag{11}
\]

The last equality is an equality of closed quadratic forms. Define its domain by
M^(-1/2)x in Dom(K^(1/2)); no unsupported assertion that an arbitrary bounded
metric square root preserves the electric domain is needed. All operators
in (9) intertwine the inherited boundary gauge action. Thus J is an actual
gauge-equivariant physical isometry, not a similarity with a discarded metric.

Let psi be the block ground. Its hidden equation gives Q psi=-Y P psi, so
psi=W P psi lies in Ran J. Its retained normalized ground is
omega_ret=M^(1/2)P psi. Since J preserves both the norm and orthogonality to
this ground, the inherited block floor transfers directly:

\[
 \boxed{h^{\rm ret}\ge\tfrac12
        (I-|\omega_{\rm ret}\rangle\langle\omega_{\rm ret}|).} \tag{12}
\]

No fresh onsite spectral approximation is hidden in this statement. The
new factor retains the old onsite floor, but the interactions between these
factors still have to be handled.

There is a useful low-energy residual identity. Let Pi_J=JJ*, let i_P embed
the original retained channel, and set L=(I-Pi_J)i_P M^(1/2). On the retained
operator domain, hJ=i_P M^(1/2)h_ret, while

\[
 L^*L=M-I,\qquad
 \boxed{\|(I-\Pi_J)hJ\,\mathbf1_{[0,E]}(h^{\rm ret})\|
                      \le\beta E.}                          \tag{13}
\]

This is an energy-restricted source estimate; it is not a full operator-norm
bound on an unbounded residual. For exact Feshbach identities we retain the
original domain-splitting P,Q,D of (7). No unproved domain property of the
rotated complementary projection I-Pi_J is substituted for that split.

## 4. The energy-memory term decreases after its metric is retained

The exact retained spectral pencil is

\[
 F_B(z)=A_B-zI-T_B(D_B-zI)^{-1}T_B^*,\qquad |z|<cb.           \tag{14}
\]

Its zero reading is K_B, and -F_B'(0)=M_B. Resolvent expansion gives the
exact identity

\[
 F_B(z)=K_B-zM_B-z^2T_BD_B^{-2}(D_B-zI)^{-1}T_B^*.           \tag{15}
\]

Since ||T_B D_B^-2 T_B*||<=beta^2, (15) proves (1). For a fixed z_* and
b>=2z_*/c the remainder is at most 2 beta^2 |z|^2/(cb). For real z<cb the
last term is positive before its minus sign. Whitening by M_B^-1/2 gives
a pencil h_B^ret-zI with the same or smaller remainder norm.

More explicitly, the spectral response has

\[
 0\le-F_B''(0)=2T_BD_B^{-3}T_B^*
                     \le\frac{2\beta^2}{cb}I.               \tag{15a}
\]

This is curvature of the energy-dependent operator pencil, not an automatic
identification with spatial Riemann curvature. If z in the stated window is
an actual block eigenvalue, its eigenvector has a nonzero P_B component and
F_B(z)P_B psi=0. The whitened remainder therefore implies
dist(z,spectrum(h_B^ret))<=beta^2 |z|^2/(cb-z_*). This direction of spectral
approximation does not by itself assert eigenvalue multiplicities or a
converse count theorem.

Keeping only K_B and replacing the metric by I would instead drop the
first-order term -z(M_B-I). Keeping the returned metric exposes a genuinely
quadratic energy remainder. The derivative metric itself is the existing
YC20 Schur/Feshbach mechanism. What is supplied here is the concrete spatial
count cut, the b-dependent hidden floor, and the resulting O(1/b) bound.

This does not mean that every error in a spatial RG step decreases as1/b.
The energy window is fixed, the count cap grows with b, and the boundary
interaction has not yet been included in (15). In particular this is not
a limit in which the lattice spacing shrinks.

## 5. An actual tensor-product coarse carrier and its compressed interaction

For a spatial partition define

\[
 \mathcal H^{\rm ret}_G=
       \left(\bigotimes_B P_B\mathcal K_B\right)^{G_\partial},
 \qquad \mathcal J=\bigotimes_B J_B.
                                                               \tag{16}
\]

Gauge equivariance makes this an isometry into the complete physical fine
carrier. Its source Gram is exactly identity; no product of local condition
numbers replaces the actual pairing. The fine Hamiltonian, relative to
the carried scalar sum_B E_B, is
H_join=sum_B h_B+V_boundary. Its form compression is

\[
 \mathcal J^*H_{\rm join}\mathcal J
       =\sum_B h_B^{\rm ret}
           +\sum_{p\ {m crossing}}\widetilde V_p,
 \quad \widetilde V_p=\mathcal J^*(-\ell\xi_pW_p)\mathcal J,
 \qquad \|\widetilde V_p\|\le\ell\xi_p.                     \tag{17}
\]

Each transformed face still acts only on the coarse factors visited by that
face and is gauge invariant under their joint action. It need not be a Wilson
multiplication operator on a single SU(2) coarse link. The correct admitted
coarse interaction class is therefore broader than the original plaquette
ansatz. It includes bounded gauge-invariant operators with the specified
boundary representation indices and the local unbounded onsite generators.

Equation(17) is an actual spatial compression with complete source pairing.
It is not the full effective Hamiltonian of the joined system: internal
residuals (13), boundary excursions and their cross terms can leave Ran J.
They must return through the complementary channel. The following construction
keeps the data for that return instead of relabelling (17) as an exact closure.

## 6. Even independent joins require the whole spectral response

For a lower-bounded self-adjoint block h and an isometry j into its full
carrier, retain the operator-valued spectral measure

\[
 \mu_j(dE)=j^*\mathsf E_h(dE)j,\qquad
 C_j(t)=j^*e^{-th}j=\int e^{-tE}\mu_j(dE),\quad t\ge0.
                                                               \tag{18}
\]

For z below the spectrum, its resolvent is
G_j(z)=integral (E-z)^-1 mu_j(dE). If j is the original P inclusion, this
is F_B(z)^-1. If j=J_B from (9), it is the response of the metric-isometric
source instead; these two source conventions must not be silently equated.

For independent blocks h_12=h_1 tensor I+I tensor h_2 and j_12=j_1 tensor j_2,

\[
 C_{12}(t)=C_1(t)\otimes C_2(t),\qquad
 G_{12}(z)=\iint\frac{\mu_1(dE_1)\otimes\mu_2(dE_2)}
                              {E_1+E_2-z}.                  \tag{19}
\]

The measure joins by energy convolution; this rule is associative by the
joint spectral theorem. Scalar energy shifts must be included in the sums.
A value of G or F at a single energy is insufficient to execute (19).

For example take h=[[2,1],[1,2]] and j=e1. Its spectral weights are1/2 at
energies1 and3, so G(0)=2/3 and F(0)=3/2. For two independent copies the
energies are2,4,6 with weights1/4,1/2,1/4. Consequently

\[
 G_{12}(0)=7/24,\qquad F_{12}(0)=24/7\ne3=F_1(0)+F_2(0).
                                                               \tag{20}
\]

This is a displayed algebraic illustration, not a numerical YM certificate.
The discrepancy is already present without an interblock interaction.

Projected propagation also keeps its missing channel explicitly:

\[
 C_j(t+s)-C_j(t)C_j(s)
       =j^*e^{-th}(I-jj^*)e^{-sh}j.                           \tag{21}
\]

At t=s this is a positive Gram operator. It vanishes for every t precisely
when the retained subspace reduces h. Thus projected propagation is generally
not a semigroup on the retained factor. Replacing it by an autonomous
generator loses this return unless an approximation estimate is supplied.

## 7. Boundary histories make the interacting join exact

Each crossing Wilson interaction has a finite expansion into bounded block
operators with contracted SU(2) matrix indices:

\[
 V_\partial=\sum_\nu g_\nu\bigotimes_B O_{B,\nu}.            \tag{22}
\]

Only blocks touched by the face have nonidentity entries. Internal vertex
indices are contracted within each block; external indices are matched in
the global gauge invariant sum. Individual boundary operators can carry
boundary charge. Keeping only their singlet components would change V.

For every ordered boundary word retain

\[
 \mathcal Z_B^{a_1\ldots a_n}(t_0,\ldots,t_n)
 =J_B^*e^{-t_0h_B}O_{B,a_1}e^{-t_1h_B}\cdots
                         O_{B,a_n}e^{-t_nh_B}J_B.             \tag{23}
\]

There is no retained projection between insertions. Excursions through
omitted local channels are therefore still in this data. Include adjoints
and every ordered word; the zero-word is (18) with j=J_B. Each factor has
the elementary norm bound <=product_i||O_(B,a_i)|| because h_B>=0.

The norm-convergent finite-volume Dyson expansion reads

\[
 \mathcal J^*e^{-t(\sum h_B+V_\partial)}\mathcal J
 =\sum_{n\ge0}(-1)^n
   \int_{\substack{t_i\ge0\\\sum_{i=0}^nt_i=t}}
   \sum_{\nu_1,\ldots,\nu_n}g_{\nu_1}\cdots g_{\nu_n}
   \bigotimes_B\mathcal Z_B^{\nu_1\ldots\nu_n}
                                      (t_0,\ldots,t_n)\,d\mathbf t.
                                                               \tag{24}
\]

An identity entry means no insertion in that block, with adjacent heat
intervals merged. Here dmathbf t means dt1...dtn with t0=t-sum_(i=1..n)ti;
the ordered simplex has volume t^n/n!. At finite volume
bounded perturbation theory gives the bound (t||V_boundary||)^n/n! for the
unexpanded Dyson term. This proves convergence and the equality (24).
Boundary gauge projection can be taken afterward because the full sum is
equivariant. Additional outer boundary insertions give the same construction
for the joined block's new Z words. This is an exact recursive composition
law for the full boundary-history family.

A further retained isometry L acts on the complete joined word by endpoint
compression L* Z_join L. It must not be inserted separately between every
boundary operation. Selecting successive isometries that keep the actual
joined vacuum and its estimates uniformly is still a dynamical obligation.

In general the zero-word C_B, or even F_B(z) at all energies, does not contain
the matrix elements of boundary operators between hidden states. The word
family is the further data needed for an interacting join. Its ordered
insertions are not replaced by products of separately projected ports.

This closes the **algebraic** join rule, not its uniform estimates: the
displayed finite-volume Dyson norm involves an extensive boundary interaction.
No volume-independent convergence constant or contraction across infinitely
many spatial scales is inferred from that norm bound.

## 8. Connection to the earlier local-coefficient route and remaining target

YC36 controls specified finite source counts and local coefficient errors
inside the admitted parent windows. It can supply local ingredients for a
fixed finite boundary word. It does not yet give a bound uniform in both the
word length and the number of joined blocks. The count constants depend on
j, so an all-word sum cannot simply use the fixed-j estimate as a constant.

The concrete next estimate is a bound on the connected boundary-return
kernels generated by (23)-(24), after vacuum energy is centered and the
source metric is restored. It must control all three of the following at
once: internal energy residuals (13), crossing boundary excursions, and
their mixed histories. Only then can one decide whether a spatial iteration
contracts in an admitted interaction norm.

There is a simple scale warning even with (8): for geometrically regular
blocks the number of crossing faces can grow like b^(2/3). A crude squared
boundary norm divided by a hidden floor of order b then grows like b^(1/3).
The O(1/b) local energy-memory estimate does not cancel that cost by itself.
Connectedness and source structure, rather than multiplication of onsite
gap floors, are the remaining issue. This is a diagnostic, not a no-go.

## 9. Zero specialization, domains and status

At ell=0, h_B=H0,B, E_B=0 and [P_B,H0,B]=0. Thus T_B=0, M_B=I,
J_B is the ordinary retained inclusion, and F_B(z)=A_B-zI exactly. Local
support energies add before the response is taken, so (19) recovers that
additive zero reading. The boundary interaction coefficients vanish; their
operators and ordered source jets remain in (23). A projected boundary
word can still visit Q between ports, so zero interaction does not authorize
erasing its response derivatives. The tensor regrouping, Gauss matching,
source, adjoint and physical pairing specialize consistently.

The auxiliary t is inverse-energy Euclidean semigroup time. Its ordered
composition follows from the represented generator; it has not been derived
as an external physical clock from tower depth. Under h_phys=c_E h and
z_phys=c_E z, use t_phys=t/c_E. M and J are dimensionless, while every term
of F and its remainder scales by c_E. No physical energy conversion or
continuum trajectory is selected by this consistency check.

The identities use standard spectral calculus and bounded perturbation
theory. The energy-dependent Schur map and derivative metric are already
in YC20; a primary reference for the classical Feshbach mechanism is
[Dusson--Sigal--Stamm, arXiv:2105.02058](https://arxiv.org/abs/2105.02058).
Projected memory is also studied in
[Widder--Zimmer--Schilling, arXiv:2503.20457v1](https://arxiv.org/abs/2503.20457v1);
no general orthogonal-dynamics existence theorem from that literature is
assumed here. Equations(18)-(24) use bounded compressed heat kernels directly.

This note contains the theoretical derivations. No scripts or certificate
packet were added or run for YC37. Practical harmonic computation, a uniform
connected-history norm and its repeated spatial bound remain outstanding.
The existing anisotropic lattice windows are unchanged; the4D continuum
mass gap is still open.
