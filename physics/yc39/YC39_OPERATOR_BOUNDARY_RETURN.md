# YC39 — a vacuum-aligned cut and the complete operator boundary excursion

10 October 2026. Research owner: Monty Dabas.
Published parent: `4be87b04122890040720dc3acada9d71e7dae3ce` (YC38
snapshot, including the recoverable original local history).

**Result.** The next intermediate construction between YC38's scalar
source bound and a full spatial iteration has three parts:

1. Justify YC37's vacuum-aligned graph cut on the actual smooth Wilson
   operator domain. Its orthogonal hidden block retains the floor cb.
2. Lift YC38's source Gram to arbitrary bounded, noncommuting operators
   on an environment, including arbitrary environment entanglement.
3. Sum every boundary insertion that remains in this hidden block. A
   boundary/volume ratio q=O(b^(-1/3)) controls the entire inverse.

For the exact vacuum-anchored compression of the high-channel Schur
return, the resulting operator inequality is

\[
 0\preceq\mathfrak R_B(z)
 \preceq {C_B\over cb-v-z}\sum_p A_p^*A_p,
 \quad C_B=\min\{N_\partial,C_\Gamma\},\quad z<cb-v.       \tag{1}
\]

Here A_p are environment operators, v bounds the complete boundary
interaction, and C_Gamma is YC38's block-size-independent covariance
constant. There is no gap assumption on the environment. The block
input in (1) is its actual isolated vacuum; the environment input is
arbitrary. This is an operator-valued, all-hidden-insertion result, not
an all-retained-input theorem or a bound for histories that repeatedly
visit retained excited channels. Those distinctions remain essential.

For regular 3D blocks and fixed port strength, (1) is O(b^(-1/3)); its
hidden lifted-response norm squared is O(b^(-4/3)). The complete hidden
inverse has a convergent series with an explicit remainder for sufficiently
large b. No numerical crossover, new global gap window or continuum limit
is asserted. This is a theory-first derivation, with additional executable
certification deferred as requested.

## 1. Existing carrier and the precise additional domain argument

Retain YC27/37's complete b-factor block, b>=8, all harmonics, boundary
charges, Haar pairing and internal Gauss action. Its actual isolated
Hamiltonian h=H0+V_int-E_B is nonnegative and has ground psi. H0 is the
sum of complete reference-factor Hamiltonians. V_int is a bounded smooth
Wilson multiplier on a finite compact product of SU(2) groups.

Let `P=1_[0,floor(b/8)](N)`, Q=I-P. The reference count commutes strongly
with H0: its individual vacuum projections commute with the corresponding
factor Hamiltonians. Thus P and Q preserve the domain of H0 and its form
domain. The original block decomposition is

\[
 h=\begin{pmatrix}A&T\\T^*&D\end{pmatrix},\quad
 D\ge dI,\quad d=cb,\quad T=PV_{\rm int}Q,\quad
 Y=D^{-1}T^*,\quad K=A-TD^{-1}T^*\ge0.                    \tag{2}
\]

YC37 gives c=47/800 or 1337/22800 at the inherited profile endpoints.
For smaller couplings use their corresponding c or these lower bounds.
All statements are in the inherited lattice electric units.

YC37 deliberately did not assume that an arbitrary graph projection
preserves an unbounded Hamiltonian's domain. Here that point is supplied
for this particular operator. On the complete compact finite block,
Dom(H0) is the second Sobolev space, with equivalent graph norm: the
electric Laplacian is elliptic and the internal potentials are smooth
and bounded. Multiplication by V_int maps this domain continuously into
itself. Commuting P,Q restrictions have the same graph-domain property.
Since D differs from QH0Q by a bounded operator, D^-1 maps QH boundedly
into Q Dom(H0). Consequently

\[
 Y:\mathcal H_P\longrightarrow\operatorname{Dom}D,
 \qquad Y^*=TD^{-1}:\mathcal H_Q\longrightarrow\operatorname{Dom}A,
 \qquad DY=T^*,\quad AY^*\text{ is bounded}.              \tag{3}
\]

The last boundedness is in the fixed finite block, by graph-norm
composition. No uniform-in-b bound on AY* is being imported into (1).
The domain argument can be performed before internal Gauss restriction;
all the operators then intertwine that action, so restriction preserves
it. For a generic abstract bounded T, (3) need not hold. The smooth
Wilson multiplication and reference spectral cut are actual hypotheses.

Define

\[
 M=I+Y^*Y,\quad N=I+YY^*,\quad
 J=\binom{I}{-Y}M^{-1/2},\qquad
 L=\binom{Y^*}{I}N^{-1/2}.                               \tag{4}
\]

J is YC37's retained isometry. L here denotes the complementary
isometry, not the differently typed residual map named L in YC37(13).
The usual graph identities give

\[
 J^*J=L^*L=I,\quad J^*L=0,\quad JJ^*+LL^*=I.             \tag{5}
\]

These are exact Hilbert-space pairings. To check domains, for either
sign write (I+Z)^(+/-1/2)-I=Z f_+/- (Z), where f is bounded and
continuous on the compact spectrum of the bounded positive Z. Equation
(3) implies that M^(+/-1/2)-I maps into Dom(A), and
N^(+/-1/2)-I maps into Dom(D). Together with (3), the unitary (J,L)
and its adjoint preserve the respective full operator domains. They
also preserve form domains, either by interpolation or the same bounded
graph-norm maps. This justifies the orthogonal cut and closed forms
used below, instead of assuming those properties from a picture.

## 2. The rotated hidden floor survives with no condition-number loss

YC37's completed-square identity is

\[
 h[(u,v)]=\langle u,Ku\rangle
          +\langle v+Yu,D(v+Yu)\rangle.                  \tag{6}
\]

For a complementary vector Lw, write v=N^(-1/2)w, so u=Y*v and
v+Yu=Nv. Applying (6) gives the exact closed-form expression

\[
 \widehat D:=L^*hL
   =N^{-1/2}YKY^*N^{-1/2}+N^{1/2}DN^{1/2}
   \succeq dN\succeq dI.                               \tag{7}
\]

The first term is a bounded positive operator by (3). The second is a
positive closed form on its transported domain. In fact the smoothing
identities above identify its operator domain with Dom(D): writing
N^(1/2)=I+R gives DR bounded, hence RD has a bounded extension as well,
and N^(1/2)DN^(1/2)=D+DR+RD+RDR on Dom(D). Equations
(3)-(7) therefore give a self-adjoint hidden generator with floor cb
on the actual orthogonal complement of Ran(J).

The isolated ground satisfies psi=J omega, omega=J*psi. Hence
L*psi=0 and L*hJ omega=0 exactly. The retained generator is
h_ret=J*hJ=M^(-1/2)KM^(-1/2), with its old onsite gap and metric intact.
The internal retained-hidden coupling on other inputs is still present;
we have not made Ran(J) an invariant subspace of h.

More precisely the transformed off-diagonal is

\[
 C:=J^*hL=M^{-1/2}KY^*N^{-1/2},\qquad
 C^*=Yh_{\rm ret}\quad\hbox{on }\operatorname{Dom}h_{\rm ret}.
                                                               \tag{7a}
\]

KY* is bounded by (3), so C is bounded at every fixed block. The second
identity follows from N^(-1/2)Y=YM^(-1/2). It also explains YC37's
low-energy residual without discarding it: for the spectral cut
`Pi_E=1_[0,E](h_ret)`, E>=0,

\[
 \|C^*\Pi_E\|\le\beta E,\qquad
 \beta=a/c\ge\|Y\|.                                    \tag{7b}
\]

Thus the internal source on low retained energies is already controlled
uniformly in b. The missing uniform estimate is the boundary source on
those excited inputs, rather than the existence of the orthogonal cut.

This is the new cut construction needed for the next step. YC38 used
the original count complement Q, for which Qpsi need not vanish. The
orthogonal graph complement now kills the actual vacuum while retaining
the large hidden floor. Changing the cut requires this argument; it is
not licensed by replacing Q by a similarly named symbol.

## 3. Join an arbitrary environment and retain its operator sources

Let the environment have a Hilbert space E and a self-adjoint generator
h_E>=0 on its stated domain. At finite volume its ground energy can be
subtracted and recorded separately; no positive environment gap or
factorization of its states is required. The actual boundary interaction
has a finite expansion

\[
 V_\partial=\sum_{p=1}^{N_\partial}O_p\otimes A_p,
 \quad V_\partial=V_\partial^*,\quad
 \|O_p\|\le1,\quad\|A_p\|\le g,\quad
 \|V_\partial\|\le v\le gN_\partial.                    \tag{8}
\]

Use exactly YC38's bounded-density local block-port assumptions for
O_p. Environment operators need not commute. Adjoint pairs and every
boundary matrix-index label are retained; neither each term's gauge
invariance nor its self-adjointness is required separately. The full
sum is gauge invariant, and the projections are gauge equivariant.
Physical matching is imposed jointly, rather than deleting charged
ports. A finite Wilson boundary has the required bounded local
decomposition, with its normalization constants in the A_p.

On the tensor carrier, relative to the carried isolated energies, let

\[
 H=h\otimes I+I\otimes h_E+V_\partial,
 \quad P_J=JJ^*\otimes I,\quad Q_J=LL^*\otimes I.          \tag{9}
\]

In hidden L coordinates its diagonal block is

\[
 \mathbb D=\widehat D\otimes I+I\otimes h_E+V_{QQ},
 \quad V_{QQ}=(L^*\otimes I)V_\partial(L\otimes I),
 \quad\mathbb D\succeq(d-v)I.                           \tag{10}
\]

The fixed-block off-diagonal internal coupling is bounded by (3)-(4),
although its norm has not been bounded uniformly in b. Tensoring with
the environment and adding bounded V_partial thus gives the ordinary
domain-compatible block resolvent/Schur construction for real z<d-v.

Let iota:E->P H tensor E send xi to omega tensor xi. The exact full
off-diagonal source restricted to this input is

\[
 S=(L^*\otimes I)H(\psi\otimes I)
   =(L^*\otimes I)V_\partial(\psi\otimes I).
                                                               \tag{11}
\]

The identity holds first on Dom(h_E) and extends as the bounded source
on the right. Both h psi and L*psi vanish, so internal and environment
free terms contribute no source here. These vanishings do not hold for
a general retained excited block state.

Define the exact vacuum-anchored Schur return

\[
 \mathfrak R_B(z)=S^*(\mathbb D-z)^{-1}S.                 \tag{12}
\]

If F_J(z) is the full Schur pencil on Ran(J) tensor E, its compression
to iota is exactly

\[
 \iota^*F_J(z)\iota
   =h_E+\sum_p\omega_p A_p-zI-\mathfrak R_B(z),\qquad
 \omega_p=\langle\psi,O_p\psi\rangle.                    \tag{13}
\]

The other retained block channels have not been eliminated by this
compression. Therefore (13) is not automatically the full effective
environment Hamiltonian. The means omega_p A_p are explicitly present;
centering the source does not discard physical one-sided terms.

## 4. Tensor amplification of the actual source Gram

Put v_p=(O_p-omega_p I)psi as in YC38. Since L*psi=0, the columns of
S are L*v_p with environment coefficients A_p. Let

\[
 G_L=(\langle L^*v_p,L^*v_q\rangle)_{pq},\quad
 G=(\langle v_p,v_q\rangle)_{pq},\quad
 \mathcal A\xi=(A_p\xi)_p.                              \tag{14}
\]

YC38 gives 0<=G<=C_B I with C_B=min(N_partial,C_Gamma). Orthogonal
projection yields 0<=G_L<=G without an entrywise decay assumption on
G_L. Tensoring a positive matrix inequality with the identity preserves
it, even on an infinite-dimensional environment. Consequently

\[
 S^*S=\mathcal A^*(G_L\otimes I_E)\mathcal A
       \preceq C_B\mathcal A^*\mathcal A
       =C_B\sum_p A_p^*A_p.                            \tag{15}
\]

Together with (10), this proves (1). Noncommuting coefficient products
and all cross terms are accounted for before the upper bound is taken.
No replacement of A_p by expectation values has occurred. The same
proof holds after tensoring an arbitrary additional spectator space,
so the bound covers environment states entangled with such spectators.
It does not allow entanglement of an excited block input with E while
still calling the block input its one-dimensional vacuum.

For all integers k>=0,

\[
 \mathfrak R_B^{(k)}(z)
   =k!S^*(\mathbb D-z)^{-k-1}S
   \preceq {k!C_B\over(d-v-z)^{k+1}}\sum_p A_p^*A_p.
                                                               \tag{16}
\]

At k=1 the left side is the Gram of the actual hidden lifted response
(mathbb D-z)^(-1)S. It is retained as metric information, not normalized
away. This is the operator-valued successor to YC38's scalar-source
estimate on a different, now justified, vacuum-aligned cut.

## 5. Sum every insertion within the hidden channel

Let D0=widehat D tensor I+I tensor h_E and, for real z<d, define

\[
 X_z=(D_0-z)^{-1/2}V_{QQ}(D_0-z)^{-1/2},\qquad
 q_z={v\over d-z}.
\]

If q_z<1, then ||X_z||<=q_z and the exact inverse is

\[
 (\mathbb D-z)^{-1}
   =(D_0-z)^{-1/2}\sum_{n=0}^{\infty}(-X_z)^n
                          (D_0-z)^{-1/2}.                \tag{17}
\]

This keeps the ordering and all environment operators. It sums every
boundary word staying in the orthogonal hidden channel between the two
source endpoints. It does not insert retained projections and drop the
resulting omitted terms: Q_J is the specified excursion channel itself.

If R_B^[m] keeps n=0,...,m, the remainder has the two-sided form bound

\[
 -E_m\preceq\mathfrak R_B(z)-\mathfrak R_B^{[m]}(z)
       \preceq E_m,\quad
 E_m={C_B\over d-z}{q_z^{m+1}\over1-q_z}\sum_p A_p^*A_p.
                                                               \tag{18}
\]

Indeed the remainder of (I+X_z)^(-1) has norm at most
q_z^(m+1)/(1-q_z); sandwich it with (D0-z)^(-1/2)S and use (15).
The finite truncation is not asserted to be positive. Positivity belongs
to the full inverse, while (18) records the actual truncation error.

## 6. A scale that makes this entire excursion small

For regular 3D blocks assume N_partial<=C_partial b^(2/3), retaining the
finite Wilson port multiplicity in C_partial. Hold g, c, C_partial,
C_Gamma and an energy window |z|<=z_* fixed. If

\[
 b\ge {4z_*\over c},\qquad
 b\ge\left({4gC_\partial\over c}\right)^3,               \tag{19}
\]

then z_*<=cb/4, v<=cb/4, q_z<=1/3 and
mathbb D-z>=cb/2. The following uniform environment-operator bounds
therefore follow:

\[
 \boxed{\|\mathfrak R_B(z)\|
       \le {2g^2C_\Gamma C_\partial\over c}\,b^{-1/3},}
 \qquad
 \boxed{\|\mathfrak R'_B(z)\|
       \le {4g^2C_\Gamma C_\partial\over c^2}\,b^{-4/3}.}
                                                               \tag{20}
\]

Thus arbitrarily many **hidden** boundary insertions no longer carry an
uncontrolled extensive Dyson factor in this one-block excursion. The
denominator grows as volume, while the complete boundary perturbation
grows as area. The endpoint Gram converts its squared source cost to
an area cost before this ratio is used.

No numerical C_Gamma or practical block size has been supplied. Its
inherited locality constants are conservative. The count cap grows
with b and its retained space still contains infinitely many harmonic
states. Equation (20) is not a fixed-dimensional computational scheme,
a small-lattice numerical result or a physical continuum scaling law.

## 7. Exact control, transport and zero reading

A two-channel control makes the cut and operator ordering explicit:

\[
 h=\begin{pmatrix}1/4&-1/2\\-1/2&1\end{pmatrix},\quad
 Y=-1/2,\quad\psi=J={1\over\sqrt5}\binom2{1},\quad
 L={1\over\sqrt5}\binom{-1}{2}.
\]

Here K=0, D=1 and L*hL=5/4>=1. For the block port
Z=diag(1,-1), omega=3/5, L*Zpsi=-4/5 and L*ZL=-3/5.
With any bounded self-adjoint environment A and h_E>=0,

\[
 \mathfrak R(z)={16\over25}
 A\left({5\over4}I+h_E-{3\over5}A-zI\right)^{-1}A.        \tag{21}
\]

The A factors cannot generally be commuted through h_E. This displayed
algebra verifies the domain-free finite control and the source/mean/
hidden-operator placement. It is not a Yang-Mills simulation or a new
executed numerical certificate.

All graph maps, cuts, ports, states and pairings transform together
under a unitary frame change. For a change of source labels the exact
Gram congruence from YC38 remains; environment coefficients transform
dually so the represented interaction is unchanged. The general bound
with sum A_p*A_p is stated in the normalized local-port presentation;
an arbitrary nonunitary source change must carry its metric cost.

At zero internal mixed coupling, Y=0, J and L are the original P and Q
inclusions, psi=Omega and widehat D=D. Thus the new construction recovers
the original cut and source. If floor(b/8)>=k_0, a single centered local
port on the product reference cannot enter Q, so S and this high-channel
return vanish exactly. Histories first passing through retained excited
channels can still reach Q after later insertions; they are not claimed
to vanish. At zero boundary strength, S=0 and R=0 for every internal
coupling; internal excited-channel residuals and source derivatives
remain in the full family.

There is also a stronger reference check for the next low-energy target.
At zero internal mixed coupling, h_ret is the restriction of H0 and
H0>=N_count/2, with commuting count. Thus Pi_E has count at most
floor(2E). A k_0-factor port can increase that count by at most k_0.
If floor(b/8)>=floor(2E)+k_0, its high-channel source B_E in section 8
vanishes identically, including operator environment coefficients.
The internal source is also zero because Y=0. This establishes the
reference value of that construction; uniform control away from zero
still requires the excited-input estimate stated there.

Under h_phys=e h, h_E,phys=e h_E, V_phys=e V and z_phys=e z,
the dimensionless graph maps are unchanged, q is unchanged and
R_phys(ez)=e R(z). The response derivative metric is dimensionless.
All isolated ground-energy shifts in (9),(13) must be carried when
comparing physical spectra. Block size is an actual cut parameter,
not an observer coordinate change or a continuum spacing by itself.

## 8. What this opens, and the exact next circle

YC38 did not justify replacing its scalar coefficients by arbitrary
environment operators, did not use the graph-orthogonal hidden inverse,
and did not sum hidden boundary interactions. Sections 1-6 provide
these three constructions, reusing its Gram rather than restarting.

The outstanding part is now precisely the retained-channel dynamics:
matrix elements with arbitrary excited block inputs, transitions among
them, and histories leaving and reentering the hidden channel through
those retained states. For such inputs the internal residual in YC37(13)
does not vanish. Neither (15) nor the vacuum-port covariance proof supplies
a uniform bound on those inputs. The block Z-port counterexample in
YC38 remains a valid warning against that promotion.

The next candidate is an energy/source-count resolved extension on the
retained block, carrying both its actual interaction and returned metric.
There is a precise intermediate target before a full iteration. Define
the boundary source on the low retained space by

\[
 B_E=(L^*\otimes I)V_\partial(J\Pi_E\otimes I).
\]

Seek an estimate, in the same normalized port presentation,

\[
 B_E^*B_E\preceq C(E)\left(\Pi_E\otimes\sum_p A_p^*A_p\right),
 \qquad C(E)\text{ independent of }b.                   \tag{22}
\]

This is an open estimate, not a consequence of the vacuum covariance
bound. Below the inherited retained gap, Pi_E is just the vacuum and
(15) supplies C(E)=C_Gamma, with the sharper pointwise C_B. For excited
Pi_E, scalar vacuum expectations
cannot replace the operator source in (22).

If (22) is established, the complete source on this input is
S_E=C*Pi_E tensor I+B_E. Equation (7b), together with the squared-norm
inequality for a sum, then gives, for every eta>0,

\[
 0\preceq S_E^*(\mathbb D-z)^{-1}S_E
 \preceq\frac{
 (1+\eta)\beta^2E^2(\Pi_E\otimes I)
 +(1+\eta^{-1})C(E)\left(\Pi_E\otimes\sum_p A_p^*A_p\right)
 }{cb-v-z}.                                             \tag{23}
\]

This conditional estimate includes the mixed internal/boundary cross
terms. At fixed E it would have an internal O(b^(-1)) contribution and
boundary O(b^(-1/3)) contribution under (19). Both sources and their
returned metrics remain available; (23) does not erase either one.

The resulting weighted norm must still control repeated retained/hidden
alternation and composition across multiple blocks. Smallness of the
one excursion (20), or even a proof of (22), alone does not prove that
the low retained dynamics is gapped or that its interactions contract
at each scale. Existing anisotropic lattice gap windows are inputs and
remain unchanged; the 4D continuum task is open.

**Verification performed:** written domain, completed-square, tensor-Gram,
resolvent-series, zero-specialization and energy-unit derivations, with
the exact finite control (21). No new executable certification script,
routine predecessor suite or numerical crossover calculation was run.
Source pins are in [SOURCES.json](SOURCES.json). The spectral/Schur and
positive-matrix methods are standard; no external priority claim is made.
