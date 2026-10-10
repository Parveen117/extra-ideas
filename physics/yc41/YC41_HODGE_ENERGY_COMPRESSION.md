# YC41 — a Hodge-type split and controlled energy compression

10 October 2026. Research owner: Monty Dabas.
Parent: c0ea33bc26bb9bf8e6452cb628e29469b0d84514 (YC40).

**Result.** The first closure problem in YC40 can be addressed without
requiring the compressed ports to remain scalar multiplication
operators. Their actual energy-form bounds survive compression:

\[
 \boxed{\|h_B^{1/2}J_BJ_B^*x\|
       \le\Lambda_B\|h_B^{1/2}x\|,\qquad \Lambda_B\le7/6.}       \tag{1}
\]

This holds on the entire form domain in the inherited 28/64-link
windows, at b>=8. The proof supplies a sharper block-dependent
Lambda_B below. The bound concerns the energy norm, whose square is
the energy form.

Three consequences follow.

1. Compression carries an adjoint-stable algebra of energy-controlled
   operators into the corresponding retained algebra. Generated ports
   can leave the scalar Wilson class while keeping an explicit bound.
2. A tensor layer of any number of independent blocks uses the maximum
   of their Lambda_B, not their product.
3. The failure of compression to preserve ordered products is an
   explicit hidden-source pairing. At fixed low-energy endpoints,
   each fixed-length word has O(b^(-1)) error from inserting retained
   projections. A stated range of logarithmically growing word lengths
   is controlled too.

The construction has an exact Hodge-type interpretation: the hidden
gradient range is closed, its orthogonal projection removes the
eliminable gradient, and the remainder is the Schur energy. This is a
two-term closed-range decomposition, not a claim about de Rham
cohomology or a continuum Yang--Mills Hodge theorem.

This stage proves a controlled spatial compression layer and an
operator class stable under its stated algebra operations. It does
not yet prove uniform contraction through arbitrarily many layers,
unrestricted time-dependent histories, or a continuum mass gap.
Frozen predecessor notes and certificates remain unchanged.

## 1. The carrier and the two different projections

Retain YC37-40's complete block, charge matching, Haar pairing and
actual isolated Hamiltonian h=H0+V-E_B>=0. Here H0 is the sum of
reference-factor Hamiltonians above their respective ground energies,
and V is the internal mixed Wilson potential. The inherited data are

\[
 \|V\|\le ab,\quad -ab\le E_B\le0,\quad
 d=cb,\quad \beta=a/c,\quad D=QhQ\ge dI.
\]

P is the original reference count cut, Q=I-P; both commute with H0.
Use YC39's domain-compatible graph maps

\[
 Y=D^{-1}T^*,\quad K=A-TD^{-1}T^*,\quad
 W=\binom I{-Y},\quad M=I+Y^*Y,\quad N_g=I+YY^*,
\]
\[
 J=WM^{-1/2},\quad L=\binom{Y^*}I N_g^{-1/2},\quad
 P_J=JJ^*,\quad Q_J=LL^*,\quad h_r=J^*hJ.                    \tag{2}
\]

Here A=PhP, T=PVQ and K>=0. The original count complement Q and
the orthogonal graph complement Q_J are distinct. In particular,
the gradient projection in section 2 uses Q, while (1) concerns P_J.
The complete-square identity and YC40's source-energy inequality are

\[
 h[(u,w)]=K[u]+\|D^{1/2}(w+Yu)\|^2,\qquad h\ge dQ_J.        \tag{3}
\]

The old domains and boundary gauge actions are carried throughout.
All form statements extend to arbitrary tensor spectators and hence
to entangled inputs. Physical restrictions are imposed by the
commuting gauge actions, not by replacing each boundary factor by
its singlet sector.

These represented spaces carry the admitted operators, gradient map,
adjoints, composition, source cuts, gauge action and returned pairing.
The mixed coupling, spatial block size, energy cap and algebra weight
tau are distinct adapters or choices on that carrier, not primitive
lambda-space or elapsed time.

## 2. The hidden gradient and its harmonic remainder

Let psi be the actual positive block ground and dnu=psi^2 dmu.
Define the closed ground-state gradient from the original Hilbert
space to the gradient Hilbert space by

\[
 {\cal G}x=\bigl(X_\mu(x/\psi)\bigr)_\mu,\qquad
 {\cal F}_1=\bigoplus_\mu L^2(d\nu),\qquad
 h={\cal G}^*{\cal G}.
\]

Its domain is the form domain of h, by YC22/40. For the original
hidden inclusion i_Q put G_Q=G i_Q. It is closed, and
G_Q*G_Q=D>=dI. Therefore

\[
 U_Q=G_QD^{-1/2},\quad U_Q^*U_Q=I_Q,\qquad
 {\cal P}_Q=U_QU_Q^*
\]

is the orthogonal projection onto the closed range of G_Q. In particular

\[
 {\cal F}_1=\operatorname{Ran}G_Q\ \mathbin{\oplus}\ \ker G_Q^*.
                                                               \tag{4}
\]

The expression G_Q D^-1 G_Q* for this projection is understood as
the bounded extension U_Q U_Q*. Closed range follows directly from
the lower bound on D and closedness of G_Q.

For u in the retained form domain, the off-diagonal form gives
G_Q*G_P u=T*u, where G_P=G i_P. Thus

\[
 {\cal G}Wu=(I-{\cal P}_Q)G_Pu,\qquad
 K[u]=\|(I-{\cal P}_Q)G_Pu\|^2,\qquad
 {\cal G}_r={\cal G}J,\quad {\cal G}_r^*{\cal G}_r=h_r.       \tag{5}
\]

The word harmonic here means that the lifted gradient is orthogonal
to all gradients from the chosen hidden channel. It does not mean
that every lifted vector has zero physical energy. The returned
Hilbert metric remains M=W*W; J restores the physical pairing exactly.

The related oblique projection

\[
 {\cal E}=WP
\]

satisfies E^2=E and h[E x]=K[Px]<=h[x] by (3). Hence it is an exact
contraction in the homogeneous energy seminorm. Its orthogonal
counterpart P_J requires the additional quantitative estimate below.
Neither projection is substituted for the other without its metric.

## 3. A local derivative bound for the internal mixed potential

In YC22's unit-S^3 electric convention put

\[
 \Gamma_V=\left\|\sum_\mu|X_\mu V|^2\right\|_\infty.
\]

For the inherited simple cubic plaquette family, each link belongs
to at most four distinct unoriented plaquettes, each mixed face
meets four reference factors, and the face-factor incidence is
I_face=48 or88. A normalized fundamental plaquette has gradient
cost at most12, as in YC40. Operator/scalar Cauchy--Schwarz gives

\[
 n_{\rm int}\le I_{\rm face}b/4,\qquad
 \boxed{\Gamma_V\le\zeta b,\qquad
        \zeta=12\ell^2 I_{\rm face}.}                       \tag{6}
\]

All coefficients xi_p are in [0,1]. The estimate uses the admitted
torus sizes, without repeated-link plaquettes. A different interaction
list must carry its actual incidence and derivative budget.

Let H_ref=H0-E_B>=0. P,Q commute with H_ref. The product of the
reference-factor positive grounds supplies the exact H0 ground form.
Adding the nonnegative constant -E_B to that form and using the
product rule gives

\[
 \|H_{\rm ref}^{1/2}Vx\|
 \le ab\,\|H_{\rm ref}^{1/2}x\|+\sqrt{\Gamma_V}\,\|x\|.
                                                               \tag{7}
\]

One can regard the constant shift as an extra zeroth-order component
of the form map; multiplication by V commutes with that component.
No bound on a ground's logarithmic derivative is used.

For w in the hidden Hilbert space,

\[
 H_{\rm ref}[D^{-1}w]\le
     \left(\frac1d+\frac{ab}{d^2}\right)\|w\|^2
     =\frac{1+\beta}{d}\|w\|^2.
\]

Since Y*=PVQ D^-1, P contracts the H_ref form and
A<=P H_ref P+ab I. Equations (6)-(7) imply

\[
 A[Y^*w]\le d\,\kappa_b^2\|w\|^2,\qquad
 \boxed{\kappa_b^2=
 \left(\beta\sqrt{1+\beta}
       +\frac{\sqrt\zeta}{c^{3/2}b}\right)^2+\beta^3.}       \tag{8}
\]

For clarity, the first squared term comes from (7) applied to D^-1w.
The additional ab||Y*w||^2 is at most beta^3 d||w||^2.
The metric identity M^(-1/2)Y*=Y*N_g^(-1/2), together with K<=A,
then gives the useful transported estimate

\[
 \boxed{\|h_r^{1/2}Y^*\|\le\kappa_b\sqrt d.}                \tag{9}
\]

This is a form-regularity bound for the actual graph correction.
A bound on ||Y|| alone would not imply it.

## 4. The orthogonal graph projection has a uniform energy bound

The exact graph algebra yields

\[
 {\cal E}=P_J+JY^*L^*.
\]

Combine its energy contraction with (9) and h>=dQ_J:

\[
 \begin{aligned}
 \|h^{1/2}P_Jx\|
 &\le\|h^{1/2}{\cal E}x\|+\|h_r^{1/2}Y^*L^*x\|\\
 &\le\|h^{1/2}x\|+\kappa_b\sqrt d\,\|L^*x\|\\
 &\le(1+\kappa_b)\|h^{1/2}x\|.
 \end{aligned}                                             \tag{10}
\]

Thus Lambda_b=1+kappa_b is an explicit admissible constant in (1).
For both old profile windows, the elementary common estimates are

\[
 \beta\le1/15,\quad c\ge1/18,\quad
 \sqrt\zeta\le3/400,\quad b\ge8.
\]

Using sqrt(18)<5 gives

\[
 \kappa_b^2
 \le(16/225+27/320)^2+1/3375<1/36,\qquad \Lambda_b\le7/6.
                                                               \tag{11}
\]

These are written rational bounds from the old parameters, not a new
numerical certificate. They are deliberately loose. At zero internal
mixed coupling beta=zeta=0, the exact constant becomes Lambda=1.

## 5. Generated operators retain their energy budget

Say a bounded O has a two-sided energy budget (m_O,a_O,b_O) when
||O||<=m_O and both O and O* preserve the form domain and satisfy

\[
 \|h^{1/2}Ox\|\le a_O\|h^{1/2}x\|+b_O\|x\|,               \tag{12}
\]

with the same bound for O*. A normalized smooth Wilson multiplication
port has m=a=1 and b=sqrt(gamma_O), by the ground-state product rule.
Its adjoint has the same gradient cost.

For O_r=J*OJ, no scalar-multiplier interpretation is necessary:

\[
 \|O_r\|\le m_O,\qquad
 \boxed{\|h_r^{1/2}O_ru\|
       \le\Lambda_b a_O\|h_r^{1/2}u\|
          +\Lambda_b b_O\|u\|.}                            \tag{13}
\]

Indeed J O_r=P_J OJ and (10) applies. The adjoint is exactly J*O*J
and obeys the same estimate. Physical pairing, rather than an
unrecorded similarity, is used here.

The collective operator-coefficient form in YC40 also survives.
For V_partial=sum f_p tensor A_p and
V_r=(J* tensor I)V_partial(J tensor I),

\[
 \| (h_r^{1/2}\otimes I)V_r\xi\|
 \le\Lambda_b\left[
 \sqrt{n\,(h_r\otimes\Sigma)[\xi]}
 +\sqrt{\rho\,\langle\xi,(I\otimes\Sigma_\gamma)\xi\rangle}
 \right].                                                 \tag{14}
\]

Sigma and Sigma_gamma are YC40's actual coefficient columns. Their
ordering and entangled input pairing are preserved. Equation (14)
is an energy budget for generated operators, not a claim that their
new literal derivative incidence still equals rho.

For products, the forward budget of AB is

\[
 m_{AB}\le m_A m_B,\quad a_{AB}=a_Aa_B,\quad
 b_{AB}^{(+)}=a_A b_B+b_A m_B.
\]

The adjoint uses the reversed order, so a common two-sided b is
max(a_A b_B+b_A m_B, a_B b_A+b_B m_A). Sums add the corresponding
budgets. This gives explicit closure under finite products, sums
and adjoints, including the generated coarse ports.

A convenient complete normed version uses tau>0:

\[
 S_\tau=(h+\tau)^{1/2},\qquad
 \|O\|_{{\cal A}_\tau(h)}
 =\max\{\|O\|,\|S_\tau O S_\tau^{-1}\|,
                  \|S_\tau O^* S_\tau^{-1}\|\}.             \tag{15}
\]

Require O,O* to preserve Dom(S_tau), and use the bounded extensions
of the displayed conjugates. This is a Banach star algebra: the
product conjugates multiply in their stated order; adjoints exchange
the last two entries; completeness follows from closedness of S_tau
and convergence of all three bounded operators. Equation (12) gives
the bound max(a_O,m_O)+b_O/sqrt(tau).

Since P_J is also a contraction in Hilbert norm, (10) applies to
the shifted energy norm with the same Lambda_b. Hence

\[
 \boxed{\|J^*OJ\|_{{\cal A}_\tau(h_r)}
          \le\Lambda_b\|O\|_{{\cal A}_\tau(h)}.}             \tag{16}
\]

If a dimensionless T has this algebra norm below one, its full
Neumann series converges in the same algebra. This supplies a
resolvent rule with an explicit condition; it does not assert that
every physical interaction satisfies that condition.

## 6. An arbitrary tensor layer uses a maximum, not a product

For finitely many blocks set

\[
 H_\oplus=\sum_B h_B,\quad {\cal J}=\bigotimes_B J_B,\quad
 {\cal P}={\cal J}{\cal J}^*=\bigotimes_B P_{J,B},\quad
 H_r=\sum_B h_{r,B}.
\]

For the B energy term, every other block projection is a Hilbert
contraction commuting with h_B. Apply (10) only on B, then sum:

\[
 H_\oplus[{\cal P}x]\le
       \sum_B\Lambda_B^2 h_B[x]
       \le\Lambda_{\max}^2 H_\oplus[x],\qquad
 \Lambda_{\max}=\max_B\Lambda_B\le7/6.                      \tag{17}
\]

Therefore (13)-(16) hold for this full compression with Lambda_max,
independently of the number of blocks. The constants of an extensive
interaction may still depend on its own incidence and coefficient
column; (17) removes an additional spurious product of projection
constants. It does not make an extensive interaction small.

The compression of the fine joined Hamiltonian remains the exact
form compression described in YC37, including its transformed
boundary terms. The full complementary return is carried separately.
The reference energy in (17) is the independent-block sum; changing
to the interacting joined generator requires its own comparison.

## 7. Quantify the source lost by multiplying compressed ports

Compression is adjoint preserving, but it is not multiplicative.
For any bounded O,F, its exact defect is

\[
 (OF)_r-O_rF_r=J^*OQ_JFJ=(L^*O^*J)^*(L^*FJ).              \tag{18}
\]

For normalized ports with their energy budgets (a=1,
b=sqrt(gamma)), and Pi_E the h_r spectral cut onto [0,E], (3)
and the product rule give

\[
 \|\Pi_E[(OF)_r-O_rF_r]\Pi_E\|
 \le\frac{(\sqrt E+\sqrt{\gamma_O})
              (\sqrt E+\sqrt{\gamma_F})}{d}.               \tag{19}
\]

The pairing in (18) is kept before the norm bound. If the fine
ports commute, their retained commutator is determined by the
difference of these two omitted-channel pairings, rather than
being declared zero.

More generally let O_1,...,O_m have norm at most one and the same
two-sided fine energy budget (1,1,sqrt(gamma)). Put

\[
 e_0=\sqrt E,\qquad e_j=\Lambda_b(e_{j-1}+\sqrt\gamma).
\]

Every retained prefix has norm at most one on a unit input and
energy norm at most e_j. The exact telescoping defect is

\[
 \begin{split}
 J^*O_m\cdots O_1J-O_{m,r}\cdots O_{1,r}
 =\sum_{j=1}^{m-1}
 J^*O_m\cdots O_{j+1}Q_J O_jJ
            O_{j-1,r}\cdots O_{1,r}.
 \end{split}                                               \tag{20}
\]

On low-energy endpoints, the uncompressed left suffix has energy
norm at most sqrt(E)+(m-j)sqrt(gamma). Apply h>=dQ_J on both
sides of each hidden pairing:

\[
 \boxed{\left\|\Pi_E
 [J^*O_m\cdots O_1J-O_{m,r}\cdots O_{1,r}]
 \Pi_E\right\|
 \le\frac1d\sum_{j=1}^{m-1}
 [\sqrt E+(m-j)\sqrt\gamma]\,[e_{j-1}+\sqrt\gamma].}          \tag{21}
\]

For fixed E,gamma,m this is O(b^(-1)). Using Lambda_b<=7/6,
a sufficient condition for growing m=m_b is
m_b^2 (7/6)^(m_b)/b ->0. In particular m_b<=alpha log b
with fixed alpha<1/log(7/6) suffices.

This estimates ordered port products with specified low-energy
endpoints. It does not insert spectral cuts after every port, and
intermediate energies are tracked by the e_j recurrence. It does
not yet estimate arbitrary propagator insertions, a sum over all
word lengths, or all retained/hidden dynamical histories.

## 8. Returned operators on the stated energy window

Let R_E(z) be YC40's complete internal-plus-boundary Schur return,
already compressed on both sides to Pi_E tensor I. Use the retained
block energy h_r tensor I in (15); an environment energy weight is
not silently included. Since the range and source of R_E lie in
that energy window, for real z below its hidden floor,

\[
 \|R_E(z)\|_{{\cal A}_\tau(h_r\otimes I)}
       \le\sqrt{1+E/\tau}\,\|R_E(z)\|.                      \tag{22}
\]

The same bound holds for each self-adjoint derivative R_E^(k)(z).
Consequently YC40's O(b^(-2/3)) return and O(b^(-5/3)) derivative
metric remain small in this algebra at fixed E,tau. The stronger
vacuum bounds remain valid too.

This is an actual window-restricted member of the transported
operator class. It does not bound the unrestricted return in the
same norm or eliminate the other retained energies. Those channels
must be carried when composing the full physical problem.

## 9. Exact controls and the zero reading

An abstract two-channel control keeps the gradient and metric visible:

\[
 h=\begin{pmatrix}k+\beta^2d&\beta d\\\beta d&d\end{pmatrix},
 \quad
 G=\begin{pmatrix}\sqrt k&0\\\beta\sqrt d&\sqrt d\end{pmatrix},
 \quad k,d>0.
\]

Here h=G*G, Y=beta, W=(1,-beta)^T and K=k. The hidden gradient
projection is diag(0,1), so removing it from G i_P leaves
(sqrt(k),0)^T exactly. This verifies (5) in a finite model.

The same example shows why the extra derivative argument matters.
For x=(0,1)^T,

\[
 \frac{h[P_Jx]}{h[x]}
        =\frac{k\beta^2}{d(1+\beta^2)^2}.                   \tag{23}
\]

At fixed nonzero beta and fixed d this grows without bound as k increases. Thus
a hidden floor and a small Hilbert graph norm alone do not establish
(1). The Wilson derivative budget in (6)-(9) supplies the missing
control. This abstract matrix check does not claim that the
counterexample family obeys that Wilson budget.

At zero internal mixed coupling Y=0, M=I, J=P and Lambda=1.
The reference cut commutes with H0, the two gradient ranges are
orthogonal, and (5) becomes the ordinary retained reference energy.
For local ports of support at most k_0, every intermediate reference
word on Pi_E stays below the count cut when

\[
 \lfloor b/8\rfloor\ge\lfloor2E\rfloor+m k_0.
\]

Then the exact defect (20) is zero. Its general upper bound need
not vanish at that specialization. Source derivatives and every
complementary channel remain in the full family.

With boundary interactions switched off, internal graph residuals
can still occur on excited inputs, as recorded in YC39/40. They are
not removed by the algebra construction. All gradient projections,
Hilbert metrics, gauge actions and source pairings transform with
their declared unitary frame changes.

Under an overall energy conversion by e, gradients scale by sqrt(e),
Gamma_V by e^3 (V itself scales by e), d by e and
H_ref by e. Thus kappa_b and Lambda_b are unchanged. In the energy
algebra tau scales by e, and b_O in (12) scales by sqrt(e) for
dimensionless O. The cut, retained pairing and word-defect bounds
are unchanged. Ground-energy shifts are still carried explicitly.

## 10. Scope of the closure and the next construction

YC40 required generated operators to retain a usable energy budget.
Equations (13)-(17) now supply that budget for one complete tensor
compression layer and give a concrete adjoint-stable algebra for
the generated ports. Equations (18)-(21) quantify the source error
of multiplying those ports, including some increasing word lengths.
Window-restricted complete returns and metrics enter that same
algebra by (22).

Three obligations remain for unrestricted spatial iteration:

1. Construct subsequent coarse cuts with energy-projection constants
   controlled across depth. A bound of 7/6 at each layer can accumulate;
   it is not a contraction.
2. Carry unrestricted retained-energy transitions and returned
   operators, rather than only their fixed-window compressions.
3. Bound time-dependent ordered histories and their connected sums,
   including the change from independent-block to interacting energy.

The Hodge-type split is exact and provides a way to retain the
eliminable gradient and its orthogonal remainder at each proposed
cut. It does not supply these quantitative iteration estimates
merely by naming the decomposition. Existing lattice gap windows
and the open continuum target are unchanged.

**Verification performed:** written gradient-range and harmonic-lift
identities; internal derivative and reference-form bounds; the rational
7/6 estimate; tensor, algebra and ordered-word derivations; exact
matrix and zero controls; energy-unit transport. Source hashes and
new links were checked. No executable mathematical certificate,
numerical crossover or routine predecessor suite was run.
See [SOURCES.json](SOURCES.json).

Standard background: Leopardi--Stern,
[The abstract Hodge--Dirac operator and its stable discretization](https://arxiv.org/abs/1401.1576),
and Corach--Maestripieri--Stojanoff,
[Oblique projections and Schur complements](https://arxiv.org/abs/math/0006120).
The former develops abstract Hodge theory; the latter studies
energy-compatible projections and shorted operators in a bounded
setting. The unbounded domains and the specific Wilson estimates
used here are supplied explicitly above and in YC39. These papers
are background, not evidence for a published proof of this YM stage.
