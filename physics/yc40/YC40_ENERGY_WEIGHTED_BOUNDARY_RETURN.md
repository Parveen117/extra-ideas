# YC40 — energy-weighted Wilson sources beyond the vacuum

10 October 2026. Research owner: Monty Dabas.
Parent: 26ddc5ffaf9f1cbe257c91f25155e1bd3237c2aa (YC39).

**Result.** YC39(22)'s missing excited-input boundary estimate follows
for the actual smooth Wilson multiplication ports. The ingredients
already existed: YC22's positive-ground-state energy form and YC39's
graph cut. Their combination gives the stronger form inequality

\[
                     h_B\succeq cb\,L_BL_B^*.                 \tag{1}
\]

Thus entering this hidden channel has an energy cost on every input,
not only on vectors already lying inside the hidden channel. Apply that
inequality after the complete boundary interaction, retaining its
derivative source and the ordering of its environment coefficients.

Let n be the number of normalized boundary ports, E a retained energy
cap, rho their maximal derivative incidence and gamma_* their maximal
electric gradient cost, defined precisely below. The boundary source
obeys

\[
 B_E^*B_E\preceq C_b(E)
       \left(\Pi_E\otimes\sum_{p=1}^n A_p^*A_p\right),\qquad
 \boxed{C_b(E)=
   \frac{(\sqrt{nE}+\sqrt{\rho\gamma_*})^2}{cb}.}               \tag{2}
\]

This includes arbitrary retained states of energy at most E, their
entanglement with the environment, and noncommuting bounded A_p.
Bounded-density ports make C_b(E) uniformly bounded in block size.
For regular three-dimensional boundaries and fixed E, the **complete
return through the graph-hidden channel**, including the internal
source, is O(b^(-2/3)). On the actual block vacuum it is O(b^(-4/3)),
with derivative metric O(b^(-7/3)).

These are written results under the inherited lattice hypotheses and
the stated smooth-port class. They do not promote an arbitrary bounded
observable to a smooth Wilson multiplier. All hidden insertions are
included; repeated visits to retained excited channels and closure of
the interaction class under spatial iteration remain separate tasks.
No numerical certificate, enlarged gap window or continuum result is
claimed. Frozen predecessor notes are unchanged.

## 1. Full carrier, exact cut and the additional source structure

Keep YC27/37's complete spatial block of b>=8 reference factors, its
full link harmonics, internal Gauss action and boundary representations.
Before gauge restriction its Hamiltonian is the compact scalar
Schrodinger operator

\[
 h=H_B-E_B=-\sum_\mu X_\mu^2+V_B-E_B\ge0.
                                                               \tag{3}
\]

Here X_mu are the orthonormal real link vector fields in YC22's
unit-S^3 normalization; V_B includes all reference and mixed Wilson
potentials. If reference ground energies were already subtracted,
carry their scalar sum so that (3) still means energy above the actual
isolated block ground. The smooth real scalar potential has a unique
strictly positive normalized ground psi.

Use exactly YC39's original count cut P,Q and graph coordinates:

\[
 h=\begin{pmatrix}A&T\\T^*&D\end{pmatrix},\quad D\ge dI,\quad
 d=cb,\quad Y=D^{-1}T^*,\quad K=A-TD^{-1}T^*\ge0,
\]
\[
 M=I+Y^*Y,\quad N_g=I+YY^*,\quad
 J=\binom I{-Y}M^{-1/2},\qquad
 L=\binom{Y^*}I N_g^{-1/2}.                                  \tag{4}
\]

N_g is a graph metric, not the port count n or the reference excitation
count. YC39 established the required operator and form domains for these
maps. J,L are complementary isometries, psi=J omega and

\[
 h_{\rm ret}=J^*hJ,\quad
 \widehat D=L^*hL\ge dI,\quad
 C^*=L^*hJ=Yh_{\rm ret},\quad \|Y\|\le\beta=a/c.
                                                               \tag{5}
\]

The inherited endpoint pairs are c=47/800, beta=3/47 for the 28-link
profile, and c=1337/22800, beta=88/1337 for the 64-link profile.
No new lattice coupling window is introduced.

Take a complete boundary expansion on an environment Hilbert space:

\[
 V_\partial=\sum_{p=1}^n f_p\otimes A_p=V_\partial^*,\qquad
 \|f_p\|_\infty\le1,\quad \|A_p\|\le g,\quad
 \|V_\partial\|\le v\le gn.                                  \tag{6}
\]

The block ports are now smooth, possibly complex scalar multiplication
functions of a uniformly bounded set of block links. Their environment
coefficients can be arbitrary bounded noncommuting operators. No term
is assumed self-adjoint separately. Adjoint partners and boundary
matrix indices are retained in the sum.

Define the actual derivative budgets in this local port presentation:

\[
 \gamma_p=\left\|\sum_\mu|X_\mu f_p|^2\right\|_\infty,\qquad
 \gamma_*=\max_p\gamma_p,\qquad
 \rho=\max_\mu\#\{p:X_\mu f_p\not\equiv0\}.                    \tag{7}
\]

These are finite and independent of b for a fixed finite list of
normalized local Wilson shapes with bounded incidence. Constants are
selected from these shapes, not fitted to the desired conclusion.

A Wilson face can be expanded across the block boundary into products
of matrix entries of its open path segments, contracting internal
indices first and retaining the boundary indices. Each normalized
block factor has norm at most one and uses at most four distinct links.
In the unit-S^3 normalization each fundamental directional generator
has norm one. Differentiating a matrix entry of such a product once at
one link therefore has absolute value at most one. There are three
directions per link, so gamma_p<=12 is a safe bound for these individual
plaquette ports. Trace factors, couplings and any sums not expanded
into these normalized ports belong in the coefficients and label count.
Repeated-link loops or newly generated coarse operators require their
own derivative budgets; the number 12 is not asserted for them.

Internal Gauss invariance is preserved by contracting a path's internal
indices. The full sum preserves the joint boundary matching. All
identities may equivalently be derived before gauge restriction and
then restricted to invariant subspaces. No charged representation or
boundary multiplicity is removed.

## 2. Strengthen the hidden floor to an inequality on all inputs

For an arbitrary form-domain vector x=(u,w) in the original P,Q
coordinates, YC37's completed square is

\[
 h[x]=K[u]+\|D^{1/2}(w+Yu)\|^2.
\]

But the orthogonal hidden coordinate of that same vector is
L*x=N_g^(-1/2)(Yu+w). Since N_g>=I and K>=0,

\[
 d\|L^*x\|^2
 \le d\|Yu+w\|^2
 \le \|D^{1/2}(Yu+w)\|^2
 \le h[x].                                                  \tag{8}
\]

This proves (1) on the full form domain and after tensoring any
spectator Hilbert space. It is stronger than merely knowing
L*hL>=dI. The latter statement alone would not control mixed
retained/hidden vectors for a general block matrix. Positivity of K
and the actual graph construction are what supply (8).

The estimate is compatible with h psi=0 because L*psi=0. It does not
assert h>=dI or assign energy cb to the retained excitations.

## 3. Reuse the ground-state energy form before estimating the source

YC22(1) applies to the actual isolated block. Put dnu=psi^2 dmu.
Multiplication by psi is unitary from L2(dnu) to the Haar Hilbert space,
and, also for environment-valued functions F,

\[
 \| (h^{1/2}\otimes I)(\psi F)\|^2
       =\sum_\mu\int\|X_\mu F\|^2\,d\nu.                    \tag{9}
\]

The identity follows by integration by parts and h psi=0 on smooth
functions, then form closure. At each finite block psi is smooth and
positive. Neither a uniform lower bound on psi nor a bound on its
logarithmic derivative is used.

Write the boundary operator as the environment-operator-valued
multiplication function

\[
 {\cal V}(x)=\sum_p f_p(x)A_p,\qquad
 \Sigma=\sum_p A_p^*A_p,\qquad
 \Sigma_\gamma=\sum_p\gamma_p A_p^*A_p.
\]

Operator Cauchy--Schwarz at each configuration gives

\[
 {\cal V}(x)^*{\cal V}(x)\preceq n\Sigma,\qquad
 \sum_\mu (X_\mu{\cal V})^*(X_\mu{\cal V})
                  \preceq\rho\Sigma_\gamma
                  \preceq\rho\gamma_*\Sigma.                 \tag{10}
\]

For the derivative inequality, each fixed mu has at most rho summands;
bound the squared norm of their sum by rho times their squared norms,
then sum mu and use (7). No A_p is commuted past A_q or replaced by an
expectation value.

The exact product rule retains two contributions:

\[
 X_\mu({\cal V}F)={\cal V}X_\mu F+(X_\mu{\cal V})F.           \tag{11}
\]

The first transports the input energy; the second is the derivative
source. Treating the full operator-valued multiplier as a real scalar
in a double-commutator identity would require extra justification.
Equations (9)-(11) avoid that substitution entirely.

## 4. An all-retained-input relative form bound

Let

\[
 B=(L^*\otimes I)V_\partial(J\otimes I).
\]

For xi in the form domain of h_ret tensor I, put
F=psi^(-1)(J tensor I)xi. Combine (8), the triangle inequality for
the gradient Hilbert space, and (10)-(11). The result is

\[
 \boxed{\sqrt d\,\|B\xi\|
 \le \sqrt{\,n\,(h_{\rm ret}\otimes\Sigma)[\xi]\,}
       +\sqrt{\,\rho\,\langle\xi,(I\otimes\Sigma_\gamma)\xi\rangle\,}.}
                                                               \tag{12}
\]

In the first term, A_p acts only on the environment and commutes with
the block gradient, ground multiplication and J. Thus its integral is
exactly the h_ret tensor A_p*A_p form, including entangled xi. In the
second term J is isometric. These facts justify both tensor pairings
in (12); they are not product-state assumptions on xi.

Equivalently, for every t>0,

\[
 B^*B\preceq\frac{(1+t)n}{d}\,h_{\rm ret}\otimes\Sigma
       +\frac{(1+t^{-1})\rho}{d}\,I\otimes\Sigma_\gamma
                                                               \tag{13}
\]

as forms. This bound covers all retained inputs with the displayed
energy weight. It is not an energy-independent norm bound on the
entire unbounded retained spectrum.

Now let Pi_E be the spectral projection of h_ret onto [0,E], E>=0,
and B_E=B(Pi_E tensor I). Since

\[
 (\Pi_E h_{\rm ret}\Pi_E)\otimes\Sigma
             \preceq E\,\Pi_E\otimes\Sigma,
\]

equation (12) proves (2). For bounded port density n<=nu_0 b and b>=8,

\[
 C_b(E)\le C(E):=\frac1c
       \left(\sqrt{\nu_0E}+\sqrt{\rho\gamma_*/8}\right)^2.
                                                               \tag{14}
\]

This supplies the block-size-independent C(E) sought in YC39(22)
for the stated Wilson multiplication class. The sharper C_b(E) is
retained in all further estimates. The broad class of merely bounded
ports in YC38/39 need not satisfy (7), so its excited-input extension
is not claimed here.

On the actual vacuum the input gradient in (11) is zero. The bound is
then sharper without a two-term squared-norm loss:

\[
 B_0^*B_0\preceq\frac{\rho}{d}\,
                 |\omega\rangle\langle\omega|\otimes\Sigma_\gamma.
                                                               \tag{15}
\]

Constants in f_p have zero derivative and cause no hidden vacuum
source. This is the energy-form counterpart of centering; the
one-sided physical boundary means remain in the retained Hamiltonian.

## 5. Complete low-energy Schur return, including the internal source

Use YC39's arbitrary environment h_env>=0 and its exact hidden
diagonal, with no positive environment gap assumption:

\[
 {\mathbb D}=\widehat D\otimes I+I\otimes h_{\rm env}
          +(L^*\otimes I)V_\partial(L\otimes I),\qquad
 {\mathbb D}\ge(d-v)I.
\]

For real z<d-v write delta_z=d-v-z. The complete off-diagonal source
on the low retained input is

\[
 S_E=(Yh_{\rm ret}\Pi_E)\otimes I+B_E,\qquad
 R_E(z)=S_E^*({\mathbb D}-z)^{-1}S_E.                        \tag{16}
\]

The environment free term has no off-diagonal part because L*J=0.
The internal term is retained and satisfies
norm(Y h_ret Pi_E)<=beta E. Therefore, for every eta>0,

\[
 0\preceq R_E(z)\preceq\frac{W_E(\eta)}{\delta_z},
\]
\[
 W_E(\eta)=(1+\eta)\beta^2E^2(\Pi_E\otimes I)
       +(1+\eta^{-1})C_b(E)(\Pi_E\otimes\Sigma).             \tag{17}
\]

This establishes YC39(23) for the Wilson class, including its mixed
internal/boundary source terms. A useful sharper norm form is

\[
 \|R_E(z)\|\le
     \frac{\bigl(\beta E+g\sqrt{n C_b(E)}\bigr)^2}{\delta_z}.
                                                               \tag{18}
\]

For every integer k>=0, keeping E fixed when differentiating z,

\[
 0\preceq R_E^{(k)}(z)
       =k!S_E^*({\mathbb D}-z)^{-k-1}S_E
       \preceq\frac{k!\,W_E(\eta)}{\delta_z^{k+1}}.          \tag{19}
\]

In particular R_E' is the Gram of the actual hidden lifted response.
It is retained as metric information.

The inverse in (16) is the full hidden inverse. With
D0=widehat D tensor I+I tensor h_env and q_z=v/(d-z)<1,
YC39's series sums every boundary insertion staying hidden. Replacing
that inverse by terms of orders 0,...,m has two-sided return error
bounded by

\[
 \frac{q_z^{m+1}}{(d-z)(1-q_z)}\,W_E(\eta).                 \tag{20}
\]

Thus the complete hidden series now applies to excited retained
endpoints and to the internal/boundary mixed source as well.

R_E is a compression of the complete retained-space Schur return.
It has not eliminated retained energies above E. Transitions through
those retained states and subsequent reentry into the hidden space
are still present in the full model, not discarded from it.

## 6. Quantitative improvement at regular block size

Assume n<=C_partial b^(2/3), with g, c, beta, rho, gamma_* and
C_partial fixed. For a fixed spectral window |z|<=z_*, use YC39's
conditions

\[
 b\ge8,\qquad b\ge4z_*/c,\qquad
 b\ge(4gC_\partial/c)^3.                                  \tag{21}
\]

They give delta_z>=cb/2 and q_z<=1/3. From
C_b(E)<=2(nE+rho gamma_*)/(cb), equation (17) with eta=1 yields

\[
 \boxed{\begin{aligned}
 \|R_E(z)\|\le{}&
 \frac{4\beta^2E^2}{c}\,b^{-1}\\
 &+\frac{8g^2C_\partial^2E}{c^2}\,b^{-2/3}
 +\frac{8g^2\rho\gamma_*C_\partial}{c^2}\,b^{-4/3}.
 \end{aligned}}                                           \tag{22}
\]

At each fixed E>0 this is O(b^(-2/3)); the response metric is
O(b^(-5/3)). The constants here no longer contain YC38's unevaluated
quasi-local covariance constant. They use the local derivative budget
instead. This does not by itself certify a practical block size.

There is also a growing-energy window: if E=E_b=o(sqrt(b)), every
term in (22) tends to zero. The statement concerns the displayed
return on that increasing input space, not a claim that its retained
dynamics is uniformly contractive under joining blocks.

On the actual vacuum the internal source vanishes. Equations (15)-(16)
give the stronger estimates

\[
 \boxed{\|R_0(z)\|
       \le\frac{2g^2\rho\gamma_*C_\partial}{c^2}\,b^{-4/3},}
 \qquad
 \boxed{\|R_0'(z)\|
       \le\frac{4g^2\rho\gamma_*C_\partial}{c^3}\,b^{-7/3}.}
                                                               \tag{23}
\]

For this same graph cut and the smooth Wilson port subclass, (23)
improves YC39's O(b^(-1/3)) vacuum-return bound. YC38/39's earlier
Gram result remains useful for nonsmooth bounded local ports and is
not invalidated.

## 7. Exact reference control: one step around a Fourier circle

As an algebraic control of (2), take h=-d^2/dtheta^2 on the unit
circle, with its constant ground and Fourier basis e_k. Retain
|k|<=r, r>=1, so Y=0, J=P, L=Q and d=(r+1)^2. For the single smooth
port f=e^(i theta),

\[
 \|f\|_\infty=1,\quad\gamma_*=\rho=n=1,\quad
 QfP e_r=e_{r+1}.
\]

At input energy E=r^2, the bound in (2) is exactly

\[
 C_b(E)\ \hbox{with general floor }d
       =\frac{(\sqrt E+1)^2}{d}=1,
\]

and the source norm squared is one: the estimate is attained. If
the input Fourier support stays below r, the source is exactly zero.
This control concerns the source inequality, which does not require
the single port to be self-adjoint. A self-adjoint complete interaction
can include its adjoint partner before using the return formulas.
It is a written circle-model identity, not a Yang--Mills numerical run.

## 8. Zero reading, pairing and units

At zero internal mixed coupling, h=H0, Y=0 and the graph cut becomes
the original reference count cut. The retained energy cap has count
at most floor(2E), since H0>=N_count/2 and the two commute. If each
port acts on at most k_0 reference factors and

\[
 \lfloor b/8\rfloor\ge\lfloor2E\rfloor+k_0,
\]

then a single boundary insertion cannot enter Q. Hence B_E=0, the
internal source is zero, and the exact R_E is zero even with the
complete hidden boundary inverse present. Equations (2) and (22)
are upper bounds and need not vanish at this specialization; the
exact constructed source does. Its coupling derivatives remain.

At zero boundary strength B_E=0, but the internal source
Y h_ret Pi_E generally remains for excited inputs. The corresponding
internal return must not be set to zero. It vanishes on the vacuum.
Histories first visiting other retained states may reach the hidden
channel later, even when the single low-energy reference source is zero.

Under a unitary frame change transport h, J, L, the source operators,
the ground-form map and its derivatives together. The exact form and
Schur return are then intertwined. For a source-label change transform
coefficients dually. The actual V and its derivative quadratic form
are preserved; the convenient n,rho,gamma_* bounds belong to the
chosen normalized local presentation and carry its metric cost.

If physical energies equal e times lattice electric energies, then
h, h_ret, E, d, z and gamma_p scale by e, while A_p and V scale by e.
The graph maps and C_b(E) are unchanged under this simultaneous
conversion. R_phys(ez)=e R(z); its first energy derivative is
dimensionless. Scalar ground-energy offsets remain recorded.

The represented operator family, domains, gauge action, pairing,
source cuts and ordered composition are the carrier. Coupling,
block size, count cap, energy cap and spectral variable have distinct
roles; this calculation does not identify any one of them with
primitive lambda-space or elapsed time.

## 9. The next obstruction is now a closure question

YC39's excited-input source obligation is supplied for the actual
smooth Wilson ports, and the complete graph-hidden excursion has a
vanishing bound on fixed and certain growing energy windows. There
is also the all-input relative form estimate (13).

For spatial iteration the new retained interaction includes
J* f_p J and energy-dependent returned operators. These are not
automatically scalar Wilson multiplication functions with the same
derivative incidence. Their energy-form representation and source
budgets must be transported or reconstructed. One must also control
transitions among retained energy windows and repeated retained/hidden
alternation. Equation (22) alone does not supply either closure.

The next candidate is therefore an energy-weighted class of retained
interactions and returned metrics that is stable under those joins.
Preserve the exact source form in (12), rather than immediately
reducing every generated interaction to its full operator norm.
Existing lattice floors remain inherited inputs; continuum scaling
and a continuum mass gap remain open.

**Verification performed:** written completed-square, ground-form,
operator Cauchy--Schwarz, product-rule, source and return derivations;
the exact Fourier-circle control; zero-specialization and dimensional
checks. No numerical crossover, executable certificate or predecessor
test suite was run. Source revisions and hashes are recorded in
[SOURCES.json](SOURCES.json).

The ground-state form is inherited explicitly from
[YC22](../yc22/YC22_GAUGE_PROTECTED_CHANNELS.md). Standard analytic
background includes Phan Thanh Nam's
[Mathematical Quantum Mechanics II, section 2.4](https://www.math.lmu.de/~nam/LectureNotesMQM2020.pdf)
for localization energy identities, and Dusson--Sigal--Stamm's
[Feshbach--Schur map paper](https://arxiv.org/abs/2105.02058)
for Schur reduction. The present derivation is written above; these
references are not evidence of a published proof of this YM stage,
and no external priority claim is made.
