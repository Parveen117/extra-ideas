# YC16 — absorb the complete one-cube return and keep the connected operator

9 October 2026. Continues YC15 at `8a976df`. Frozen packets are unchanged.
The user authorized the next step: include the returned contribution in the
reference cube, then bound the remaining interaction on the actual operator.

**Result.** The complete local second-order return has a unique decomposition
relative to the actual cube grounds into scalar, two one-cube operators and
a connected two-cube operator. The one-cube operators are absorbed into a
new gapped reference. An exact add/subtract identity retains all compensating
terms in the original Hamiltonian. The full cluster proof on this new reference
still gives **gap>9/40 for eta<=1/4320**, theta_c<=1, uniformly on even tori
with sides>=4. This preserves YC15's window; it does not enlarge it.

The connected **vacuum source** has norm<=1/8658 on the original reference,
and <1/8657 after transporting it to the new reference. It vanishes if
either incident cube is free, and starts at mixed internal order theta_A theta_B
with coefficient19/86261760. Its **operator** is retained on all excited states:
it is nonzero even at the free point. Thus absorption is not a replacement
of the whole return by a scalar Wilson-coupling shift.

## 1. Declared source, metric and full carrier

Use YC12/YC15's full cube Hilbert spaces, normalized Haar measure and
unit-S³ electric normalization. Let h_c=H_c−E_c>=0 have its actual normalized
ground Omega_c. For 0<=theta_c<=1 the charged and uncharged excited states
all have the certified floor

\[
g_* =1279511/4672512.
\]

The original reference is the tensor product of these cubes and constant
bridge states. For a two-bridge external face p, YC15 supplies the bounded
operator on its complete two-cube carrier

\[
R_p=\frac14\int_0^\infty e^{-6t}M_A(t)\otimes M_B(t)\,dt,
\quad M_c(t)=\sum_{\alpha=0}^3 u_{e,\alpha}e^{-th_c}u_{e,\alpha}.
                                                               \tag{1}
\]

Its source is Wp Omega, not an arbitrary vector confined to a finite harmonic
span. For omega_c(X)=<Omega_c,X Omega_c>, put

\[
m_c(t)=\omega_c(M_c(t)),\qquad N_c(t)=M_c(t)-m_c(t)I.
\]

These are bounded-operator conditional expectations with the actual ground
pairing. Their use here is an explicit native retained/returned split; there
is no identification of R_p with ST or of theta_c with observation lambda.

The local formula (1) is exact for the four active factors. A full-lattice
resolvent on excited spectators includes their heat factor as YC15 explains.
Only the reference-source coefficient has all spectators in their grounds.
The all-order proof below works directly with the original Hamiltonian and
does not discard spectator energies by replacing it with (1).

## 2. Canonical one-cube and connected pieces

There is exactly one decomposition with zero reference means and zero
connected conditional expectations:

\[
\boxed{R_p=\alpha_p I+A_{p,A}\otimes I+I\otimes A_{p,B}+C_p,}
                                                               \tag{2}
\]
\[
\alpha_p=(\omega_A\otimes\omega_B)(R_p),\quad
A_{p,A}=(\mathrm{id}\otimes\omega_B)(R_p)-\alpha_p I,
\quad A_{p,B}=(\omega_A\otimes\mathrm{id})(R_p)-\alpha_p I.
\]

The definition of C_p is then forced by (2). In particular
omega_A(A_p,A)=omega_B(A_p,B)=0 and
(id tensor omega_B)(C_p)=(omega_A tensor id)(C_p)=0.
Applying these expectations to any proposed split proves uniqueness.

Equation (1) gives the full operator formulas

\[
\alpha_p=\frac14\int e^{-6t}m_A(t)m_B(t)dt,
\quad A_{p,A}=\frac14\int e^{-6t}m_B(t)N_A(t)dt,
\]
\[
\boxed{C_p=\frac14\int_0^\infty e^{-6t}N_A(t)\otimes N_B(t)dt.} \tag{3}
\]

All cube electric modes remain in these integrals. No rank-one replacement
of M_c, A_p,c or C_p has been made.

From 0<=M_c(t)<=I, 0<=m_c(t)<=1, ||N_c(t)||<=1. Therefore

\[
\|C_p\|\le1/24.                                          \tag{4}
\]

This uses the heat-product representation; an arbitrary positive two-factor
matrix need not have connected-part norm bounded by its own norm.
Also 0<=R_p<=I/24 and YC15 gives 1/48<=alpha_p<=a_*:=235/11232.
Each conditional expectation of R_p lies in[0,I/24]. Since alpha_p>=1/48,

\[
\boxed{\|A_{p,c}\|\le\alpha_p\le a_*.}                  \tag{5}
\]

Gauge transformations commute with R_p and its conditional expectations:
the bridge vacuum projections and the other cube's ground are invariant.
Thus these are admissible gauge-invariant cube operators on the full carrier.

## 3. A substantially smaller connected source, not a discarded operator

Let T_(c,e)=integral |grad_e Omega_c|². YC15 proves
T_(c,e)<=theta_c²/[96(1−theta_c/2)²]<=1/24. Define
s_alpha=(h_c−3)u_(e,alpha) Omega_c. The ground product rule gives

\[
\sum_\alpha\|s_\alpha\|^2=4T_{(c,e)}.
\]

Both u_alpha Omega_c and s_alpha lie in a charged sector orthogonal to the
cube ground. Use the conservative full floor g=1/4 on that entire sector.
Duhamel's identity and the row-contraction norm of (u_0,...,u_3) give

\[
\|N_c(t)\Omega_c\|\le2\sqrt{T_{(c,e)}}J_g(t),\quad
J_g(t)=\frac{e^{-gt}-e^{-3t}}{3-g}.                        \tag{6}
\]

Indeed first compare M_c(t)Omega_c with exp(−3t)Omega_c, write the difference
as the integral of exp(−(t-s)h_c)s_alpha exp(−3s), and then project off the
ground. This cannot increase the norm. The scalar m_c(t) is a positive
spectral average on the same charged sector, hence m_c(t)<=exp(−gt).

Integrate (6) in (3), retaining both decay rates:

\[
\boxed{\|A_{p,A}\Omega_A\|\le\frac4{481}\sqrt{T_{(A,e_A)}},
\qquad
\|C_p\Omega_A\Omega_B\|\le\frac4{1443}
                \sqrt{T_{(A,e_A)}T_{(B,e_B)}}\le\frac1{8658}.} \tag{7}
\]

The constants follow from
(1/2) integral exp(−(6+g)t)J_g(t)dt=4/481 and
integral exp(−6t)J_g(t)²dt=4/1443 at g=1/4.
The connected source is in q_A H_A tensor q_B H_B, where
q_c=I−|Omega_c><Omega_c|: either partial vacuum projection kills it.
If either internal coupling is zero its ground is constant, T_e=0, and
this connected source vanishes exactly, even with the other cube interacting.

The leading mixed derivative can also be read exactly. YC15's
f(t)=exp(−3t)/84−exp(−9t)/48+exp(−17t)/112 gives

\[
\frac{C_p\Omega}{\Omega}
=\boxed{\frac{19\theta_A\theta_B}{86261760}
     \left(\sum_{f\ni e_A}W_f\right)
     \left(\sum_{f\ni e_B}W_f\right)}
 +O(\theta_A^2\theta_B+\theta_A\theta_B^2).               \tag{8}
\]

Here the quotient denotes the returned reference function, not the complete
operator. The coefficient is (1/4) integral exp(−6t)f(t)²dt. Each sum has two
faces. Equation (8) is an exact derivative at the internal free point; (7)
is the separate full-window bound, which uses no linear or quadratic truncation.

**An actual free-point negative control.** For a gauge-invariant internal
face Wf containing e, free kinetics gives

\[
M_e(t)W_f=(e^{-9t}/4+3e^{-17t}/4)W_f,
\qquad m_e(t)=e^{-3t}.
\]

Consequently at theta_A=theta_B=0,

\[
A_{p,A}W_f=-\frac{19}{1872}W_f,\qquad
C_p(W_f\otimes W_{f'})=\frac{3931}{599040}(W_f\otimes W_{f'}), \tag{9}
\]

while A_p,A Omega_A=C_p Omega=0. Both witnesses are in gauge-invariant
cube sectors. They prove directly why a silent vacuum reading cannot justify
removing an operator from the excitation problem or replacing it by a scalar
potential: a nonzero multiplication operator could not annihilate the constant
free ground everywhere and still act this way on Wf.

## 4. The returned one-cube operator is now in the reference

For each cube define the bounded self-adjoint operator

\[
\boxed{B_c=\sum_{p\ni c}\eta_p^2 A_{p,c},\qquad
        \widetilde h_c=h_c-B_c.}                          \tag{10}
\]

The sum is over its24 two-bridge interface faces, with their actual possibly
nonuniform couplings. Centering ensures <Omega_c,B_c Omega_c>=0. From (5)
and (7), for eta=max eta_p,

\[
\|B_c\|\le\epsilon:=\frac{235}{468}\eta^2,
\quad\|B_c\Omega_c\|\le\frac{24}{481}\eta^2<\frac{\eta^2}{20}=:s.
                                                               \tag{11}
\]

For the source bound use sqrt(T_e)<=1/4 and sum24 faces; no cancellation is
assumed. Let e_c be the bottom of h_c−B_c. Min–max gives its next energy
at least g_*−epsilon, and the old ground trial gives e_c<=0. Thus its
ground is simple and its complete charged gap is at least

\[
\widehat g=g_*-\epsilon.                                 \tag{12}
\]

Choose the normalized new ground hatOmega_c with positive overlap with Omega_c,
and set hat h_c=h_c−B_c−e_c. This reference remains gauge invariant; a simple
ground under the connected SU(2) vertex group has the trivial character.
No pointwise-positivity or differential-operator assumption for this new
bounded-operator perturbation is needed by the cluster argument.

The original-ground complement of h_c−B_c has floor hat g. Its complete
Schur inverse proves

\[
|e_c|\le s^2/\widehat g,\qquad
\|\widehat\Omega_c-\Omega_c\|\le2s/\widehat g.             \tag{13}
\]

Throughout eta<=1/4320, hat g>27/100. Hence the convenient complete envelopes
are |e_c|<=eta⁴/108 and ||hatOmega_c−Omega_c||<=10eta²/27. These control
the actual reference ground, not merely its perturbative coefficient.

The exact original Hamiltonian, after its usual scalar Wilson constants,
can now be written

\[
\boxed{H=\sum_c\widehat h_c+H_{\rm bridges}
       +\underbrace{\left[-\sum_p\eta_pW_p+\sum_cB_c\right]}_{\widehat\Phi}
       +\sum_c e_c.}                                    \tag{14}
\]

All compensating one-cube operators are present in hatPhi. This is an exact
identity on the original domain, not a Hamiltonian obtained by silently
subtracting its second-order self-energy. Scalar alpha_p remains in the
return/energy ledger; it is not subtracted again from the original H.

## 5. What this absorption cancels, precisely

Fix the internal theta_c and write eta_p=zeta t_p to count interface order.
At a finite volume the simple isolated ground has its usual analytic
perturbation coefficients near zeta=0. The first interface source excites
four factors, so it has no component on one cube alone. At second order,
all bridges can return to their ground only if the two face sources are
the same face (YC14/YC15's signature theorem). Projecting onto one cube c
therefore gives the exact ground-vector coefficient

\[
q_c\Psi^{(2)}_{\text{one cube}}
  =h_c^{-1}\sum_{p\ni c}t_p^2 A_{p,c}\Omega_c.             \tag{15}
\]

This inverse is on the whole excited cube carrier. Differentiating the new
reference ground of h_c−zeta² sum t_p² A_p,c gives exactly the same coefficient.
Thus in the creation expansion relative to hatOmega, the one-cube ground
creator is removed through interface order two. The connected part C_p
remains as a two-cube source, alongside the scalar shift alpha_p and all
bridge-excited channels. Products of the first four-factor creators cannot
create a one-cube component, so the same cancellation holds for the cluster
logarithm as for the ground vector.

This is the content of absorption. It does not remove all one-cube creators
at higher orders, and it does not identify the complete induced operator
with a Wilson-coupling change. The higher orders are bounded next without
using a truncated expansion or replacing full spectator energies.

The source must also be transported when the reference changes. With
delta=10eta²/27, (4) and the product-vector distance bound imply

\[
\|C_p\widehat\Omega_A\widehat\Omega_B\|
\le\frac1{8658}+\frac\delta{12}
=\frac1{8658}+\frac5{162}\eta^2<\frac1{8657}.              \tag{15a}
\]

It is not generally doubly excited relative to the new ground cuts.
For unit vectors the norm of their expectation-functional difference is
at most twice their vector distance. Since both old conditional expectations
of C_p vanish, recentering this same C_p in the new grounds gives a raw
one-cube conditional norm<=delta/12, and a scalar norm<=delta²/6. After
centering again and summing the24 incident eta_p² C_p contributions,

\[
\|\text{regenerated centred one-cube operator}\|
\le\frac{20}{27}\eta^4+\frac{400}{729}\eta^6,
\qquad
|\text{regenerated scalar per face}|
\le\frac{50}{2187}\eta^6.                                \tag{15b}
\]

These are complete norm bounds for recentering C_p, not claims that every
higher bridge interaction starts at those orders. Such interactions remain
in the full cluster estimate below. If a cube is free, its A_p,c operators
annihilate its constant ground; the corrected ground stays that same vector.
The connected source therefore still vanishes exactly when either cube is free.

## 6. Certifying the full remaining interaction on the new reference

The new two-bridge source still excites four factors: the cube grounds are
gauge invariant, and averaging an open edge against their invariant density
kills it. Its Haar norm remains1/2. To compare its complete inverse with
YC15's, use the common two-bridge fundamental sector before projecting
either cube onto its old or new excited space. Both operators have floor6
there, irrespective of cube excitations.

Write delta=10eta²/27 and d_E=eta⁴/108. The change of the two-cube product
ground is at most2delta; multiplication by Wp on constant bridges has
norm1/2. The two cube operators change by at most2(epsilon+d_E).
The full resolvent identity therefore gives

\[
\|\widehat H^{-1}W_p\widehat\Omega-H^{-1}W_p\Omega\|
\le\frac\delta6+\frac{\epsilon+d_E}{36}
\le\frac{\eta^2}{13}.                                   \tag{16}
\]

No inverse is compressed to the source span. Combine with YC15 to obtain
the new inverse-vector bound27/640+eta²/13. The four-bridge value remains1/24.
For the compensating one-cube source in (14),

\[
\|\widehat h_c^{-1}\widehat q_c B_c\widehat\Omega_c\|
\le(s+\epsilon\delta)/\widehat g\le3\eta^2/16.            \tag{17}
\]

These inequalities use only the full reference gap and bounded perturbation.
The exact rational replay checks the endpoint reserves; every envelope is
monotone on the declared interval.

For the full cluster map, support k=4 and nonlinear minimum support s_min=1
remain unchanged. With r=1/128, M=4, exp(8r)<16/15, the incidence and source
budgets are now

\[
\widehat\beta\le24\eta+\epsilon,\qquad
\|\widehat{\mathcal T}(0)\|_*
\le\frac{81}{20}\eta+\frac3{16}\eta^2+\frac{96}{13}\eta^3=:s_0.
                                                               \tag{18}
\]

The bridge-root source is smaller. Import YC10's complete derivative and
excitation estimates on these actual factors and interactions:

\[
q\le\frac{4\widehat\beta}{\widehat g}\frac{16}{15}
                 [2+8(1+2r)],\quad
b\le\frac{8\widehat\beta}{\widehat g}\frac{16}{15},\quad
\|\widehat{\mathcal T}(c)\|_*\le s_0+qr.                 \tag{19}
\]

At eta=1/4320 the exact bounds satisfy

| Quantity | Outward bound |
|---|---:|
| New charged cube gap | 6218422849/22708408320 >27/100 |
| Complete nonlinear contraction q | 5450044392/6218422849 <877/1000 |
| Excitation relative return b | 3229655936/18655268547 <174/1000 |
| Ball reserve r−s0−qr | >1/40000 |
| Original Hamiltonian gap via new reference | 15425612611/68125224960 >9/40 |

All source and incidence costs increase with eta, while hat g decreases;
the endpoint proves the whole window. The same full-domain creation
similarity and Neumann-inverse argument as YC9/YC10 applies to the new
compact-resolvent factors plus bounded local interactions. It includes
every electric harmonic, boundary charge and cluster size. It proves the
constructed eigenvector is the unique actual ground of the original H.
Gauge and centre symmetry therefore give the physical vacuum conclusion too.

The explicit floor from this reorganized proof is

\[
\boxed{\Delta_{\rm full}\ge
g_*-\frac{1024}{5}\eta-\frac{517}{108}\eta^2>9/40
\quad(0\le\eta\le1/4320).}                              \tag{20}
\]

YC15's slightly stronger numerical floor for the same H is retained.
The quadratic cost in (20) pays for the complete compensating operators;
it is not an assertion that the physical gap changed under a relabelling.
The new result is a justified operator-level reference update and an
independently controlled remaining interaction, with the certified window intact.

## 7. Evidence and next target

Thirty exact checks and fifteen new tests pass, along with the twenty-seven
predecessor tests in YC10 and YC15. The certificate records exact conditional-expectation identities, full-source
Duhamel integrals, free gauge-sector negative controls and rational all-order
budgets. Its code, note, tests and frozen inputs have SHA256 pins. Analytic
operator/domain arguments are written above; finite checks are not a formal
proof-assistant verification.

Run `python physics/yc16/yc16_absorbed_cube_return.py --check` and
`python -m unittest discover -s physics/yc16 -p 'test_*.py'` from the repository.

Next target: use the centred connected channel in a repeated reference update
or a larger genuinely coupled block, with its spectator-energy and full
interaction errors controlled. A single compensated reorganization is not
yet a spatial coarse-graining law or a 4D continuum mass-gap theorem.
