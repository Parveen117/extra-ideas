# YC15 — the bridge source returns into specific cube channels

9 October 2026. Continues YC14 at `a7ae437`. Frozen packets are unchanged.
Unit-S³ link Laplacian, normalized SU(2) Haar measure, and the actual correlated
cube reference of YC12/YC14 are retained throughout.

**Result.** For internal cube couplings 0<=theta_c<=1 and external plaquette
couplings 0<=eta_p<=eta, on every even periodic spatial lattice with sides>=4,

\[
\boxed{\eta\le\frac1{4320}\quad\Longrightarrow\quad
\Delta_{\rm full}\ge\frac{1279511}{4672512}-\frac{1024}{5}\eta
\ge\frac{979629}{4326400}>\frac9{40}.}                    \tag{1}
\]

The interface window grows by 32/27 relative to YC14's 1/5120, while the
endpoint floor improves from 21/100 to more than 9/40. All factor harmonics,
boundary charges and nonlinear cluster orders are included. The physical
gauge-invariant, centre-even vacuum sector inherits the same bound. This
remains a correlated, nonuniform lattice theorem, not an enlarged isotropic
window, spatial renormalization map or continuum mass-gap construction.

The main structural result is more specific than the numerical gain:
the complete returning source has a small energy variance, and its first
nonconstant bridge-free return feeds the two faces incident on each participating
cube edge, with exact coefficient **1/22464**. Single-cube return channels
really occur; the first source's four-factor support does not persist.

## 1. Carrier, reference and the return being computed

YC12 partitions all links into disjoint twelve-edge cubes and individual
bridges. A cube Hamiltonian is

\[
H_c=H_{0,c}+\theta_c(6-\mathcal S_c),\qquad
H_{0,c}=-\sum_{e\in c}\Delta_e,\qquad \mathcal S_c=\sum_{f\subset c}W_f.
\]

Its actual normalized positive ground is Omega_c, with energy E_c. Set
h_c=H_c−E_c>=0, Omega=product_c Omega_c times the constant bridge readings,
and H_ref=sum_c h_c+sum_bridge(−Delta_e). The pairing is the full Haar L²
pairing, or equivalently Omega² dmu after the exact ground-state transform.
No trial comparison function replaces Omega_c.

For an external face p, the source r_p=Wp Omega has norm²=1/4 and
<r_p,H_ref r_p>=3. Its normalized mean energy is exactly 12. A two-bridge
face excites both bridges and both incident cubes, so its complete reducing
source sector has floor d=13/2. Four-bridge faces have exact source energy12.
These are YC14's full-carrier facts.

We compute two related objects:

1. The full reference inverse D0^−1 acting on r_p, where D0 is H_ref on the
   complement of Omega. Its norm controls the initial cluster creator.
2. The part of Wp D0^−1 Wp Omega returned to the bridge vacuum. It acts on
   cube readings even though the first source excited four factors.

The spectral intertwiner remains the complete creation similarity of YC9/YC10.
The information carried back is a vector/operator with a declared inverse;
it is not identified with the thermodynamic scalar ST or with a physical
observer-induced energy. YC14's local-zero-cut interpretation is preserved.

## 2. Two bridge integrals retain a product of two cube readings

Write the internal edge quaternions as a,b and the two bridge quaternions
as u,v, with Wp=Re(a u b^−1 v^−1). With a second pair a',b', normalized
Haar integration gives the exact kernel

\[
\boxed{\int W_p(a,b,u,v)W_p(a',b',u,v)\,du\,dv
        =\tfrac14(a\cdot a')(b\cdot b').}                 \tag{2}
\]

One v integration gives Re(a u b^−1 b' u^−1 a'^−1)/4. Averaging conjugation
by u replaces b^−1 b' by its scalar part, proving (2). The code also expands
all quaternion coefficients and contracts both second-moment tensors exactly.

This kernel is already a return with memory: it records both internal links,
not just the norm1/4 at a=a', b=b'. Differentiating it in unit-quaternion
tangent directions A at a and B at b yields

\[
\int (\partial_A W_p)^2=|A|^2/4,\quad
\int (\partial_B W_p)^2=|B|^2/4,\quad
\int (\partial_A W_p)(\partial_B W_p)=0,                  \tag{3}
\]

because a·A=b·B=0. These identities require no assumption of independent
links inside either cube; bridge Haar variables alone are integrated.

## 3. Exact source variance equals the two internal edge kinetic energies

Let s_p=(H_ref−12)r_p. Acting on the product Wp Omega and using
H_ref Omega=0 and H0 Wp=12Wp gives

\[
s_p=-2\sum_{e\in p\,\mathrm{internal}}
             \nabla_e W_p\cdot\nabla_e\Omega.             \tag{4}
\]

There are two internal edges, one in each cube. Apply (3) conditionally and
then integrate the cube density. For
T_(c,e)=integral |grad_e Omega_c|² dmu_c, this proves

\[
\boxed{\|s_p\|^2=T_{(A,e_A)}+T_{(B,e_B)},\qquad
       \langle r_p,s_p\rangle=0.}                        \tag{5}
\]

The second identity is also YC14's mean-energy sum rule. Thus changes in
the actual cube ground appear as an exact variance of the boundary source.
It is the first moment together with this variance that controls the inverse.

Here is a full-window bound for the kinetic terms, not a perturbative
ground approximation. Work on the gauge-invariant cube carrier for its
ground only. Let Omega_c=a 1+q, q perpendicular to1, and
epsilon=E_c−6theta<=0 by the constant trial. H0 on this complement has
floor12, ||S||<=6, and ||H0^−1/2 S1||²=(3/2)/12=1/8. The ground equation
projected onto q and Cauchy–Schwarz in the H0 form give

\[
(1-\theta/2)T_c
\le\theta a\langle q,\mathcal S1\rangle
\le\theta |a|\sqrt{T_c/8},
\quad
T_c=\langle\Omega_c,H_0\Omega_c\rangle
\le\frac{\theta^2}{8(1-\theta/2)^2}.                     \tag{6}
\]

The scalar a has |a|<=1. Cube automorphisms and uniqueness of the actual
ground make the twelve edge kinetic expectations equal. Hence, for
0<=theta_A,theta_B<=1,

\[
\boxed{\|s_p\|^2\le v(\theta_A,\theta_B)
 :=\sum_{c=A,B}\frac{\theta_c^2}{96(1-\theta_c/2)^2}
 \le\frac1{12}.}                                        \tag{7}
\]

Only this ground-state estimate uses the gauge free floor12. The boundary
source and the eventual excitation theorem retain **all charged states**;
no gauge-sector gap is substituted for a charged-cube gap.

## 4. The entire inverse, with its quadratic residual retained

For x>0, the exact scalar identities at mean mu are

\[
\frac1x=\frac1\mu-\frac{x-\mu}{\mu^2}
                    +\frac{(x-\mu)^2}{\mu^2x},
\]
\[
\frac1{x^2}=\frac1{\mu^2}-\frac{2(x-\mu)}{\mu^3}
             +(x-\mu)^2\frac{2x+\mu}{\mu^3x^2}.          \tag{8}
\]

Use their spectral calculus on the complete source sector with mu=12,
norm²(r)=1/4, <r,s>=0 and x>=d=13/2. The inverse remainder coefficient
is positive and decreasing. Consequently

\[
\langle r,D_0^{-1}r\rangle
=\frac1{48}+\frac1{144}\langle s,D_0^{-1}s\rangle,
\quad
\boxed{\frac1{48}\le\langle r,D_0^{-1}r\rangle
\le\frac1{48}+\frac{v}{144d}\le\frac{235}{11232}.}        \tag{9}
\]

The final upper/lower ratio is 235/234. This improves YC14's upper1/26
without assuming the inverse stays in the one-dimensional source span.
Similarly,

\[
\boxed{\|D_0^{-1}r\|^2
\le\frac1{576}+v\left(\frac1{864d}+\frac1{144d^2}\right)
\le\frac{773}{438048}<\left(\frac{27}{640}\right)^2.}      \tag{10}
\]

Four-bridge sources still have inverse-vector norm exactly1/24. The source
variance vanishes at zero internal coupling and both inequalities recover
the exact free result. No expansion in the interface coupling is used in
(7)–(10). They hold over the whole declared internal-coupling window.

## 5. A full return operator and its nonconstant cube part

First restrict to the four active factors of a two-bridge face: two full
cubes and two full bridges. Let P_b project only the bridges onto constants,
leaving the complete two-cube carrier. Set H_p=h_A+h_B+H_bridges and Q_b=1−P_b.
The complete local reference return at energy zero is

\[
\mathcal R_p=P_bW_pQ_b(Q_bH_pQ_b)^{-1}Q_bW_pP_b.
\]

The bridge complement has a strictly positive floor. For internal-edge
quaternion coordinates a_alpha, define the positive operator

\[
M_{c,e}(t)=\sum_{\alpha=0}^3
          a_\alpha e^{-t h_c}a_\alpha.
\]

Equation (2), the exact Laplace representation of the inverse and the two
bridge energy6 give

\[
\boxed{\mathcal R_p=\frac14\int_0^\infty
       e^{-6t}M_{A,e_A}(t)\otimes M_{B,e_B}(t)\,dt.}       \tag{11}
\]

All cube eigenstates enter the heat operators. Since h_c>=0 and sum a_alpha²=1,
0<=M_c(t)<=I, so 0<=R_p<=I/24 on the full two-cube carrier. It is generally
an integral operator, not multiplication by a scalar potential.

On the full lattice, spectator cube energies must also be kept in the inverse:
the corresponding expression includes exp(−t H_rest). They reduce to1 only
when acting on their reference grounds, as in the boundary-source return
computed here. Equation (11) alone is not a local effective Hamiltonian for
arbitrary excited spectators or a spatial RG map.

For the actual reference Omega, P_b Wp² Omega=Omega/4 and
||P_b Wp||=1/2 (the adjoint square is multiplication by the bridge average
Wp²=1/4). From D0^−1 r=r/12−D0^−1 s/12,

\[
\boxed{\left\|\mathcal R_p\Omega-\frac{\Omega}{48}\right\|^2
\le\frac{v}{576d^2}\le\frac1{292032}.}                   \tag{12}
\]

The same bound controls its projection perpendicular to Omega, hence all
nonconstant cube readings together. The scalar part is the narrower (9).
For distinct p,q, P_b Wq D0^−1 Wp Omega=0: their different bridge signatures
leave at least one unmatched degree-one bridge, as established in YC14.
Thus at second order in interfaces only a self-paired face can return to
the complete bridge vacuum. Later interacting returns are not set to zero.

## 6. Which cube readings return first: an exact derivative at the free point

The actual normalized cube ground has Omega_theta=1+theta S/12+O(theta²),
since its centered energy derivative vanishes and H0 S=12S. This is a
derivative of the exact ground branch, not the trial ground used in YC8.
For a face f incident on the internal edge e, the product a_alpha Wf has
degree0 or2 on e and degree1 on the other three face edges. Its free
energies are9 and17, with contractions of weights1/4 and3/4, respectively.
For a nonincident face its free energy is15.

Differentiate the exact heat operator and its ground factor. If
F_(c,e)(t)=Omega_c^−1 M_(c,e)(t) Omega_c, the result is

\[
F_{c,e}(t)=e^{-3t}
 +\theta_c f(t)\sum_{f\ni e}W_f+O(\theta_c^2),
\]
\[
\boxed{f(t)=\frac{e^{-3t}}{84}-\frac{e^{-9t}}{48}
                                  +\frac{e^{-17t}}{112}.} \tag{13}
\]

Nonincident faces cancel against the reference-ground change. For completeness,
the contribution of an incident face before simplification is
sum_(E=9,17) w_E [exp(−Et)/12+(exp(−3t)−exp(−Et))/(E−3)]−exp(−3t)/12.
This includes the Duhamel derivative and the derivative of Omega itself.
Dropping either would give a different, incorrect coefficient.

The response has f(0)=f'(0)=0, f''(0)=1 and f(t)>0 for every t>0:
with z=exp(−2t),
4−7z³+3z⁷=(1−z)²(3z⁵+6z⁴+9z³+12z²+8z+4), and
f(t)=exp(−3t)(4−7z³+3z⁷)/336.

Insert (13) into (11). The returned reference reading is

\[
\boxed{\Omega^{-1}\mathcal R_p\Omega
=\frac1{48}+\frac1{22464}
 \left(\theta_A\sum_{f\subset A,\ f\ni e_A}W_f
      +\theta_B\sum_{f\subset B,\ f\ni e_B}W_f\right)
 +O((|\theta_A|+|\theta_B|)^2).}                         \tag{14}
\]

Here Omega^−1 R_p Omega denotes the function obtained by applying the operator
to Omega and dividing by that positive function, not the entire transformed
operator. The coefficient is exactly
(1/4) integral exp(−9t) f(t) dt=1/22464.
An independent stationary-resolvent derivation uses YC12's adjacent product
energies18 and26 and gives
(1/6)[1/(16·18)+3/(16·26)]−1/576=1/22464. A nonincident face instead gives
(1/6)/(4·24)−1/576=0.

Each sum has two faces. To first order in the internal couplings these are
single-cube excitations; there is no two-cube-excited component at that order.
The scalar susceptibility (9) changes only at second order, while this
nonconstant returned reading already changes at first order. A scalar energy
record can therefore miss an active spatial return channel. This is an exact
first-visible-order distinction on the same source, not a comparison of
unrelated quantities carrying the same name.
The Schur effective restoring operator subtracts eta_p² times this return,
so the displayed positive returned face reading carries a negative sign there.
It is not a proved replacement of the entire effective operator by a Wilson
potential with a new coupling.

Equations (13)–(14) specify the exact derivative at theta_A=theta_B=0.
No numerical quadratic Taylor remainder is supplied for that expansion at
theta=1. The full-window bound (12) and the all-order gap proof below remain
valid independently of a linear truncation. Analyticity near zero follows
from the isolated simple cube ground and bounded potential perturbation.

## 7. Full nonlinear return: use the source variance and the unused cube reserve

For a two-bridge interface, (10) bounds the inverse-vector norm by27/640;
its first creator has four factors. A cube meets24 such faces, while each
bridge meets two of this type and two four-bridge faces. Therefore

\[
\boxed{\|\mathcal T(0)\|_*\le\frac{81}{20}\eta,}
\qquad \text{bridge-root cost}\le\frac{161}{240}\eta.     \tag{15}
\]

This improves YC14's96eta/13. Keep its exact source-centred map estimate
||T(c)||_*<=||T(0)||_*+q r, including every subsequent creator support.
Equation (14) explicitly demonstrates why the nonlinear minimum support
must remain **s=1**, even though the first creator has four factors.

Also retain YC12's entire certified charged-cube reserve instead of rounding
it down to1/4:

\[
g_*:=\frac{1279511}{4672512}>1/4.
\]

This is a full-space floor for all theta_c in[0,1]; bridges have the larger
gap3. In YC10's complete cluster criterion take k=4, M=4, a=8,
beta<=24eta, r=1/128 and exp(8r)<16/15. Then

\[
q\le\frac{5184}{5g_*}\eta,\qquad
b\le\frac{1024}{5g_*}\eta.                              \tag{16}
\]

At eta=1/4320 the exact outward budget is

| Quantity | Bound |
|---|---:|
| First creator norm | 3/3200 |
| Full nonlinear contraction q | 28035072/31987775 <1 |
| Ball reserve r−seed−qr | 11417/409443520 >0 |
| Excitation relative return b | 5537792/31987775 <1 |
| Actual full gap g_*(1−b) | 979629/4326400 >9/40 |

These left-hand costs are linear in eta, so the whole interval is covered.
The domain, triangular similarity and relative-Neumann inverse argument of
YC9/YC10/YC14 is unchanged: at each finite volume all subset creators are
included, each component belongs to Dom H_I, the similarity preserves the
complete domain, and no nonzero spectrum lies below g_*(1−b). The constructed
eigenvector is therefore the unique true ground. This proves (1) uniformly
over the allowed finite volumes, not just on enumerated sample boxes.

Both improvements matter for this endpoint. Replacing (15) by YC14's old
seed fails the ball inequality, even with g_*. Replacing g_* by1/4 also
fails, even with the new seed. The heat-coefficient expansion (14) is not
used as an uncontrolled approximation in this gap proof.

## 8. Evidence and next gate

There are 31 exact checks and 15 new tests; the 29 predecessor tests in YC12
and YC14 also pass. The executable certificate uses exact quaternion/Haar contractions, spectral
residual identities, independent heat and resolvent coefficients, rational
whole-window budgets and frozen-input SHA256 pins. Tests include signed Haar
cubature, source variance factors, an inverse that leaves its source span,
and failure controls for omitted reference-ground and nonlinear support terms.
The analytic full-carrier proof is written above; numerical checks are not a
formal proof-assistant verification.

Run `python physics/yc15/yc15_boundary_return_channels.py --check` and
`python -m unittest discover -s physics/yc15 -p 'test_*.py'` from the repository.

The classical tools are the spectral theorem, ground-state transform,
Duhamel differentiation, Schur complement and the Kirkwood–Thomas full-return
construction already derived and credited in YC9/YC10. These are applied here
to the actual correlated cube source with new complete residual bounds.

Next target: feed the now-resolved single-cube return into a revised reference
block and control the remaining connected two-cube/interacting return. Keep
spectator energy, source dependence and the induced operator rather than
equating a ground-vector reading with a local scalar potential. This is the
next possible scale step; a continuum interpretation still needs its proof.
