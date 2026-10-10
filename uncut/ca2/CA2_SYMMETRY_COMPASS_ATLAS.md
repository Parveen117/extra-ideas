# CA2 — Literal transport of the TVSP compass to another symmetry

Monty Dabas. 9 October 2026. Exact symbolic algebra; Python 3.12.

The framework's four-place diagram is a derivative template. For another
valid symmetry sector, do not call its properties heat capacities and do not
invent a fresh analogy. Replace the four typed readings in the existing TVSP
formulas, while carrying the potential, held-fixed labels, response metric and
observation map through the same substitution.

For the electromagnetic constitutive sector the owner's substitution is

\[
\boxed{T\mapsto E,\qquad V\mapsto B,\qquad
S\mapsto D,\qquad P\mapsto H.}                              \tag{1}
\]

This stage proves that the old thermodynamic ratio identities then become
exact conditional permittivity and permeability identities.

## 1. The literal four-property transport

The thermodynamic templates retained from CA1 are

\[
\begin{aligned}
C_P&=T\left(\frac{\partial S}{\partial T}\right)_P,&
C_V&=T\left(\frac{\partial S}{\partial T}\right)_V,\\
C_S&=-P\left(\frac{\partial V}{\partial P}\right)_S,&
C_T&=-P\left(\frac{\partial V}{\partial P}\right)_T.
\end{aligned}                                               \tag{2}
\]

Apply (1) position by position. After choosing the positive orientation of
the passive electromagnetic response, the four transported properties are

\[
\begin{aligned}
C_H&=E\left(\frac{\partial D}{\partial E}\right)_H,&
C_B&=E\left(\frac{\partial D}{\partial E}\right)_B,\\
C_D&=H\left(\frac{\partial B}{\partial H}\right)_D,&
C_E&=H\left(\frac{\partial B}{\partial H}\right)_E.
\end{aligned}                                               \tag{3}
\]

The pressure minus sign in (2) belongs to the thermodynamic orientation
\(P=-U_V\). Electromagnetic energy uses \(H=+u_B\). Carrying that orientation
with the label gives the positive definitions in (3). Keeping a common minus
sign on the second pair would change both \(C_D,C_E\), but not their ratio.

The names under the amplitude factors are the conditional differential
constitutive coefficients:

\[
\epsilon_B=\left(\frac{\partial D}{\partial E}\right)_B,\quad
\epsilon_H=\left(\frac{\partial D}{\partial E}\right)_H,\quad
\mu_D=\left(\frac{\partial B}{\partial H}\right)_D,\quad
\mu_E=\left(\frac{\partial B}{\partial H}\right)_E.          \tag{4}
\]

Thus \(C_B=E\epsilon_B,\ C_H=E\epsilon_H,\ C_D=H\mu_D,\)
and \(C_E=H\mu_E\). These are the exact electromagnetic equivalents of the
four old properties. The dimensionless objects are their ratios.

## 2. Exact electromagnetic calculation

Take one scalar or fixed-polarization reciprocal constitutive sector with
energy density

\[
u(D,B)=\frac12(aD^2+2bDB+cB^2),\qquad
E=u_D,\quad H=u_B,                                          \tag{5}
\]

where

\[
M=\frac{\partial(E,H)}{\partial(D,B)}
=\begin{pmatrix}a&b\\b&c\end{pmatrix}>0,\qquad
\Delta=ac-b^2>0.                                            \tag{6}
\]

Holding the second source fixed gives

\[
\epsilon_B=\frac1a,\qquad \mu_D=\frac1c.                    \tag{7}
\]

Holding the conjugate response fixed performs the same Schur complement that
turns \(C_V\) into \(C_P\):

\[
\epsilon_H=\frac{c}{\Delta},\qquad
\mu_E=\frac{a}{\Delta}.                                     \tag{8}
\]

Consequently the old TVSP law is transported literally:

\[
\boxed{
\frac{C_B}{C_H}
=\frac{\epsilon_B}{\epsilon_H}
=\frac{C_D}{C_E}
=\frac{\mu_D}{\mu_E}
=\frac{\Delta}{ac}
=1-\frac{b^2}{ac}
=\chi_{\rm EM}.}                                            \tag{9}
\]

The inverse form is

\[
\frac{C_H}{C_B}
=\frac{\epsilon_H}{\epsilon_B}
=\frac{C_E}{C_D}
=\frac{\mu_E}{\mu_D}
=\chi_{\rm EM}^{-1}.                                        \tag{10}
\]

For the exact witness \(a=4,b=1,c=3\),

\[
\chi_{\rm EM}=\frac{11}{12},\qquad
1-\chi_{\rm EM}=\frac1{12}.                                 \tag{11}
\]

For an uncoupled vacuum constitutive block \(b=0\), so
\(\chi_{\rm EM}=1\). After axis/impedance normalization the block is the
identity. Reciprocal magnetoelectric coupling is precisely the lost share
\(b^2/(ac)\) on this carrier.

These \(\epsilon\) and \(\mu\) are static or quasistatic conditional
differential constitutive coefficients. If “transport coefficient” means a
frequency-dependent dissipative coefficient, the response is generally
complex or nonreciprocal and section 6 applies.

## 3. Why the same algebra keeps appearing

For any two-source positive reciprocal potential \(\Phi(x,y)\), put

\[
p=\Phi_x,\qquad q=\Phi_y,\qquad
\nabla^2\Phi=\begin{pmatrix}A&B\\B&C\end{pmatrix}>0.         \tag{12}
\]

Its pre-observation four-place diagram is

\[
\boxed{p-y-x-q}.                                             \tag{13}
\]

The four transported scale responses are

\[
\begin{aligned}
\mathcal C_{x|y}&=p_0(\partial x/\partial p)_y=p_0/A,&
\mathcal C_{x|q}&=p_0(\partial x/\partial p)_q=p_0C/\Delta,\\
\mathcal C_{y|x}&=q_0(\partial y/\partial q)_x=q_0/C,&
\mathcal C_{y|p}&=q_0(\partial y/\partial q)_p=q_0A/\Delta.
\end{aligned}                                               \tag{14}
\]

Therefore

\[
\frac{\mathcal C_{x|y}}{\mathcal C_{x|q}}
=\frac{\mathcal C_{y|x}}{\mathcal C_{y|p}}
=\frac{\det\nabla^2\Phi}{\Phi_{xx}\Phi_{yy}}.                \tag{15}
\]

The repetition is forced by one operation: fixing the conjugate response
eliminates the other source and returns the Schur complement. The same
two-by-two determinant ratio must therefore appear in every valid literal
replacement.

## 4. Atlas of replacements

| Sector | Literal four-place diagram | Pure ratio |
|---|---|---|
| thermodynamics | \(T-V-S-P\) | \(C_V/C_P=C_S/C_T\) |
| electromagnetic constitutive channel | \(E-B-D-H\) | \(\epsilon_B/\epsilon_H=\mu_D/\mu_E\) |
| two elastic channels | \(\sigma_1-\epsilon_2-\epsilon_1-\sigma_2\) | conditional compliance ratio |
| two source observables | \(\langle O_1\rangle-j_2-j_1-\langle O_2\rangle\) | \(1-\operatorname{Corr}(O_1,O_2)^2\) |

For elasticity, replace
\((T,V,S,P)\) by
\((\sigma_1,\epsilon_2,\epsilon_1,\sigma_2)\) in the same templates.
For a source-deformed partition function

\[
W(j_1,j_2)=\log Z,\qquad
\partial_a\partial_b W=\operatorname{Cov}(O_a,O_b),         \tag{16}
\]

the same ratio is

\[
\chi=1-\frac{\operatorname{Cov}(O_1,O_2)^2}
{\operatorname{Var}(O_1)\operatorname{Var}(O_2)}.           \tag{17}
\]

For Yang--Mills, choose gauge-invariant \(O_1,O_2\), then make the literal
replacement

\[
(T,V,S,P)\mapsto
(\langle O_1\rangle,j_2,j_1,\langle O_2\rangle).             \tag{18}
\]

The source pair remains part of the physical statement. Gauge symmetry says
which sources are admissible; it does not choose one pair by itself.

## 5. Symmetry restriction and observation

If a group acts on the two source channels by \(R_g\), at a symmetry-fixed
background

\[
R_g^T H R_g=H.                                               \tag{19}
\]

A full quarter-turn symmetry on a real symmetric two-channel block forces
\(H=aI\), hence \(\chi=1\). A nontrivial \(\chi<1\) means that the chosen
background or channel split retains coupling instead of the full isotropic
symmetry.

When more than two admissible channels exist, the source pair must still be
declared. The positive witness

\[
G=\begin{pmatrix}4&1&0\\1&3&1\\0&1&2\end{pmatrix}            \tag{20}
\]

has \(\chi_{12}=11/12\) and \(\chi_{23}=5/6\). Literal
substitution works for each pair, but “the symmetry” has not selected which
one an experiment reads.

The owner defines TVSP before observation and VTSP after observation.
Accordingly every analogue uses

\[
(p,y,x,q)\longmapsto(y,p,x,q),                              \tag{21}
\]

with the source, metric, held-fixed conditions and response frame pushed
through the chart. For electromagnetism this is

\[
E-B-D-H\longmapsto B-E-D-H.                                 \tag{22}
\]

## 6. Reciprocity boundary

The exact single-ratio formulas above assume a positive reciprocal Hessian.
If the mixed responses differ,

\[
L=\begin{pmatrix}A&B_1\\B_2&C\end{pmatrix},\qquad B_1\ne B_2,
\]

then

\[
\chi=1-\beta^2+\omega^2,\qquad
\beta=\frac{B_1+B_2}{2\sqrt{AC}},\quad
\omega=\frac{B_2-B_1}{2\sqrt{AC}}.                         \tag{23}
\]

The turn share \(\omega\) is extra retained memory. This is the correct
extension for nonreciprocal media, Hall-type response, dissipation or
path-dependent cuts. The literal replacement remains useful, but one
dimensionless ratio no longer classifies the response.

The propagating vacuum-wave carrier of EM1 is also separate. There Maxwell's
equations constrain \(E,B\); they are not independent constitutive sources.
EM1 reads the field invariant and energy flow. CA2 reads the constitutive
Hessian with \(D,B\) independently varied. A further theorem may connect
those carriers; matching their letters is insufficient.

## 7. Claim boundary

\[
\begin{array}{ll}
T\to E,\ V\to B,\ S\to D,\ P\to H
\text{ transports the four derivative templates} & \textbf{PROVED},\\
\epsilon_B/\epsilon_H=\mu_D/\mu_E=\chi_{\rm EM}
\text{ for a positive reciprocal block} & \textbf{PROVED},\\
\text{every positive reciprocal two-source potential obeys (15)}
& \textbf{PROVED},\\
\text{one symmetry name selects a unique source pair}
& \textbf{FALSE; exact witness},\\
\text{one }\chi\text{ classifies full tensor electromagnetism}
& \textbf{NOT CLAIMED},\\
\text{a dissipative/nonreciprocal response needs no turn number}
& \textbf{FALSE}.
\end{array}
\]

Reproduce from the CA2 directory with
python ca2_symmetry_compass_atlas.py and
python -m unittest -v test_ca2.

The certificate has nine exact checks and eight focused tests. It pins CA1,
UP3, UP4 and EM1 without changing those packets.
