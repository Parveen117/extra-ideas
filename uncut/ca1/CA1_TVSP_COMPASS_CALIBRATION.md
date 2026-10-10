# CA1 — The TVSP compass as a calibration law

Monty Dabas. 9 October 2026. Exact symbolic algebra; Python 3.12.

The framework supplies pure numbers. This stage proves the precise part of the owner's proposed bridge: a TVSP
response compass turns one pure number into a **ratio of physical response coefficients**, and one measured
coefficient then supplies its conjugate coefficient. The pure number fixes response shape. A dimensional anchor
fixes size.

The four properties are retained separately. (C_P,C_V) are heat capacities. (C_S,C_T) below are
pressure-volume compliances, defined in a different way and carrying volume units. They are not renamed heat
capacities.

## 1. Carrier, cut, target and hypotheses

Take a twice differentiable potential (U(S,V)) at a stable point and its positive Hessian response

\[
H=\begin{pmatrix}A&B\\B&C\end{pmatrix},\qquad
A=U_{SS}>0,\quad B=U_{SV},\quad C=U_{VV}>0,\quad
\Delta=AC-B^2>0.                                           \tag{1}
\]

The carrier is the cone of these positive two-channel response blocks. The source is a displacement
((dS,dV)); the target is the pair of first-response changes ((dT,-dP)). Changing the units of (S,V)
acts by positive diagonal congruence (H\mapsto DHD). The observer keeps the ratio of the determinant to the
two diagonal responses:

\[
\boxed{\chi(H)=\frac{\det H}{H_{11}H_{22}}
=1-\frac{B^2}{AC}},\qquad 0<\chi\le1.                       \tag{2}
\]

This is TD1's (F=R-D), UP3's cross-ratio and CI1's
(exp(-2I_{\rm cut})). Write also

\[
m=1-\chi,\qquad \gamma=\chi^{-1},\qquad
I_{\rm cut}=-\tfrac12\log\chi.                              \tag{3}
\]

These are four charts of one invariant, not four freely matchable constants.

## 2. The four (C)-properties

Define the thermal pair in the usual way,

\[
C_P=T\left(\frac{\partial S}{\partial T}\right)_P,qquad
C_V=T\left(\frac{\partial S}{\partial T}\right)_V,          \tag{4}
\]

and define the second pair as dimensionful pressure-volume compliances,

\[
C_S=-P\left(\frac{\partial V}{\partial P}\right)_S,qquad
C_T=-P\left(\frac{\partial V}{\partial P}\right)_T.        \tag{5}
\]

For comparison keep their reciprocal stiffnesses

\[
K_S=-V\left(\frac{\partial P}{\partial V}\right)_S,qquad
K_T=-V\left(\frac{\partial P}{\partial V}\right)_T.        \tag{6}
\]

Equation (1) gives, without an equation-of-state model,

\[
C_V=\frac{T}{A},\quad C_P=\frac{TC}{\Delta},\qquad
C_S=\frac{P}{C},\quad C_T=\frac{PA}{\Delta},                \tag{7}
\]

and

\[
K_S=VC,\qquad K_T=\frac{V\Delta}{A},\qquad
C_SK_S=C_TK_T=PV.                                           \tag{8}
\]

Therefore the whole compass law is

\[
\boxed{
\chi=\frac{C_V}{C_P}=\frac{C_S}{C_T}=\frac{K_T}{K_S}},
\qquad
\boxed{
\gamma=\frac{C_P}{C_V}=\frac{C_T}{C_S}=\frac{K_S}{K_T}}.   \tag{9}
\]

This is the requested separation: (C_S,C_T) are defined independently of (C_P,C_V), yet the same
dimensionless compass ratio intertwines both pairs. In UP5/UP6 language the same statement is
(lambda_p/\lambda_v=z_t/\lambda_s=\chi); the reciprocal operations remain distinct.

## 3. Exact calibration theorem

For positive diagonal unit changes (D=\operatorname{diag}(a,b)), equation (2) is unchanged. Conversely, two
positive symmetric (2\times2) blocks (H,H') are related by such a unit change exactly when they have the
same (chi) and the same sign of (B). Indeed choose

\[
a=\sqrt{A'/A},\qquad b=\sqrt{C'/C}.
\]

Equality of (chi) gives (|B'|=ab|B|), and the orientation sign gives (B'=abB). Thus
((\chi,\operatorname{sign}B)) is the complete unit-orbit address. If reversal of one oriented axis is also
allowed, (chi) alone is complete.

The practical conversion table follows immediately.

| Measured anchor | Compass conversion | Recovered partner |
|---|---:|---|
| (C_V) | divide by (chi) | (C_P=C_V/\chi) |
| (C_P) | multiply by (chi) | (C_V=\chi C_P) |
| (C_S) | divide by (chi) | (C_T=C_S/\chi) |
| (C_T) | multiply by (chi) | (C_S=\chi C_T) |
| (K_S) | multiply by (chi) | (K_T=\chi K_S) |
| (K_T) | divide by (chi) | (K_S=K_T/\chi) |

One anchor is enough **within each dimensional family**. Reconstructing all four (C)-properties needs one
thermal anchor and one pressure-volume anchor. Reconstructing the full block in fixed physical coordinates
needs two axis scales and the orientation of (B).

The executable witness

\[
H=\begin{pmatrix}4&2\\2&3\end{pmatrix},\qquad
H'=\begin{pmatrix}36&30\\30&75\end{pmatrix}
\]

has (chi=2/3) in both unit systems even though every dimensional entry changes.

## 4. Comparing volume or another coefficient

Absolute volume is not invariant under a change of units. The compass compares the reduced volume

\[
v=V/V_0                                                       \tag{10}
\]

after a reference (V_0) and the remaining reduced controls have been supplied. If a material has a calibrated
branch (chi_i(v)) that is one-to-one on the region used, then the point carrying compass value (q) is

\[
v_i(q)=\chi_i^{-1}(q),\qquad
\frac{V_1(q)}{V_2(q)}=\frac{V_{0,1}}{V_{0,2}}
\frac{\chi_1^{-1}(q)}{\chi_2^{-1}(q)}.                       \tag{11}
\]

For the exact witness curves (chi_1(v)=1/(1+v)) and
(chi_2(v)=1/(1+2v)), the same (q=2/3) identifies (v_1=1/2) and (v_2=1/4). Their absolute volumes still
need (V_{0,1},V_{0,2}). If a curve is not one-to-one, a phase or branch label is part of the address.

The same rule applies to any other dimensional coefficient: compare its reduced form directly, or use one
measured member of a conjugate pair and equation (9) to recover the other.

## 5. TVSP before observation and VTSP after observation

The owner's chart involution is

\[
\Pi(T,V,S,P)=(V,T,S,P).                                      \tag{12}
\]

TVSP is the diagram before observation; VTSP is its observed chart. Equation (9) survives when the response
frame, source, metric and units are pushed forward together by (Pi). The labels in (4)-(6) carry their
definitions through the chart. A bare exchange of the letters (T,V) while leaving the derivative conditions
fixed is not the observation map.

## 6. What changes when the cuts do not commute

For a response with two different mixed readings

\[
L=\begin{pmatrix}A&B_1\\B_2&C\end{pmatrix},\qquad
\bar B=\frac{B_1+B_2}{2},\quad w=\frac{B_2-B_1}{2},          \tag{13}
\]

positive unit changes preserve two normalized numbers

\[
\beta=\frac{\bar B}{\sqrt{AC}},\qquad
\omega=\frac{w}{\sqrt{AC}},
\quad
\chi=1-\beta^2+\omega^2.                                   \tag{14}
\]

UP4 identifies (2w=[D_S,D_V]U), and a small loop reads (2w) times its area. The ratio of the two ends is
still (chi), but (chi) alone no longer classifies the response: different ((\beta,\omega)) can give the
same (chi), and (chi) can exceed one. The certificate gives an exact pair with the same (chi=3/4) and
different (omega^2). Thus the commuting compass needs one shape number; the uncut/turning compass retains an
additional signed memory number.

## 7. Adapter to another symmetry, including Yang--Mills

For any system with two selected source channels (j_1,j_2), form a positive response Gram/Hessian

\[
G_{ab}=\frac{\partial^2\log Z(j)}{\partial j_a\partial j_b}
=\operatorname{Cov}(O_a,O_b),                               \tag{15}
\]

on the declared ensemble and carrier. If (G>0), then

\[
\chi_G=\frac{\det G}{G_{11}G_{22}}                          \tag{16}
\]

is its TVSP compass address, and equation (9) defines the corresponding effective response ratios. This is the
structural adapter: carrier (G), target its positive diagonal-unit orbit, metric the covariance pairing,
source the two declared observables, intertwiner (G\mapsto\chi_G).

For Yang--Mills, the next legitimate calculation is therefore to choose two gauge-invariant source observables
and compute or bound (15). CM2's dimensionless internal-gap coefficient (>1.0112) is a spectral number, not
the determinant ratio of a declared two-source response block. Calling it (C_P/C_V) merely because both are
dimensionless would skip the intertwiner. It may be compared only after a theorem relates that spectral target
to (15), or after a response Hessian is computed independently.

## 8. Claim boundary

\[
\begin{array}{ll}
\chi=C_V/C_P=C_S/C_T=K_T/K_S & \textbf{PROVED}\\
(\chi,\operatorname{sign}B)\text{ classifies positive symmetric blocks up to positive axis units}
& \textbf{PROVED}\\
\text{one anchor recovers its partner within each response family} & \textbf{PROVED}\\
\text{reduced-volume matching on a supplied injective calibration branch} & \textbf{PROVED}\\
\chi\text{ alone classifies a noncommuting response} & \textbf{FALSE; exact witness}\\
\chi\text{ supplies an absolute coefficient or absolute volume} & \textbf{FALSE}\\
\text{an arbitrary dimensionless result equals a heat-capacity ratio} & \textbf{NOT ESTABLISHED}\\
\text{CM2's Yang--Mills gap coefficient is a response ratio} & \textbf{NOT ESTABLISHED}.
\end{array}
\]

## Reproduce

```text
cd uncut/ca1
python ca1_compass_calibration.py
python -m unittest -v test_ca1
```

The certificate has twelve exact checks and nine focused tests. Its source hashes pin UP3-UP6, CI1 and TD1.
