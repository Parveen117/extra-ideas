# R46: an operational port for the existing native return law

Monty Dabas · 2 October 2026

**Result.** The already proved RKF R2 return coefficient has an exact positive
finite-path readout. An explicitly labelled electrical interface realizes that
readout as the input admittance of a passive resistor ladder. The native probe
has a definite component-control law, and the existing Publications third-probe
prediction receives finite-length and component-error bounds.

The new contribution is this operational adapter and its error budget. The
R2 return theorem, NI parameter reconstruction and CR third-probe polynomial
are consumed unchanged. Resistor ladders and series/parallel reduction are
established circuit theory. No new universal circuit law or external priority
is claimed. The numerical fixture below is constructed, not measured.

## Native inputs, constructed target and physical interpretation

Use the canonical cut field and R2 algebra

\[
R^2=-I,\quad K^2=I,\quad KR=-RK,\quad L=KR,\quad L^2=I.
\]

R2 derives, from paired native edge products, the return
\(F=I+xL\) and update \(f_a(x)=1/(a^{-1}+x)\), for positive radial
\(a\) and nonnegative tail \(x\). An alternating cell profile
\((a,b,a,b,\ldots)\) has the unique completed coefficient
\(bx^2+x-a=0\), selected by the native finite-aperture limit.

We construct a finite readout target with positive radial edge weights. Its
flow and balance relations are definitions on that target; their elimination
is proved below using native scalar arithmetic. They are not inferred to be
the uniquely selected material dynamics. Only **after** that derivation do we
identify the target with an ideal steady-DC resistor network. Ohmic component
semantics, the DC regime and the measurement units are physical interface
assumptions, not new primitive-only consequences. No classical circuit premise
is used to re-prove the upstream native R2 theorem.

## Written results

### R46.1 — A positive finite-path target realizes the paired native return

Take nodes \(0,1,\ldots,M\), source value \(v_0\), positive radial
series weights \(s_j\) and shunt weights \(g_j\), \(1\le j\le M\).
Define edge flow \(q_j=(v_{j-1}-v_j)/s_j\), shunt flow \(g_jv_j\),
and an optional nonnegative terminal weight \(\ell\). Internal balance is
\(q_j=g_jv_j+q_{j+1}\), with \(q_{M+1}=\ell v_M\).

Let \(y_j\) be the downstream flow/value ratio seen just after section
\(j\), with \(y_M=\ell\). Eliminating node \(j\) gives

\[
T_j(y)=\frac{g_j+y}{1+s_j(g_j+y)},\qquad
y_0=T_1\circ\cdots\circ T_M(\ell).
\tag{46.1}
\]

For native cells \(a_{2j-2}=1/s_j\), \(a_{2j-1}=1/g_j\),

\[
T_j=f_{a_{2j-2}}\circ f_{a_{2j-1}}.
\tag{46.2}
\]

Thus \(q_1/v_0\) equals R2's coefficient for the corresponding **2M native
cells** and terminal coefficient \(\ell\). The intermediate reciprocal
coordinate alternates its role; it must not be read as the same physical
quantity at every half-step.

**Proof.** Balance at one node gives
\((v_{j-1}-v_j)/s_j=(g_j+y)v_j\). Solving it proves (46.1);
the denominator is positive. Direct substitution in the two native updates
proves (46.2). Iteration gives the finite equality. Independently, the internal
coefficient matrix has quadratic form, for a zero source variation \(h_0=0\),
\[
\sum_{j=1}^M\frac{(h_{j-1}-h_j)^2}{s_j}
+\sum_{j=1}^Mg_jh_j^2+\ell h_M^2.
\]
It is strictly positive for a nonzero radial variation, hence the finite
linear equations are nonsingular. Multiplying the node balances by \(v_j\)
and telescoping yields
\[
v_0q_1=\sum_{j=1}^M\frac{(v_{j-1}-v_j)^2}{s_j}
+\sum_{j=1}^Mg_jv_j^2+\ell v_M^2.
\tag{46.3}
\]
This identity supplies positivity and the later power-balance interpretation.
Neither a primitive Hilbert space nor a physical energy postulate is needed
for these finite scalar identities. ∎

### R46.2 — The native common probe has a definite electrical control law

Introduce a positive reference resistance \(R_*\), and interpret
\[
s_j=R_{s,j}/R_*,\quad g_j=R_*/R_{p,j},\quad
v_j=V_j/V_{\rm in},\quad q_1=R_*I_{\rm in}/V_{\rm in}.
\]
Here each section is one series resistor \(R_s\), followed by a shunt
resistor \(R_p\) to the reference terminal. Later sections continue from
that junction. Under the stated ideal ohmic interface, these definitions
give exactly the target equations in R46.1, and its quadratic identity is
electrical input power equals total resistor dissipation after restoring units.

For the existing cut family, choose \(b>0\) and set
\[
\boxed{R_s(z)=\frac{R_*}{(b+1)z},\qquad
R_p(z)=R_*bz,\qquad z>0.}
\tag{46.4}
\]
Every section then has native cells \(((b+1)z,bz)\). Its completed port obeys
\[
\boxed{x(z)=R_*Y_{\rm in}(z),\qquad
bz x(z)^2+x(z)-(b+1)z=0.}
\tag{46.5}
\]

**Proof.** The inverse series coordinate is \(R_*/R_s=(b+1)z\);
the inverse shunt-conductance coordinate is \(R_p/R_*=bz\).
R46.1 and R2 SY1 now give (46.5). Conversely these two cell assignments
force (46.4), so the control law is necessary as well as sufficient for
this particular component realization. Equivalently, after fixing the
baseline, verify both
\[
R_s(z)R_p(z)=R_*^2\frac b{b+1},\qquad
R_p(z)/R_p(1)=z.
\tag{46.6}
\]
The product alone does not establish the dial calibration. Multiplying every
resistor by \(z\) instead produces cells \(((b+1)/z,bz)\), and multiplying
the input voltage only leaves the admittance unchanged. Neither is the
common native probe. ∎

### R46.3 — Finite length and component intervals remain explicit

For \(n=2M\) native cells, use R2's transfer matrix
\[
\prod_{j=0}^{n-1}\begin{pmatrix}0&1\\1&a_j^{-1}\end{pmatrix}
=\begin{pmatrix}A&B\\C&D\end{pmatrix}.
\]
Since \(n\) is even, \(AD-BC=1\), and
\[
L_M=B/D\le x\le U_M=A/C,\qquad U_M-L_M=1/(CD).
\tag{46.7}
\]
The open-ended M-section port is exactly \(L_M\); a terminal short gives
\(U_M\). Every nonnegative finite terminal load lies between them. The
infinite positive periodic response lies in the same interval by R2, not
by treating the finite device as infinite.

For uncertain positive component intervals, calculate a lower response with
all \(s_j\) at their upper endpoints and all \(g_j\) at their lower
endpoints, and an upper response with the opposite choices, retaining the
declared terminal-load bounds. This also covers nonuniform sections.

**Proof.** R2 RI1 gives the fractional-linear map and determinant formula;
even parity fixes the endpoint order. Direct differentiation of (46.1) gives
\(T_y=T_g=[1+s(g+y)]^{-2}>0\) and
\(T_s=-(g+y)^2/[1+s(g+y)]^2<0\). Composition proves the component bounds.
Increasing every physical resistor decreases the admittance. If all resistors
of an open-ended ladder are within relative \(\epsilon<1\) of their
nominal values, uniform rescaling at the two corners therefore gives
\[
\frac{x_M^{\rm nom}}{1+\epsilon}\le x_M^{\rm actual}
\le\frac{x_M^{\rm nom}}{1-\epsilon}.
\tag{46.8}
\]
For a normalized ratio between two configurations, a conservative bound
multiplies the nominal ratio by \((1-\epsilon)/(1+\epsilon)\) and its
reciprocal. This allows independently varying component errors; correlation
can improve it but is not silently assumed. ∎

### R46.4 — The existing third prediction becomes a finite port prediction

Use the previously committed NI/CR fixture, without another fit:
\[
b=5/7,\quad\delta=3/4,\quad
z_-=1/4,\quad z_+=7/4,\quad z_3=1+\delta/2=11/8.
\]
NI already proves that completed readings \(x_-=2/5\), \(x_+=6/5\)
determine this \(b,\delta\) uniquely under its symmetric-probe contract.
CR already gives
\[
55x_3^2+56x_3-132=0,\qquad
x_3=(-28+2\sqrt{2011})/55.
\tag{46.9}
\]
R46.2 now supplies the missing component/readout realization of those controls.

For the illustrative \(R_*=1000\,\Omega\), the settings are:

| Setting | z | Each series resistor, ohm | Each shunt resistor, ohm |
|---|---:|---:|---:|
| Reference | 1 | 1750/3 | 5000/7 |
| Minus | 1/4 | 7000/3 | 1250/7 |
| Plus | 7/4 | 1000/3 | 1250 |
| Third | 11/8 | 14000/33 | 6875/7 |

With 20 sections, the ideal open-ended third response and the completed
response both lie in \([1.1216063508,1.1216063510]\); the sharper R2 tail
width is below \(3\times10^{-11}\). At \(V_{\rm in}=0.1\,\mathrm V\)
the corresponding ideal current is enclosed by
\([112.16063508,112.16063510]\,\mu\mathrm A\).
These are computed design values, not laboratory readings. With 0.1% component
uncertainty, the component interval (46.8) dominates that tiny truncation error;
instrument and reference-standard errors must be included separately.

**Proof.** The component table is exact substitution into (46.4). R2's unchanged
solver encloses the third root, and the positive-root polynomial changes sign
across its rational endpoints. Independent finite node elimination gives the
open-ended endpoint. Conversion to current is \(I=V_{\rm in}x_M/R_*\).
For gain-cancelling readout \(m_z=A x_M(z)\), \(m_z/m_1=x_M(z)/x_M(1)\).
It is essential to retain the finite reference error: if
\(|x_M(z)-x(z)|\le e_z\), \(|x_M(1)-1|\le e_0<1\), then
\[
\left|\frac{x_M(z)}{x_M(1)}-x(z)\right|
\le\frac{e_z+x(z)e_0}{1-e_0}.
\tag{46.10}
\]
Subtract the two ratios and apply the triangle inequality. One may replace
\(x(z)\) on the right by its upper enclosure. For uncertain finite components
and readouts use the interval ratio \([L_z/U_1,U_z/L_1]\), with \(L_1>0\).
The finite ratio is not inserted into NI as an exact completed reading. A
common multiplicative gain cancels; offsets and configuration-dependent gains
do not. ∎

### R46.5 — Target faithfulness does not turn an inverse pole into port divergence

The port determines \(x\), and hence \(I+xL\), **within the declared unit
response line**. It does not identify arbitrary EMK coefficients, the individual
native opening/closing amplitudes, or an unknown deep profile. In particular,
the completed port at \(z=1\) has \(Y=1/R_*\), finite and positive.
The CR inverse \((I+xL)^{-1}\) nevertheless has a complementary-channel pole
as \(x\to1\). This network measures the former return coefficient, not that
inverse-channel observable.

**Proof.** The map \(x\mapsto I+xL\) is injective because \(L^2=I\).
The scalar readout extracting its L coefficient is not injective on the whole
four-channel algebra: changing the I, R or K coefficient leaves that scalar
unchanged. Also, native products factor as \(b_hc_h=a_j\), and reciprocal
rescaling of each pair preserves the product and return. Finally R2 gives
\(x(1)=1\), while native complementary idempotents give
\[
(I+xL)^{-1}=\frac{P_+}{1+x}+\frac{P_-}{1-x},\quad
P_\pm=(I\pm L)/2.
\]
Thus the pole belongs to a different observation. This is an explicit physical
readout boundary, consistent with CR-1/CR-2, rather than a particle or resonance
claim. Positive finite networks also have unique static responses: this adapter
does not supply R13's dynamical noise, a physical clock or hysteresis. ∎

## Evidence and next experimental use

[Protocol](../04-operator-evolution/r46-return-port/PROTOCOL.json),
[verifier](../04-operator-evolution/r46-return-port/verify.cjs),
[verification](../04-operator-evolution/R46_VERIFICATION.json) and
[source audit](../04-operator-evolution/r46-return-port/SOURCE_AUDIT.md) bind this
application. General proofs, exact finite checks and a future physical run are
separate. The canonical arithmetic and R2 solver are imported unchanged.

A physical run should record actual series/shunt values at each setting, input
voltage and current, terminal loading and uncertainty. Check (46.6) before
using NI's calibration. Freeze the two-reading inference and its intervals,
then compare the third reading with the finite-component prediction. A mismatch
tests the declared network/control/readout model; this classical electrical
realization is not a discriminator against circuit theory or a derivation of
alpha. No hardware run, fundamental constant or whole-framework certification
is claimed.
