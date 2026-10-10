# HC1 — the signed horizon reading and its native curvature

9 October 2026. Continues the owner's TVSP/tower/centre lead after YC10
at `b8a4805`. This is a conditional gravity/response result; it changes
no Yang–Mills spectral bound. Units are c=1.

**Result.** On the existing static spherical gravity profile, the native
curvature-square diagnostic is exactly the geometric curvature square
times `1/(2 N^8)`. It diverges at the horizon although the geometric
curvature stays finite. The folded response radius alone loses the sign
of outward propagation; one signed coordinate restores it. Thus the
owner's increasing-native-curvature lead has an exact realization, but
neither infinite matter energy nor a zero-area physical horizon follows.

## 1. Carrier, premises, and notation

Read and preserved: MO1 (falling frame), MA1 (free radial stretch under the
vacuum law), CV1 (geometric curvature), NC1 (native response curvature),
GR1/GB1 (clock and unit block), LT1 (response generation), TS1 (mirror),
HM1 (spacing-dependent memory), HX1 (algebraic cut continuation), YC10
(typed scale audit).

MA1 derives the MO1 frame **within its stationary spherical ansatz**, from
CV1/TP1's vacuum law and frames at rest at infinity. We inherit those
hypotheses; the TVSP drawing alone is not the metric construction. Choose
the positive-mass branch r_s>0, inward orientation, and time orientation
dt>0. Its regular falling coframe on r>0 is

\[
 \vartheta^0=dt,\quad \vartheta^1=dr+\beta dt,\quad
 \vartheta^2=r\,d\theta,\quad\vartheta^3=r\sin\theta\,d\phi,
 \qquad \beta=\sqrt{r_s/r}.
\]

With frame pairing diag(1,-1,-1,-1), this gives

\[
 ds^2=dt^2-(dr+\beta dt)^2-r^2d\Omega^2.                 \tag{1}
\]

For the response construction we initially restrict to **r>r_s**:

\[
 \beta=\tanh\eta,\quad N=\operatorname{sech}\eta,
 \quad N^2=1-\beta^2,\quad G=e^{\psi n\cdot C},\quad\psi=2\eta.
                                                               \tag{2}
\]

Correction to the preceding chat notation: LT1's response-disc radius is
**u=tanh(psi)=2 beta/(1+beta^2)**, not beta=tanh(eta). NC1 also uses the
letter rho for a different curvature ratio, -1/2. HC1 avoids that collision.
Here beta is a fall parameter, not YC9/YC10's incident-coupling budget.

The transfer data are explicit: the source is G(r,n) on the exterior
spatial carrier; the response connection is A=G^-1 dG/2 with scalar pairing
sc(B)=tr(B)/2. Its target is the already declared Lorentzian coframe (1).
The map uses beta^2=r_s/r and the signed coordinate below. There is **no
claimed intertwiner of the two curvature tensors or their dynamics**.

## 2. A signed circle recovers the causal side

For any beta>0 define

\[
 u=\frac{2\beta}{1+\beta^2},\qquad
 v=\frac{1-\beta^2}{1+\beta^2}.
                                                               \tag{3}
\]

**H1.** The map is a bijection from beta>0 to the open right semicircle
u>0, u^2+v^2=1. Its inverses and mirror are

\[
 \beta=\frac{u}{1+v},\quad
 \frac r{r_s}=\frac{1+v}{1-v},\quad
 \beta\mapsto\beta^{-1}\Longleftrightarrow (u,v)\mapsto(u,-v).
                                                               \tag{4}
\]

Proof: substitute (3); conversely u^2=(1-v)(1+v) and -1<v<1 give
the positive beta and the two identities. Equivalently beta=tan(alpha/2),
(u,v)=(sin(alpha),cos(alpha)), 0<alpha<pi. On the exterior, v=sech(psi)>0.
Extending (3) algebraically inside does **not** extend (2) as a positive
real static response. There is no real static observer there.

For every non-horizon u in (0,1), a sheet sign in addition to u selects
one of exactly two beta values. The nonnegative invariant 1-u^2=v^2
forgets precisely this sign. TS1 already proves this mirror for the tower;
H1 records the missing signed coordinate and its causal role.

**H2.** Radial null propagation in (1) obeys

\[
 c_+=\frac{dr}{dt}=1-\beta
       =v\frac{1+\beta^2}{1+\beta},\qquad
 c_-=-1-\beta.                                                \tag{5}
\]

For a general future causal curve, divide ds^2>=0 by dt^2:

\[
 (\dot r+\beta)^2+r^2(\dot\theta^2+
                   \sin^2\theta\,\dot\phi^2)\le1.
\]

Therefore dot r<=1-beta, including angular motion. This proves:

| Signed side | Outgoing radial direction | Causal conclusion |
|---|---|---|
| v>0, r>r_s | c_+>0 | an outward null ray can escape |
| v=0, r=r_s | c_+=0 | the outgoing generator stays on the surface |
| v<0, 0<r<r_s | c_+<0 and c_-<0 | every future causal curve moves inward |

In the stationary asymptotically flat ingoing extension (1), every
exterior radius connects by an outgoing null ray to arbitrarily large
radii, while no interior curve can do so. Hence r=r_s is its event
horizon. Without that global extension H2 only supplies the local
trapping test, not a general event-horizon theorem.

A concrete lost-sign witness: beta=1/2 and beta=2 both give u=4/5, but
v=3/5 and -3/5, and c_+=1/2 and -1, respectively. The identical folded
reading cannot tell whether escape is allowed.

## 3. The tower is not a trajectory through the horizon

LT1's infinite-lambda generation (response squaring) gives

\[
 T(\beta)=\frac{2\beta}{1+\beta^2},\qquad
 f(r)=\frac{(r+r_s)^2}{4r}.
\]

**H3.** On beta>0, T<=1, with equality only at beta=1. Also

\[
 f(r)-r_s=\frac{(r-r_s)^2}{4r}\ge0,\quad
 f(r)=f(r_s^2/r),\quad
 f'(r)=\frac14(1-r_s^2/r^2).                              \tag{6}
\]

These are the TS1 fold, now compared with H2. An interior radius is
mapped outside by this algebraic update. No future causal worldline in
(1) does that. Thus the response generation cannot be physical infalling
time evolution with r interpreted as the same areal radius at each step.
Its known reference-radius interpretation in LT1 remains intact. A
signed chart records the side; it does not change the folding map into a
causal evolution or prove traversability of TS1's mirror sheet.

On the exterior set delta=1-beta. The exact escape/clock step is

\[
 \delta_{j+1}=\frac{\delta_j^2}{1+\beta_j^2},\qquad
 \frac{\delta_j^2}{2}\le\delta_{j+1}\le\delta_j^2,
 \qquad N_{j+1}^2=v_j^2.                                  \tag{7}
\]

Do not confuse the last identity with the normalized response invariant
at the same level, 1-u_j^2=v_j^2; the static clock N_j is a different
reading. For 0<beta_0<1, beta_j=tanh(2^j artanh(beta_0)) tends to 1.

For actual outgoing propagation, put w=sqrt(r/r_s)>1. An exact primitive
of dt/dr=1/(1-beta) is

\[
 t_*(r)=r_s[w^2+2w+2\log(w-1)].                         \tag{8}
\]

The delay from r_s+epsilon to a fixed exterior radius grows as
2 r_s log(r_s/epsilon)+O(1). The outgoing speed has a simple zero,
c_+=(r-r_s)/(2r_s)+O((r-r_s)^2). This is propagation in the declared
metric, not a derivation of quantum temperature or a clock per tower step.

## 4. Native divergence and finite geometric curvature

Use NC1's **spatial Cartesian** connection A=X/2, X=G^-1 dG. The
curvature is F=-[X_i,X_j]/4. The angular dependence of n is retained;
freezing the principal axis would erase this curvature. On an orthonormal
radial/transverse basis of the flat spatial slice, its signed squares are

\[
 I_{ra}=\operatorname{sc}(F_{ra}^2)
       =-\frac{(\psi'\sinh\psi)^2}{4r^2}\quad\text{(twice)},
 \qquad I_{aa}=-\frac{\sinh^4\psi}{4r^4}.
\]

**H4.** Since psi'=-beta/[r(1-beta^2)] and
sinh(psi)=2 beta/(1-beta^2), substitution gives

\[
 I_{ra}=-\frac{r_s^2}{r^2(r-r_s)^4},\qquad I_{aa}=4I_{ra}.
                                                               \tag{9}
\]

Define a diagnostic with a stated contraction, positive **on this
family**: D=-(2 I_ra+I_aa). It is not a universal positive norm of every
native matrix field, nor a declared matter-energy density. CV1's
Levi-Civita Kretschmann scalar is K=12 r_s^2/r^6. Thus

\[
 \boxed{\mathcal D=\frac{6r_s^2}{r^2(r-r_s)^4}
                    =\frac{K}{2N^8}},\qquad
 \boxed{N^8\mathcal D=\frac K2}.                           \tag{10}
\]

This is an exact scalar relation **for this profile and these pairings**;
no general information/curvature equivalence is inferred. In particular,
the contractions on the two sides are different. A vanishing factor
N^8 regularizes this diagnostic, but is not an invertible map at the
horizon or a smooth gauge transformation removing native curvature.

As r decreases to r_s, D tends to infinity while K tends to 12/r_s^4.
Indeed `(r-r_s)^4 D -> 6`. The native connection/response has a genuine
singular limit in its declared invariant pairing, even though the
spacetime geometry is regular at that surface. Sampling this same profile
at the decreasing exterior radii of LT1, D increases without bound. Its
logarithmic derivative is -2/r-4/(r-r_s)<0; the geometric K increases only
to its finite limit. This is a profile-sampling statement. Squaring the
entire response field at a fixed base radius is a different operation,
treated next.

The diagnostic has a useful exact integrability test. The MO1 spatial
slice dt=0 has volume 4 pi r^2 dr. For r_s+epsilon<R,

\[
 \int_{r_s+\epsilon}^{R}\mathcal D\,4\pi r^2dr
 =8\pi r_s^2\big[\epsilon^{-3}-(R-r_s)^{-3}\big].        \tag{11}
\]

So a proposed functional assigning this unweighted quadratic diagnostic
as a density on that slice has an infinite horizon contribution. This
is a constraint on **that functional**, not a proof of infinite physical
energy. NC1 already finds that simple native curvature-square
stationarity does not select this gravity profile. Multiplying by N^8
changes the functional and requires a separate reason; (10) does not
authorize doing so silently.

### Fixed-base squaring also increases native curvature, but changes the profile

**H5.** Hold the base coordinates, their flat spatial pairing and the
principal-axis field n fixed, and actually square G. Then psi becomes
2 psi and psi' becomes 2 psi'. Using sinh(2 psi)=2 sinh(psi) cosh(psi),
the two native squares change by

\[
 I_{ra}[G^2]=16\cosh^2\psi\ I_{ra}[G],\qquad
 I_{aa}[G^2]=16\cosh^4\psi\ I_{aa}[G].                  \tag{12}
\]

Both are nonzero and negative on the exterior gravity family. Their
positive magnitudes grow strictly. Iterating gives D_j>=16^j D_0 at
each fixed finite exterior radius. More explicitly, for q=2^j,

\[
 \mathcal D_j=
 \frac{q^2\psi'^2\sinh^2(q\psi)}{2r^2}
 +\frac{\sinh^4(q\psi)}{4r^4}.                          \tag{13}
\]

This directly realizes increasing native curvature under the squaring
tower, **with a turning principal axis**. A fixed principal axis would
give commuting X_i and zero native curvature regardless of the growing
response eigenvalues. It is not a theorem for every response tower.

The native curvature ratio changes from -1/2 to
`r (2 psi')/sinh(2 psi) = -1/(2 cosh psi) = -v/2` after one square.
If the new response is interpreted as a new fall field on the **same**
radial base, its memory is m_new=4 r_s r/(r+r_s)^2, and

\[
 (r m_{\rm new})'=\frac{8r_s^2r}{(r+r_s)^3}>0.          \tag{14}
\]

It fails CV1's empty-space condition `(r m)'=0`. Thus (10) applies to
the original gravity profile and its radial samples; it cannot be reused
unchanged for G_j at the original fixed r. Conversely, changing radius
to f(r) as a coordinate map requires transporting derivatives, metric and
measure; declaring the new radius flat afresh is a new physical reading.
This separates three operations: fixed-base response generation,
sampling a field at another place, and transporting its coordinate chart.

## 5. What shrinks, and what an observer can do

The determinant of the radial metric block in (1) is -1, including at
beta=1. Its coefficients and inverse are smooth there. A falling clock
dr/dt=-beta has ds=dt and crosses the horizon in finite proper time;
from r_0 to r_1 it takes
`2(r_0^(3/2)-r_1^(3/2))/(3 sqrt(r_s))`. The physical horizon area is
4 pi r_s^2, not zero. The geometric singular limit of this solution is
r=0, where K diverges.

The static exterior observer has clock ds=N dt and ceases to be a
timelike static observer at the horizon. GR1's static reading factor
1/N diverges there; this is not an identification with matter density.
At dt=0 radial proper length is dr, whereas on a static-time slice it
is dr/N. Length comparisons need an observer and simultaneity convention.
There is no derived rule that every physical length contracts with the
drawn response circle.

Even N alone cannot determine a cone: on the exterior,
`ds^2=N^2 dT^2-A^2 dr^2-r^2 dOmega^2` has radial speed N/A. Choosing
A=N or A=1/N gives 1 or N^2 with the same static clock. These are two
different geometries, not interchangeable gauges of (1). GB1's spatial
frame input is essential; no regular interior extension of the first
alternative is asserted.

In the signed circle, u->0,v->+1 is spatial infinity, whereas u->0,v->-1
is the r->0 singular limit. The horizon is u=1,v=0. Thus a **folded
disc centre alone is ambiguous**, and in LT1's chart the horizon is the
boundary, not the centre. Another centre chart is possible only with its
map stated. Nothing here identifies the thermodynamic ST-centre potential
itself with a black hole.

## 6. Status, classical identification, and Yang–Mills handoff

The classical target is Schwarzschild in ingoing Painleve–Gullstrand
coordinates. Hamilton and Lisle, *The river model of black holes*,
arXiv:gr-qc/0411060v2, equations (1)–(2), use the same falling form with
the opposite overall metric-sign convention:
<https://arxiv.org/pdf/gr-qc/0411060>. H1–H5 are the native dictionary,
causal-sign completion and contraction audit here, not a new black-hole
solution or a new empirical prediction.

| Statement | Status |
|---|---|
| Signed response reconstructs side and null-direction sign | Written exact proof and symbolic/rational checks |
| Tower folds the two sides and cannot itself be infalling time | TS1 mirror plus the explicit causal comparison |
| Native D=K/(2 N^8), with divergent unweighted spatial integral | Exact on the declared exterior profile and contractions |
| Fixed-base squaring amplifies native curvature and fails the same vacuum profile | Exact factors (12), with the nonzero residual (14) |
| Horizon finite area and finite geometric curvature | Known result of the inherited spherical metric |
| Native response smoothly continues as the same positive static field inside | Not established; its exterior chart degenerates |
| ST centre is physically a black hole; matter energy diverges there | Not established |
| Native lambda is YM scale flow; horizon yields a quantum mass gap | Not established |

YM keeps YC10's certified small-coupling, volume-uniform finite-lattice
vacuum bounds. The next spectral target is its correlated reference block:
certify the block's actual local gap and the full interface source/inverse
budget. HC1's useful discipline for that step is to preserve sheet/sign
information and distinguish a singular response normalization from the
operator's physical rate. No gravity-to-YM operator intertwiner is supplied.

## Reproduce

From this directory:

```sh
python hc1_signed_horizon.py --check
python -m unittest test_hc1
```

Thirty-eight exact checks and nine focused tests pass. The result JSON pins
the source packets, this note, code and tests by SHA256.
Exact computations audit identities and counterexamples; the causal and
limit arguments are written above. This is not machine-formal verification.
