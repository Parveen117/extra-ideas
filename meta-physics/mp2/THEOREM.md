# MP-2: evolving response with protected curvature and retained energy records

Research owner: Monty Dabas. Development: 3 October 2026. Python 3.12 only.

## Advance and contract

This packet extends the **frozen-response free protocol** in Publications
YM54 to an explicit evolving Hessian-response family. A selected relaxation
law conserves its native curvature marker, decreases its response budget,
accounts for that decrease in paired records, and retains a uniform decay
bound for the time-ordered free heat process. The moving process has an
additional, signed driving-work term in its derivative-energy balance.

The two main conclusions are

\[
 d(t)=d_0\ne0,\quad B(t)\le B_0,
 \qquad \|U(t,s)f\|_\Phi\le e^{-\gamma_*(t-s)}\|f\|_\Phi,
 \quad \Phi(f)=0,
\]
\[
 \gamma_* = \min\left\{\frac{|d_0|}{4},\frac{d_0^2}{2B_0}\right\}>0.
\]

The law is derived from a **declared selection criterion**: least change
from radial response relaxation while conserving oriented shape area.
This criterion, an optional frame rotation, the initial source and the heat
clock are inputs. Nature's choice of them is not derived. This is a concrete
admissible dynamics, not a universal law for every thermodynamic material.

Reuse NT-1--NT-4, YM50--YM52 and YM54 at Publications commit
`b7160ae47c76342088ca37334775a6e8458ed2c6`. In particular,

\[
 K^2=L^2=1,\quad R^2=-1,\quad L=RK,\quad[K,L]=-2R,
 \quad Y_i=p_iK+q_iL.
\]

Let V have columns `(p_i,q_i)`, `d=det V`, `B=||V||_F^2`, and use the
existing native matching norm. The symbol B in this packet is a scalar
budget, **not** NT's positive factor. At the response origin that factor
will be I. YM54 supplies, without a new algebra,

\[
 f=d/2,\quad G=2V^TV,\quad\tau=2B,\qquad
 C=\frac12\begin{pmatrix}VV^T&0\\0&0\end{pmatrix}.
 \tag{1}
\]

Its marker f is for the selected connection `A_i=H^-1 partial_i H/2`.
The rank-two compact tensor C is a counted-protocol tensor, not spacetime
curvature. Let J be the coefficient-plane turn `[[0,-1],[1,0]]`.
Conjugation of the native shape words by `Exp_Sigma(theta R/2)` induces
`Exp_Sigma(theta J)` on their coefficients. Thus the frame rotation below
uses an existing native operation.

We use the earned scalar completion, roots and factorial exponential from
Extra Ideas R16/R42. Time is a declared control/heat parameter. All heat
operators act on the **same** YM50 compact carrier with the **same** Phi.
Their coefficients depend on time only, not on the compact state. No
primitive Hilbert space or replacement native engine is introduced.

## T1 — Every evolving shape pair has an actual local Hessian realization

Write V=`[[p1,p2],[q1,q2]]`. At each time choose the dimensionless potential

\[
 u_t(x,y)=\frac{x^2+y^2}{2}+\frac16\left[
 (2p_1+q_2)x^3+3q_1x^2y+3q_2xy^2+(q_1-2p_2)y^3\right].
 \tag{2}
\]

Then H=`D^2 u_t` has H(0)=I and, at the origin,

\[
 H_x=\begin{pmatrix}2p_1+q_2&q_1\\q_1&q_2\end{pmatrix},\quad
 H_y=\begin{pmatrix}q_1&q_2\\q_2&q_1-2p_2\end{pmatrix}.
\]

Removing their scalar parts gives exactly Y1 and Y2. Thus (1) is realized
by integrable response jets, not two unrelated assigned matrices. If
`B(t)<=Bmax` and a positive supplied bound r satisfies `r^2>=Bmax`, then
on the common patch `|x|+|y|<=1/(6r)` one has

\[
                       \tfrac12 I\preceq H_t\preceq\tfrac32 I.
 \tag{3}
\]

**Proof.** Differentiate (2). The scalar parts are
`sigma_x=p1+q2`, `sigma_y=q1-p2`. Each is at most `sqrt(2Bmax)` in
magnitude. A trace-free shape word has operator norm at most `sqrt(Bmax)`,
so each Hessian derivative has norm at most `3r`. Since H is affine in
(x,y), `||H-I||<=3r(|x|+|y|)<=1/2`, proving (3).
With H(0)=I, NT gives `F_xy=-[Y1,Y2]/4=(d/2)R`; hence f=d/2.
This realization is a time-dependent constitutive family evaluated at a
fixed chart origin. It is not claimed to be a trajectory in a single
autonomous material equation of state, nor to control all time-space
curvature components away from that origin.

## T2 — A selected flow preserves curvature and reduces response budget

For a two-by-two V put

\[
 V^\#=\operatorname{cof}V
 =\begin{pmatrix}V_{22}&-V_{21}\\-V_{12}&V_{11}\end{pmatrix}.
\]

This cofactor map is linear, self-adjoint and involutive for the matching
pairing. Its norm square is B and `<V,V#>=2d`. Fix `kappa>=0` and a
continuous, locally bounded angular rate omega(t). Define

\[
 \boxed{\dot V=\omega J V-\kappa\left(V-\frac{2d}{B}V^\#\right).}
 \tag{4}
\]

For d0 nonzero it has a global solution, with

\[
 \dot d=0,\qquad
 \dot B=-2\kappa\frac{B^2-4d_0^2}{B}\le0,
 \qquad
 B(t)^2=4d_0^2+(B_0^2-4d_0^2)e^{-4\kappa t}.
 \tag{5}
\]

**Selection proof.** Among velocities Z satisfying `<V#,Z>=0`, minimize
`||Z+kappa V||^2`. The unique orthogonal projection is
`Z=-kappa(V-(2d/B)V#)`. This is the declared selection rule. Adding
`omega J V` does not change d or B, because determinant is unchanged
by a determinant-one rotation and `<V,JV>=0`.

**Balance and existence proof.** Differentiate determinant by its cofactor
and differentiate B by `2<V,Vdot>`. The cofactor identities give (5).
Also `B>=2|d|`, since `||V-V#||^2` and `||V+V#||^2` are nonnegative.
For explicit existence, remove the frame rotation by O(t) solving
`Odot=omega J O`, namely the native exponential of `J integral omega`.
Set `V_+=(V+V#)/2`, `V_-=(V-V#)/2` and initial weights
`a_+0=||V_+0||^2`, `a_-0=||V_-0||^2`. They are orthogonal and
`a_+-a_-=2d0`, `a_++a_-=B`. Let
`a_+(t)=(B(t)+2d0)/2`, `a_-(t)=(B(t)-2d0)/2`. Then

\[
 V(t)=O(t)\left[
 \sqrt{a_+(t)/a_{+0}}\,V_{+0}
 +\sqrt{a_-(t)/a_{-0}}\,V_{-0}\right].                 \tag{6}
\]

An initially zero summand stays zero and its fraction is omitted.
Differentiate (6) to check (4). The weights are positive for every finite
time when initially positive; B stays above `2|d0|`. This constructs a
global solution directly in the earned scalar completion. If kappa>0,
the smaller protocol eigenvalue increases toward `|d0|/2`; frame rotation
can continue. Area preservation is not an assertion that arbitrary
dissipative evolution preserves native curvature.

## T3 — Equal paired loss has an exact retained-record account

In the co-rotating frame, both branch weights obey

\[
 \dot a_+=\dot a_-=-\frac{4\kappa a_+a_-}{a_++a_-}.
\]

Therefore the same nonnegative record

\[
 m(t)=\frac{B_0-B(t)}2,\qquad
 a_+(t)+m(t)=a_{+0},\quad a_-(t)+m(t)=a_{-0},
 \qquad \frac{B(t)}2+m(t)=\frac{B_0}2                 \tag{7}
\]

accounts for the removed amount from each branch. The conserved difference
is exactly twice d0. At finite times the current V, the cumulative frame
O and m recover the initial V by undoing O, splitting into its two
cofactor parts, and multiplying each nonzero part by
`sqrt((a_+(t)+m)/a_+(t))` or `sqrt((a_-(t)+m)/a_-(t))`.

**Proof.** Project (4) into the two cofactor eigenspaces, after removing O:
`Vdot_+=-2kappa a_-/B V_+`, `Vdot_-=-2kappa a_+/B V_-`.
Taking squared norms proves the common loss and then (7). Formula (6)
gives the inverse. An initially absent branch has m=0 and is omitted.
At infinite relaxation one branch can vanish despite nonzero initial
weight: then m retains its weight but not its direction. Recovering that
direction needs an additional record. We do not call scalar energy
accounting complete memory, nor identify m with physical entropy or a
retarded memory kernel. The distinct observable variance in T5 is derived
from actual counted turn records.

## T4 — The evolving free protocol retains a uniform contraction bound

At each t, apply the four equally counted native turns of YM54 with the
current two columns of V(t). Write U(t,s) for their time-ordered heat
limit; later intervals act on the left. The tensor C(t) in (1) has positive
eigenvalues

\[
 \beta(t)=\frac{B-\sqrt{B^2-4d_0^2}}4,
 \quad c(t)=\frac{B+\sqrt{B^2-4d_0^2}}4.
\]

YM52's exact frozen-generator rates therefore give

\[
 \gamma_{\rm full}(t)=\min\{B(t)/8,\beta(t)\},\qquad
 \gamma_{\rm even}(t)=\beta(t),\qquad
 \beta(t)\ge\frac{d_0^2}{2B_0}.
\tag{8}
\]

For centered f in the existing recognition completion,

\[
 \boxed{\|U(t,s)f\|_\Phi\le
 \exp\left(-\int_s^t\gamma_{\rm full}(r)\,dr\right)\|f\|_\Phi
 \le e^{-\gamma_*(t-s)}\|f\|_\Phi,}
\quad
 \gamma_* =\min\{|d_0|/4,d_0^2/(2B_0)\}.                 \tag{9}
\]

The even sector has the analogous integral bound with beta. The tensor
may rotate: the generators at distinct times need not commute. There is
no single stationary generator whose spectral gap is asserted in (9).

**Construction and proof.** On polynomials of degree at most n, native
derivatives preserve a finite-dimensional space and have matching gain
at most n/2 per unit direction. Thus
`||L_C-L_C'||<=n^2 sum_ab |C_ab-C'_ab|/4`. Freeze C on a partition of
[s,t]. Its ordered products of the previously constructed contractions
are Cauchy; two coefficient paths differ by at most
`(t-s)n^2 sup_r sum_ab |C_ab(r)-C'_ab(r)|/4`, by telescoping/Duhamel on
this finite space. Continuity of C makes the bound vanish on refinement.
Each frozen heat factor is itself the YM54 native counted-turn limit.
For a partition interval of length h_j, its finite-count error is at most
`5 n^4 h_j^2 (2B0)^2/(1536 N_j)`; summing those errors gives a genuine
ordered-word approximation. These are degree-dependent bounds, not
all-degree operator-norm errors at time zero.

The limit is positive, unital and Phi-preserving in the uniform completion,
and contractive in the recognition completion. The existing polynomial
density extends it to both, consistently. For centered polynomial f_t,
`d||f_t||_Phi^2/dt=-2 Phi(f_t L_C(t) f_t)`. Apply the existing common-Phi
coercivity (8), integrate, and then complete. There is no derivative of a
moving vacuum/norm because Phi and its constant vacuum are fixed. Using
`B>=2|d0|` and `B<=B0` proves (9).

This removes freezing for this **free, time-driven** model only. The
interacting chain has a changing ground source when coefficients change;
YM53's fixed-source interacting bound cannot be transferred by this proof.

## T5 — Moving response contributes signed work; noise keeps an exact balance

Retain YM52's native product form Gamma_C and derivative energy
`E_C(f)=Phi(Gamma_C(f))`. For f_t=U(t,0)f and differentiable C,

\[
 \frac{d}{dt}\Phi(f_t^2)=-2\mathcal E_{C(t)}(f_t),
\]
\[
 \boxed{\frac{d}{dt}\mathcal E_{C(t)}(f_t)
 =\underbrace{\Phi(\Gamma_{\dot C(t)}(f_t))}_{\text{driving work}}
   -2\Phi((L_{C(t)}f_t)^2).}                              \tag{10}
\]

Gamma for dot C is its linear extension, not a claimed positive form.
Consequently derivative energy need not decrease when the response is
driven, even though retained squared distinction decreases.

For the actually discarded turn records define

\[
 N_t(f)=U(t,0)(f^2)-(U(t,0)f)^2
       =2\int_0^t U(t,r)\Gamma_{C(r)}(f_r)\,dr\ge0,
\]
\[
 \boxed{\Phi(f^2)=\Phi(f_t^2)+\Phi(N_t(f)).}             \tag{11}
\]

**Proof.** Differentiate the first two expressions using the product rule,
symmetry of the instantaneous L in the same Phi pairing, and
`L(fg)=fLg+gLf-2Gamma(f,g)`. This gives (10) and
`Ndot=-L_C(t)N+2Gamma_C(t)(f_t)`, with N0=0. T4 variation of constants
gives (11); positivity and Phi invariance give its sign and mean balance.
For polynomial f all terms live in finite degree (use degree 2n for f²).
At finite approximation, the ordered turns have a joint labelled word
record; its nonnegative conditional squared deviation is exactly the
difference between the mean squared record and squared mean. Its limit is
N. Thus the discarded record term is derived, not inserted as free noise.
Derivative identities are stated on the polynomial/common smooth core;
the integrated squared-norm balance extends by completion where defined.

An exact counterexample to omitting work: V=diag(4,1), kappa=10, omega=0
gives `C=diag(8,1/2,0)` and `Cdot=diag(-2400/17,150/17,0)`.
For `f=q0^2+q1^2-(sum qj^2)/2`, Phi(f²)=1/12, E_C(f)=1/24,
the work is 25/34, and the derivative in (10) is **283/408>0**.
The retained squared-norm derivative is still -1/12. The executable
control uses the unchanged native polynomial/reference engine.

## T6 — At fixed curvature the full and sign-blind observers prefer different shapes

For fixed nonzero d let b<=c be the two nonzero eigenvalues of C.
Their product is d²/4. The sign-blind rate is b, maximized at b=c=|d|/2.
The full-observer rate is `min((b+c)/4,b)`, with maximum

\[
 \boxed{\max\gamma_{\rm full}=\frac{|d|}{2\sqrt3},
 \qquad c=3b,\quad b=\frac{|d|}{2\sqrt3}.}              \tag{12}
\]

Thus full-observer optimization at fixed curvature does not select equal
shape strengths. The ratio three in (12) is a dimensionless optimum for
this declared objective, not a coupling constant of nature.

**Proof.** Put p=|d|/2 and c=p²/b, with `0<b<=p`. The two arguments
of the minimum cross at `b=p/sqrt3`. Below it the minimum is b and
increases; above it the minimum is `(b+p²/b)/4`, which decreases up to
p. Hence the maximum occurs at the crossing. The balanced-shape rate
is |d|/4. A relaxing trajectory starting sufficiently anisotropic passes
the full-observer optimum before reaching the even-observer optimum.
The underlying eigenvalue formula is credited to YM52/Lauret; (12) is
an elementary fixed-area optimization of that formula.

## T7 — Integrated source error has a quantitative survival budget

Add a continuous perturbation E(t) to (4), using the **current** d and B
inside its feedback. Suppose `||E(t)||_F<=e(t)`, choose `r0^2>=B0`, and
let `epsilon(t)=integral_0^t e(r) dr`. Then

\[
 \sqrt{B(t)}\le r_0+\epsilon(t),\qquad
 |d(t)-d_0|\le r_0\epsilon(t)+\epsilon(t)^2/2.            \tag{13}
\]

For a horizon with epsilon<=epsilon_* define

\[
 D_*=|d_0|-r_0\epsilon_*-\epsilon_*^2/2>0,\quad
 B_*=(r_0+\epsilon_*)^2.
\]

The source stays full rank with the same oriented sign, T1 has a common
stable patch, and T4 holds on that horizon with

\[
 \boxed{\gamma_{\rm robust}\ge
       \min\{D_*/4,D_*^2/(2B_*)\}>0.}                  \tag{14}
\]

Its open record balances are
`ddot=<V#,E>` and `(B/2+m)dot=<V,E>`, where
`mdot=kappa(B-4d²/B)>=0`. The error acts as an explicit area/energy source.

**Proof.** The unperturbed terms do not increase B, so
`(sqrt B)dot<=e`. Determinant differentiation gives
`|ddot|<=sqrt B e<=(r0+epsilon)e`. Integrate both inequalities to obtain
(13). The strict D_* gate prevents a zero determinant, justifying the
estimates throughout the horizon by continuation; bounded B and the
locally bounded coefficients prevent finite-time escape. Apply the same
instantaneous common-Phi inequalities with `|d|>=D_*` and `B<=B_*`.
The balance identities follow by differentiating exactly as in T2/T3.
An all-time conclusion requires an all-time integrated error budget; a
small fixed error rate alone does not supply one.

## Source repair, duplication boundary and unresolved physics

UAL sources **033** (gauge/spectral mixing), **053** (dissipative splitting)
and **093** (sectoral curvature balance) motivate this shared development.
Their original formulas are not promoted wholesale. In particular:

- Conservation of d is selected by (4), not inferred from action stationarity
  or a sum of differently dimensioned physical constants.
- Positive response loss alone does not imply a curvature floor:
  V=diag(1,exp(-t)) has decreasing B but d tends to zero. Fixed d alone
  is also insufficient: V=diag(exp(t),exp(-t)) has unbounded B and a
  vanishing smaller protocol eigenvalue. These obstructions are already
present in YM54 and are reused as controls, not new discoveries.
- R10 already proves a sourced native Bianchi identity for a declared
  connection. T2 is a different statement: conservation of one constitutive
  response marker along a selected control-time flow, with a budget bound.
  It neither replaces R10 nor infers time conservation from Bianchi alone.
- An Euler step need not preserve determinant even for tangent velocity.
  The exact pair-step checker enforces the invariant; it is not falsely
  labelled a generic structure-preserving numerical integrator.
- The scalar record m, observed turn variance N, and a retarded memory
  kernel are different objects. UAL 051's cross-time/passivity problem
  remains open. MP-2 does not claim to solve it by adding a bookkeeping term.
- Response energy B/2, observable derivative energy E_C and physical
  energy in joules are different typed quantities. Calibration remains open.
- State-dependent diffusion coefficients, physical sector selection,
  interacting moving-vacuum estimates and four-dimensional Yang--Mills
  remain open. No Yang--Mills branch is advanced by this application packet.

The cofactor projection uses standard constrained-gradient mathematics;
see Absil, Mahony and Sepulchre, *Optimization Algorithms on Matrix
Manifolds*, chapter 3 ([authors' book page](https://sites.uclouvain.be/absil/amsbook/)).
Lauret's [homogeneous-sphere eigenvalue paper](https://arxiv.org/abs/1801.04259)
is credited through YM52 with its normalization and singular-extension
contract. These are antecedents, not claims of new general methods.
The advance here is the explicit native Hessian source, its chosen evolving
law and record account, and the proof that the existing free protocol and
its work/noise balances extend together under quantified errors.

## Replay and certification scope

`SOURCE_PINS.json` pins the exact native dependencies and prior evidence.
The verifier imports the unchanged Publications engines; none is copied
into this folder. It replays MP-1 on its frozen inputs and checks R1--R46.
The original MP-1 folder/index/audit stay byte-identical; the live sequence
is listed in [INDEX.md](../INDEX.md).

From the Extra Ideas root, with Publications checked out at the commit
above as a sibling directory named `Publications`:

```sh
python3.12 -B meta-physics/mp2/verify.py --publications-root ../Publications --check
```

The seven general statements have written proofs; exact finite controls
test implementations and boundary examples. This is not formal
proof-assistant verification, independent peer review or experimental
confirmation. The packet's selection assumptions remain visible.
