# YC30 — sign cuts and a vacuum-centered Fredholm family

10 October 2026. Research owner: Monty Dabas. Continues YC29 at `0b0909c`.
Keep unit-S³ electric normalization and YC27's complete physical carriers,
internal profiles and mixed-coupling windows. No predecessor packet changes.

**Decision and result.** Use a native cut-return construction, represented
by an exact Schur map, and retain its source and pairing. The sign labels
are the two objects/channels of that cut; their involution is not the
operator that exchanges them, and neither determines the energy scale.

The first proposed global normalization needs a correction. Removing only
the reference vacuum from YC28's raw odd return gives a Birman–Schwinger
operator whose norm satisfies

\[
 \boxed{\|\mathcal B_{\rm bare}(0)\|
       \ge \frac{\lambda^2}{1300}\sum_p\xi_p^2.}       \tag{1}
\]

This holds in both YC27 windows. Thus a fixed positive density of mixed
faces defeats a uniform norm-below-one test at this uncentered energy,
even though YC27 already proves a positive physical gap. Equation(1)
diagnoses the extensive reference-energy burden, not closure of the gap.

Center at the **actual interacting vacuum** and let A be the complete
physical excitation operator in its actual Hilbert pairing. Construct

\[
 \boxed{\mathcal K_\tau=\tau(A+\tau)^{-1},\qquad
 \mathfrak F_\tau=I-\mathcal K_\tau=A(A+\tau)^{-1},\quad\tau>0.} \tag{2}
\]

At each finite lattice K_tau is compact and F_tau is invertible Fredholm
of index0. More significantly, its positive contraction class is closed
under **every exact admissible nested cut**. If 0<=K_tau<=qI with q<1,
the reduced return also satisfies0<=K_eff<=qI. The exact lift and source
have bounds independent of the number of eliminations. The normalization
commutes with the energy Schur map; its derivative metric is explicitly
converted back to the physical one.

Using YC27's existing gaps, tau=4 gives full-family bounds q=175/323 or
3325/6122. These are uniform in admitted lattice volume and mixed coupling.
They re-express the existing gap theorem; they are not an independent new
gap proof. On the retained energy window <=4, the physical lift metric is
bounded by175/148 or3325/2797, again without multiplying losses per cut.

This supplies a controlled operator class and an exact cut-composition
law. It does not show that the resulting operators are local coarse
Yang–Mills Hamiltonians, enlarge the coupling windows, or construct the
continuum. The missing estimate is now localization/approximation of the
centered operator under actual spatial refinement and physical calibration.

## 1. Construct on the full carrier; then read the signs

Start with the full represented family from YC27–29,

\[
 H_\lambda(u)=H_0(u)-\lambda\sum_{p\ {\rm mixed}}\xi_pW_p,
 \qquad 0\le\xi_p\le1.                              \tag{3}
\]

All internal factor interactions remain in H0, above their actual product
ground Omega. The physical carrier includes all global centre sectors.
The domain is the inherited finite-product Laplacian domain, with bounded
Wilson potentials and normalized Haar pairing. The scalar coordinates
are adapters into this operator carrier, not a definition of the primitive
lambda-space. Cuts, ordered composition, source and pairing are carried
before taking any zero evaluation.

YC28 supplies the actual centre involution J, satisfying

\[
 J^*=J,\quad J^2=I,\quad
 P_\pm=(I\pm J)/2,\quad JH_0J=H_0,\quad JVJ=-V.       \tag{4}
\]

The labels + and - can therefore name the two channel objects. J reads
their signs; the off-diagonal source T=Pminus V Pplus crosses between
them. It is generally a map between large spaces, not a scalar sign flip.
For its polar decomposition T=U|T|, let E+=support(T*T) and
E-=support(TT*). On E+ plus E- define

\[
 X=\begin{pmatrix}0&U^*\\U&0\end{pmatrix},\qquad R=JX.
 \quad X^2=I,\quad JX=-XJ,\quad R^2=-I.               \tag{5}
\]

Extended by zero outside these supports, X² is the support projection,
not the full identity. The kernels and unmatched channels must remain.
This is T21's polar cut-bridge construction realized on the actual mixed
source. The positive modulus |T| carries information that signs alone
do not contain: replacing T by cT, c>0, keeps the polar signs but changes
every quadratic return by c². No energy or curvature magnitude is derived
from the two labels alone.

Use the retained source T when evaluating lambda=0; taking the polar
part of lambda*T and then setting lambda=0 would erase that source and
introduce an avoidable discontinuity. The physical spectral sign of a
Hamiltonian is also distinct from this independently defined cut grading.

## 2. Test the uncentered return before choosing its norm criterion

Let Pi=|Omega><Omega|, Qe=Pplus-Pi and He=Qe H0 Qe. In YC27's windows
He>=4. Let Hminus=Pminus H0 Pminus and T=Pminus V Pplus. The natural
bare reference-ground-excised operator, for real z<4, is

\[
 \mathcal B_{\rm bare}(z)=\lambda^2( H_e-z)^{-1/2}
 Q_eT^*(H_--z)^{-1}TQ_e( H_e-z)^{-1/2}.              \tag{6}
\]

It is positive and compact at each finite lattice. Compactness follows
from the reference compact resolvents and bounded finite-volume source.
Finite-volume Fredholm structure does not give a uniform reserve.

Take any normalized factor-neutral f perpendicular to Omega with
mu=<f,H0 f>. YC28 gives normalized sources Up=2WpP with exact mean
<Up f,H0 Up f>=mu+12. Different p sources are resolvent-orthogonal.
For positive Hminus-z, Cauchy–Schwarz in its square-root pairing gives

\[
 \langle U_pf,(H_--z)^{-1}U_pf\rangle
 \ge\frac1{\mu+12-z}.
\]

Test(6) on (He-z)^(1/2)f divided by sqrt(mu-z). With M2=sum xi_p²,

\[
 \boxed{\|\mathcal B_{\rm bare}(z)\|
 \ge\frac{\lambda^2 M_2}{4(\mu-z)(\mu+12-z)}.}        \tag{7}
\]

Each admitted tiling contains complementary square components. Choose
f to excite the first neutral level on one such component and leave all
other factors in their true grounds. For its single coupling u<=1/2,
the potential norm is<=u. The free neutral first excitation is12, so
min–max gives its ground-to-first separation mu<=12+2u<=13. This is
an upper bound used only to test the norm, not a new gap lower bound.
At z=0, (7) proves(1). At the completely free internal point one may use
f=2Wr Omega0 for a pure face, mu=12, giving the stronger
lambda² M2/1152.

At any fixed nonzero lambda and fixed positive density of nonzero xi,
these lower bounds grow with volume. In particular merely dropping the
reference ground does not subtract the vacuum returns on spectator blocks.
The spectral origin must be the actual energy E_lambda; the physical
vacuum must also be removed. This rejects a uniform norm test for(6)
at fixed z=0. It does not reject the exact Schur identity or tests with
a properly shifted energy and a correctly transported vacuum.

## 3. Actual-vacuum centering and its original pairing

Let Psi_lambda be the normalized true physical ground constructed by
YC27, E_lambda its energy above H0's carried reference, and
Q_lambda=I-|Psi_lambda><Psi_lambda|. On the complete excitation carrier
E_lambda_space=Q_lambda H_G put

\[
 A_\lambda=(H_\lambda-E_\lambda)|_{\mathcal E_\lambda}.
                                                               \tag{8}
\]

The symbol E_lambda_space denotes a space, not the scalar energy. A is
self-adjoint on the inherited restricted domain. YC27 supplies

\[
 A_\lambda\succeq mI,\qquad
 m=592/175\ \text{or}\ 11188/3325                     \tag{9}
\]

through mixed caps1/3200 or1/5700 in the respective internal windows.
There is no harmonic truncation, sector omission or approximate ground
projection in(8).

One can use YC19's fixed reference-quotient coordinates without discarding
their metric. Write S=exp(C), M=S*S and
G=MQQ-MQ0 M00^-1 M0Q on Q0=I-Pi. Define

\[
 \mathcal W f=S\bigl(f-\Omega M_{00}^{-1}M_{0Q}f\bigr),
 \qquad \mathcal W^*\mathcal W=G,
 \qquad \mathcal U=\mathcal W G^{-1/2}.               \tag{10}
\]

W is onto the orthogonal complement of S Omega and U is unitary from
Q0 onto the actual physical excitation carrier. If A_q is YC19's quotient
generator, then W A_q=A W and U* A U=G^(1/2) A_q G^(-1/2). This follows
also from the vanished ground row in the physical Hermitian pairing.
No uniform condition number, tensor factorization or spatial locality of
G^(1/2) is asserted. Equations(8),(10) specify precisely which metric all
norms below use.

## 4. The selected bounded family and its Fredholm status

Construct(2) by functional calculus of(8), for a declared positive
energy reference tau. Equivalently,

\[
 \mathcal K_\tau=Z_\tau^*Z_\tau,\qquad
 Z_\tau=\sqrt\tau\,(A+\tau)^{-1/2},\qquad
 \mathcal K_\tau=\tau\int_0^\infty e^{-\tau t}e^{-tA}\,dt. \tag{11}
\]

The integral is norm convergent. Its variable t is an auxiliary spectral
parameter; a physical clock requires the separate energy calibration.
This is a constructed return from the actual operator family and pairing,
not a freely fitted interaction kernel. It is a scoped YM realization;
it does not construct a universal primitive lambda-space from two signs.

By the spectral theorem, if A>=mI then

\[
 0\preceq\mathcal K_\tau\preceq qI,\quad
 \delta I\preceq\mathfrak F_\tau\preceq I,\qquad
 q=\frac\tau{m+\tau},\quad\delta=\frac m{m+\tau}.      \tag{12}
\]

The positive finite-volume elliptic excitation operator has compact
resolvent, so K_tau is compact. F_tau is invertible with inverse norm
<=1/delta and is Fredholm of index0. No trace-class claim or ordinary
Fredholm determinant is needed. Compactness need not survive a volume
limit; a surviving norm bound q<1 would still give invertibility directly.

On a real spectral interval z<m one can instead use
K_tau(z)=tau(A-z+tau)^-1 and F_tau(z)=(A-z)(A-z+tau)^-1.
All statements apply with m replaced by m-z. No exclusion above the
already known m is inferred by this substitution.

## 5. The contraction class is closed under exact cuts

Here write K for any bounded positive operator with K<=qI, q<1, and
F=I-K. Choose an orthogonal cut P,Q in its actual pairing. The hidden
block D=QFQ=I_Q-K_QQ obeys D>=delta I_Q. Exact elimination gives

\[
 \boxed{S_P(F)=I_P-K_{\rm eff},\qquad
 K_{\rm eff}=K_{PP}+K_{PQ}(I_Q-K_{QQ})^{-1}K_{QP}.}    \tag{13}
\]

The complete lift and inhomogeneous source are

\[
 W=P+Q D^{-1}K_{QP},\qquad
 f_{\rm eff}=W^*f=P f+K_{PQ}D^{-1}Q f.               \tag{14}
\]

If S_P(F)u=f_eff, the full solution is W u+Q D^-1 Q f. Keep this
particular hidden-source term; the homogeneous lift alone is insufficient.

The returned term in(13) is positive, so S_P(F)<=I_P. Completing the
square or minimizing over the hidden component gives S_P(F)>=delta I_P.
Consequently

\[
 \boxed{0\preceq K_{\rm eff}\preceq qI_P,\qquad
 I_P\preceq W^*W\preceq\delta^{-1}I_P,\qquad
 \|f_{\rm eff}\|\le\delta^{-1/2}\|f\|.}              \tag{15}
\]

For the lift bound, W* F W=S_P(F)<=I and F>=delta I. If K is compact,
the block products in(13) show K_eff is compact too.

For nested cuts the Schur map and source compose exactly, and
W_total=W_first W_second ... equals the direct lift from the final
retained space. Apply(15) to that direct cut: its metric bound is still
1/delta, **not a product of bounds from all previous steps**. This is
closure of a bounded operator class, not strict contraction q_next<q.
The cuts must be nested in one consistently identified carrier; changing
the operator or its metric between steps needs its own comparison.

## 6. Normalization commutes with the energy cut

For a positive A and an energy-admissible cut, let
A_eff=S_P(A), with the full hidden inverse and inherited form domains.
Its compressed inverse is A_eff^-1=P A^-1 P. Since
F_tau(A)^-1=I+tau A^-1, the inverse identity for the bounded Schur map gives

\[
 \boxed{S_P\bigl(A(A+\tau)^{-1}\bigr)
       =A_{\rm eff}(A_{\rm eff}+\tau)^{-1}.}          \tag{16}
\]

The effective return in(13) is therefore precisely
tau(A_eff+tau)^-1. This makes the bounded step compatible with the
original energy reduction, rather than an unrelated positive operator.
For general bounded cuts the inverse-compression formula defines the
positive reduced operator; the physical energy lift below additionally
requires the closed-form/domain admissibility specified in YC20.

The bounded-pencil lift W_F and physical energy lift W_A are different.
For an admissible cut, A W_A=P A_eff and P W_A=I. Direct substitution gives

\[
 W_F=(A+\tau)W_A(A_{\rm eff}+\tau)^{-1},\qquad
 QW_F=\tau QW_A(A_{\rm eff}+\tau)^{-1}.              \tag{17}
\]

Thus the bounded lift metric in(15) must not be called the physical
energy graph metric. The latter has its own strong, energy-relative bound:

\[
 M_A=W_A^*W_A,\qquad
 A_{\rm eff}=W_A^*A W_A\succeq m M_A.
\]

On the spectral subspace of A_eff with energies<=L,

\[
 \boxed{I\preceq M_A\preceq (L/m)I}                 \tag{18}
\]

as compressed quadratic forms. This holds for a direct cut and hence
for its exact nested composition. It is not a bound on all high-energy
retained states or on the full coordinate metric G.

For completeness, let S_A(z)=S_P(A-z) and M_A(z)=-S_A'(z).
Differentiating the actual identity(16), rather than an inequality, gives

\[
 -\partial_z S_P(F_\tau(z))
 =\tau(S_A(z)+\tau)^{-1}M_A(z)(S_A(z)+\tau)^{-1}.    \tag{19}
\]

The parent derivative is tau(A-z+tau)^-2, not I. This records the full
metric correction for an energy-dependent normalized pencil. At a
spectral root S_A(z)u=0 where the inverses remain valid, the right-hand
quadratic form is ||W_A(z)u||²/tau. The zero-vector correspondence is
therefore compatible with the original physical norm.

## 7. A complete error bound and a conditional limit passage

Suppose F=I-K and Ftilde=I-Ktilde both belong to the same positive
class delta I<=F,Ftilde<=I, with ||K-Ktilde||<=epsilon. For the same cut,
write their lifts W,Wtilde and reductions S,Stilde. Exact hidden
annihilation on each side gives

\[
 S-\widetilde S=\widetilde W^*(F-\widetilde F)W.       \tag{20}
\]

Using(15) yields

\[
 \boxed{\|K_{\rm eff}-\widetilde K_{\rm eff}\|
          \le\epsilon/\delta.}                      \tag{21}
\]

Also W-Wtilde=-Q D^-1 Q(F-Ftilde)Wtilde, so

\[
 \|W-\widetilde W\|\le\epsilon\delta^{-3/2},\qquad
 \|W^*W-\widetilde W^*\widetilde W\|
      \le2\epsilon\delta^{-2}.                       \tag{22}
\]

These constants concern a comparison of the full inputs followed by any
exact nested elimination. Separate approximations inserted at successive
steps still need a cumulative error budget; they cannot all be charged
once as epsilon. Equation(22) concerns the bounded-pencil lift metric;
physical energy transport uses(17)–(19).

If transported full operators K_n converge in norm on a common recovered
carrier with the same q<1, their limit has that reserve, and each fixed
compatible cut, source and bounded lift converges with(21),(22). This is
a conditional completion result. No such norm convergence across YM
volumes or lattice spacings has been proved here. The source, carrier
and metric identifications required in T28 remain substantive inputs.

## 8. Apply the bounds at the existing YM windows

Choose tau=4, the certified physical reference floor used by YC27, in
the present electric units. This choice is an energy reference, not a
derived dimensional mass or a clock. Equations(9),(12),(15),(18) give:

| Old block | Mixed cap | Return norm q | Reserve delta | Bounded lift metric ceiling | Physical lift metric on A_eff<=4 |
|---|---:|---:|---:|---:|---:|
| 28 links | 1/3200 | 175/323 | 148/323 | 323/148 | 175/148 |
| 64 links | 1/5700 | 3325/6122 | 2797/6122 | 6122/2797 | 3325/2797 |

All bounds retain the full excitation carrier and all allowed internal
profiles. They are uniform in the admitted volumes and the number of
exact nested eliminations. The gap input remains YC27; this table does
not enlarge its interval or prove an additional independent mass gap.

If c_n is the physical energy conversion at a refinement level and
A_phys,n=c_n A_n, a fixed physical reference tau_phys must be represented
by tau_n=tau_phys/c_n. Then

\[
 \tau_{\rm phys}(A_{{\rm phys},n}+\tau_{\rm phys})^{-1}
 =\tau_n(A_n+\tau_n)^{-1}.                            \tag{23}
\]

Keeping tau=4 in changing lattice units does not prove a uniform physical
gap. Equation(23), construction of a nontrivial physical limit, and
localized spatial transport are the remaining scale obligations.

## 9. Zero reading, verification and the next measurable task

At lambda=0, Psi=Omega, E=0, S=I, G=I and A=H0 restricted to Omega's
physical orthogonal complement. Hence K_tau recovers
tau(H0_exc+tau)^-1 exactly. It does not vanish: internal reference
dynamics remain. Its sign-cut/source history also remains. YC28's unitary
J carries Psi_lambda to Psi_-lambda and maps the actual excitation
carriers accordingly, giving K_tau(-lambda)=J K_tau(lambda)J with those
typed domains. Both the exact return and normalization commute with this
conjugation. Fixed compatible cuts have the same zero specialization
before and after(13),(16).

```sh
python physics/yc30/yc30_centered_cut.py
```

The new direct checks are exact rational identities: polar exchange with
a retained zero mode; Schur compatibility of the bounded normalization;
direct versus nested lifts and inhomogeneous sources; the derivative
metric conversion; the full-input error identities; and the table's
fractions. The all-volume lower bound(1), operator-class closure and
physical low-energy metric bounds have written proofs above. No old
suite, harmonic diagonalization or frozen certificate was rerun.

**Next task:** obtain a local representation or controlled local
approximation of this centered kernel and its physical cut maps under
actual spatial blocking. The present positive class prevents artificial
reserve loss from exact repeated elimination; it does not prove locality,
strict contraction, or preservation of the microscopic Yang–Mills form.
The gap floor used here cannot also be presented as a newly derived
consequence of the normalization. Any larger-window or continuum step
must supply fresh source/scale estimates.

Source route: T21/T23/T27 for cut-generated forms, polar support and the
derived Fredholm representation; YC19/20 for the exact ground quotient,
metric and associative Schur lift; YC27/28 for the physical windows,
four-factor mean and source orthogonality; YC29 for the retained full
neutral/charged return. The identities here are standard Schur and
spectral calculus applied with the stated native carrier and metric;
no new general Fredholm theorem or priority claim is made.

The RKF source statements were read at
`e8e3089745bf7dc0532d01c62a5f660f5658c6f0` in
`theorum/21_cut_memory_spectral_isomorphism.md`,
`theorum/23_native_cut_generated_object_theorem.md`, and
`theorum/27_rsc_primitive_to_completion_hierarchy.md`.
