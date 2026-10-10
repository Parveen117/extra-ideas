# YC23 — primitive carrier, cut energy and the vacuum gap

10 October 2026 (India). Continues YC22 and the observation discussion at
`0c6bdb0`. Earlier packets are unchanged. Electric normalization remains
the unit-S³ Laplacian, Wp=Tr(Up)/2. This is a written analytic result with
exact computational controls, not a proof-assistant certificate.

**Result.** A specified cut instrument on the actual finite-lattice YM
vacuum has an exact energy cost. For smooth two-outcome cuts it is one
quarter of the vacuum-averaged Fisher metric, or one eighth of the
squared derivative of a lifted involution. The physical vacuum gap is
exactly the infimum of energy cost divided by the probability of leaving
the vacuum, over these calibrated physical probes. This supplies an
observation-to-spectrum bridge, not a larger gap window or continuum proof.

## 0. Primitive order and what is being proposed

The owner now proposes **lambda-space as primitive**, with cuts acting on
it and a first agitation after a cut. This supersedes treating a scalar
potential or a real interval of lambda values as the primitive object.
Use a provisional symbol A_Lambda for the carrier of admissible operations
and records. A selected cut need not already be present; the meaning of
noncommutation depends on which operations are compared. Do not identify
this carrier with the earlier source amplitude, tower order, observation
zero cut, time or RG coordinate.

Research decision: represent a vacuum by a positive normalized reference
state omega_0 on the operational carrier, not by equating one state with
the entire carrier. A product, dagger, positivity and a dynamics/energy
form must be supplied or derived before spectral words apply. No theorem
here uniquely constructs these from the word "lambda-space". R15's typed
cut-record/topology construction remains valid on its declared carrier;
it is not replaced by an assumed smooth connected lambda manifold.

For the concrete YM adapter these objects already exist. A Hilbert-space
representation is used to calculate; it is not asserted to be foundational.
Three cut uses must remain distinct: a relabelling, an algebraic retained
channel, and an implemented state-changing instrument. Only the last is
assigned an injected energy below. The formulas describe system energy
after that instrument; they do not derive apparatus dynamics or a
spontaneous first event from the vacuum.

## 1. What a sharp represented cut does to a ground state

Let h=H-E0>=0 have normalized unique vacuum Omega. Let K=K*=K^-1,
P±=(I±K)/2 and assume K Omega is in the quadratic-form domain of h.
Implement the nonselective projective instrument

    rho_cut = P+ |Omega><Omega| P+ + P- |Omega><Omega| P-.

The energy increase, a sum of quadratic forms, is

\[
 E_{\rm cut}=\tfrac12\|h^{1/2}K\Omega\|^2.                 \tag{1}
\]

Proof: P±Omega=(Omega±KOmega)/2 and h^(1/2)Omega=0. Both branches
contribute one quarter of the same norm. If K Omega is in D(h), this is
also omega_0([K,[h,K]])/4. The first energy moment is (1); the squared
norm of [h,K]Omega is the **second** moment, not the same energy.

Write k=omega_0(K), Var(K)=1-k². Then

\[
 q_{\rm exc}:=1-\langle\Omega,\rho_{\rm cut}\Omega\rangle
       =\tfrac12(1-k^2),\qquad
 E_{\rm cut}\ge\Delta\,q_{\rm exc}.                     \tag{2}
\]

This follows by applying the spectral gap inequality to (K-k)Omega.
A balanced sharp cut k=0 has q_exc=1/2. Equation (2) uses a proved gap;
it does not produce a positive gap from balance. Zero energy is determined
on the vacuum: [h,K] can be nonzero elsewhere while K Omega=±Omega.
For nonunique vacua replace the rank-one vacuum projection accordingly;
that case is not used below. Arbitrary discontinuous multiplication cuts
need not satisfy the form-domain assumption.

## 2. An admissible smooth cut on the complete YM carrier

On a finite connected SU(2) link graph take H=-sum_e Delta_e+V, with real
smooth gauge-invariant V. YC22 supplies the unique positive smooth vacuum
Omega, dnu=Omega²dmu, and the exact ground-state form

\[
 \langle\Omega f,h\Omega f\rangle
 =\mathcal E_\nu(f)=\sum_e\int|\nabla_e f|^2d\nu.        \tag{3}
\]

Let f be real smooth gauge invariant, omega_0(f)=0 and Var_nu(f)>0.
Fix eta with 0<|eta| ||f||_infinity<1. This **probe strength eta** is not
identified with primitive lambda-space or the owner's observation lambda.
Take positive multiplication operators

\[
 M_\pm=\sqrt{(1\pm\eta f)/2},\qquad
 \rho_\eta=\sum_{s=\pm}|M_s\Omega\rangle\langle M_s\Omega|.
                                                               \tag{4}
\]

They are smooth and preserve the physical gauge sector. The two outcome
probabilities are exactly 1/2, at every strength, since f is centered.
Their conditional states and energy need not agree. Using (3) separately
on both branches, including the full interacting vacuum, gives

\[
 \boxed{E_\eta=\frac{\eta^2}{4}\int
       \frac{\sum_e|\nabla_e f|^2}{1-\eta^2f^2}\,d\nu.}  \tag{5}
\]

The potential has not been dropped: it determines Omega and nu, and cancels
explicitly only through the exact ground transform. The assumptions are
finite-volume ellipticity and smoothness, not a truncation of link harmonics.
In particular E_eta>0 for every nonconstant f, and E_eta tends to zero as
eta tends to zero. Balanced outcomes alone impose no strength-independent
minimum agitation. At eta=0 both branches are I/sqrt(2): this is an
uninformative probe, **not** a proposed model of pure lambda-zero observation.

## 3. The energy is the variation of a cut, with its full metric

Pointwise set u=eta f, v=sqrt(1-u²), W=(M+,M-)^T and Pi=WW^T. The lifted
cut and the fixed outcome readout are

\[
 K_\eta=2\Pi-I=\begin{pmatrix}u&v\\v&-u\end{pmatrix},
 \quad K_\eta^2=I,\qquad J=\operatorname{diag}(1,-1),
 \quad W^T J W=\eta f.                                  \tag{6}
\]

K_eta is the retained-image cut of the lift; J is the outcome cut. The
compressed readout eta f is not itself an involution. This distinction
prevents replacing the smooth instrument by a discontinuous sharp cut.

For configuration derivatives, the complete horizontal tensor is
G_ij=(partial_i W)^T(I-Pi)partial_j W. Here W^T dW=0, so

\[
 G_{ij}=\frac{\eta^2\partial_i f\partial_j f}
                    {4(1-\eta^2f^2)},\qquad
 E_\eta=\int\operatorname{tr}_{\rm electric}G\,d\nu
       =\frac18\int\sum_e\operatorname{tr}
                    ((\nabla_e K_\eta)^T\nabla_e K_\eta)\,d\nu.
                                                               \tag{7}
\]

All vector-field components of each link gradient are summed. The electric
metric and its energy unit are inputs from the actual YM Hamiltonian.
No dimensional mass is supplied merely by K_eta²=I.

For outcome probabilities p±=(1±eta f)/2 the classical Fisher tensor
I_ij=sum_s(partial_i p_s)(partial_j p_s)/p_s equals 4G_ij. Thus (7) is an
explicit information-to-energy equality for this protocol, with factor1/4.
Fisher information here means local sensitivity of the two outcome laws;
it is not their entropy or a claim of complete state reconstruction.

**Curvature control.** The projected connection W^T dW is zero and its
curvature is G_ij-G_ji=0, while (7) is positive. Also, for the native real
turn R, (R K_eta)²=I: the UP8 return-metric half-share connection is flat
for this self-adjoint cut family. This is an exact example of agitation
with zero curvature of these specified connections. The energy is the
symmetric variation; the ordered curvature is a different reading of the
full response. This does not assert that every other geometry on the
primitive carrier is flat, or equate this protocol with UP8's deformed
cross-corner source. Keep the connection type attached to each assertion.

## 4. Calibrated agitation is exactly the physical gap target

The probability of leaving the unique vacuum is

\[
 q_\eta=1-\sum_s\omega_0(M_s)^2
       =\sum_s\operatorname{Var}_\nu(M_s)>0.              \tag{8}
\]

Applying the physical spectral inequality to each centered M_s gives
E_eta>=Delta_G q_eta. Conversely, the uniformly convergent square-root
expansion on |eta| ||f||<1 gives

\[
 E_\eta=\tfrac{\eta^2}{4}\mathcal E_\nu(f)+O(\eta^4),
 \qquad q_\eta=\tfrac{\eta^2}{4}\operatorname{Var}_\nu(f)
                         +O(\eta^4),
\]
\[
 \boxed{\Delta_G
   =\inf_{f,\eta}\frac{E_\eta}{q_\eta}
   =\inf_f\frac{\mathcal E_\nu(f)}{\operatorname{Var}_\nu(f)}.} \tag{9}
\]

The infimum is over the admissible centered real smooth gauge-invariant
nonconstant f and nonzero probe strengths. Proof of the reverse bound:
fix f and let eta tend to zero; then use density of smooth gauge-invariant
functions in the physical form domain. On each finite compact carrier
Omega and Omega^-1 are smooth bounded, and gauge averaging preserves the
form core. No volume-uniform bound on either multiplication is needed.
Real functions suffice because the real and imaginary Rayleigh forms add.

An individual probe ratio is therefore an **upper** variational bound on
Delta_G. A lower bound requires the inequality for every physical probe.
YC22's strong rectangular window already gives E_eta>=1619 q_eta/1000
uniformly on its declared volumes. This is a translation of that earlier
bound, not a new estimate or wider window. Scaling h by c>0 scales E and
Delta by c and leaves q and the dimensionless cut unchanged.

## 5. Actual plaquette calculation, not an abstract two-level example

Let Wp be an ordinary plaquette with four distinct link variables. Each
link appears in its fundamental quaternion representation. On each link,
-Delta_e Wp=3Wp and |nabla_e Wp|²=1-Wp². Consequently, at arbitrary V,

\[
 \mathcal E_\nu(W_p)=4\omega_0(1-W_p^2),\qquad
 \mathcal R_p=\frac{4(1-\omega_0(W_p^2))}
                  {\omega_0(W_p^2)-\omega_0(W_p)^2}.    \tag{10}
\]

Use f=Wp-omega_0(Wp) and a sufficiently small eta to obtain balanced
probes; subtracting the mean changes no gradient. R_p is their weak-probe
limit and an upper bound, not a mass-gap lower bound.

At V=0, Omega=1, omega(Wp)=0 and omega(Wp²)=1/4. The source Wp has
exact free electric energy12. Haar multiplication gives Wp the semicircle
density (2/pi)sqrt(1-c²), whose even moments are Catalan_n/4^n. For f=Wp,
t=eta² and s=sqrt(1-t), summing its convergent moment series gives

\[
 \boxed{E_\eta=1-\frac{2s^2}{1+s},\quad 0\le t<1,}
 \qquad E_\eta=\tfrac34\eta^2+\tfrac18\eta^4+O(\eta^6),
\]
\[
 q_\eta=\tfrac1{16}\eta^2+\tfrac9{1024}\eta^4+O(\eta^6),
 \qquad E_\eta/q_\eta\longrightarrow12.                 \tag{11}
\]

Derivation: sum_n omega(c^(2n))t^n=2/(1+sqrt(1-t)); insert this into
t omega((1-c²)/(1-t c²)). Also 0<=E_eta<=eta² pointwise. The limit
eta→1 of this expression is1; no endpoint smoothness claim is used.
This calculation proves a local physical excitation cost at the free
point; it neither identifies the lowest excitation in every volume nor
controls all neutral loop combinations at interacting coupling.

## 6. Lineage, analytic obligations and the next actual YM task

Inputs are YC22's exact true-ground form, YC20's cut-crossing source and
YC21's projected-return metric/curvature distinction. The nonselective
measurement update and measurement energy are classical quantum results;
see Yi, Talkner and Kim, *Single-Temperature Quantum Engine Without
Feedback Control*, equations5–9, https://arxiv.org/abs/1703.04359.
The metric/curvature split belongs to established state geometry; see
Provost and Vallee, *Riemannian structure on manifolds of quantum states*,
Commun. Math. Phys.76(1980),289–301,
https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-76/issue-3/Riemannian-structure-on-manifolds-of-quantum-states/cmp/1103908308.full .
YC23's role is the typed native-cut adapter, exact YM normalization,
balanced smooth probe and calibrated spectral target, not a claim to have
invented these classical identities.

The general primitive carrier, its state-selection law, connected topology,
analytic completion and first dynamical actualization remain proposals.
They are not prerequisites silently supplied by the scalar letter lambda.
For YM the carrier, measure and electric form are declared and all formulas
above are on their full finite-lattice domains. The next substantive target
is to transport both E and q through an actual spatial block, keeping all
generated neutral readings. Proving a uniform positive ratio along the
physical continuum trajectory remains the difficult step. A single
plaquette, a balanced cut, or nonzero curvature does not settle it.

## Reproduce

Run with Python3.12 from this directory:

    python yc23_cut_agitation.py --write
    python -m unittest test_yc23

The written proofs cover all smooth physical probes. Exact matrix,
polynomial, moment-series and normalization controls exercise the stated
identities; they do not replace those proofs by a finite harmonic sample.
