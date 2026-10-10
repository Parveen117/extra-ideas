# YC35 — excitation-count control of the inverse and its hidden return

10 October 2026. Research owner: Monty Dabas. Continues YC34 and the
target-reference audit at `1db89bc`. Keep the complete YC27 carrier,
internal profiles, admitted volumes and existing physical gap floors.
Frozen predecessor packets are unchanged.

**Result.** The exact inverse response to an at-most-k-factor source has
a volume-independent excitation-count tail in the moving vacuum frame:

\[
 \boxed{\|N^{1/2}\widehat A^{-1}I_k\|\le\mathcal B\sqrt{k},
 \qquad
 \|(I-P_j)\widehat A^{-1}I_k\|
       \le\mathcal B\sqrt{\frac{k}{j+1}}.}                   \tag{1}
\]

Here P_j retains the excited physical states with at most j excited
reference factors, I_k embeds its k-factor source space, and j>=k>=1.
All local harmonics, physical centre sectors and matched charged factors
remain. B is an explicit finite filter integral below, independent of the
number of factors. This is a count of reference excitations, not a
conserved physical particle number.

For the exact bounded Fredholm family F=Ahat(Ahat+tau)^-1, let S_k be its
full hidden-space Schur reduction onto P_k. First compress F to P_j,
then eliminate P_j-P_k, obtaining S_(k,j). This is a genuine reduction
of that bounded pencil, and

\[
 \boxed{0\le S_{k,j}-S_k
       \le\frac{\tau^2\mathcal B^2 k}{j+1}I.}                \tag{2}
\]

It retains the existing reserve m/(m+tau), with source and lift-metric
bounds. The constant is not numerically optimized; no practical j or
radius has been certified. The input P_j F P_j still comes from the
complete interacting operator. Replacing it by a local computable family
is a separate approximation task, explicitly stated in section7.

## 1. Carrier, domains, source, units and zero-reading

Use H(s)=H0-s ell sum_p xi_p Wp on the full tensor carrier of complete
YC27 factors. At a fixed admitted s, let U=U_s be YC31's exact unitary
transport of Omega=product_f Omega_f to the actual vacuum Psi_s. Put

\[
 \mathsf A=U^*(H(s)-E_s)U,\quad
 \mathsf A\Omega=0,\qquad
 \widehat A=\mathsf A|_{Q_0\mathcal H_G},\quad
 Q_0=I-|\Omega\rangle\langle\Omega|,
 \qquad \widehat A\ge mI.                                   \tag{3}
\]

H_G denotes the globally gauge-invariant space. The gap inputs are
m=592/175 or11188/3325 in the respective28/64-link profiles. The full
tensor gap>=1/2 underlies the unitary locality construction; it is not
replaced by the larger physical gap on nonphysical states.

The unbounded domain is the actual unitary pullback of the electric
domain. This stage only takes bounded functional calculus and bounded
orthogonal cuts; it does not need a new assertion about energy-domain
invariance of the count cut. The physical pairing is the original Haar
Hilbert pairing in these unitary coordinates.

Set q_f=I-|Omega_f><Omega_f| and N=sum_f q_f. Define P_j=1_[1,j](N)
on the physical space, and I_k:P_k H_G -> Q0 H_G as inclusion. The vacuum
and every omitted support are retained in the full construction before
these cuts are chosen. The corresponding physical source is U I_k, and
its exact physical response is U Y_k, where Y_k=Ahat^-1 I_k.

The primitive data are the represented operator family, adjoint, product,
pairing, source, complementary cuts and ordered transport history. The
coupling path s is an adapter into them. The variable t below is auxiliary
inverse-electric-energy time; count k, cutoff j and path depth are not
physical time or a spatial renormalization scale.

At ell=0, U=I, Ahat=H0 on Q0 H_G, and H0 commutes with N. The inverse
preserves each exact excitation support. The tail in (1) and the error
in (2) are then exactly zero for j>=k. Internal correlated energies and
source jets have not been set to zero. Further compatibility is recorded
in section8.

## 2. The useful positive local readings

Define the vacuum-preserving evolution and its observable action

\[
 \Gamma_t=e^{it\mathsf A},\qquad
 \beta_t(B)=\Gamma_tB\Gamma_t^*
          =\alpha^{-1}\tau_t^s\alpha(B),\qquad
 \alpha(B)=UBU^*.                                            \tag{4}
\]

The scalar ground-energy phases cancel in the observable action, and
Gamma_t Omega=Omega. Therefore O_(f,t)=beta_t(q_f) is a positive projection
annihilating Omega. Let E_Z be expectation outside Z in the product
reference, as in YC31. Then

\[
 0\le E_Z(O_{f,t})\le I,\qquad
 \langle\Omega_Z,E_Z(O_{f,t})\Omega_Z\rangle=0
 \quad\Longrightarrow\quad E_Z(O_{f,t})\Omega_Z=0.           \tag{5}
\]

The last implication is positivity: the expectation is the squared norm
of the positive square root applied to Omega_Z. Thus these local readings
annihilate their local vacua **before** summing over f. This is stronger
than merely centering an arbitrary shell as in YC34. No creation lift or
count-dependent binomial factor is needed for these positive observables.

If a self-adjoint local z_Z annihilates Omega_Z, put
Q_Z=I-|Omega_Z><Omega_Z|, extended by identity outside. Then

\[
 -\|z_Z\|Q_Z\le z_Z\le\|z_Z\|Q_Z,
 \qquad Q_Z\le\sum_{g\in Z}q_g.                             \tag{6}
\]

Consequently a local family with incident norm
J=sup_g sum_(Z containing g)||z_Z|| satisfies
-JN<=sum_Z z_Z<=JN as quadratic forms. Multiplicities of identical supports
are included. This estimate acts on all vectors, not only fixed-count
inputs; it uses local vacuum annihilation and the actual Hilbert pairing.

## 3. A polynomial-in-time number propagation bound

Use kappa=6750, I=48 or88, and v=8e I ell from YC31. Choose flow exponent16
and write l(r)=min{2,C16(r+1)^(-12)} for its one-factor norm-locality bound,
valid for alpha and alpha^-1. Here C16 is the coefficient in YC31(12)
after removing ||A|| |X|, not its convolution constant C_n. All constants
are uniform along the admitted
coupling path. For r>=3 set u=floor(r/3). Three successive localizations
in (4) give

\[
 \begin{split}
 \|O_{f,t}-E_{B(f,r)}O_{f,t}\|&\le h_t(r),\\
 h_t(r)&=\min\big\{2,\;2[(1+\kappa(2u+1)^3)l(u)\\
 &\hspace{37mm}+\min\{2,4\kappa(u+1)^3e^{-u+v|t|}\}]\big\}.
 \end{split}                                                  \tag{7}
\]

The first error is l(u), from alpha(q_f) to its radius-u approximation.
Physical dynamics then grows that support by u with error bounded by
the second term in (7), using YC31(4). Finally alpha^-1 grows the radius-2u
support by u with error <=kappa(2u+1)^3 l(u). Each intermediate localizing
map is contractive. Replacing the resulting supported approximation by
E_(B(f,r)) of the exact operator costs at most another factor two.

An exponential in |t| cannot be integrated against a merely
faster-than-polynomial filter tail. Instead grow the core radius with t.
Choose

\[
 R_t=6\lceil v|t|\rceil+6.
\]

For r>=R_t, u>=2v|t|, so the physical exponential in (7) is bounded by
exp(-u/2). Define, for r>=3,

\[
 \bar h(r)=\min\{2,\;2[(1+\kappa(2u+1)^3)l(u)
                  +4\kappa(u+1)^3e^{-u/2}]\},\quad
 u=\lfloor r/3\rfloor,
\]
\[
 J_*=\kappa\sum_{r\ge7}(r+1)^3[\bar h(r)+\bar h(r-1)]<\infty,
\qquad
 C_N(t)=\kappa(R_t+1)^3+J_* .                                \tag{8}
\]

Indeed bar h(r)=O(r^-9), making the series summable. Expand O_(f,t) as
the positive core E_(B(f,Rt))O_(f,t), of norm<=1, plus shells
(E_(B(f,r))-E_(B(f,r-1)))O_(f,t) for r>Rt. All these operators annihilate
their local vacua by (5); each shell has norm <=bar h(r)+bar h(r-1).
A factor lies in at most kappa(r+1)^3 balls B(f,r). Applying (6) to
the sum over f and all shells proves

\[
 \boxed{\Gamma_t^*N\Gamma_t\le C_N(t)N,
 \qquad C_N(t)\le\kappa(6v|t|+13)^3+J_* .}                   \tag{9}
\]

The shell argument first gives beta_t(N)<=C_N(t)N; replace t by -t to
obtain the displayed orientation. At finite volume N is bounded and all
expressions are well defined; the displayed constants are independent of
volume. A radius exceeding a finite component simply makes its tail zero.
This proof neither constructs an infinite-volume Hamiltonian nor asserts
physical conservation of factor number.

The auxiliary locality inputs are those of YC31 and sections3.2,4,5.3,6
of [Nachtergaele--Sims--Young, arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2).
Equations(5)-(9) give the particular number estimate used here. Its large
constants have not been numerically optimized. At ell=0 the stronger
identity Gamma_t* N Gamma_t=N holds directly.

## 4. Integrate the complete inverse, keeping its count weight

Use YC32's smooth even function f(omega), zero on |omega|<=m/2 and equal
to1/|omega| for |omega|>=m. Its inverse Fourier weight k_inv is integrable
with every absolute time moment. Fix the convention

\[
 f(\omega)=\int k_{\rm inv}(t)e^{it\omega}dt,\qquad
 f(\mathsf A)|_{Q_0\mathcal H_G}=\widehat A^{-1}.
\]

In particular no inverse is taken on the vacuum. Equation(9) and the
triangle inequality give the stronger weighted mapping estimate

\[
 \boxed{\|N^{1/2}\widehat A^{-1}N^{-1/2}\|\le\mathcal B,
 \qquad\mathcal B=\int|k_{\rm inv}(t)|\sqrt{C_N(t)}dt<\infty.} \tag{10}
\]

The inverse N^-1/2 here is on Q0 H_G, where N>=1. At every finite volume
the weighted integral identity is legitimate; its bound is uniform in
volume. The time-growth power3 in (9) requires only an absolute moment
of order3/2 of k_inv, among the moments already available. This avoids
integrating an uncontrolled exponential time bound.

For input I_k, ||N^1/2 I_k||<=sqrt(k). On I-P_j in the excited physical
space, N>=j+1. These facts prove (1). One can use the sharper of this
tail bound and ||Y_k||<=1/m. Also the existing physical energy estimate is

\[
 \|\widehat A^{1/2}Y_k\|\le1/\sqrt m.                       \tag{11}
\]

The count weight in (10) and energy weight in (11) are different. No
inequality N<=constant Ahat is asserted. Nor does (10) bound the effect
of an arbitrary unbounded error operator on the high-count tail.

A gap alone would not prove (10). For example, on two excited vectors
assign N_L=diag(1,L) and
A=(3/2,1/2;1/2,3/2), with eigenvalues1,2. Then
||N_L^1/2 A^-1 e1||^2=(9+L)/16. The fixed gap coexists with a diverging
number-weighted response. The uniform locality and vacuum-preserving
structure in (5)-(9) are doing essential work.

## 5. A genuine hidden-channel cut and its error

On Q0 H_G use the exact positive bounded pencil

\[
 F=\widehat A(\widehat A+\tau)^{-1},\qquad
 \delta I\le F\le I,\quad\delta={m\over m+\tau}.
\]

Its full Schur reduction and minimizing lift onto P_k are

\[
 T_k=I_k^*\widehat A^{-1}I_k,\quad S_k=(I+\tau T_k)^{-1},
 \quad W_k=F^{-1}I_kS_k
 =I_k+\tau(I-P_k)Y_kS_k.                                    \tag{12}
\]

The identities are on the excited physical carrier; I-P_k denotes its
full hidden complement. They give P_k W_k=I_k, F W_k=I_k S_k,
and I<=W_k*W_k<=delta^-1 I. All source channels and their pairing remain.

For j>=k, let F_j=P_j F P_j restricted to P_j H_G. It obeys the same
bounds delta I<=F_j<=I. Define S_(k,j) as its exact Schur complement onto
P_k, eliminating P_j-P_k, and embed its minimizing lift W_(k,j) into the
full space by zero outside P_j. The variational characterization gives

\[
 \delta I\le S_k\le S_{k,j+1}\le S_{k,j}\le I.                \tag{13}
\]

Let L_j=(I-P_j)W_k. From (1) and (12),

\[
 \|L_j\|\le\tau\eta_{k,j},\qquad
 \eta_{k,j}=\mathcal B\sqrt{k/(j+1)}.
\]

The vector P_j W_k is an admissible trial lift for the j-cut. Its retained
reading is unchanged. Exact stationarity F W_k=I_k S_k makes both cross
terms with L_j vanish, giving

\[
 (P_jW_k)^*F(P_jW_k)=S_k+L_j^*F L_j.
\]

Minimizing within P_j and using F<=I proves (2). Thus an amplitude tail
of order j^-1/2 gives an operator error of order j^-1. This is a positive
variational error, not the difference of two independent spectral bounds.
For requested return error epsilon, any j>=k with
j+1>=tau^2 B^2 k/epsilon suffices formally, independently of total volume.
The displayed constants may require an impractically large j.

For the actual optimizing lift, D_j=W_(k,j)-W_k has zero retained component.
The same stationarity yields the exact identity

\[
 \boxed{S_{k,j}-S_k=D_j^*F D_j.}                             \tag{14}
\]

Consequently the complete lift, retained source and physical Gram errors
obey

\[
 \begin{split}
 \|D_j\|&\le{\tau\eta_{k,j}\over\sqrt\delta}=:w_{k,j},\\
 \|(W_{k,j}^*-W_k^*)b\|&\le w_{k,j}\|b\|,\\
 \|W_{k,j}^*W_{k,j}-W_k^*W_k\|
      &\le2\delta^{-1/2}w_{k,j}.
 \end{split}                                                  \tag{15}
\]

Both lifts have norm<=delta^-1/2, which proves the last line. A general
inhomogeneous full solution additionally has its particular hidden-source
solution. That entire hidden inverse is not approximated here. These are
the bounded-pencil lifts with their actual Hilbert norms, not an automatic
energy-Schur lift certificate for the unbounded Hamiltonian.

## 6. Why cutting before normalizing is a different construction

F_j in section5 is the compression of the **full normalized** operator.
It is generally not A_j(A_j+tau)^-1 for A_j=P_j Ahat P_j. Hidden channels
contribute to the inverse even when their source begins in P_j.

This distinction is why (2) now controls a genuine full hidden-return
comparison while YC34(21), by itself, did not. It does not make F_j
computable from the ordinarily compressed Hamiltonian. Count j also does
not impose a spatial support diameter or truncate local harmonic energies;
the retained dimension still grows with volume and is generally infinite.

In particular, inserting the count tail (1) directly into YC34's unbounded
interaction error on the omitted space would be invalid. That error was
controlled only at fixed count, with constants growing in count. The
bounded variational pencil in (12)-(14) is the reason the present transfer
needs no such uncontrolled multiplication.

## 7. The combined budget for a future local coefficient construction

Suppose a further construction supplies an operator F_tilde_j on P_j with
delta I<=F_tilde_j<=I and full retained-input norm error
||F_tilde_j-F_j||<=epsilon_loc. Let S_tilde_(k,j) be its exact reduction
onto P_k. The minimizing-lift comparison, using lift norms<=delta^-1/2,
gives

\[
 \boxed{\|\widetilde S_{k,j}-S_k\|
   \le{\tau^2\mathcal B^2 k\over j+1}
                      +{\epsilon_{\rm loc}\over\delta}.}    \tag{16}
\]

The two j-cut lifts differ in norm by at most epsilon_loc/delta^(3/2):
subtract their hidden stationarity equations, invert the hidden block
with norm<=1/delta, and apply the source lift bound. Adding (15) gives
a full lift/source error <=w_(k,j)+epsilon_loc/delta^(3/2); the Gram error
is at most2 delta^-1/2 times that sum.

These local-approximation hypotheses have not been established here for
the complete P_j carrier. YC33 supplies a specific one-factor coefficient
construction; it cannot be silently substituted for epsilon_loc at
arbitrary j. YC34 supplies an unbounded-interaction bound on fixed-count
inputs; translating it into this bounded, reserve-preserving coefficient
estimate remains work. Equation(16) tells that work its precise target.

## 8. Zero compatibility, calibration and verification

At ell=0 every h_f preserves its ground and excited spaces. Gamma_t and
H0^-1 preserve N, U=I, and Y_k has support in P_k. Consequently L_j=0,
W_k=I_k, S_(k,j)=S_k, and all errors in (2),(15) vanish exactly for j>=k.
On each physical support J, S_k is H_J/(H_J+tau) with H_J=sum_(f in J)h_f,
as in YC34. This proves typed zero compatibility for source, inverse,
count cut, reduction, lift and pairing. Use this exact zero statement
rather than the loose positive-coupling upper envelope C_N(t).

The constants are selected from the actual gap m, the spectral-flow filter,
the physical propagation bound and the inverse filter. They are not fitted
to a desired mass. B has inverse-energy units. If A_phys=c Ahat and
tau_phys=c tau, then Y_phys=Y/c, B_phys=B/c and the dimensionless return
budget tau^2 B^2 k/(j+1) is unchanged. This consistency does not determine
c or construct a continuum trajectory.

Written proofs establish the locality, positivity, count estimate and
variational transfer on the full family. A focused exact script checks
positive localization with a fixed vacuum, the gap-alone counterexample,
the stationary-lift identity, monotone Schur error, the square-tail bound
and the noncommutation of compression with normalization on finite
matrices. These controls are not Yang–Mills harmonic truncations or
numerical certificates for B,j,R. No predecessor suites were rerun.

**Next:** construct the required bounded local/count-restricted coefficients
with reserve and full-input error epsilon_loc; combine (16) with spatial
and box budgets. A repeated spatial blocking law and physical continuum
gap remain open. YC27's volume-uniform lattice windows are unchanged.
