# YC19 — local excitation dressing, its metric and a complete truncation error

9 October 2026. Continues YC18 at `9a63c80`. Frozen packets are unchanged.

**Result.** The exact ground similarity gives a full excitation operator
expressible as a summable family of bounded local terms, each annihilating
the product reference vacuum. The incident interaction norm is <=130|zeta|
through |zeta|<=rho=1/4320. Its Taylor coefficients have bounded range; the
complete post-fourth-order tail is controlled both in local interaction norm
and relative to the full electric/reference energy. A resolvent test transfers
certified information from that truncation to the actual excitation operator.

The similarity is not unitary. Its full positive metric and the induced
excitation-quotient metric are explicit and retained. Finite block regrouping
has a controlled cost, but does not automatically improve the coupling budget
or produce independent block metrics. The actual gap window remains YC15's;
this is not an iterated spatial RG or continuum mass-gap theorem.

## 1. Source, carrier and the operator to be transported

Use YC15/YC18's product of actual correlated cubes and bridge factors on even
tori with sides>=4. The cube couplings remain real in[0,1], with full charged
floor g=1279511/4672512. Each factor has its ground Omega_i; H_0=sum_i h_i
and H_I>=g|I| on the full excitation component I. No electric harmonic is
discarded. Write H(zeta)=H_0+Phi(zeta), Phi=-zeta sum_p v_p W_p, |v_p|<=1.
The local incidence budget is beta<=24|zeta|. Scalars removed in earlier
normalizations may be restored without changing excitation differences.

YC18 supplies the holomorphic full creation collection c(zeta) in
||c||_*=max_i sum_(I containing i)|I| ||c_I||<=r=1/128 on |zeta|<=rho.
Let C=sum chat_I, S=exp(C), and E be the eigenvalue with S Omega as its
intermediate-normalized ground for real physical couplings. Define

\[
 \widetilde H=S^{-1}(H-E)S,\qquad W=S^{-1}\Phi S.
                                                               \tag{1}
\]

Finite-volume S and S^-1 preserve Dom(H_0): the full fixed-point equation
places every c_I in Dom(H_I), and [H_0,C] is bounded at fixed volume,
as established in YC9. The estimates below do not use an extensive bound
on ||S|| or a uniform Hilbert-space condition number for S.

## 2. Subtract the complete vacuum creator locally

For a bounded A supported on factors Y, decompose A Omega_Y into its exact
vacuum/excitation components a_I, I subset Y. Define

\[
 \Gamma_Y(A\Omega_Y)=a_\emptyset I+\sum_{\emptyset\ne I\subset Y}
       |a_I\rangle\langle\Omega_I|\otimes I_{Y\setminus I},
 \qquad \mathcal N_Y(A)=A-\Gamma_Y(A\Omega_Y).             \tag{2}
\]

Gamma is a *creation lift*, not the rank-one operator |A Omega><Omega| on
the whole Y space. It includes the scalar. N_Y(A) is supported on Y and
annihilates Omega_Y. Extending Y with idle factors extends N by the identity.
For any excitation creator uhat_J, creation commutativity gives

\[
 \mathcal N(A)u_J=[A,\widehat u_J]\Omega.                 \tag{3}
\]

The ground equation and [H_0,C]=sum widehat(H_I c_I) imply the exact identity

\[
 \boxed{\widetilde H=H_0+\mathcal N(W).}                 \tag{4}
\]

In particular this transports the **entire operator on excited states**.
It does not replace the excitation problem by W Omega or by a scalar return.
The spectator issue in YC18 is resolved here by the exact similarity and
creation subtraction, not by setting excited spectator energies to zero.

## 3. A summable local interaction family with an explicit norm

Expand W= sum_p exp(-C)Phi_p exp(C) into nested commutators and then operator
words. A nonzero commutator tuple only uses creators I_j meeting the original
four-factor face X_p. Assign each word the support label
Y=X_p union I_1 union ... union I_l. Apply N_Y separately to each word,
retaining its sign and 1/l! coefficient. Call this labelled family Phi_tilde.
Labels may overestimate minimal operator support, which is harmless.

For a fixed nonzero word A, all factors outside X_p have fixed excitation
status when A acts on Omega: only the four factors in X_p can vary. Thus
there are at most16 mutually orthogonal source components. Their norms sum
to at most4||A Omega||, even on infinite-dimensional excited factor spaces.
Consequently ||N_Y(A)||<=5||A||. This factor5 is specific to these words;
it is not claimed for an arbitrary operator on an arbitrarily large support.

Use the norm of a labelled interaction family

\[
 \|\Phi_{\rm loc}\|_{\rm inc}
  =\max_i\sum_{\tau:i\in Y_\tau}\|\Phi_\tau\|.           \tag{5}
\]

YC9's root count now has no inverse to insert. Roots in X cost beta;
roots in a marked creator cost beta|I|, paid by ||c||_*. The 2^l expanded
words and unmarked creator sum <=4r give

\[
 \|\widetilde\Phi\|_{\rm inc}
 \le5\beta(1+2r)e^{8r}
 \le5(24|\zeta|)(65/64)(16/15)
 =\boxed{130|\zeta|}.                                  \tag{6}
\]

Every term kills its supported reference vacuum. These generally
non-self-adjoint terms are not separately positive Hamiltonians. Equation(6)
is a uniform full local-operator bound, not merely a vacuum-source bound.

## 4. Connected finite-range coefficients and the full tail

At fixed volume the above labelled family is a holomorphic Banach-valued
function on the same complex disk as c. Supports are fixed by the index
tuple, independent of zeta. Uniform summability in(6) justifies Cauchy
estimates in the incident norm. Put

\[
 M=130\rho=13/432,\quad x=|\zeta|/\rho<1,
 \quad\widetilde\Phi(\zeta)=\sum_{n\ge1}\zeta^n\Phi_n.
\]

Then ||Phi_n||_inc<=M rho^-n. A nonzero coefficient at order n contains
one face and creator coefficients of total order n-1. By YC18 each creator
coefficient has a connected face witness, and it meets the original face.
Thus its support has a connected witness of at most n faces, at most4n
factors, and diameter at most n in the graph joining factors that share an
external face. The output support itself need not be connected.

\[
 \boxed{\left\|\widetilde\Phi-
       \sum_{n=1}^4\zeta^n\Phi_n\right\|_{\rm inc}
       \le M\frac{x^5}{1-x}.}                           \tag{7}
\]

At half radius the bound is13/6912, and at quarter radius13/331776. There
is also a summable size bound and an exponential range bound:

\[
 \max_i\sum_\tau |Y_\tau|\|\widetilde\Phi_\tau\|
       \mathbf1_{i\in Y_\tau}\le\frac{4Mx}{(1-x)^2},
 \qquad
 \max_i\sum_\tau b^{\operatorname{diam}Y_\tau}
       \|\widetilde\Phi_\tau\|\mathbf1_{i\in Y_\tau}
       \le\frac{Mbx}{1-bx},\quad b\ge1,\ bx<1.          \tag{8}
\]

For these weighted bounds the interaction family is refined by Taylor
order, so its terms are zeta^n Phi_(tau,n); the unrefined word need not
have finite range. Absolute convergence of this refinement follows from
the Cauchy bounds for x<1. For x<=1/4 one may choose b=2, obtaining a
range-weighted norm <=M. This makes the distance decay explicit; it does
not say the operator has exactly finite range. At x=1 the full unrefined
bound(6) survives, whereas these particular Taylor estimates do not.

## 5. Excitation coefficients and an error measured against their energy

Let Q=I-|Omega><Omega|, and take the quotient by the ground line. Its
operator is A=H_0+Q U Q, U=N(W). On nonempty components use
||u||_1=sum_I ||u_I||, retaining all excitation vectors in every component.
Equation(3) is exactly the YC9/YC10 full excitation commutator. It gives

\[
 \|QUQH_0^{-1}\|_{1\to1}
 \le\frac{8\beta}{g}e^{8r},\qquad
 B:=\frac{8(24\rho)}g\frac{16}{15}
       =\frac{5537792}{31987775}<1.                     \tag{9}
\]

Both this bounded relative operator and the local family are holomorphic.
For A_4=H_0+Q sum_(n=1)^4 zeta^n Phi_n Q,

\[
 \boxed{\|(A-A_4)H_0^{-1}\|_{1\to1}
          \le\delta_4(x):=B\frac{x^5}{1-x}.}            \tag{10}
\]

At half radius delta_4=346112/31987775<0.01083. This controls the full
excitation operator, not an error inferred from the ground tail alone.

To compute the coefficients let C_j be the creators from YC18. Define
L0=V, L1=[V,C1],
L2=[V,C2]+[[V,C1],C1]/2, and
L3=[V,C3]+[[V,C1],C2]+[[[V,C1],C1],C1]/6. Then

\[
 \boxed{\Phi_n=-\mathcal N(L_{n-1}),\quad1\le n\le4.}  \tag{11}
\]

The scalar energy and [H_0,C] contributions are included by N, not forgotten.
No C4 is required explicitly in(11); its ground equation provides the
cancellation. Exact finite controls compare(11) with an independently
expanded full similarity, including the energy through fourth order.

**Resolvent transfer.** If a computation or analytic bound gives an inverse
of A4-z with ||H_0(A4-z)^-1||_1<=R4(z), and delta4 R4(z)<1, then

\[
 \boxed{A-z\text{ is invertible},\qquad
 \|H_0(A-z)^{-1}\|_1\le\frac{R4(z)}{1-\delta4 R4(z)}.} \tag{12}
\]

This is the complete relative-error Neumann argument on Dom(H_0).
For example set b4=B sum_(n=1)^4 x^n; for real 0<=z<g(1-b4),
R4(z)<=1/[1-z/g-b4]. At x=1/2 and z=g/2 the transfer product is at most
692224/21604415<1. This is a conservative consistency certificate, not a
new sharper gap. A4 is generally nonnormal; eigenvalue differences or a
Hermitian variational principle for A4 would not establish(12).

## 6. The metric travels with the observer change

For real physical couplings H is self-adjoint. Set Mcal=S^dagger S>0.
Then, with the transported domains,

\[
 \boxed{\widetilde H^\dagger\mathcal M
             =\mathcal M\widetilde H.}                 \tag{13}
\]

Relative to Omega plus its original orthogonal complement, write metric
blocks m00,m0Q,mQ0,MQQ. The full transformed operator is upper triangular
with blocks (0,ell;0,A). Eliminating the metric cross terms gives

\[
 \boxed{G=MQQ-mQ0\,m00^{-1}m0Q>0,\qquad A^\dagger G=GA.} \tag{14}
\]

G is the metric on the quotient by the ground line. Thus the *exact*
excitation operator has the correct real physical spectrum although it
need not be Hermitian in the original component pairing. Neither Mcal nor
G is asserted to factor over blocks, have a volume-uniform condition number,
or coincide with the TVSP cut metric. Equation(13) is for real couplings;
the complex disk is used only for analytic estimates.

This is the typed counterpart of transporting the metric with an observer
reset: the state, operator and pairing change together. It supplies no new
identification of observation lambda with interface coupling or RG time.

## 7. Finite block regrouping and the remaining iteration gate

Partition the factor set into blocks, each containing at most m factors.
Replace each labelled support Y by the set of blocks it meets, preserving
the operator and all its coefficients. The block reference is sum_(i in B)h_i
with product ground and gap>=g. Exact counting gives

\[
 \|\widetilde\Phi\|_{\rm inc,blocks}\le130m|\zeta|,
 \qquad \|\text{post-fourth tail}\|_{\rm inc,blocks}
       \le mM\frac{x^5}{1-x}.                           \tag{15}
\]

No copy of a shared factor is introduced. If larger correlated blocks are
then *absorbed* into the reference, one must prove their spectral/metric
properties and a new interaction budget. The factor m in(15) does not
automatically contract; long interactions and metric correlations survive.
Consequently one cannot reapply a self-adjoint product-block theorem to
these non-self-adjoint local terms while throwing away Mcal or G.

The next measurable target is a block/metric construction that improves
this budget while keeping the excitation resolvent reserve. YC19 closes
the full-operator locality and controlled truncation step; it does not
close that repeated-block contraction or the physical continuum trajectory.
YC15's gap window is unchanged. Proofs use the complete infinite-dimensional
factor spaces; finite controls and exact rational checks are not formal
proof-assistant verification.

Run `python3 -B physics/yc19/yc19_local_excitation_operator.py --check` and
`python3 -B -m unittest discover -s physics/yc19 -p 'test_*.py'`.
