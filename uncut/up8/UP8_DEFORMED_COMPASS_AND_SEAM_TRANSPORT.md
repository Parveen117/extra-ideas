# UP8 — a deformed compass, with its seam transport sourced by the centre

9 October 2026. Continues UP4–UP7, TT1, DW1 and LC1; repository base
`d1ac633`. The owner's direction is to build the general theory with a
recoverable flat TVSP case, retaining the nonzero-lambda compass rather
than imposing a closed diagram at every level.

**What is added.** An explicit potential deformation closes UP6's open
parameter-family gate. Its mixed-response defect determines a cut, the
cut determines a positive cycle-return metric, and the existing native
half-share transport of that metric gives

\[
 \boxed{\mathfrak a=\ell R\,d\phi,\qquad
        \mathfrak F=R\,d\ell\wedge d\phi,\qquad
        \Theta(\gamma)=-\oint_\gamma\ell\,d\phi},
 \qquad
 \boxed{\ell=\frac{w^2}{p^2+q^2-w^2}}.                    \tag{1}
\]

Here w is the antisymmetric mixed response, p,q are its symmetric
anisotropy, and phi=arg(p+i q). Thus the turn is computed from the centre
and the declared transport, not supplied as an independent angular speed.
The finite-dimensional geometry is general on its cut sector. It is not
a derivation of a physical time law, mass, or a new Yang–Mills gap.

## 1. Source map and hypotheses

Sources read and preserved: UP1 (centre and degrees), UP4/UP5 (frame and
weight defects), UP6 (cross-corner scale operations), UP7 (all primitive
cuts), ON1/CP1/SS1/DW1 (running ratios), TT1 (cycle as a boost), PH3/NC1
(native half-share connection), QC4/QC5 (return and record memory), LC1
(connection equals the root's lost part), HC1 (typed tower/metric limits).

Let U_lambda(S,V) be a C^5 family of scalar centre potentials on a
positive two-variable patch. Fix the coordinate units and the real
two-dimensional coefficient pairing. T=U_S, P=-U_V. Require the derivatives
used as denominators below to be nonzero. Set

\[
 D_1=S\partial_S,\qquad
 D_2=V\partial_V+mD_1,\qquad
 m=-\frac{VT_V}{ST_S}.                                   \tag{2}
\]

D_1 is the S-at-V scale reading and D_2 the V-at-T reading. There is no
extra rotation postulate in (2). UP6's exact identities give

\[
 [D_1,D_2]=D_1(m)D_1,\qquad
 L=\begin{pmatrix}D_1D_1U&D_2D_1U\\D_1D_2U&D_2D_2U\end{pmatrix},
 \quad 2w=L_{21}-L_{12}=D_1(m)\,ST.                       \tag{3}
\]

Ordinary coordinate Maxwell reciprocity U_SV=U_VS still holds. Equation
(3) compares two **different corner operations**, not two commuting
coordinate derivatives. L is not the thermodynamic Hessian of U in the
(S,V) chart and is not assumed positive or symmetric.

Use the primitive matrices

\[
 R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
 K=\operatorname{diag}(1,-1),\quad J=RK,
 \quad L=aI+pK+qJ+wR.
\]

Write rho^2=p^2+q^2 and Delta=rho^2-w^2. Assume **Delta>0** on the patch.
All statements below are invariant under a common orthogonal rotation of
this fixed carrier. Changing units or the pairing requires transporting
the construction; this is not an untyped invariant of all scalings.

## 2. The open compass and its positive return

**D1 — source to cut.** UP7/TT1 give, now with the source (3),

\[
 \kappa=\frac{L-aI}{\sqrt\Delta},\qquad \kappa^2=I,
 \quad \kappa=\cosh\eta\,k(\phi)+\sinh\eta\,R,
 \quad k(\phi)=\cos\phi K+\sin\phi J,
\]
\[
 \sinh\eta=\frac w{\sqrt\Delta},\qquad
 \tanh\eta=\frac w\rho,\qquad \phi=\arg(p+iq).          \tag{4}
\]

The four-step compass return is

\[
 B=(R\kappa)^2
   =\exp[-2\eta n(\phi)],\qquad
 n=Rk=-\sin\phi K+\cos\phi J.                           \tag{5}
\]

Proof: R k is a symmetric involution, R kappa=cosh(eta)n-sinh(eta)I;
squaring gives cosh(2eta)I-sinh(2eta)n. Hence B is symmetric positive,
det B=1, with eigenvalues exp(±2eta). It is I iff w=0. No new result is
claimed for this fixed-cut part: it is TT1's return with the normalization
and source retained.

There is a contracting direction and an expanding direction for w!=0.
This is **not isotropic inward contraction**. A fixed B repeated is a
fixed-axis boost, not a Euclidean vortex. Moving axes are needed for an
angular return. UP7's earlier word "spiral" should be read with this
distinction from TT1.

The symmetric balanced seam is the pair of rays
`z^T(pK+qJ)z=0`, at angles phi/2±pi/4. The +1 and -1 eigenlines of the
actual oblique cut kappa are different: if delta=atan(sinh eta), their
angles are `(phi+delta)/2` and `(phi+pi-delta)/2`. Thus the general compass
retains both its orientation and its failure of orthogonality. Neither
angle is silently identified with motion in physical space.

## 3. Centre to native connection to return-angle

**D2 — exact seam-transport law.** Choose the existing PH3/NC1 protocol
for the positive return metric B:

\[
 \nabla=d+\mathcal A,\qquad
 \mathcal A=\tfrac12B^{-1}dB,\qquad
 \mathcal F=-\tfrac14[B^{-1}\partial_iB,B^{-1}\partial_jB]
                     \,dx^i\wedge dx^j\quad(i<j).        \tag{6}
\]

The unambiguous component convention is
`F_ij=-[X_i,X_j]/4`, X=B^-1 dB. This is the equal share of two flat
connections, 0 and X; it is a stated transport choice, not the only
possible connection compatible with B. It preserves the metric because
`dB=A^T B+B A`.

Put h=B^(1/2), and express a transported vector in its orthonormal
reading z=h v. The connection becomes

\[
 \mathfrak a=h\mathcal A h^{-1}-dh\,h^{-1}.
\]

For B=exp(t n(phi)), t=-2eta, direct multiplication gives

\[
 \mathfrak a=\tfrac12(\cosh t-1)R\,d\phi
            =\sinh^2\eta\ R\,d\phi=\ell R\,d\phi.       \tag{7}
\]

For example, n'=dn/dphi obeys n n'=-R. Then
`B^-1 dB=n dt+(cosh(t)sinh(t)n'+sinh²(t)R)dphi`.
Conjugation by h and subtraction of dh h^-1 cancel the n and n' terms,
leaving (7). This proves (1), since every coefficient of mathfrak a is
in the same so(2) direction R. Parallel transport is dz=-mathfrak a z.
For a closed loop in the patch, the returned line angle is Theta in (1),
modulo pi for an unoriented line. Ordered compass labels retain more
information than that quotient.

This reuses LC1's law exactly: h=exp(-eta n) loses sinh²eta in a cut
across n, so its lost part is ell. The new handoff is that ell and phi
are both reconstructed from the noncommuting **centre response** (3).
No claim of inventing the underlying connection/holonomy identity is made.

On a surface with coordinates x,y its signed curvature coefficient is

\[
 f_{xy}=\ell_x\phi_y-\ell_y\phi_x,\qquad
 \phi_i=\frac{p q_i-q p_i}{p^2+q^2}.                     \tag{8}
\]

Thus a nonzero defect alone is insufficient: its squared normalized
magnitude and its orientation must vary independently. If phi is fixed,
or if both depend on only one scalar z(x,y), this curvature vanishes.
The signed w must still be kept in the return ledger: ell alone loses
the sign that distinguishes B from B^-1.

For clarity about types, (3) is frame anholonomy. For a connection on a
noncoordinate frame the curvature is
`[nabla_D1,nabla_D2]-nabla_[D1,D2]`. A scalar frame commutator is not by
itself Riemann curvature. Equations (6)–(8) now supply a declared native
matrix connection with its actual curvature; no spacetime identification
is used.

## 4. A genuine lambda-family, flat at zero and curved for every positive value

Use fixed dimensionless S,V and an arbitrary common energy unit. Choose

\[
 \boxed{U_\lambda(S,V)=\frac{S^3}{V}
            +\lambda\frac{S^4}{V}
            +\lambda\frac{S^3}{V^2}},\qquad S,V>0,\ \lambda\ge0.
                                                               \tag{9}
\]

This is a constructive example, not a selected equation of state for a
substance. Every term has a positive-semidefinite ordinary Hessian, and
the first is positive definite. Indeed for S^a V^b the Hessian determinant
is `a b (1-a-b) S^(2a-2) V^(2b-2)`; exponents (3,-1),(4,-1),(3,-2)
give positive, positive, zero determinants with positive first diagonal.
Thus U_lambda is strictly convex for the whole domain. T,P,T_S and the
standard corner Jacobians have their required nonzero signs.

Set x=lambda S, y=lambda/V, h0=1+2x+y, z0=S^3/V. Then

\[
 m=\frac{3+4x+6y}{6h_0},\qquad
 D_1m=-\frac{x(1+4y)}{3h_0^2},\qquad
 \frac w{z_0}=-\frac{x(1+4y)(4x+3y+3)}{6h_0^2}.          \tag{10}
\]

**D3 — closed flat limit.** At lambda=0 the potential is the pure power
S^3/V, all eight UP6 scale operations have constant coefficients in
(log S, log V) and commute. Its cross response is

\[
 L_0=z_0\begin{pmatrix}9&3/2\\3/2&1/4\end{pmatrix},
 \quad w_0=0,\quad \Delta_0=z_0^2\,1369/64,
 \quad B_0=I,\quad\mathfrak a_0=\mathfrak F_0=0.         \tag{11}
\]

Ordinary Hessian positivity in this example must not be confused with
positivity of this cross-corner matrix: L_0 has rank one, and L_lambda
need not have a positive symmetric part. The return metric B is positive
by (5), independently of that distinction.

**D4 — complete positivity certificate for the cut sector and curvature.**
Using L/z0 in (4), the discriminant is

\[
 \Delta/z_0^2=\frac{H(x,y)}{5184h_0^6},
\]

where H has 45 monomials, all with strictly positive coefficients and
H(0,0)=110889. Hence the cut never crosses Delta=0 on x,y>=0. The
numerators of p/z0 and q/z0 also have positive coefficients, so phi has
a single smooth choice on this patch.

The exact curvature coefficient simplifies to

\[
 f_{xy}=\frac{10368x(4y+1)h_0^3(4x+3y+3)^2G(x,y)}{H(x,y)^2},       \tag{12}
\]

where G has 28 positive-coefficient monomials and G(0,0)=945. Complete
exponent/coefficient lists for H and G are in `UP8_RESULT.json`; the
script reconstructs them from the derivatives and checks the exact
factorization. These polynomial signs prove **f_xy>0 for x>0,y>=0**;
this is not a claim extrapolated from a sample grid.

For each fixed lambda>0, dx wedge dy=−lambda²/V² dS wedge dV. Therefore

\[
 \boxed{\mathfrak F_{SV}
       =-\frac{\lambda^2}{V^2} f_{xy}(\lambda S,\lambda/V)R\ne0}
 \quad\text{everywhere on }S,V>0.                         \tag{13}
\]

This supplies a parameter family with an actual flat limit and actual
nonzero native transport curvature, not merely a tilted drawing.

Two rational controls:

| x,y | ell | f_xy |
|---|---|---|
| 1,1 | 40000/45653449 | 76489856000/2084237405595601 |
| 1/2,1/3 | 81/186664 | 483327/4355431112 |

**D5 — first visible orders.** At fixed S,V near lambda=0,

\[
 \frac{w}{z_0}=-\frac{\lambda S}{2}+O(\lambda^2),\quad
 \ell=\frac{16}{1369}\lambda^2 S^2+O(\lambda^3),
\]
\[
 \phi=\arg(35+12i)
       -\frac{652}{1369}\lambda S
       +\frac{420}{1369}\frac\lambda V+O(\lambda^2),
\]
\[
 \boxed{\mathfrak F_{SV}
   =-\frac{13440}{1874161}\frac{S}{V^2}\lambda^3 R
      +O(\lambda^4)}.                                    \tag{14}
\]

The seam direction changes at first order, the lost weight at second,
and the native transport curvature first appears at third. These orders
belong to (9); they are not a universal assignment to every lambda tower.
At S=V=1 the derivative of phi with respect to lambda is
`-232/1369` at zero and `-1852126/137080347` at one.

## 5. Ordered cycles can turn without putting in a rotation

**D6 — two-cycle rotation.** Write B_i=exp(t_i n(phi_i)), with t_i=-2eta_i
from (4). Put c_i=cosh(t_i), s_i=sinh(t_i). The product B_2 B_1 has
scalar and skew coefficients

\[
 M=c_2c_1+s_2s_1\cos(\phi_2-\phi_1),\qquad
 W=s_2s_1\sin(\phi_2-\phi_1).
\]

Its polar rotation is exp(Omega R), with tan(Omega)=W/M and M>0.
This follows by multiplying the primitive cuts and requiring the
remaining factor to be symmetric positive. Equal axes give W=0;
reversing the two cycles changes the sign of W. No independent angular
generator was added.

For the centre (9), at S=V=1 and lambda_1=1, lambda_2=2, the exact
coefficient is

\[
 \tan\Omega=-\frac{273483623347200}{8726504489924666083}\ne0.
\]

This is an **ordered product of two pointwise compass returns**, distinct
from the connection holonomy (1) around a closed path of states. Their
experimental protocols are not interchangeable. The mechanism is the
classical rotation from noncollinear boosts (Thomas–Wigner/polar
decomposition); the native contribution here is its explicit centre-source
and lambda-family adapter. See Simon et al., *Hamilton's Turns for the
Lorentz Group*, <https://arxiv.org/abs/quant-ph/0601060>, and Moretti,
*The interplay of the polar decomposition theorem and the Lorentz group*,
<https://arxiv.org/abs/math-ph/0211047>.

## 6. What the generalization does and does not say

Lambda in (9) is a **potential deformation**. It is not automatically
the polynomial coefficient in L+tau L², the derivative order of the
covariant Jacobian tower, the PH3 sharing fraction, physical time, or the
YM coupling theta. In particular, L+tau L² only rescales L's traceless
part and preserves its normalized cut on a branch where that scalar is
positive (TT1). It cannot generate (10) from a reciprocal fixed L.

The following three meanings of "flat" are kept separate:

| Object | Flat/closed criterion |
|---|---|
| Cross-corner operations | Their frame commutator vanishes |
| One algebraic compass | w=0, hence (R kappa)^2=I |
| Native transport of the return metric | d ell wedge d phi=0 locally |

The example (9) has all three at lambda=0 and nonzero defects in all
three senses for every finite lambda>0. The criteria are not equivalent
for arbitrary fields. For instance, a one-scalar response can have w!=0
but d ell wedge d phi=0. Likewise an identically fixed seam gives no
curvature even if the lost weight changes.

Nor does curvature have to grow monotonically with lambda: for (9), at
S=V=1 its magnitude tends to zero both as lambda->0 and lambda->infinity,
while staying nonzero at every finite positive lambda. A scale crossover
can give a peak; an infinite spike needs an additional degenerating
denominator or unbounded derivative. Neither an irreversible arrow nor
isotropic contraction follows from the determinant-one returns.

This is a generalization of the primitive two-cut response/transport
theory. It explicitly contains the closed TVSP compass. Earlier nonflat
packets (UP4, QC4, NC1, etc.) are not retroactively relabelled as flat.
The four-corner potential identities and ordinary coordinate Maxwell
identities can remain exact while cross-corner transport is curved.
No theorem here embeds every previous physical result in one uniquely
selected lambda-family.

The normalized connection is so(2) on this real two-component carrier;
it is not already the full nonabelian SU(2) Yang–Mills connection. A
larger-carrier/source adapter and the actual gap/physical scale remain
separate obligations. YC10's spectral baseline is unchanged. A concrete
next theory gate is to derive the admissible U_lambda or its flow from a
selected physical source, rather than fitting the deformation after the fact.

## Reproduce and evidence

```sh
python up8_deformed_compass.py --check
python -m unittest test_up8
```

The source-pinned JSON contains every polynomial coefficient used in the
global sign proof, exact point/jet controls, and the ordered-return
witness. Symbolic checks audit the identities; their general proofs and
connection choice are stated above. This is not machine-formal verification.
Thirty-seven exact checks and ten focused tests pass.
