# YC2 — centre response on the actual compact three-link Yang–Mills carrier

9 October 2026. Exact product-Haar moments and finite centre-character
algebra. The executable certificate is `yc2_compact_centre_response.py`; the
focused tests are in `test_yc2.py`.

YC1 constructed the TVSP potential of one SU(2) centre record. YC2 embeds
that record into OL1/CM2's actual one-site compact Yang–Mills carrier and
separates the responses that can see a centre character from those that are
forced to be blind to it.

Let

\[
U_i=(a_i,u_i)\in S^3\simeq SU(2),\qquad i=1,2,3,
\]

and use normalized product Haar measure. Simultaneous conjugation is the
one-site gauge action. The three independent centre flips are

\[
Z_i:U_i\longmapsto-U_i.
\]

The gauge-invariant centre-odd holonomies and centre-even plaquette energies
are

\[
c_i=\frac12\operatorname{Tr}U_i=a_i,
\qquad
X_{ij}=|u_i\times u_j|^2,
\qquad
Q=X_{12}+X_{23}+X_{31}.                                    \tag{1}
\]

OL1's exact commutator identity is

\[
1-W_{ij}=2X_{ij},
\]

and the compact Hamiltonian retained by CM2 is

\[
H_\theta=-\sum_{i=1}^3\Delta_{S^3,i}+2\theta Q.             \tag{2}
\]

## 1. The multi-link source potential

Introduce odd sources \(\kappa_i\), even self-sources \(h_i\), and even
plaquette sources \(j_{ij}\):

\[
\boxed{
\mathcal W_\theta(\boldsymbol\kappa,\mathbf h,\mathbf j)
=\log\int_{SU(2)^3}
\exp\!\left[-2\theta Q+
\sum_i(\kappa_i c_i+h_i c_i^2)+
\sum_{i<j}j_{ij}X_{ij}\right]dU.}                            \tag{3}
\]

This is a positive reciprocal response potential on the compact
configuration carrier. Its Hessian is the covariance matrix of the declared
observables. It is the configuration response of the multiplication
potential, not the quantum heat trace of \(H_\theta\). When \(\theta=0\) and
all cross-link sources vanish, equation (3) factorizes into three copies of
YC1.

The source and symmetry types are load-bearing:

| channel | gauge action | centre parity |
|---|---|---|
| \(c_i\) | invariant under simultaneous conjugation | odd under \(Z_i\), even under the other flips |
| \(c_i^2\) | invariant | even under all flips |
| \(X_{ij}\) | invariant | even under all flips |
| \(Q\) | invariant | even under all flips |

## 2. Centre symmetry block theorem

At the symmetry seam \(\boldsymbol\kappa=0\), the weight in (3) is even
under every \(Z_i\), for every finite \(\theta,\mathbf h,\mathbf j\). Hence

\[
\langle c_i\rangle=0,
\qquad
\operatorname{Cov}(c_i,c_j)=0\quad(i\ne j),                 \tag{4}
\]

and for every centre-even observable \(E\),

\[
\operatorname{Cov}(c_i,E)=0.                                \tag{5}
\]

Therefore the complete response Hessian has the exact form

\[
\boxed{
\nabla^2\mathcal W_\theta\big|_{\boldsymbol\kappa=0}
=\begin{pmatrix}
D_{\rm odd}&0\\
0&H_{\rm even}
\end{pmatrix},
\qquad
D_{\rm odd}=\operatorname{diag}(\nu_1,\nu_2,\nu_3)>0.}       \tag{6}
\]

On the isotropic line, link-permutation symmetry gives

\[
D_{\rm odd}=\nu I_3.                                       \tag{7}
\]

Thus every selected pair of centre-odd holonomies has

\[
\boxed{\chi_{\rm odd}=1}                                   \tag{8}
\]

on the centre-symmetric seam, at every coupling. Equation (8) is stronger
than the one-face result in one direction: the interaction may couple all
three links, but independent centre parities still prohibit mixed odd
response.

This does not mean that the even Wilson channels are uncorrelated. Centre
symmetry only separates the odd and even blocks.

## 3. Exact shared-link response at the Haar point

At \(\theta=0\) and zero sources, the quaternion vectors are independent and

\[
\mathbb E[u_x^{2k_x}u_y^{2k_y}u_z^{2k_z}]
=\frac{\prod_{r\in\{x,y,z\}}(2k_r-1)!!}
{4\cdot6\cdots(4+2K-2)},
\qquad K=k_x+k_y+k_z.                                       \tag{9}
\]

Expanding \(X_{ij}=|u_i|^2|u_j|^2-(u_i\cdot u_j)^2\) gives

\[
\mathbb E X_{ij}=\frac38,
\qquad
\operatorname{Var}(X_{ij})=\frac{13}{192},
\qquad
\operatorname{Cov}(X_{ij},X_{ik})=\frac1{64}.               \tag{10}
\]

The last number is the response memory contributed by the shared link. The
three-channel Hessian is

\[
\boxed{
H_X=\frac1{192}
\begin{pmatrix}
13&3&3\\
3&13&3\\
3&3&13
\end{pmatrix}.}                                             \tag{11}
\]

For any selected pair of Wilson-energy channels, CA2's compass is therefore

\[
\boxed{
\chi_X=1-\left(\frac{3}{13}\right)^2
=\frac{160}{169},
\qquad
1-\chi_X=\frac9{169}.}                                      \tag{12}
\]

The full symmetric mode and the two anisotropy modes have eigenvalues

\[
\lambda_{\rm sym}=\frac{19}{192},
\qquad
\lambda_{\rm aniso}=\frac5{96},                             \tag{13}
\]

so no response direction is missing. The normalized three-channel
determinant is

\[
\frac{\det H_X}{(H_{11}H_{22}H_{33})}
=\frac{1900}{2197}.                                         \tag{14}
\]

For the actual Wilson-core channel \(Q\),

\[
\boxed{\mathbb E Q=\frac98,qquad
\operatorname{Var}(Q)=\frac{19}{64}.}                       \tag{15}
\]

These are exact compact-carrier numbers. They are not fitted TC1 spectral
coefficients.

## 4. TVSP before observation and VTSP after observation

For two selected centre-odd channels, the literal pre-observation chart is

\[
(T,V,S,P)
=(\langle c_i\rangle,\kappa_j,\kappa_i,\langle c_j\rangle).
\]

The observed chart is

\[
\boxed{(V,T,S,P)
=(\kappa_j,\langle c_i\rangle,\kappa_i,\langle c_j\rangle).} \tag{16}
\]

For two Wilson-energy channels it is similarly

\[
(\langle X_\alpha\rangle,j_\beta,j_\alpha,
\langle X_\beta\rangle)
\longmapsto
(j_\beta,\langle X_\alpha\rangle,j_\alpha,
\langle X_\beta\rangle).                                   \tag{17}
\]

Writing the observation permutation as

\[
\Pi=\begin{pmatrix}
0&1&0&0\\1&0&0&0\\0&0&1&0\\0&0&0&1
\end{pmatrix},
\]

an arrow \(A\), source covector, and response metric must all be pushed:

\[
A_{\rm obs}=\Pi A,
\qquad
H_{\rm obs}=\Pi H\Pi^T.                                    \tag{18}
\]

This is the exact rule needed before reassigning Maxwell arrows. A pure
first-two-coordinate swap reflects the arrow across the diagonal: the first
and third quadrants remain themselves, while the second and fourth exchange.
Moving an arrow from the third quadrant to the first requires an additional
orientation reversal, such as a physical conjugate-sign convention; it does
not follow from \(T\leftrightarrow V\) alone. Mixed response can rotate the
arrow, but its sign must be derived from the transported constitutive block.

## 5. Why the Wilson potential cannot determine centre-sector splitting

Every \(X_{ij}\) and therefore the full multiplication potential \(2\theta Q\)
is unchanged by all three centre flips. Any untwisted configuration integral
built only from these even channels has no centre-character label. Thus

\[
\boxed{\text{centre-even Wilson response alone cannot distinguish }
\sigma\in\{\pm1\}^3.}                                      \tag{19}
\]

This is the structural reason CM2 finds the same leading folded core in all
eight sectors. The character resides in how the wavefunction is joined across
the eight centre sheets, not in the local folded potential.

The correct sector-sensitive observer is the character Fourier transform of
twisted heat traces. Let \(Z_z\) be the unitary implementing
\(z\in(\mathbb Z_2)^3\). For a character \(\sigma\),

\[
P_\sigma=\frac18\sum_z\chi_\sigma(z)Z_z,
\qquad
K_z(t)=\operatorname{Tr}_{\rm gauge}
\left(Z_z e^{-tH_\theta}\right),                            \tag{20}
\]

and

\[
\boxed{
Z_\sigma(t)=\operatorname{Tr}_{\rm gauge}
(P_\sigma e^{-tH_\theta})
=\frac18\sum_z\chi_\sigma(z)K_z(t).}                       \tag{21}
\]

Character orthogonality makes the eight projectors exact and mutually
orthogonal. Since the compact operator has discrete spectrum,

\[
E_{0,\sigma}
=-\lim_{t\to\infty}\frac1t\log Z_\sigma(t),                 \tag{22}
\]

and relative to the all-even sector,

\[
\boxed{
E_{0,\sigma}-E_{0,+}
=-\lim_{t\to\infty}\frac1t
\log\frac{Z_\sigma(t)}{Z_+(t)}.}                            \tag{23}
\]

Equation (23) is the new computational target for CM2's missing absolute
sector-splitting rate. It inserts the observer that the even potential lacks.

## 6. What advanced and what remains

YC2 completes YC1's first multi-link obligation at the symmetry level:

- YC1 embeds exactly into each link when the interaction and cross sources
  are removed.
- For the interacting compact carrier, centre parity proves an all-coupling
  odd/even response decomposition.
- The shared-link covariance is computed exactly at the Haar point, giving
  the pure response numbers \(160/169\) and \(9/169\).
- The missing centre-sector observer is identified as the twisted heat-trace
  transform (21), rather than another derivative of the Wilson potential.

No absolute value or asymptotic rate for (23) is computed here. The result is
one-site and SU(2). It does not address nonconstant spatial modes, volume
uniformity, the continuum limit, or the four-dimensional mass gap.

The next useful stage is to bound the seven nontrivial twisted traces relative
to \(K_+(t)\) at weak coupling. CM2's eight-patch folding gives their local
character phases; CB1's positive heat-kernel estimates give the untwisted
return. A successful bound must retain paths that cross between centre sheets,
because paths confined to one folded core cancel under the nontrivial
character sum.

## Sources used unchanged

- YC1 for the one-link centre compass and its typed TVSP source potential.
- CA2 for the two-source response compass and observation chart.
- OL1 for the exact one-site compact operator and plaquette identity.
- CM2 for its eight centre characters, folded core, and open absolute
  splitting rate.

The result JSON pins these four source files by SHA-256.
