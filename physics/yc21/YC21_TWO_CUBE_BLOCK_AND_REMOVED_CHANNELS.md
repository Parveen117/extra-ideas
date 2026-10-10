# YC21 — a coupled two-cube block retains its removed channels

9 October 2026. Continues YC20 at `88e40f7`. Frozen packets are unchanged.
Unit-S³ electric normalization, Wp=Tr(Up)/2, as in YC12–YC20.

**Actual spatial result.** Join two twelve-link cubes by all four bridges
and all four interface faces. This is the full 28-link rectangular block
with 4×2×2 vertices (3×1×1 elementary cells). For internal cube couplings
theta_A,theta_B in [0,1] and four interface couplings eta_p in [0,1/2],

\[
\boxed{\Delta_{B,\mathrm{full}}\ge1/5.}                 \tag{1}
\]

All boundary charges and harmonics are included. At eta_p<=1/4 the stronger
floor 13/50 holds. The bridge-vacuum Schur reduction retains **all cube
states**, and its complete returned metric satisfies, for real z<=1/2,

\[
\boxed{I\le M(z)=-F'(z)\le(50/49)I.}                  \tag{2}
\]

This is a concrete block and energy-dependent metric bound, not just YC20's
conditional adapter. The four-face interaction is not truncated in order.

Tiling periodic lattices by these blocks, with sides Lx divisible by4 and
at least8, Ly,Lz even and at least4, gives a new **anisotropic** uniform result:
the above internal couplings and remaining plaquette couplings
0<=epsilon_p<=epsilon<=1/12800 imply

\[
\boxed{\Delta_{\mathrm{full}}\ge\frac15-\frac{2048}{5}\epsilon
                                  \ge\frac{21}{125}.} \tag{3}
\]

The physical gauge-invariant vacuum sector inherits this floor. This does
not improve YC10's isotropic window. The next-block boundary cost increases;
no iterated RG contraction or continuum mass gap is claimed.

## 1. The owner's dimension-removal reminder

The intended principle is that an unseen direction's contribution remains
in the law read on the retained carrier. DM1 gives the exact native example
n²−r1²−r2²=r3²; its three-cut mode becomes the two-cut turn with g=−c k3.
DM1 does not select a physical k3 or a dimensionful mass. Preserve this
earlier result instead of treating the idea as merely fewer coordinates.

Here we eliminate **bridge Hilbert channels**, not a spacetime coordinate.
Their contribution is retained in an energy return, a norm correction and,
when a parameter-dependent connection is declared, a curvature contribution.
These three objects have different types. In particular, an inverse norm or
an energy Hessian is not itself spacetime curvature. Section5 supplies the
explicit projected-connection identity and its actual block bound.

The primitive order remains cut -> complementary channels -> source and
balance -> response/potential. No scalar potential is postulated at the
centre, and observation lambda is not identified with eta or z.

## 2. Full carrier, parity cut and exact source norm

Let h_A=H_A−E_A and h_B=H_B−E_B be the actual full cube Hamiltonians.
YC12/YC15 prove unique positive grounds and the charged floor

\[
g_* =1279511/4672512>1/4.
\]

On L²(SU(2)^28), with no boundary Gauss quotient imposed, set

\[
H_B=h_A+h_B+\sum_{e=1}^4(-\Delta_e)-\sum_{p=1}^4\eta_pW_p.
\]

Omitted cube and plaquette constants are scalar shifts only. The compact
elliptic operator plus bounded real potential is self-adjoint with compact
resolvent and a unique positive ground. All stated Schur operations are
valid on its closed form domain; their cross blocks are bounded.

Label the four bridge links cyclically. Each interface face contains two
consecutive bridges, with signatures {12},{23},{34},{41}. Let J be simultaneous
centre flip of all four bridges. Every Wp commutes with J. The retained
projector P fixes every bridge to its Haar constant and is identity on both
full cubes. Split its complement into Qe (even J, excluding P) and Qo (odd J).

A link harmonic of degree n has energy n(n+2) and centre parity (−1)^n.
Every nonconstant even four-link state costs at least6: either two odd links
cost3+3, or a nonconstant even link costs at least8. Odd states cost at least3.
Thus, with s=sum_p |eta_p|,

\[
D:=Q_eH_BQ_e\ge(6-s)Q_e,\qquad Q_oH_BQ_o\ge(3-s)Q_o. \tag{4}
\]

There is no coupling between the odd carrier and P or Qe. We do not use6
as a floor for the entire hidden carrier; the odd sector is kept separately.

Write r_p=Qe Wp P, R=sum eta_p r_p. Haar integration gives

\[
PW_pP=0,\quad r_p^\dagger r_p=I/4,\quad
r_p^\dagger r_q=0\ (p\ne q),\quad
\boxed{R^\dagger R=vI,\quad v=\tfrac14\sum_p\eta_p^2.} \tag{5}
\]

For the cross term, different signatures leave an unmatched odd bridge.
For a square, integrating a bridge in a fundamental normalized trace gives
1/4 independently of the cube edge values. This is an operator identity
on all cube states, not a vacuum average. At |eta_p|<=eta, v<=eta², s<=4eta.
The interacting D does mix individual signatures: (5) does **not** discard
cross terms R†(D−z)^−1R or the full YC17 tube return.

## 3. Full returned energy, metric and a genuine gap

For z<d:=6−s the exact even Schur pencil and lift are

\[
F(z)=h_A+h_B-z-R^\dagger(D-z)^{-1}R,
\qquad W_z=\binom{I}{(D-z)^{-1}R}.
\]

The positive sign in the hidden component follows from the off-diagonal
Hamiltonian entry −R. The odd carrier is omitted from W only because it
is exactly decoupled and has the separate floor (4).

\[
0\le\Sigma(z)\le\frac{v}{d-z}I,\qquad
M(z)=W_z^\dagger W_z=I+R^\dagger(D-z)^{-2}R,
\]
\[
\boxed{I\le M(z)\le\left(1+\frac{v}{(d-z)^2}\right)I.} \tag{6}
\]

This returns every hidden component, including all excited cube propagation.
No inverse is replaced by its action on the four source vectors. YC20's
two-energy lift identity and complete-residual enclosure apply directly.
Here the starting Hilbert pairing is I: its metric-only Schur complement
is I, whereas the energy-dependent lifted metric is (6).

Let Omega0 be the product of the two true cube grounds and four bridge
constants. Its energy expectation is zero, so E0(H_B)<=0. The retained
operator h_A+h_B has a simple zero and all other levels at least g_*.
If 0<z<min(d,3−s) and

\[
\boxed{g_*-z-\frac{v}{d-z}>0,}                         \tag{7}
\]

F(z) is positive on the codimension-one space perpendicular to Omega0's
retained component. Full Schur inertia therefore permits at most one
eigenvalue below z, and that eigenvalue is the ground. Hence E1>=z and
Delta>=z. Compact-resolvent form inertia, not a finite-matrix approximation,
is used here. Endpoint estimates increase in s,v, so the entire box follows
from the largest eta.

| eta cap | even hidden floor d | odd floor | z | exact reserve (7) |
|---|---:|---:|---:|---:|
| 1/4 | 5 | 2 | 13/50 | 6019313/9228211200 |
| 1/2 | 4 | 1 | 1/5 | 3572617/443888640 |

At eta<=1/2 and z<=1/2, (6) gives v/(d−z)²<=1/49, proving (2).
The ground energy itself is enclosed by −v/d<=E0<=0: at negative z,
F(z)>=−z−v/(d−z), excluding spectrum below −v/d. In particular E0>=−1/16.
This is a scalar normalization bound, not a physical mass.

The bounded positive square root M(z)^(1/2) and its inverse exist inside
this one block. Its condition ratio is <=50/49, with ||M^−1/2||<=1.
This does not identify tensor products of these energy-dependent metrics
with YC19's global quotient metric or prove a volume-independent global
condition number. Wz is a Schur lift, not an isospectral unitary at every z.

## 4. Joining the actual new blocks

Partition vertices into translates of {0,1,2,3}×{0,1}×{0,1}. Internal links
form precisely the above 28-link blocks, with 12 old cube faces and four
joining faces. All other links are single bridge factors. Every remaining
plaquette has either two blocks and two bridges or four bridges; supports
have four factors. A block's external incidence is

\[
2[4+10+10]=48.                                         \tag{8}
\]

These are the internal-edge counts on the two boundary planes of each
orientation. A remaining bridge meets four plaquettes. Unlike the cube
tiling, some bridges meet three two-bridge faces and one four-bridge face;
the proof uses the upper count four, not the previous 2+2 assumption.

Take each block's actual positive normalized ground as its reference and
use its full floor g=1/5. Local gauge invariance makes a single internal-edge
quaternion have zero ground expectation. Consequently a two-bridge face
source excites both blocks and both bridges, has norm1/2 and complete support
floor 6+2g=32/5. Its inverse-vector norm is at most5/64. A four-bridge source
has exact energy12 and inverse-vector norm1/24. This covers all harmonics.

Using YC10/YC14's same weighted full-creator norm and all subsequent support
sizes (s_min=1), the seed is at most

\[
\max(48\cdot4\cdot5/64,\ 4\cdot4\cdot5/64)\epsilon
=15\epsilon.                                         \tag{9}
\]

With beta<=48epsilon, k=M_source=4, radius r=1/128 and
exp(8r)<16/15, the established complete cluster estimates give

\[
q\le10368\epsilon,\qquad b\le2048\epsilon,
\qquad\|\mathcal T(c)\|_*\le15\epsilon+qr.             \tag{10}
\]

At epsilon=1/12800, q<=81/100, b<=4/25 and the ball reserve is1/3200.
The domain-compatible full creation similarity and relative resolvent
argument of YC9/YC10/YC14 then prove the actual gap g(1−b), giving (3).
The exact positive block ground, not a truncated Schur eigenvector, is the
reference. No inverse spectator energy is dropped and no induced interaction
is subtracted from the original Hamiltonian. This step changes its tensor
grouping and the allowed anisotropic parameter profile.

The boundary budget per new factor doubles from24epsilon to48epsilon,
while the convenient charged floor falls from g_* to1/5. Thus this certificate
**does not contract the boundary budget**. It nevertheless permits complete
internal tube couplings up to1/2, previously constrained to the tiny external
window in YC15. Next requires improving the boundary source or a multiscale
norm, not silently repeating (10) and claiming convergence.

## 5. Where the removed directions enter curvature

Fix theta_A,theta_B and real z below the hidden floor. Vary the four actual
interface couplings. On the fixed full even carrier the Schur graph W=(I,Z)^T,
Z=(D−z)^−1R, is smooth in operator norm. Set

\[
M=W^\dagger W,\quad \Pi=WM^{-1}W^\dagger,\quad
\mathcal A_i=M^{-1}W^\dagger\partial_i W.
\]

This is the metric-compatible connection induced by projecting the ambient
flat derivative onto the graph. It is a **declared coupling-parameter
connection**, not an unproved identification with UP8 or physical YM gauge
curvature. Direct differentiation of M and Pi gives

\[
\boxed{M\mathcal F_{ij}
=(\partial_iW)^\dagger(I-\Pi)\partial_jW
 -(\partial_jW)^\dagger(I-\Pi)\partial_iW.}             \tag{11}
\]

The difference of the two ordered departures into the complementary carrier
is exactly the curvature contribution. This is the projected-connection
(Gauss/Berry) identity, not a new general differential-geometric theorem.
For a general compatible ambient connection there is also the projected
ambient curvature term. Here that ambient term is zero by construction.

Dimension removal alone does not fix (11)'s sign, size or nonvanishing.
A single parameter path has no curvature two-form. Constant eliminated
directions, or cancellation of the two terms, can leave a flat connection
even with a nonzero energy return. Thus DM1's native turn identity is kept,
while this distinct target gets its own explicit intertwiner and pairing.

For this actual block, letting V_i=Qe W_i Qe,

\[
\partial_iZ=(D-z)^{-1}(r_i+V_i Z),\qquad
\|\partial_iZ\|\le\frac{1/2+\sqrt v/(d-z)}{d-z}.
\]

In the box eta_p in[0,1/2], z<=1/2 this is <=9/49. Since M>=I and
I−Pi is an orthogonal projector,

\[
\boxed{\|M^{1/2}\mathcal F_{ij}M^{-1/2}\|
       \le162/2401.}                                  \tag{12}
\]

At eta=0, D0 preserves each individual bridge signature. Hence
r_i†(D0−z)^−2 r_j=0 for i≠j and the curvature (11) is exactly zero there,
for all allowed internal cube couplings. This does not say curvature stays
zero in the coupled block; neither a nonzero value there nor a curvature
lower bound is computed. The gap proof (7) does not assume either.

## 6. Evidence, classical tools and next obligation

There are 34 exact checks and 14 new tests; the 19 YC20 and15 YC15 tests
also pass under Python 3.12.14. The executable packet checks exact rational reserves, actual lattice
incidence, Haar/source contractions, bridge-parity floors and independent
finite graph-curvature and spectral-count controls. These support the
written full-space proofs; they are not formal proof-assistant verification
or a finite harmonic replacement of the actual YM operator.

The classical spectral tool is the Feshbach–Schur map; see Dusson, Sigal and
Stamm, [arXiv:2105.02058](https://arxiv.org/abs/2105.02058). The connection
identity is derived above in a nonorthonormal frame; its rank-one complex
version is the familiar Berry curvature, see
[Berry (1984)](https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/07/berry120.pdf).
YC9/YC10 supply the full-creator contraction proof used without change.

Reproduce from the repository:

```sh
python physics/yc21/yc21_two_cube_block.py --check
python -m unittest discover -s physics/yc21 -p 'test_*.py'
```

Next: quantify the new block's boundary-edge kinetic/source variance, then
test whether its full returned boundary cost can improve under a second
spatial update. The current finite-lattice anisotropic theorem leaves the
physical scale trajectory, repeated contraction and 4D continuum construction
open. Removing a carrier direction is not silently removing physical volume.
