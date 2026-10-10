# YC24 — the actual bridge response tower retains its full metric

10 October 2026. Research owner: Monty Dabas. Continues LS2 at `7648401`,
using the actual 28-link block and lattice joins of YC21/YC22. Frozen
predecessor packets are unchanged. Unit-S³ electric normalization,
normalized Haar measure and Wp=Tr(Up)/2 throughout.

**Result.** A second bridge-parity cut makes the complete hidden operator
bipartite. Its exact return is a positive operator tower in the square of
the interface scale, with full energy and graph-metric remainder bounds.
All bridge harmonics and cube excitations remain in that tower. Resolving
the eight even parity sectors improves the hidden floor from2 to5/2 at
interface couplings up to1 and improves the actual physical block gaps:

| Internal cube coupling cap | Four interface coupling cap | YC22 physical block floor | YC24 physical block floor |
|---|---:|---:|---:|
| theta_A,theta_B<=2 | eta_p<=1 | 17/10 | **7/3** |
| theta_A,theta_B<=1 | eta_p<=1 | 15/8 | **12/5** |

The full, including charged, block floor remains3/4. The existing all-volume
join therefore gives, on the same rectangular tori and coupling profiles,

| Cube cap | Tube cap | External cap | Physical lattice gap lower bound |
|---:|---:|---:|---:|
| 2 | 1 | 1/4480 | **461/225** |
| 1 | 1 | 1/4000 | **6476/3125** |

The second bound was1619/1000 in YC22. These enlarge the certified floors,
not the external coupling windows. At z<=1/2 the complete block graph
metric now obeys I<=M(z)<=8I/7, improving YC22's13I/9 in this same tube
window. This is a genuine spatial block result and its lattice consequence;
repeated scale contraction and the 4D continuum remain open.

## 1. Actual carrier and the extra cut

Take two disjoint twelve-link cubes, joined along their facing squares by
four bridges b1,...,b4. The resulting open rectangle has28 links and16
faces. The complete cube Hamiltonians are
H_c=-sum_(e in c) Delta_e+theta_c(6-sum_(f in c)Wf),
with one theta_c per cube as in the predecessor. Let h_c=H_c-E_c>=0.
No boundary Gauss condition is imposed until the physical-sector count.

Write eta_i=t xi_i, with fixed real |xi_i|<=1 and 0<=t<=1. This ray notation
covers every box |eta_i|<=1 by choosing t=max |eta_i|; at t=0 any xi may be
used. The spectral results and lattice profiles above use nonnegative
couplings. Define

\[
 H_t=h_A+h_B+\sum_{j=1}^4(-\Delta_{b_j})
                    -t\sum_{i=1}^4\xi_iW_i.             \tag{1}
\]

The four face signatures are12,23,34,41. The projector P fixes all bridges
to Haar constants and is identity on the full two-cube Hilbert space.
The simultaneous flip of all bridges splits off the total-odd sector.
The physical gauge-invariant carrier is total-even: this full bridge flip
is the vertex-centre gauge transformation on the left cube. YC22 separately
controls the charged sector; no positive odd-sector floor is invented.

Let Z_j be the centre flip of bridge j. Within the total-even carrier,
introduce the alternating cut J_alt=Z1 Z3. The free/reference H0 commutes
with it, and each W_i anticommutes with it. P is in its plus side.
After removing P, the eight even bridge signatures split as follows:

| Hidden carrier | Signatures | Reference electric floor |
|---|---|---:|
| E, alternating-cut minus | 12,23,34,41 | 6 each |
| B0, alternating-cut plus | empty, with P removed | 8 |
| B4, alternating-cut plus | 1234 | 12 |
| Bd1,Bd2, alternating-cut plus | 13 and24 | 6 each |

An odd link costs at least3; a nonconstant even link costs at least8.
In the empty signature P removes every bridge-constant state, including
all cube excitations with constant bridges. Thus its remaining floor is8.
Each row is an **infinite-dimensional sector**, not a single harmonic.
The cube energies remain nonnegative operators in every row.

Every E signature connects to each of the four B signatures by exactly
one face flip. There are no within-E or within-B interaction entries.
Let P_E,P_B project onto these two groups. Set A=P_E H0 P_E,
B=P_B H0 P_B, C=P_E(sum xi_i W_i)P_B and R=P_E(sum xi_i W_i)P.
In E,B order the entire even hidden operator is

\[
 D_t=\begin{pmatrix}A&-tC\\-tC^\dagger&B\end{pmatrix},
 \qquad R^\dagger R=v_0I,\quad v_0=\tfrac14\sum_i\xi_i^2\le1.
                                                               \tag{2}
\]

The source identity is YC21's full conditional Haar identity. Resolvents
do mix signatures through C. No off-diagonal neutral return is deleted.

## 2. Two eliminations, one exact retained return

For real z<6 write A_z=A-z, B_z=B-z and

\[
 T_z=A_z^{-1/2}CB_z^{-1}C^\dagger A_z^{-1/2}\succeq0,
 \qquad S_z=A_z^{-1/2}R,\qquad q=t^2.
\]

Whenever q||T_z||<1, nested Schur elimination gives

\[
 \boxed{\Sigma_t(z)
    =qS_z^\dagger(I-qT_z)^{-1}S_z,\qquad
 F_t(z)=h_A+h_B-z-\Sigma_t(z).}                          \tag{3}
\]

This is the complete return t²R†(D_t-z)^-1_EE R. The operators A_z,B_z
include all cube states, bridge harmonics and their energies. Their bounded
inverses and the bounded interaction C make the form Schur operations
legitimate on the inherited closed form domains. P has infinite rank;
no finite source span is assumed invariant.

This is the promised actual YM tower, rather than LS2's scalar source
recurrence imported by analogy. Its parameter q is the square of a declared
interface amplitude at fixed cubes and fixed spectral z. It is not the
source lambda in UP8, observation-zero cut, lattice spacing or an RG step.

**Y1 — all-order response and a signed remainder.** At fixed z,

\[
 \Sigma^{[N]}_t=q\sum_{n=0}^{N}q^nS_z^\dagger T_z^nS_z,
\]
\[
 \boxed{\Sigma_t-\Sigma^{[N]}_t
 =q^{N+2}S_z^\dagger T_z^{N+1}(I-qT_z)^{-1}S_z\succeq0.} \tag{4}
\]

For n>=1 the exact derivative is

\[
 \partial_q^n\Sigma_t
 =n!S_z^\dagger T_z^{n-1}(I-qT_z)^{-n-1}S_z\succeq0.    \tag{5}
\]

These follow by the finite geometric identity and functional calculus of
the single positive T_z. A and C need not commute. The complete coefficient
on one fixed ray is positive; an individual mixed four-face word need not
be. YC17's nonzero connected tube return is included, together with repeated
faces and every later order. No finite recurrence or closure of this
infinite-dimensional tower is asserted.

## 3. Sector-resolved bounds for the whole hidden operator

Set a=6-z and, for z<6,

\[
 c(z)=4\left(\frac1{8-z}+\frac1{12-z}+\frac2{6-z}\right),
 \quad b(z)=4\left(\frac1{(8-z)^2}+\frac1{(12-z)^2}
                              +\frac2{(6-z)^2}\right),
\]
\[
 \tau(z)=c(z)/a,\qquad \alpha(q,z)=a-qc(z).             \tag{6}
\]

For a column C_gamma from one B signature into the four E signatures,
every block has norm<=1, so ||C_gamma||<=2. B_z is block diagonal in these
four signatures. Its floors in the table imply the full operator bounds

\[
 CB_z^{-1}C^\dagger\preceq c(z)I,\qquad
 CB_z^{-2}C^\dagger\preceq b(z)I,\qquad
 \|T_z\|\le\tau(z),\quad\|S_z\|^2\le v_0/a.           \tag{7}
\]

This is an upper comparison of complete operators, not a replacement by
eight basis vectors. If alpha>0, then D_t-z is positive by block-form
congruence, and

\[
 \boxed{0\preceq\Sigma_t(z)\preceq\frac{qv_0}{\alpha(q,z)}I.} \tag{8}
\]

At q=1,z=5/2, alpha=193/2926>0. Thus D_t>=5I/2 throughout the tube box.
At z=0, tau=13/36, alpha(1,0)=23/6 and return<=6I/23. At z=1/2,
alpha(1,1/2)=24017/7590 and tau=17728/41745<1. All these bounds are
monotone worst cases in q and z, so they cover the stated intervals.

From(4),(7), with rho=q tau<1,

\[
 0\preceq\Sigma_t-\Sigma_t^{[N]}
 \preceq\frac{qv_0}{a}\frac{\rho^{N+1}}{1-\rho}I.      \tag{9}
\]

In particular the energy tail is controlled at full interface strength,
not only as a formal expansion about zero.

## 4. The complete metric travels with the tower

The exact hidden lift Z=(Z_E,Z_B)^T is

\[
 Z_E=tA_z^{-1/2}(I-qT_z)^{-1}S_z,
 \qquad Z_B=tB_z^{-1}C^\dagger Z_E.
\]

With W=(I,Z_E,Z_B)^T, the original Haar pairing returns as

\[
 \boxed{M_t(z)=W^\dagger W=-\partial_zF_t(z)
                  =I+Z_E^\dagger Z_E+Z_B^\dagger Z_B.} \tag{10}
\]

The sign in Z_B follows from the negative off-diagonal in(2). Define

\[
 L^2(q,z)=\frac{qv_0[1+qb(z)]}{\alpha(q,z)^2}.
\]

By(7),(8), ||Z||²<=L² and I<=M<= (1+L²)I. At q=1,z=0,
L²=189/2116. For all q<=1,z<=1/2,

\[
 \boxed{I\preceq M_t(z)\preceq
 \left(1+\frac{78682276}{576816289}\right)I
                           \preceq\frac87I.}           \tag{11}
\]

This directly bounds the whole lifted norm; it is not obtained by
differentiating an inequality between two self-energies.
It is a one-block graph metric bound, not a volume-independent condition
number for a tensor product or for YC19's global excitation-quotient metric.

For a finite tower use

\[
 Z_{E,N}=tA_z^{-1/2}\sum_{n=0}^N(qT_z)^nS_z,
 \quad Z_{B,N}=tB_z^{-1}C^\dagger Z_{E,N},
 \quad M_N=I+Z_N^\dagger Z_N.
\]

Let r_N=rho^(N+1). Geometric summation gives
||Z-Z_N||<=L r_N and ||Z_N||<=L(1-r_N). Consequently

\[
 \boxed{\|M-M_N\|\le L^2(2r_N-r_N^2).}                \tag{12}
\]

Thus the full omitted norm is retained along with the omitted energy.
At q=1,z=0,v0<=1:

| Retained return terms | Energy-return tail cap | Graph-metric error cap |
|---|---:|---:|
| t² through t⁶ (N=2) | 2197/178848 (<0.01229) | 1401257585/170595237888 (<0.00822) |
| t² through t⁸ (N=3) | 28561/6438528 (<0.00444) | 665891061017/221091428302848 (<0.00302) |

The metric enclosure is two-sided around M_N. Even though(4) is positive,
M-M_N need not be positive: the independent noncommuting finite control
has negative determinant for that difference. Nor is M_N automatically
-partial_z of a truncated Schur pencil. Using the graph norm with(12)
avoids both errors. This is the full metric required by YC20's retained
cut, not a Gaussian curvature or an assumed physical spacetime metric.

## 5. Improved block and physical lattice gaps

On the gauge-invariant bridge-vacuum carrier P, the two cubes are separately
gauge invariant. The retained reference has one zero state; its remaining
spectrum is at least delta=27/5 for theta_A,theta_B<=2, and at least11
for theta_A,theta_B<=1, by YC13. At a positive z below the hidden floor,
equation(8) gives the codimension-one Schur reserve

\[
 \delta-z-\frac{qv_0}{\alpha(q,z)}.                    \tag{13}
\]

At the endpoint q=v0=1 the two choices have strictly positive reserves:

| Cube floor delta | Trial z | alpha(1,z) | Reserve(13) |
|---:|---:|---:|---:|
| 27/5 | 7/3 | 5941/16269 | 29251/89115 |
| 11 | 12/5 | 311/1260 | 7073/1555 |

Since both z<5/2, the entire hidden block is positive there. The compact
elliptic Hamiltonian has a simple positive ground, and its reference
product has expectation0, so E0<=0. Full form Schur inertia permits at
most one physical level below z. Therefore E1>=z and the physical gap
E1-E0>=z. This proves the opening block table with all modes retained.
The charged complement separately has gap>=3/4 above the true ground
by YC22, so the full block floor is still3/4.

For the existing rectangular tiling, Lx is divisible by4 and>=8; Ly,Lz
are even and>=4. Internal cube faces, the four tube faces and external
faces keep exactly the previous partition. YC22's full creator proof
uses the unchanged full factor floor g=3/4, the same true positive
block ground,48 external incidences, source norm caps and full nonlinear
support bound. Nothing in those estimates is worsened by improving the
physical gap. Its force profile uses only the same coupling caps and
degrees, so it remains valid.

The physical reference floor is min(block physical gap,6). Both improved
block floors are below6. Retain YC22's relative return b (which still uses
g, not delta) and apply its physical excitation estimate delta(1-b):

\[
 \boxed{\begin{aligned}
  \theta_c\le2,\ \eta_p\le1,\ \epsilon_p\le1/4480:
   &\quad\Delta_G\ge\frac73\left(1-\frac{64}{525}\right)
                         =\frac{461}{225},\\
  \theta_c\le1,\ \eta_p\le1,\ \epsilon_p\le1/4000:
   &\quad\Delta_G\ge\frac{12}{5}\left(1-\frac{256}{1875}\right)
                         =\frac{6476}{3125}.
 \end{aligned}}                                       \tag{14}
\]

The corresponding full lattice gaps remain461/700 and1619/2500.
Uniformity is over these finite volumes and anisotropic profiles. The
positive local tower does not by itself construct an infinite-volume QFT.

## 6. What this contributes to the framework and next scale

The extra cut has a concrete operational role: a bridge interaction crosses
it once, so a return crosses it an even number of times. Eliminating its
opposite side produces a positive q=t² response operator. Its derivative
tower, source, norm and remainder are computed on the actual YM carrier.
The diagram-to-operator map has therefore survived this transfer without
assuming that the compass's two-mode algebra exhausts the hidden space.

The method reuses Schur/Feshbach elimination, parity decomposition and
positive-resolvent functional calculus. YC13 already has a positive tower
for the one-cube cut; the new content is this eight-sector,28-link adapter,
its full metric-tail control, and the stronger block/lattice bounds.
Reference for the classical spectral tool: Dusson, Sigal and Stamm,
[The Feshbach–Schur map and perturbation theory](https://arxiv.org/abs/2105.02058).
All identities and numerical comparisons needed here are proved above.

The remaining scale issue is concrete: cube incidence/full-floor is24/1;
the rectangle still has48/(3/4)=64. This budget has not contracted. The
geometric ratio13/36 in the **local coupling expansion** is not an RG
contraction factor. A second spatial join must control its complete
connected neutral source, spectators and lifted metric with a smaller
scale-normalized boundary cost. Physical running and the continuum limit
remain additional obligations. No arbitrary count of remaining steps is
assigned, and no dimensionful mass is extracted from these lattice units.

## 7. Reproduce

```sh
python physics/yc24/yc24_parity_response_tower.py --check
python -m unittest discover -s physics/yc24 -p 'test_*.py'
```

The certificate checks the entire finite parity graph, exact comparison
floors, spectral and lattice reserves, and independent finite controls of
the noncommuting inverse, retained metric and signed energy tail. Written
proofs extend these controls to all harmonics of the actual block. The
finite matrices are controls, not harmonic replacements or proof-assistant
verification. Source hashes preserve the exact predecessor boundary.

Verification: **35 exact checks, 15 new tests and 29 predecessor tests**
(YC21:14, YC22:15) pass on Python3.12.14 / SymPy1.14.0. No predecessor
proof, code or evidence packet was modified.
