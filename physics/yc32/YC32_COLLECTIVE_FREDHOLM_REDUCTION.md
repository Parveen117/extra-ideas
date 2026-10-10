# YC32 — collective source control and a local Fredholm reduction

10 October 2026. Research owner: Monty Dabas. Continues YC31 at `5af1305`.
Keep YC27's full electric carrier, internal profiles, torus restrictions
and all-sector physical gap inputs. Frozen predecessors are unchanged.

**Clarification and result.** The volume-uniform lattice gap was already
established in the declared YC27 windows: m=592/175 at mixed cap1/3200,
or m=11188/3325 at mixed cap1/5700, in unit-S³ lattice electric units.
YC31's remaining volume issue concerned
the collective approximation error, not a missing proof of those gaps.
Physical continuum refinement remains a distinct task.

The centered Fredholm operator remains the selected candidate. Here its
exact reduction onto a spatially repeated source space is

\[
 \boxed{S=(I+\tau T)^{-1},\qquad T=V^*A^{-1}V.}       \tag{1}
\]

A is the actual physical excitation operator, and V is a constructed
isometry, including every neutral harmonic on each individual factor.
The hidden space still contains all other physical excitations. It is
eliminated exactly in (1), not discarded.

The new result is a volume-independent operator-norm approximation on
this whole retained source space. Local source tails combine in an l2
norm using product-vacuum orthogonality and overlap counts. There is no
factor equal to the number of blocks. The reduced Fredholm operator then
has positive finite-range approximants with an explicit error budget.
Exact transported sources have zero vacuum leakage. No new parent gap,
closed coarse Yang–Mills interaction, or continuum limit is asserted.

## 1. Construct the retained carrier before taking its zero reading

Use the full family H(s), actual ground Psi=U Omega and exact unitary U
of YC31, at a fixed point of its admitted coupling path. Let Q remove
Psi in the full gauge-invariant Hilbert space, and let
A=(H-E)|Q. Thus A>=m. The unbounded electric domains and all physical
centre sectors remain in the parent operator.

For each correlated reference factor f, define
E_f = H_f^(factor-gauge-invariant) orthogonal to Omega_f. This space may
be infinite dimensional; an acyclic factor can have E_f={0}. Put

\[
 \mathcal E=\bigoplus_f E_f,\qquad
 C_f(\phi)=|\phi\rangle\langle\Omega_f|\otimes I_{f^c},
 \quad
 V_0(\{\phi_f\})=\sum_f C_f(\phi_f)\Omega,
 \qquad V=UV_0 .                                   \tag{2}
\]

Different summands in V0 have exactly one excited factor in different
positions. Hence

\[
 V_0^*V_0=V^*V=I_{\mathcal E},\qquad QV=V.           \tag{3}
\]

Each creator is bounded with norm ||phi|| and commutes with the global
gauge action. E is the retained one-factor neutral source space, not the
entire physical excitation space. Multi-factor neutral excitations,
matched charged factors, winding sectors and all other states remain
in the full complementary cut Q-P, where P=VV*. No finite harmonic
truncation or assumption that a lowest excitation lies in E is used.

The source metric is exactly I without a condition-number estimate.
This follows from the native source choice and its retained unitary
history, not from flattening YC19's older nonunitary quotient metric.
The scalar mixed coupling is the same scoped adapter as in YC31;
neither its path nor the localization radius defines primitive lambda-space
or physical time. Composition, pairing and the complement remain explicit.

## 2. Pull the complete inverse source back to the product reference

Choose a smooth even frequency function f(omega), zero for |omega|<=m/2
and equal to 1/|omega| for |omega|>=m. Use the cutoff constructed in
YC31(15). Its inverse Fourier transform k is in L1, has every absolute
time moment, and obeys YC31(16) with this f. Write

\[
 \mathcal R(B)=\int k(t)\tau_t(B)dt,\quad
 L_k=\int|k|,\quad
 \alpha(B)=UBU^*,\quad
 \mathcal D_f(\phi)=\alpha^{-1}\mathcal R\alpha(C_f(\phi)).       \tag{4}
\]

All integrals are strong operator integrals. On a physical source,
f(H-E)=A^-1 Q. Consequently the complete inverse response Y=A^-1 V
has columns

\[
 U^*Y\phi_f=\mathcal D_f(\phi_f)\Omega,
 \quad \langle\Omega,\mathcal D_f(\phi_f)\Omega\rangle=0,
 \quad \|Y\|\le M:=1/m .                           \tag{5}
\]

The source is transformed by alpha before filtering, and transformed
back afterwards. Filtering a reference creator in the actual Hamiltonian
without these two maps would source a different vector.

Here is one explicit uniform bound for the columns in (5). Let kappa=6750
be YC31's factor-graph ball bound, and choose its spectral-flow exponent
n=16. Let
l(r)=min{2,C16(r+1)^(-12)}, with C16 the constant in YC31(12) after removing
||B|| |support B|. For ell>0 and r>=3, put u=floor(r/3) and

\[
 h(r)=\min\left\{2L_k,\;
 2\left[L_k\{1+\kappa(2u+1)^3\}l(u)
       +\epsilon_k(u,\kappa(u+1)^3)\right]\right\},              \tag{6}
\]

where epsilon_k is YC31(5); set h(0)=h(1)=h(2)=2L_k. Then

\[
 \|\mathcal D_f(\phi)-E_{B(f,r)}\mathcal D_f(\phi)\|
       \le h(r)\|\phi\| .                          \tag{7}
\]

To prove this, first localize alpha(C_f) in B(f,u), apply R and localize
in B(f,2u), then apply alpha^-1 and localize in B(f,3u). The three errors
are bounded by L_k l(u), epsilon_k(u,kappa(u+1)^3), and
L_k kappa(2u+1)^3 l(u), times ||phi||. Replacing that supported approximant
by E_{B(f,r)} of the exact operator costs at most another factor two.
This proves (6)-(7). With, for example, the twelfth derivative tail for k,
h(r)=O((r+1)^(-9)). Higher flow exponents give any desired inverse power.
At ell=0 the three maps preserve the single factor exactly; set h(r)=0
there rather than using the loose positive-coupling bound.

## 3. Product-vacuum shells remove the extensive source count

Let L_Z=I_Z tensor |Omega_(Z^c)><Omega_(Z^c)| be the projection onto states
whose excitations are confined to Z. The product-state expectation satisfies

\[
 E_Z(D)\Omega=L_ZD\Omega.                            \tag{8}
\]

For each column use d_(f,0)=L_{B(f,0)} D_f(phi) Omega and, for r>=1,
d_(f,r)=(L_{B(f,r)}-L_{B(f,r-1)}) D_f(phi) Omega. These are local vectors,
linear in phi, and orthogonal to Omega. Their bounds are

\[
 \|d_{f,0}\|\le L_k\|\phi\|,\qquad
 \|d_{f,r}\|\le h(r-1)\|\phi\|\quad(r\ge1).          \tag{9}
\]

The second inequality uses nested *orthogonal vector projections*, so
no sum of two operator-localization errors is necessary. If d(f,g)>2r,
the r-shell vectors based at f and g are supported on disjoint factors
and both have zero vacuum component. Their inner product is exactly zero.
There are at most N_r=kappa(2r+1)^3 possibly overlapping centres per f.
For a coefficient vector phi in E, Cauchy–Schwarz and the adjacency row
bound give

\[
 \left\|\sum_f d_{f,r}(\phi_f)\right\|^2
 \le N_r h(r-1)^2\sum_f\|\phi_f\|^2.                \tag{10}
\]

This also holds when every E_f is infinite dimensional: it uses the norm
of the local linear map, not a count of its basis vectors.

Define the truncated synthesis map and collective tail

\[
 Y_R\phi=U\sum_f L_{B(f,R)}\mathcal D_f(\phi_f)\Omega,
 \qquad
 e_R=\sum_{r>R}\sqrt\kappa(2r+1)^{3/2}h(r-1).        \tag{11}
\]

Summing (10) in norm establishes

\[
 \boxed{\|Y-Y_R\|_{\mathcal E\to Q\mathcal H}\le e_R,
       \qquad QY_R=Y_R.}                           \tag{12}
\]

With the concrete exponent in (6), e_R=O((R+1)^(-13/2)); stronger polynomial
rates follow by increasing that exponent. All constants are independent
of lattice volume and of the dimension of each neutral factor space.
This closes YC31's collective-source problem for the declared E. It is
not a norm approximation of A^-1 on arbitrary physical vectors.

The exact U is retained at the last step of (11). Thus vacuum leakage is
zero, rather than merely small. Replacing U itself by a strictly local
rounding would be a different approximation requiring a new leakage bound.

## 4. The inverse compression is an exact spatially reduced operator

Set T=V*Y and T_R=V*Y_R. Equations(3),(12) imply

\[
 0\preceq T\preceq M I,\quad \ker T=\{0\},\qquad
 \|T-T_R\|\le e_R.                                 \tag{13}
\]

In fact T_R is precisely the block-distance truncation of T. To see this,
pair a column in (11) with a one-factor excitation at g. If g is in
B(f,R), L_{B(f,R)} fixes the bra; if not, the pairing vanishes. Therefore

\[
 (T_R)_{gf}=\begin{cases}T_{gf}&d(g,f)\le R,\\0&d(g,f)>R.
              \end{cases}                          \tag{14}
\]

In particular T_R is self-adjoint and has range R on the retained factor
graph. It need not be positive. Use the explicit positive correction

\[
 T_R^+=T_R+e_R I,
 \qquad 0\preceq T\preceq T_R^+\preceq T+2e_RI
                              \preceq(M+2e_R)I.     \tag{15}
\]

This retains the range R. It does not claim that entrywise truncation
preserves positivity or that the correction is an actual YM coupling.

For F=A(A+tau)^-1=I-K_tau on the *full* excitation carrier,
F^-1=I+tau A^-1. The bounded Schur inverse identity on P=VV* now gives

\[
 S_P(F)\ \hbox{in E coordinates}
      =(V^*F^{-1}V)^{-1}=(I+\tau T)^{-1}=S.          \tag{16}
\]

This is not V*F V: it includes the full hidden return. Put
S_R=(I+tau T_R^+)^-1. Inverse order and the resolvent identity prove

\[
 \begin{split}
 0&\preceq S-S_R\preceq 2\tau e_R I,\\
 \delta_R I&\preceq S_R\preceq I,
 \quad \delta_R=\frac1{1+\tau(M+2e_R)} .             \tag{17}
 \end{split}
\]

Thus I-S_R is a positive contraction with a volume-independent reserve.
As R grows the reserve approaches YC30's existing delta=m/(m+tau).
The finite-volume exact T is compact because A^-1 is compact; S is
identity minus a compact operator and invertible Fredholm of index zero.
Adding e_R I need not preserve compactness on infinite local channel
spaces. S_R is still invertible and hence Fredholm of index zero; no
ordinary determinant or trace-class assertion is being made.

## 5. Preserve positivity while making the reduced map finite range

The inverse in S_R need not have finite range even though T_R^+ does.
It has a controlled finite-range polynomial approximation. Define

\[
 b_R=1+\tau(M+2e_R),\quad
 Z_R=I-\frac{I+\tau T_R^+}{b_R},\quad
 q_R=1-b_R^{-1},\quad
 S_{R,p}=\frac1{b_R}\sum_{j=0}^p Z_R^j.              \tag{18}
\]

Here 0<=Z_R<=q_R I. Each term has range at most jR, so S_(R,p) has range
pR. Summing the geometric tail, without any lattice-volume factor, gives

\[
 \boxed{\delta_R I\preceq S_{R,p}\preceq S_R\preceq S\preceq I,
 \qquad\|S-S_{R,p}\|\le2\tau e_R+q_R^{p+1}.}         \tag{19}
\]

For illustration choose tau=4 and any R with e_R<=1/800; such an R exists
uniformly in volume by (11). Use the conservative correction 1/800 if the
exact tail is smaller. Degree p=8 then gives the following exact bounds:

| Old factors | q_R upper bound | delta_R lower bound | Operator error |
|---|---:|---:|---:|
| 28 links | 1103/2028 | 925/2028 | <3/200 |
| 64 links | 335297/614997 | 279700/614997 | <3/200 |

Both ninth powers are <1/200 and 2tau e_R<=1/100. These are exact
conditional arithmetic bounds, not a numerically selected radius: YC31's
coarse filter constants have not been optimized to give a practical R.
The finite range is in the retained factor index; local neutral spaces
are still complete. The local coefficients themselves are constructed
from the exact transported inverse response, not yet from an autonomous
small-box numerical solve.

## 6. Retain the lift, source and physical pairing

Let Qh=Q-P denote the full hidden excitation space. The exact bounded
Fredholm lift and a constraint-preserving approximation are

\[
 W_F=V+\tau Q_hYS,\qquad
 W_{F,R}=V+\tau Q_hY_RS_R .                          \tag{20}
\]

Both obey V*W=I and QW=W exactly. From (12),(17), ||Y||<=M and ||S_R||<=1,

\[
 \|W_F-W_{F,R}\|\le w_R:=\tau e_R(1+2\tau M).
 \quad
 \|W_F^*W_F-W_{F,R}^*W_{F,R}\|
 \le2\delta^{-1/2}w_R+w_R^2 .                       \tag{21}
\]

Here ||W_F||<=delta^-1/2 is YC30's direct-cut bound. For a complete input
source b in QH, the retained source error is
||W_F*b-W_(F,R)*b||<=w_R||b||. The separate particular hidden solution
(Qh F Qh)^-1 Qh b remains in the full inhomogeneous solution; (21) does
not provide an approximation of that whole hidden inverse.

For the physical energy pairing, define A_eff=T^-1 on its natural dense
domain and use its low-energy projection J_L=1_[1/L,M](T), with L>=m.
On this subspace T_L>=1/L and

\[
 W_A=YT^{-1}J_L=VJ_L+Q_hYJ_LT_L^{-1},\quad
 W_A^*AW_A=J_LT^{-1}J_L,\quad
 J_L\preceq W_A^*W_A\preceq(L/m)J_L.                \tag{22}
\]

The exact lift lies in Dom(A): AY=V and T^-1 is bounded on this window.
It is the minimum-energy extension with the specified retained reading.
Let T_(R,L)^+=J_L T_R^+ J_L on ran J_L and set
W_(A,R)=VJ_L+Qh Y_R J_L (T_(R,L)^+)^-1. Its retained reading and vacuum
orthogonality are exact. Since T_(R,L)^+>=T_L>=1/L,

\[
 \begin{split}
 \|W_A-W_{A,R}\|&\le v_R:=L e_R(1+2L/m),\\
 \|W_A^*W_A-W_{A,R}^*W_{A,R}\|
       &\le2\sqrt{L/m}\,v_R+v_R^2 .                \tag{23}
 \end{split}
\]

This controls the actual Hilbert pairing of the physical lift. The
approximate W_(A,R) has not been shown to preserve Dom(A), so (23) is
not an energy-residual or variational certificate for that approximation.
The spectral window J_L itself need not be spatially local.

## 7. Zero recovery, direct verification and the next decision

At zero mixed coupling U=I and H=H0=sum h_f. The one-factor carrier is
invariant and

\[
 Y_0\phi_f=h_f^{-1}\phi_f\otimes\Omega_{f^c},\qquad
 T_0=\bigoplus_f(h_f|E_f)^{-1},\qquad
 S_0=\bigoplus_f h_f(h_f+\tau)^{-1}|E_f .             \tag{24}
\]

Localization then has zero error already at R=0. The exact construction
and its typed zero specialization agree; internal dynamics are retained.
The finite polynomial specialization is the same polynomial of (24),
with its declared geometric error, not falsely equal to the exact inverse
at a finite degree. YC31's nonzero first mixed source jet remains retained.

```sh
python physics/yc32/yc32_collective_return.py
```

Direct exact checks cover full inverse compression versus the Schur
complement, the distinction from simple compression, positive correction,
the source/lift identities, and the table's rational error budgets. The
all-volume synthesis bound is proved by (8)-(12), not inferred from a
finite lattice calculation. No generic tests, harmonic diagonalization,
or predecessor reruns were performed.

The classical names remain useful: operator Schur/Feshbach reduction and
quasi-local spectral filtering. Background sources are
[Dusson–Sigal–Stamm, arXiv:2105.02058](https://arxiv.org/abs/2105.02058)
and [Nachtergaele–Sims–Young, arXiv:1810.02428](https://arxiv.org/abs/1810.02428).
The source isometry, collective shell estimate and approximation budget
above are written for the actual declared YM carrier.

**Next task:** obtain the band coefficients from controlled local-box
data rather than the exact global U/response, and compare the resulting
effective interaction with an admissible next-scale family. In particular,
inverse compression has not proved that the coarse operator has the YM
form or that successive spatial budgets contract. No finer-spacing energy
conversion or continuum trajectory is supplied by (19). The already
proved volume-independent lattice gap and these new approximation bounds
must remain separate from those remaining scale questions.
