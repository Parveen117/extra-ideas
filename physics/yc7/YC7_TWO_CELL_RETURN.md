# YC7 — the first spatial join: a two-cell vacuum gap with the full source

9 October 2026. Continues YC6 (`0b3b9d2`) and the signed-compass handoff
(`e7282d5`). Python 3.12, exact rational Haar integrals and inertia.

**Result.** On the actual periodic **2 × 1 × 1 spatial lattice**, with all
six SU(2) links, both local Gauss laws and all three even global centre
characters, the internal vacuum-sector gap obeys

\[
\boxed{\Delta_{\rm vac}^{(2,1,1)}(\theta)\ge\frac65
\qquad(0\le\theta\le2).}                                      \tag{1}
\]

The minimum certified cell reserve is 12431/10000. The first nonconstant
spatial degree of freedom is retained: this is not two copies of the
one-site calculation. In particular the free gap is **6**, versus 8 on the
one-site vacuum carrier. Two centre-odd winding readings can combine into
an allowed centre-even excitation on different sites.

This remains a finite spatial lattice theorem, with continuous Hamiltonian
time. The transverse directions have one cell each. It is neither a
volume-uniform bound nor a four-dimensional continuum construction.

## 1. Carrier, gauge constraints and centre characters

Use two sites x=0,1 (modulo 2), with y,z periodic modulo 1. Label links

| Link | Direction | Source → target |
|---|---|---|
| A0 | x | 0 → 1 |
| A1 | x | 1 → 0 |
| B0, C0 | y, z | 0 → 0 |
| B1, C1 | y, z | 1 → 1 |

A0 and A1 are independent links; A1 is not A0 inverse. With normalized
product Haar measure dμ, the unreduced carrier is L²(SU(2)^6,dμ). At each
vertex impose

\[
U_e\longmapsto g_{s(e)}U_e g_{t(e)}^{-1}.
\]

The three global centre operations may be represented by

\[
Z_x:A_0\mapsto-A_0,\qquad
Z_y:(B_0,B_1)\mapsto(-B_0,-B_1),\qquad
Z_z:(C_0,C_1)\mapsto(-C_0,-C_1).
\]

They preserve every plaquette. A gauge-centre operation at one vertex flips
both A links, so the alternative plane choice for Z_x is equivalent on
physical states. Flipping B0 alone is **not** a symmetry. We restrict to
the gauge-invariant subspace with Z_x=Z_y=Z_z=+1, with no further reflection
or spatial-exchange restriction.

Write W(U)=Tr(U)/2. The six plaquette words are

\[
\begin{split}
&A_0B_1A_0^{-1}B_0^{-1},\quad A_1B_0A_1^{-1}B_1^{-1},\\
&A_0C_1A_0^{-1}C_0^{-1},\quad A_1C_0A_1^{-1}C_1^{-1},\\
&B_0C_0B_0^{-1}C_0^{-1},\quad B_1C_1B_1^{-1}C_1^{-1}.
\end{split}
\]

The operator, in precisely the link-Laplacian convention of CB1/YC6, is

\[
H_\theta=H_0+\theta V,\qquad
H_0=-\sum_{e=1}^{6}\Delta_{S^3,e},\quad
V=\sum_{p=1}^{6}(1-W_p),\qquad0\le V\le12.                   \tag{2}
\]

Each edge is counted once in H0, however many plaquettes share it. No tree
gauge fixing replaces its electric energy by a guessed independent-loop
kinetic operator. The quadratic form on the compact product, restricted to
the stated invariant sector, defines a self-adjoint operator with compact
resolvent. Its positive full-space ground is unique and invariant under
the gauge and centre actions, so it lies in this vacuum sector.

Classical input: the electric-Casimir/plaquette Hamiltonian is the
Kogut–Susskind construction (Phys. Rev. D 11, 395 (1975),
https://doi.org/10.1103/PhysRevD.11.395). We declare our normalization in (2),
rather than importing physical coupling factors from another convention.
The local gauge/loop distinction is also discussed by Mathur and Sreeraj,
https://arxiv.org/abs/1509.04033. No continuum conclusion from either source
is used in this certificate.

## 2. A complete free cut, including the new winding pair

Let b_x=Tr(B_x)/2 and c_x=Tr(C_x)/2. The following seven vectors are
orthonormal in product Haar measure:

\[
\mathcal B=\big(1,\ 4b_0b_1,\ 4c_0c_1,\
4b_0^2-1,\ 4c_0^2-1,\ 4b_1^2-1,\ 4c_1^2-1\big).
\tag{3}
\]

Their electric energies are (0,6,6,8,8,8,8). Let P be their orthogonal
projection and Q=I−P inside the vacuum sector. Then

\[
\boxed{QH_0Q\ge12Q,\qquad QH_\theta Q\ge12Q.}              \tag{4}
\]

**Completeness proof.** On an SU(2) edge, degree n=2j contributes energy
n(n+2)=4j(j+1) and has endpoint representation j at each end. Peter–Weyl
decomposition followed by the vertex singlet constraints is complete.
Below 12 only n=0,1,2 can occur. The centre conditions are n_A0 even,
n_B0+n_B1 even, n_C0+n_C1 even. A vertex-centre gauge constraint additionally
forces n_A0+n_A1 even. Thus neither A link can have n=1 in this sector.
A single n=2 A link has no vertex singlet at either endpoint, and adding
enough other representations to make it invariant costs at least 12.

Among the four loop links, two n=1 readings with the required parity must
be the B pair or the C pair. At each vertex its single loop has a unique
invariant character. A single n=2 loop likewise has one invariant
character. Three n=1 loops cannot satisfy both B and C parities. Together
with the constant, these give exactly (3). The next allowed energy is 12,
for example the product of all four fundamental loop characters. The
certificate independently enumerates the endpoint singlet multiplicities:
one at 0, two at 6, four at 8, none else below 12.

Consequently a cut retaining only independently centre-even cell readings
would miss actual vacuum excitations. The new products do not assert a
positive string tension or a long-string energy law; they are explicit
finite-lattice winding observables with known electric energy.

## 3. Exact full source and interface cross terms

Let M=PVP, r=QVP and B=r* r. In the ordered basis (3),

\[
M=\begin{pmatrix}
11/2&-1/2&-1/2&-1/4&-1/4&-1/4&-1/4\\
-1/2&5&0&-1/2&0&-1/2&0\\
-1/2&0&5&0&-1/2&0&-1/2\\
-1/4&-1/2&0&21/4&1/12&0&0\\
-1/4&0&-1/2&1/12&21/4&0&0\\
-1/4&-1/2&0&0&0&21/4&1/12\\
-1/4&0&-1/2&0&0&1/12&21/4
\end{pmatrix}.                                               \tag{5}
\]

Every entry is an exact six-link Haar integral. For a single S³ link,
odd moments vanish and

\[
\mathbb E\prod_{a=0}^{3}q_a^{2k_a}
=\frac{\prod_a(2k_a-1)!!}{\prod_{m=0}^{\sum_a k_a-1}(4+2m)}.
\]

The polynomial generator is the complete spherical Laplacian. Projecting
each full residual r_i onto its finite harmonic support gives

\[
r_i=\sum_\lambda r_{\lambda,i},\qquad
B_\lambda=(\langle r_{\lambda,i},r_{\lambda,j}\rangle)_{ij}
\succeq0,\qquad \sum_\lambda B_\lambda=B,                     \tag{6}
\]

with λ in {12,14,16,18,20,22,24,26,32}. The generator eigen-equations and
source reconstruction are checked exactly, including the vanishing of
every possible component below 12. This finite source support is **not**
assumed invariant under V; the hidden dynamics remains infinite-dimensional.

Separate V=V_L+V_I into the two local yz plaquettes and the four interface
x plaquettes. Put r_L=QV_LP and r_I=QV_IP. The complete source square is

\[
B=r_L^*r_L+r_I^*r_I+
\underbrace{r_L^*r_I+r_I^*r_L}_{B_{LI}}.                       \tag{7}
\]

This cross matrix is not zero. For example, with zero-based indices in
(3), (B_LI)_{1,3}=1/4 and (B_LI)_{1,4}=1/12. These terms are retained in
both the source Gram and its resolved matrices. Keeping only separate
local/interface norms would not give the same Schur form.

For the constant retained reading, specifically,

\[
\langle V\rangle=11/2,\quad
\operatorname{Var}_\mu V=43/24,\quad
B_{00}=25/24,\quad
(B_{14})_{00}=3/4,\quad(B_{16})_{00}=7/24.                    \tag{8}
\]

The other entries (B_λ)_{00} vanish. The 3/4 part is the interface source
after subtracting its retained winding product: its horizontal adjoint
degree costs 8, and its two fundamental loop degrees cost 3+3, totaling
14. The local residual lies at 16. Their different kinetic denominators
are explicit spatial information.

As an additional check, retaining **only** the constant for perturbation
theory gives the full spectral source weights
\((\lambda,w)=(6,1/2),(8,1/4),(14,3/4),(16,7/24)\). Hence

\[
E_0(\theta)=\frac{11}{2}\theta-\frac{167}{896}\theta^2+R_3(\theta),
\quad
|R_3(\theta)|\le
\frac{12(43/24)\theta^3}{6(6-11\theta/2)}
\quad(0\le\theta<12/11).                                    \tag{9}
\]

For the error, the rank-one complement is at least 6, 0≤E0≤11θ/2,
and the resolvent identity bounds the difference from its free inverse by
12θ/[6(6−11θ/2)]. Equation (9) does not replace the seven-reading
certificate. It prices the complete kinetic source, not only its static
variance. The constant Haar reading is not assumed to be the interacting
vacuum.

## 4. Full hidden return and certified gap

For z<12, write D0=QH0Q, W=QVQ and

\[
\Sigma_\theta(z)=r^*(D_0+\theta W-z)^{-1}r,\quad
S_\theta(z)=K+\theta M-zI-\theta^2\Sigma_\theta(z),
\quad K=\operatorname{diag}(0,6,6,8,8,8,8).
\]

Since 0≤W≤12Q, inverse order gives the full-space sandwich

\[
\boxed{\sum_\lambda\frac{B_\lambda}{\lambda+12\theta-z}
\preceq\Sigma_\theta(z)\preceq
\sum_\lambda\frac{B_\lambda}{\lambda-z}.}                    \tag{10}
\]

No commutativity between D0 and W is needed. Let L and U be the resulting
lower and upper Schur forms, respectively. Their finite inertias bound the
inertia of S; the positive entire hidden block supplies the full-space
Schur congruence. At most j negative directions in L imply E_j≥z. At least
j+1 nonpositive directions in U imply E_j≤z. If an upper lies at or above
12, use the retained Ritz form instead; no hidden inverse above its floor
is used.

`YC7_ENDPOINTS.json` specifies rational thresholds at all 17 nodes θ=k/8,
0≤k≤16. The independent replay uses exact rational elimination for every
sign. Floating roots only propose these thresholds in `explore_yc7.py`.

| θ | Vacuum gap lower | Vacuum gap upper |
|---:|---:|---:|
| 0 | 6 | 6 |
| 1/2 | 5.6989 | 5.7262 |
| 1 | 5.1133 | 5.4684 |
| 3/2 | 3.5078 | 5.4722 |
| 2 | 1.2753 | 5.7570 |

All decimals in this table are exact rationals. For every cell [a,b], the
ordered eigenvalues are nondecreasing because V≥0. Thus

\[
\Delta(\theta)\ge l_1(a)-u_0(b)\ge12431/10000>6/5.
\]

This proves (1) on the entire interval. It does not assume gap monotonicity.
The wide upper/lower separation at θ=2 describes this comparison's error,
not a measured collapse of the spectrum. The first excited spatial symmetry
and its multiplicity at nonzero coupling are not identified here.

## 5. An exact spatial-join law with a changing source

To expose what changes on joining cells, keep this six-link carrier and set

\[
H_{\theta,\eta}=H_0+\theta(V_L+\eta V_I),\qquad0\le\eta\le1.
\]

At η=0 the magnetic terms are the two local yz blocks; the A links retain
their electric energies. The positive ground factors into those two local
grounds and constants on A0,A1. This is not two copies of YC6's three-link
one-site operator. At η=1 the full spatial lattice (2) is recovered.

Put r_η=Q(V_L+ηV_I)P, s=QV_IP, D_η=QH_{θ,η}Q,
W_I=QV_IQ, R_η=(D_η−z)^{-1}, and Σ_η=r_η*R_ηr_η. For an increment h
with 0≤η≤η+h≤1,

\[
\boxed{\begin{split}
\Sigma_{\eta+h}-\Sigma_\eta
={}&h\big(s^*R_{\eta+h}r_\eta+r_\eta^*R_{\eta+h}s\big)
  +h^2s^*R_{\eta+h}s\\
 &-\theta h\,r_\eta^*R_{\eta+h}W_I R_\eta r_\eta.
\end{split}}                                                 \tag{11}
\]

Proof: r_{η+h}=r_η+hs and
R_{η+h}−R_η=−θh R_{η+h}W_I R_η. Expand the new source on both sides of
its inverse. The resulting Schur step is
S_{η+h}−S_η=θh PV_IP−θ²(Σ_{η+h}−Σ_η).
The source-cross terms in the first line are the part that would be missed
by simply reusing YC6's fixed-source coupling identity. Their sign is not
discarded. A noncommuting exact-matrix control checks (11) independently.

Equation (11) is an interface-strength interpolation on a fixed spatial
carrier. It is not a coarse-graining map from the two-site Hilbert space
onto YC6's one-site space, nor a closed scalar recursion for gaps.

## 6. Actual vacuum pairing, signed compass and scale

Let ψθ>0 be the normalized actual six-link ground and dνθ=ψθ²dμ. The
ground equation and integration by parts give, by form closure,

\[
\langle\psi_\theta f,(H_\theta-E_0)\psi_\theta f\rangle_\mu
=\int\sum_{e=1}^{6}|\nabla_e f|^2d\nu_\theta.
\]

Multiplication by ψθ is unitary and preserves the gauge/centre restrictions.
Consequently (1) is equivalently the **actual vacuum** inequality

\[
\int\sum_{e=1}^{6}|\nabla_e f|^2d\nu_\theta
\ge\frac65\operatorname{Var}_{\nu_\theta}(f),\qquad0\le\theta\le2,
\tag{12}
\]

for vacuum-sector form-domain readings f. No product or Haar approximation
to νθ at η=1 is assumed. The shared-link metric is the original six-link
electric metric in (2).

The signed compass rule J(x,y)=(−y,−x) transports source and response
frames together. A unitary change of retained/hidden frames conjugates the
whole Schur return, including (7) and (11); the inertia and (12) are
unchanged. An orientation change cannot remove the nonzero cross source.
No claim is made that this chart itself is a physical observation operator.

With a declared calibration Ephys=(γg²/a)H_{κ/g⁴}+cI, (1) would read
gap(Ephys)≥6γg²/(5a), only when κ/g⁴≤2. This is a strong-coupling window,
not a bound along the weak-coupling continuum trajectory. At fixed a the
spatial box here is (2a,a,a). Increasing one cell count at fixed spacing
is distinct from shrinking a or from refining a fixed physical box. No
renormalized coupling matching has been supplied by this calculation.

**Next obligation:** enlarge the spatial comparison or retained cut while
controlling the interface cross sources and their actual-vacuum form cost
as cells are added. Equation (11) and the resolved source data make that
obligation concrete. A lower bound independent of lattice size, a physical
scale trajectory and a nontrivial continuum field remain unproved.

## 7. Replay and claim boundary

```text
python physics/yc7/yc7_two_cell_return.py --check
python -m unittest discover -s physics/yc7 -p 'test_*.py' -v
```

Eighteen exact certificate checks and thirteen focused tests cover the
complete free cut, both local gauge actions, global centre parity, direct
SU(2) plaquette evaluation, full harmonic source, cross terms, noncommuting
resolvent comparisons, changed-source join, signed transport and every
endpoint/cell. The certificate uses standard-library rationals; the
independent matrix controls use SymPy. The analytic completeness, Schur
and ground-transform arguments are written above, not claimed as formal
proof-assistant verification.

Earlier stage packets are unchanged. Classical tools retained by name:
Peter–Weyl decomposition, SU(2) tensor-product singlet counting,
Rayleigh–Ritz/min–max, Schur/Feshbach return, inverse operator order, the
second resolvent identity and the positive ground-state transform. The new
stage supplies their explicit six-link source data, a finite two-cell gap
certificate and the interface identity (11); it does not claim new general
versions of those classical theorems.
