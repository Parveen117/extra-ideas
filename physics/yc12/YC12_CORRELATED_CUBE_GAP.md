# YC12 — actual cube kinetics and weakly joined correlated blocks

9 October 2026. Continues YC11 at `8f15a7d`, using YC8's local dressing
and YC10's complete product-reference return criterion. Frozen packets
are unchanged. Unit S³ link metric and normalized SU(2) Haar measure
throughout; every electric harmonic is retained.

**Results.** The ordinary twelve-link cube has an actual full-space
quantum gap greater than 1/4 for 0≤θ≤1, including every boundary charge
sector. Its rank-seven gauge-invariant cut has a complete hidden-source
resolution at energies 18, 24, 26 and 32, with rigorous bounds on the
interacting hidden inverse. Disjoint correlated cubes can then be joined
on every even periodic spatial lattice with sides ≥4:

\[
\boxed{0\le\theta_c\le1,\quad0\le\eta_p\le1/16384
\quad\Longrightarrow\quad\Delta_{\rm full}\ge55/256.}       \tag{1}
\]

Here θ_c is the common coupling of the six internal faces of cube c;
η_p is the coupling of an external elementary plaquette. This is an
actual nonuniform Wilson Hamiltonian, not a comparison Hamiltonian or
a tensor product of gauge-reduced cubes. The bound also holds in the
physical gauge-invariant, global-centre-even vacuum sector. It is uniform
in the finite lattice volume. It does **not** enlarge YC10's isotropic
coupling window, establish an RG step or give a continuum mass.

## 1. The owner's “ST plus missing term” has a typed quantum counterpart

Keep the intuition that the retained centre reading needs the contribution
of what its cut misses. On a declared quantum cut P+Q=I write

\[
A=PHP,\quad D=QHQ,\quad PHQ=\theta r^*,\qquad
S(z)=A-z-\theta^2\underbrace{r^*(D-z)^{-1}r}_{\Sigma(z)}.    \tag{2}
\]

For z below the full hidden spectrum, Σ(z) is positive. Solving the hidden
equation and substituting it into the retained equation derives the minus
sign: the omitted response **reduces the effective restoring operator**.
The complete reading includes the return, but its algebraic sign cannot
be chosen by writing “plus some term.” It depends on the target equation.
The inverse in (2) is the full compression; its range is not assumed to
lie in the finite source span.

This is a precise retained/missing-channel correspondence. There is no
proved equality between the thermodynamic product ST, UP9/YC11's centre
scalar U_mid and the retained quantum operator A. They have different
carriers and definitions. YC11 supplies the actual cube configuration
source; the present packet supplies its electric dynamics. The source
parameter λ remains distinct from the Hamiltonian coupling θ.

## 2. Complete low electric cut on the ordinary cube

Let E_c be the 12 edges of the cube graph and P_c its 6 faces. Set

\[
H_c(\theta)=H_0+\theta V,\quad H_0=-\sum_{e\in E_c}\Delta_e,
\quad V=6-\mathcal S,\quad\mathcal S=\sum_{p\in P_c}W_p,
\quad W_p=\tfrac12\operatorname{Tr}U_p.                    \tag{3}
\]

For this section impose Gauss invariance at its eight vertices. Constants
have energy zero. A link of degree n has electric energy n(n+2), with
fundamental degree one costing 3. A nonzero invariant spin-network support
has no degree-one vertex. The cube graph has girth four. A nonempty
support with total energy <18 uses at most five edges and must be one
of the six four-cycles: a connected minimum-degree-two graph with ≤5
edges is a cycle, unless it contains two cycles; two distinct cycles of
the cube require at least seven edges. The cube has no five-cycles.
There is also no room for two nonempty components. Along a degree-two
cycle the singlet condition forces the same representation on all edges.
Only the fundamental four-cycle has energy below 18, and its invariant
space is one-dimensional.

Consequently the **complete** subspace below 18 is

\[
P=\operatorname{span}\{1,b_p=2W_p:p\in P_c\},\quad
PH_0P=\operatorname{diag}(0,12I_6),\quad QH_0Q\ge18.        \tag{4}
\]

The basis is orthonormal. Independent Haar integration of a free boundary
link proves face orthogonality; each holonomy is Haar with E[W²]=1/4.
YC11's proper-subset Haar independence applies to face products here.
At energy 18 there are sixteen six-cycles, including four skew cycles
that are not boundaries of adjacent face pairs. All are in Q, together
with every higher spin network. We do not claim face polynomials are a
complete description of Q.

Let J_6 be the all-ones face matrix, O the permutation exchanging each
opposite pair, and A_o=J_6−I−O the face adjacency matrix. Then

\[
PVP=6I_7-\tfrac12
\begin{pmatrix}0&\mathbf1^T\\\mathbf1&0_{6\times6}\end{pmatrix}.
                                                               \tag{5}
\]

All cubic face moments vanish. With C_p=4W_p²−1 and B_pq=W_pW_q,
the complete columns of r=QVP are

\[
r_0=0,\qquad r_p=-\tfrac12C_p-2\sum_{q\ne p}B_{pq}.       \tag{6}
\]

For adjacent p,q, write Wp=a·u and Wq=b·u in the common unit quaternion,
including traversal signs, and set A_pq=a·b. YC8 proves

\[
H_0C_p=32C_p,\quad H_0A_{pq}=18A_{pq},\quad
H_0B_{pq}=26B_{pq}-2A_{pq}.                               \tag{7}
\]

Thus each adjacent term −2B splits into −A/2 at energy 18 and
−2(B−A/4) at 26. An opposite-face product has energy 24. The Haar norms
are ||C||²=1, ||A||²=1/4, ||B||²=1/16 and ⟨A,B⟩=1/16; conditioning on
the common link gives E_u B=A/4. Hence ||B−A/4||²=3/64.
Different features at the same energy have different per-link harmonic
degrees, so are orthogonal. For adjacent pairs the six-link boundaries
are all distinct; the degree-two common link distinguishes the 26 modes.
This resolves the whole source, with face-block Gram matrices

\[
\boxed{B_{18}=(4I+A_o)/16,\quad B_{26}=3B_{18},\quad
B_{24}=(I+O)/4,\quad B_{32}=I/4.}                         \tag{8}
\]

Their vacuum row/column is zero, and their sum is (5I+J_6)/4.
In particular the face columns are not treated as independent returns.

| Face subspace | Dimension | B18 | B24 | B26 | B32 |
|---|---:|---:|---:|---:|---:|
| Uniform | 1 | 1/2 | 1/2 | 3/2 | 1/4 |
| Opposite-even, zero sum | 2 | 1/8 | 1/2 | 3/8 | 1/4 |
| Opposite-odd | 3 | 1/4 | 0 | 3/4 | 1/4 |

These are source eigenvalues, not quantum excitation energies.

## 3. Full hidden inverse, including the interaction

The elementary bounds 0≤V≤12 imply, for θ≥0 and z<18,

\[
\sum_E\frac{B_E}{E+12\theta-z}
\preceq\Sigma(z)\preceq\sum_E\frac{B_E}{E-z}.             \tag{9}
\]

A sharper centred envelope uses a cube-specific parity. There is a link
centre-sign assignment flipping every Wp: equivalently the six face
parity equations have an even total on this closed cube. A concrete
twelve-bit assignment is supplied by the replay. Its unitary J commutes
with H0, preserves P,Q and sends S to −S. The vacuum is J-even, each bp
is J-odd, and every nonzero source column in (6) is J-even.

Write z=6θ+ζ, D0=QH0Q and T=QSQ. Then ||T||≤6 and JTJ=−T.
Put A0=D0−ζ and C=θ A0^(−1/2) T A0^(−1/2). For

\[
\zeta<18,\quad u=6\theta/(18-\zeta)<1,\qquad
R_0(\zeta)=\sum_E B_E/(E-\zeta),
\]

the normalized source A0^(−1/2)r has even parity. All odd powers of C
vanish between these source vectors; every even power is positive.
The Neumann series is norm convergent, so

\[
\boxed{R_0(\zeta)\preceq\Sigma(6\theta+\zeta)
\preceq\frac{R_0(\zeta)}{1-u^2}.}                         \tag{10}
\]

This includes all repeated hidden interactions, even though T sends
the source outside its span. It uses the full hidden floor 18 and a
bounded perturbation, not invariance of the listed four source shells.

Substitute the upper and lower sides into (2). Below the hidden floor,
the number of negative eigenvalues of S(z) equals the number of full
gauge-carrier eigenvalues below z. Ordering the two pencils brackets
these counts. Exact rational congruence (YC4) and rational bisection give:

| θ | Actual gauge gap E1−E0, outward decimal enclosure |
|---:|---:|
| 1/4 | [11.9993, 12.0004] |
| 1/2 | [11.9819, 12.0017] |
| 1 | [11.2838, 12.0096] |

The result file retains exact fractions and ground/excitation enclosures
separately. These are point certificates, not an interpolated full-window
gauge-gap claim. For joining blocks we need the full carrier, including
charged states which already have free energy 3 rather than 12. The next
argument supplies that distinct requirement on an entire interval.

## 4. A full-space cube gap: finite residual makes the transfer possible

Return to L²(SU(2)^12) without imposing independent block boundary Gauss
constraints. Use YC8's exact positive reading

\[
f_1=-\mathcal S/12,\qquad
f_2=\sum_p C_p/4608-\sum_{p\sim q}A_{pq}/1404
                         +\sum_{p\sim q}B_{pq}/1872,
\quad F=\theta f_1+\theta^2f_2,\quad\psi=e^{-F}.            \tag{11}
\]

For the six-face cube its identities give exactly

\[
H_c\psi/\psi=6(\theta-\theta^2/48)+R,\qquad
R=-2\theta^3\Gamma(f_1,f_2)-\theta^4\Gamma(f_2).            \tag{12}
\]

The coefficients in (11) do not require four faces per edge. Only the
incidence bounds change: a cube edge belongs to two faces, each face has
four adjacent faces, and each edge belongs to seven adjacent-pair unions.
Exactly one pair shares that edge, while six contain it noncommonly.
The joined A loops therefore use it six times, the B derivatives eight
times when the common edge is counted twice. YC8's quaternion derivative
bounds imply

\[
|\nabla_e f_1|\le1/6,\qquad
|\nabla_e f_2|\le k:=8/4608+6/1404+8/1872=77/7488.         \tag{13}
\]

For the Hessian block-row sums, W has bound 1 per block; C has bound
8 on the same link and 16 on the other three links; A has bound 1 on
its six links; B has bound 4 on its seven links. Consequently

\[
\|\mathrm{Hess} f_1\|\le2/3,\qquad
\|\mathrm{Hess} f_2\|\le
112/4608+36/1404+196/1872=193/1248.                        \tag{14}
\]

The comparison K=Hc−6(θ−θ²/48)−R has ground ψ at zero, and its
ground transform is the weighted diffusion for ψ²dμ. Product Ricci
curvature is 2. Integrated weighted Bochner, as derived in YC8, gives

\[
\operatorname{gap}K\ge\rho(\theta)
=2-4\theta/3-(193/624)\theta^2.                           \tag{15}
\]

Summing (13) over **twelve** edges and retaining the nonpositive quartic
term in (12) gives

\[
-4k\theta^3-12k^2\theta^4\le R\le4k\theta^3,\qquad
\operatorname{osc}R\le8k\theta^3+12k^2\theta^4.            \tag{16}
\]

Min–max yields E1(Hc)≥scalar+ρ+inf R and E0(Hc)≤scalar+sup R.
Thus the **actual** full-space gap satisfies

\[
\boxed{\operatorname{gap}H_c\ge G(\theta)
=2-4\theta/3-(193/624)\theta^2-8k\theta^3-12k^2\theta^4.} \tag{17}
\]

All nonconstant coefficients of G are negative. On [0,1] its minimum is

\[
G(1)=1279511/4672512>1/4,
\quad G(1)-1/4=111383/4672512.                            \tag{18}
\]

The actual cube ground is unique and strictly positive by the compact
elliptic Schrödinger positivity argument; it is invariant under every
vertex gauge transformation. It is not ψ unless R is constant. The
bound (17) applies to all states above it, including nonsinglet boundary
charges. The gauge gap in section 3 cannot replace this bound in a block
join. Conversely the conservative 1/4 does not claim that the physical
gauge excitation is so low.

YC8 could not transfer its global comparison gap by subtracting an
extensive residual on a growing lattice. Here the residual is paid on
a fixed twelve-link block before the all-volume argument. That is the
reason this transfer is available.

## 5. Disjoint correlated cubes cover an actual lattice

Take a rectangular periodic lattice with every side an even integer ≥4.
Partition vertices by the blocks {2m_i,2m_i+1} in each direction. A link
with both endpoints in one block is an internal cube edge. There are
twelve per block. Every other link is a bridge, retained as its own
factor. These disjoint edge sets cover the lattice exactly, hence

\[
\mathcal H=\bigotimes_c L^2(SU(2)^{12})
             \otimes\bigotimes_{e\ {m bridge}}L^2(SU(2)). \tag{19}
\]

This is a true factorization before gauge reduction. No shared physical
link is duplicated, and no singlet constraint is independently imposed
at a cube's boundary. The reference operator is

\[
H_{\rm ref}=\sum_c H_c(\theta_c)
                 +\sum_{e\ {m bridge}}(-\Delta_e).
\]

Subtract each reference ground energy. Its unique ground is the product
of the actual correlated cube grounds and the constant bridge readings.
Each factor has gap at least g=1/4; bridge factors have gap 3. Each
ground is invariant under the relevant endpoint gauge actions, so the
product respects all lattice Gauss symmetries. Keep the complete excited
space of each factor, not just its gauge-singlet excitations.

Each remaining elementary face is one of two types, classified by the
parity of its two coordinate intervals:

| External face | Factor support | Count per vertex block |
|---|---|---:|
| One interval crosses a block boundary | two cubes, two bridges | 12 |
| Both intervals cross block boundaries | four bridges | 6 |

If neither interval crosses, the face is internal (six per block).
The sides ≥4 ensure the displayed factors are distinct even at periodic
seams. A cube factor meets 24 external faces: each of its twelve edges
has two noninternal incident faces, with no repeated cube incidence
within a face. A bridge meets four faces. Thus for Φ=−Σ_external ηp Wp,

\[
k_{\rm support}=4,\qquad
\beta:=\max_{\rm factor\ i}\sum_{p\ni i}\|\eta_p W_p\|
\le24\eta,\qquad \eta=\max_p|\eta_p|.                    \tag{20}
\]

These are all-size parity/incidence proofs. The replay checks three
different boxes as controls; finite enumeration is not the proof of
uniformity. In particular the old single-link incidence β=4η cannot be
used for a cube factor. A cube already supports a nonconstant invariant
excitation on one factor. Use **s=1**, not YC10's single-link s=4.

## 6. Full cluster return and uniform gap

Apply YC10 section 5 with the now-established correlated reference:

\[
g=1/4,\quad k=4,\quad s=1,\quad M=2^{k/2}=4,
\quad r=1/8,\quad a=2k/s=8,\quad e^{ar}=e<3.              \tag{21}
\]

The proof does not need a finite factor dimension or a known explicit
formula for its ground. If Ωi is that unique ground and qi its orthogonal
complement, the excited support I has H_I≥g|I|. The exact creation map
is c_I=−H_I^−1 P_I exp(−C)Φexp(C)Ω, C=Σ ĉ_I. The orthogonal binary
projections over a four-factor interaction still cost at most √16=4.
The weighted norm max_i Σ_(I containing i)|I| ||c_I|| cancels the free
inverse. These are precisely YC9's complete cluster construction and
YC10's refined estimates, with every excited block state retained.

With β≤24η, their three bounds are

\[
\|\mathcal T(c)\|_*\le60\beta\le1440\eta,\qquad
\operatorname{Lip}\mathcal T\le576\beta\le13824\eta,
\qquad b\le96\beta\le2304\eta.                           \tag{22}
\]

At η=1/16384 they respectively equal 45/512<1/8, 27/32<1 and 9/64<1.
For the exponential bound a rational series estimate is
e≤1+1+1/2+(1/6)/(1−1/4)=49/18<3. Banach contraction constructs the
full dressed reference eigenvector; all cluster sizes are included.
At finite volume the creation sum is bounded and nilpotent; each fixed
point coefficient is in Dom H_I, so the similarity preserves the full
reference domain. The quotient excitation operator is H_ref+QUQ with

\[
\|QUQH_{\rm ref}^{-1}\|_{1\to1}\le b\le9/64.
\]

For z<g(1−b), the relative Neumann inverse exists just as in YC9.
The triangular similarity then proves that the constructed eigenvector
is the unique true ground and that no other full-space eigenvalue lies
below g(1−b). In fact the displayed estimates give the useful formula

\[
\Delta_{\rm full}\ge\tfrac14(1-2304\eta)
\ge55/256.                                               \tag{23}
\]

The unique positive true ground is gauge invariant and even under the
three global centre flips. Restricting to that sector cannot decrease
its excitation gap. This proves (1) on every allowed finite volume;
construction of a particular infinite-volume representation or
continuum limit is not claimed here.

The classical method is weak-interaction / Kirkwood–Thomas dressing.
A primary reference allowing infinite-dimensional single-site spaces is
D. A. Yarotsky, *Quasi-particles in weak perturbations of non-interacting
quantum lattice systems*, section 2, https://arxiv.org/abs/math-ph/0411042.
The numerical constants here come from the reproduced YC9/YC10 estimates
and the present block geometry, not a quotation of that paper's constants.
For the positive local-energy method see A. Mouchet,
https://arxiv.org/abs/math/0505541; the min–max and weighted Bochner steps
used for the gap are explicitly displayed above and in YC8.

## 7. What this closes, and what remains

YC11 ended with an actual configuration source but no electric block
inverse or gap. Equations (6)–(10) return that block's complete kinetic
source. Equations (17)–(23) establish the correlated reference and its
first volume-uniform join, including the boundary charges and interface
norms required by YC10. The gauge cut and the full-space block are kept
distinct throughout.

The internal coupling may now be of order one, while the interface must
remain ≤1/16384 in this certificate. Setting every coupling equal gives
no improvement over YC10's isotropic 3/160 window. A genuine scale step
must control interactions between strongly correlated blocks and the
ratio β/g after elimination. Neither the centre intuition nor the small
interface theorem performs that step automatically. The next target is
to improve the **boundary source/inverse and its interface norm**, rather
than replace all boundary response by ||Wp||∞=1. Arbitrarily strong
interfaces, an isotropic weak-coupling path and physical scale matching
remain open.

## Replay

Run `python physics/yc12/yc12_correlated_cube_gap.py --check` and
`python -m unittest discover -s physics/yc12 -p 'test_*.py'` from the repo
root. The result stores exact fractions, source-shell Grams, rational
inertia enclosures, geometric controls and SHA256 source pins. Tests
include free-spectrum completeness, skew channels, a full-hidden parity
model that leaves the source span, invalid geometries and the charged
boundary distinction. This is a written analytic proof with exact
computational certificates, not a proof-assistant formalization.
