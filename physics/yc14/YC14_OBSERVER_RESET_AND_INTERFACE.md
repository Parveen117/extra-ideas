# YC14 — a locally aligned observer retains memory; cube interfaces cost less

9 October 2026. Continues YC13 at `c730847`, with the owner's subsequent
correction to the interpretation of lambda zero. Frozen packets are unchanged.
The source, metric, complete carrier and transfer hypotheses are specified below.

**New Yang–Mills result.** In YC12's unit-S³ electric normalization, on every
periodic spatial lattice with even sides at least four,

\[
0\le\theta_c\le1,\qquad 0\le\eta_p\le\eta\le1/5120
\quad\Longrightarrow\quad
\boxed{\Delta_{\rm full}\ge\frac14-\frac{1024}{5}\eta
                         \ge\frac{21}{100}.}                 \tag{1}
\]

Here theta_c couples the six internal faces of correlated cube c and eta_p
couples each external face. The complete Hilbert space includes all boundary
charges, electric harmonics and cluster sizes. The conclusion also holds in
the physical gauge-invariant, centre-even vacuum sector. The allowed interface
coupling is **16/5 times** YC12's 1/16384; at that old endpoint the new floor is
19/80 rather than 55/256. This is not a larger isotropic window, an RG step,
or a four-dimensional continuum mass-gap theorem.

The gain comes from resolving the *first* boundary source before applying the
full nonlinear return estimate. It does not assume that later returns keep
the first source's four-factor support.

## 1. Corrected observation convention

The owner clarified after YC13:

- Lambda zero is the balanced cut at which pure observation is possible:
  symmetry/equilibrium and balanced TVSP-to-VTSP exchange, with distinction
  still present. It is not merely the first member of a movie-potential family.
- A new frame changes the observer's relative orientation. The observer
  realigns to the diagonal/local zero cut for the next reading.
- Retained ratio-memory distinguishes frames after this local realignment.
- Radial motion is undecided. No contraction, black-hole motion or physical
  clock is imposed. The proposed connection to the missing YM middle is a
  question to test, not an identification of an existing gap theorem.

YC13's potential-family and transport theorems remain true for their declared
parameters. They do not yet identify their source amplitude with this
observation condition. Below **local reset** means an explicit isometry of
readings and their pairing. We do not define the owner's lambda by fiat as
w, eta, the Hamiltonian coupling, or a numerical observer angle. TVSP before
observation and VTSP after remain the declared convention; the signed arrow
swap and the unsigned four-slot permutation retain their distinct actions.

## 2. Exact diagonal alignment: retain the metric and the turn

Use UP8's real two-component carrier, R²=−I, K²=I, KR=−RK. On Delta>0,

\[
\kappa=\cosh\eta\,k(\phi)+\sinh\eta\,R,\quad
k=\cos\phi K+\sin\phi RK,\quad n=Rk,
\quad B=(R\kappa)^2=e^{-2\eta n}.                         \tag{2}
\]

Define the following **cut metric** and observer map:

\[
G=B^{1/2}=e^{-\eta n},\qquad
C=O B^{1/4},\qquad O=e^{-\phi R/2}.
\]
\[
\boxed{C\kappa C^{-1}=K,\quad C^TC=G,\quad
       \kappa^TG=G\kappa,\quad\det C=1.}                  \tag{3}
\]

Proof: n²=I and [n,k]=2R, so
kappa=exp(eta n/2) k exp(−eta n/2); conjugation by B^(1/4) removes
that factor, and O k O^T=K. The remaining identities follow directly.
G is positive because it is the square of the nonsingular symmetric
positive matrix B^(1/4). This G is **not** the metric B used for UP8's
half-share connection. The different roots serve different purposes.

The cut projections P±=(I±kappa)/2 are orthogonal in G. For v on this carrier,

\[
\|v\|_G^2=\|P_+v\|_G^2+\|P_-v\|_G^2=\|Cv\|^2,
\]
\[
v^TG\kappa v=\|P_+v\|_G^2-\|P_-v\|_G^2
              =(Cv)^TK(Cv).                              \tag{4}
\]

Thus the balanced seam becomes the ordinary diagonals of the reset reading.
Both total information norm and signed contrast are retained. This is an
exact finite-dimensional realization of an undistorted local *reading*, not
a claim that a physical measurement of every observable is nondisturbing.

On the genuinely flat cut eta=0, G=I and C=O: **rotation alone suffices**.
For eta≠0 an orthogonal rotation preserves the antisymmetric coefficient
w, so it cannot remove obliqueness. The positive factor has reciprocal
singular values exp(±eta/2); no isotropic radial motion follows from it.

Crucially, also carry the primitive turn:

\[
R'=CRC^{-1},\qquad
\boxed{(R'K)^2=CBC^{-1}.}                                \tag{5}
\]

An aligned local cut does **not** erase the full return. If eta≠0, keeping
R numerically fixed after the nonorthogonal reset would incorrectly replace
the genuine return by (RK)²=I. The eigenvalues exp(±2eta) of the actual
cycle remain. This is why local pure reading and retained memory can coexist.
The signed observation chart J=[[0,−1],[−1,0]] transports G, kappa, R and
the source together. It does not alter (3)–(5) as geometric identities.

## 3. Reset each frame, keep the transition ledger

Let C_i reset reading i, and let T_(j<-i) be a **specified physical or native
transport on that same carrier**. The transition seen by the reset observer is

\[
\widetilde T_{j\leftarrow i}=C_jT_{j\leftarrow i}C_i^{-1},
\quad
\prod_{i=0}^{N-1}\widetilde T_{i+1\leftarrow i}
 =C_N\left(\prod_{i=0}^{N-1}T_{i+1\leftarrow i}\right)C_0^{-1}. \tag{6}
\]

The products are ordered with later steps to the left. Every intermediate
reset cancels, while the transport history remains. An O(lambda,x)-dependent
observer must retain −dO O^T in its connection, as in YC13. The same formula
with −dC C^−1 applies to a general invertible frame on the appropriate carrier.
UP8's orthonormal transport and (2)'s original coefficient carrier must be
intertwined by UP8's h=B^(1/2) before applying (6); their matrices are not
silently substituted for one another.

A simple flat example takes C0=I and rotations C1,C2 with cosine/sine pairs
(3/5,4/5) and (5/13,12/13), and kappa_i=C_i^T K C_i. Every local reset
shows K, yet M10=C1 C0^T and M21=C2 C1^T have slope ratios 4/3 and 16/63.
Their complete ratio-memory composes as

\[
\frac{4/3+16/63}{1-(4/3)(16/63)}=12/5.                   \tag{7}
\]

This is the tangent addition law, with the full rotation matrix retained
at ratio poles and for angle branches. These open-path ratios refer to
specified endpoint frames; they are not frame-independent observables.
If T=I, a closed ledger with C_N=C_0 has identity product. Relabelling alone
does not manufacture curvature or energy. For a genuine closed transport,
(6) instead preserves its conjugacy class. A loop winding an angular patch
must also retain the sign/transition of its chosen half-angle frame.

## 4. The equilibrium observation map on the actual YM reference

Use YC12's exact factorization into disjoint twelve-link cubes and individual
bridge links. Each cube has its **actual**, unique positive normalized ground
Omega_c and a full excitation floor g=1/4 for theta_c in [0,1]. Its explicit
comparison reading from YC8 is not substituted for Omega_c. Put

\[
\Omega=\bigotimes_c\Omega_c\otimes\bigotimes_{e\,\mathrm{bridge}}1,
\quad H_{\rm ref}=\sum_c(H_c-E_c)+\sum_{e\,\mathrm{bridge}}(-\Delta_e),
\quad d\nu=\Omega^2d\mu.
\]

The classical ground-state transform U_Omega f=Omega f is a unitary from
L²(nu) to L²(mu), with

\[
\langle\Omega,F\Omega\rangle=\int F\,d\nu,\qquad
\langle\Omega f,H_{\rm ref}\Omega f\rangle
       =\sum_e\int|\nabla_e f|^2d\nu.                    \tag{8}
\]

Here F is a bounded multiplication reading. Integration by parts and
H_ref Omega=0 prove the form identity, first for smooth f and then by form
closure. It retains both the equilibrium measure and all kinetic directions.
Equation (8) supplies a precise quantum/probability equality for this class
of readings; it does not turn arbitrary real-time quantum measurements into
classical ones. The transformed operator has the same excitation gap.

| Transfer | Carrier/source | Pairing | Intertwiner and target |
|---|---|---|---|
| Diagonal reset | UP8 admissible two-component cut | G=B^(1/2) | C carries norm and contrast to the Euclidean diagonal |
| Equilibrium reading | Complete cube/bridge Hilbert space, actual Omega | Haar mu ↔ Omega² mu | U_Omega preserves the exact quadratic form and gap |
| Interface return | r_p=W_p Omega, complete excited factor spaces | Haar Hilbert norm | H_ref inverse and the full cluster similarity retain the source's energy cost |

There is no claimed operator map from the two-component compass into all
Yang–Mills degrees of freedom. The useful shared prescription is explicit:
realign the reference, transport the metric, and keep the returning source.

## 5. Every first interface source excites all four factors

Discard the additive scalar sum eta_p in the Wilson Hamiltonian; it does
not change any gap. The interaction is Phi=−sum_p eta_p Wp. For each external
face let r_p=Wp Omega. YC12 gives two possible supports:

| Face type | Excited factors in r_p | Source energy floor | Norm² | Inverse-vector norm bound |
|---|---|---:|---:|---:|
| One crossing interval | two cubes, two bridges | 6+2g=13/2 | 1/4 | 1/13 |
| Two crossing intervals | four bridges | exactly 12 | 1/4 | exactly 1/24 |

These are statements about the **whole source vector**, not selected low
harmonics. A bridge occurs once in its face holonomy, so the source has
fundamental degree one and electric energy 3 on that bridge, and is orthogonal
to its constant ground. For a cube factor the face uses one open internal edge.
A centre gauge transformation at either endpoint changes that edge's sign,
leaves every internal cube Wilson loop unchanged, and fixes Omega_c. Averaging
against Omega_c² therefore kills every matrix entry of that open edge.
Projection onto the cube ground vanishes. Both cube factors are excited.
Their complete excited spaces, including boundary charges, each cost at least g.

Condition on all links except one bridge. Its holonomy is Haar, irrespective
of the other links and their correlated cube density. Thus

\[
\|r_p\|^2=\mathbb E_\nu W_p^2=1/4,
\qquad \mathbb E_\nu W_p=0.                              \tag{9}
\]

H_ref preserves every per-bridge harmonic sector and every cube vacuum/excited
projection. Its inverse on the source's four-factor sector consequently has
the energy floors displayed in the table. No unknown explicit ground function
is required for these estimates.

There is also an exact energy sum rule, including interacting cubes. On each
of the four unit-S³ links, |grad_e Wp|²=1−Wp². Equation (8) gives

\[
\boxed{\langle r_p,H_{\rm ref}r_p\rangle=3,
\qquad \langle H_{\rm ref}\rangle_{r_p}=12.}              \tag{10}
\]

Let D0 be H_ref compressed to the full complement of Omega. Cauchy–Schwarz
in its spectral measure, followed by the complete source-sector floor, yields

\[
\boxed{\frac1{48}\le\langle r_p,D_0^{-1}r_p\rangle
                   \le\frac1{26}}                      \tag{11}
\]

for a two-bridge face, and exactly 1/48 for a four-bridge face. The lower
bound is ||r||⁴/<r,Hr>. At theta_c=0 the two-bridge source too has exact
energy 12. At nonzero theta_c, (10) fixes its mean, not its whole spectral
distribution; replacing the inverse by division by the mean is invalid.

Distinct external faces have distinct sets of bridge edges. For two-bridge
faces the parallel pair and its neighboring endpoints uniquely reconstruct
the face; for four-bridge faces the four edges reconstruct its elementary
square. Side lengths at least four exclude the two-edge periodic ambiguity.
There is therefore a bridge on which one source has degree one and the
other degree zero. D0 inverse preserves this distinction, proving

\[
\langle r_p,r_q\rangle=\tfrac14\delta_{pq},\quad
\langle r_p,H_{\rm ref}r_q\rangle=3\delta_{pq},\quad
\langle r_p,D_0^{-1}r_q\rangle=0\quad(p\ne q).            \tag{12}
\]

This is the complete reference inverse. The interacting inverse does not
preserve those bridge signatures; later returns are handled next, not set
to zero by (12). The order-eta² ground shift is −sum eta_p² times the
susceptibilities (11); this perturbative coefficient alone is not a gap bound.

## 6. A source-centred full-return ball gives the larger lattice window

Keep YC9/YC10/YC12's full creation map and its weighted local norm:

\[
\mathcal T(c)_I=-H_I^{-1}P_Ie^{-\mathcal C}\Phi e^{\mathcal C}\Omega,
\quad \|c\|_* =\max_i\sum_{I\ni i}|I|\,\|c_I\|.
\]

The first source has |I|=4, so the preceding inverse-vector bounds give
weighted costs 4/13 per two-bridge face and 1/6 per four-bridge face.
A cube meets 24 two-bridge faces; a bridge meets two of each type. The latter
count follows by inspecting each of its two transverse directions: one of
the positive/negative intervals stays inside a cube and the other crosses.
Hence the initial full-map seed obeys

\[
\boxed{s_0:=\|\mathcal T(0)\|_*\le\frac{96}{13}\eta,}
\qquad s_{0,\mathrm{bridge}}\le\frac{37}{39}\eta.           \tag{13}
\]

Later creators may occupy a single cube. Retain **s=1**, g=1/4, support k=4,
M=sqrt(2^k)=4, a=2k/s=8, and beta<=24eta. YC10's proven derivative and
excitation estimates, on a ball of radius r, are

\[
q\le\frac{M\beta}{g}e^{ar}[2+a(1+2r)],\qquad
b\le\frac{2M\beta}{g}e^{ar}.                             \tag{14}
\]

Replace only the coarse absolute map bound by the elementary consequence
of the same Lipschitz estimate:

\[
\boxed{\|\mathcal T(c)\|_*\le s_0+q r.}                  \tag{15}
\]

This is valid on the same complete Banach space; it does not require T(0)
to be the fixed point. Choose r=1/128. The entire exponential series satisfies
exp(8r)=exp(1/16)<sum_(n>=0)(1/16)^n=16/15. At eta=1/5120 the exact bounds are

\[
s_0\le3/2080,\quad q\le81/100,\quad b\le4/25,
\quad s_0+qr\le1293/166400<1300/166400=r.                  \tag{16}
\]

All left sides decrease with eta, so the entire interval is covered. The
Banach fixed point retains all cluster orders and all factor harmonics. Its
iteration remainder is at most q^m s0/(1−q), and its norm at the endpoint
is at most 15/1976<1/128. This is a local cluster norm, not a volume-uniform
Hilbert norm approximation to an explicitly computed vacuum.

For completeness, the remaining spectral transfer is exactly YC9's full
domain argument. At each finite volume there are finitely many subset
creators, even though each excited factor is infinite dimensional; the fixed
point components belong to Dom H_I. Their bounded creation sum preserves
Dom H_ref. The similarity exp(−C)(H−E)exp(C) is triangular relative to Omega,
with quotient H_ref+QUQ. Its component-sum relative norm is at most b in (14).
For real z<g(1−b), the relative Neumann inverse excludes every nonzero
eigenvalue there, including negative z. The zero eigenvalue is simple, so E
is the actual ground, not merely an assumed branch. Thus

\[
\Delta_{\rm full}\ge g(1-b)
\ge\tfrac14(1-\tfrac{4096}{5}\eta),                       \tag{17}
\]

which is (1). Positivity and symmetry put the ground in the physical vacuum
sector; restricting to that sector cannot lower its excitation gap.
No lower-dimensional truncation or separate block singlet restriction enters.

The classical method is the ground-state transform and weak-interaction
Kirkwood–Thomas expansion. D. A. Yarotsky's primary reference
[Quasi-particles in weak perturbations of non-interacting quantum lattice
systems](https://arxiv.org/abs/math-ph/0411042) treats general gapped factors;
YC9/YC10 give the estimates reproduced here. The constants (13)–(17) are
derived here from actual cube boundary sources, not quoted from that paper.

## 7. Evidence, status and next target

Run `python physics/yc14/yc14_observer_reset_and_interface.py --check` and
`python -m unittest discover -s physics/yc14 -p 'test_*.py'` from the repo root.
The replay has 44 exact checks and 15 new tests; the 36 predecessor tests
in UP8, YC10 and YC12 also pass. The result pins its note, code, tests and
frozen inputs with SHA256. Exact
symbolic/rational checks accompany the written analytic proofs; they are not
a formal proof-assistant verification. Three periodic boxes check incidence
and source signatures; the parity/uniqueness arguments above prove all sizes.

The owner's reset-and-ratio distinction has a precise local representation.
Its useful YM translation is to keep the exact equilibrium pairing and pay
the actual source inverse before estimating the remaining return. That produces
a stronger certified interface theorem here. It does not establish that pure
observation is the missing continuum mass mechanism.

Next target: resolve how the first bridge source returns into single-cube and
neighboring-cube channels, and sharpen the nonlinear inverse without falsely
retaining four-factor support. The physical scale/energy conversion and a
continuum construction remain separate obligations. Radial motion stays open.
