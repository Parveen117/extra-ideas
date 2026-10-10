# YC9 — shared-plane joining and an actual volume-uniform vacuum gap

9 October 2026. Base: YC8, `15ae97f`. Carrier and metric are unchanged:
periodic rectangular three-dimensional spatial SU(2) lattices with every
side length at least 3, unit S³ link metrics, normalized Haar measure,
all local Gauss laws, and the even character of each global centre flip.

**Result.** Write P for the number of plaquettes. For the actual operator

\[
H_\theta=-\sum_e\Delta_e+\theta\sum_p(1-W_p),
\qquad W_p=\tfrac12\operatorname{Tr}U_p,
\]

the vacuum-sector excitation gap satisfies

\[
\boxed{\Delta_{\rm vac}(\theta)\ge
12-\frac{2560}{3}\theta\ge\frac{32}{3},
\qquad 0\le\theta\le\frac1{640}.}                         \tag{1}
\]

The constants do not depend on any side length. This is a theorem for the
original finite-volume YM Hamiltonians, not YC8's comparison operator.
Its window is deliberately small. The proof covers the full, untruncated
link Hilbert spaces and returns all excitation clusters by a convergent
fixed point. It does not give a weak-coupling continuum limit, identify the
first excited multiplet, or assert a constructed infinite-volume field.

The general weak-interaction mechanism is classical: Kirkwood–Thomas /
coupled-cluster ground-state dressing and stability of product vacua.
The packet supplies its own estimates and explicit constants for this
four-link interaction, including the physical-sector restriction. It is
not presented as a newly discovered general existence of strong-coupling
lattice gaps.

## 1. What the owner's three intersecting planes supply

The owner proposes joining two more TVSP planes to one, at a common centre,
keeping shared edges and their signs. On a spatial lattice the precise
counterpart is the three oriented face families xy, yz and zx meeting at
common vertices and links. Three transverse coordinate planes have eight
local octants; two have four. Octants are spatial regions, not eight
independent quantum sectors or the eight global centre characters.

| Link direction | Incident plane families | Plaquette incidences |
|---|---|---:|
| x | xy and xz | 2 + 2 |
| y | xy and yz | 2 + 2 |
| z | xz and yz | 2 + 2 |

A shared link is one SU(2) variable, with one kinetic term. Three
independent planar tensor products would duplicate it. Before Gauss
restriction the correct factorization is over **links**,
\(\mathcal H=\bigotimes_e L^2(SU(2))\). Gauge transformations at their
common vertex act on all incident links together. The physical subspace
does not inherit an independent tensor factor for each plane.

An oriented boundary uses inverse transport on a reversed link:

\[
U_{x,ij}=U_{x,i}U_{x+\hat i,j}
 U_{x+\hat j,i}^{-1}U_{x,j}^{-1}.
\]

Reversing a loop inverts its holonomy; its normalized SU(2) trace is
unchanged. Noncommuting products and YC8's cross sources are retained.
In a declared smooth small-loop limit, plane-indexed curvature is the
antisymmetric tensor/2-form \(F_{ij}\); the finite lattice object here is
the ordered group holonomy. The compass chart is not substituted for a
link orientation or an assumed continuum limit.

The useful consequence for the gap is a **local join budget**, rather
than the total number of faces. Allow nonuniform real couplings t_p and set

\[
H=H_0+\sum_p t_p(1-W_p),\quad
\beta=\max_e\sum_{p\ni e}|t_p|.                            \tag{2}
\]

The more general conclusion proved below is

\[
\boxed{\beta\le\frac1{160}\quad\Longrightarrow\quad
\Delta_{\rm vac}\ge12\left(1-\frac{160}{9}\beta\right)
\ge\frac{32}{3}.}                                         \tag{3}
\]

The scalar \(\sum_p t_p\) is removed during the proof and restored at the
end; it changes no gap. Uniform coupling gives β=4θ and (1).

For plane-dependent couplings a=t_xy, b=t_yz, c=t_zx,

\[
\beta=2\max(|a|+|c|,|a|+|b|,|b|+|c|).                    \tag{4}
\]

Thus adding planes along (a,b,c)=(t,st,st), 0≤s≤1 and 0≤t≤1/640,
has β=2t(1+s) and keeps the actual vacuum gap at least 32/3 throughout.
This is a spatial join path at fixed normalization, not an RG step.

## 2. Excitation cuts and their full source

For each link let Ω_e=1, p_e=|Ω_e><Ω_e| and q_e=1−p_e. The Hilbert
space has the exact orthogonal decomposition

\[
\mathcal H=\bigoplus_{I\subset E}\mathcal H'_I\otimes\Omega_{E\setminus I},
\qquad \mathcal H'_I=\bigotimes_{e\in I}q_e L^2(SU(2)).       \tag{5}
\]

Let P_I extract this component, with P_empty the scalar vacuum component.
Every nonempty I has

\[
H_I:=\left.H_0\right|_{\mathcal H'_I}\ge3|I|,
\qquad\|H_I^{-1}\|\le\frac1{3|I|}.                         \tag{6}
\]

This uses the entire single-link spectrum n(n+2), n≥1, not a harmonic
truncation. For c_I∈H'_I define the bounded creation operator
\(\widehat c_I=|c_I\rangle\langle\Omega_I|\), extended by the identity
outside I. Its norm is ||c_I||. Two such operators commute: their product
is zero when their supports overlap, and a tensor product otherwise.

Put C=Σ_(I nonempty) ĉ_I and Ψ=exp(C)Ω. This is a creation-operator
exponential, not multiplication by YC8's scalar exp(−F). At each finite
volume C is bounded and nilpotent; exp(C) is invertible and
<Ω,Ψ>=1. Use the weighted local norm

\[
\|c\|_*:=\max_e\sum_{I\ni e}|I|\,\|c_I\|.                 \tag{7}
\]

Its weight will cancel the inverse in (6) and the number of places where
a cluster can meet another face. No fixed maximum cluster size is imposed.

Let Φ=−Σ_p t_p W_p, so H−Σ_p t_p=H0+Φ. Set
\(W(c)=e^{-C}\Phi e^C\). Since
\([H_0,\widehat c_I]=\widehat{H_Ic_I}\) and creation operators commute,
the eigenvector equation is exactly the fixed point

\[
\boxed{\mathcal T(c)_I=-H_I^{-1}P_I W(c)\Omega,\quad I\ne\varnothing.}
                                                               \tag{8}
\]

At a solution the eigenvalue, with its scalar restored, is
\(E=\sum_p t_p+P_\varnothing W(c)\Omega\).

## 3. Local commutators and the projection-count lemma

We state the estimates for interactions Φ_X on sets X of size at most
k, with \(\max_e\sum_{X\ni e}\|\Phi_X\|\le\beta\). Here k=4 and
||−t_p W_p||=|t_p|.

In the commutator expansion of e^(−C) Φ_X e^C, **each** participating
I_j must intersect X. Indeed a creation operator missing X commutes with
Φ_X and with every other creation operator, so the Jacobi identity lets
its zero commutator be taken first. It need not be bounded by the size of
an expanding union. Consequently

\[
\sum_{I:I\cap X\ne\varnothing}\|c_I\|\le k\|c\|_*.        \tag{9}
\]

For a fixed ordered tuple I_1,…,I_l, expand a nested commutator into its
2^l operator words and apply a word to Ω. Every nonzero word has all links
in (I_1∪…∪I_l)\X excited, and every link outside X∪I_1∪…∪I_l in its
vacuum. Only the at most k links of X can change their binary vacuum /
excited status. There are therefore at most 2^k nonzero P_J components.
Each projection is a contraction, so their sum of Hilbert norms is bounded
by 2^k times the word norm. This deliberately loose bound does not count
the internal dimension of an excited link space.

The same statement holds for a word containing an additional creator
\(\widehat u_J\): outside X the excited set is now
(J∪I_1∪…∪I_l)\X, or the word vanishes. This second version controls the
entire excitation operator, not only its action on Ω.

There is no factor 2^(total number of links) in either estimate.
At finite volume these expansions are finite; summing them under an
infinite exponential majorant is an upper bound, not a convergence
assumption about an infinite spatial volume.

## 4. Explicit contraction: every hidden cluster is returned

For ||c||_*≤r, (6), the 2^k projection bound and (9) give

\[
\boxed{\|\mathcal T(c)\|_*\le
\frac{2^k\beta}{3}(1+2r)e^{2kr}.}                         \tag{10}
\]

Here is the root count, including the term that removes the volume cost.
After |J| cancels against H_J inverse, an output containing a selected
link e must come from either e∈X or e∈I_a for at least one creator.
In the first case sum ||Φ_X|| over X containing e to get β, and sum the
remaining creators using (9), giving β exp(2kr). In the second case mark
one of the l creators containing e. For that I_a,

\[
\sum_{X:X\cap I_a\ne\varnothing}\|\Phi_X\|
\le\beta|I_a|.
\]

Its sum is at most βr by (7), while the other l−1 creators contribute
(kr)^(l−1). The sum of l·2^l/l! is 2 exp(2kr), giving the second term
2βr exp(2kr). Multiple markings only overcount and are safe.

For c,d in the same ball, telescope each ordered product by replacing
one creator at a time. If the replaced creator also carries the marked
root, its norm sum is ||c−d||_*; otherwise it contributes at most
k||c−d||_*. Summing the two cases gives

\[
\boxed{\|\mathcal T(c)-\mathcal T(d)\|_*\le
\frac{2^k\beta}{3}e^{2kr}
\big[2+2k(1+2r)\big]\|c-d\|_*.}                          \tag{11}
\]

Take k=4, r=1/16. The elementary rational estimate exp(1/2)<5/3
suffices. For example n!≥2·3^(n−2) for n≥2 implies
exp(1/2)≤3/2+(1/8)/(1−1/6)=33/20<5/3.
Equations (10)–(11) become

\[
\|\mathcal T(c)\|_*\le10\beta\le\frac1{16},\qquad
q:=\operatorname{Lip}(\mathcal T)\le\frac{880}{9}\beta
\le\frac{11}{18}<1                                        \tag{12}
\]

whenever β≤1/160. The finite direct sum in (7) is a Banach space even
though its component Hilbert spaces are infinite dimensional. Banach's
fixed-point theorem gives the full solution c and hence an exact
eigenvector Ψ. For domain control, each final c_I in (8) lies in Dom H_I;
there are finitely many subsets at fixed volume. Thus [H0,C] is bounded,
C and exp(±C) preserve Dom H0, and the displayed similarity identities
hold on that domain. A bounded perturbation defines the self-adjoint H.

The existence of Ψ as an eigenvector is not yet being used to assume it
is the ground. The spectral argument next proves this.

Starting at c^(0)=0, iteration has the certified local remainder

\[
\|c-c^{(m)}\|_*\le\frac{q^m}{1-q}
\|c^{(1)}-c^{(0)}\|_*\le\frac{q^m}{1-q}\frac\beta6.       \tag{13}
\]

The last estimate is specific to ordinary YM plaquettes: W_p belongs
entirely to I=∂p, H0 W_p=12 W_p and ||W_p||_2=1/2, so
c^(1)_∂p=t_p W_p/12 and its local weighted sum is β/6.
Formula (13) certifies convergence; this packet does not enumerate a
large-volume numerical vacuum or use a finite cluster truncation as exact.

## 5. The actual excitation operator has a relative local bound

Let \(\widehat H=e^{-C}(H-E)e^C\). It kills Ω. Acting on a nonempty
excitation u_J=û_J Ω gives

\[
\widehat H u_J=H_Ju_J+Uu_J,\qquad
Uu_J=[W(c),\widehat u_J]\Omega.                            \tag{14}
\]

The [H0,C] term commutes with every creator and cancels in this identity.
Also

\[
[W(c),\widehat u_J]
=\sum_X e^{-C}[\Phi_X,\widehat u_J]e^C.
\]

Only X meeting J contribute. All dressing creators must still meet X,
because they commute with û_J as well as one another. Applying the second
projection-count lemma, the extra commutator contributes a factor 2:

\[
\sum_I\|P_I Uu_J\|
\le2^{k+1}\beta |J|e^{2kr}\|u_J\|
\le\frac{160}{3}\beta |J|\|u_J\|.                         \tag{15}
\]

Use on the nonempty direct sum the norm ||u||_1=Σ_(I nonempty)||u_I||.
At fixed volume it is equivalent to the Hilbert norm; its equivalence
constant may grow with volume, but no bound below uses that constant.
The quotient by Ω has operator

\[
A=H_0+QUQ,\qquad Q=1-P_\varnothing,
\quad\|QUQH_0^{-1}\|_{1\to1}\le b:=\frac{160}{9}\beta
\le\frac19.                                               \tag{16}
\]

More generally restrict to a symmetry subspace preserved by H0, U and
all P_I, containing Ω, whose free nonempty bottom is δ. For 0≤z<δ,
the bound H_I≥3|I| and the restricted free energies λ≥δ imply

\[
\|QUQ(H_0-z)^{-1}\|_{1\to1}
\le\frac{b}{1-z/\delta}.                                 \tag{17}
\]

Indeed |I|/(λ−z)≤(1/3)λ/(λ−z)≤[3(1−z/δ)]^−1. For z<0 the
bound is at most b. Thus the Neumann inverse of
I+QUQ(H0−z)^−1 exists for every real z<δ(1−b).

In the decomposition Ω⊕QH, the full transformed operator is upper
triangular with diagonal blocks 0 and A. Its off-diagonal row is bounded
relative to H0 at fixed volume. The same inverse therefore excludes all
nonzero real spectrum below δ(1−b), while the kernel is exactly the
one-dimensional Ω line. Similarity preserves the spectrum of the
original self-adjoint operator. On the full product space δ=3, proving
that E is its unique ground energy and

\[
\Delta_{\rm full}\ge3(1-b)\ge\frac83.                    \tag{18}
\]

This removes the extensive scalar *and* controls the nonconstant operator.
No global estimate ||Φ||≤Σ_p|t_p| or ||R_YC8||≤Pε is substituted into
a gap difference. The price of a large excitation is paid by its free
energy 3|I|, matching its local interaction incidence.

## 6. Gauge and centre restriction raises the free floor to 12

Every p_e and q_e commutes with left and right SU(2) actions and with the
centre flip. H0, every plaquette interaction, and the fixed-point iteration
preserve gauge and global-centre symmetries. Starting at zero, each c_I
is invariant under the transformations acting on its supported links;
Ω_I is invariant too. Hence every ĉ_I, C, exp(±C) and U commutes with
these symmetries. The argument (17) applies to the vacuum sector itself,
not to a gauge quotient with a guessed metric.

For completeness its free nonconstant bottom is exactly δ_vac=12.
By Peter–Weyl, nontrivial link spin j has electric rate 4j(j+1)≥3.
Gauss invariance forbids a vertex incident to exactly one nontrivial
representation. A nonempty support therefore has no degree-one vertex.
On the ordinary torus graph with every N_i≥3, one or two distinct links
cannot meet this condition. A three-link support can only be a triangle;
in this graph such a triangle winds once around a direction of length 3.
Degree-two singlet constraints force the same spin on all three edges.
Its only energy below 12 would be j=1/2, at energy 9, and that winding
trace is odd under the corresponding global centre flip. It is excluded
from the vacuum sector. Integer j≥1 instead gives energy at least 24.
Every support with at least four links has energy at least 12.
Finally an elementary fundamental plaquette W_p is gauge invariant,
centre even, orthogonal to Ω, and has energy exactly 12.

Applying (17) with δ=12 proves (3), and β=4θ gives (1). The true
ground is gauge and centre invariant: its unique positive elliptic
ground function is unchanged by these commuting symmetries. Alternatively
the invariant fixed point already constructs it.

The graph argument covers all allowed sizes. The certificate checks its
short-cycle and centre-parity ingredients on representative finite boxes;
those samples are not the proof for arbitrarily many cells.

## 7. Continuity with YC8 and the remaining scale problem

The first cluster return gives Ψ=1+θS/12+O(θ²), consistent with YC8's
first logarithmic dressing. Its energy return is

\[
E_0=P\theta-\theta^2\langle S,H_0^{-1}S\rangle+O(\theta^3)
=P\theta-\frac{P}{48}\theta^2+O(\theta^3).
\]

YC8 remains the explicit second-order local-energy enclosure, on its
larger stated range. YC9 sums the entire creation-cluster response in a
small window and establishes an actual uniform excitation gap there.
It does not prove that YC8's comparison gap transfers through θ=1/2.
Nor does the creator C equal the multiplication operator −F_YC8.

Under YC6's declared physical normalization
H_phys=(γ g0²/a) H_(κ/g0⁴) plus a scalar, (1) supplies only a
strong-coupling lattice window κ/g0⁴≤1/640. The weak-coupling continuum
trajectory g0(a)→0 leaves it. The next quantitative task is to enlarge
the local convergence domain using correlated blocks / the YC8 dressing,
while preserving the local inverse and source bound at every join.
The continuum also needs scale matching and a constructed limit.

## 8. Replay, provenance and scope

```text
python physics/yc9/yc9_shared_plane_gap.py --check
python -m unittest discover -s physics/yc9 -p 'test_*.py' -v
```

The executable checks rational constants, local geometry, centre parity,
and finite exact creator/commutator witnesses. Tests separately exercise
an excited space of dimension greater than one, including an entangled
cluster. The analytic projection-count, contraction, domain and resolvent
proofs above supply the all-volume and untruncated claims. This is not
machine formal verification or independent peer review.

Primary context for the named classical methods:

- D. A. Yarotsky, *Quasi-particles in weak perturbations of non-interacting
  quantum lattice systems*, https://arxiv.org/abs/math-ph/0411042,
  section 2: creation operators, fixed-point dressing and relative
  spectral estimates, including infinite-dimensional local Hilbert spaces.
- S. Bravyi, D. DiVincenzo and D. Loss, *Polynomial-time algorithm for
  simulation of weakly interacting quantum spin systems*,
  https://arxiv.org/abs/0707.1894: explicit Kirkwood–Thomas estimates for
  bounded-degree qubit interactions. Its qubit/2-local theorem is not
  silently applied to SU(2); sections 2–6 above give the needed proof.
- D. A. Yarotsky, *Ground states in relatively bounded quantum
  perturbations of classical lattice systems*,
  https://arxiv.org/abs/math-ph/0412040: broader stability background.
  No unknown smallness constant from that theorem is used in (1).

The owner's shared-plane intuition chooses the overlap geometry and the
local budget. The kinetic metric, gauge constraints, full excitation
inverse and contraction establish the spectral statement. Each role is
kept explicit.
