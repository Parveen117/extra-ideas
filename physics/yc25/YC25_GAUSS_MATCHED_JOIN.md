# YC25 — Gauss matching pays for the complete spatial join

10 October 2026. Research owner: Monty Dabas. Continues YC24 at `6552e55`
and the construction policy at `7e07ed6`. Unit-S³ electric normalization,
normalized Haar measure and Wp=Tr(Up)/2. Frozen predecessor packets remain
unchanged. The result concerns the actual lattice Hamiltonian.

**Result.** On physical excitation supports the inverse can use a stronger
*energy per excited factor*, because a charged reference block must be
matched by excited bridges. The full, symmetry-preserving creator family
and its excitation operator therefore admit larger spatial join windows:

| Reference | Internal couplings | Previous external cap | New external cap | Physical gap at new cap |
|---|---|---:|---:|---:|
| Twelve-link cubes | cube theta<=2 | 1/1920 | **1/840** | **5643/1225** |
| 28-link rectangles | cube theta<=2, tube eta<=1 | 1/4480 | **1/1792** | **209/105** |
| 28-link rectangles, resolved boundary forces | cube theta<=1, tube eta<=1 | 1/4000 | **1/1680** | **12396/6125** |

The last cap is 50/21 times its predecessor. On the old external window
epsilon<=1/4000 its new physical floor is **49044/21875**, improving
YC24's6476/3125. These are all-volume anisotropic lattice bounds. The
nonphysical complement has the separate full-lattice charged floor1/2.
No repeated blocking contraction or 4D continuum limit is claimed.

## 1. Construct the family before its zero-reading

Use precisely YC22/YC24's tilings: cubes on even tori with sides>=4;
rectangles on Lx divisible by4 and>=8, with Ly,Lz even and>=4. Every
vertex belongs to one reference block. Remaining links are individual
bridge factors between distinct blocks. Retain each block's full Hilbert
space, including all boundary charges, and every bridge harmonic.

Let Omega_b be each actual positive block ground, h_b=H_b-E_b, and

\[
 \Omega=\bigotimes_b\Omega_b\otimes\bigotimes_e1,
 \qquad H_0=\sum_b h_b+\sum_{e\ {\rm bridge}}(-\Delta_e).
\]

For a fixed external profile 0<=xi_p<=1, construct the full family

\[
 H_\lambda=H_0+\lambda\Phi,\qquad
 \Phi=-\sum_{p\ {\rm external}}\xi_pW_p.                 \tag{1}
\]

The omitted Wilson constants and block ground energies are scalar shifts.
Here lambda is the **declared interface coordinate of this YM adapter**.
It is not an identification of the entire primitive lambda-space with a
real interval, an observation-reset parameter, tower index or RG scale.
The complete creator family below exists holomorphically in a complex
disk and gives the physical family on its real interval. No derivative
truncation defines that family.

The reference pairing is the full Haar Hilbert pairing. Denote the exact
factor excitation projections by P_I, and their number operator by
N=sum_i(1-|Omega_i><Omega_i|). Thus N=|I| on P_I. Each P_I and H_0
commutes with the full vertex gauge action. They need not split the
physical carrier into a tensor product of independent block singlets.

## 2. The missing energy count

**G1 — matched physical support floor.** Suppose every reference block
has charged floor g and neutral excited floor delta, with
0<g<=3 and delta>=g. On the globally gauge-invariant carrier,

\[
 \boxed{H_0\succeq\mu N,\qquad
        \mu=\min\{\delta,(g+3)/2\}.}                   \tag{2}
\]

The estimate is extensive in the number of excited *factors*, not a
replacement of a full-factor gap by the physical one-level gap.

**Proof.** Decompose into each block's vertex-gauge isotypic sectors and
each bridge's full Peter–Weyl harmonic degree. These decompositions
commute with H_0, all factor excitation cuts and global gauge averaging.
They retain all multiplicities and internal excited energies. In one
such component let

- n0 count neutral excited blocks;
- nc count charged blocks;
- bo count bridges of odd harmonic degree;
- be count nonconstant bridges of even harmonic degree.

Then N=n0+nc+bo+be and

\[
 H_0\succeq \delta n_0+g n_c+3b_o+8b_e.                \tag{3}
\]

Flip the gauge centre simultaneously at every vertex of one block.
Every internal link is unchanged. On a physical state this forces an
even number of incident odd bridge ends at that block. Consequently a
charged block incident to an odd bridge has at least two such ends.
There are 2bo odd bridge ends in total, so at most bo charged blocks
can be of this kind.

A charged block without odd bridge ends needs at least one nonconstant
even bridge: if all its bridges were Haar constants, gauge invariance
would require its internal state to be neutral at every vertex. There
are 2be even bridge ends, so at most 2be charged blocks are of this
second kind. Thus the exact resource inequality is

\[
 \boxed{n_c\le b_o+2b_e.}                              \tag{4}
\]

It permits parallel coarse bridges, repeated representations and arbitrary
neutral loops; no minimum four-edge coarse cycle is assumed. Because
g<=mu<=3 and mu<=delta, subtracting mu N from(3), then using(4), gives

\[
 H_0-\mu N\succeq
 (\delta-\mu)n_0+(g+3-2\mu)b_o+(2g+8-3\mu)b_e\succeq0.
\]

The last even-bridge coefficient is at least(g+7)/2. Summing the
orthogonal components and closing the positive forms proves(2). QED.

YC22 gives g=1 for cubes and g=3/4 for rectangles. Their physical block
floors are27/5,7/3,12/5 in the three table rows, respectively. Hence

\[
 \boxed{\mu_{\rm cube}=2,\qquad
        \mu_{\rm rectangle}=15/8.}                     \tag{5}
\]

For every nonempty *physical* support I, therefore,
||H_I^-1||<=1/(mu|I|). The independent physical reference gap remains
delta_ref=min(delta,6), as in YC22; here it equals the displayed block
floor in each row. Delta_ref controls the first physical level, whereas
mu controls the cost per factor of arbitrarily large supports.

## 3. Apply the stronger inverse only where symmetry allows it

Let c_I be globally invariant excitation vectors, including the
centre-even restriction for the vacuum sector, and define

\[
 \widehat c_I=|c_I\rangle\langle\Omega_I|,\quad
 C=\sum_{I\ne\varnothing}\widehat c_I,\quad
 \|c\|_* =\max_i\sum_{I\ni i}|I|\|c_I\|.
\]

Every creator commutes with global gauge transformations: both its
source vacuum and its target vector are invariant on its support.
Each individual plaquette, excitation projection and H_I inverse also
preserves that symmetry. The complete map

\[
 \boxed{\mathcal T_\lambda(c)_I
 =-\lambda H_I^{-1}P_Ie^{-C}\Phi e^C\Omega}             \tag{6}
\]

therefore acts on the closed invariant creator space. Every word to
which the inverse is applied lies in that space. This proves the use
of(2) throughout the nonlinear map, not just on its first source.

YC9/YC10's root count now applies with mu in the inverse denominator.
The support size is k=4, the orthogonal output factor is sqrt(2^4)=4,
and the minimum nonlinear support remains1: a single neutral excited
block is allowed. Write beta=max_i sum_(p incident i)|lambda xi_p|.
For ||c||_*<=r,

\[
 q\le\frac{4\beta}{\mu}e^{8r}[2+8(1+2r)],\qquad
 \|\mathcal T_\lambda(c)\|_*\le s_0+qr.                \tag{7}
\]

The unchanged exact first-source estimates of YC22 give

| Reference row | beta upper | s0 upper |
|---|---:|---:|
| Cubes | 24 abs(lambda) | 6 abs(lambda) |
| Rectangles, cube theta<=2 | 48 abs(lambda) | (64/5) abs(lambda) |
| Rectangles, cube theta<=1 | 48 abs(lambda) | (228/25) abs(lambda) |

In particular no later source is assumed to retain the first source's
four-factor support. All cluster orders, spectators, charged exchanges
and neutral multiplicity mixing remain in(6).

The same improvement applies to the **physical excitation quotient**.
For physical u_J, H_0^-1 u_J is still physical and has inverse bound
1/(mu|J|). The full commutator/projection estimate is unchanged:

\[
 \sum_I\|P_I Uu_J\|\le8\beta e^{8r}|J|\|u_J\|,
 \qquad
 \boxed{b_G:=\|QUQH_0^{-1}\|_{1\to1,G}
                  \le\frac{8\beta}{\mu}e^{8r}.}        \tag{8}
\]

Here Uu_J=[e^-C(lambda Phi)e^C,widehat u_J]Omega and Q removes the
reference vacuum. Physical P_J components remain orthogonal and invariant,
so the component-sum norm proof of YC9 applies without change. For real
z below delta_ref(1-b_G), the physical reference resolvent satisfies

\[
 \|H_0(H_0-z)^{-1}\|_{1\to1,G}
 \le\begin{cases}1,&z\le0,\\
 (1-z/\delta_{\rm ref})^{-1},&0<z<\delta_{\rm ref}.
 \end{cases}
\]

Thus the relative Neumann argument excludes every nonzero physical
eigenvalue of the dressed operator below delta_ref(1-b_G). The zero
eigenvalue is simple. This establishes

\[
 \boxed{\Delta_G(H_\lambda)\ge\delta_{\rm ref}(1-b_G).} \tag{9}
\]

For completeness, the actual compact elliptic Hamiltonian has a unique
positive ground, invariant under all gauge and centre symmetries. Its
lowest eigenvalue is therefore in this same physical vacuum sector.
The branch constructed by(6) is that ground, since the physical quotient
has no lower level. We do not need to apply(2) to nonphysical vectors.
YC22 separately gives charged gap>=3/6=1/2 on the full torus above this
true ground. Full gap>=min(1/2,delta_ref(1-b_G)) follows.

At finite volume there are finitely many supports, each with its complete
infinite-dimensional excited carrier. The fixed point belongs to Dom H_I;
bounded creators and their commutators preserve Dom H_0 as in YC9. All
similarity and resolvent assertions use that inherited domain argument.

## 4. Exact new windows

Choose r=1/32. The elementary exponential majorant is

\[
 e^{1/4}\le\sum_{j=0}^3\frac{(1/4)^j}{j!}
 +\frac{(1/4)^4}{4!}\frac1{1-1/20}
 =\frac{12491}{9728}<\frac97.                           \tag{10}
\]

Insert it into(7),(8). At the three new endpoints:

| External cap | q upper | b_G upper | r-s0-qr lower | Physical gap lower |
|---:|---:|---:|---:|---:|
| 1/840 | 27/35 | 36/245 | 0 | 5643/1225 |
| 1/1792 | 27/35 | 36/245 | 0 | 209/105 |
| 1/1680 | 144/175 | 192/1225 | 3/28000 | 12396/6125 |

A closed-ball self-map with q<1 suffices; the zero displayed margins
in the first two rows are also strictly positive with(10)'s sharper
rational value. The estimates are monotone in abs(lambda), so every
smaller coupling and each allowed nonuniform profile is included.

For rectangles the whole-window physical relative bound is

\[
 b_G\le\frac{9216}{35}|\lambda|.                       \tag{11}
\]

At theta,eta<=1 and lambda<=1/4000, equation(9) becomes
(12/5)(1-288/4375)=49044/21875. At the new endpoint1/1680 it is
12396/6125, approximately2.02384. The old stronger full-space bounds
on their old domains are not withdrawn; this proof's independent
charged estimate1/2 is simply less sharp there.

## 5. Zero-reading, complete history and metric

Uniform contraction in each complex disk of the above radius constructs
the unique holomorphic fixed-point family c(lambda) in the declared ball.
It is the uniform limit of the holomorphic iterates of(6). This is a
construction of the complete interface tower, with finite-volume domain
control and volume-uniform local-norm constants; it is not a finite
Taylor polynomial or an infinite-volume state construction.

At lambda=0, equation(6) is identically zero, giving exactly

\[
 c(0)=0,\quad S(0)=e^{C(0)}=I,\quad
 \Psi(0)=\Omega,\quad E(0)=0,\quad A(0)=H_0|_Q.         \tag{12}
\]

Here E is measured relative to the omitted reference scalar, and A is
the exact excitation quotient of S^-1(H_lambda-E)S. The first jet remains

\[
 c'_I(0)=-H_I^{-1}P_I\Phi\Omega,\qquad E'(0)=0,         \tag{13}
\]

with complete source inverses. The zero-value reading does not erase this
retained derivative/history. For real lambda carry the exact pairing

\[
 M_\lambda=S^\dagger S,\qquad
 G_\lambda=M_{QQ}-M_{Q0}M_{00}^{-1}M_{0Q}.              \tag{14}
\]

As in YC19, A^dagger G=GA on the excitation quotient. At zero both M and G
are identity on their respective carriers. The maps in

\[
 \boxed{\pi_0(\operatorname{Construct}_\lambda)
              =\operatorname{Construct}_0}             \tag{15}
\]

are evaluation at zero on this operator/state/metric family; equations
(6),(12),(14) prove compatibility. The zero model retains the interacting
blocks. It is the established decoupled-interface model, not a claim
that all of Yang–Mills becomes classical or that internal curvature
vanishes. No volume-uniform condition number of S or G is asserted.

## 6. What changed, and what the next scale still needs

The charged block and its compensating bridge were previously bounded
as separate possible excitations. On a physical source they are linked
by the gauge cut. Returning that constraint to the energy count raises
the extensive inverse floor from1 to2 for cubes and from3/4 to15/8
for rectangles. This controls the complete source and excitation towers,
and enlarges their actual spatial join intervals.

The new scale comparison is nevertheless24/2=12 for cubes against
48/(15/8)=128/5 for rectangles. Both are improved budgets; the second
still exceeds the first. This stage therefore does not manufacture a
contracting sequence by comparing one improved level with an old coarse
bound. Next is a second block update retaining this matching constraint
and controlling its connected neutral return. A physical continuum
trajectory and the continuum observable construction are still needed.

The classical analytic mechanism is Kirkwood–Thomas creation dressing
and a relative Neumann resolvent estimate. Background:
[Yarotsky, Quasi-particles in weak perturbations of non-interacting quantum
lattice systems](https://arxiv.org/abs/math-ph/0411042). The concrete
matched-support inequality and rational windows are derived here from
YC22/YC24, not attributed to that reference.

## 7. Minimal verification and source route

The all-support proof is(3)–(4); the all-order proof is the invariant
fixed point and relative resolvent above. Zero specialization is proved
in(12)–(15). The small calculator `yc25_join_bounds.py` reproduces only
the necessary exact rational budgets at zero and the stated endpoints:

```sh
python physics/yc25/yc25_join_bounds.py
```

It is arithmetic for this proof, not a numerical lattice spectrum,
independent formal verification or a replacement proof of the source
theorems. No new unit-test suite, frozen certificate regeneration or
predecessor regression run accompanies this documentation-led stage.

Read the frozen inputs at commit `6552e55`:

- `physics/yc9/YC9_SHARED_PLANE_GAP.md`, sections2–6: complete creator,
  invariant-domain and relative-resolvent construction.
- `physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md`, sections2–5:
  invariant creator restriction and orthogonal projection factor4.
- `physics/yc22/YC22_GAUGE_PROTECTED_CHANNELS.md`, sections2,5–7:
  true-ground charged floors, exact first-source budgets, tilings and
  physical reference gap.
- `physics/yc24/YC24_PARITY_RESOLVED_RESPONSE_TOWER.md`, section5:
  full-harmonic physical rectangle floors7/3 and12/5.
- `physics/yc19/YC19_LOCAL_EXCITATION_OPERATOR.md`: the complete
  similarity and physical excitation-quotient metric.
