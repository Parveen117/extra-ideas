# YC22 — symmetry protects the complete charged return

10 October 2026 (India). Continues YC21 at `c45bc50`. Earlier packets
remain frozen. Electric normalization is the unit-S³ Laplacian, Wp=Tr(Up)/2.

**Result.** A ground-state form inequality gives a coupling-independent
floor for every non-Gauss-invariant sector of a finite SU(2) link graph.
It also proves which boundary sources remain orthogonal under the entire
interacting block inverse. This removes a weak charged-sector estimate from
the cube/block join; it does not replace the physical vacuum problem by
the charged one.

Combining this with YC13's actual gauge gap and YC21's complete bridge cut:

- A twelve-link cube has **full gap >=1** for all theta in [0,2].
- A 28-link two-cube block has **full gap >=3/4**, cube couplings <=2,
  and all four joining couplings <=1. Its **gauge gap is >=17/10**.
- If its cube couplings are <=1, its gauge gap is **>=15/8**.
- At joining couplings <=1/2, a separate full-space count improves
  YC21's 1/5 to **9/10**, now allowing cube couplings <=2.

The actual finite-lattice joins are uniform over all the declared volumes:

| Reference factors | Internal couplings | Remaining external cap | Full gap | Physical vacuum gap |
|---|---|---:|---:|---:|
| Cubes | cube theta<=2 | 1/1920 | 67/75 | 603/125 |
| Two-cube blocks | cube theta<=2, tube eta<=1 | 1/4480 | 461/700 | 7837/5250 |
| Two-cube blocks, resolved boundary forces | cube theta<=1, tube eta<=1 | 1/4000 | 1619/2500 | 1619/1000 |

Cube tilings require even side lengths >=4; two-cube tilings require Lx a
multiple of4 >=8 and Ly,Lz even >=4. All couplings in this table are
nonnegative. Full link harmonics, all boundary charges and every cluster
order remain. These are anisotropic finite-lattice results. YC10's larger
isotropic window remains the isotropic baseline. No continuum theorem is
claimed, and the improved reference does not yet give repeated blocking
contraction.

## 1. A symmetry cut whose channels really do not mix

The owner proposes perpendicular information readings with no overlap or
mixing. The relevant exact condition is an invariant channel decomposition,
not orthogonality at just one input vector. Here the symmetry supplies it.

Let Gamma be a finite connected simple graph with at least one edge, and
put H=−sum_e Delta_e+V on L²(SU(2)^E,dmu). The real smooth potential V
is invariant under the full vertex gauge group SU(2)^vertices. There are
no matter degrees of freedom. The compact elliptic Hamiltonian has a
unique positive normalized ground Omega, necessarily gauge invariant.
Write h=H−E0 and dnu=Omega²dmu. Multiplication by Omega is a unitary
from L²(dnu) to the original Hilbert space, and its exact form is

\[
\langle\Omega f,h\Omega f\rangle
=\mathcal E_\nu(f):=\sum_e\int|\nabla_e f|^2d\nu.       \tag{1}
\]

All derivatives are first checked on smooth functions and extend by form
closure. Smooth positivity on this finite compact carrier suffices; no
volume-uniform pointwise bound on Omega is assumed.

Let P_G be Haar averaging over the vertex gauge group. This is an exact
orthogonal projection in both pairings. K_G=2P_G−I is a represented cut,
and [H,K_G]=0. Thus Q_G H P_G=0 and the charged complement creates no
Schur source in the physical sector. This is a concrete realization of
the owner's no-mixing condition, not a conclusion that all physical
closed-loop channels are independent.

## 2. Gauge-orbit energy from the ground-state form

At vertex v of degree d_v, let G_v^a be the infinitesimal gauge generator:
the signed sum of the left/right link vector fields incident to v, in the
same normalization as −Delta. Define C_v=−sum_a(G_v^a)². With n_v=2j_v,
its irreducible eigenvalues are n_v(n_v+2). In particular the fundamental
has Casimir3 and the adjoint has Casimir8.

Pointwise Cauchy–Schwarz and left/right equality of the single-link metric
give

\[
\sum_a|G_v^a f|^2\le d_v\sum_{e\ni v}|\nabla_e f|^2.
\]

Since nu is gauge invariant, the generators are skew-adjoint in nu.
Each edge belongs to two stars; after integration,

\[
\boxed{h\succeq\frac12\sum_v\frac{C_v}{d_v},
\qquad h\succeq\frac{C_v}{d_v}\quad\hbox{for each fixed }v.} \tag{2}
\]

These are form inequalities in the original Hilbert space too, because
Omega commutes with gauge action. The second inequality is separate; it
must not be summed without the factor1/2 in the first one.

All simultaneous vertex-centre flips act trivially on pure link variables.
Consequently every occurring charge pattern obeys sum_v n_v even. A
nonzero allowed pattern has either two or more odd n_v, costing at least
3+3, or a nonzero even n_v, costing at least8. With dmax=max d_v,

\[
\boxed{h|_{\operatorname{ran}Q_G}\succeq(3/d_{\max})I.} \tag{3}
\]

This is independent of V's size and of vertex count. It is a spectral
floor **above the true ground**, not a free-kinetic estimate with an
uncontrolled potential subtraction. It applies, for example, to charged
sectors on ordinary three-dimensional tori with dmax=6, giving1/2 at
every finite coupling. The Casimirs vanish on the physical carrier, so
(2) supplies zero there. Claiming a physical mass gap from (3) would be
invalid.

For constant positive electric coefficients a_e the same proof uses
kappa_v=sum_(e incident v) a_e^−1 in place of d_v. This is a weighted
link-product form statement; it does not automatically cover an arbitrary
nonlocal metric after blocking.

## 3. What survives the full interacting inverse

For an oriented edge e=(u,v), let u_e,alpha, alpha=0,...,3, be its quaternion
coordinates and xi_e,alpha=u_e,alpha Omega. Their charge pattern is
fundamental at u and v and trivial elsewhere. Different edges of a simple
graph have different endpoint patterns. Functional calculus of h preserves
each full isotypic carrier, so for any bounded Borel b,

\[
\boxed{\langle\xi_{e,\alpha},b(h)\xi_{f,\beta}\rangle=0
       \quad(e\ne f).}                                \tag{4}
\]

For a fixed e, SU(2)_u×SU(2)_v irreducibility gives a scalar response in
its four coordinate indices. Its common norm is1/4, since a single-edge
marginal of a gauge-invariant probability is Haar. The source mean is3:
the ground-form identity gives sum_alpha <xi_alpha,h xi_alpha>=3, and
symmetry divides it equally. Resolvents may leave the four source vectors
and enter infinitely many copies of the same representation. No invariant
four-dimensional source span is assumed.

Thus the full inverse preserves orthogonality **between distinct charge
patterns**, even for strong internal interactions. Inside a repeated
representation it can mix vectors. In particular all gauge-invariant
Wilson-loop readings have charge pattern zero and can mix freely. The
finite certificate includes both the protected and same-charge mixing
controls. External plaquettes also exchange the individual reference
blocks' boundary charges: (4) cannot be used to delete those later cluster
returns. The full join below retains them.

## 4. A local force bound for the actual boundary variance

Let L_e^a be a left link generator, and let u be that endpoint. Because
Omega is gauge invariant, the three functions L_e^a Omega transform as
one adjoint at u and as singlets elsewhere. The second inequality (2)
therefore gives an energy floor8/d_u on this carrier. The exact ground
equation, and [−Delta,L_e^a]=0, give

\[
h L_e^a\Omega=-(L_e^a V)\Omega.
\]

Use the complete inverse on this adjoint carrier and sum the three
components. The right generators give the same gradient norm and allow
the better endpoint degree:

\[
\boxed{T_e:=\int|\nabla_e\Omega|^2d\mu
\le\frac{\min(d_u,d_v)^2}{64}
      \int|\nabla_e V|^2\Omega^2d\mu.}                 \tag{5}
\]

For V=−sum_p t_p Wp, |grad_e Wp|²=1−Wp²<=1 on an ordinary plaquette.
Hence, with c_e=sum_(p containing e)|t_p|,

\[
T_e\le\min(d_u,d_v)^2 c_e^2/64.                        \tag{6}
\]

This concerns the true interacting ground and is local in its incident
couplings. It does not need a ground-vector truncation or the total number
of plaquettes. It need not beat YC15's special symmetric-cube estimate;
its use here is the nonsymmetric rectangular block.

For an external plaquette joining two true blocks and two Haar bridges,
YC15's exact conditional Haar calculation still gives
||r||²=1/4, <r,H_ref r>=3 and

\[
s=(H_{\rm ref}-12)r,\quad\langle r,s\rangle=0,
\quad\|s\|^2=T_{A,e_A}+T_{B,e_B}.                      \tag{7}
\]

If the complete source floor is d>0, the unchanged full residual bound is

\[
\|H_{\rm ref}^{-1}r\|^2
\le\frac1{576}+\|s\|^2
       \left(\frac1{864d}+\frac1{144d^2}\right).         \tag{8}
\]

This retains all electric harmonics and all hidden propagation.

## 5. Actual cube and stronger tube blocks

The cube has dmax=3. Its charged floor is1 for every gauge-invariant
potential. YC13's frozen full gauge certificate gives Delta_G>=27/5 for
theta<=2 (and>=11 for theta<=1). The unique true ground is gauge
invariant, so the orthogonal union of those two estimates proves

\[
\boxed{\Delta_{\rm cube,full}\ge1\quad(0\le\theta\le2).} \tag{9}
\]

For the 28-link block, dmax=4 gives charged floor3/4 at every coupling.
On its **gauge-invariant** carrier the simultaneous flip of all four
joining bridges is a gauge transformation: flip every vertex of the left
cube by the centre. Thus the odd bridge sector is absent here. The retained
bridge-vacuum carrier consists of the two separate gauge-invariant cubes.
Its first nonzero energy is at least27/5, or11 when theta_A,theta_B<=1.

YC21's entire even hidden inverse and source identity remain valid for
eta_p<=1: d=6−sum eta_p>=2 and R†R=(sum eta_p²/4)I<=I. In the physical
Schur pencil, the codimension-one reserve is therefore

\[
\delta_{\rm cube,G}-z-\frac1{2-z}.
\]

At z=17/10 and delta_cube,G=27/5 it equals11/30. At z=15/8 and
delta_cube,G=11 it equals9/8. The actual ground expectation in the
reference is zero, so both are actual gap lower bounds. Combine them with
the independent charged floor3/4 to obtain the stated full block gap.
No false positive odd-hidden floor is used at eta=1.

At eta<=1/2 the full retained cube gap1 can instead be used in YC21's
original full-carrier count: d>=4, v<=1/4, odd floor>=1, and z=9/10
has reserve3/155. This proves the stronger full gap9/10 in that subwindow.

In the enlarged eta<=1 window the physical graph metric still has a full
bound I<=M(z)<=13I/9 for z<=1/2. YC21's tighter50/49 bound remains valid
in its original smaller parameter window. A metric bound alone is not
substituted for any of the gap counts above.

## 6. Boundary costs and the actual lattice joins

All joins use the true positive ground of each reference factor. Put
r=1/128, exp(8r)<16/15, k=M_source=4 and retain nonlinear support s_min=1.
If each full factor has gap>=g and the maximal incident interaction budget
is beta, the established complete creator bounds are

\[
q\le\frac{4\beta}{g}\frac{16}{15}
           [2+8(1+2r)],\qquad
b\le\frac{8\beta}{g}\frac{16}{15}.                    \tag{10}
\]

The source-centred map costs seed+qr; gap_full>=g(1−b). These are
YC10/YC14's full-return estimates, not a new truncation.

For cubes, g=1, beta<=24epsilon. A two-bridge source has floor8 and norm1/2,
so inverse norm<=1/16 and seed<=6epsilon. At epsilon=1/1920,
q<=27/50, b<=8/75 and ball reserve3/6400.

For strong two-cube blocks, g=3/4, beta<=48epsilon. The source floor is
15/2, giving inverse norm<=1/15 and seed<=64epsilon/5. At epsilon=1/4480,
q<=108/175, b<=64/525 and ball reserve3/22400.

There is an additional whole-window improvement when theta_A,theta_B<=1,
eta_p<=1. Enumerate the block's external face incidences by their internal
edge type. Matching edges on neighbouring identical rectangular tilings
have the same type, even when their couplings differ:

| Boundary edge class | Incidences per block | T_e upper from (6) | Inverse norm² upper from (8) | Outward inverse norm cap |
|---|---:|---:|---:|---:|
| Endpoint degree minimum3, two incident faces | 32 | 9/16 | 59/28800 | 91/2000 |
| Middle joining edge, degree4, two faces | 8 | 1 | 11/4800 | 6/125 |
| Inner transverse edge, degree4, three faces | 8 | 9/4 | 43/14400 | 11/200 |

Four-bridge sources retain inverse norm1/24. A single remaining bridge has
at most four external incidences, whose seed cost is less than the block
cost. Thus

\[
\mathrm{seed}\le4\epsilon
 [32(91/2000)+8(6/125)+8(11/200)]=(228/25)\epsilon.
\]

At epsilon=1/4000, q<=432/625, b<=256/1875 and ball reserve53/400000.
All endpoint bounds are monotone majorants, proving their whole coupling
windows. The no-mixing statement is used where symmetry proves it;
the existing sum norm and complete nonlinear return are kept elsewhere.

## 7. The physical free floor is kept separate from the charged floor

In each reference product, block ground projections and bridge electric
projections commute with global Gauss action. For physical states with no
excited bridges, each excited block is separately gauge invariant, so the
reference floor is its gauge gap delta_G.

If a bridge has an odd harmonic, centre invariance for the gauge transform
that flips all vertices in one block requires an even number of odd bridge
ends at that block. One odd bridge alone is impossible. Two odd bridges
cost at least6; a nonconstant even bridge costs at least8. Thus any physical
reference state with excited bridges has bridge energy at least6. Parallel
coarse edges are allowed: do not assume four distinct coarse edges are
necessary. Therefore the reference physical floor is at least

\[
\delta=\min(\delta_G,6).
\]

The full creator fixed point preserves Gauss and global centre symmetry.
The energy-relative excitation estimate (10) still uses **g**, the full
factor floor. On the restricted physical carrier, YC9's exact relative
resolvent argument then gives

\[
\boxed{\Delta_{\rm vac}\ge\delta(1-b).}               \tag{11}
\]

Using delta=27/5,17/10,15/8 gives the three physical columns in the opening
table. This does not replace g by delta in the inverse/source/contraction
budget. All later charged exchanges and all neutral multiplicity mixing
are retained. The original operator, not a ground-response approximation,
is bounded, uniformly in the declared finite volumes.

## 8. Continuum gate and evidence

This stage supplies an all-coupling, size-independent **charged** floor and
larger **physical** anisotropic lattice windows. The floor in (3) vanishes
as an inequality on the physical subspace because all C_v are zero there.
Closed-loop interactions can mix inside that same neutral carrier; they
are the remaining physical problem. Orthogonal initial readings do not
justify replacing the full neutral inverse by a diagonal one.

The current blocking cost has not become contractive: incidence/g changes
from24/1 for cubes to48/(3/4) for the strong rectangular factors. A new
scale step needs control of the neutral connected interaction and its full
excitation inverse. A physical continuum trajectory also needs the electric
normalization, running coupling, spatial limit and construction of continuum
observables. None follows by renaming this lattice gap a mass.

Classical ingredients are the ground-state transform, Cauchy–Schwarz,
Peter–Weyl charge decomposition, Schur's lemma, full Feshbach count and the
previously proved creator contraction. The normalization and exact gauge
action follow the Hamiltonian link formulation of
[Kogut and Susskind (1975)](https://doi.org/10.1103/PhysRevD.11.395).
The new inequalities are proved above rather than attributed to that paper.

There are 34 exact checks and 15 new tests; the 14 YC21, 14 YC13 and15 YC15
tests also pass under Python 3.12.14. The executable packet checks quaternion generator normalization, a direct
weighted ground-transform witness, charge/parity floors, full-source inverse
controls, graph incidence and all rational reserves. Tests include mixing
within one charge pattern and rejection of a charged-to-physical inference.
Finite algebra checks support the written full-carrier arguments; they are
not a harmonic truncation or formal proof-assistant verification.

```sh
python physics/yc22/yc22_gauge_protected_channels.py --check
python -m unittest discover -s physics/yc22 -p 'test_*.py'
```
