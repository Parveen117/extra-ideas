# YC26 — the next spatial block retains the charge carried by its bridges

10 October 2026. Research owner: Monty Dabas. Continues YC25 at `759f8f5`.
Unit-S³ electric normalization, normalized Haar measure, Wp=Tr(Up)/2.
All predecessor proof and evidence packets remain unchanged.

**Result.** Gauss matching improves the complete hidden inverse before
the next spatial join. The actual28-link rectangle now has physical gap
at least **4** for cube theta<=2, or **9/2** for theta<=1, with every tube
coupling<=1. Its complete physical returned metric is at most26I/25
through z<=1/2.

Join two such rectangles in the perpendicular y direction. The resulting
4×4×2-vertex block has **64 links, 42 faces, eight new bridges and ten
new joining faces**. The entire physical hidden carrier is retained:

| Cube cap | First x-tube cap | New y-join cap | Physical64-link gap lower bound |
|---:|---:|---:|---:|
| 2 | 1 | 1 | **1** |
| 1 | 1 | 1 | **11/10** |
| 2 | 1 | 3/4 | **5/2** |
| 2 | 1 | 1/2 | **7/2** |
| 1 | 1 | 1/2 | **4** |

Each row permits nonuniform nonnegative joining profiles. At full new
join strength the returned metric obeys M(0)<=10I/7 and M(z)<=101I/50
through z<=1/2. The charged complement has the separate floor3/5 above
the true block ground. These are original-operator bounds, not a
replacement by finitely many electric harmonics.

All improved Schur and metric bounds in this stage refer to the globally
gauge-invariant retained/hidden carriers. Hidden parents may carry matched
nonzero charges; they are not independently restricted to singlets.

The stronger28-link parent bounds also improve YC25's existing lattice
floors without changing its external windows:

| Cube/tube caps | External cap | YC25 physical lattice floor | New physical lattice floor |
|---|---:|---:|---:|
| theta<=2, tube<=1 | 1/1792 | 209/105 | **836/245** |
| theta<=1, tube<=1 | 1/1680 | 12396/6125 | **9297/2450** |

The last floor is approximately3.79469. This stage supplies a perpendicular
second block and its complete return. Its remaining boundary budget still
does not contract, as quantified in section7.

## 1. Operator family and the two zero-readings

The owner's latest clarification treats lambda-space as an operator space
and leaves the mathematical construction choices to the researcher.
Here the admitted YM representation is explicit: the full64-link Hilbert
carrier, the inherited Haar pairing, the elliptic electric domains, the
vertex gauge action and the following operator-valued family. The two
interface coordinates are coordinates of this adapter, not a claim to
have constructed or identified the universal primitive operator space.

Start with four disjoint cubes, the eight bridges making two x rectangles,
and the eight bridges joining these rectangles in y. With fixed profiles
0<=xi_p,zeta_p<=1, the unshifted joint family is

\[
 H(t_x,t_y)=\sum_{c=1}^4h_c+\sum_{e\ {\rm new}}(-\Delta_e)
       -t_x\sum_{p\in X}\xi_pW_p-t_y\sum_{p\in Y}\zeta_pW_p,
 \quad |X|=8,\quad |Y|=10.                              \tag{1}
\]

There are16 new bridge links in total. Internal cube couplings are fixed
in their declared windows. Scalar Wilson constants can be restored at
any stage; they change no gap. The family is defined on the full inherited
domain before any specialization or channel elimination.

Evaluation at t_y=0 recovers the two complete x blocks and the eight free
y bridges. Setting t_x=0 too recovers the four cubes and sixteen free
bridges. The two parameter evaluations commute on(1). They do not erase
the corresponding source derivatives or set the internal cube couplings
to zero. Actual block grounds are used when shifting an intermediate
block to zero energy; these scalar shifts are carried separately.

## 2. A bridge carries its endpoint block charge

**J1 — the one-new-link-per-vertex gluing law.** Join two disjoint gauge
blocks A,B by bridges such that each interface vertex meets exactly one
new bridge. Let d_A(v),d_B(w) be its degrees inside the old blocks, and
let h_A,h_B be their full Hamiltonians above their true grounds. On the
gauge-invariant carrier of the joined graph,

\[
 \boxed{h_A+h_B+\sum_e(-\Delta_e)
  \succeq\sum_{e=(v,w)} a_e(-\Delta_e),\qquad
  a_e=1+\frac1{2d_A(v)}+\frac1{2d_B(w)}.}                \tag{2}
\]

**Proof.** YC22's true-ground form gives
h_b >= (1/2)sum_v C_(b,v)/d_b(v). At an interface vertex the physical
Gauss equation is (G_(b,v)+G_(e,v))psi=0. The generators on block and
bridge factors commute, so their squared Casimirs agree on that
physical carrier: C_(b,v)psi=C_(e,v)psi. In unit-S³ normalization the
bridge Casimir is exactly -Delta_e, at either endpoint. Keep the
nonnegative contributions at the interface vertices in the two block
inequalities and add the bridge kinetic energy. This proves(2) as a
closed quadratic-form inequality, with every harmonic and every block
charge multiplicity retained. QED.

The hypothesis matters. If several new bridges meet the same vertex,
the block charge matches their **combined** generator, not the sum of
their individual Casimirs. Equation(2) must not be applied to all
external links of a tiled lattice without that additional fusion analysis.

## 3. Repair the first rectangle's hidden floor

For two cubes every old interface degree is3. Thus a_e=4/3. In YC24's
complete even hidden decomposition, replace only its reference floors:

| Sector group | Bridge signatures | New physical reference floor |
|---|---|---:|
| E | 12,23,34,41 | 8 |
| B0 | empty, excluding all bridge constants | 32/3 |
| B4 | 1234 | 16 |
| Bd | 13,24 | 8 each |

An odd bridge costs at least(4/3)3=4; a nonconstant even bridge costs
at least(4/3)8=32/3. The total-odd bridge carrier is absent by the
whole-cube centre gauge transformation. No charged block has been
projected independently to its singlet ground.

Even bridge parity alone is insufficient for this improvement: two odd
bridges with both parents in their grounds have reference energy6 on
the unrestricted charged carrier. The new floor8 uses the physical
Gauss matching in addition to parity.

Retain YC24's complete bipartite Schur formula and its source identity
R^dagger R<=I. For interface amplitudes at most1 put

\[
 c_1(z)=4\left((32/3-z)^{-1}+(16-z)^{-1}+2(8-z)^{-1}\right),
\]
\[
 b_1(z)=4\left((32/3-z)^{-2}+(16-z)^{-2}+2(8-z)^{-2}\right),
 \quad\alpha_1(z)=8-z-c_1(z).                           \tag{3}
\]

Whenever alpha_1>0, the entire physical hidden operator is above z,
the retained return is at most1/alpha_1, and

\[
 M(z)\preceq\left[1+\frac{1+b_1(z)}{\alpha_1(z)^2}\right]I.
                                                                  \tag{4}
\]

All arguments are YC24's full-sector norm and form bounds, now using(2).
In particular no lower-floor assertion is made on an independently
charged carrier where the matching Gauss equation would be unavailable.

At z=4, alpha_1=16/15 and the theta<=2 retained reserve is
27/5-4-15/16=37/80>0. At z=9/2, alpha_1=2595/11914 and the theta<=1
reserve is11-9/2-11914/2595=9907/5190>0. The product reference has
expectation zero, so the true ground is at most zero. Full Schur inertia
therefore gives the physical gaps4 and9/2 above that true ground.

At z=0 the complete return cap is8/51 and the metric cap10705/10404.
At z=1/2 the metric cap is111440101677/107585968009<26/25. Monotonicity
in coupling and z proves the claimed whole-window norm bound.

## 4. The full transverse hidden carrier

Use the two actual28-link Hamiltonians from section3 as the parent
blocks, each shifted by its own true ground energy. They remain full
operators on their original link spaces. We do not substitute an
energy-dependent parent pencil into a tensor product or multiply two
local metric condition numbers.

The new bridges lie between y=1 and y=2, labelled by (x,z), with
x=0,1,2,3 and z=0,1. A new face corresponds to a nearest-neighbour
edge of this4×2 interface grid: six x edges and four z edges, ten in all.
Its two bridge endpoints determine its centre-flip signature uniquely.

The parent internal endpoint degree is3 at x=0,3 and4 at x=1,2.
Equation(2) therefore gives a_e=4/3 at the four end bridges, and5/4
at the four middle bridges. Let m be the subset of odd bridges. The
global physical carrier has even |m|, so it has128 parity sectors.
Each is still an infinite-dimensional harmonic and block-state space.
For the bridge-vacuum cut P, its hidden reference floors are

\[
 d_m=\begin{cases}
 10,&m=\varnothing\text{ with }P\text{ removed},\\
 \sum_{e\in m}w_e,&m\ne\varnothing,
 \end{cases}\qquad
 w_e=\begin{cases}4,&x=0,3,\\15/4,&x=1,2.\end{cases}   \tag{5}
\]

The10 keeps the first nonconstant even middle bridge. Discarding the
empty-signature hidden sector would drop genuine returning states.

Fix a new-face profile zeta_p and write t=t_y. For the complete even
hidden operator D_t and source R_t=Q(t sum zeta_p Wp)P, introduce the
scalar comparison on all128 signatures:

\[
 K_{t,z}=\operatorname{diag}(d_m-z)-t A,\quad
 A_{mn}=\#\{p:m\mathbin\triangle\partial_b p=n\},\quad
 u_m(t)=\begin{cases}t/2,&m=\partial_b p\text{ for a face }p,\\
 0,&\text{otherwise}.\end{cases}                        \tag{6}
\]

Here t>=0 bounds the magnitudes of the actual ten couplings. The source
factor1/2 follows from complete conditional Haar integration
P Wp^2 P=I/4. Distinct face signatures are different sectors. Every
interaction entry has operator norm at most t. These bounds hold on
arbitrary retained parent states; they do not assume a closed source span.

For eta_m=||psi_m|| the full hidden quadratic form satisfies

\[
 \langle\psi,(D_t-z)\psi\rangle\ge\eta^T K_{t,z}\eta. \tag{7}
\]

**J2 — full inverse and norm domination.** If K_(t,z) is positive,
then D_t-z is positive. Its inverse has the complete bounds

\[
 \boxed{0\preceq\Sigma_t(z)=R_t^\dagger(D_t-z)^{-1}R_t
       \preceq s(t,z)I,\quad s=u^T K_{t,z}^{-1}u,}      \tag{8}
\]
\[
 \boxed{I\preceq M_t(z)=I+R_t^\dagger(D_t-z)^{-2}R_t
       \preceq[1+\|K_{t,z}^{-1}u\|_2^2]I.}            \tag{9}
\]

For(8), use(7) in the variational inverse formula and bound the source
pairing componentwise by u_m||psi_m||||f||. For(9), expand the inverse
around its block-diagonal reference. Positivity of the symmetric
nonnegative-adjacency comparison implies its scaled adjacency has
spectral radius<1. The full block Neumann series is componentwise
dominated by the nonnegative scalar series, giving
||[(D_t-z)^-1 R_t f]_m||<=(K^-1u)_m||f||. Sum the squares. This proves
the original graph norm; it is not obtained by differentiating an
inequality for the returned energy.

The actual Schur pencil is F_t=h_A+h_B-z-Sigma_t and -partial_z F_t=M_t.
Both source and metric thus retain all neutral mixing and every harmonic.

## 5. Exact finite comparison, not a harmonic approximation

The comparison in(6) has two reflections, x->3-x and z->1-z. They give
44 orbits of the128 even masks. This symmetry belongs to the uniform
upper comparison; the actual ten couplings and the two parent blocks
need not have reflection symmetry.

For each required (t,z), the accompanying short exact calculator solves
K y=1 and K x=u using this orbit reduction, then lifts the answers and
checks **every original128-row equation**. All entries of y are strictly
positive and all entries of x are nonnegative. This proves positivity
on the entire128-dimensional comparison, including nonsymmetric vectors:

\[
 v^\dagger K v=\sum_m\frac{(Ky)_m}{y_m}|v_m|^2
 +\sum_{m<n}(-K_{mn})y_my_n
       |v_m/y_m-v_n/y_n|^2>0\quad(v\ne0).              \tag{10}
\]

The return is u^T x and the metric excess is sum_m x_m^2. Orbit
multiplicities are retained in both sums. Only exact rational arithmetic
is used for the following outward bounds:

| t | z | Return upper bound | Metric upper bound when needed |
|---:|---:|---:|---:|
| 1 | 0 | 77/100 | 10/7 |
| 1 | 1/2 | 109/100 | 101/50 |
| 1 | 1 | 113/50 | — |
| 1 | 11/10 | 3 | — |
| 3/4 | 5/2 | 4/5 | — |
| 1/2 | 7/2 | 6/25 | — |
| 1/2 | 4 | 1/3 | — |

The scalar inverse series has nonnegative coefficients. Its components,
source return and norm bound increase with t and z while the comparison
stays positive. The endpoint calculations therefore cover every smaller
nonnegative profile and all smaller spectral parameters. Arbitrary signs
with the same magnitude caps are also dominated, although the physical
tables use nonnegative couplings.

On P the two parents are individually gauge invariant. Their reference
has one zero vector and remaining spectrum at least delta=4 or9/2.
Subtracting z and the return bound in the table gives the respective
codimension-one reserves

\[
 37/50,\quad2/5,\quad7/10,\quad13/50,\quad1/6>0.        \tag{11}
\]

Hidden positivity and full form Schur inertia allow at most one physical
level below each trial z. Its true ground is nonpositive by the parent
product trial, hence the actual gaps are the five headline values. The
finite compact elliptic ground is unique and positive. The new graph has
maximum degree5, so YC22 gives charged gap>=3/5 separately above it.

## 6. Complete tower and zero compatibility

Let J be the product of centre flips over the four black vertices of
the checkerboard4×2 interface grid. Every joining face crosses J once.
The bridge-vacuum cut is on its plus side. Splitting the complete hidden
carrier into its two J sides gives the same exact all-order structure as
YC24, with the present larger operators:

\[
 D_t-z=\begin{pmatrix}A_z&-tC\\-tC^\dagger&B_z\end{pmatrix},
 \quad T_z=A_z^{-1/2}CB_z^{-1}C^\dagger A_z^{-1/2}\succeq0,
 \quad S_z=A_z^{-1/2}R,
\]
\[
 \boxed{\Sigma_t=t^2 S_z^\dagger(I-t^2T_z)^{-1}S_z.}     \tag{12}
\]

Here R is the source at unit amplitude for the fixed profile. The
positive comparison proves coercivity of the full hidden block; its
off-diagonal scaled norm is strictly below one, so the positive return
tower converges on these windows. This also follows from the dominated
full Neumann series. Complete operators replace neither of their infinite
carriers by one vector per parity sector.

At t_y=0 the source and its full energy return vanish exactly, the
hidden lift is zero, and M=I on the original full parent carrier. Thus

\[
 \pi_{t_y=0} F(t_x,t_y;z)=h_A(t_x)+h_B(t_x)-z,\qquad
 \pi_{t_y=0}M(t_x,t_y;z)=I.                             \tag{13}
\]

The first source derivative and higher returned coefficients remain in
the family. Applying the first-join zero-reading too leaves the four
cube operators on the retained sixteen-bridge-vacuum cut, with identity
metric. This agrees with specializing the original joint family(1)
first and then constructing that retained reading. It is the typed
construction/zero compatibility; no claim that both parameters exhaust
lambda-space or that pure reading deletes earlier source history is used.

## 7. Lattice consequence and the next concrete obstruction

For the28-link reference, the actual charged floor3/4, YC25 matched
support floor15/8, incidence48, boundary-force bounds and creator ball
all remain valid. Only the physical block floor improves. Consequently
YC25's complete relative physical argument immediately gives

\[
 \Delta_G\ge4(1-36/245)=836/245
       \quad(\theta_c\le2,\ \eta_p\le1,\ \epsilon_p\le1/1792),
\]
\[
 \Delta_G\ge(9/2)(1-192/1225)=9297/2450
       \quad(\theta_c\le1,\ \eta_p\le1,\ \epsilon_p\le1/1680).
                                                                  \tag{14}
\]

These use the same rectangular tori as YC25. At the older1/4000 external
cap the second floor is36783/8750. The additional64-link block is a
distinct result, not the reference used in(14).

For a4×4×2 block the external incidence is
2(10+10+24)=88. Even in its half-strength transverse window, where its
physical floor is4, YC25's matched per-factor floor is only
min(4,(3/5+3)/2)=9/5. Its consistently estimated incidence/floor is
440/9, exceeding the parent's128/5. At full transverse strength and
theta<=1, the present certified physical floor11/10 instead limits
the extensive floor, giving ratio80. Thus stronger local gap estimates
and a successful second join have not yet made repeated spatial blocking
contractive.

The next source to control is now explicit: several exterior bridge
generators can meet one vertex and fuse into a smaller total charge,
including a singlet. The one-bridge law(2) then retains a combined
Casimir with cross terms. A useful next boundary estimate must retain
that fusion and its connected neutral return, rather than replacing
the block by one weak charged floor per arbitrarily many boundary faces.
Physical running and the4D continuum observable construction remain
additional obligations.

## 8. Focused verification and provenance

```sh
python physics/yc26/yc26_block_comparison.py
```

This performs the new argument's exact positive-vector/source equations,
outward return and metric arithmetic, spectral reserves and zero return.
It also reproduces the inherited lattice arithmetic. There is no new
unit-test suite, predecessor regression run or certificate-hash campaign.
The128 scalar rows bound full harmonic sectors by(7)–(9); the program
does not compute or diagonalize a truncated Yang–Mills spectrum.

Source route at `759f8f5`: YC22 sections1–3 for the true-ground gauge
Casimir forms; YC24 sections1–5 for complete bipartite elimination and
metric transport; YC25 for the matched-support and lattice-join theorem.
Equations(2),(7)–(10) give the new bridge/parent matching and comparison
arguments explicitly. They are written proofs with exact finite comparison
arithmetic, not a claim of independent peer or proof-assistant review.
