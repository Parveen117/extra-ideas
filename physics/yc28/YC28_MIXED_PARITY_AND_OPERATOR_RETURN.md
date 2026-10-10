# YC28 — the mixed-interface cut and its full operator return

10 October 2026. Research owner: Monty Dabas. Continues YC27 at `b67d93e`.
Unit-S³ electric normalization, normalized Haar, Wp=Tr(Up)/2. Keep YC27's
geometry and full factor carriers. The algebraic and charged-return bounds
below allow arbitrary finite nonnegative internal Wilson couplings in both
partitions. The numerical susceptibility improvements, analytic vacuum disk
and lattice gap still use YC27's stated windows. No predecessor packet
is changed.

**Result.** The remaining mixed interaction has an exact global grading:
a product of old-block vertex-centre transformations reverses every mixed
face and leaves both interacting reference partitions unchanged. Eliminating
its odd side gives an exact operator pencil **quadratic in lambda**, with
the original graph metric. This identity holds at every real mixed strength
on each admitted finite lattice for spectral z<4; it is not a small-coupling
Taylor replacement.

On the fully factor-neutral target, distinct first face sources are
orthogonal under the complete reference resolvent. For each face its
normalized source has the exact **operator** mean and variance

\[
 \boxed{U_p^\dagger H_0U_p=H_P+12,
 \qquad U_p^\dagger H_0^2U_p-(H_P+12)^2=4K_p,}         \tag{1}
\]

where Kp is the electric energy of its four selected links. The returned
operator is therefore controlled on excited neutral states as well as on
the vacuum. Spectator energies remain in HP.

At fourth order, four distinct mixed faces return only as the four sides
of an elementary cube. Each mixed face belongs to exactly two such tubes.
Their two ends can be old loops **or complementary loops**. The full
24-order tube return has norm<=5/64 at z=0, at arbitrary finite internal
couplings. Repeated-face channels also remain; section5 gives their bounds.

The charged inverse floor4, the operator identities and the tube norm
bound are independent of the internal coupling sizes. They do not give
a lower bound on a neutral block's excitation gap at those arbitrary
couplings.

These results resolve the next source and its locality pattern. The physical
lattice gap windows remain YC27's592/175 and11188/3325 at mixed caps1/3200
and1/5700. No wider gap interval, contracting spatial iteration or4D
continuum theorem is claimed.

## 1. Construct the family and the grading before the return

Use YC27's old factors O and complementary cube/square/link factors D on
the4×2×2 or4×4×2 tilings. Each shared vertex belongs to exactly one factor
of each kind. On the full link Hilbert space with its Haar pairing, let

\[
 H_0(u)=\sum_{f\in O\cup D}h_f(u),\qquad
 V=\sum_{p\ {\rm mixed}}\xi_pW_p,
 \qquad H_\lambda(u)=H_0(u)-\lambda V.                 \tag{2}
\]

Here hf is above its actual positive ground. The carried scalar reference
energy can be restored. The common domain is the inherited finite-product
Laplacian domain; V is bounded. All statements below about physical operators
use the globally gauge-invariant carrier. Factor singlets are imposed only
at the explicitly named projection P. This is the same operator-space
adapter as YC27, not a new identification of primitive lambda-space with
the scalar coordinate lambda.

For a lattice vertex v=(x,y,z), define the binary weight

\[
 \omega(v)=xy+xz+yz\pmod2.                            \tag{3}
\]

Every side length is even, so this is well-defined periodically. Let Zv^O
apply the centre element -I at v to its **old factor's links only**, and set

\[
 \boxed{J=\prod_v(Z_v^O)^{\omega(v)},\qquad J^2=I.}    \tag{4}
\]

This is a unitary self-adjoint cut. It commutes with the full global gauge
action and every hf. It is generally not a global gauge transformation:
the complementary links at the same vertex are not transformed.

A mixed face has two opposite old edges. Their four endpoints are precisely
its four corners. The second finite difference of(3) in any two distinct
coordinate directions is1. Hence

\[
 JW_pJ=-W_p\quad(p\ {\rm mixed}),\qquad
 JH_0J=H_0,\qquad \boxed{JH_\lambda J=H_{-\lambda}.}   \tag{5}
\]

Old internal faces are unchanged by old-factor gauge invariance, and
complementary internal faces have no transformed links. This proves(5)
for arbitrary allowed nonuniform internal profiles and every finite volume
in YC27, including the shortest periodic sides4.

The actual ground energy is even in real lambda by unitary equivalence.
On YC27's analytic vacuum branch it has only even Taylor coefficients.
The full creator satisfies c(-lambda)=Jc(lambda), since J fixes the
reference ground and respects its factor cuts. Its odd source jets do
not vanish; they live on J's odd side. No quadratic estimate for the
entire ordered excitation spectrum follows merely from this symmetry.

## 2. An exact square-coupling pencil, with its metric

Put Pplus=(I+J)/2, Pminus=(I-J)/2 and Hplus/minus=Pplus/minus H0 Pplus/minus.
Write T=Pminus V Pplus. A J-odd physical state cannot be neutral in every
factor. YC27's matched-charge proof gives Hminus>=4. Only its charged part
is needed: the paired coefficient is>=1/3; odd charges form an even-degree
subgraph of the simple bipartite incidence graph and cost>=4; a nonzero
even charge costs>=16/3 by the two matched single-vertex bounds. YC22's
true-ground Casimir inequalities hold at every finite internal coupling.
No internal neutral excitation floor is used here. Consequently for
every real z<4 the complete odd inverse exists, independent of lambda.
Exact elimination gives

\[
 \boxed{F_+(\lambda,z)=H_+-z-\lambda^2T^\dagger(H_--z)^{-1}T,}
                                                               \tag{6}
\]
\[
 \boxed{M_+(\lambda,z)=I+\lambda^2T^\dagger(H_--z)^{-2}T
                       =-\partial_zF_+.}             \tag{7}
\]

The lifted state is (f,lambda(Hminus-z)^-1 T f);(7) is its actual norm
metric. These formulas keep the entire J-even carrier, including charged
factor states. They are exact for arbitrary real lambda at this finite
volume, not assertions that all neutral interactions have been integrated
out or that the full spectrum is known. In particular T's global norm can
grow with volume; the identity alone is not a volume-uniform stability
estimate. The true ground has energy<=0 by the reference-vacuum trial,
so its odd elimination is always within z<4.

Now let P=product_f PG,f average the gauge group of each factor separately,
and let Q=I-P on the physical carrier. P is contained in Pplus. Set HP=PH0P
and, on Q, Gz=[Q(H0-z)Q]^-1. Again QH0Q>=4. For one mixed face put

\[
 U_p=2W_pP,\qquad U_p^\dagger U_p=I_P,\qquad
 R_p(z)=PW_pG_zW_pP.                                 \tag{8}
\]

Independent factor gauge averaging makes each selected open-edge marginal
Haar, proving PWp²P=I_P/4. WpP has fundamental endpoint charge at the four
corners in both matched factors. Distinct elementary faces have different
corner sets on these tori. H0 preserves all factor vertex-charge isotypic
spaces, including their multiplicities. Thus for p!=q,

\[
 PW_pG_z^nW_qP=0\quad(n=1,2).                         \tag{9}
\]

This gives the exact neutral **compression** of(6),(7):

\[
 PF_+P=H_P-z-\lambda^2\sum_p\xi_p^2R_p(z),\qquad
 PM_+P=I_P+\lambda^2\sum_p\xi_p^2J_p(z),
 \quad J_p(z)=PW_pG_z^2W_pP.                          \tag{10}
\]

The remaining coupling to Pplus-P starts at lambda²; eliminating that
carrier creates further returns beginning at fourth order. Equation(10)
is not the fully reduced operator on P. Even its same-face second-order
kernel can couple several neutral factors; diagonal in face labels does
not mean scalar or independent block energies.

## 3. The complete four-factor kernel and operator variance

For a selected edge e in factor f, let u_alpha be its four quaternion
coordinates. On that factor's neutral carrier define

\[
 A_{f,e}(t)=\sum_{\alpha=0}^3P_{G,f}u_\alpha
                 e^{-t h_f}u_\alpha P_{G,f}.           \tag{11}
\]

Endpoint SU(2)×SU(2) covariance makes the coordinate-index compression
scalar: PG u_alpha exp(-t hf) u_beta PG=delta_alpha,beta A_(f,e)(t)/4.
This follows by gauge averaging the irreducible coordinate representation;
the scalar here is an operator on the full neutral multiplicity space.

The tensor contraction of Wp therefore gives the exact heat representation

\[
 R_p(z)=\frac14\int_0^\infty e^{tz}
 \left[\bigotimes_{(f,e)\in p}A_{f,e}(t)\right]
 \otimes e^{-tH_{\rm spectators}}\,dt.                \tag{12}
\]

There are four active factors. Each endpoint-charge carrier has its full
positive Casimir floor; the four floors sum to at least4 by YC27's paired
degree argument. Thus(12) converges in operator norm for z<4. The spectator
heat factor is retained exactly and cannot be replaced by I on excited
spectators.

Let Ke=-Delta_e. On the smooth neutral core, the identities sum u_alpha²=1,
sum u_alpha grad u_alpha=0 and the spherical coordinate metric give

\[
 A_{f,e}(0)=I,\qquad -A'_{f,e}(0)=h_{f,G}+3,
\]
\[
 A''_{f,e}(0)=(h_{f,G}+3)^2+4K_{e,G}.                 \tag{13}
\]

For the second derivative, use
hf(u_alpha f)=u_alpha(hf+3)f-2 grad u_alpha dot grad_e f.
Summing its squared norms kills the cross terms and adds4||grad_e f||².
This proof works for the actual gauge-invariant potential and shifted
ground energy; no free-state replacement is made.

Set Kp=sum_(four selected e) P Ke P. Differentiating the tensor kernel on
the common smooth core proves(1), with forms extended by closure. It
retains the entire neutral excitation operator HP, not just the constant12
seen on Omega. The variance4Kp is positive and local in its selected links.

This supplies an exact residual bound for the full return. Put
Az=HP+12-z and Lp=(QH0Q)Up-Up(HP+12). Then

\[
 U_p^\dagger L_p=0,\qquad L_p^\dagger L_p=4K_p,
\]
\[
 \boxed{R_p(z)=\frac14A_z^{-1}
 +\frac14A_z^{-1}L_p^\dagger G_zL_pA_z^{-1}.}          \tag{14}
\]

Indeed Gz Up=Up Az^-1-Gz Lp Az^-1; substitution and Up†Lp=0 prove(14).
Since Gz<=1/(4-z), it follows as an operator/form inequality that

\[
 \boxed{0\preceq R_p(z)-\tfrac14A_z^{-1}
 \preceq\frac1{4-z}A_z^{-1}K_pA_z^{-1}.}              \tag{15}
\]

All displayed products define bounded forms. For example Kp is bounded
above by the sum of the four active hf plus a finite scalar potential
bound, so Kp^(1/2) Az^-1 is bounded. The scalar is local to those factors;
no total-volume potential norm enters(15).

For the metric use the same exact lift, not the derivative of an inequality.
For f in P define a=||Az^-1 f|| and
b=2||(Kp)^(1/2)Az^-1 f||/(4-z). The two triangle inequalities give

\[
 \boxed{\tfrac14\max(0,a-b)^2
 \le\langle f,J_p(z)f\rangle\le\tfrac14(a+b)^2.}      \tag{16}
\]

The simpler bound Jp(z)<=I/[4(4-z)²] also holds. No positivity of
Jp(z)-Az^-2/4 is inferred from the positive energy remainder in(15).

On the product reference ground, HP Omega=0 and <Omega,Kp Omega> is
exactly YC27's four-edge residual variance vp. Thus the vacuum face
susceptibility obeys the following sharper bounds in YC27's internal
windows:

\[
 \frac1{48}\le\chi_p:=\langle\Omega,R_p(0)\Omega\rangle
 \le\frac1{48}+\frac{v_p}{576}
 \le\begin{cases}179/6144,&28\hbox{-link tiling},\\
                  2197/73728,&64\hbox{-link tiling}.
 \end{cases}                                         \tag{17}
\]

At arbitrary finite internal couplings the uniform weaker interval
1/48<=chi_p<=1/16 still follows from(14) and ||G0||<=1/4. No uniform
vacuum analytic radius at those arbitrary couplings is asserted.

In YC27's windows, E(lambda)=-lambda² sum_p xi_p² chi_p+O(lambda⁴). The full
analytic family is already supplied by YC27. A loose but volume-uniform
per-face remainder can also be stated: with its radius R=1/3200 or1/5700,
x=|lambda|/R<1 and M1=sum_p xi_p,

\[
 \left|E(\lambda)+\lambda^2\sum_p\xi_p^2\chi_p\right|
 \le\frac{M_1R}{14}\frac{x^4}{1-x^2}.                 \tag{18}
\]

To prove this, use the creator root ball r=1/32. Only creators supported
inside a face can contribute to <Wp Omega,exp(C)Omega>. Their norm sum is
at most4r, so |E(zeta)|<=|zeta| M1(exp(4r)-1)/2<=|zeta| M1/14 on the
complex radius R. Here exp(1/8)<8/7. Cauchy's estimate and evenness sum
the omitted even powers to(18). This is an error bound, not a numerically
sharp prediction of the interacting vacuum energy or an enlarged gap.

## 4. The complete returning face patterns through fourth order

Relative factor-centre parity at each physical vertex is preserved by H0,
P, Q and Gz. A mixed face toggles its four corner parities. Consequently a
word returning to P must have zero XOR of its corner sets. Repetitions
are counted modulo2. The grading(4) additionally kills every odd-length
neutral return, at all orders.

Through fourth order the allowed multisets are exactly:

| Order | Possible face multiset |
|---:|---|
| 1 | none |
| 2 | p,p |
| 3 | none |
| 4 | p,p,p,p; p,p,q,q; four distinct sides of one elementary cube |

Here "possible" means not excluded by this selection rule, not that every
matrix element is nonzero. For four distinct faces, the cube must have
exactly four mixed faces and two pure opposite end faces.

**All-volume proof of the new four-face classification.** Take one of
four distinct faces whose corner XOR vanishes. Its four vertices must
be covered by the other three faces, so some pair shares at least two
vertices. Distinct ordinary elementary faces share at most two; such
a pair shares an edge. The other pair has the same four-vertex XOR and
must also share an edge.

If the first pair is coplanar, its outer vertices form a2×1 rectangle.
There is a unique two-face filling of those corners, except for the
alternative two-face path around a periodic side of length4. In the latter
case that coordinate has block side2, so the two adjacent faces have
opposite boundary-crossing status there; with their other coordinate
fixed they cannot both be mixed. Thus this exceptional relation is absent.

If the first pair is perpendicular, its XOR consists of two parallel unit
edges at two perpendicular offsets from their shared edge. The only other
two-face filling uses the opposite parallel edge of that unit cube.
These are precisely the four side faces of the cube. This proves the
classification, including the minimal periodic volumes. Removing even
face multiplicities proves the other rows of the table.

Write bi=1 if the unit cube crosses an old boundary in coordinate i.
When b has one entry1 the pure ends are in two old blocks. When it has
two entries1 they are in two complementary blocks. With zero or three
entries1 there are no mixed faces. Every mixed face borders two unit
cubes, both of which have one or two boundary crossings, so it belongs
to **exactly two returning tubes**, independently of block size and volume.

## 5. Full fourth-order bounds and the old/complementary exchange

All following words use Gz on the complete physical Q carrier. Define
the tube coefficient as the sum over its24 distinct orders,

\[
 \mathcal K_t(z)=\sum_{\pi\in S_4}
 PW_{p_{\pi(4)}}G_zW_{p_{\pi(3)}}G_z
                 W_{p_{\pi(2)}}G_zW_{p_{\pi(1)}}P.     \tag{19}
\]

This is the coefficient in the neutral Schur **return**; that return is
subtracted from the effective Hamiltonian. Reverse orders are adjoints,
but no positivity of the separate tube coefficient is assumed. The full
neutral Schur family exists near zero at each finite volume by the bounded
perturbation inverse. Its raw global convergence radius is not asserted
uniform in volume; equations(6)–(7) and YC27's local analytic construction
remain the distinct complete-family statements.

An odd corner charge costs at least1 in the paired-degree sum of YC27.
After one and three distinct tube insertions four corners are odd. After
two adjacent faces four remain odd; after two opposite faces eight do.
Endpoint source norms are1/2 by(8). There are16 adjacent-pair orders and
eight opposite-pair orders. Thus, for0<=z<4,

\[
 \|\mathcal K_t(z)\|\le
 \frac4{(4-z)^3}+\frac2{(4-z)^2(8-z)},\qquad
 \boxed{\|\mathcal K_t(0)\|\le\frac5{64}.}            \tag{20}
\]

For a repeated pair p,q let r be their number of common vertices, which
is0,1 or2. A hidden even-parity but charged state costs at least16/3 by
YC27's separate single-vertex bound on its matched factors. A different
initial pair has8-2r odd vertices. When r=0 the two equal-initial-pair
orders vanish: the later face cannot remove any residual charge at the
first face's corners, while Q forbids an already completely neutral state.
The same source norms and all six pair orders give:

| Fourth-order multiset | Full return norm upper bound at z=0 |
|---|---:|
| p^4 | 3/1024 |
| p²q², r=0 | 1/128 |
| p²q², r=1 | 25/1536 |
| p²q², r=2 | 11/512 |
| Four distinct tube faces | 5/64 |

No ground-state or harmonic truncation enters these bounds. Repeated
p²q² channels on overlapping factor supports can remain even if their
corners do not overlap; full neutral dynamics inside a shared factor
must not be deleted using first-source charge orthogonality.

The static tube has an exact exchange identity. Let WL,WR be the two
pure end-face Wilson traces. In the old-ended case the four transverse
links belong to four different complementary factors; in the complementary-
ended case they belong to four different old factors. Independent factor
gauge averaging integrates each selected open link with Haar measure.
The quaternion second moment I4/4 and the four-link contraction give

\[
 \boxed{P\prod_{p\in t}W_pP=\frac1{64}W_LW_RP.}       \tag{21}
\]

This is the same character-gluing normalization as YC17, now valid with
the roles of the two correlated partitions exchanged. It is not a kinetic
coefficient. At the completely free internal point, the exact kinetic
source does recover YC17's local twelve-link calculation in either case:

\[
 \mathcal K_t(0)\Omega_0
 =\frac1{64}\left(\frac{16}{12\cdot18\cdot24}
                  +\frac8{12\cdot24\cdot24}\right)W_LW_R\Omega_0
 =\frac{11}{165888}W_LW_R\Omega_0.                    \tag{22}
\]

Completed transverse links have no later insertion on their factors;
their nonconstant charge cannot survive P. The free inverse preserves
that decomposition, giving the displayed denominators. Interacting
factors instead require the full inverses in(19);(22) is not extended
unchanged to them.

## 6. Connected return, spectators and what is still expensive

For the union X of a word's factor supports, exact spectral resolution
of the spectator reference Hamiltonian gives

\[
 \mathcal K_{\rm global}(z)=
 \int\mathcal K_X(z-E)\otimes d\Pi_{\rm spectators}(E).\tag{23}
\]

Differentiating the actual three-resolvent words, with the same floors
as(20), gives the tube's relative spectator bound

\[
 \|(\mathcal K_{t,\rm global}(0)-\mathcal K_{t,X}(0)\otimes I)f\|
 \le\frac{29}{512}\|H_{\rm spectators}f\|.            \tag{24}
\]

The constant is12/4⁴+4/(4³8)+2/(4²8²). No differentiation of an operator
inequality is used to infer an operator sign. A scalar norm bound does
not permit omitting the spectator dependence in a new lattice Hamiltonian.

For the actual ground energy, disconnected factor clusters cancel exactly.
If p and q have disjoint factor supports and all other mixed couplings are
turned off, H splits into independent factors, so its mixed energy derivative
at p²q² vanishes. The same is true in the commuting creation logarithm.
This extends YC18's connected argument to the present correlated partition.
Raw Schur words can still be nonzero before their energy-dependent
cancellation; they are not induced distant forces.

The four-distinct-face mechanism now has a fixed incidence of two tubes
per face, with selected norm cost at most(5/32)lambda⁴ for uniform mixed
magnitude lambda. This is an incidence-summed coefficient estimate, not
yet a local interaction norm after eliminating spectators. The repeated
channels with overlapping factors still grow with block boundary incidence.
They are the next concrete target for a contracting boundary estimate,
together with transport of the full excitation metric. Neither parity nor
the static loop-exchange identity supplies that contraction by itself.

## 7. Zero recovery and focused verification

At lambda=0, Fplus=Hplus-z and Mplus=I exactly. The first return starts
at lambda²; the additional return through the remaining even charged
carrier starts at fourth order, while
the actual first source derivative remains nonzero on the odd side.
The source kernel, variance and graph metric all come from the full
family(2) before this evaluation. Setting complementary couplings u=0
recovers the free-bridge reference under YC27's tensor regrouping;
setting every internal coupling to zero recovers(22). Source history
and the original Haar pairing are retained in these typed readings.

```sh
python physics/yc28/yc28_mixed_return.py
```

The calculator checks the new grading on every face of the two smallest
admitted tori and enumerates their pair-XOR collisions. All returning
four-face sets are precisely their elementary tubes:96 for8×4×4 and176
for8×8×4, with two tubes per mixed face. It also checks the new rational
susceptibility, repeated-word and tube bounds. The all-volume classification,
operator identities and domain/form arguments are written above. There is
no new unit-test suite, predecessor rerun or finite-spectrum approximation.

Source route at `b67d93e`: YC27 for the correlated partition, matched-degree
floors, source variance and full creator disk; YC22 for charge-invariant
functional calculus; YC17 for the Haar tube contraction and free kinetic
normalization; YC18 for spectator resolution and connected cancellation.
The new grading, full operator variance and four-face classification on
the complementary partition are proved here. No independent peer or formal
proof-assistant verification is claimed.
