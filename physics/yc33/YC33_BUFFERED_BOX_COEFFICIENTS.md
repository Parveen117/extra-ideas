# YC33 — moving observation frame and buffered local-box coefficients

10 October 2026. Research owner: Monty Dabas. Continues YC32 at `2452c03`.
Keep YC27's declared profiles, full electric carrier and existing physical
gap floors. Frozen predecessor packets are unchanged.

**Result.** YC32's reduced coefficients can now be constructed from finite
buffered boxes, with a volume-independent error budget. The box calculation
uses bounded spectral filters and a source-restricted unitary flow. It
does **not** require an assumed spectral gap or the interacting ground
state of the open box. Complete local harmonics remain in the construction.

If R is the retained coefficient range and L>=max{4,R} the buffer parameter,
the collective error is

\[
 \boxed{\varepsilon_{R,L}=e_R+\sqrt{N_R}\,\eta_L,
       \qquad N_R=\kappa(2R+1)^3,\quad\kappa=6750.}              \tag{1}
\]

Here e_R is YC32's spatial response tail and eta_L is an explicit local-box
comparison bound below. For fixed R, eta_L tends to zero faster than every
inverse power if correspondingly high locality bounds are used. Neither
term contains the total number of factors. Symmetrization and the positive
correction epsilon I preserve the reduced Fredholm reserve.

The owner's moving-diagonal interpretation also has a precise connection
identity in section1. It applies to the transported vacuum spectral split;
a general retained source cut is not automatically an energy-diagonal cut.
The scalar YM path remains an adapter, not the whole primitive lambda-space.

## 1. What the moving observation frame actually cancels

Use YC31's full finite-lattice family H(s), its actual vacuum projection
Pi_s and the unitary U_s with U_s'=iD_s U_s. The filter conventions are

\[
 D_s=\int W_\gamma(t)\tau_t^s(H'(s))dt,\qquad
 W_\gamma'=\delta_0-w_\gamma,\qquad\gamma=1/4.        \tag{2}
\]

The real even w_gamma has integral one and frequency support in
[-gamma,gamma]. Integration by parts, including the jump at zero, gives

\[
 \boxed{B_s:=H'(s)+i[H(s),D_s]
              =\int w_\gamma(t)\tau_t^s(H'(s))dt.}              \tag{3}
\]

The commutator identity is first a weak identity on Dom(H(s)); its right
side minus H'(s) is bounded at finite volume and gives the bounded
commutator extension. The vacuum-to-excited frequencies have magnitude
at least1/2 on the full tensor carrier, whereas the filter support ends
at1/4. Thus

\[
 [B_s,\Pi_s]=0,\qquad
 \Pi_s'=i[D_s,\Pi_s],\qquad
 [H'(s),\Pi_s]+[H(s),\Pi_s']=0.                    \tag{4}
\]

In transported coordinates, the derivative of U_s* H(s) U_s is the
pullback of B_s, interpreted on the common-domain/form core. Its
vacuum/excitation cross-block is zero. This is the exact sense in which
the chosen observation frame moves to keep that spectral split diagonal.
The diagonal energy derivative need not vanish. This condition neither
equates every information-curvature functional nor proves a mass scale.

For another retained cut P_s=U_s P0 U_s*, the transport equation still
holds, but [H(s),P_s]=0 is not guaranteed. The path, source and choice of
connection determine the move. Adding a generator commuting with Pi_s
can retain the same vacuum transport while changing the excited frame.
The diagram alone does not select that remaining freedom.

## 2. A finite box with a smaller source region

Work on YC31's factor graph. Its complete factors carry h_f, including
all internal interactions. A mixed face p has four-factor support Xp,
diameter one, multiplier Wp of norm<=1, and incidence at most I=48 or88.
Let ell be the admitted mixed strength, with the same caps as YC27.

Around a source factor f choose
C=B(f,L) and B=B(f,2L). On the full local Hilbert space of B construct

\[
 H_B(s)=\sum_{g\in B}h_g-s\ell\sum_{Xp\subset B}\xi_pW_p,
 \qquad J_C=-\ell\sum_{Xp\subset C}\xi_pW_p,                    \tag{5}
\]

and the buffered flow

\[
 D_{B,C}(s)=\int W_\gamma(t)\tau_t^{B,s}(J_C)dt,
 \qquad U_{B,C}'=iD_{B,C}U_{B,C},\quad U_{B,C}(0)=I.            \tag{6}
\]

H_B has the inherited unbounded electric domain and bounded potentials.
J_C is bounded. Equations(5)-(6) exist without a gap assumption for H_B.
This is a source-restricted comparison flow: J_C is generally not H_B'.
There is no claim that U_(B,C) transports the actual open-box vacuum.
Gauge symmetry is preserved, since factors and complete plaquettes, not
fragments of Wilson loops, were selected.

Let alpha_BC be conjugation by the endpoint U_BC. At the endpoint H_B(1),
use the *same* inverse-source filter k as YC32:

\[
 \mathcal R_B(A)=\int k(t)\tau_t^{B,1}(A)dt,\qquad
 \mathcal D_f^B(\phi)=\alpha_{BC}^{-1}\mathcal R_B\alpha_{BC}
                          (C_f(\phi)).                         \tag{7}
\]

The frequency cutoff is zero near zero, so this is a bounded map even if
the open box has low excited levels. It is not being identified with an
uncertified inverse of H_B minus its ground energy. All inputs in (7) are
confined to B: complete factor operators, local couplings, the parent
filter parameters and the reference product state Omega_B.

## 3. A boundary estimate for physical dynamics and its filters

Compare the full dynamics with the Hamiltonian obtained by deleting mixed
faces crossing the boundary of B. The latter factors into H_B plus its
outside Hamiltonian. There are at most I|B| deleted faces, each of norm
<=ell. If an observable A has support X inside B, put
d0=dist(X,B^c)-1 and n=|X|. A crossing face is at distance at least d0
from X. Duhamel's formula and YC31's propagation bound give, for d0>=0,

\[
 \|\tau_t(A)-\tau_t^B(A)\|
 \le\|A\|\min\{2,\,2\ell I|B|n|t|e^{-d_0+v|t|}\},
 \qquad v=8eI\ell .                                \tag{8}
\]

Only bounded crossing interactions enter this comparison. The unbounded
on-site terms remain on both sides and preserve support in the interaction
picture. The estimate is uniform along the path and applies to either
sign of the auxiliary time t.

For an integrable weight w with first absolute moment M1(w), define

\[
 b_w(d,n,N)=\min\{2L_w,\,
    2\ell I N n M_1(w)e^{-d/2}+2T_w(d/(2v))\},\quad d\ge1,     \tag{9}
\]

where L_w=integral |w| and T_w(t)=integral_(|u|>t)|w(u)|. Set b_w=2L_w
at d=0. Splitting the time integral at d/(2v) proves

\[
 \left\|\int w(t)(\tau_t(A)-\tau_t^B(A))dt\right\|
       \le\|A\|b_w(d_0,n,|B|).                    \tag{10}
\]

At ell=0 all these differences vanish; define the corresponding bounds
as zero instead of dividing by v=0. The explicit |B| factor is allowed:
the box size depends on the requested accuracy, not on total lattice volume.

## 4. Compare the true transport with the buffered flow

Use YC31's flow-locality bound with exponent32. Write
l(r)=min{2,C32(r+1)^(-28)} for either orientation of any path segment of
length<=1. Let epsilon_W(r,4) be its single-source filtered localization
bound from YC31(5), with epsilon_W(0,4)=2L_W. For integer d>=1 define

\[
 r_d=\lfloor(d-1)/3\rfloor,\qquad
 c(d,n)=\min\{2L_W,\,
              2\epsilon_W(r_d,4)+2L_W n\,l(r_d)\}.              \tag{11}
\]

Consider one full generator term d_p(s) and a globally transported
observable initially supported on X. Localize each in radius r_d, with
d=dist(Xp,X). Their localized supports are disjoint. The two localization
errors therefore bound their commutator by ell ||A|| c(d,|X|).

First replace D_s by D_C^global(s), the sum of its original terms with
Xp subset C, still using the *full* physical dynamics. For X subset B(f,a),
a<L, every omitted face has min_(g in Xp) dist(f,g)=q>=L. At most
I kappa(q+2)^3 faces have a given q. A unitary-flow Duhamel comparison
thus has the uniform local error

\[
 C_L(a,n)=\ell I\kappa\sum_{q\ge L}(q+2)^3 c(q-a,n).            \tag{12}
\]

This comparison uses the locality of the full flow on the test observable;
it never estimates the norm of the extensive omitted generator by its
number of terms.

Next replace the physical dynamics inside each retained generator term
by H_B dynamics. Those terms are at boundary distance d0>=L. There are at
most I kappa(L+1)^3 such faces. Equation(10) gives the full bounded-generator
error

\[
 z_L=\ell I\kappa(L+1)^3
        b_W(L,4,\kappa(2L+1)^3).                   \tag{13}
\]

The resulting unitary error is <=z_L, and conjugation costs at most2z_L.
Combining both steps gives, for either alpha or its inverse,

\[
 \|\alpha^{\pm1}(A)-\alpha_{BC}^{\pm1}(A)\|
 \le\|A\|a_L(a,n),\qquad
 a_L(a,n)=\min\{2,C_L(a,n)+2z_L\}.                 \tag{14}
\]

Only the global transport uses the established global spectral gap.
The comparison flow needs bounded filtered generators and locality.
This separation is why an extra open-box gap assumption is unnecessary.
The classical comparison mechanism is Duhamel plus Lieb–Robinson;
the unbounded-site version is covered by section3.2 and Theorem3.4 of
[Nachtergaele–Sims–Young, arXiv:1810.02428v2](https://arxiv.org/abs/1810.02428v2).
Equations(8)-(14) give the particular conservative bounds used here.

## 5. Assemble the full three-map box response

Let a=floor(L/4), n_a=kappa(a+1)^3 and n_(2a)=kappa(2a+1)^3. Recall the
exact column operator D_f=alpha^-1 R alpha(C_f) from YC32. Then

\[
 \boxed{\|\mathcal D_f(\phi)-\mathcal D_f^B(\phi)\|
             \le\eta_L\|\phi\|,}                               \tag{15}
\]

where one explicit bound is

\[
 \begin{split}
 \eta_L=\min\{2L_k,\;&4L_k l(a)+2\epsilon_k(a,n_a)\\
                 &+b_k(2L-a,n_a,\kappa(2L+1)^3)\\
                 &+L_k[a_L(0,1)+a_L(2a,n_{2a})]\}.             \tag{16}
 \end{split}
\]

For clarity, the composition proof retains the following two local
intermediates. First set A1=E_(B(f,a)) alpha(C_f), with error l(a)||phi||.
Then set A2=E_(B(f,2a)) R(A1), with total error
[L_k l(a)+epsilon_k(a,n_a)]||phi|| relative to R alpha(C_f), and
||A2||<=L_k||phi||. Apply (14) to A2 for the outer inverse transport.
The physical filter difference on A1 is bounded by (10). The inner
transport difference is (14) on C_f. Replacing the two nonlocal
intermediates costs twice their respective errors. Their sum is (16).

With the displayed flow exponent and sufficiently high finite moments
of W and k, eta_L=O((L+1)^(-21)); larger exponents give higher powers.
The constants are explicit series and filter integrals, not optimized
numerical radii. The construction and its proof use complete local
operator spaces, not a finite sample of their matrix elements.

## 6. Collective assembly, centering and positivity

The box state used for coefficients is Omega_B, the known product of
factor grounds. It need not be the interacting box ground. Center each
box column in this reference:

\[
 D_f^{B,c}(\phi)=D_f^B(\phi)
       -\langle\Omega_B,D_f^B(\phi)\Omega_B\rangle I.            \tag{17}
\]

For d(g,f)<=R, define the raw coefficient block by pairing D_f^(B,c)Omega_B
with the one-factor excitation at g; set it to zero for d(g,f)>R. This
requires only the box around f since R<=L. Call the assembled operator
T_box^raw on YC32's complete retained space E.

Different columns use different boxes, so this raw operator need not be
self-adjoint. Before correcting it, control its error collectively.
Apply the vector support projection L_(B(f,R)) to each centered box
column and compare with the exact YC32 column. Centering is Q0 on the
vector, and Q0 commutes with that support projection. Since the exact
column already has zero vacuum component, (15) bounds the projected
vector error by eta_L ||phi||, without a factor two. Errors for disjoint
radius-R balls are orthogonal in the product reference. YC32's overlap
argument gives

\[
 \|T_{\rm box}^{\rm raw}-T_R\|
       \le\sqrt{N_R}\eta_L .                       \tag{18}
\]

The same bound holds for the corresponding response synthesis error
before compression. No local-channel dimension appears. Symmetrize and
use (1):

\[
 B_{R,L}=\tfrac12(T_{\rm box}^{\rm raw}
                         +(T_{\rm box}^{\rm raw})^*),\quad
 \|B_{R,L}-T\|\le\varepsilon_{R,L},\quad
 T\preceq B_{R,L}+\varepsilon_{R,L}I
                \preceq T+2\varepsilon_{R,L}I.                  \tag{19}
\]

Here T=V* A^-1 V is the exact inverse compression including the full
hidden carrier. Thus S_box=[I+tau(B_(R,L)+epsilon I)]^-1 satisfies

\[
 \|S-S_{\rm box}\|\le2\tau\varepsilon_{R,L},\qquad
 \frac{I}{1+\tau(M+2\varepsilon_{R,L})}
            \preceq S_{\rm box}\preceq S\preceq I,\quad M=1/m. \tag{20}
\]

Apply YC32's positive geometric polynomial to this box-built matrix.
Its degree-p version has retained spatial range pR and error

\[
 \boxed{\|S-S_{R,L,p}^{\rm box}\|
       \le2\tau\varepsilon_{R,L}+q^{p+1},\quad
 q=1-\frac1{1+\tau(M+2\varepsilon_{R,L})}.}          \tag{21}
\]

Its coefficients depend on boxes of radius2L around their source factors.
A degree-p row therefore needs data within radius pR+2L, not the whole
lattice. This is a local-data construction of the selected reduced
Fredholm approximation, not merely a formal band truncation of an unknown
global matrix.

For tau=4, first choose R with e_R<=1/1600. Then choose L>=max{4,R} such
that sqrt(N_R)eta_L<=1/1600. Both choices exist independently of total
volume. YC32's degree8 bound <3/200 and its reserve table now apply to
coefficients obtained from these finite boxes. No practical numerical
values of R,L are certified here.

An eventual numerical box solver has its own error. If every returned
block is accurate in *operator norm* to nu, with at most N_R blocks per
row and column, the block Schur test adds at most N_R nu to epsilon.
Entrywise accuracy on a finite harmonic sample does not supply this
operator-norm bound on the complete local channels. A certificate for
coefficient blocks alone also does not certify separately computed
response vectors used in a numerical lift.

## 7. Physical source identification and what the box result does not replace

To identify the box response as a vector in the actual physical Hilbert
space, apply the original exact U to its centered reference synthesis.
Then vacuum leakage is exactly zero and its collective error is bounded
by (1). YC32's constraint-preserving lifts and source/metric estimates
apply with this enlarged budget. The numerical coefficient construction
itself requires only (5)-(7),(17); the exact U is retained for the theorem
identifying the physical carrier.

Replacing that final U by unrelated separate local unitaries is not part
of this claim. Nor have energy-domain error bounds for a numerically
approximated physical lift been supplied. The parent volume-independent
lattice gap remains YC27's established input in lattice electric units.

The retained carrier is a direct sum of one-factor neutral excitations.
It is not yet a tensor product of coarse YM link spaces and is not closed
under multiplying two independent source creations. Consequently the
finite-range operator (21) alone is not a full spatial RG step. This is
the next structural issue, now separated from the local-data problem.

## 8. Zero specialization and direct checks

At ell=0, both true and buffered transports are identity, all boundary
comparison errors vanish, and the filtered response of a single-factor
source stays on that factor. Hence eta_L=e_R=0 and the assembled inverse
compression is exactly the direct sum of inverse neutral factor energies.
This is YC32's zero reading. The frequency filters, local domains,
pairings and retained first source jet remain; they are not erased by
the zero value of the interaction. Finite geometric polynomials specialize
to the corresponding reference polynomial with its stated residual.

```sh
python physics/yc33/yc33_box_assembly.py
```

The direct exact check verifies the observer connection cancellation on
a moving two-level model with changing diagonal energies, a deliberately
non-Hermitian raw coefficient assembly, symmetrization/positive correction,
and the two-part rational error budget. All-volume boundary and synthesis
bounds have written proofs above. No harmonic truncation, predecessor
rerun or generic test suite was used. A practical local-box computation
still requires numerical filter/energy error control.

**Next task:** extend the retained construction to simultaneous sources
on several factors, retaining their products, connected returns and
physical metric. Determine whether the resulting effective operator
admits a controlled coarse interaction family. Neither a shrinking lattice
spacing nor its physical energy conversion follows from the present
fixed-spacing box localization.
