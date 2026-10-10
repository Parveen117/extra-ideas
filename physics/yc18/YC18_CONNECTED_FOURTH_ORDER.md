# YC18 — complete fourth-order bridge selection and a local ground-return expansion

9 October 2026. Continues YC17 at `16bb4a0`. Frozen packets are unchanged.

**Result.** All bridge-vacuum return words through fourth order are now
classified, including four-bridge plaquettes and repetitions. Third order
vanishes. Fourth order contains only p^4, p^2 q^2 and YC17's four-face tubes.
Full operator bounds are given for each class. Spectator energy is retained
as an exact shifted-energy argument, with an explicit relative error bound.

The next local object is the **connected creation exponent of the actual
ground**, not the raw global bridge Schur operator. Its complete coefficients
through order four have a local recursion, and the full omitted tail has a
volume-independent bound. The tube's connected creation coefficient at free
internal coupling is **11/3981312**. Disconnected four-face histories cancel
in this exponent, even though their raw Schur return need not vanish.

This does not close the excitation operator under spatial blocking or enlarge
the gap window. YC15's actual gap>9/40 through eta<=1/4320, theta_c<=1 remains
the numerical baseline. Locality of a vacuum expansion is a different target
from an effective Hamiltonian valid on all excited states.

## 1. The declared carrier and all low-order return words

Use YC17's actual cube/bridge tiling on even spatial tori with sides >=4,
normalized SU(2) Haar measure and unit-S3 electric energies l(l+2). Let
H_0=sum h_c+sum_b(-Delta_b), with actual cube h_c>=0 and theta_c in[0,1].
P_b fixes every bridge to its constant function and retains every cube state;
Q_b=I-P_b. G_z is YC17's full Q_b inverse, defined for real z<3.
The new bounds below use z=0. This is the original reference, with no terms
removed from the original Hamiltonian by forgetting YC16 compensation.

For an external face p let B_p be its set of bridge edges and m_p=|B_p|,
which is 2 or 4. Each W_p is odd under the centre flip of every bridge in
B_p. The reference, its resolvents and its cuts preserve these parities.
Thus a necessary return condition is XOR_p B_p=empty, with multiplicity.

**Selection theorem.** A returning word of length <=4 has the following form:

| Length | Possible multisets |
|---:|---|
| 1 | none |
| 2 | p,p |
| 3 | none |
| 4 | p,p,p,p; p,p,q,q with p!=q; four distinct sides of one two-bridge tube |

Proof. Distinct ordinary elementary plaquettes on these tori share at most
one edge. If a set of at most four distinct external faces contained a
four-bridge face, each of its four bridges would need another face to cancel
its parity. At least four *other* distinct faces are required, impossible.
If all faces have two bridges, YC17's disjoint four-cycle graph shows that
the only minimal nonempty returning subset is a complete tube. Distinct
face bridge signatures are consequently different. Remove even face
multiplicities from a word: the possible odd supports now prove every entry
of the table. In particular the conclusion for order three includes mixed
two/four-bridge words, extending YC17's restricted selection statement.

The proof is for all the stated volumes. The certificate also enumerates
all pair-XOR collisions on three different tori, including periodic side4.
No claim about arbitrary higher-order odd words follows from this table.

## 2. Bounds for the repeated channels, with the full hidden spectrum

Every endpoint satisfies ||W_p P_b||=1/2. After one or three insertions,
the odd bridge signature gives floor 3m_p for its corresponding face.
After two equal faces the signature is empty, but the intermediate Q_b
requires at least one nonconstant **even** bridge harmonic. Its floor is
8, not zero and not the odd fundamental floor3. The cube and all spectator
energies remain nonnegative and are retained in the inverse.

The single p^4 word therefore obeys

\[
 \boxed{\|K_{p^4}\|\le\frac1{288m_p^2}.}                 \tag{1}
\]

For p!=q put r=|B_p intersect B_q| in{0,1}, and s=m_p+m_q-2r.
There are six distinct orders of p,p,q,q. The two with an equal initial
pair have middle floor8; the other four have middle floor3s. If r=0,
the equal-initial-pair words are **zero**: every bridge of the completed
first face must already be constant to survive the final projection,
whereas the middle Q_b forbids all bridges being constant. No later face
touches those bridges, so this projection can be commuted backwards.
Hence the sum of the six ordered words has norm bounded by

\[
 \boxed{B(m,n,r)=\frac{\mathbf1_{r=1}}{144mn}
            +\frac{(m+n)^2}{108(m+n-2r)m^2n^2}.}         \tag{2}
\]

| m,n | No common bridge | One common bridge |
|---|---:|---:|
| 2,2 | 1/432 | 11/1728 |
| 2,4 | 1/1152 | 5/2304 |
| 4,4 | 1/3456 | 17/20736 |

For four distinct returning faces, YC17 gives the complete 24-order tube
bound5/216. Together these bound **each allowed fourth-order multiset**;
they do not yet give a summable local interaction norm for their global sum.
In particular p and q may be arbitrarily far apart in (2).

## 3. Keep spectator energy as an argument instead of dropping it

For a fixed word let X contain all its active cube and bridge factors.
On the retained bridge carrier write H_0=H_X+H_ext. The word multipliers
act only on X. Spectral resolution of H_ext gives the exact identity

\[
 \boxed{K_{\rm global}(z)=
     \int K_X(z-E)\otimes d\Pi_{\rm ext}(E).}            \tag{3}
\]

Equivalently, the three-resolvent heat kernel is the local three-time
kernel times exp[-(t1+t2+t3)H_ext]. Outside bridge factors are constant
at both ends and throughout; arbitrary spectator cube energies remain.
Thus replacing K_global by K_X(z) tensor I is wrong on general excited
spectators. On a reference source with H_ext Omega_ext=0 it is exact.

For YC17's tube, differentiating the three resolvents, or integrating
(t1+t2+t3) against the same heat bounds, gives

\[
 \sup_{E\ge0}\|\partial_E K_t(-E)\|
 \le\frac{29}{2592},\qquad
 \boxed{\|(K_{t,\rm global}(0)-K_{t,X}(0)\otimes I)\psi\|
       \le\frac{29}{2592}\|H_{\rm ext}\psi\|.}          \tag{4}
\]

This holds on Dom(H_ext). The constant is the derivative at z=0 of
4/(6-z)^3+2/[(6-z)^2(12-z)]. Equation(4) is a relative spectator-energy
bound, not a small absolute operator error on all states. Summing it over
all distant tubes without a connected cancellation would reintroduce a
volume cost; we do not perform that invalid step.

## 4. Why a disconnected Schur word is not an induced long-range force

At the completely free point choose faces p,q on disjoint factor supports.
The four interlaced p,p,q,q histories returning to the vacuum have energies
12,24,12, and each has final Haar weight1/16. The other two orders vanish.
Their raw scalar return coefficient is

\[
 \frac4{16\cdot12\cdot24\cdot12}=\frac1{13824}>0.        \tag{5}
\]

Nevertheless the two uncoupled systems have additive ground energy. In the
rank-one vacuum Schur equation E=-Sigma(E), each second-order contribution
is eta_p^2/[4(12-E)]. The mixed fourth-order term from the energy argument is

\[
 \frac2{48\cdot576}=\frac1{13824},                      \tag{6}
\]

which cancels minus(5) in E. This rank-one calculation is a declared scalar
control, not a replacement of the cube-retaining P_b above. It shows the
danger of extracting a zero-energy raw return and calling it an interaction.
Keeping the eigenvalue equation or using a connected creation exponent
implements the cancellation. No distant coupling is created in an exactly
decoupled Hamiltonian.

## 5. A complete local recursion on the actual ground carrier

Now use YC9/YC10/YC15's full excitation cut. Each factor i has its actual
reference ground Omega_i; I labels the excited factors, and P_I extracts
their component. H_I=sum_(i in I) h_i restricted to this component has
H_I>=g_* |I|, with g_*=1279511/4672512. Bridge floors are larger.
For c_I define the bounded creator chat_I=|c_I><Omega_I|, identity outside
I, and C=sum chat_I. These creators commute, with overlapping products zero.
The ground is written in intermediate normalization as Psi=exp(C)Omega.

Fix any profile |v_p|<=1, V=sum v_p W_p and H(zeta)=H_0-zeta V. Here zeta
is an **interface coupling parameter**, distinct from the spectral energy z
in(3), native response lambda and observation lambda. Write

\[
 c(\zeta)=\sum_{n\ge1}\zeta^n c_n,\qquad
 \mathcal B_l(d_1,...,d_l)_I=
 H_I^{-1}P_I[...[V,\widehat D_1],...,\widehat D_l]\Omega.  \tag{7}
\]

B_0 means H_I^-1 P_I V Omega. D_j is the summed creator for collection d_j.
Since creators commute, the nested-commutator expression is symmetric in
its creator arguments. The fixed point c=zeta H^-1 P e^-C V e^C Omega
therefore gives the **complete** recursion

\[
\begin{aligned}
 c_1&=\mathcal B_0,\\
 c_2&=\mathcal B_1(c_1),\\
 c_3&=\mathcal B_1(c_2)+\tfrac12\mathcal B_2(c_1,c_1),\\
 c_4&=\mathcal B_1(c_3)+\mathcal B_2(c_1,c_2)
                   +\tfrac16\mathcal B_3(c_1,c_1,c_1).
\end{aligned}                                                   \tag{8}
\]

All faces, repeated insertions, charged cube states and all electric harmonics
remain in (7)-(8). This is an exact operator recursion, not a numerical table
of all its interacting-cube matrix elements. The finite controls compare it
with independently derived Rayleigh-Schrodinger coefficients and the creation
logarithm through order four.

**Locality theorem for these coefficients.** A nonzero order-n monomial has
a connected witness consisting of at most n face supports in the factor
overlap graph. Its output I lies in their union (at most4n factors), but I
itself need not be connected. It depends on the reference operators in
that union only. If the monomial's faces split into disjoint factor sets,
turn off every other coupling: the Hamiltonian and intermediate-normalized
ground factorize, Psi=Psi_A tensor Psi_B. The commuting creation logarithm
then has C=C_A+C_B, so every mixed derivative is zero. This proves the
connected selection rule, including repeated faces. The nested-commutator
recursion gives the same conclusion inductively.

Spectators are in their reference grounds when evaluating (7); their energy
vanishes exactly, not by an approximation. Equations(3)-(4) still apply when
the target is instead the full operator on spectator excitations. The two
different cuts and targets must not be interchanged.

The bridge-vacuum part of(8) has precisely the multiset possibilities of
section1. Disjoint-support p^2 q^2 coefficients vanish in C. Overlapping
ones are retained by(8); no assumption that only tube terms remain is made.

## 6. A uniform bound on the complete omitted ground-return tail

The classical mechanism is the Kirkwood-Thomas/coupled-cluster construction
already credited in YC9/YC10; see also Yarotsky, arXiv:math-ph/0411042,
section2. The following explicit analytic bound is derived from YC15's
pinned constants, rather than from a finite-order truncation.

Use ||c||_*=max_i sum_(I containing i) |I| ||c_I||, r_*=1/128 and
rho=1/4320. The bounds for the fixed point depend only on absolute couplings
and Hilbert norms, so hold also for complex |zeta|<=rho and complex profiles
|v_p|<=1. The reference h_c and its ground remain the fixed real ones.
YC15 gives seed<=81|zeta|/20 and contraction at the endpoint

\[
 q_*\le\frac{28035072}{31987775}<1,\qquad
 r_*-({81\rho}/{20}+q_*r_*)
        \ge\frac{11417}{409443520}>0.                  \tag{9}
\]

The finite-volume map is holomorphic in zeta and in the creation vectors
(the fixed ground bra is not conjugated with the variable vector). Uniform
Picard convergence on the disk gives a holomorphic fixed point bounded by
r_*, with the same constants for every volume. Banach-valued Cauchy estimates
therefore give ||c_n||_*<=r_* rho^-n and, for x=|zeta|/rho<1,

\[
 \boxed{\left\|c(\zeta)-\sum_{n=1}^4\zeta^n c_n\right\|_*
       \le\frac1{128}\frac{x^5}{1-x}.}                 \tag{10}
\]

For |zeta|<=1/8640 the tail is <=1/2048; for |zeta|<=1/17280 it is
<=1/98304. These are conservative complete-tail bounds in the local
creation norm. They are not full many-body Hilbert-vector error bounds,
energy errors, or new spectral gaps. At the endpoint x=1 this particular
Taylor estimate does not converge; YC15's full fixed point still does.

## 7. The tube becomes a local connected creation coefficient

At free internal coupling, select the coefficient of the product of the
four *distinct* tube couplings. Every nonempty proper subset has a nonzero
bridge parity. In a product of disjoint creators, such odd bridge supports
cannot disappear; overlapping creators multiply to zero. Thus no creation-
log subtraction contributes to the bridge-vacuum part of this multilinear
coefficient. The lower energy corrections also have no such proper-subset
vacuum monomial. It equals H_0^-1 applied to YC17's returned source.

Both endpoint faces have electric energy12, so their product has energy24.
Consequently the two-cube component is

\[
 \boxed{(c_4)_{t,\,\mathrm{bridge\ vacuum}}
        =\frac{11}{3981312}W_AW_B\Omega_0.}              \tag{11}
\]

Its ordinary Hilbert norm is11/15925248, and its contribution at either
endpoint to the weighted local creation norm is11/7962624 before multiplying
the four couplings. This refines the target: 1/64 is the static conditional
record, 11/165888 the kinetic return source, and 11/3981312 the connected
ground-creation coefficient. None is an all-state scalar-potential coupling.

## 8. What is now closed, and what is next

The fourth-order bridge parity classification is complete. The full
ground-creation recursion through that order, its connected locality and
its volume-independent higher-order tail are controlled. The coefficients
on arbitrary interacting cubes are specified by full inverses, not all
numerically evaluated. The static/kinetic/creator distinction is retained
in the framework continuity map.

The next target is the transformed **excitation** operator: its connected
interactions, energy dependence and locality under a block update. Equation
(10) controls the ground exponent and cannot alone replace that operator
bound. No repeated spatial RG contraction, weak-coupling continuum trajectory
or 4D mass theorem is claimed. Existing lattice-gap certificates remain valid.

Reference: https://arxiv.org/abs/math-ph/0411042 .
Run `python3 -B physics/yc18/yc18_connected_fourth_order.py --check` and
`python3 -B -m unittest discover -s physics/yc18 -p 'test_*.py'`.
