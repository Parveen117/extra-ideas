# YC17 — a four-face bridge return that survives at the free cube point

9 October 2026. Continues YC16 at `fa07e53`. Frozen packets are unchanged.

**Result.** The two-bridge faces in YC12's actual cube tiling partition into
four-face tubes. Bridge parity kills every nonempty proper subset of one
tube, even with the full reference resolvents inserted. Returning all four
faces produces a genuine two-cube channel. Its static Haar coefficient is
1/64, whereas its fourth-order kinetic return on the free reference vacuum
has coefficient **11/165888**. On arbitrary correlated cube states the sum
of its 24 ordered full-resolvent words has norm at most **5/216**.

These are one selected, now explicitly resolved, part of the next spatial
return. The full lattice gap remains YC15/YC16's **>9/40 through
eta<=1/4320**, internal theta_c<=1. No enlarged window, complete fourth-order
effective Hamiltonian, iterated blocking or continuum gap is claimed.

## 1. Carrier, retained cut and the next question

Use disjoint twelve-link cubes and individual bridge links on a spatial
three-torus with even sides >=4, as in YC12. The normalization is normalized
SU(2) Haar measure, W_p=Tr(U_p)/2 and fundamental electric energy 3 per link.
Let

\[
 H_0=\sum_c h_c+\sum_{b\ {\rm bridge}}(-\Delta_b),\qquad
 h_c=H_c-E_c\ge0.
\]

The actual correlated cube Hamiltonians and all their charged states are
retained. Define **P_b** to project every bridge to its constant function,
while acting as the identity on every cube; Q_b=I-P_b. This is not the
rank-one global vacuum cut. Both projections commute with the gauge action.
The comparison norm below is the full product L2 norm, not a finite trial
norm. The source and target are explicit: external W_p insertions leave
P_b, propagate through Q_b, and return an operator on the complete cube
carrier. Nonlocal spectator energies stay in H_0.

YC16 changes the reference by exact compensation. Here we resolve an
additional coefficient using its original h_c reference, which is still
available; we do not silently omit YC16's counterterms or assert that this
coefficient already equals a blocked Hamiltonian on its changed reference.

For z<3, set

\[
 G_z=Q_b\big[Q_b(H_0-z)Q_b\big]^{-1}Q_b.                 \tag{1}
\]

It exists because a nonconstant bridge costs at least 3. For one tube with
distinct faces p_0,...,p_3, define its multilinear fourth-order return

\[
 K_t(z)=\sum_{\pi\in S_4}
 P_bW_{p_{\pi(4)}}G_zW_{p_{\pi(3)}}G_z
 W_{p_{\pi(2)}}G_zW_{p_{\pi(1)}}P_b.                    \tag{2}
\]

Here pi(1),...,pi(4) list 0,...,3. Reverse orders pair as adjoints, so
K_t(z) is self-adjoint. For H=H_0-sum eta_p W_p, this is the coefficient of
prod_(p in t) eta_p in the Schur **return**; the Schur effective H subtracts
that return. This does not say K_t is a positive operator. On each finite
lattice the coefficient is an actual resolvent derivative at zero external
couplings. A uniform convergence radius for the *entire* bridge Schur
expansion does not follow from the selected coefficient estimates here.

## 2. The geometric selection rule, without a harmonic truncation

Take the four bridges joining a face of cube A to the facing face of its
neighbor B. Enumerate the vertices around the square by i modulo 4; bridge
u_i joins the corresponding vertices. The side faces have internal links
a_i in A and b_i in B, traversed in corresponding cycle orientations:

\[
 W_i=\tfrac12\operatorname{Tr}(a_i u_{i+1}b_i^{-1}u_i^{-1}). \tag{3}
\]

Their bridge-centre signatures are {u_i,u_(i+1)}. Every bridge has exactly
two two-bridge faces, the adjacent sides of this square. Hence the complete
two-bridge-face signature graph is a disjoint union of four-cycles, one for
each positive-direction neighboring cube interface. There are three tubes
per cube in total counting and six incident tubes at each cube. With only
two cubes in a periodic direction, its two interfaces remain distinct;
they are not merged because their endpoints coincide.

Let Z_b flip u_b to -u_b. H_0, G_z and P_b commute with every Z_b, while
W_i has the above odd signature. A word returning to P_b therefore vanishes
unless the XOR of its signatures is empty. In a four-cycle a subset has
even degree at every vertex iff it is empty or the whole cycle. Therefore:

* no proper nonempty subset of the four distinct side faces can return;
* among distinct **two-bridge** faces, minimal nonempty returning sets are
  exactly these four-face tubes;
* at three insertions, even allowing repetition, a two-bridge-only return
  is zero (odd multiplicity support cannot be empty or a union of tubes).

This does not classify words containing four-bridge plaquettes, repeated
fourth-order faces or higher-order returns. Parity is only a necessary rule
in general; the nonzero calculations below establish this particular channel.

## 3. Static closed return: exact boundary information

Put W_A=Tr(a_0a_1a_2a_3)/2 and W_B=Tr(b_0b_1b_2b_3)/2. Fundamental
SU(2) character gluing, with all four bridge integrations performed once,
gives the pointwise conditional identity

\[
 \boxed{\int\prod_{i=0}^3 W_i\,du_0du_1du_2du_3
                 =\frac1{64}W_AW_B.}                  \tag{4}
\]

One direct proof uses quaternion components. Write W_i=u_i^T T_i u_(i+1).
Each bridge appears exactly twice; its second moment is I_4/4. Contracting
the four moments yields Tr(T_0T_1T_2T_3)/256. Left/right quaternion
multiplication gives Tr(T_0...T_3)=4 Re(a_0...a_3) Re(b_0...b_3), proving
(4). No small-link or weak-field approximation is involved. The code checks
the full four-linear one-cell tensor identity on all 4^4 basis choices;
cyclic matrix multiplication then proves the general contraction.

As a normalization cross-check, close the two ends with W_A and W_B and
integrate the cube links at Haar. Each squared face has mean 1/4, giving
(1/64)(1/4)^2=1/1024: exactly YC11's six-face closed-cube moment. This
uses the same twelve participating links, not duplicated shared links.

Equation (4) is a configuration record, not a kinetic resolvent coefficient.

## 4. Full operator bound with all spectator energies

After one and three distinct insertions, exactly two bridges are odd.
After two, there are two odd bridges for adjacent faces and four for
opposite faces. Odd SU(2) harmonics have electric energy at least 3, so
the full H_0 floors at these steps are respectively

\[
 (6,6,6)\quad\hbox{or}\quad(6,12,6).                   \tag{5}
\]

All cube and spectator Hamiltonians are nonnegative; retaining them can
only increase these floors. This argument does not project their spectra
to a finite set. Also ||W_i||<=1 and

\[
 \|W_iP_b\|=\|P_bW_i\|=1/2,                            \tag{6}
\]

because integrating W_i^2 over its two constant bridges gives 1/4 pointwise
in its internal links. Of the 24 orders, 16 have an adjacent initial pair
and eight an opposite initial pair. Thus for 0<=z<3,

\[
 \|K_t(z)\|\le
 \frac4{(6-z)^3}+\frac2{(6-z)^2(12-z)},\qquad
 \boxed{\|K_t(0)\|\le\frac5{216}.}                     \tag{7}
\]

The zero-energy bound also holds for z<=0. The rational expression has a
formal larger domain z<6 on the selected parity sectors, but (1) has only
been defined globally for z<3; no claim crosses that distinction.

For |eta_p|<=eta, the sum of **selected tube norm costs incident to a cube**
is at most 6(5/216)eta^4=5eta^4/36, independently of the number of cubes.
With excited spectators present, K_t need not be a two-cube-supported
operator, because its resolvents include their energies. This is an
incidence-summed coefficient estimate, **not** an already justified local
interaction norm for a new all-volume cluster theorem. In particular it
must not be inserted into YC10 as beta while ignoring spectator dependence.

## 5. Exact kinetic vacuum source at free internal coupling

Now set every internal theta_c=0 and apply (2) at z=0 to the constant
reference Omega_0. At the final return, a bridge that has appeared twice
must be in its constant representation: no later face in that word acts
on it. Its final projection commutes backwards through all later free
resolvents and the later faces. It can therefore be imposed as soon as
that bridge is completed. This removes its adjoint branch exactly, rather
than assuming the intermediate product is a single harmonic.

Each of the 2k internal edges used after k distinct faces occurs once and
has fundamental energy 3. Each still odd bridge contributes another 3.
The surviving branch therefore has intermediate energies

| Initial pair | Number of orders | Three resolvent energies |
|---|---:|---|
| Adjacent | 16 | 12, 18, 24 |
| Opposite | 8 | 12, 24, 24 |

The remaining final bridge contraction is (4), so

\[
 \boxed{K_t(0)\Omega_0
 =\frac1{64}\left(\frac{16}{12\cdot18\cdot24}
                 +\frac8{12\cdot24\cdot24}\right)W_AW_B\Omega_0
 =\frac{11}{165888}W_AW_B\Omega_0.}                     \tag{8}
\]

The face functions have zero Haar mean and norm 1/2. This is a nonzero
q_A tensor q_B source, with exact norm **11/663552**. Spectator cubes and
bridges are in their free grounds for this source calculation only.

This explains why YC16's free-point silence is not closure: its connected
second-order source C_p Omega_0 is zero, but the four-distinct-face channel
(8) survives at fourth interface order. Neither that silence nor absorbing
one-cube feedback allows this channel to be discarded. A source coefficient
at theta=0 and a full-operator bound on theta_c in[0,1] are different results;
we do not extend (8) unchanged to interacting cube grounds. We also do not
replace K_t on excited states by multiplication by the right side of (8).

## 6. The owner's commutator distinction, precisely typed

The owner's lead is that off-diagonal entries alone do not make curvature.
In a coordinate frame with D_mu=partial_mu+A_mu,

\[
 [D_\mu,D_\nu]=F_{\mu\nu},\qquad
 F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu+[A_\mu,A_\nu].
                                                               \tag{9}
\]

In a moving frame, subtract the frame bracket from the commutator; the
components include -c_(mu nu)^rho A_rho. The noncommuting object in the
general statement is the **covariant transport**, not a selected pair of
off-diagonal matrix entries. For constant coordinate connections only,
(9) reduces to the connection-matrix commutator. Tong's *Gauge Theory*,
section 2.1, uses the equivalent Hermitian convention D=partial-iA.

Two exact controls keep the useful lead from becoming an incorrect rule.
First, UP8/YC13 already have a=ell R dphi: all its matrix components commute,
yet F=R d ell wedge dphi can be nonzero. Second, for a Maurer-Cartan form
L=g^-1 dg, A=L is flat even when [L_mu,L_nu]!=0; the derivative cancels
the commutator. The native constant-share connection A=pL instead has
F=p(p-1)L wedge L. These are conditional identities for that connection,
not a universal equivalence between response asymmetry and YM field strength.

For the present lattice calculation W_p are scalar multiplication operators
and commute with one another. Their commutation does not kill (8): the
electric operator does not commute with them. Already at the free point,
[H_0,W_p]1=12W_p. The kinetic denominators in (8), absent in (4), record this
distinction. Neither static covariance nor matrix off-diagonality alone
is a quantum mass-gap criterion.

Reference: David Tong, *Gauge Theory*, section 2.1,
https://davidtong.org/pdfs/teaching/gauge-theory/gauge.pdf .

## 7. Claim boundary and next measurable obligation

The exact tensor checks, parity enumeration and rational arithmetic support
the written proofs; this is not a formal proof-assistant verification.
Frozen predecessor hashes are replayed. There is no numerical Taylor
remainder for the whole bridge expansion and no stronger lattice gap here.
Mixed four-bridge words, repeated-face returns, ground-energy dependence,
transport to the absorbed reference and spectator-dependent nonlocal terms
remain in the original full problem; none is claimed to vanish.

The next bounded task is to place this tube channel together with those
other fourth-order channels in a local expansion whose spectator dependence
and remainder have a summable bound. The larger gates are controlled
repeated spatial blocking, positive physical gap on a continuum scale
trajectory, and construction of the limiting quantum field theory. These
are proof obligations, not a known count of remaining development stages.

Run `python3 -B physics/yc17/yc17_four_face_bridge_return.py --check` and
`python3 -B -m unittest discover -s physics/yc17 -p 'test_*.py'`.
