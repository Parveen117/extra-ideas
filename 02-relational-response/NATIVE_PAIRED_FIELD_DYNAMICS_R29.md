# R29: native paired fields, constitutive closure and boundary residue

Research owner: Monty Dabas. Development: 1 October 2026.

R28 connected native H/K curvature to hidden exchange and a signed wave
response. R29 extracts two finite-difference fields from that same response.
Their coupled evolution, their relation to each other, a positive conserved
form and a local current are derived from the source transport. A native
compatibility residue distinguishes the actual source from the larger class
of paired-wave solutions. Zero-field backgrounds and the conserved residue
between two different backgrounds are then classified.

No electromagnetic field equation or constitutive coefficient is an input.
The field letters below denote event and address differences; their electric
and magnetic identification is not assumed. The R18 response target remains
explicit. Selection of that target as a physical process is not implied by
the following exact results within it.

## Native source, event pairing and field target

Retain R16's source roles, pairing, H,K and P=(I+H)/2, Q=(I-H)/2.
Let \(\chi_n=U^n[H,K]v/2\) be the signed curvature response of R28.
At even events and even word addresses set

\[
\phi_m(y)=\chi_{2m}(2y),\qquad (Tf)(y)=f(y-1),\quad
\delta=T-I,\quad\delta^\dagger=T^{-1}-I.
\tag{29.0}
\]

Here m counts pairs of source events, and y is a quotient of the existing
word address. No spatial coordinate, physical duration or differential
operator is supplied. The odd-address sector, when present, is a separate
copy of this construction. Finite sums use the native matching pairing;
dagger for T is proved by reindexing the finite sum.

The source two-event map and a derived local factor are

\[
V=U^2=\tfrac12\begin{pmatrix}T+I&T-I\\I-T^{-1}&I+T^{-1}\end{pmatrix},
\quad N=\tfrac12(P+T^{-1}Q)(H+K)
       =\tfrac12\begin{pmatrix}I&I\\T^{-1}&-T^{-1}\end{pmatrix}.
\tag{29.1}
\]

The two **field readouts** are explicitly constructed differences:

\[
\boxed{\mathsf e_m=\phi_{m+1}-\phi_m,\qquad
       \mathsf b_m=\delta\phi_m.}
\tag{29.2}
\]

Both have native two-role values. They are not two physical polarizations,
and their names do not identify H/K roles as electric/magnetic components.
The construction starts on finite address fields. Where stated, we also
use sequences constant outside a finite interval, keeping their two end
values as a finite record. Only their finite field differences are summed;
no infinite background norm is assigned. A cyclic count quotient is a
separately specified finite boundary target, not a chosen physical topology.

## Written results

### R29.1 — The source derives a local constitutive relation between the fields

Exact multiplication of the source blocks gives

\[
\boxed{V-I=N\delta,\quad N^\dagger N=NN^\dagger=I/2,\quad
       N^{-1}=2N^\dagger.}
\tag{29.3}
\]

Consequently every source-derived field pair satisfies

\[
\boxed{\mathsf e_m=N\mathsf b_m,\qquad
       \|\mathsf e_m\|^2=\tfrac12\|\mathsf b_m\|^2.}
\tag{29.4}
\]

N is determined by H,K and the source address shift; it is not a supplied
medium response or a fitted constant. This norm ratio belongs to the exact
readout normalization (29.2). Rescaling either readout changes the ratio;
it is not yet a physical impedance, speed or dimensionless coupling.

**Proof.** Factor T-I from every entry of V-I, using
\(I-T^{-1}=T^{-1}(T-I)\), to obtain N. The factor P+T^{-1}Q preserves
pairing because P,Q are disjoint and T is a shift. Anticommutation gives
\((H+K)^2=2I\), proving both dagger products and the stated inverse.
Insert \(\phi_{m+1}=V\phi_m\) into (29.2). This proves the constitutive
relation and norm ratio. In this paper "constitutive" names the derived
relation between the two constructed readouts; it is not a physical-medium
postulate. All identities are finite Laurent identities.

### R29.2 — The two native fields obey closed reversible induction equations

The source-derived fields obey

\[
\boxed{\mathsf b_{m+1}-\mathsf b_m=\delta\mathsf e_m,
\qquad \mathsf e_{m+1}-\mathsf e_m
             =-\tfrac12\delta^\dagger\mathsf b_{m+1}.}
\tag{29.5}
\]

These equations have no external forcing term. Their coefficient 1/2 is
the already derived R18 normalization. Once written as an update on any
finite pair, they define the reversible local map

\[
\mathcal M=
\begin{pmatrix}I-\tfrac12\delta^\dagger\delta&-\tfrac12\delta^\dagger\\
                \delta&I\end{pmatrix}.
\tag{29.6}
\]

Each field separately obeys

\[
f_{m+1}-2f_m+f_{m-1}=-\tfrac12\delta^\dagger\delta f_m.
\tag{29.7}
\]

The larger map (29.6) contains pairs which do not arise from the original
source. R29.3 gives the exact selection condition; a wave equation alone
does not certify source provenance.

**Proof.** Taking the address difference of the first field in (29.2)
gives the first equation. R28's signed wave identity at the middle event
\(m+1\) gives the second, since
\(-\delta^\dagger\delta=T+T^{-1}-2I\). Substitution yields (29.6).
The inverse is obtained without dividing by a difference operator:
\(\mathsf e=\mathsf e'+\delta^\dagger\mathsf b'/2\), then
\(\mathsf b=\mathsf b'-\delta\mathsf e\). Both directions have finite
range. Taking one further event difference in (29.5) gives (29.7) for
each field. No classical wave or induction law has been used as a premise.

### R29.3 — A conserved compatibility residue selects the native source sector

For any pair evolved by (29.6), define

\[
\mathsf c_m=\mathsf e_m-N\mathsf b_m.
\]

Then the complete evolution takes the triangular form

\[
\boxed{\mathsf c_{m+1}=V^\dagger\mathsf c_m,\qquad
       \mathsf b_{m+1}=V\mathsf b_m+\delta\mathsf c_m.}
\tag{29.8}
\]

Thus \(\|\mathsf c_m\|^2\) is conserved. The constitutive condition
\(\mathsf c=0\) is preserved exactly, and a nonzero compatibility residue
does not decay into that sector. Its influence on b is a calculable signed
source \(\delta\mathsf c\), not a stochastic noise assumption.

For finite-support fields on the unwrapped address line, a pair is derived
from a finite-support potential \(\phi\) exactly when

\[
\boxed{\mathsf c=0,\qquad\sum_y\mathsf b(y)=0.}
\tag{29.9}
\]

The potential is unique there. On a cyclic quotient of L addresses, the
same criterion holds, and the potential is unique up to one uniform native
two-role value. The image has rank 2L-2. For sequences with potentially
different constant end values, the zero-sum condition is replaced by the
boundary record in R29.7.

**Proof.** Unitarity of V and (29.3) give
\(V+V^\dagger=2I-\delta^\dagger\delta/2\) and
\(NV^\dagger=N+\delta^\dagger/2\). Substituting (29.6) into c proves
the first identity in (29.8); the second follows from
\(I+\delta N=V\). This proves norm conservation and sector invariance.

For a finite-support b with zero sum, define
\(\phi(y)=-\sum_{j\le y}\mathsf b(j)\). Both tails vanish and
\(\delta\phi=\mathsf b\). If c=0, (29.3) also gives
\((V-I)\phi=\mathsf e\). Conversely a finite difference telescopes to
zero total. A zero difference is a constant, and a finite-support constant
is zero, proving uniqueness. On a cycle the identical finite recursion
closes exactly when the sum is zero; its initial two-role value is free.
There are 2L potential coordinates and precisely two constant kernel
coordinates, giving the image rank. This is a native integration/recovery
constraint, not an assumed electromagnetic Gauss law.

### R29.4 — The paired update has an exact positive conserved field form

For arbitrary finite paired fields, define

\[
\mathcal E(\mathsf e,\mathsf b)=\|\mathsf e\|^2
 +\tfrac12\|\mathsf b\|^2+\tfrac12\langle\mathsf b,\delta\mathsf e\rangle.
\tag{29.10}
\]

It is invariant under (29.6) and has the positive native decomposition

\[
\boxed{\mathcal E=
 \tfrac12\|\mathsf b+\delta\mathsf e/2\|^2
 +\tfrac12\|\mathsf e\|^2
 +\tfrac18\|(T+I)\mathsf e\|^2.}
\tag{29.11}
\]

It vanishes only for the zero pair. On the source sector c=0 there are also
the conserved quantities

\[
\mathcal E_{\rm src}=\|\mathsf e\|^2+\tfrac12\|\mathsf b\|^2
                   =\|\mathsf b\|^2,
\quad
\mathcal E=\|\mathsf b\|^2-\tfrac18\|\delta\mathsf b\|^2,
\tag{29.12}
\]

with \(\|\mathsf b\|^2/2\le\mathcal E\le\|\mathsf b\|^2\).
The two forms are distinguished; neither is renamed physical EM energy.
The ordinary uncoupled sum in \(\mathcal E_{\rm src}\) is not generally
conserved off the source sector.

**Proof.** Put \(d=\delta\mathsf e\),
\(\mathsf b'=\mathsf b+d\),
\(\mathsf e'=\mathsf e-\delta^\dagger\mathsf b'/2\).
Expanding (29.10) at the primed pair cancels both squared
\(\delta^\dagger\mathsf b'\) terms and leaves
\(\|\mathsf e\|^2+\|\mathsf b'\|^2/2-\langle d,\mathsf b'\rangle/2\),
which equals (29.10). Completing a square and expanding
\(\|(T-I)\mathsf e\|^2+\|(T+I)\mathsf e\|^2=4\|\mathsf e\|^2\)
proves (29.11). Its terms are native matching squares.

For c=0, R29.1 and \(\mathsf b'=V\mathsf b\) prove conservation of
\(\mathcal E_{\rm src}\). The cross term is
\(\langle\mathsf b,(V-I)\mathsf b\rangle/2
=-\|\delta\mathsf b\|^2/8\), giving (29.12). Since the shift preserves
norm, \(0\le\|\delta\mathsf b\|^2\le4\|\mathsf b\|^2\), proving
the bounds. Commutation of V with delta gives conservation directly as well.

### R29.5 — The source field has an exact local current and a visible defect source

Write \(\mathsf b(y)=(u_y,v_y)\),
\(h_y=u_y+v_y\), \(k_y=u_y-v_y\), and
\(\rho_y=u_y^2+v_y^2\). Define the native link current

\[
\boxed{j_{y+1/2}=\tfrac14(h_y^2-k_{y+1}^2)
                         +\tfrac12h_y k_{y+1}.}
\tag{29.13}
\]

On the source sector it obeys exactly

\[
\boxed{\rho_{m+1}(y)-\rho_m(y)
 =j_m(y-1/2)-j_m(y+1/2).}
\tag{29.14}
\]

For a general paired field, the right side acquires the explicit term

\[
2\langle(V\mathsf b)_y,(\delta\mathsf c)_y\rangle
       +\|(\delta\mathsf c)_y\|^2.
\tag{29.15}
\]

Thus dropping the compatibility record can make the b reading look locally
sourced. The complete paired form (29.10) still conserves its energy. This
source term is a response-ledger effect, not an identified electric charge.

**Proof.** The two source components are
\(u'_y=(h_{y-1}+k_y)/2\), \(v'_y=(h_y-k_{y+1})/2\).
Expand their two squares and subtract
\(\rho_y=(h_y^2+k_y^2)/2\). Group the terms by their left and right
links to obtain (29.13)–(29.14). This includes the cross-address signed
interference term. Equation (29.8) gives (29.15) by expanding the square
of \(V\mathsf b+\delta\mathsf c\). Finite-support or periodic sums
telescope, yielding the corresponding global balance. No field-energy
current from classical electrodynamics has been supplied.

### R29.6 — Zero-field backgrounds are exactly the stationary uniform sector

On the address line or a finite cyclic count quotient,

\[
\boxed{\mathsf e=\mathsf b=0
 \ \Longleftrightarrow\ \phi(y)=q\text{ for every }y
 \ \Longleftrightarrow\ V\phi=\phi.}
\tag{29.16}
\]

Distinct q values have identical field readings. The field map forgets
exactly this uniform background and nothing else within the stated sequence
class. A uniform nonzero q can still have a nonzero local native curvature
response:

\[
\|[H,K]q\|^2=4\|q\|^2.
\tag{29.17}
\]

For a finite perturbation \(\eta\), the fields of q+eta equal those of
eta, and \(V^m(q+\eta)=q+V^m\eta\). Disturbances therefore propagate
on a stationary zero-field background without selecting its amplitude.
No infinite background energy is used. On finite-support potentials alone,
the sole uniform background is zero. Stationarity here is after two source
events; a single-event role transformation is a different target.

**Proof.** A zero address difference means \(\phi(y-1)=\phi(y)\), hence
a constant q. For a constant, T acts as I, and (29.1) gives V=I. Conversely
\((V-I)\phi=N\delta\phi=0\); the local inverse of N forces
\(\delta\phi=0\). Applying this argument to the difference of two
potentials identifies the field-map kernel. R28.1 proves (29.17).
Linearity and stationarity prove the perturbation statement. This derives
uniform-background blindness of this observer. It is not a claim of an
arbitrary local gauge symmetry or a physically identified vacuum state.

### R29.7 — An interface between zero-field backgrounds has a conserved boundary residue

Suppose phi is constant outside a finite interval, with end records
\(q_-\) and \(q_+\). Its b field is finite and

\[
\boxed{\mathcal Q:=\sum_y\mathsf b(y)=q_--q_+,
       \qquad\mathcal Q_{m+1}=\mathcal Q_m.}
\tag{29.18}
\]

Every finite b and left end record reconstruct such a source potential:

\[
\phi(y)=q_--\sum_{j\le y}\mathsf b(j),\qquad
q_+=q_--\mathcal Q,\qquad\mathsf e=N\mathsf b.
\tag{29.19}
\]

The native continuation changes only a finite interval at every finite
event, so both end records persist. A nonzero boundary residue cannot be
erased by this continuation. Zero total residue still allows nonzero
propagating fields: \(\phi=\delta_0 e_0\) has
\(\mathsf b=(\delta_1-\delta_0)e_0\), total zero and field energy two.

For a single interface \(q_-=e_0,q_+=0\), take
\(\mathsf b_0=\delta_0e_0\). One paired source event gives

\[
\mathsf b_1=\tfrac12\delta_1e_0
 +\tfrac12\delta_0(e_0+e_1)-\tfrac12\delta_{-1}e_1.
\tag{29.20}
\]

Its total residue is still e0 and its field energy is still one.

If b occupies at most M addresses, its source field energy satisfies
\(\|\mathsf b\|^2\ge\|\mathcal Q\|^2/M\). Equality is attained
by M equal consecutive values \(\mathcal Q/M\). Hence a fixed boundary
residue alone does not select a nonzero universal energy or particle mass;
different source preparations may spread the interface over different
numbers of native addresses. This is not a loss of energy during any one
conserved evolution.

**Proof.** Telescoping \(\phi(y-1)-\phi(y)\) proves (29.18), and finite
partial sums prove (29.19). Summing the first equation in (29.5) shows that
the total b is conserved even for the larger paired-field map. In the
source sector the compact increment N b leaves both tails unchanged.
The examples follow by the source Laurent matrix (29.1). Finally expand
\(\sum_{j=1}^M\|\mathsf b_j-\mathcal Q/M\|^2
=\sum_j\|\mathsf b_j\|^2-\|\mathcal Q\|^2/M\).
Native positivity gives the bound and its equality case. The conserved
two-role quantity Q is a boundary residue; quantized electric charge and
its interaction law have not been identified by this calculation.

## Evidence, lineage and physical interpretation

Seven written proofs are accompanied by exact source-derived Laurent
identities, finite evolution and reconstruction checks, native word replay,
negative controls and pinned source hashes. The proof graph excludes
imported/admitted/comparison premises. Its metadata gate is not a semantic
proof assistant. The earlier R28 and foundation chain are replayed without
rewriting their frozen evidence.

| Pinned source | Premise status and use |
| --- | --- |
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native-derived:** signed/refined coefficients, roles, matching pairing and H/K. Finite equality, counting and induction remain the explicitly declared proof infrastructure. |
| [R18.1–R18.5](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | **Native-derived within an explicit target:** source-word address, signed aggregation, normalization, finite transport and wave identity. Physical selection of the target is not supplied. |
| [R28.1–R28.6](https://github.com/Parveen117/extra-ideas/blob/9ed6ae27d3ddb005ed7f4f3f27417036729698fd/02-relational-response/NATIVE_CUT_CURVATURE_PROPAGATION_R28.md) | **Native-derived:** curvature response, transition coupling, hidden memory, signed propagation and the boundary against inferring unique dynamics from algebra alone. R29 changes the readout, not that physical-selection boundary. |
| [Canonical RKF algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | **Native exact arithmetic and equality replay:** unchanged engine; scoped use, no whole-engine or formal-proof-assistant claim. |
| [R28 comparison audit](https://github.com/Parveen117/extra-ideas/blob/9ed6ae27d3ddb005ed7f4f3f27417036729698fd/04-operator-evolution/R28_DERIVATION_LEDGER.json) | **Comparison only:** prior declared zero-cut topology, declared-carrier information response, and the admitted smooth/tangent/metric Riemann adapter remain excluded from the new proofs. Source inclusion does not promote their imported inputs. |

The completed advance is a source-derived paired field system with its own
constitutive operator, positive conservation, local current, source-sector
constraint, stationary zero-field backgrounds and conserved interface
residue. This makes the proposed cut-medium interpretation mathematically
more specific. Three physical spatial directions, transverse EM response,
an electromagnetic source/charge identification, local gauge structure and
physical units remain to be derived or operationally justified. Calling
the two differences electric and magnetic would not establish those results.
No value of physical c, alpha or particle mass is certified here.
