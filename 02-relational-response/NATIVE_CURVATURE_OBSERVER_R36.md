# R36: response-derived geometry, a curvature-complete observer and information current

Research owner: Monty Dabas. Development: 2 October 2026 (India).

R35 supplies an isotropic signed-envelope sector of the retained R32
transport. R36 now constructs a way to read that sector using the native
curvature lineage of R28–R30. A source response determines a
preparation-independent phase metric and its dual count length. The same
response operators determine an oriented curvature. Its cut contrast
recovers exactly the information hidden by one native cut.

The complete paired reading has the exact source dynamics, not just a
scalar wave equation with extra initial states. In its controlled leading
flow, the local information density and current obey

\[
\partial_t\rho+\partial_xj_x+\partial_yj_y=0,\qquad
\boxed{j_x^2+j_y^2\le\rho^2/2.}
\]

More precisely, the deficit from this bound is a derived mismatch of the
two retained readings. No free coupling or statistical noise is supplied.
This connects the R35 phase coefficient to a native information-current
bound and a complete signed observer.

These fields remain readouts of one constructed source. The current law
is for the derived leading flow; it is not an exact finite-lattice current
law. R35's faster nonzero support front remains. Neither a physical field,
physical rods/clocks, energy/action units nor c, h or alpha is identified
by these results.

## Source order and the existing curvature-to-Riemann bridge

Use the unchanged R16 role/count construction, R20 tuple matching and R32
retained interface. R35 bijectively groups its signed carrier into sixteen
roles per component address. One component step is two R30 coarse steps;
one Z block is eight original source events. With component shifts S_x,S_y,

\[
Z^\dagger Z=I,\qquad Z+Z^\dagger=2I-L_*,\qquad
L_*=\tfrac12(\lambda_x+\lambda_y)-\lambda_x\lambda_y/16,
\quad\lambda_i=(S_i-I)^\dagger(S_i-I).
\tag{36.0}
\]

The exact Laurent coefficients of Z are constructed by R35.2 from the
original H/K source. Its first component-shift jets give self-dagger maps
A_x,A_y with

\[
A_x^2=A_y^2=I/2,\qquad A_xA_y+A_yA_x=0.
\tag{36.1}
\]

Scalar phase evaluations use only the iota and completion earned in R16,
and the native factorial Euler theorem F00-E. They are finite coefficient
probes, not an assumed Fourier basis. Matching norms and real parts are
those of the derived cut coefficients; ordinary complex/Hilbert structures
are not primitive inputs.

The earlier [R10](NATIVE_CURVATURE_DESCENT_R10.md) remains the precise
native-to-Riemann bridge. It proves connection descent after a smooth
tangent observer and a compatible torsion-free metric connection are
specified. Its smooth/tangent/Levi-Civita adapter is explicitly admitted.
R36 uses the native operator and information side of the programme to
construct a new observer and response metric. It does not silently import
R10's adapter or claim those new objects already satisfy its tangent
connection contract. Native protocol curvature and a Riemann tensor retain
their distinct types.

## Written results

### R36.1 — Actual source distinguishability derives a phase metric and dual count length

Let \(Z(q)=Z(\operatorname{Exp}_\Sigma(\iota q_x),
\operatorname{Exp}_\Sigma(\iota q_y))\). For any nonzero finite
sixteen-role coefficient v, define its relative response cost

\[
\mathcal I(q;v)=\frac{\|(Z(q)-I)v\|^2}{\|v\|^2}.
\]

This is squared signed distinguishability in the native matching pairing;
it is not assumed Shannon entropy, a probability or Fisher information.
Writing \(r=(1-\operatorname{Cos}_\Sigma q_x)/2\) and
\(s=(1-\operatorname{Cos}_\Sigma q_y)/2\), one has exactly

\[
\boxed{\mathcal I(q;v)=2(r+s)-rs,}
\qquad
\boxed{\gamma(a,b)=\tfrac12(a_xb_x+a_yb_y).}
\tag{36.2}
\]

Here gamma is the quadratic response at the zero carrier:
\(\mathcal I(\epsilon a;v)=\epsilon^2\gamma(a,a)+O(\epsilon^4)\),
with its bilinear value obtained by polarization. It is independent of v.
It is a metric on the phase-probe directions, not yet a spacetime metric
or a claim that the Hessian stays positive at every finite phase.

For a dual component displacement d, native optimization gives

\[
\boxed{\mathscr L(d)^2
 =\sup_{a\ne0}\frac{(a_xd_x+a_yd_y)^2}{\gamma(a,a)}
 =2(d_x^2+d_y^2).}
\tag{36.3}
\]

Thus the dual count-length form is \(\gamma^{-1}=2I_2\). The phase
form and its dual must not be interchanged. This is a constructed
information length from the specified response, not a material rod.

**Proof.** R35.1 gives \((Z-I)^\dagger(Z-I)=L_*\) coefficientwise.
After the native phase substitution, L_* is the scalar
\(\ell=2(r+s)-rs\); multiplying by v proves the exact cost and
its preparation independence. The native factorial cosine series gives
r=epsilon^2 a_x^2/4+O(epsilon^4), and similarly for s. Its tail is
controlled as in R35.4–R35.5. This proves the quadratic value, and finite
polarization proves gamma. No statistical weights enter.

The finite algebraic identity
\(|a|^2|d|^2-(a\cdot d)^2=(a_xd_y-a_yd_x)^2\ge0\)
proves the bound in (36.3); a=d attains it when d is nonzero. When d=0
every ratio vanishes. This proves the dual form without importing a
classical variational or representation theorem.

### R36.2 — The same response products derive an oriented curvature and its area law

Put \(A(a)=a_xA_x+a_yA_y\) and define

\[
J=[A_x,A_y]=2A_xA_y,\qquad H_c=A_x+A_y,
\qquad K_c=A_x-A_y.
\]

Then

\[
J^\dagger=-J,\quad J^2=-I,\quad
H_c^2=K_c^2=I,\quad H_cK_c=-K_cH_c,\quad J=K_cH_c.
\tag{36.4}
\]

The complete operator-valued response product is

\[
\boxed{A(a)^\dagger A(b)=\gamma(a,b)I
               +\tfrac12\det(a,b)J.}
\tag{36.5}
\]

Its symmetric part gives the metric and its antisymmetric part retains
ordered response. If \(\mathcal F(a,b)=[A(a),A(b)]\), then

\[
\boxed{\mathcal F(a,b)=\det(a,b)J,\qquad
\mathcal F(a,b)^\dagger\mathcal F(a,b)
 =4\{\gamma(a,a)\gamma(b,b)-\gamma(a,b)^2\}I.}
\tag{36.6}
\]

This is an exact native area/curvature relation for the derived response
target. Orientation reversal changes F while preserving its square.
It is not a Riemann sectional-curvature equation.

For the derived involutions, the role curvature and actual protocol loop are

\[
F_c=[H_c,K_c]=-2J,\qquad F_c^\dagger F_c=4I,
\qquad H_cK_cH_cK_c=-I.
\tag{36.7}
\]

Thus R28's cut-curvature law reappears inside the R35 propagation response;
no new H/K or iota is installed as a primitive.

**Proof.** Apply (36.1) to the displayed sums and products. Reversing
the order of A_x,A_y changes the sign, giving J-dagger=-J and
J-squared=-4A_x^2A_y^2=-I. The H_c,K_c relations and J=K_cH_c follow
by expansion. Expanding A(a)A(b) separates its dot and determinant
coefficients, proving (36.5). Subtract the reverse product and square it
with the dagger to get (36.6), using the two-coordinate determinant
identity from R36.1. Finally (36.4) proves (36.7) by finite multiplication.
These are actual ordered response maps and a four-arrow source-law
protocol; no spatial connection or small-loop approximation is presumed.

### R36.3 — Curvature supplies exactly the missing cut and a minimal complete observer

Let \(P=(I+H_c)/2\), \(Q=I-P\). Define two signed readouts
and their decoder by

\[
\boxed{a=Pu,\qquad b=PF_cu/2=-PJu,\qquad
       u=a+Jb.}
\tag{36.8}
\]

Both a and b belong to the same eight-role plus-cut image. The map
\(\mathcal O u=(a,b)\) is a bijective pairing isometry onto
\(\operatorname{im}P\oplus\operatorname{im}P\):

\[
\boxed{\|a\|^2+\|b\|^2=\|u\|^2.}
\tag{36.9}
\]

Two readings are minimal among linear readouts of this fixed rank-eight
cut when the target is the complete sixteen-role coefficient. This is
a target-relative minimum, not a universal physical observer count or
a replacement of the earlier spectral papers' different targets.

**Proof.** H_c is a self-dagger involution, so P,Q are disjoint
pairing cuts. JH_c=-H_cJ gives PJ=JQ and J-dagger PJ=Q.
Consequently the second readout has dagger square Q, proving (36.9).
Also J P J=-Q, and hence a+Jb=Pu-JPJu=u. Conversely, for any a,b
in im(P), the two readings of a+Jb are exactly a,b. This proves the
inverse, onto property and isometry without an additional norm choice.

J is an isomorphism between the plus and minus cut images; their
dimensions therefore agree and add to sixteen. Each is eight. One
rank-eight linear readout cannot inject a sixteen-role carrier; two
readings attain the bound. In particular Pu=0 can hold for a nonzero u
while its curvature reading b is nonzero. Reading only the two squared
sizes cannot perform this reconstruction: u and -u have the same two
intensities and different signed states. Signed curvature information is
retained, not replaced by an independently sampled noise channel.

### R36.4 — The complete observer has exact closed dynamics and the omitted cut has memory

Apply O pointwise to a finite component field. Its inverse is
\(\mathcal D(a,b)=a+Jb\), and the observed transport is

\[
\boxed{\mathcal M=\mathcal O Z\mathcal D
 =\begin{pmatrix}
 PZP&PZJP\\
 -PJZP&-PJZJP
 \end{pmatrix}}
\tag{36.10}
\]

on the two plus-cut targets. It is local, reversible and pairing
preserving, with inverse \(\mathcal O Z^\dagger\mathcal D\).
Since L_* is a scalar shift polynomial,

\[
\mathcal M+\mathcal M^\dagger=2I-L_*,\qquad
\boxed{\binom{a_{n+1}}{b_{n+1}}-2\binom{a_n}{b_n}
       +\binom{a_{n-1}}{b_{n-1}}
       =-L_*\binom{a_n}{b_n}.}
\tag{36.11}
\]

Both readouts are therefore actual compatible source fields with the
same finite wave operator. The full curvature field w=Ju instead has
\(Z_J=JZJ^\dagger\); it has the same exact scalar wave relation.
These are equivalent readings of one source, not independently selected
physical fields whose universality has been established.

One cut alone does not have a closed one-step law for general source
inputs. Its exact retained-memory equation is

\[
\boxed{a_{n+1}=M_{11}a_n+
 \sum_{k=0}^{n-1}M_{12}M_{22}^kM_{21}a_{n-1-k}
 +M_{12}M_{22}^n b_0.}
\tag{36.12}
\]

All kernel coefficients and the initial residue come from the original
transport. No random force or additional event law is assumed.

**Proof.** The identities D O=I and O D=I on the stated targets give
exact intertwining, the inverse and preservation of the matching form.
The matrix in (36.10) is direct composition of (36.8); all its entries
have the original finite support. O and D commute with scalar shifts.
Transport (36.0) through them to get (36.11). Conjugating by J proves
the curvature-field assertion. The first-order compatibility law always
remains; arbitrary initial pairs from a larger scalar second-order wave
solution space have not been admitted.

For a direct failure of one-cut closure, P A_x Q=P K_c Q/2 is nonzero,
because K_c maps the nonzero minus cut bijectively to the plus cut. It is
the first phase jet of PZQ, so PZQ is not the zero Laurent operator.
There are finite hidden inputs with Pu_0=0 but PZu_0 nonzero; the zero
input has the same initial readout and a different next one. Thus even
a putative function of the current a alone cannot describe all source
inputs. Iterating b_{n+1}=M_{21}a_n+M_{22}b_n and substituting into the
first row gives (36.12) by finite induction, applying R28's elimination
mechanism to the newly constructed observer. Hidden memory is derived
here from a specific omitted cut, not postulated as unexplained noise.

### R36.5 — Native information current has the same sharp cone and an exact retained-mismatch deficit

For a coefficient u define the principal density and currents

\[
\rho=\langle u,u\rangle,\qquad
j_i=\langle u,A_i u\rangle.
\]

The j_i are real by self-dagger A_i. In the complete observer coordinates
of (36.8), put \(r_a=\|a\|^2\), \(r_b=\|b\|^2\) and
\(h=\operatorname{Re}\langle a,b\rangle\). Then

\[
\boxed{\rho=r_a+r_b,\qquad
j_x=(r_a-r_b+2h)/2,\quad j_y=(r_a-r_b-2h)/2,}
\tag{36.13}
\]
\[
\boxed{\rho^2/2-j_x^2-j_y^2
       =2(r_a r_b-h^2)\ge0.}
\tag{36.14}
\]

For nonzero density this proves

\[
\boxed{|j|/\rho\le1/\sqrt2,
       \qquad\mathscr L(j/\rho)\le1.}
\tag{36.15}
\]

The second expression uses the response-dual count length from R36.1;
it is the same bound in that constructed length convention. It does not
fix a physical metre or second. The bound is sharp: proportional native
real readings attain it, including b=0. Independent orthogonal readings
of equal nonzero norm give zero principal current with positive density.
Curvature J reverses both currents while preserving density.

There is an actual local conservation law for the native leading-flow
fields constructed as finite sums

\[
u(t,x,y)=\sum_{\nu=1}^{N}
 \operatorname{Exp}_\Sigma[-\iota(k_{\nu x}x+k_{\nu y}y)]
 \operatorname{Exp}_\Sigma[\iota t A(k_\nu)]v_\nu.
\tag{36.16}
\]

They satisfy

\[
\partial_tu=-A_x\partial_xu-A_y\partial_yu,
\qquad
\boxed{\partial_t\rho+\partial_xj_x+\partial_yj_y=0.}
\tag{36.17}
\]

The parameters are native completed count/phase variables. This is a
constructed finite-phase limit class, not an assumed complete Fourier
expansion or a continuum theorem for every lattice field. No finite
global norm is assigned to a plane phase. The density is matching
information, not an identified physical energy or probability density.

**Proof.** Under O, H_c acts by (a,b)->(a,-b), K_c by (a,b)->(b,a),
and J by (a,b)->(-b,a). Since A_x=(H_c+K_c)/2 and
A_y=(H_c-K_c)/2, pairing these actions with (a,b) gives (36.13).
Squaring the two currents and subtracting from rho-squared/2 gives
(36.14). For a nonzero a, expand the native positive square
\(\|b-(h/r_a)a\|^2=r_b-h^2/r_a\). This proves its sign; if a=0
the assertion is immediate. Thus no imported Cauchy-Schwarz theorem is
needed. Equality holds exactly when this real matching residual is zero,
with the evident zero-cut cases. Proportional real readings attain the
bound and realize every current direction by the native circle
parametrization, with the zero-first-reading endpoint included. The
eight-dimensional plus space contains orthogonal equal-norm nonzero
readings, proving the zero-current example. The displayed action of J
proves current reversal.

F00-E and R35.5 justify derivatives of each scalar/matrix factorial series
in the finite sum (36.16). Each term has time derivative iota A(k)u and
spatial derivatives -iota k_i u, giving the first equation of (36.17).
Differentiate the finite matching square. Self-dagger, constant A_i give
\(\partial_t\rho=-\sum_i\partial_i\langle u,A_i u\rangle\).
This proves the second equation and its local interpretation without
assuming a classical field equation. Equations (36.14)–(36.17) are not
asserted to be the exact one-block current of Z.

### R36.6 — The paired and curvature fields inherit a controlled common limit, with a finite correction

On the paired plus-cut space write

\[
H_o(a,b)=(a,-b),\qquad K_o(a,b)=(b,a),\qquad
\widetilde A(k)=\tfrac12[(k_x+k_y)H_o+(k_x-k_y)K_o].
\tag{36.18}
\]

This generator is derived by the exact observer intertwining. It obeys
\(\widetilde A(k)^2=|k|^2I/2\). For
\(K=|k_x|+|k_y|\), \(0\le\epsilon K\le1/2\),

\[
\boxed{\|\mathcal M(\epsilon k)^n-
 \operatorname{Exp}_\Sigma[\iota n\epsilon\widetilde A(k)]\|
 \le89n\epsilon^2K^2.}
\tag{36.19}
\]

The curvature field has leading generator -A(k), the same squared
coefficient and the same inherited error. Hence R35's phase cone and
R36's information-current cone agree in this controlled limit.
Both exact phase bands remain in the complete observer and curvature
field. No independent physical-field universality is inferred.

One must not replace the full curvature continuation by reversed source
time: in general \(JZJ^\dagger\ne Z^\dagger\). The exact source
defect is nonzero, but its constant and first jets vanish. Quantitatively,

\[
\boxed{\|JZ(\epsilon k)J^\dagger-Z(\epsilon k)^\dagger\|
          \le176\epsilon^2K^2.}
\tag{36.20}
\]

**Proof.** R36.3 makes O,D mutual pairing-isometric inverses.
The three native actions in the proof of R36.5 give (36.18).
Conjugating each polynomial partial sum and its completed limit by these
maps transports R35.5's 89 n epsilon-squared K-squared bound unchanged,
proving (36.19). Conjugation by J works the same way, with
J A_i J-dagger=-A_i. It preserves the exact source characteristic
polynomial and band ranks as well.

For (36.20), R35's complete nine-coefficient moment gives
\(Z(\epsilon k)=I+\iota\epsilon A(k)+E\),
\(\|E\|\le88\epsilon^2K^2\). Apply J and dagger separately;
their linear terms agree and their two error norms add, proving the bound.
Nonzero exact Laurent coefficients in the source-built defect verify
that equality with Z-dagger is false. The certificate checks that defect
coefficientwise, as well as its vanishing zeroth and first jets. Thus
current reversal and the common leading cone do not erase the finite
retained-memory correction.

### R36.7 — Frame transport preserves the reading; the Riemann and physical identifications remain typed

Let G be a native pairing-preserving constant role frame. Transform
A_i,H_c,K_c,J,P and source values by the same frame. The response metric,
curvature-area law, complete observer, current bound and decoder are
unchanged as native readings. For a local frame on the count addresses,
transport Z and the shifts as well; the finite identities and locality
are preserved. The transported phase probes, rather than unchanged bare
coordinate formulas, express the same limiting statements.

The scalar component-address shifts commute and have identity elementary
loop. They coexist with the nonzero role protocol curvature (36.7) and
the constant response-dual metric (36.3). These are different mathematical
objects. Thus the nonzero protocol curvature does not by itself certify
a Riemann curvature tensor of the information-length form, and a flat
address reading cannot discard the retained minus-identity role loop.

The full native source has now supplied a response metric, its dual
length, a curvature-sensitive complete observer, compatible paired
dynamics and a bounded leading information current. A physical reading
still requires selection of this source/interface, physical field/source
and rod/clock identifications, and control of finite-front observables.
For a later adapter assigning length L per component step and duration T
per block, the count-current bound converts to \(L/(\sqrt2T)\).
No coefficient here selects L/T, physical energy/action, charge or alpha.

**Proof.** Pairing-preserving conjugation preserves every finite product,
sum and dagger. In particular P^G G u=G P u and
-P^G J^G G u=-G P J u; the reconstructed value is G u. The matching
forms, signed current pairings and error bounds are therefore preserved.
For local frames, adjacent inverse frames cancel in each transported
operator word, exactly as in R30.6. A pure transformed shift link between
addresses z,w is G_z G_w^dagger. Both orders along a scalar address cell
telescope to the same endpoint frame, so their loop stays identity.
In contrast the role protocol loop is conjugate to -I and remains -I.
This proves their separation by explicit source words.

The R10 smooth tangent descent and Levi-Civita conditions have not been
derived for this response target; they retain their existing admitted
adapter status. They are recorded as comparison, not used to prove any
new equation here. Finally source count equations contain neither L nor
T; substituting positive choices leaves them unchanged. The conversion
and the absence of a physical normalization follow. R35's exact faster
support witness also persists under the complete observer, so the new
leading-current law cannot upgrade its cone to a microscopic causal
front. These statements retain the exact earlier results while specifying
which part of the physical bridge remains to construct.

## Certification and pinned reuse

The [native application](../04-operator-evolution/native_curvature_observer.cjs)
constructs Z and its generators again from the pinned R17 role tags using
the unchanged canonical RKF engine. It checks the complete response Gram,
curvature products, observer/decoder, exact observed wave and memory,
current/deficit identities, native first-jet continuity, finite curvature
correction, frames and signed-readout counterexamples.

The [verification](../04-operator-evolution/R36_VERIFICATION.json) binds
seven written proofs, nine exact check groups, fourteen native word
replays, fifteen false alternatives and nine graph mutation controls.
It replays the frozen R35 source chain and preserves all 278 earlier
non-navigation files. The [ledger](../04-operator-evolution/R36_DERIVATION_LEDGER.json)
keeps the R1–R33 synthesis as inherited context and binds the specific
curvature and observer results reused here.

Written mathematical arguments, finite execution, provenance, formal
verification and physical identification remain separate. The graph
gate audits declared dependencies and hashes; it does not semantically
prove arbitrary prose. No whole-engine, Lean, external-review or physical
experiment PASS is claimed.

| Pinned prior result | Status and use in R36 |
|---|---|
| [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | **Native construction/derivation:** roles, matching, H/K, earned iota and completion. Finite equality, counting and induction remain explicit infrastructure. |
| [R20 native tuples](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | **Native construction:** retained role tuples; no physical tensor or measurement axiom supplied. |
| [R28.1–R28.5](https://github.com/Parveen117/extra-ideas/blob/9ed6ae27d3ddb005ed7f4f3f27417036729698fd/02-relational-response/NATIVE_CUT_CURVATURE_PROPAGATION_R28.md) | **Native-derived:** cut/curvature distinction, hidden-channel elimination and signed propagation. R36 applies that mechanism to the new response-derived cut. |
| [R30.3, R30.6–R30.8](https://github.com/Parveen117/extra-ideas/blob/7ac5b3e5baba4dafb51f4330bfe84d001a9ab1d4/02-relational-response/NATIVE_DIRECTIONAL_LOOP_INTERACTION_R30.md) | **Native-derived within stated targets:** response intertwiners, simultaneous frame transport and observer-kernel closure. The source target remains explicit. |
| [R32](https://github.com/Parveen117/extra-ideas/blob/1105a43209b144bd34a81bfbb99acc7876948e95/02-relational-response/NATIVE_RETAINED_LOOP_GAP_R32.md) | **Native-derived retained interface:** source events, nonflat memory links and central wave polynomial. Physical vacuum/field selection is not an input. |
| [R35.1–R35.7](https://github.com/Parveen117/extra-ideas/blob/9c13e1752263c756d356dcf5851ac873d2e0d102/02-relational-response/NATIVE_SIGNED_ENVELOPE_R35.md) | **Native-derived:** signed envelope, full-rank generators, native phase bands, sharp phase slopes, modewise error and exact-front boundary. |
| [RKF F00-E §§1–5](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md) | **Native analytic theorem on the earned field:** factorial flows, derivatives and circular identities. Written completion arguments and inherited finite audit are distinct evidence. |
| [R10 curvature-to-Riemann descent](https://github.com/Parveen117/extra-ideas/blob/6035e66bd40fe4016eeded41a4d704ec4ff4b342/02-relational-response/NATIVE_CURVATURE_DESCENT_R10.md) | **Comparison only; admitted smooth/tangent/metric/Levi-Civita adapter.** Its conditional quotient theorem is preserved. None of those inputs enters the new native proofs. |
| [Publications U30–U32 response tensor](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/NATIVE_RESPONSE_TENSOR_AND_CUT_LEDGER.md) | **Comparison only:** prior declared-carrier response geometry and separate admitted ensemble weights. R36 derives its required product, metric and curvature from R35 source coefficients. |
| [Canonical operator engine](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | **Unchanged native exact arithmetic/replay:** scoped application evidence; no second engine or whole-engine certification. |

```bash
python3.12 -B 04-operator-evolution/verify_r36.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```

The next physical task can now use an explicit observer and information
length. It must identify an operational source/field and a clock/rod
process, then decide which finite-front and long-scale observables it
reads. Curvature is already part of that observer construction; no
separate classical curvature assumption needs to be added for it.
