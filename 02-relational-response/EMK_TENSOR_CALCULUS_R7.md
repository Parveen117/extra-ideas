# R7 — A typed tensor foundation for EMK

Research owner: Monty Dabas. Development date: 30 September 2026.

This development fixes the tensor layer before choosing physical equations. Its starting point is the existing EMK KIR carrier, ordered recognition transport and lawful memory ledger. It supplies precise tensor types, transformation laws, finite derivatives, curvature and metric compatibility. It also distinguishes these constructions from proposals in older Vault manuscripts that require additional hypotheses.

The executable core uses exact rational arithmetic. The smooth extension below is conditional on separately admitted differentiable structure. Neither a spacetime, a Hilbert space, a universal lambda, nor a field equation is derived here.

## 1. What belongs to EMK, and what is established mathematics

| Layer | Classical tensor/bundle calculus | EMK realization in this development |
| --- | --- | --- |
| Tensor | Multilinear object on a specified carrier and its dual | Ordered, named E/M or recognition-carrier slots |
| Transport | Linear maps; dual maps are inverse transposes | Existing KIR maps lifted separately to each slot |
| Connection | Rule for comparing neighboring fibers | Ordered recognition transitions first; smooth RTC connection only in an admitted smooth sector |
| Curvature | Failure of two-direction transport to agree | KIR/channel order defect, retaining cross-commutators |
| Metric | Additional nondegenerate symmetric bilinear form | Must be selected or transported; the full two-mode R/K flow does not preserve a fixed symmetric metric |
| Return | Holonomy depends on the chosen representation | Audit carrier return, tensor-type return and sheet ledger separately |
| Recognition | May be modeled by observation/projection maps | Active recognition and lawful ledger determine the declared closure decision |

Classical geometry already includes rotations, holonomy, vector bundles and covering-space monodromy. Sheet rotation alone does not distinguish EMK from all classical mathematics. The framework-specific work is specifying its native operators, admitted state sectors, recognition targets and lawful memory policy, then proving which representations preserve those distinctions.

The KIR relations and mixed RTC curvature already exist in Vault and public EMK certificates. Tensor products, induced connections, exterior algebra and Bianchi identities are established mathematics. R7's contribution to this repository is their consistent typed implementation on the pinned native carrier, together with explicit rank-parity, metric and closure diagnostics. This is not a claim of historical priority for those general identities.

## 2. Admitted scalar and carrier layer

Fix a central field F of characteristic zero. Use F=Q for exact finite certificates and F=R for the optional smooth sector. At each description x, admit finite-dimensional carriers V_e(x), V_m(x), or other explicitly named V_f(x). A slot is `(fiber, variance, dimension)`; variance +1 means V_f and −1 means V_f*.

For an ordered slot list s, define

\[
\mathcal T_s(x)=\bigotimes_{a=1}^{r} W_a(x),\qquad
W_a=V_{f_a}\ \text{or}\ V_{f_a}^{*}.
\]

The rank-zero carrier is F. Addition requires identical slot lists. Tensor multiplication concatenates lists. Contraction evaluates a carrier's own dual against that carrier. Permutation changes slot order; symmetrization or alternation averages only identical slot types.

These admissions matter. A set of nonlinear states is not automatically a vector space. A partial UGD phase/scale/sheet product is not automatically a globally distributive scalar ring. KIR is an operator algebra over F, not the field of scalar coefficients. A morphic overlay operation denoted by a tensor symbol is not automatically this multilinear tensor product: identifying them requires a representation preserving their declared operations.

The source master object `(e,m,kappa; Rec,E,tau,Kernel; Phi,Omega,Lambda,varphi,Z; Smriti,Seam)` is a structured packet. Assign each component a carrier, type and transformation rule. Merely listing the components does not make the entire packet one multilinear tensor of a single rank. A gap functional kappa may select coefficients at a base description; it cannot depend arbitrarily on each input argument if multilinearity is asserted.

## 3. Frame changes and the E/M bridge

Write component frame changes as v'_f=G_f v_f. Evaluation invariance forces

\[
\omega'_f=G_f^{-T}\omega_f.
\]

Consequently a tensor transforms by G on every vector slot and G^{-T} on every dual slot. For example,

\[
T'{}^{i_1\cdots i_p}_{j_1\cdots j_q}
=G^{i_1}_{a_1}\cdots G^{i_p}_{a_p}
(G^{-1})^{b_1}_{j_1}\cdots(G^{-1})^{b_q}_{j_q}
T^{a_1\cdots a_p}_{b_1\cdots b_q}.
\]

Different named fibers may change frames independently. An E-covector cannot evaluate an M-vector without an admitted map B:V_m→V_e. Its scalar pairing is

\[
\omega_e^T Bv_m,\qquad B'=G_e B G_m^{-1}.
\]

**Proof of invariance.** Substitute the three transformed expressions; adjacent inverse matrices cancel. Holding E fixed is the special gauge choice G_e=I, not a universal rule that upper indices never transform. Without B, the cross-fiber contraction is undefined and the implementation rejects it.

For an endomorphism T∈V⊗V*, the law is T'=GTG^{-1}. For a metric H∈V*⊗V*, it is H'=G^{-T}HG^{-1}. They are different tensor types even when both are displayed as square matrices.

## 4. Native KIR action and tensor-rank parity

Use the unchanged public EMK two-mode realization

\[
I=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
K=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

It obeys R²=−I, K²=I and RK=−KR. Its linear span with RK is all real 2×2 matrices; this is inherited EMK algebra, not a new R7 construction. The old cut swap is represented by K in this realization. An older expression involving complex conjugation must be treated as real-linear when used here; it is not generally complex-linear.

For invertible U, let rho_s(U) apply U on vector slots and U^{-T} on dual slots. Multiple fibers use independent U_f.

**Theorem 1 — transport lift and parity.** For a single carrier and rank r,

\[
\rho_s(UV)=\rho_s(U)\rho_s(V),\qquad
\rho_s(-I)=(-1)^r I_{\mathcal T_s},
\]
\[
\rho_s(R)^2=(-1)^r I_{\mathcal T_s},\quad
\rho_s(K)^2=I_{\mathcal T_s},\quad
\rho_s(R)\rho_s(K)=(-1)^r\rho_s(K)\rho_s(R).
\]

**Proof.** The vector and dual assignments are group representations: `(UV)^(-T)=U^(-T)V^(-T)`. Their tensor products preserve composition. Every slot receives −I under central transport −I, including dual slots. Apply the representation to each native relation. This proves the identities for every rank and variance pattern, not just sampled components.

Thus all even-rank types are blind to I versus −I. A rank-two metric and a rank-two endomorphism both return under −I while vectors change sign. The two ordered native paths with transports KR and RK differ by central −I; their even-rank lifts agree. A metric-only closure check loses that distinction. This is representation blindness, not evidence that classical tensor calculus lacks rotation.

The finite tensor lift is generally **not linear in U**: rho(U+V) need not equal rho(U)+rho(V). Therefore it is a group representation, not a linear homomorphism of the entire associative KIR algebra. Requiring primitive R²=−I on every tensor rank would contradict Theorem 1.

## 5. Infinitesimal lift and the actual derivative rule

For an endomorphism A_f on each carrier define

\[
d\rho_s(A)T
=\sum_{\text{vector slots }a} A_{f_a}^{(a)}T
-\sum_{\text{dual slots }a}(A_{f_a}^{T})^{(a)}T.
\]

The superscript means action on that slot alone. This action is linear in A.

**Theorem 2 — induced derivative.** The slot action is a derivation of tensor multiplication, commutes with contraction and identical-type symmetrization, and obeys

\[
[d\rho_s(A),d\rho_s(B)]=d\rho_s([A,B]).
\]

**Proof.** Actions on different slots commute. On a vector slot the commutator is [A,B]. On a dual slot it is `[−A^T,−B^T]=−[A,B]^T`. In a contraction, A on the vector cancels −A^T on its paired dual. Dividing slots into those belonging to the first and second factors proves the product rule; identical slot permutations commute with the sum.

The lift preserves the Lie bracket, not arbitrary associative products or the primitive square relation. For type (1,1), d rho(A)T=[A,T]. For a vector it is Av; for a metric it is −A^T H−HA. A universal formula `[A,T]` for every tensor field is therefore incorrect.

No clock is needed for this algebraic operation. Turning it into a rate requires an admitted parameter and a limit that exists. If a normalization clock has zero increment, division is undefined; a zero clock must not silently set every derivative to zero.

## 6. Seam decomposition on all slots

For an admitted involution K on a carrier, let P_±=(I±K)/2. These are complementary idempotents. On dual slots use P_±^T. For a sign list eta∈{+1,−1}^r, define the corresponding tensor block by applying its projector to each slot.

**Theorem 3 — typed block resolution.** Every rank-r tensor is the sum of its 2^r seam blocks. Under frame change, transform K to GKG^{-1}; the entire block resolution transforms covariantly.

**Proof.** Expand the product of `(P_++P_-)=I` on every slot. The projectors are complementary on each slot. Conjugation transports vector projectors, and inverse-transpose frame change transports their dual transposes.

For rank two, `++`, `--`, `+-`, `-+` are explicit parallel/transverse/mixed sector labels after declaring this projector convention. They are not functions supported only on the literal diagonal of an unrelated manifold product. The latter interpretation would force a smooth function supported only on a proper diagonal to vanish.

## 7. Metric obstruction and the alternating invariant

The operator K has both a finite involution interpretation and a continuous generator interpretation. These must not be confused.

**Theorem 4 — no shared fixed symmetric metric on the standard carrier.** On the real two-mode carrier, no nonzero symmetric H is preserved by both continuous flows exp(tR) and exp(sK) for all real t,s.

**Proof.** Differentiating preservation at zero gives

\[
R^T H+HR=0,\qquad K^T H+HK=0.
\]

Write H=[[x,y],[y,z]]. The R condition forces y=0 and z=x, so H=xI. The K condition then gives 2xK=0, hence x=0. The argument excludes indefinite forms as well as positive ones. Conversely these infinitesimal equations suffice for each constant-generator flow, by differentiating exp(tA)^T H exp(tA).

This theorem is restricted to this representation and these continuous flows. Finite R and K themselves both preserve H=I; exp(sK) for nonzero real s does not. It does not prohibit metrics on larger representations, restricted transport families or a moving metric field. It also does not invalidate R3's metric for its separately restricted circular generator.

Set epsilon=[[0,1],[-1,0]]=−R. Direct multiplication proves, for every real 2×2 matrix U and every A,

\[
U^T\epsilon U=\det(U)\epsilon,\qquad
A^T\epsilon+\epsilon A=\operatorname{tr}(A)\epsilon.
\]

Therefore traceless R, K and RK generate transport preserving a nondegenerate **alternating area form**. Scalar dilation generated by I gives a conformal weight instead. An alternating form is not a symmetric distance metric and supplies no positive length. The natural shared invariant in this two-mode sector is area; selecting a physical metric requires further structure.

For a chosen nondegenerate symmetric H, raising and lowering use H^{-1} and H. They commute with frame change when H is transformed with its own covariant law. Exterior products need alternation but no metric. Hodge star, norm, divergence and metric Laplacian need additional metric and, where applicable, orientation and signature conventions.

## 8. Moving metrics and global compatibility

For an arrow e:x→y with carrier transport U_e define its metric defect

\[
Q_e=U_e^T H_y U_e-H_x.
\]

Q_e=0 is exactly compatibility. One arrow can always transport an admitted metric by `H_y=U_e^(-T)H_x U_e^(-1)`. Under reframe, Q'_e=G_x^{-T}Q_e G_x^{-1}; the zero decision is invariant.

**Theorem 5 — path consistency of a metric.** In a connected finite graph with invertible reciprocal transports, a chosen root metric H_0 extends compatibly to every edge if and only if every based loop holonomy W satisfies W^T H_0 W=H_0.

**Proof.** Necessity follows by composing metric compatibility around a loop. For sufficiency, choose a root-to-x path with transport U_x and set H_x=U_x^{-T}H_0U_x^{-1}. Two choices U_x and V_x differ by the root loop V_x^{-1}U_x. Its preservation of H_0 makes the two transported metrics equal. Appending an edge gives compatibility. A spanning tree and its chord-loop generators suffice to check all loops because preservation is closed under products and inverses.

This is the existing holonomy criterion applied to the admitted EMK transport. It separates choosing a metric from asserting that it glues globally. It also shows why metric consistency alone does not prove full carrier or sheet return: central −I preserves every symmetric H.

## 9. Clock-free tensor increments and lawful ledgers

Let tensor fields at x and y have the same named slot type. For a typed tensor ledger L_e in the target fiber use the backward comparison convention

\[
\delta_e T=\rho_s(U_e)T_x-T_y-L_e.
\]

It transforms in the target tensor type. This is a finite transported field difference, not an automatically derived smooth derivative. The general RSC naturality defect compares actual recognition with recognized transport; it has additional Rec maps. In a path-field specialization its forward convention has the opposite sign to this backward comparison, with the corresponding ledger sign changed. R7 does not redefine that upstream primitive.

For zero ledgers, tensor multiplication obeys the exact twisted product identity

\[
\delta_e(T\otimes S)
=(\delta_eT)\otimes\rho(U_e)S_x+T_y\otimes(\delta_eS).
\]

With declared ledgers L_T, L_S and L_TS, add the coupling term

\[
L_T\otimes\rho(U_e)S_x+T_y\otimes L_S-L_{TS}.
\]

**Proof.** Add and subtract `T_y tensor rho(U)S_x`, then substitute each ledger-subtracted increment. This works for arbitrary tensor types whose concatenation is admitted.

The zero-ledger difference is linear; a fixed nonzero subtractive ledger makes the comparison affine. More general state-dependent ledgers require their own declared linearity policy. A recognition map need not preserve tensor products. A Smriti correction need not be a derivation. Additional product defects must be retained unless the recognition and declared-output maps are proved multiplicative. The current RSC algebra-calculus section already states this qualification; R7 instantiates it with explicit tensor types.

Ledgers are fixed by the model or protocol independently of the residual being audited. Choosing a ledger equal to a failed residual afterward proves no predictive closure.

## 10. Ordered integration and the separate sheet register

For successive arrows e_j, write E_j=delta_ej T and B_j=L_j+E_j. The update is

\[
T_j=\rho(U_j)T_{j-1}-B_j.
\]

**Theorem 6 — finite covariant integration.** With W=U_n⋯U_1,

\[
T_n=\rho(W)T_0-
\sum_{j=1}^{n}\rho(U_n\cdots U_{j+1})B_j.
\]

**Proof.** Substitute the update successively, using Theorem 1. The empty downstream product is identity. This is an exact finite fundamental theorem along a declared path, with no approximation or external clock.

In the implemented **pure winding sector**, a move also carries an integer n_e. Moves compose as `(U_b,n_b)(U_a,n_a)=(U_bU_a,n_b+n_a)` and reverse as `(U^(-1),−n)`. Under integer sheet-reference changes b_x,

\[
n'_e=n_e+b_y-b_x.
\]

The path sums telescope, so closed-loop winding is reference invariant. This direct-product sector is an explicit admission, not a derivation of every UGD sheet module; rational/scale-coupled memory sectors need different rules.

There are three distinct questions:

1. Does a tested tensor state or an entire tensor type return?
2. Does the full native carrier return?
3. Does the active sheet register return, or match its preregistered lawful ledger?

A loop with U=−I returns on every even-rank tensor type but fails carrier return. A loop with U=I and n=2 returns on every tensor type but fails unledgered sheet return. The checker audits these separately. It does not assert that n is automatically an integral of matrix curvature: such an identification requires a separately specified readout, normalization and quantization theorem.

## 11. Finite exterior calculus and curvature

Admit an ordered simplicial complex with a carrier at each vertex and invertible edge transport U_ij:V_i→V_j. Tensor-valued edges use rho_s(U_ij). For a vertex section s and an edge one-cochain alpha stored in its terminal fiber, define

\[
(D_0s)_{ij}=s_j-U_{ij}s_i,\qquad
(D_1\alpha)_{012}=\alpha_{12}-\alpha_{02}+U_{12}\alpha_{01}.
\]

For a triangle set

\[
F_{012}=U_{02}-U_{12}U_{01}:V_0\to V_2.
\]

**Theorem 7 — curved square and finite Bianchi.**

\[
(D_1D_0s)_{012}=F_{012}s_0,
\]
\[
F_{123}U_{01}-F_{023}+F_{013}-U_{23}F_{012}=0.
\]

**Proof.** Substitute D_0 in D_1; the s_1 and s_2 terms cancel. For Bianchi, expand every F. The four direct-transport terms cancel, and associativity cancels the two copies of `U_23 U_12 U_01`. Flatness is not needed for the Bianchi identity.

Under vertex frame changes F'_012=G_2F_012G_0^{-1}; its zero status is invariant. For tensor coefficients first lift the edges, then subtract:

\[
F^{s}_{012}=\rho_s(U_{02})-
\rho_s(U_{12})\rho_s(U_{01}).
\]

It is generally incorrect to write rho_s(F_012), because the finite lift is nonlinear and F_012 may be singular. For a curved connection D² is not zero. Consequently `kernel(D)/image(D)` is not automatically cohomology: one must first establish square-zero, or specify an invariant restricted module on which the curvature acts trivially. The full higher-degree calculus requires the corresponding typed cochain definitions; R7 implements the vertex/edge/triangle layer and its tetrahedron identity.

This is established discrete bundle calculus specialized to EMK. See the primary discrete-vector-bundle reference below; nonzero Bianchi is not a new universal EMK law.

## 12. Optional smooth connection and curvature

Now separately admit a smooth base B, smooth carrier bundles and sufficiently differentiable coefficient fields. Use the convention

\[
\nabla_X T=\partial_XT+d\rho_s(A(X))T.
\]

For component frame changes v'=Gv the connection law is

\[
A'=GAG^{-1}-(dG)G^{-1}.
\]

It makes `(nabla T)'=rho(G)nabla T`. The tensor derivative has the product and contraction rules of Theorem 2. If short-edge parallel transport is `U_h=I−h A(X)+o(h)`, the zero-ledger backward increment divided by h tends to −nabla_X T. Different forward conventions change that sign; it must be declared.

For vector fields X,Y,

\[
F(X,Y)=[\nabla_X,\nabla_Y]-\nabla_{[X,Y]},
\]
\[
F_{ij}=\partial_iA_j-\partial_jA_i+[A_i,A_j]-A([X_i,X_j]).
\]

Omitting the last term is valid in a commuting coordinate basis, not for arbitrary noncoordinate directions. The induced tensor curvature is d rho_s(F), by Theorem 2. In the implementation, derivative coefficients are supplied exact inputs; omitting them declares constant coefficients rather than estimating derivatives.

If the existing RTC channels all act on the same admitted carrier, set A=sum_a A^a. Then

\[
F=\sum_a(dA^a+A^a\wedge A^a)
+\sum_{a<b}(A^a\wedge A^b+A^b\wedge A^a).
\]

Mixed cross-commutators cannot be dropped. A sum of maps with mismatched domains is not a connection until common-carrier embeddings have been supplied. This cross-channel requirement is already present in EMKT1/EMKT2.

For a genuine associative smooth connection,

\[
d_A F=dF+A\wedge F-F\wedge A=0.
\]

**Proof.** Expand F=dA+A∧A, use d²=0 and the graded product rule; all terms cancel. Sufficient differentiability is required for these operations. A reported nonzero Bianchi residual therefore concerns a projected/truncated quantity, omitted channels, a different derivative or an invalid connection declaration. It does not overturn the identity for the full connection.

The smooth metric defect is

\[
Q(X)=\partial_XH-A(X)^T H-H A(X).
\]

Compatibility is Q=0. Theorem 4 is its constant-H specialization. A scalar density or an integral of raw matrix entries is not automatically a topological charge. Nonabelian Stokes formulas require ordered transport and a defined readout; matrix flux, Wilson holonomy and integer winding are different typed objects.

## 13. Why Christoffel, Riemann, torsion and Ricci need more data

An End(V)-valued curvature two-form is defined without identifying V with tangent directions. To define torsion, admit a direction bundle D, its anchor/bracket and a soldering map S:D→V. Then

\[
\operatorname{Tor}(X,Y)=\nabla_X(SY)-\nabla_Y(SX)-S[X,Y].
\]

Without such S, a “native torsion” expression silently pairs different types. Ricci contraction likewise requires identifying the appropriate direction and carrier indices. A general operator-valued two-form cannot be contracted as a four-index Riemann tensor without that identification.

If B is a smooth manifold, V=TB, S is identity, g is nondegenerate, and the connection is metric-compatible and torsion-free, one recovers the Levi-Civita connection. In coordinate components,

\[
\Gamma^a_{bc}=\tfrac12 g^{ad}
(\partial_b g_{dc}+\partial_c g_{db}-\partial_d g_{bc}).
\]

Its curvature is the usual Riemann tensor, with a declared sign convention. Thus Christoffel and Riemann are valid conditional realizations, not primitives that must be renamed or universally rejected.

The current Vault warped metric `g_A=A(v)^2 du^2+dv^2` already chooses this smooth metric sector. It has Gamma^u_uv=A'/A, Gamma^v_uu=−AA' and Gaussian curvature −A''/A. Those are classical Levi-Civita calculations for a chosen EMK realization. They do not prove that arbitrary recognition curvature equals Gaussian curvature, or that every gap kappa is universally curvature. R7 credits the existing metric work rather than reproducing it as a new theorem.

## 14. Exact implementation and reproduction

Files:

- [emk_tensor_calculus.py](../04-operator-evolution/emk_tensor_calculus.py): ordered slots, finite and infinitesimal lifts, contraction, symmetry, wedge, index conversion, increments, integration, curvature, metric defects and return audits.
- [test_emk_tensor_calculus.py](../04-operator-evolution/test_emk_tensor_calculus.py): 16 focused tests, including all 30 rank/variance signatures at ranks 1–4, mixed-frame pairing, seam blocks, product-ledger defects, curved D², finite Bianchi, projected false closure and the metric obstruction.
- [R7_SOURCE_PINS.json](../04-operator-evolution/R7_SOURCE_PINS.json): exact upstream snapshots, file hashes and source-reading depth.
- [R7_VERIFICATION.json](../04-operator-evolution/R7_VERIFICATION.json): actual test result, concrete witnesses and preserved R1–R6 evidence.
- [EMK_VAULT_SOURCE_AUDIT_R7.md](EMK_VAULT_SOURCE_AUDIT_R7.md): source coverage, inherited work and necessary corrections.

Use Python 3.11 or 3.12 and a separate checkout of Publications at the commit recorded in the pins:

```bash
python3.12 -B 04-operator-evolution/verify_r7.py --publications-root ../Publications
```

The verifier checks unchanged public EMK2/EMKT1/EMKT2 source bytes before importing them. It reuses R2's exact matrix utilities without editing them, and checks every hash in the R1–R6 verification records. No private Vault source is copied into this public repository. Tests certify the finite examples; the general statements above have separate human-readable proofs.

## 15. Foundation status before physics

| Item | R7 status | Remaining admission or development |
| --- | --- | --- |
| Typed algebra and tensor transport | Implemented with general proofs and exact checks | Nonlinear states require a separate representation theorem |
| Seam blocks and rank-dependent return | Implemented | Which blocks each target must observe is a model choice |
| Finite difference, ledger product rule and path integration | Implemented | Nonlinear recognition requires its own typed defect maps |
| Triangle curvature, curved D² and tetrahedron Bianchi | Implemented | Full higher-degree tensor-valued cochains |
| Fixed metric obstruction and moving-metric compatibility | Proved; local defects and invariant forms implemented | Select a native metric policy and audit global holonomy |
| Smooth connection, curvature and Bianchi | Conditional algebraic derivation | Smooth realization, regularity and analytic existence |
| Torsion, Ricci, Hodge, divergence and Laplacian | Typed prerequisites specified | Direction carrier, soldering, metric/orientation and analytic domains |
| Physical equations or experiments | Not supplied | Only after the preceding choices and consistency checks |

The next foundation step is to choose and prove the native direction/soldering and metric policy, then build torsion/Ricci and higher exterior operators on that declared sector. An action, physical units, constitutive laws and experiment cannot be inferred merely from the existence of tensor notation.

## 16. Primary classical references and native lineage

- Dominic Joyce, [Introduction to Differential Geometry](https://people.maths.ox.ac.uk/joyce/Nairobi2019/IntroDiffGeom.html), Oxford course notes (2019): tensors, connections, curvature, torsion, Levi-Civita and bundles.
- Richard Melrose, [Differential Geometry](https://math.mit.edu/~rbm/18.157-F05.pdf), MIT notes: induced tensor/dual connections and Bianchi.
- Daniel Berwick-Evans, Anil Hirani and Mark Schubel, [Discrete Vector Bundles with Connection](https://arxiv.org/abs/2104.10277), consulted v3 revised 23 April 2026: discrete bundle connections and algebraic differential identities.
- Publications [EMK–UGD algebra](https://github.com/Parveen117/Publications/tree/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra): existing native KIR and RTC certificates; exact file pins are in the source index.
- Vault current math papers and recovered EMK/UGD/RSC foundations: private source paths and immutable commits are listed in the audit and pins. Their declarations are evaluated under the explicit assumptions above.
