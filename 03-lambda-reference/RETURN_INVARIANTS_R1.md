# R1: Reference-independent return multipliers

Research program: Monty Dabas / Extra Ideas. Development date: 30 September 2026.

**Result.** In an explicitly defined finite scalar response network, products around closed walks are unchanged by independent local reference rescalings. A spanning tree gives a complete set of such invariants. Their constancy during evolution has a separate necessary and sufficient condition.

This constructs a restricted, information-preserving four-way translation. It does not determine the draft's universal reference `lambda_*`. Write the new quantities as `Lambda_C`, indexed by a specified closed walk `C`, to keep the two roles distinct.

## 1. Model and added assumptions

Let `G=(V,E)` be a finite connected labelled simple graph, with `n` vertices and `m` edges. Choose one orientation for each edge. Work over a specified scalar field `F`; the implementation uses exact rational numbers. At vertex `i`, attach a one-dimensional space `L_i`.

Assign each oriented edge `e:i->j` an invertible linear transport `L_i -> L_j`. In chosen local coordinates it has a nonzero coefficient `a_ij`. Reverse traversal uses `a_ji=a_ij^{-1}`.

These are assumptions of this development: finiteness, scalar invertible transport, reciprocal reversal, and a fixed labelled graph. They are not all established by the supplied manuscript. No metric, inner product, spacetime, or Hilbert completion is required.

An edge is a **process from one labelled state space to another**, not a demand that all equations `x_j=a_ij x_i` hold simultaneously for a single nonzero global state. The latter demand would force every closed return to be the identity.

Change the local coordinates by `x_i'=g_i x_i`, with `g_i` nonzero. Then

$$a_{ij}'=g_j a_{ij}g_i^{-1}.$$

Call two networks equivalent under a reference change when they are related this way, keeping the graph and labels fixed. This is the precise transformation group for every invariance claim below.

Independent vertex rescaling is an explicit choice for this model, stronger than a single shared-reference change. A common scaling of the form `g_i=c^{k_i}` is a subgroup, so the return invariants survive it too. If native laws restrict the allowed `g_i`, the completeness statement and the tree conclusion must be recomputed for that smaller group; the manuscript does not by itself establish unrestricted local reference freedom.

For a walk `P=(i_0,...,i_k)`, define

$$A(P)=a_{i_{k-1}i_k}\cdots a_{i_0i_1}.$$

The empty-length walk at a vertex has gain one. For a closed walk `C`, write `Lambda_C=A(C)`. If coefficients carry input/output units, those units cancel around the closed walk.

## 2. The return invariant

**Theorem 1.** A reference change transforms a path gain by

$$A'(P)=g_{i_k}A(P)g_{i_0}^{-1}.$$

Consequently `Lambda_C'=Lambda_C` for every closed walk. Also, for two paths `P,Q` with the same starting and ending vertices,

$$R(P,Q)=A(P)/A(Q)$$

is invariant and equals the gain of the walk that follows `P` and then reverses `Q`.

**Proof.** Substitute the transformation law into the ordered product. Every intermediate factor `g_i^{-1}g_i` cancels. Only the endpoint factors remain. For a closed walk, the endpoint factors cancel because the coefficients are scalars in a commutative field. Equal-endpoint path gains acquire the same multiplicative factor, so their ratio is unchanged. Reciprocal reversal proves the final assertion. QED.

The walk matters: reversing it replaces `Lambda_C` by `Lambda_C^{-1}`. Different cycles can have different values. An invariant of the complete model is the collection of return values, with its composition relations; it is not automatically one distinguished number.

## 3. Complete classification by a spanning tree

Fix a root `r` and a spanning tree `T`. Let `p_i` be the transport gain along the unique tree path from `r` to `i`, with `p_r=1`. For an oriented edge `e:i->j` outside the tree, define

$$h_e=p_i a_{ij}/p_j.$$

This is the return gain for the closed walk `r -> i` along the tree, then `i -> j`, then the reverse tree path `j -> r`.

**Theorem 2.** Set `g_i=p_i^{-1}`. All tree-edge gains become one, and each non-tree edge becomes `h_e`. For the fixed graph, tree, and orientations, two networks are equivalent under a reference change if and only if their `m-n+1` values `h_e` agree. Every assignment of nonzero values to these non-tree edges occurs.

**Proof.** Along a tree edge `i->j`, the rooted path gains obey `p_j=a_ij p_i`, including the case where the root path traverses the chosen orientation backwards. Thus the normalized coefficient `p_j^{-1}a_ij p_i` is one. On a non-tree edge it is `h_e`.

Theorem 1 proves that reference-equivalent networks have the same `h_e`. Conversely, equal `h_e` give identical tree-normalized edge data, so reversing one normalization and applying the other supplies the required reference change. For existence, put gain one on the tree and the prescribed gain on each remaining edge. A spanning tree has `n-1` edges, leaving `m-n+1` free nonzero gains. QED.

Equivalently, for this scalar model the reference quotient has coordinates

$$\{\text{edge transports}\}/\{\text{vertex reference changes}\}
\ \cong\ (F^\times)^{m-n+1}.$$

These coordinates depend on the chosen tree. The complete return assignment does not: each closed-walk gain can be reconstructed as a product of normalized non-tree gains and their inverses. Changing the tree changes the invariant coordinates, not the return values of specified walks.

**Consequences.** A tree has no nontrivial return parameter. A graph with one independent cycle has one such parameter. More cycles generally require more parameters; reducing them to a single lambda would need an additional structural relation or a distinguished cycle.

## 4. When an invariant is also constant during evolution

Assume now that `a_e(s)>0` are continuously differentiable functions on an interval. The graph and its edge orientations are fixed. Define the logarithmic rate

$$b_e(s)=\frac{d}{ds}\log a_e(s).$$

For a closed walk `C`, let `epsilon_e(C)` be its signed net traversal count of edge `e`.

**Theorem 3.**

$$\frac{d}{ds}\log\Lambda_C(s)=\sum_e\epsilon_e(C)b_e(s).$$

Thus a specified return is constant exactly when this sum vanishes throughout the interval. All closed returns are constant if and only if, at every `s`, there are vertex rates `eta_i(s)` such that

$$b_{ij}(s)=\eta_j(s)-\eta_i(s).$$

Equivalently, the evolution stays within a single reference-equivalence class:

$$a_{ij}(s)=g_j(s)\,a_{ij}(s_0)\,g_i(s)^{-1},\qquad g_i(s_0)=1.$$

**Proof.** Differentiate the logarithm of the product to obtain the first identity. Vertex differences telescope to zero around every walk. For the converse, define `eta_r=0` and let `eta_i` be the sum of the signed rates along the tree path `r->i`. This gives the stated difference on tree edges. Vanishing of the sum on each fundamental cycle gives it on every non-tree edge. Finally take `g_i(s)=exp(integral_{s_0}^s eta_i(u) du)` and integrate the logarithmic edge equation. QED.

Reference independence therefore does not by itself imply a conserved quantity. A genuine change in an independent cycle parameter can change `Lambda_C` even though its value remains independent of local units at each instant.

## 5. Ordinary common-parameter derivatives are flat

**Proposition 4.** Let `X_i(s)` be differentiable observables at the same parameter value, with all `dX_i/ds` nonzero. Set

$$a_{ij}=\frac{dX_j/ds}{dX_i/ds}.$$

Every closed return is one.

**Proof.** Put `v_i=dX_i/ds`. Then `a_ij=v_j/v_i`; the cycle product telescopes. Equivalently, normalize using `g_i=v_i^{-1}` to make every edge one. QED.

Consequently, a nontrivial scalar return cannot be obtained by relabelling ratios of a single common tangent vector. It requires genuinely path-dependent transport or other additional structure. This proposition does not assert that derivatives under different held-fixed constraints, or at different states, must give one. Those require their own typed model.

## 6. The four-way translation in this model

Let `W` be the labelled direct sum of the spaces `L_i`. In coordinates, `W=F^n`, with basis vector `e_i` and projector `P_i`. Write `E_ji` for the matrix unit that maps `e_i` to `e_j`, killing the other coordinate vectors. Define the edge operator

$$T_{ji}=a_{ij}E_{ji}.$$

**Proposition 5.** A composable path gives

$$T(P)=A(P)E_{i_k i_0},\qquad T(C)=\Lambda_C P_{i_0}.$$

With `D_g=diag(g_i)`, reference changes act by `T_ji'=D_g T_ji D_g^{-1}`. The labelled operators recover every edge gain, and tree normalization recovers the entire reference class from the fundamental return values.

**Proof.** Use `E_kj E_ji=E_ki` successively. Diagonal conjugation multiplies `E_ji` by `g_j/g_i`. Reading the unique nonzero entry of a labelled edge operator recovers its scalar gain. Theorem 2 gives the final reconstruction. QED.

| Way | Explicit realization | Data retained |
| --- | --- | --- |
| 1: intrinsic coordinates | `W=F^n` with labelled coordinate lines and projectors | State coordinates and source/target labels |
| 2: response relations | Nonzero edge coefficients and reciprocal path transport | The full response network |
| 3: reference normalization | Tree gains one; non-tree gains `h_e` | The full reference-equivalence class |
| 4: operator composition | Labelled `T_ji=a_ij E_ji` | All edge gains and every ordered path action |

The labels and operators are essential. A bare vector does not encode a network. Way-3 intentionally removes local reference choices; include the factors `p_i` if an exact original representative is required. This is a constructive correspondence for the stated model, not a proof of equivalence of all mathematical structures in `4ways.tex`.

## 7. Exact example and observable interpretation

Choose `a_AB=2`, `a_BC=3`, and `a_CA=1/7`. Then

$$\Lambda_{A\to B\to C\to A}=2\cdot3\cdot\frac17=\frac67.$$

Under `(g_A,g_B,g_C)=(5,11,13)`, the edge gains become `22/5`, `39/11`, and `5/91`; their product remains `6/7`. With tree `AB,BC`, the normalized gains are `1,1,6/7`. The closed operator is `(6/7)P_A`.

For a nonzero input on `L_A`, this model predicts an output/input ratio `6/7` after that particular process cycle. Comparing two processes with the same endpoints gives the same ratio without requiring a common calibration between those endpoints. This is a conditional prediction of the supplied example data, not an experimental result.

Changing the actual `BC` coefficient from `3` to `4` changes the return to `8/7`; such a change cannot be removed by a reference rescaling. The example values were chosen for transparent arithmetic and have not been derived from the framework.

## 8. What this says about the manuscript's six-stage chain

The response chain `OM -> NA -> MA -> SHI -> VA -> YA`, interpreted as five nonzero scalar transports, has `n=6`, `m=5`, and hence zero independent cycle parameters. Every edge coefficient can be normalized to one. The open chain alone therefore selects no nontrivial invariant scalar from its edge gains under the admitted local reference changes.

Adding a return edge would produce one independent cycle. Its multiplier would be the product of the six gains. But that return coefficient would be a new input unless the native rules determine it. Indeed, for any desired nonzero real `r`, choose the return gain as `r/(a_1 a_2 a_3 a_4 a_5)`; the cycle value is then `r`. This proves that closure topology alone does not select a universal numerical constant.

## 9. Progress and remaining selection problem

R1 provides an explicit reference-free candidate `Lambda_C`, a complete scalar classification, a conservation criterion, an operator realization, and exact reproducible checks. It replaces the unrestricted expression `x=n lambda_*` as the sole source of invariant information with a specified return process.

The next selection problem is concrete: derive the admissible return process and its coefficient from additional native laws, then show that all admissible models give the same distinguished invariant. If multiple independent cycles remain, prove the relations that reduce their free parameters. A scale-free invariant can be computed from a given network without being universal across networks.

The present reciprocal scalar assumption should also be tested against any intended cut or information-loss process. A noninvertible cut cannot be inserted as an invertible edge without changing the model.

## 10. Verification and mathematical provenance

Run from the repository root with Python 3.11 or 3.12:

```bash
python3.12 04-operator-evolution/verify_r1.py
```

The standard-library implementation and its focused rational checks are in [response_transport.py](../04-operator-evolution/response_transport.py) and [test_response_transport.py](../04-operator-evolution/test_response_transport.py). The generated [verification record](../04-operator-evolution/R1_VERIFICATION.json) records source hashes and the scope of the run. Exact finite checks support the examples; the general results rely on the proofs above.

The underlying construction belongs to established gain-graph mathematics: invertible edge labels, reversal, reference changes (switching), and cycle balance. See Thomas Zaslavsky, *Biased graphs. I. Bias, balance, and gains*, Journal of Combinatorial Theory, Series B 47 (1989), 32–52, [author-hosted paper](https://people.math.binghamton.edu/zaslav/Tpapers/bg1.jctb1989.pdf), and [the author's description of gain graphs](https://people.math.binghamton.edu/zaslav/Bsg/index.html). This development supplies explicit definitions, proofs, implementation, and a selection test for the four-way proposal; it does not claim discovery of gain-graph invariants.
