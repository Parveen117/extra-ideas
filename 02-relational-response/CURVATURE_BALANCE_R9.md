# R9: Native curvature balance and the information ledger

The research question is conditional: when does recognition balance information-bearing, order-sensitive transport, and what remains after that balance? Observation is not assumed to flatten every curvature. Off-diagonal matrix entries are not the definition of curvature.

R9 supplies one explicit native mechanism. The existing derived cut-swap `K` acts on two related branches. Requiring the recognition output to be independent of their exchange uniquely selects equal branch weights within the declared mixing family. The resulting observer removes cut-odd curvature, retains cut-even curvature, and has an exact information/branch ledger. The elementary averaging operation and grading are existing mathematics; the linked EMK selection, dynamics and source checks are this research adapter's contribution.

## 1. Recognition curvature and Riemann curvature are different typed objects

Classical Riemann curvature acts on tangent vectors through a tangent-bundle connection, commonly the metric's Levi–Civita connection. EMK recognition curvature acts on a declared information/seam carrier and tracks the mismatch of lawfully comparable transports. Master-channel curvature, recognized residue and the closure ledger are additional declared targets.

| Quantity | Carrier and required data | Certificate here |
|---|---|---|
| Riemann curvature | Tangent bundle and a declared connection; a metric selects Levi–Civita in that special case | No identification with the recognition carrier is imposed |
| Native operator curvature | Common information/seam carrier, declared actions and lawful composition | Exact constant-coefficient commutators and their cut grades |
| Recognized curvature | A declared observer acting on the operator or state | Exact branch balance and retained even sector |
| Closure residue | Active sectors and independently declared compensation/ledger | Existing master policy is retained; a zero scalar mean does not replace it |

The commutator form is shared general connection mathematics. The carrier, transport law and recognition policy determine which tensor it describes. For example one may use a flat Euclidean base with zero Levi–Civita curvature and an independent information-fibre connection with constant coefficients `K,R`, whose curvature is `[K,R] != 0`. Equating the two curvatures would require an explicit bundle/connection adapter.

The distinction is therefore in the mathematical objects and framework contract; the general formula `[nabla_i,nabla_j]` is not claimed as a new discovery. No metric or spacetime is primitive in the algebraic results below.

## 2. Native mechanism and its selection condition

Use the existing [EMK-T1 cut-swap](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/emkt1_master_tensor_and_time.py): it is derived as `K`, with `K^2=I`. Define

\[
\alpha(A)=KAK,\qquad A_e=(A+\alpha A)/2,\qquad A_o=(A-\alpha A)/2.
\]

The cut scalars and lawful composition precede this chart. These are the existing native even/odd grades. R9 declares a convex two-branch recognition protocol

\[
\Phi_\theta(A)=(1-\theta)A+\theta\alpha(A)
=A_e+(1-2\theta)A_o,\qquad 0\leq\theta\leq1.
\]

For generic algebra results only the involution is required. For density-matrix/channel statements we use the admitted Hermitian/unitary two-mode `K`; a non-unitary involution alone would not guarantee positivity.

**Theorem R9.1 — unique exchange-balanced endpoint.** Suppose the cut action is nontrivial, and require the output to agree for the two exchanged inputs:

\[
\Phi_\theta\alpha=\Phi_\theta.
\]

Within this declared family the unique solution is `theta=1/2`, giving `Phi=A_e`. It satisfies `Phi^2=Phi`, `alpha Phi=Phi`, and `Phi(A_o)=0`.

**Proof.** The difference `Phi_theta alpha - Phi_theta` is `(1-2 theta)(alpha-id)`. Since `alpha-id` is nonzero and the central field has characteristic zero, exchange invariance forces `1-2 theta=0`. The other identities follow from `alpha^2=id`. □

Exchange invariance is the selection hypothesis. The half weight follows from it; a universal physical recognition law, a rate, and a clock do not follow from `K^2=I` alone. This result selects the balanced endpoint rather than post-selecting a value because one readout happened to vanish.

## 3. Theorem R9.2: balance selects curvature sectors

For `F=[A,B]`, cut parity and bilinearity give

\[
\boxed{\Phi(F)=[A_e,B_e]+[A_o,B_o].}
\]

The mixed even/odd terms are cut-odd and cancel. The odd/odd term is cut-even and remains. Consequently

\[
\boxed{\Phi([A,B])-[\Phi(A),\Phi(B)]=[A_o,B_o].}
\]

**Proof.** Expand both factors into their two grades. The cut action is an algebra automorphism, so even/even and odd/odd products are even, while the two mixed products are odd. Applying `Phi` retains exactly the two even commutators. □

This is the antisymmetric form of the product-memory identity

\[
\Phi(AB)-\Phi(A)\Phi(B)=A_oB_o.
\]

Thus recomputing curvature after separately averaging its generators can erase a retained even curvature. `Phi` is a projection onto the even algebra, but is generally not multiplicative. R8's compression identity is a different adapter; neither permits a universal observation/flatness claim.

On the existing chart

\[
K=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
RK=\begin{pmatrix}-1&0\\0&1\end{pmatrix},
\]

take an even radial action `K` and two odd tangential actions `R,RK`:

| Pair | Raw curvature | Cut parity | Recognized operator curvature |
|---|---|---|---|
| `K,R` | `-2 RK` | Odd | `0` |
| `K,RK` | `-2 R` | Odd | `0` |
| `R,RK` | `-2 K` | Even | `-2 K` |

These radial/tangential labels specify this adapter, not every physical direction. The table makes the conditional flatness precise: this balanced observer removes the mixed-parity sectors while retaining the tangential odd/odd memory. For a differentiable connection and a constant cut, the same commutator correction applies after linearly transforming the derivative terms. A varying cut requires its derivative contributions.

The unchanged EMK-T1 channel engine also reproduces this result with channels `(K,0)` and `(0,R)`: individual channel curvatures are zero, the mixed term is `-2 RK`, and its balanced readout vanishes. Replacing the pair by `(R,0),(0,RK)` retains the even mixed curvature. The raw master and mixed term remain in the audit.

## 4. Theorem R9.3: clock-free balance dynamics

For an ordered sequence of recognition events with declared weights `theta_j`, set `lambda_j=1-2 theta_j`. Then

\[
\boxed{A_n=A_e+\Lambda_nA_o,\qquad
\Lambda_n=\prod_{j=1}^n\lambda_j.}
\]

**Proof.** Each event is identity on the even subspace and scalar multiplication by `lambda_j` on the odd subspace. Compose those actions. □

If initial odd content is nonzero, exact finite balance occurs precisely when an event has weight `1/2`. An event above `1/2` reverses the odd response sign. A fixed weight strictly between zero and one contracts odd content asymptotically, except that `1/2` balances in one event. Boundary weights zero and one respectively retain or repeatedly flip the imbalance. No physical time variable is needed.

Example: initial axis imbalance `m=3/5`, followed by weights `1/4,3/4,1/2`, has odd multipliers `1/2,-1/4,0` and imbalances `3/10,-3/20,0`. Six events of weight `1/4` give multiplier `1/64`; that is an exact finite tolerance bound, not exact balance.

Trace cyclicity gives the dual state/readout law

\[
\operatorname{Tr}(\Phi_\theta(\rho)F)
=\operatorname{Tr}(\rho\Phi_\theta(F)).
\]

For an odd fixed probe its mean therefore scales by `lambda`, including sign reversal. A state is fully cut-balanced only when its entire odd part vanishes; a single zero mean need not establish that. Conversely a cut-balanced state can read even curvature: `rho=(I+K/2)/2` is balanced, yet `Tr(rho [R,RK])=-1`.

## 5. Theorem R9.4: mean-curvature and fluctuation budget

In the admitted positive-state chart let a Hermitian probe `H` be cut-odd. Since `H^2` is even, recognition preserves its second moment while scaling its mean:

\[
\mu'=\lambda\mu,\qquad s'=s,\qquad
v'=v+(1-\lambda^2)\mu^2.
\]

Here `mu=Tr(rho H)`, `s=Tr(rho H^2)`, and `v=s-mu^2`. The identity follows from the trace-dual law and `alpha(H^2)=H^2`; positivity makes the variance nonnegative.

For `H=[K,R]/2=-RK` and `rho_m=diag((1+m)/2,(1-m)/2)`, this is

\[
\boxed{\mu=m,\qquad s=1,\qquad \mu^2+v=1.}
\]

At the balance point the signed mean is zero and the variance is one. The operator's noncommuting origin and second moment remain present. This is a mean-versus-fluctuation curvature budget for the declared probe; it is not an identification of variance with every notion of Fisher information.

## 6. Theorem R9.5: transform the whole information family

The existing [QTH-1 adapter](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/qth1_quantum_recognition_information.py) has `sigma_x=K`, `sigma_y=iota R`, `sigma_z=-RK` in its admitted complex chart. No primitive imaginary unit is added to the native engine. Consider the full-rank local family

\[
\rho_m(x,y)=\tfrac12(I+x\sigma_x+y\sigma_y+m\sigma_z),
\qquad |m|<1,
\]

evaluated at `x=y=0`. Its SLDs are `L_x=sigma_x`, `L_y=sigma_y`. The source tensor has symmetric information `g=I_2` and antisymmetric pairing

\[
B_{xy}=\frac1{2\iota}\operatorname{Tr}(\rho[L_x,L_y])=m.
\]

This recovers the earlier static witness: different `m` values can reverse or balance the mean while the current probe pair remains noncommuting. Applying a recognition channel to the same family is a further operation. Both state and derivatives transform:

\[
(x,y,m)\longmapsto(x,\lambda y,\lambda m),\qquad
\partial_x\rho\mapsto\sigma_x/2,\quad
\partial_y\rho\mapsto\lambda\sigma_y/2.
\]

The recomputed SLDs are `L'_x=sigma_x`, `L'_y=lambda sigma_y`, so

\[
\boxed{g'=\operatorname{diag}(1,\lambda^2),\qquad
B'_{xy}=\lambda^2m.}
\]

**Proof.** Substitute these matrices into the source SLD equation. Their anticommutators with the output state give exactly the transformed derivatives. Then the symmetric products give the displayed metric, and the commutator contributes one `lambda` from `L'_y` and one from the output state imbalance. □

This yields the linked identity `B'_{xy}=m g'_{yy}` for this adapter. A fixed normalized `sigma_z` probe instead has mean `lambda m`. Their different transformation laws are certified explicitly. Branch exchange at `theta=1` flips that fixed mean while the co-transformed information tensor keeps its orientation pairing. Freezing derivatives would incorrectly report preserved information during averaging.

### Branch ledger and exact recovery

With the branch label `b` retained, each realized output is `K^b rho K^b`, and applying `K^b` again recovers the input. State-independent branch weights make the flagged information the weighted sum of the branch information metrics; each branch is a unitary cut action, so this is `I_2`.

If the label is discarded, the visible channel has `g'=diag(1,lambda^2)`. Thus

\[
\boxed{g_{\rm flagged}=g_{\rm visible}+g_{\rm discarded},\qquad
g_{\rm discarded}=\operatorname{diag}(0,1-\lambda^2).}
\]

Equal branch weights balance the output but discard one tangent information direction from its visible state. Keeping the actual branch record permits exact inverse recovery. At that endpoint distinct inputs `rho_(3/5)` and `rho_(-3/5)` have the same unflagged output, so the raw state/history cannot be reconstructed from that average alone. This consumes the existing master ledger idea with independently specified branch data; it does not invent a compensating ledger after seeing the desired result.

## 7. Theorem R9.6: zero mean curvature and joint recognition

A further control distinguishes mean-curvature balance from complete simultaneous parameter recognition. For any finite one-copy qubit effects

\[
E_a=w_a(I+\mathbf n_a\cdot\boldsymbol\sigma),\quad
w_a\geq0,\quad |\mathbf n_a|\leq1,\quad
\sum_aw_a=1,\quad\sum_aw_a\mathbf n_a=0,
\]

the `x,y` classical information at `(0,0,m)` is

\[
C=\sum_a\frac{w_a}{1+mn_{az}}
\begin{pmatrix}n_{ax}^2&n_{ax}n_{ay}\\n_{ax}n_{ay}&n_{ay}^2\end{pmatrix},
\qquad \boxed{\operatorname{tr}C\leq1.}
\]

**Proof.** Outcome probabilities are `p_a=w_a(1+m n_az)` and their tangent derivatives are `w_a n_ax,w_a n_ay`, giving the matrix formula. Since `|m|<1`, the denominators are positive. Also

\[
\frac{n_{ax}^2+n_{ay}^2}{1+mn_{az}}
\leq\frac{1-n_{az}^2}{1+mn_{az}}
\leq1-mn_{az}.
\]

The final inequality uses `1-n_z^2 <= 1-m^2 n_z^2`. Multiply by the weights, sum, and use completeness. □

At `m=0`, mean antisymmetric curvature is zero, while the original information metric is still `I_2`, of trace two. Thus this finite one-copy observation cannot retain both tangent information directions completely. An `x` projective readout has `C=diag(1,0)`; equal random choice between `x` and `y`, with its choice/outcome recorded, gives `C=I_2/2`. Both have an explicit positive information-loss ledger and trace gap one. The statement is for the admitted measurement model; no multi-copy attainability claim is made.

The bound is established quantum-estimation mathematics specialized here. R9 adds the exact comparison to the native balance mechanism and its residual ledger, rather than claiming a new universal uncertainty law.

## 8. Evidence, lineage and reproduction

The single active operator engine remains [RKF/operator_foundation](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation). R9 calls it unchanged. Its packet replays six identities in a free associative presentation with only `K^2=I`, and three KIR certificates: nonzero raw curvature, cancelled radial/tangential curvature, and surviving tangential/tangential curvature. No finite carrier is required for the six generic rewrite proofs. The grading, conditional-expectation identities and information methods are credited as existing; the EMK source integration and linked balance ledger are the extension.

EMK-T1 supplies the derived cut and mixed-channel master. QTH-1 supplies SLD/tensor and directional measurement operations. CID-1 supplies the prior recognized/discarded covariance interpretation. R8 remains byte-identical and supplies the independent compression adapter. No upstream certificate or private Vault body is copied or edited here.

From the repository root, use Python 3.11 or 3.12, Node, and separate checkouts at the pinned commits:

```bash
python3.12 -B 04-operator-evolution/verify_r9.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
```

[R9_SOURCE_PINS.json](../04-operator-evolution/R9_SOURCE_PINS.json) fixes both upstream sources and the symbolic input contract. The verifier checks source SHA256 and Git blobs before execution, runs 24 focused exact tests, compares 25 full-family tensor pushforwards and 75 directional information values with unchanged QTH-1, replays nine symbolic certificates, and requires six deliberate mathematical mutations and two native alterations to fail. Recorded R1–R8 and master-review hashes are checked without rerunning their unchanged suites.

[R9_NATIVE_CERTIFICATE.json](../04-operator-evolution/R9_NATIVE_CERTIFICATE.json) contains the symbolic inputs and rewrite witnesses. [R9_VERIFICATION.json](../04-operator-evolution/R9_VERIFICATION.json) binds proof, code, tests, pins and earlier evidence. Written general proofs, rewrite certificates and exact source computations are explicit evidence levels; new Lean formalization, external peer review and physical experiments are separate and are not asserted by the finite PASS.
