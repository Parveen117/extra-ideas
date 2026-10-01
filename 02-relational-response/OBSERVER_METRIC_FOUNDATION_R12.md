# R12: derive the observer and metric from a native response law

Author: Monty Dabas. Development: 1 October 2026.

R11 supplied spectral repair tools. R12 supplies a conditional foundation: a specified native event experiment determines its future observer quotient and its operational metric together. The observer is the quotient by every direction invisible to future responses; the metric is the Fisher information of the complete tagged experiment. Neither is inserted as an independent metric matrix.

A balanced, square-zero EMK sheet record gives an explicit geometric result:

\[
\boxed{g= (1+\kappa v^2)\,du^2+dv^2,\qquad
\kappa=\frac{r}{1-r}c^2.}
\tag{1}
\]

The quadratic warp of the existing EMK-G1 family is therefore derived within this event law. The coupling c, continuation probability r, record calibration and tangent interpretation are admitted physical inputs. R12 does not derive those inputs from the bare KIR algebra, identify information distance with a Lorentzian spacetime interval, or select a universal constant.

## 1. What the experiment supplies

Work in an admitted finite real representation V of lawful native path composition. Supply:

| Input | Meaning |
|---|---|
| Invertible maps \(T_a:V\to V\) | Permitted native moves |
| A linear record map \(D:V\to\mathbb R^m\) | Calibrated terminal responses |
| A positive precision \(W=\Sigma^{-1}\) | Measured record-noise covariance |
| \(p_a>0,\ \sum_a p_a=1\) | State-independent probabilities of the permitted moves |
| \(0<r<1\) | Probability of another recognition event; stopping probability is \(t=1-r\) |

Draw a word w of k moves with probability \(t r^k p_w\), where \(p_w\) is the product of its letter probabilities. Apply its ordered transport \(T_w\), then record

\[
Z=D T_w x+\epsilon,\qquad \epsilon\sim N(0,\Sigma).
\tag{2}
\]

Record the word, including its length and branch labels, alongside Z. These tags are part of this experiment. A permitted marker acting on the state is an intervention, as R11 already distinguished from reading a marked response.

Gaussian calibration is a declared likelihood model, not a consequence of noncommutativity. A different likelihood can produce a different information metric. Event count is not a physical clock.

For convergence, exhibit \(H\succ0\) and \(0\le\beta<1\) such that

\[
r\sum_a p_a T_a^T H T_a\preceq\beta H.
\tag{3}
\]

H is a checkable stability witness, not the physical metric being derived.

## 2. The metric selection theorem

Let \(Q=D^T W D\) and define the positive operator on symmetric forms

\[
\mathcal L(B)=r\sum_a p_a T_a^T B T_a.
\]

**Theorem 1 — tagged response metric.** Under (2)–(3), the expected Fisher information exists and is the unique symmetric solution of

\[
\boxed{G=tQ+\mathcal L(G).}
\tag{4}
\]

It is positive semidefinite, and

\[
G=t\sum_{k=0}^{\infty}\mathcal L^k(Q)
 =\sum_w t r^{|w|}p_w\,T_w^TQ T_w.
\tag{5}
\]

**Proof.** For a fixed word, the Gaussian score in direction h is
\(h^T T_w^T D^T W(Z-DT_wx)\). Its expected product with the score in direction l is \(h^T T_w^T Q T_w l\). Word probabilities do not depend on x, so the joint tagged Fisher information is their weighted sum.

Choose \(\alpha\ge0\) with \(Q\preceq\alpha H\). Positivity and (3) imply \(\mathcal L^k(Q)\preceq\alpha\beta^kH\); the series converges. Splitting off the empty word gives (4). On symmetric forms, the order norm induced by H makes \(\mathcal L\) a contraction of factor at most beta. A homogeneous solution \(B=\mathcal L(B)\) must consequently vanish. This proves uniqueness. ∎

The expected recognition balance is

\[
G-r\sum_a p_a T_a^T G T_a=tQ\succeq0.
\tag{6}
\]

Individual transports need not preserve G. Equation (6) is an experiment budget, not a claim that continuous R and K flows share a fixed symmetric metric.

A finite event aperture M has the exact unresolved remainder

\[
G-G_{\le M}=\mathcal L^{M+1}(G)\succeq0,\qquad
G-G_{\le M}\preceq
\frac{t\alpha\beta^{M+1}}{1-\beta}H.
\tag{7}
\]

This connects finite recovery tools to a completed response metric without replacing a tail by zero. The exact implementation solves (4) on symmetric rational forms and checks both its residual and its separate stability contract.

## 3. The observer is selected by the same form

**Theorem 2 — observer/metric co-selection.**

\[
\boxed{\ker G=\bigcap_w\ker(DT_w)=:\mathcal N_\infty.}
\tag{8}
\]

This is the largest subspace of \(\ker D\) invariant under all admitted transports. The quotient \(V/\mathcal N_\infty\) is the minimum future-complete linear refinement of D, and G induces a positive definite metric on it.

**Proof.** Every weight in (5) is positive and W is positive definite. Thus \(x^TGx=0\) exactly when every \(DT_wx=0\). Appending a move to any word proves invariance of this kernel. Conversely, any invariant subspace inside \(\ker D\) is killed by every future word and lies in (8). RKF N03 then supplies the minimum-rank statement. ∎

Choose any independent row basis C of G. Its kernel is (8). If B is any section with \(CB=I\), then

\[
\bar G=B^TGB,\qquad G=C^T\bar G C,\qquad
\bar T_a=CT_aB,\qquad CT_a=\bar T_a C.
\tag{9}
\]

Changing B adds only invisible directions, so the quotient metric and descended transports are independent of that choice. A change of quotient coordinates changes their matrices in the ordinary covariant way. The quotient itself, rather than a particular row basis, is selected.

Unlike G, the future invisible kernel is independent of positive noise weights and positive event probabilities. It depends on the response map and the support of the permitted catalogue. Removing an action or recording fewer responses can change it. State-observer rank is not R11's count of scalar channels reconstructing an operator-valued curvature target.

For example, \(D=(1,1)\), \(T=\operatorname{diag}(1,1/2)\), \(r=1/2\), \(W=1\) give

\[
G=\begin{pmatrix}1&2/3\\2/3&4/7\end{pmatrix}.
\]

One current scalar response therefore selects a two-dimensional future observer. With R11's two reverse markers admitted as moves, two visible seed responses select all four state modes; the exact metric is

\[
\operatorname{diag}(64/63,64/63,8/63,8/63).
\]

Without those reverse actions, the two hidden state directions remain invisible. The unchanged canonical N03 engine independently certifies these observer completions.

## 4. Deriving a metric field and its derivatives

Under a constant native reference change \(x'=Sx\),

\[
D'=DS^{-1},\quad T'_a=ST_aS^{-1},\quad
G'=S^{-T}GS^{-1}.
\tag{10}
\]

A change of record units \(Z'=UZ\) requires
\(D'=UD,\ W'=U^{-T}WU^{-1}\) and leaves G unchanged. Both statements follow by substitution in the uniquely solved equation (4).

For smooth state-dependent transport data, keep r and p fixed. Differentiate (4), rather than supply metric derivatives independently:

\[
(I-\mathcal L)G_i
=tQ_i+r\sum_a p_a
\left(T_{a,i}^TGT_a+T_a^TGT_{a,i}\right).
\tag{11}
\]

The second derivative satisfies

\[
\begin{aligned}
(I-\mathcal L)G_{ij}=tQ_{ij}+r\sum_a p_a\big(&
T_{a,ij}^TGT_a+T_a^TGT_{a,ij}\\
&+T_{a,i}^TG_jT_a+T_a^TG_jT_{a,i}\\
&+T_{a,j}^TG_iT_a+T_a^TG_iT_{a,j}\\
&+T_{a,i}^TGT_{a,j}+T_{a,j}^TGT_{a,i}\big).
\end{aligned}
\tag{12}
\]

An admitted constant tangent map \(E:TU\to V\) now gives
\(g=E^TGE,\ g_i=E^TG_iE,\ g_{ij}=E^TG_{ij}E\). It is a tangent metric only where this pullback is nondegenerate. A varying tangent map needs its own derivative terms; they are not omitted by interpreting this constant-map implementation as general.

Constant response rank on a region is needed for a smooth observer bundle. The tests include a rank jump at v=0, so pointwise quotient selection is not silently promoted to a global bundle theorem.

## 5. The EMK sheet family selects the quadratic seam warp

Use four native coordinates, with two calibrated tangent modes, a recorded sheet mode and one invisible mode. Set

\[
N=E_{31},\quad N^2=0,\quad
D=(I_3\;0),\quad W=I_3,\quad
T_\pm(v)=I\pm cvN.
\tag{13}
\]

The moves are inverse lawful transports. The sector cut J has \(JNJ=-N\). Exchange symmetry of the two allowed branches selects \(p_+=p_-=1/2\), the existing R9 mechanism.

The first moment is I, but the quadratic response satisfies

\[
\frac{T_+^TBT_++T_-^TBT_-}{2}
=B+c^2v^2N^TBN.
\tag{14}
\]

The square-zero relation implies
\(N^T(N^TQN)N=0\). Consequently (4) has the explicit solution

\[
\boxed{G=Q+\frac{r}{1-r}c^2v^2N^TQN
=\operatorname{diag}(1+\kappa v^2,1,1,0).}
\tag{15}
\]

The stability requirement is genuinely met: take
\(H=\operatorname{diag}(1+2\kappa v^2,1,1,1)\) and
\(\beta=(1+r)/2\). The matrix margin in (3) is nonnegative everywhere for this family. Thus (15) is a convergent completed experiment metric, rather than merely a formal fixed point.

Its derived future observer is \(C_\infty=(I_3\;0)\). The coordinate tangent map \(E=(e_1,e_2)\) pulls (15) back to (1). The physical identification of these two directions is supplied by this chart; the observer quotient itself has three retained state modes.

This is a smooth family of local, anchored contrast experiments: x is an infinitesimal displacement at a base point with normal coordinate v. It is not asserted to be one global statistical family whose sufficient statistics have differential \(cv\,du\). That one-form is nonclosed:
\(d(cv\,du)=c\,dv\wedge du\). Its native directional response is permitted to carry memory.

Equations (11)–(12) derive

\[
g_v=\operatorname{diag}(2\kappa v,0),\qquad
g_{vv}=\operatorname{diag}(2\kappa,0),
\]

with every u and mixed partial zero. The R10 metric-jet construction then gives the existing EMK-G1 formulas

\[
\Gamma^u_{uv}=\frac{\kappa v}{1+\kappa v^2},\qquad
\Gamma^v_{uu}=-\kappa v,\qquad
\boxed{K_g=-\frac{\kappa}{(1+\kappa v^2)^2}.}
\tag{16}
\]

The seam v=0 is geodesic by reflection symmetry. Positive record cost derives \(\kappa\ge0\) here. EMK-G1's admissible \(\kappa<0\) strip is not derived by this experiment; it needs a different native response mechanism.

## 6. Balance, forgetting and curvature are different operations

For the R11 connection \(A_u=P_{\rm hidden}, A_v=cN\), raw curvature is cN. Its spectrum is blind to this nonzero square-zero operator. The balanced linear observation removes it because it is cut-odd. The squared cost \(N^TQN\) is cut-even when Q commutes with J, and survives in (15).

**Zero balanced first-moment curvature therefore does not force zero metric curvature.** At v=0 and c nonzero, the metric value is Euclidean but its derived second jet gives \(K_g=-\kappa\ne0\). Values alone cannot certify flatness.

Now actually forget every event tag in experiment (2). At x=0 all conditional Gaussian densities coincide. The marginal score is the average conditional score, so

\[
M=t\left(I-r\sum_a p_aT_a\right)^{-1},\qquad
G_{\rm untagged,0}=(DM)^TW(DM).
\tag{17}
\]

The score covariance gives the exact local information ledger

\[
G=G_{\rm untagged,0}+\Delta,\qquad\Delta\succeq0.
\tag{18}
\]

For (13), \(M=I\), so \(G_{\rm untagged,0}=Q\) and
\(\Delta=\kappa v^2E_{11}\). The anchored tangent information metric becomes I after this tag erasure. This particular experiment loses the warp; this is not an all-observation flatness law. At nonzero displacement x, the conditional means differ, and (17) is not the Fisher information of the Gaussian mixture. Replacing that mixture by a Gaussian with the averaged mean is also a different experiment.

CID-1 already supplies a classical recognized/discarded covariance ledger; QTH-1 already distinguishes a quantum information metric and its antisymmetric commutator part. R12's tagged Gaussian score ledger and event-law selection are a new adapter to those existing foundations, not a proof of uniqueness of all quantum Fisher metrics.

## 7. Relation to the classical curvature special case

The information-derived g does not automatically turn the original native connection into its Levi–Civita connection. On the two-mode corner of \(A_u=P_{\rm hidden}, A_v=cN\), the visible native connection is zero. For the derived metric, R10 certifies

\[
0=R(g)+\mathcal D(-\Gamma(g))+0,
\tag{19}
\]

with nonzero Riemann curvature and an exactly cancelling distortion term. Native-to-Levi–Civita descent fails, as it should.

Conversely, the derived metric can be used in R10's complete descent contract. Extending its computed Levi–Civita jet by an independent K/R hidden sector gives a certified Riemann special sector and retains additional native curvature. That extension is a supplied compatible connection law; information metric selection alone did not force it.

## 8. What is now derived, and what remains physical input

| Structure | R12 result | Remaining choice |
|---|---|---|
| Future observer | Minimum quotient determined by all future response kernels | Permitted moves and terminal response preparation |
| Operational metric | Unique expected tagged Fisher form under the stable event law | Likelihood/noise calibration, event probabilities and stopping law |
| Seam warp | \(\kappa=c^2r/(1-r)\), from balanced nilpotent sheet cost | Coupling, event law and tangent chart |
| Metric derivatives | Derived twice from the selection equation | Smooth source data and constant tangent map in this implementation |
| Classical curvature | Computed from the derived g; exact EMK-G1 agreement | A native connection needs the full R10 descent law |
| Signature and dynamics | Positive information distance | Lorentz signature, physical clock and field equations are unselected |

The residual freedom is testable, not hidden. With the same c=1, v=1 and the same observer, \(r=1/2\) gives \(g=\operatorname{diag}(2,1)\), while \(r=1/3\) gives \(g=\operatorname{diag}(3/2,1)\). These are not related by one overall scale. Bare transport algebra therefore cannot select a unique physical metric without a response law.

This is a foundation derivation conditional on concrete experimental structure. Physical selection of that structure remains the next substantial question; adding another named metric matrix would not answer it.

## 9. Certification and lineage

The [implementation](../04-operator-evolution/observer_metric_foundation.py) and [28 tests](../04-operator-evolution/test_observer_metric_foundation.py) check event convergence gates, exact quotient factorization, state/record covariance, tag loss, all differentiated product terms, 27 sheet parameter cases, 18 independent EMK-G1 source comparisons, feedback observer repair and rank changes. The [verification record](../04-operator-evolution/R12_VERIFICATION.json) also records six rejected mathematical mutations, eight native equality replays, one native distinctness replay, four canonical N03 completions and three rejected native alterations. R1–R11 and master-review evidence stays unchanged.

~~~bash
python3.12 -B 04-operator-evolution/verify_r12.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

Use Python 3.11 or 3.12 and the pinned canonical Node source. The [pins](../04-operator-evolution/R12_SOURCE_PINS.json) bind source bytes, commits and the externally checked native input. The [native certificate](../04-operator-evolution/R12_NATIVE_CERTIFICATE.json) separates abstract rewrite results from finite observer completions. Infinite-series convergence is a written theorem under (3); the runtime checks rational witnesses and finite cases. There is no new Lean formalization, external peer review or physical experiment.

The Fisher score, observability kernel, Stein equation and quotient constructions are established mathematics. R12's contribution is their joint EMK selection contract, its balance/second-moment mechanism, the derived positive quadratic seam family, differentiated geometry, and its certified interface with R9–R11.

- [RKF native algebra, N03](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md): minimum future state observer; unchanged canonical engine.
- [Native metric graph pairing](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/35_native_metric_graph_boundary_pairing.md): existing graph form \(I+C^*C\), boundary pairing and quotient decoder. It admits its Hilbert/graph data; it does not select this experiment.
- [CID-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/cid1_curvature_information_duality.py): existing classical information and metric-free obstruction separation.
- [QTH-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/qth1_quantum_recognition_information.py): existing quantum information/commutator tensor; no quantum metric selection is claimed by R12.
- [EMK-G1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-recognition-geometry/certificates/emkg1_rotational_seam_metric.py): existing seam metric geometry, consumed unchanged.
- [R11](SPECTRAL_CURVATURE_OBSERVER_R11.md): spectral blindness, marked curvature repair and the difference between response reading and intervention.
