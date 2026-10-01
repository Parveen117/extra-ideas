# R13: complementary native memory generates noise, event laws and metric data

Author: Monty Dabas. Development: 1 October 2026.

R12 derived an observer and metric after specifying Gaussian record noise and an event law. R13 constructs those statistical ingredients from native complementary transport instead. It proves an exact memory/force equation, derives hidden-force covariance from a balanced unit-norm state, derives a monitored event law from native branch norms, and realizes the completed R12 metric equation with **native binary records rather than Gaussian noise**.

The central link is constructive. In a native rotation, the sine component is a complementary information channel. Eliminating it produces a visible memory kernel and an unresolved force. Monitoring its norm produces event weights. Reading balanced native amplitude pairs produces a likelihood and Fisher metric. These are related outputs of explicit protocols, rather than independent noise parameters.

## 1. Native amplitudes and the observed cut

Use the established cut amplitude \(z=x+\iota y\), with \(\iota^2=-1\) and path/scalar dagger. In its admitted positive two-mode chart,

\[
N(z)=z^\dagger z=x^2+y^2,\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

No physical Hilbert space or spacetime metric is introduced as primitive. This is the finite positive representation of the existing native cut norm.

The existing KIR basis also derives

\[
J=-RK=\operatorname{diag}(1,-1),\qquad
P=\frac{I+J}{2},\quad Q=\frac{I-J}{2}.
\tag{1}
\]

P reads x and Q retains its complementary y channel. This observer-adapted cut is distinct from treating the matrix K itself as this same coordinate projection.

For \(c^2+s^2=1\), native transport is

\[
T=cI+sR=\begin{pmatrix}c&-s\\s&c\end{pmatrix}.
\tag{2}
\]

Writing \(c=\cos\theta,\ s=\sin\theta\) is an angle chart, not a choice of which component is inherently classical. Exact certificates use rational c,s from the existing Cayley flow:
\(c=(1-a^2)/(1+a^2),\ s=2a/(1+a^2)\).

The visible and complementary channels satisfy

\[
PTP=cP,\qquad QTP=sQRP,\qquad
\boxed{[P,T]=-sK,\quad [P,T]^2=s^2I.}
\tag{3}
\]

Thus the same complementary amplitude measures cut/transport noncommutation and norm leakage:

\[
P-(PTP)^\dagger(PTP)=(QTP)^\dagger(QTP)=s^2P.
\tag{4}
\]

A cut/transport commutator is not automatically the Riemann tensor or the curvature of an unspecified connection. Two rotations generated only by R commute, yet (3) can be nonzero: the observed cut itself is the incompatible action.

## 2. Hidden elimination derives memory and an unresolved force

For any constant invertible native transport split into visible and hidden blocks,

\[
\binom{x_{n+1}}{z_{n+1}}
=\begin{pmatrix}A&B\\C&D\end{pmatrix}
\binom{x_n}{z_n},
\]

solve the hidden recurrence and substitute it back:

\[
\boxed{x_{n+1}=Ax_n+
\sum_{j=0}^{n-1}BD^{\,n-1-j}Cx_j+BD^nz_0.}
\tag{5}
\]

**Proof.** Induction gives
\(z_n=D^nz_0+\sum_{j<n}D^{n-1-j}Cx_j\). The visible row is \(Ax_n+Bz_n\), which gives (5). ∎

The memory kernel and hidden initial force are

\[
K_\ell=BD^\ell C,\qquad f_n=BD^nz_0.
\tag{6}
\]

Neither needs a random-number source. If z0 is known, fn is known. If only x0 is observed, the unresolved initial complement can appear as noise.

For the rotation (2),

\[
\boxed{x_{n+1}=cx_n
-s^2\sum_{j=0}^{n-1}c^{\,n-1-j}x_j
-sc^ny_0.}
\tag{7}
\]

The sine leg both leaves the observed sector and returns from it. Its first round-trip correction is

\[
PT^2P-(PTP)^2=PTQTP=-s^2P.
\tag{8}
\]

This is the existing native N06 memory defect made into an actual evolution law. Dropping it changes the dynamics.

Two visible time records determine the hidden initial component when s is nonzero:

\[
y_0=\frac{cx_0-x_1}{s},\qquad
x_{n+2}=2cx_{n+1}-x_n.
\tag{9}
\]

The unexplained force can therefore become reconstructed memory. With response errors \(\epsilon_0,\epsilon_1\), its sharp recovery radius is
\((|c|\epsilon_0+\epsilon_1)/|s|\). Near-zero complementary gain creates an observation-conditioning problem, consistent with R11's stability tools.

## 3. A native preparation derives the noise covariance

For an arbitrary hidden ensemble, center
\(\eta_0=z_0-\mathbb E[z_0\mid x_0]\), and let its conditional covariance be S. Equation (5) determines

\[
\xi_n=BD^n\eta_0,\qquad
\boxed{\operatorname{Cov}(\xi_n,\xi_m\mid x_0)
=BD^nS(D^m)^TB^T.}
\tag{10}
\]

A covariance is not selected by a transport matrix alone. R13 closes this remaining preparation freedom in an explicit native family.

Specify a known unit total norm and exchange balance under the derived J:

\[
x_0^2+y_0^2=1,\qquad
y_0\in\{+a,-a\},\qquad p_+=p_-=\frac12.
\]

The unit-norm state fixes \(a^2=1-x_0^2\), and exchange balance fixes its conditional mean to zero. Hence

\[
\boxed{\operatorname{Var}(y_0\mid x_0)=1-x_0^2,}
\]

\[
\boxed{\operatorname{Cov}(\xi_n,\xi_m\mid x_0)
=s^2c^{n+m}(1-x_0^2).}
\tag{11}
\]

There is no independent noise variance in (11). State normalization and the native complementary preparation determine it.

For this rotation, the same derived covariance and memory kernel obey

\[
\boxed{\operatorname{Cov}(\xi_\ell,\xi_0\mid x_0)
=-(1-x_0^2)K_\ell,\qquad K_\ell=-s^2c^\ell.}
\tag{12}
\]

This is a finite preparation-relative noise/memory identity, not a claim of thermodynamic equilibrium or a universal fluctuation–dissipation theorem.

Its noise is temporally correlated and has covariance rank one. The balanced hidden pair has fourth moment \(a^4\), rather than the Gaussian value \(3a^4\). Endogenous complement noise therefore does not justify Gaussian or independent-in-time noise. If \(a=0\), the unresolved force is zero while the kernel can remain nonzero; memory is not identical to noise variance.

A hidden record that never returns, B=0, produces neither visible memory nor visible force, even when the hidden state changes. A complementary coordinate becomes relevant noise only through the observer's actual coupling.

## 4. Deriving branch probabilities from native norm weights

To interpret split norms as probabilities, state the event rule precisely:

1. Branch weight depends only on its nonnegative native squared norm.
2. Weights add when an orthogonal branch is refined into disjoint subchannels.
3. Unit total norm has unit weight.

**Theorem — normalized additive norm weights.** These requirements select branch weight \(w(e)=e\) for \(0\le e\le1\).

**Proof.** Additivity gives \(w(m/n)=m/n\); nonnegativity gives monotonicity. Approximating any real e from below and above by rationals forces \(w(e)=e\). No Gaussian assumption enters. ∎

These axioms select a consistent operational probability law. Whether a particular physical detector implements that law is an empirical/protocol question; probability semantics do not follow from bare multiplication alone.

Under this norm readout, a visible unit preparation transported through (2) has

\[
\boxed{r=c^2,\qquad q=s^2=1-r.}
\tag{13}
\]

The variance of a binary event indicator is \(rq=c^2s^2\). This is distinct from the hidden-force covariance (11); both are now derived, with their source and meaning shown.

## 5. Monitored first exit derives the event law

Monitor the cut after each native move: continue on P, terminate on Q. The branch amplitude for k continuations followed by first exit is

\[
H_k=QT(PTP)^kP=sc^kQRP,
\]

so

\[
\boxed{\Pr(K=k)=H_k^\dagger H_k\big|_{\mathrm{visible}}
=q r^k=s^2c^{2k}.}
\tag{14}
\]

The geometric event law is derived from the repeated cut operation. The event count is clock-free. Its finite ledger is

\[
\sum_{k=0}^{M}qr^k+r^{M+1}=1.
\tag{15}
\]

If s=0, the unresolved tail is one and there is no exit. If c=0, exit is immediate. For a varying sequence of native moves,

\[
\Pr(K=k)=s_k^2\prod_{j<k}c_j^2,\qquad
\mathrm{tail}=\prod_j c_j^2,
\tag{16}
\]

so a constant geometric hazard is not imposed on arbitrary dynamics.

Coherent retention is a different policy. With no intermediate cut, two-step visible amplitude is \(c^2-s^2\), while monitored survival amplitude is \(c^2\). At c=3/5, s=4/5 their squared values are respectively 49/625 and 81/625. The geometric law does not describe the fully retained coherent trajectory.

## 6. Event information is derived without a Gaussian likelihood

For regular two-outcome weights \(r=c^2,\ q=s^2>0\), the native angular tangent satisfies \(c'=-s,\ s'=c\). The binary Fisher information is

\[
\boxed{g_{\theta\theta}^{\rm one\ event}
=\frac{(r')^2}{r}+\frac{(q')^2}{q}=4.}
\tag{17}
\]

For a complete first-exit record,

\[
\mathbb E K=\frac r q,\qquad
\operatorname{Var}K=\frac r{q^2},\qquad
\partial_r\log\Pr(K=k)=\frac{k}{r}-\frac1q,
\]

\[
\boxed{g_{rr}^{\rm exit}=\frac1{rq^2},\qquad
g_{\theta\theta}^{\rm exit}=\frac4{s^2}.}
\tag{18}
\]

Expected trial count is \(1/s^2\); information per expected trial remains four. Divergence near rare exit is an observation-cost effect. Endpoint distributions have changing support and are not covered by the regular Fisher formula.

The exact angular scores agree with QTH-1's existing projective classical Fisher function. This consumes that readout unchanged; it does not apply the source's interior-state SLD formula to a pure-state boundary.

## 7. Native binary meters realize the completed R12 metric

Prepare each meter in the balanced native amplitude pair, with density representation

\[
\rho_{\rm bal}=\frac12(I+K)
=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]

This density can be written rationally even though its normalized state-vector entries are \(1/\sqrt2\). Given local response amplitude a, rotate it with the native rational Cayley map \(U(a/4)\), and read P/Q. If its coefficients are \(c_a,s_a\), the exact probabilities are

\[
p_0(a)=\frac{(c_a-s_a)^2}{2},\qquad
p_1(a)=\frac{(c_a+s_a)^2}{2}.
\tag{19}
\]

They are nonnegative, sum to one and have, at a=0,

\[
p_0=p_1=\frac12,\qquad
p'_0=-\frac12,\quad p'_1=\frac12.
\tag{20}
\]

The signed record has zero mean and unit variance; its mean derivative is one. Its likelihood and record variance have been generated by the native amplitude split.

For response a_j = (DT_w delta)_j, separately prepared product meters have independent signed records at the anchored zero contrast. Their score form is

\[
\sum_{j,b}\frac{dp_{j,b}\otimes dp_{j,b}}{p_{j,b}}
=T_w^TD^TDT_w.
\]

Tag each word as in R12. The complete native binary experiment therefore selects

\[
\boxed{G=(1-r)D^TD+r\sum_{\pm}\frac12T_\pm^TGT_\pm,\qquad
\ker G=\bigcap_w\ker(DT_w).}
\tag{21}
\]

This reuses the R12 convergence and quotient proofs for the new likelihood. Gaussian noise and external noise precision are absent from this realization. Product meter preparation, normalized contrast coordinates and the anchored zero-contrast point are stated parts of the construction.

Now use the same monitored rotational gate to set \(r=c^2\), and identify the paired sheet-record gain with its actual complementary leg s:

\[
T_\pm(v)=I\pm svN,\quad N=E_{31},\quad D=(I_3\;0).
\]

Branch exchange fixes the two sign weights to one half. R12's square-zero solution gives

\[
\boxed{G=\operatorname{diag}(1+c^2v^2,1,1,0),\qquad
g=(1+c^2v^2)du^2+dv^2.}
\tag{22}
\]

The old independent coupling and continuation combine as
\(\kappa=s^2r/(1-r)=c^2\). R12's derivative solver supplies the exact first/second metric jets. The derived three-mode observer is \(C_\infty=(I_3\;0)\); the admitted tangent chart uses the first two modes.

Matching sheet gain to the monitored complementary leg is an explicit construction connecting two native registers, not a theorem that every transport must have this gain. All formulas use the selected rotational phase; they do not select that phase universally.

At zero contrast, all conditional meter distributions coincide. Forgetting word tags averages the score, so the same positive information ledger as R12 holds:

\[
G_{\rm untagged,0}=D^TD,\qquad
G-G_{\rm untagged,0}=c^2v^2E_{11}\succeq0.
\tag{23}
\]

This conclusion now uses actual binary records. It is not a replacement of a non-Gaussian mixture by a Gaussian. If the gate angle is also estimated, its first-exit information is (18); at zero contrast its score has zero cross-pairing with the conditional response scores.

## 8. What the complementary mechanism selects

| Quantity | Native derivation | Structure that must be specified |
|---|---|---|
| Visible memory | \(BD^\ell C\) | Transport and observed cut |
| Apparent hidden noise | \(BD^n\eta_0\) | Initial hidden state/preparation |
| Noise covariance in the balanced unit family | \(s^2c^{n+m}(1-x_0^2)\) | Known total norm and exchange balance |
| Event continuation and exit | \(c^2,s^2\) | Additive norm readout and monitoring policy |
| First-exit law | \(s^2c^{2k}\) | Repeated fixed native move |
| Meter likelihood and variance | Native P/Q probabilities; unit contrast variance | Balanced product meter preparation |
| Completed observer and seam metric | (21)–(22) | Matched complementary sheet gain and tangent chart |

Noise and event law are derived quantities in this construction. The surviving physical freedom has moved to concrete native preparation, phase, detector policy and channel identification. Treating those as inspectable native structure is stronger than leaving an arbitrary covariance and stopping probability unrelated.

Squared event norms and noise covariances still lose orientation: s and -s give identical weights and covariance while changing signed trajectories. R11's phase/marker channels remain useful. A complementary state coordinate is also not the independent integer sheet register of R7; neither scalar sine nor covariance reconstructs erased winding.

## 9. Evidence and lineage

The [implementation](../04-operator-evolution/complement_memory_noise.py) is an adapter to the existing native carrier and R12 solver. The [verification record](../04-operator-evolution/R13_VERIFICATION.json) records **31 exact tests**, including 36 full/reduced rotation histories, 108 hidden-component recoveries, 12 independent QTH score comparisons, 72 first-exit branch-norm checks, multi-hidden feedback, non-Gaussian binary meter probabilities and the derived seam jets. It also records ten native equality replays, one distinctness replay, three canonical N03 completions, six rejected mathematical mutations and three rejected native alterations.

The first-exit native replay explicitly reuses the separately replayed projector and single-step identities before squaring its amplitude. The [native packet](../04-operator-evolution/R13_NATIVE_CERTIFICATE.json) records those dependencies; it does not hide them as unproved simplifications.

~~~bash
python3.12 -B 04-operator-evolution/verify_r13.py \
  --publications-root ../Publications \
  --rkf-root ../Recognition-Kernel-Framework
~~~

Use Python 3.11 or 3.12 and the separately pinned RKF source. The [pins](../04-operator-evolution/R13_SOURCE_PINS.json) check source bytes before execution. All R1–R12 and master-review evidence is preserved, including R12's original Gaussian experiment. No canonical engine code is copied or edited. Written derivations, finite tests and native replays are separate evidence levels; physical experiments and a new Lean formalization are not claimed.

Hidden-state elimination and the broad memory/random-force link are established projection mathematics. [Mori's original paper](https://doi.org/10.1143/PTP.33.423) derives projected random forces linked to damping/memory; [Gouasmi, Parish and Duraisamy](https://arxiv.org/abs/1611.06277v2) discuss the resolved memory/noise split and its exact linear case. R13 does not claim priority for that general mechanism.

The added EMK result is the explicit complementary-cut construction, shared sine/cut-commutator strength, norm-selected first-exit law, preparation-derived colored covariance, and non-Gaussian realization of the observer/metric foundation with \(\kappa=c^2\).

- [EMK-2](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/emk-ugd-algebra/certificates/emk2_native_carrier.py): existing native even/odd carrier and Cayley composition.
- [QTH-1](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/papers/curvature-information-duality/certificates/qth1_quantum_recognition_information.py): unchanged projective classical information and quantum/classical distinction.
- [RKF N03/N06](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md): future observers and returning compression memory.
- [R12](OBSERVER_METRIC_FOUNDATION_R12.md): completed positive selection equation, derivative solver and exact metric/connection boundary.

