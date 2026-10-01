# R17: cut-history counting derives erased-record noise and signed memory return

Research programme: Monty Dabas. Development: 1 October 2026.

The [constructive cut foundation](emk-topology-foundation/README.md) is now included in extra-ideas as the frozen R16 foundation packet. R17 develops its actual cut-role actions. It does not start with a complex scalar, a selected angle, an external noise distribution or a supplied balanced preparation.

The new source is the complete finite tree of the foundation's two cut constructors. Counting each distinct record once derives a consistent record content. Forgetting those records produces a mean operator, a covariance and an exact curvature/spread identity. The first return of the hidden component has a count-derived waiting law and a signed amplitude. The sign is stored in the intervening cut history.

The principal results are

\[
M=\frac{H+K}{2},\qquad M^2=\frac12 I,
\qquad \mathcal V_1(v)=\frac12\|v\|^2
=\frac18\|[K,H]v\|^2,
\tag{1}
\]

\[
\mathcal V_n(v)=(1-2^{-n})\|v\|^2,
\qquad c(L=m)=2^{-(m+1)},
\qquad P H K^m v=(-1)^m P H v.
\tag{2}
\]

Here `c` is normalized **record count**, and the norm is the already constructed role pairing. The statements concern this explicit source tree, not empirical detector frequencies or every possible UGD representation. Defining the whole source tree and counting its records does not prove that nature samples it uniformly. Every identity below is unconditional for the constructed census, with no adjustable angle, covariance or stopping parameter.

## 1. Foundation inputs and the single source

Use C1, C3, C4, C6–C8 of the [R16 manuscript source](emk-topology-foundation/emk_topology_foundation.tex):

- the two typed event tags and their full free records;
- signed/refined cut counts and their completed ordered radial field;
- the free two-role carrier with basis `e_k,e_c`;
- role parity `H`, role exchange `K`, and their derived quarter-turn `R=KH`;
- the coefficient-matching pairing `B` and its squared norm.

The active carrier is the signed/radial role carrier
\(V=\mathbb R_\Sigma e_k\oplus\mathbb R_\Sigma e_c\), not the full UGD state. In particular,

\[
H^2=K^2=I,\quad H^\dagger=H,\quad K^\dagger=K,
\quad HK=-KH,\quad R^2=-I.
\tag{3}
\]

These are proved on the free role generators in C4 and C6. The array representations used by the checker are generated from those same signed generator images. They are not a second definition of the source.

A chronological record \(w=(g_1,\ldots,g_n)\), with each \(g_j\in\{H,K\}\), acts by \(U_w=g_n\cdots g_1\). Matrix/operator products act rightmost first. The interpretation of the original event tags is the explicit R16 representation `k → H`, `c → K`. It is not an assertion about arbitrary raw primitive maps named chi and kappa.

Let \(W_n\) contain every length-\(n\) record, once each. Define, directly by the native counting construction,

\[
c_n(A)=\frac{\#A}{\#W_n},\qquad
\langle f\rangle_n=\frac1{\#W_n}\sum_{w\in W_n}f(w).
\tag{4}
\]

These are definitions on actual finite records. No stochastic library, Born rule or probability density is used to obtain them.

### R17.1 — Consistent census and record-preserving lift

There are \(2^n\) source records. A length-\(j\) prefix has exactly \(2^{n-j}\) extensions to depth \(n\), so its content is \(2^{-j}\), independently of the final depth. Disjoint-event additivity and normalization follow from finite counting.

For any \(v\in V\), the tagged lift \(v\mapsto(U_wv)_{w\in W_n}\), with the average of the native pairings on its record sectors, preserves \(\|v\|^2\). Each individual tagged output determines the initial role state by the reversed sequence of the same involutions.

**Proof.** Appending either tag doubles the record count. Fixing a prefix leaves precisely its suffix choices, proving the cylinder count and consistency. Both H and K preserve the pairing; induction gives \(U_w^\dagger U_w=I\). Averaging their identical squared norms proves the isometry. The inverse is \(U_w^{-1}=g_1\cdots g_n\). ∎

The tagged lift is a mathematical record construction. Erasing its labels is a separate operation. A scalar average never reconstructs a deleted word by itself.

## 2. Curvature, mean retention and erased-record spread

For the one-step source census, define its mean map and its difference map:

\[
M=\frac{H+K}{2},\qquad D=\frac{H-K}{2}.
\tag{5}
\]

Thus \(Hv=Mv+Dv\) and \(Kv=Mv-Dv\). The coefficients one half follow from the two actual records in \(W_1\).

### R17.2 — Exact curvature/spread identity

For every role state,

\[
M^2=D^2=\frac12I,\qquad MD=\frac12R,\qquad DM=-\frac12R.
\tag{6}
\]

Let \(\Omega=[K,H]=KH-HK=2R\). The mean energy and the centered record spread satisfy

\[
\|Mv\|^2=\frac12\|v\|^2,
\qquad
\mathcal V_1(v):=\frac{\|Hv-Mv\|^2+\|Kv-Mv\|^2}{2}
=\|Dv\|^2
=\frac18\|\Omega v\|^2.
\tag{7}
\]

Consequently \(\|Mv\|^2+\mathcal V_1(v)=\|v\|^2\).

**Proof.** Expand (5), use \(H^2=K^2=I\) and \(HK+KH=0\). This gives (6). H and K are self-adjoint for the constructed pairing, so M and D are too; their squares give the two norm factors. Also \(R^\dagger R=I\), so \(\Omega^\dagger\Omega=4I\). Substitution gives (7). ∎

This is a precise relation between native ordered mismatch and variance caused by erasing the two cut records. It is not a Riemann tensor identity or a universal thermodynamic fluctuation law. The spread is completely fixed by the source and readout; it is not an independently fitted noise strength.

## 3. Entire source trees, covariance and retained count information

Write \(v_w=U_wv\), \(m_n=\langle v_w\rangle_n\), \(E=\|v\|^2\), and

\[
C_n=\left\langle(v_w-m_n)(v_w-m_n)^\dagger\right\rangle_n.
\tag{8}
\]

The outer product is the native rank-one map \(u\mapsto v_w B(v_w,u)\). It introduces no physical density operator.

### R17.3 — Exact finite-depth mean and covariance

For every \(n\ge0\),

\[
m_n=M^nv,\qquad
\|m_n\|^2=2^{-n}E,\qquad
\mathcal V_n:=\langle\|v_w-m_n\|^2\rangle_n=(1-2^{-n})E.
\tag{9}
\]

For every \(n\ge1\),

\[
\boxed{C_n=\frac E2 I-m_nm_n^\dagger.}
\tag{10}
\]

If the two depths are read along the same source record, then for \(n,r\ge0\),

\[
C_{n+r,n}:=\left\langle(v_{n+r}-m_{n+r})(v_n-m_n)^\dagger\right\rangle_{n+r}
=M^r C_n.
\tag{11}
\]

Here \(C_0=0\). Thus erased-state uncertainty has a derived lag structure; it is not assumed white or Gaussian.

**Proof.** Each prefix has exactly its two one-step extensions. Averaging those gives M applied to the prefix state, hence the mean recursion. Since \(M^\dagger M=I/2\), norm induction gives the second identity in (9). Expanding the centered square and using the mean gives the last identity.

For a radial state \(u=ae_k+be_c\), the two branch columns are \((a,-b)\) and \((b,a)\). Expanding their native outer products gives

\[
\frac12(Huu^\dagger H+Kuu^\dagger K)
=\frac{a^2+b^2}{2}I.
\tag{12}
\]

Every prefix has norm E, so its next-step second moment is \(EI/2\). Subtracting the mean outer product proves (10). For (11), fix the length-n prefix and average its length-r suffixes. Their mean endpoint is \(M^r v_n\) by the same count recursion. Finite distributivity over prefixes proves (11). ∎

For nonzero v and \(n\ge1\), Cn has eigenvalues \(E/2\) and \(E(1/2-2^{-n})\), on the direction perpendicular to mn and the direction mn respectively. This follows by applying (10) to those two explicit directions. At depth one the covariance has rank one; from depth two it has rank two. The covariance is always a finite source sum, so a continuous noise distribution has not been inserted.

Using the existing native positive logarithm [N2 below], the mean information per full record is also exact:

\[
\mathcal I_n:=\left\langle\operatorname{Log}_\Sigma(1/c_n(\{w\}))\right\rangle_n
=n\operatorname{Log}_\Sigma2
=\operatorname{Log}_\Sigma\frac{E}{\|m_n\|^2}\quad(E>0).
\tag{13}
\]

The proof is substitution of the finite word count and the native logarithm product law. This is a defined record-count information, not a derivation of physical heat or thermodynamic entropy units. Equation (13) quantifies the information omitted by the mean in this census; that information stays in the tagged records.

M itself is invertible: \(M^{-1}=2M\). At a known finite depth the exact mean therefore still determines the initial role state. The omitted information in this section is the realized cut record and its endpoint fluctuation, not a claim that the noiseless mean makes the initial state unidentifiable. Retaining a state and retaining its path are different obligations.

## 4. The derived seam observer and exact memory recovery

Use the C7 seam aperture

\[
P=\frac{I+K}{2},\qquad Q=I-P.
\tag{14}
\]

Let \(u=e_k+e_c\), \(z=e_k-e_c\). Write a role state as \(v=a u+b z\). These are constructed role sums, not a new primitive complex chart. Then

\[
P(a u+b z)=a u,\quad H(a,b)=(b,a),\quad K(a,b)=(a,-b).
\tag{15}
\]

### R17.4 — One returning cut completes the observer

The following identities hold on every role state:

\[
HPH=Q,\quad P+HPH=I,\quad
v=Pv+H(PHv),\quad
\|Pv\|^2+\|PHv\|^2=\|v\|^2.
\tag{16}
\]

Thus the pair of responses \((Pv,PHv)\) reconstructs the state isometrically. One additional scalar response is necessary and sufficient beyond P to recover all role states and all their future H/K responses.

**Proof.** Anticommutation gives \(HKH=-K\), so conjugating (14) by H gives Q. The reconstruction formula and norm identity follow from the complementary orthogonal projectors and H's isometry. P has a one-dimensional kernel, while the two coefficient rows of P and PH are independent. A single scalar refinement of P cannot distinguish its nonzero kernel; these two rows do. The canonical native N03 closure reproduces ranks \(1\to2\). ∎

This is an exact source-derived recovery channel, with no inverse small sine factor. The response PH means applying the known H action before the same P readout. It is a specified mathematical response catalogue, not a guarantee that an arbitrary laboratory measurement is nondestructive.

For a known chronological cut sequence, (15) already gives the full memory equation: K preserves the visible coefficient and reverses the hidden one; H exchanges them. Discarding b changes the next H response. The returning component is therefore a retained state coordinate before it becomes an unresolved observed fluctuation.

## 5. A balanced hidden preparation derived from the first cut histories

### R17.5 — The first source layer generates the hidden pair

Take the distinguished first-role state \(v_0=e_k\). The two records in \(W_1\) give

\[
Hv_0=e_k=\frac12u+\frac12z,
\qquad Kv_0=e_c=\frac12u-\frac12z.
\tag{17}
\]

They have the same P response, distinct opposite Q responses, unit native total norm, and exactly equal count content. The hidden coefficient has mean zero, variance \(1/4\), and fourth centered moment \(1/16\). Its native squared-norm spread is \(1/2\). Applying H makes that previously hidden signed coefficient visible.

**Proof.** Apply the already derived generator maps to e_k and use the two role sums. There is one record of each type, so every displayed moment is its two-term native count sum. H swaps the two coordinates by (15). ∎

This constructs the balanced unit hidden pair which R13 had to specify. The construction follows the first layer of the free cut source; it does not assert that all preparations are balanced. Its fourth moment is the square of its variance, so replacing it by a Gaussian law with the same variance changes its moments. Retaining the first tag, or completing the state response through (16), removes this unresolved ambiguity.

## 6. First return: a derived event law with an active sign ledger

Start with any \(v=a u+b z\). Let L be the number of K tags before the first H tag. At a finite depth N, keep all first-return records with \(L<N\), together with the unresolved `K`-only tail. This is a finite prefix partition of the same source tree, not a new stopping distribution.

### R17.6 — First-return count law and exact tail

For \(0\le m<N\),

\[
c(L=m)=2^{-(m+1)},\qquad c(L\ge N)=2^{-N},
\tag{18}
\]

and

\[
\sum_{m=0}^{N-1}2^{-(m+1)}+2^{-N}=1.
\tag{19}
\]

Native completion sends the tail to zero. The mean number of tags through first return is 2.

**Proof.** The first-return cylinder consists of m fixed K tags followed by one H; R17.1 gives its content. Its remaining suffix is unrestricted. The tail has N fixed K tags. Disjoint prefix counting gives (19). For the mean,
\(\sum_{m=0}^{N-1}(m+1)2^{-(m+1)}=2-(N+2)2^{-N}\), obtained by subtracting twice-shifted finite sums. Its explicit remainder tends to zero in the completed radial field. ∎

H is a source tag that returns the hidden sector. It can occur even when b=0; an H event and a nonzero measured signal are different predicates. This count law is also different from R13's norm-weighted monitored P/Q first-exit law.

### R17.7 — Returning memory remembers the intervening cuts

The visible first-return response is

\[
\rho_m(v)=P H K^m v=(-1)^m P H v=(-1)^m b u.
\tag{20}
\]

Keeping m, or just its parity for this response, recovers b exactly. If the event depth is erased, the limiting count contents are

\[
c(m\text{ even})=\frac23,\qquad c(m\text{ odd})=\frac13,
\tag{21}
\]

and the mean response, second moment and native centered spread are

\[
\bar\rho=\frac13 b u,\qquad
\langle\|\rho\|^2\rangle=\|b u\|^2,\qquad
\langle\|\rho-\bar\rho\|^2\rangle=\frac89\|b u\|^2.
\tag{22}
\]

At finite depth N the signed return sum has the exact remainder

\[
\sum_{m=0}^{N-1}(-1)^m2^{-(m+1)}P H v
=\frac{1-(-1/2)^N}{3}P H v.
\tag{23}
\]

**Proof.** Since K fixes u and reverses z, every intervening K reverses the hidden coefficient. H then brings it into P, proving (20). Alternatively \(PHK=-PH\), and induction gives the same result. Summing even and odd finite geometric packets and taking their explicit vanishing remainders gives (21). Their difference is 1/3 and their sum is one, proving (22). Finite geometric multiplication proves (23). ∎

The factors 1/3 and 8/9 are exact consequences of this source-return census. They are not fitted constants and are not identified with a physical coupling. In particular, small mean response can coexist with full retained return amplitude; the omitted parity record accounts for the difference.

## 7. State completeness still leaves event history active

### R17.8 — Recovering the role state does not recover an erased history

Distinct source records can have the same evaluated operator and the same completed observer response. For example, the empty record and `HH` differ as free records but both act as identity. At any fixed depth \(n>3\), \(2^n\) records map into at most eight role operators, so some distinct equal-length records also share an evaluated operator.

**Proof.** The first example uses \(H^2=I\) and free-word length. The equal-depth conclusion follows by counting the eight operators constructed in C4. Every state-only observer factors through the evaluated action, so it cannot distinguish such records for any initial role state. ∎

The complete cut history, the present hidden role component, and a return-parity summary are different memory objects. R17 does not replace UGD seam memory by a scalar, a pair of state coordinates, or the parity bit. C11/C12's completed history and truncation hierarchy remain the record carrier.

## 8. How much a finite history aperture recovers

At a fixed endpoint depth N, keep only the first j tags, with \(0\le j\le N\). This is the actual C12 truncation aperture on the length-N source histories. For a retained prefix p, define its endpoint estimate by the count-average over all suffixes consistent with p:

\[
\widehat v_{N|p}:=2^{-(N-j)}\sum_{s\in W_{N-j}}U_sU_pv.
\tag{24}
\]

### R17.9 — Exact recovery budget of the native history tower

Every prefix fibre has

\[
\widehat v_{N|p}=M^{N-j}U_pv,\qquad
\|\widehat v_{N|p}\|^2=2^{-(N-j)}E,
\tag{25}
\]

and its remaining endpoint spread is

\[
\mathcal V_{N|j}=(1-2^{-(N-j)})E.
\tag{26}
\]

The spread resolved by keeping the prefix, relative to the completely erased record, is

\[
\left\langle\|\widehat v_{N|p}-m_N\|^2\right\rangle_j
=\big(2^{-(N-j)}-2^{-N}\big)E.
\tag{27}
\]

These two amounts add exactly to \((1-2^{-N})E\). Within every fibre, the count-mean estimate uniquely minimizes the average squared native error among constant role-state estimates. The remaining suffix count information is \((N-j)\operatorname{Log}_\Sigma2\), using that same suffix count.

**Proof.** The suffix census is the R17.3 construction with initial state \(U_pv\), whose energy is E. This proves (25) and (26). Averaging the fibre estimates over prefixes gives \(m_N\). Expanding their centered squares gives (27), and adding proves the budget. For any alternative estimate a, the finite identity
\(\langle\|v_{ps}-a\|^2\rangle_s=
\langle\|v_{ps}-\widehat v_{N|p}\|^2\rangle_s+
\|\widehat v_{N|p}-a\|^2\)
follows because the centered suffix sum is zero. Positive definiteness proves uniqueness. The information count follows from the number of suffix records and the native log product law. ∎

This is an optimal estimate on each actual source fibre, derived without a Gaussian likelihood. It quantifies endpoint uncertainty; the full word still has distinctions which no role-state estimate can recover. A prefix is a canonical aperture, not necessarily the smallest target-specific record summary.

## 9. Count content fixes the boundary recognition metric

For an infinite continuation h, let `[p]` denote the cylinder of all infinite histories with finite prefix p. Its content is already fixed by R17.1: \(c([p])=2^{-|p|}\). This content is consistent across all finite depths; no general measure-extension theorem is needed for these cylinder statements.

### R17.10 — The prefix metric is the count of the common recognition cylinder

For distinct infinite histories h and g, let p be their longest common prefix. Define

\[
d_c(h,g):=c([p]),\qquad d_c(h,h)=0.
\tag{28}
\]

Then \(d_c(h,g)=2^{-|p|}\), exactly the boundary restriction of C10's prefix ultrametric. It generates the cylinder recognition topology and obeys
\(d_c(h,u)\le\max\{d_c(h,g),d_c(g,u)\}\).

**Proof.** Prefix counting gives the numerical equality. Agreement of h and g on n tags and of g and u on n tags forces agreement of h and u on those n tags. The common-prefix lengths therefore satisfy the minimum inequality, which becomes the displayed ultrametric inequality after applying the count formula. Fixing a cylinder is the same finite-prefix condition as the corresponding balls, proving the topology. ∎

Thus native halving on the infinite continuation boundary is fixed once distance is defined as common-prefix cylinder content in this source census. C10 allowed other compatible metric normalizations; R17.10 identifies the count-calibrated one. It does not claim that topology alone selects this distance. Termination is not a third stochastic branch of \(W_n\); the metric on finite terminated histories retains its C10 extension contract. No physical length unit is selected.

## 10. What this advances, and what remains a later identification

R13 proved generic hidden-state elimination but supplied a balanced preparation and an observation protocol. R14 used a scalar specialization and an arbitrary source refinement ledger. R17 supplies a concrete earlier source: the full cut-history tree, its derived H/K actions, the C7 seam aperture, and the native count of its own records.

Within that source, the hidden pair, its moments, mean-loss covariance, the returning event law, its signed memory, the finite-aperture recovery budget and the count-calibrated boundary metric all follow. No independent theta, Gaussian law, hidden variance or geometric stopping parameter is introduced. The same noncommuting operators give the curvature/spread equality (7), while the same retained tag gives the recovery rule (20).

The choice to study the full source census is explicit. This is not a theorem selecting every actual physical trajectory, a Born detector law, a thermodynamic equilibrium preparation or a spacetime signature. The original raw-closure claim remains outside the proved theory exactly as recorded in R16. These boundaries are not extra inputs to any displayed result.

The next mathematical target is the smallest target-specific record quotient: determine which word distinctions must remain for a selected response, beyond the canonical prefix aperture just quantified. This can connect to the existing observer/metric machinery without supplying an external noise model, while preserving the difference between state completeness and full path memory.

## 11. Verification and premise-labeled citations

Run:

```bash
python3.12 -B 04-operator-evolution/verify_r17.py \
  --rkf-root /absolute/path/to/Recognition-Kernel-Framework
```

[R17_VERIFICATION.json](../04-operator-evolution/R17_VERIFICATION.json) records the finite source-tree checks, all coefficient cases of the quadratic identities, symbolic replays, minimum observer completion, return-tail and signed-memory checks, and rejected alterations. [R17_SOURCE_PINS.json](../04-operator-evolution/R17_SOURCE_PINS.json) binds the unchanged foundation packet and canonical engine before the adapter loads. Written induction and geometric-tail proofs cover arbitrary depths; finite enumeration is not described as an infinite proof. No Lean verification or physical experiment is claimed.

- **R16 C1/C3/C4/C6–C8 — native constructive input.** [Frozen foundation at commit `4cc46d1`](https://github.com/Parveen117/extra-ideas/tree/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation). Typed recursion, role completion and pairing are declared constructions with written proofs; arbitrary raw chi/kappa closure is not a premise.
- **RKF N03/N06 — native algebraic input, no Hilbert primitive.** [Canonical algebra theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md). N03 supplies the minimum future-state observer; N06 supplies the cut-corner interpretation. Their action and readout inputs are constructed explicitly here.
- **RKF N09 — directly proved block-memory lineage.** [Completion and memory theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/02_COMPLETION_AND_MEMORY.md). R17 uses its algebraic idea only; its optional Hilbert/GNS and later analytic representation machinery are not dependencies.
- **N2: RKF F00G — native analytic citation.** [Positive logarithm](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00G_NATIVE_LOGARITHM_AND_POWERS.md). Used only for (13), after C8's ordered completion. No ordinary logarithm is imported; the existing finite F00GHI audit is replayed through the foundation verifier.
- **R13/R14 — scoped comparison, not runtime foundations.** Their generic memory/preparation/readout results explain the gap advanced here. No physical likelihood, Gaussian, Fisher metric or scalar-source premise from those editions is imported into (1)–(28).

Finite averaging, covariance expansion and geometric sums are standard mathematical mechanisms, proved here from native counts and the derived operator relations. No priority claim is made for those general mechanisms.
