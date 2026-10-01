# R19: native current retention, returning pair memory and its exact horizon

Research owner: Monty Dabas. Development: 1 October 2026.

**Certification rule:** every mathematical input to a derived claim must trace
to the native source or to an explicitly constructed object over that source.
A classical evolution equation, probability rule, Hilbert carrier or physical
metric cannot enter as an unproved premise. A construction is not silently
upgraded into a law selected by nature. The machine-readable
[derivation ledger](../04-operator-evolution/R19_DERIVATION_LEDGER.json) enforces
these dependency categories; the written proofs below supply the mathematics.
The ledger checker audits the declared graph, not the truth of arbitrary prose.

R18 derived outgoing energy and current from cut histories. Its local pair
(rho,j) suffices for the next outgoing energies, but it does not determine the
next current. R19 derives the missing pair record, its returning-memory law,
and the exact separation needed to predict a finite horizon. Repeated removal
of that record reproduces R18's count transport. Keeping it produces a current
return two events after the first seed and changes the energy distribution on
the third event. No noise strength, coupling or reset probability is fitted.

This is the next native retention result toward physical selection. The
response propagation is still linear. Its quadratic cross terms are not being
renamed a new nonlinear force or an already selected physical interaction.

## Single source and definitions

We consume unchanged [R16 C3–C8](emk-topology-foundation/README.md),
[R17](CUT_HISTORY_MEMORY_R17.md) and [R18](CUT_TRANSPORT_METRIC_R18.md).
Their native two-role module has pairing B, radial refinement coefficients,
and the derived operators

\[
H^2=K^2=I,\quad HK=-KH,\quad R=KH,\quad R^2=-I.
\]

Use R18's existing address x and its finite address-response module. Set

\[
F=H+K,\quad \Pi_\pm=(I\pm H)/2,\quad
L_\pm=\Pi_\pm F,\quad C=F/\sqrt2.
\]

R18 already derived the update

\[
\psi_{n+1}(x)=\frac1{\sqrt2}
 \{L_+\psi_n(x-1)+L_-\psi_n(x+1)\}.
\tag{19.0}
\]

In its ordered roles, F sends (a,b) to (a+b,a-b). The positive square root
comes from R16's completed refinement field. Pair expressions below contain
only its square, so the finite verifier uses rational native coefficients.
All sums have finite support at each finite event. No infinite completion or
continuum equation is needed in R19.

The preparation class for the general observer statements is the finite
signed/refined address-response module already used in R18.3. It contains
finite sums and translations of the two source roles. The single original
seed psi_0(0)=e_k is treated separately wherever a source-specific value is
claimed. Counterexamples on the larger module do not assert that the fixed
single-seed trajectory has an unknown initial state.

R18's address target, addition within an address and finite sum pairing are
inherited explicit constructions. R19 does not claim those choices have now
been physically selected. The full cut-word source remains upstream of this
observation target, as R19.9 makes explicit.

## Written results

### R19.1 — Pair records and their update are derived by native multiplication

For a finite response psi define the pair record by its action on a role v:

\[
G(x,y)v=\psi(x)B(\psi(y),v).
\tag{19.1}
\]

This uses only the already derived B and native scalar multiplication. In the
marked-role coordinates its entries are G_ab(x,y)=psi_a(x)psi_b(y). In
particular G(y,x)=G(x,y)^dagger. Substituting (19.0) gives the exact law

\[
\boxed{(\mathcal T G)(x,y)=\tfrac12
 \sum_{s,t\in\{+1,-1\}}L_sG(x-s,y-t)L_t^\dagger.}
\tag{19.2}
\]

For the coordinate sum tr A=A_kk+A_cc define

\[
\rho(x)=\operatorname{tr}G(x,x),\qquad
j(x)=\operatorname{tr}(K G(x,x)).
\]

These are R18's a²+b² and 2ab. Their total energy is preserved by T.
Equation (19.2) extends by finite linearity to the module of finite pair
records. No independent state-matrix or measurement postulate is introduced.

**Proof.** Expand the two factors of (19.1) using (19.0) and distribute their
finite sums. The two factors 1/√2 give 1/2, and reversing the second factor
gives its native dagger. This proves (19.2), including its order. Directly,
F^dagger F=2I, Pi_+Pi_-=0 and Pi_++Pi_-=I. In the sum of diagonal traces after
the update the s≠t terms have zero trace: a product mapping between distinct
roles has zero diagonal coefficient. The s=t terms sum to tr G(y,y) after
shifting the address and using (1/2)Σ_s L_s^dagger L_s=I. This proves energy
preservation even on the finite linear pair module. Positivity is asserted
only for records actually built from (19.1), or sums of such records with
nonnegative native coefficients; it is not an axiom for every pair array.

### R19.2 — Local current requires a joint signed record at two addresses

Writing psi_n(x)=a_n(x)e_k+b_n(x)e_c, the next current is

\[
\boxed{j_{n+1}(x)=
 [a_n(x-1)+b_n(x-1)]
 [a_n(x+1)-b_n(x+1)].}
\tag{19.3}
\]

Thus even the full family of local quadratic records G_n(x,x), at every x,
does not determine the next current for the finite preparation module.

**Proof.** The first output role in (19.0) comes from x-1 and the second
from x+1. Multiplying them and applying j=2ab cancels the factor 1/2 and
gives (19.3). For an exact witness take

\[
f_+(-1)=e_k,\quad f_+(1)=e_k,\qquad
f_-(-1)=e_k,\quad f_-(1)=-e_k,
\]

with zero values elsewhere. All their local pair records coincide, so every
local quadratic response coincides. After one source update their values at
zero are (e_k+e_c)/√2 and (e_k-e_c)/√2. The currents are respectively +1 and
-1. The initial total energy is two in either case; equal normalization of
both preparations leaves the separation intact. Their two-address record
G(-1,1) distinguishes them before transport. R18.8's minimum two local
responses concerns outgoing energies at that event; it did not establish a
closed autonomous evolution of those two responses.

### R19.3 — Erased pair memory has an exact return equation

Define the local-record cut and its complement by existing address equality:

\[
(\mathcal P G)(x,y)=\mathbf1_{x=y}G(x,y),\quad
\mathcal Q=I-\mathcal P.
\]

Then P²=P, Q²=Q and PQ=QP=0. Split the already derived T into maps between
their images,

\[
A=\mathcal P\mathcal T\mathcal P,\quad
B=\mathcal P\mathcal T\mathcal Q,\quad
C_m=\mathcal Q\mathcal T\mathcal P,\quad
D=\mathcal Q\mathcal T\mathcal Q.
\]

With d_n=P G_n and q_n=Q G_n, the complete local evolution is

\[
\boxed{d_{n+1}=A d_n+B D^nq_0+
 \sum_{r=0}^{n-1} B D^{n-1-r}C_m d_r.}
\tag{19.4}
\]

Every coefficient is fixed by the source update and equality cut. There is
no freely supplied memory kernel. In particular

\[
\mathcal P\mathcal T^2\mathcal P
 -(\mathcal P\mathcal T\mathcal P)^2
=\mathcal P\mathcal T\mathcal Q\mathcal T\mathcal P=B C_m.
\tag{19.5}
\]

**Proof.** Testing equality twice proves the cut identities. Split
G_n=d_n+q_n in (19.2): d_(n+1)=A d_n+B q_n and
q_(n+1)=C_m d_n+D q_n. Repeated substitution gives
q_n=D^n q_0+Σ_(r<n)D^(n-1-r)C_m d_r. Substituting this into the first equation
proves (19.4). Inserting P+Q between two T factors proves (19.5). These are
the native N09 and N06 identities instantiated on derived pair records;
their proof is reproduced here by finite composition. The initial term cannot
be removed for preparations with q_0≠0, as R19.2 already witnesses. P is a
record cut, not the canonical proper depth seam or the physical action of a
detector selected by nature.

### R19.4 — The original single seed generates a returning current on event two

For psi_0(0)=e_k, G_0 is local and q_0=0. The first omitted pair records are

\[
q_1(1,-1)=\tfrac12 e_k e_c^\dagger,\qquad
q_1(-1,1)=\tfrac12 e_c e_k^\dagger.
\]

Their returned local record is exactly

\[
\boxed{(B C_m d_0)(0,0)=K/4,}
\tag{19.6}
\]

with zero other local blocks. Its energy trace is zero but its current
trace is 1/2. The coherent and repeatedly cut sources consequently have
the same energy distribution through event two, while on event three

\[
\rho_3(x)-p_3(x)=
 \begin{cases}+1/4,&x=1,\\-1/4,&x=-1,\\0,&\text{otherwise}.\end{cases}
\tag{19.7}
\]

**Proof.** Equation (19.0) gives psi_1(1)=e_k/√2 and
psi_1(-1)=e_c/√2, proving the omitted records. Under (19.2) only the inward
role shifts can return each of those pairs to a diagonal address, and both
return at zero. Their contributions are e_c e_k^dagger/4 and its dagger;
the sum is K/4. Since tr K=0 and tr K²=2, the energy/current traces follow.
The second-event full local record at zero is
(e_k+e_c)(e_k+e_c)^dagger/4; the repeatedly cut one has diagonal entries
1/4,1/4 and zero cross-role entries. Their energy agrees and their current
differs by 1/2. R18's outgoing formula (rho±j)/2 then produces the ±1/4
third-event difference, in agreement with its direct eight-word census.
Thus an initially absent hidden record is generated, returns visibly as
current, and only then changes the energy readout.

### R19.5 — Cutting the joint record every event derives the count transport

Let d^cut_0=G_0 for the single seed and recursively apply the available cut
after each event:

\[
d^{\rm cut}_{n+1}=\mathcal P\mathcal T d^{\rm cut}_n.
\tag{19.8}
\]

The cut output has only diagonal address blocks and, after every event, only
diagonal role entries. Its total local energy obeys

\[
\rho^{\rm cut}_{n+1}(x)=\tfrac12[
 \rho^{\rm cut}_n(x-1)+\rho^{\rm cut}_n(x+1)],
\quad \rho^{\rm cut}_n(x)=p_n(x).
\tag{19.9}
\]

This constructs R18's count transport from a definite native record operation.
Choosing to execute that cut at every event is explicit; its occurrence in a
physical interaction has not been inferred from the identity.

**Proof.** For a local input, G(x-s,x-t) in (19.2) is zero unless s=t when
the output address pair is (x,x). Each surviving term is projected into a
single role, so no cross-role local entry remains. If a local block is
diag(u,v), native multiplication gives
Pi_±F diag(u,v) F Pi_±=(u+v)Pi_±. The factor 1/2 in (19.2) makes its two
outgoing energies each (u+v)/2. The initial e_k record is already diagonal,
so induction proves (19.9) from event zero. R18.2 derives the same recurrence
and initial value from literal source counts, hence induction makes the two
values identical at every event. No stochastic law, collapse postulate or
external diffusion equation is used as a premise.

### R19.6 — The retained current exactly drives the count-versus-response residue

Define eta_n(x)=rho_n(x)-p_n(x) using the two outputs already constructed.
The coherent source satisfies

\[
\boxed{\eta_{n+1}(x)=\tfrac12[\eta_n(x-1)+\eta_n(x+1)]
 +\tfrac12[j_n(x-1)-j_n(x+1)].}
\tag{19.10}
\]

Initially eta_0=0, and Σ_x eta_n(x)=0 for every n. The apparent excess or
deficit relative to source counts has an exact signed-current source.

**Proof.** Substitute r=(rho+j)/2 and l=(rho-j)/2 into R18.4's local update.
Subtract (19.9). The current terms are a difference of two shifted finite
sums, so their total vanishes. Induction from eta_0=0 proves the sum rule.
The current itself is fixed by (19.3) and the retained pair evolution; it is
not a fitted noise term. Equations (19.6)–(19.7) exhibit its first nonzero
single-seed contribution. This residue is deterministic native response
information. A random-noise or measured-frequency interpretation would need
its own derivation and is not certified by renaming eta.

### R19.7 — Finite-horizon current needs a sharp pair-separation radius

For h≥1, the current j_h(x) depends only on initial pair records with addresses
in [x-h,x+h], of the required parity, and separation at most 2h. The bound
2h is sharp uniformly over the finite preparation module: retaining all
records of separation strictly less than 2h cannot determine that current.
For energy rho_N(x), N≥1, separation at most 2(N-1) suffices; for N≥2 this
bound is also sharp. These are bounds on a record's pair separation, not
counts of scalar channels or a minimality result for every nonlinear decoder.

**Proof.** Each source step shifts a role by one. After h events a value at
x receives only initial roles from x-h,x-h+2,...,x+h. Expanding its quadratic
current gives the asserted sufficient pair window. For necessity let f_±
be supported at -h and h, with values e_k and ±e_k. Their initial pair
records agree at every separation <2h. At address zero, only the entirely
positive path from -h and entirely negative path from h can arrive. Repeated
source multiplication gives

\[
\psi_h^\pm(0)=2^{-h/2}
 [e_k\pm(-1)^{h-1}e_c],\qquad
j_h^\pm(0)=\pm(-1)^{h-1}2^{1-h}.
\tag{19.11}
\]

This proves sharpness at every h, not just checked finite examples. The two
packet supports first meet at zero at event h. Their local energies agree
everywhere through that event, while the currents at zero differ. At event
h+1 the outgoing energy at address 1 differs by
(-1)^(h-1)2^(1-h), with its opposite at -1. This proves energy sharpness
with N=h+1. Energy sufficiency follows from the previous event's outgoing
rho and j, whose pair separations are at most 2(N-1); at N=1 it depends only
on local records. Thus no fixed finite pair-separation cutoff supports
arbitrarily long exact prediction for this whole preparation module. A
known special source or a different compressed encoding is not excluded.

### R19.8 — Odd separation is an invariant invisible sector for these local targets

Define E to retain pairs with x-y even, and O=I-E. Then

\[
E\mathcal T=\mathcal T E,\qquad
\mathcal P O=0.
\tag{19.12}
\]

All local energy and current targets at every future event are therefore
unchanged by discarding the odd-separation sector. The even off-diagonal
sector cannot be discarded in the same way, as (19.6) and (19.11) show.

**Proof.** Each term of (19.2) changes a pair separation by s-t, which is
0 or ±2. Parity is preserved. Diagonal pairs have zero, hence even,
separation. It follows by induction that P T^n O=0 for all n. Both local
targets factor through P, proving the claim. This is a derived, target-faithful
quotient for those local quadratic outputs. It is not claimed to be the
largest invisible subspace, faithful to all possible source targets, or a
quotient of the full UGD algebra. In particular off-diagonal position alone
does not decide whether a record is necessary, curved, or safely erasable.

### R19.9 — The pair calculation lifts to the unmerged cut-word register

For the complete depth-n word source W_n, retain ordered pair tags (w,u)
and the native coefficient v_n(w)v_n(u)^dagger/2^n. Address evaluation gives

\[
\boxed{G_n(x,y)=2^{-n}\!
 \sum_{\substack{w,u\in W_n\\X_n(w)=x,\ X_n(u)=y}}
 v_n(w)v_n(u)^\dagger.}
\tag{19.13}
\]

Thus the pair record is a derived readout of source words, rather than a
replacement primitive. Keeping its ordered tags retains all original words.
Merging equal addresses forgets information; an equal evaluated role state
does not identify two source words or their independent sheet records.

**Proof.** Multiply R18's two finite signed sums for psi_n(x) and
psi_n(y)^dagger and distribute. Their normalizations multiply to 2^-n,
proving (19.13). Every source state is a signed unit role, so every outer
pair coefficient is nonzero. Every original word occurs in the first and
second coordinates of the ordered tags; those tags reconstruct W_n even
when evaluated coefficients cancel after address merging. An additional
native label already attached to a word can be carried on the same tag
without modifying its coefficient or identifying that label with x. No
new law for such a label is introduced. R16 C13 and Publications NC retain
the integer-sheet and proper-seam obligations; R19 does not identify either
with the address shift or impose a four-step sheet closure.

## Derivation and citation audit

| Input or earlier result | Source | Status on R19's proof path |
|---|---|---|
| Signed/refined scalars, role operators and native pairing | [R16 C3–C8, commit 4cc46d1](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native construction and written derivation, consumed in its stated role-module scope |
| Free source census and signed history | [R17, commit f561681](https://github.com/Parveen117/extra-ideas/blob/f561681959f5af93389db6bc77763c6794d76a87/02-relational-response/CUT_HISTORY_MEMORY_R17.md) | Native derived source counts; no detector-frequency premise |
| Address target, signed aggregation, normalization and transport | [R18, commit 7afb4f7](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Explicit native target constructions followed by derived identities; physical selection remains open |
| Compression defect and future observation | [RKF N03/N06](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Native finite composition proofs; R19 instantiates and rederives the required identity |
| Exact hidden-record elimination | [RKF N09](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/02_COMPLETION_AND_MEMORY.md) | Native substitution/induction proof, rederived in (19.4); no classical memory equation imported |
| Sheet, clock and source-to-physics obligations | [Publications NC](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-ugd-propagation/THEOREM.md) | Scope/provenance reference; its declared product carrier, arrow and clock are not inserted as R19's law |

The Publications reference is a research-branch development, not a merged or
externally reviewed release. Its classical-sector CP/NP/SC ancestors and
R18's historical Hadamard-walk comparison are not dependencies of an R19
proof. Native pair multiplication and finite block substitution can coincide
with established matrix methods; no general historical-priority claim is made.
Their mathematical identities have been derived here without importing an
external evolution, metric, probability or Hilbert-space axiom.

The certification categories are: native source, explicit native construction,
written native derivation, comparison-only reference, and open physical
selection. Any imported classical premise or unproved claim on a certified
proof path makes the R19 derivation gate fail. The gate also rejects promoting
an inherited target definition to an unconditional physical-selection result.
This is stricter bookkeeping supported by the proofs, not an automated proof
assistant or an assertion that axioms and definitions cease to exist.

## Exact evidence and next development

The [native adapter](../04-operator-evolution/current_memory_retention.cjs)
uses the unchanged canonical RKF arithmetic and proof replayer. It checks
the pair update against independent response propagation and literal ordered
word pairs, checks the exact memory expansion including initial residue,
reproduces the first return and every tested reset census, and tests the
sharp-horizon witnesses and invariant parity sector. Altered source pins,
algebraic claims and forbidden derivation paths must be rejected.

```bash
python3.12 -B 04-operator-evolution/verify_r19.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```

The [verification report](../04-operator-evolution/R19_VERIFICATION.json)
separates written arbitrary-depth proofs, finite exact checks, source hashes
and formal status. R18/R17/R16 are replayed into temporary outputs. Earlier
evidence stays byte-identical; no broad archive or canonical-master PASS is
claimed by this application certificate.

**Historical replay repair.** R17's local-input pins include its three README
navigation files. R18 subsequently changed navigation, but its nested replay
attempted to run R17 in that updated checkout, causing a source-hash failure.
R19 now assembles the frozen mathematical inputs in a temporary snapshot and
restores only those three historical READMEs from their exact R17 commit.
[The replay-input archive](../04-operator-evolution/R19_HISTORICAL_REPLAY_INPUTS.json)
is checked against both Git blob hashes and the original R17 SHA-256 pins.
This repairs reproduction without loosening the gate, changing mathematical
results, replacing active navigation, or refreshing an old certificate.
Use the R19 command for the complete source-chain replay on the current tree;
the old in-place R18 command does not perform this isolation.

**Next target:** a native interaction/readout must account for the relevant
even-separation pair records, either by retaining them or by deriving its
record-cut operation and the resulting return terms. That supplies a concrete
criterion for testing the physical-selection candidates already in the
framework. A shared event address, physical clock/length process, continuum
limit and gauge normalization still require their own native derivations.
R19 derives current retention and its memory cost; it does not select c,
alpha, a detector mechanism or a new force by importing their classical laws.
