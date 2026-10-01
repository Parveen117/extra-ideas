# R18: cut-address transport before physical metric selection

Research owner: Monty Dabas. Development: 1 October 2026.

**Next target toward c and alpha:** obtain a propagation law from the existing
cut records before selecting a spacetime metric or fitting a physical constant.
The existing Publications NC correction already derives sheet-retaining
propagation for its chosen arrow RZ. R18 develops the different, fixed R17
full-census source and its role-address target; it is not the first native
propagation result. R18 constructs an event address, retains the source signs, and derives an exact
finite transport recurrence. Its two-event transfer satisfies a discrete wave
identity. This is a source-defined mathematical process, not yet a selected
physical interaction. The response aggregation and the metric interpretation
must remain visible in the derivation.

The [R16 foundation](emk-topology-foundation/README.md) and
[R17](CUT_HISTORY_MEMORY_R17.md) remain unchanged. No complex number field,
spacetime continuum, wave equation, stochastic transition parameter, or physical
speed is supplied to the finite construction. In particular the full UGD source
is not replaced by its scalar sector. We work in R16's constructed two-role
module, a particular target of that source, with its specified first seed.

## What is reused and what is constructed

R16 C3–C8 supplies the radial cut field, its positive square-root completion,
the ordered roles e_k,e_c, coefficient pairing B, and

\[
H e_k=e_k,\quad H e_c=-e_c,\qquad
K e_k=e_c,\quad K e_c=e_k,
\]
\[
H^2=K^2=I,\quad HK=-KH,\quad R=KH,\quad R^2=-I.
\tag{18.0}
\]

These maps are read from the verified R16/R17 source packet, not entered as a
new physical ansatz. H and K are self-adjoint isometries for B. R17 supplies
the complete free word census and M=(H+K)/2, M²=I/2. A word is chronological:
v_j=G_j v_(j-1), G_j in {H,K}, v_0=e_k. Its operator product is consequently
G_n ... G_1. Distinct words stay distinct at the source.

The new **address target** records event count n and accumulated role contrast
x. The new **coherent response target** adds signed role responses with the
same address, then uses the positive scalar normalization that preserves the
finite sum of B-energies. These are explicit definitions of an observation
protocol. No assertion says that nature must implement this target or erase
records in this way. A scalar normalization is selected within this protocol;
selection among all possible protocols is a separate question.

## Written results

### R18.1 — Source roles construct an event address and its reachability cone

Every prefix state is ±e_k or ±e_c. Define

\[
s_j=B(v_j,Hv_j)\in\{+1,-1\},\qquad
X_n=\sum_{j=1}^n s_j,\qquad N_n=n.
\tag{18.1}
\]

For this unit seed each cut word has exactly one sign sequence (s_1,...,s_n),
and every sign sequence comes from exactly one cut word. The address set is

\[
\mathcal A_n=\{-n,-n+2,\ldots,n\}.
\]

From any realized prefix a continuation of length m can produce an address
increment d exactly when |d|≤m and d≡m mod 2. The quadratic address expression

\[
q_{\rm addr}(m,d)=m^2-d^2=4n_+n_-\ge0
\tag{18.2}
\]

is thus derived from the two continuation counts n_++n_-=m and n_+-n_-=d.
It is zero on a unidirectional continuation. This is a reachability relation
on record classes, not a spacetime light-cone claim.

**Proof.** H preserves the occupied role and K exchanges it; both preserve
its unit magnitude. If the preceding role sign is s_(j-1), the next tag is H
when s_j=s_(j-1) and K otherwise, starting with s_0=+1. This constructs the
inverse word map and proves bijectivity. A continuation therefore allows every
binary sign string. Counting its two signs proves the address, parity and
reachability statements, and expanding (n_++n_-)²-(n_+-n_-)² proves (18.2).
No spatial lattice was presumed: its integer addresses are images of the
source words under (18.1). The address forgets order, sign and final role in
general; it does not replace the full recognition history.

### R18.2 — Erasing signs gives an exact census transport

The number of words at x in A_n is

\[
C_n(x)=\binom{n}{(n+x)/2},\qquad p_n(x)=2^{-n}C_n(x),
\tag{18.3}
\]

with zero values outside A_n. Here the binomial symbol means the count of
subsets of the indicated size, not an imported probability law. It obeys

\[
p_{n+1}(x)=\tfrac12[p_n(x-1)+p_n(x+1)],\quad
\sum_x xp_n(x)=0,\quad\sum_x x^2p_n(x)=n.
\tag{18.4}
\]

**Proof.** R18.1 identifies words at x with choices of (n+x)/2 positive
positions. Removing the last sign partitions this set into the two preceding
addresses. Dividing the exact counts by 2^(n+1) gives the recurrence. At each
extension, every prefix has one + and one - continuation, so the total first
moment is unchanged and the second moment gains one. Their initial values
are zero. This proves the moments for all n by induction. The finite census
has spreading width squared n. No continuum heat equation or measured random
frequency is needed for these identities.

### R18.3 — Retaining source signs selects a norm-preserving transport normalization

For each address define the unnormalized response

\[
A_n(x)=\sum_{w:\,|w|=n,\,X_n(w)=x} v_n(w),\quad
A_0(x)=\mathbf1_{x=0}e_k.
\]

Let Pi_+=(I+H)/2 and Pi_-=(I-H)/2 be the native role projectors. On finite
address fields f define (Tf)(x)=f(x-1) and

\[
S=\Pi_+T+\Pi_-T^{-1},\qquad W=S(H+K).
\tag{18.5}
\]

Then A_(n+1)=W A_n. For the finite address pairing
<f,g>_A=Σ_x B(f(x),g(x)), S is an invertible isometry and

\[
W^\dagger W=WW^\dagger=2I.
\]

There is exactly one positive scalar a for which aW preserves this pairing:
a=1/√2. Consequently

\[
\psi_n=2^{-n/2}A_n,\quad
\psi_{n+1}=U\psi_n,\quad U=S C,\quad C=(H+K)/\sqrt2,
\quad\sum_x B(\psi_n(x),\psi_n(x))=1.
\tag{18.6}
\]

**Proof.** From an occupied role, each next tag contributes H v or K v; its
output role determines whether its new address is x+1 or x-1. Projecting
H+K onto the two output roles gives (18.5). Linear extension gives the same
law on arbitrary finite address fields. Pi_± are orthogonal idempotents and
Pi_++Pi_-=I; shifting each projected component preserves the finite energy
sum. The inverse is Pi_+T^-1+Pi_-T. Anticommutation gives
(H+K)†(H+K)=(H+K)(H+K)†=2I. Thus aW is an isometry precisely when 2a²=1.
R16 C8 constructs the unique positive square root used in (18.6); an ordinary
complex scalar is not a premise. Induction gives the response and energy laws.
The pairing of different addresses by a direct finite sum is part of this
target's definition, not an independently proved physical superposition or
measurement rule. No infinite Hilbert-space representation is used.

### R18.4 — The exact current retains interference discarded by counting

Write ψ_n(x)=a e_k+b e_c and define its outgoing response energies

\[
r_n(x)=(a+b)^2/2,\quad l_n(x)=(a-b)^2/2,
\quad\rho_n(x)=a^2+b^2,\quad j_n(x)=2ab.
\tag{18.7}
\]

Then r+l=ρ, r-l=j, ρ≥|j|, and

\[
\rho_{n+1}(x)=r_n(x-1)+l_n(x+1).
\]

With J_n(x+1/2)=r_n(x)-l_n(x+1), there is the exact finite conservation law

\[
\rho_{n+1}(x)-\rho_n(x)
=J_n(x-1/2)-J_n(x+1/2).
\tag{18.8}
\]

This normalized native energy is different from record-count content. At
n=3, addresses (-3,-1,1,3) have

\[
p_3=(1,3,3,1)/8,\qquad \rho_3=(1,1,5,1)/8.
\tag{18.9}
\]

**Proof.** C sends (a,b) to ((a+b)/√2,(a-b)/√2), after which S shifts the two
roles in opposite directions. Their squares give (18.7). Nonnegative r,l
give ρ≥|j|; subtracting the outgoing energy at x gives (18.8). Direct source
expansion gives A_3(-3)=e_c, A_3(-1)=-e_k, A_3(1)=2e_k+e_c and A_3(3)=e_k.
R18.2 gives the four census counts independently. This proves (18.9) without
sampling. ρ is a defined quadratic response; interpreting it as detector
probability would be an additional physical claim.

### R18.5 — Two events obey an exact rational discrete wave identity

Put V=U²=W²/2. In the finite Laurent algebra let z denote only the address
shift T. It commutes with role matrices because they act on different indices;
z is not iota and is not a supplied complex phase. In the ordered role basis,

\[
W(z)=\begin{pmatrix}z&z\\z^{-1}&-z^{-1}\end{pmatrix},\qquad
V(z)=\tfrac12\begin{pmatrix}z^2+1&z^2-1\\1-z^{-2}&1+z^{-2}\end{pmatrix}.
\tag{18.10}
\]

The exact polynomial identity is

\[
V^2-\tfrac12(z^2+2+z^{-2})V+I=0.
\tag{18.11}
\]

Every source response therefore satisfies, for n≥2,

\[
\boxed{\psi_{n+2}-2\psi_n+\psi_{n-2}
=\tfrac12(T^2-2I+T^{-2})\psi_n.}
\tag{18.12}
\]

**Proof.** Insert the source-derived H,K and projectors in (18.5). Multiplying
W by itself gives (18.10). Direct multiplication of the four Laurent entries
gives (18.11); equivalently its determinant is one and its trace is the displayed
scalar, with the corresponding two-by-two identity verified by expansion.
Apply it to ψ_(n-2) and use ψ_n=Vψ_(n-2) to obtain (18.12). In entirely integral
unnormalized responses the same statement is

\[
A_{n+4}=(T^2+2I+T^{-2})A_{n+2}-4A_n.
\tag{18.13}
\]

The canonical native Laurent runtime checks (18.11) coefficient by coefficient,
and an independent word census checks (18.13). The wave stencil is a consequence
of the signed update. Its unrestricted second-order solutions can have extra
initial data; the stencil alone is not equivalent to U or a physical field law.

### R18.6 — The first transport symbol is selected within this address protocol

Let Q_±=(I±K)/2. Native multiplication gives C H C=K and hence

\[
V=(\Pi_+T+\Pi_-T^{-1})(Q_+T+Q_-T^{-1}).
\tag{18.14}
\]

Introduce a formal bookkeeping symbol d, with d³=0, and replace
T by 1-εd+ε²d²/2 and T^-1 by 1+εd+ε²d²/2. This is a finite jet substitution,
not a spacetime or convergence assumption. It gives

\[
V=I-\epsilon(H+K)d+\epsilon^2(I+HK)d^2\pmod{d^3}.
\tag{18.15}
\]

Per two source events, the first transport coefficient is therefore exactly
the R17 census operator M=(H+K)/2. Its scalar characteristic polynomial is

\[
\det(\omega I-kM)=\omega^2-\tfrac12 k^2.
\tag{18.16}
\]

This selects a (1+1)-dimensional *candidate* co-metric diag(1,-1/2) for a
smooth first-order reading of this protocol; its inverse is diag(1,-2).
It is not the address reachability form diag(1,-1) of R18.1.

**Proof.** C²=I and C H C=K follow by expanding H²=K²=I and HK=-KH.
Conjugating the first factor in S C S C proves (18.14). Its two factors have
jets I-εHd+ε²d²/2 and I-εKd+ε²d²/2; multiplying them proves (18.15). Division
of the first coefficient by the block length two gives M. Its trace is zero
and M²=I/2, or directly its entries are [[1,1],[1,-1]]/2, which proves (18.16).
The quadratic forms have the stated signatures and inverses by inspection
of their derived diagonal coefficients. If smooth t,x are later assigned,
the formal leading equation is ∂_tψ+M∂_xψ=0 and its scalar square is
∂_t²ψ-(1/2)∂_x²ψ=0. **No continuum-limit convergence theorem, arbitrary-mode
dispersion law, or observed spacetime identification is claimed here.** The
certificate proves the finite symbol and stencil, not those further statements.

### R18.7 — The exact front differs from the slow-variation characteristic slope

For every n≥1, both extreme addresses are nonzero and

\[
\rho_n(n)=\rho_n(-n)=2^{-n},\qquad
\operatorname{supp}\psi_n\subseteq[-n,n].
\tag{18.17}
\]

The exact support front is one address increment per source event. The formal
characteristic slope of (18.16) is 1/√2. Their ratio is 1/√2; these are different
notions of propagation in this construction.

**Proof.** Repeated shifts in S change an address by at most one each event.
Only H^n has every output role e_k, so A_n(n)=e_k. Only the chronological
word K H^(n-1) has every output role e_c, so A_n(-n)=(-1)^(n-1)e_c. Squaring
and applying the normalization gives (18.17); no cancellation is possible at
either endpoint because each has one word. The slopes of (18.16) have squared
magnitude 1/2. Thus a coefficient read from the long-scale symbol must not be
labelled the exact microscopic front, and neither is yet physical c.

### R18.8 — A minimal local response distinguishes directed transport

The local quadratic pair (ρ,j)=(B(v,v),B(v,Kv)) determines both outgoing
energies by r=(ρ+j)/2, l=(ρ-j)/2. The scalar ρ alone cannot determine them.
Within linear spans of quadratic forms on the two-role module the minimum
number of independent response scalars for both outgoing energies is two.

**Proof.** The formulas follow from R18.4. States e_k+e_c and e_k-e_c both
have ρ=2, but their currents are +2 and -2 and their outgoing energies are
(2,0) and (0,2). Moreover the forms r and l are independent: evaluating a
linear relation ar+bl=0 on those two states separately forces a=b=0. A one-row
linear quadratic-response bank cannot span both independent targets. The pair
(ρ,j) has rank two and suffices. This is a target-specific minimum, not a claim
about arbitrary nonlinear encodings of two numbers into one. R17's complete
P/PH state observer also suffices to compute the pair. Directed energy flow
therefore needs retained cross-role information; local energy or plain record
count alone does not supply it.

### R18.9 — Physical calibration remains outside the derived address law

Let an eventual rod/clock adapter assign x_phys=L_0 x and t_phys=T_0 n with
positive L_0,T_0. R18's exact support speed and candidate characteristic speed
would then be

\[
c_{\rm front}=L_0/T_0,\qquad
c_{\rm symbol}=L_0/(\sqrt2 T_0).
\tag{18.18}
\]

All preceding finite identities remain unchanged when L_0/T_0 is changed.
Consequently this packet selects neither a dimensional c nor a universal
physical metric. Its ratio is a property of its defined transport protocol,
not yet an experimentally assigned prediction. It supplies one address
direction, not a derivation of three spatial directions or a common cone for
different physical fields.

**Proof.** A displacement of one address per event becomes L_0/T_0 by the
adapter definition; the other slope is multiplied by the same factor. Neither
constant occurs in (18.0)–(18.17), so substituting any other positive adapter
preserves those identities. This is a boundary of the present specified source
target, not a universal impossibility theorem for the generalized recognition
framework. Deriving a physical rod/clock process or a stronger interaction
selection law would add evidence not present in this packet.

## Existing results are referenced with their premise status

All GitHub links below are commit-pinned. The Publications developments cited
here are on the `research/uncut-cut-measurement-2026-09-25` branch (PR #5),
not its main branch or an independently reviewed release. Exact local bytes are bound by
[R18_SOURCE_PINS.json](../04-operator-evolution/R18_SOURCE_PINS.json).

| Result used or compared | Source and evidence | Premise status in this development |
|---|---|---|
| Signed roles, iota, coefficient pairing, radial completion | [R16 C3–C8](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex); written proofs plus finite exact certificate | Consumed native construction in its stated two-role scope; neither full UGD reduction nor a physical module-selection theorem |
| Full word census and M²=I/2 | [R17.1–R17.3](https://github.com/Parveen117/extra-ideas/blob/f561681959f5af93389db6bc77763c6794d76a87/02-relational-response/CUT_HISTORY_MEMORY_R17.md) | Consumed native count derivation; physical equal-frequency law not asserted |
| Exact operator, Laurent and proof replay routines | [canonical RKF operator home](https://github.com/Parveen117/Recognition-Kernel-Framework/tree/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation) | Unchanged native finite arithmetic; a computational checker, not a proof assistant |
| Full native propagation and retained sheets | [Publications NC-1–NC-7](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-ugd-propagation/THEOREM.md) | Native normal-form, returning-memory and clock identities under an explicit product carrier and chosen arrow RZ. The sheet, proper seam, UGD numeral operations and clock are distinct types; no universal physical propagator is selected |
| Positive response Lorentz cone | [Publications CP-1–CP-2](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/ugd-kahler-propagation/THERMO_GAUGE_COMPLETION.md) | Comparison only: selects a two-block J-compatible self-adjoint sector; spacetime interpretation and extra relative phase are not derived there |
| Clifford wave symbol | [Publications NP-1–NP-2](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/ugd-kahler-propagation/NATIVE_PROPAGATION_ACTION.md) | Comparison only: admitted tensor representation, smooth continuum and propagation coefficients; not imported as R18's evolution law |
| Common vacuum speed and calibration boundary | [Publications SC-1–SC-6](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/ugd-kahler-propagation/SPEED_AND_CALIBRATION.md) | Comparison only: supplied common-coframe Einstein–Maxwell–matter action and physical adapters; no primitive-only c selection |
| Gauge coupling normalization | Publications CP-3–CP-4, same pinned file as CP-1 | Supplied Z,g,e and gauge kinetic law remain explicit; R18 does not use them to manufacture alpha |

The normalized matrix C and conditional shift coincide with the established
real Hadamard walk. Historical comparison: A. Ambainis, E. Bach, A. Nayak,
A. Vishwanath and J. Watrous, [One-Dimensional Quantum Walks](https://cs.uwaterloo.ca/~watrous/Papers/OneDimensionalQuantumWalks.pdf),
STOC 2001. **External/classical literature, attribution only:** its complex
Hilbert/probability setup and continuum/asymptotic results are not proof
dependencies here. Equations (18.1)–(18.17) are derived above from the pinned
cut construction. Reproducing a known walk in native notation establishes
lineage and source closure for this target, not historical novelty, quantum
theory from no premises, or experimental validation.

## Certification and the next gate

[Exact checks](../04-operator-evolution/cut_transport.cjs) use the canonical
native arithmetic and Laurent runtime. They compare literal source-word
enumeration with the derived response recurrence, verify the full polynomial
identity and its finite jets, check conservation and the minimal quadratic
observer, and reject incorrect normalizations and metric interpretations.
[The verifier](../04-operator-evolution/verify_r18.py) first checks source
hashes and reproduces R17/R16 in temporary outputs. Existing evidence stays
byte-identical. The nine written proofs have separately hashed sections in
the [claim ledger](../04-operator-evolution/R18_CLAIM_LEDGER.json).

```bash
python3.12 -B 04-operator-evolution/verify_r18.py \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --publications-root /path/to/Publications
```

The full word register is retained upstream of the address aggregation. R18
never asserts that an equal address or a returned two-role state closes that
register. Its address shift T is bilateral and is not the canonical proper
depth seam with TS=1 but ST≠1, nor is x identified with a UGD sheet. NC-2
already proves why those identifications fail. R18 does not discharge the
full source-to-physics interface simply by exhibiting a two-role stencil.

This is written mathematics with a scoped exact computational certificate.
It is not formal proof-assistant verification or external review. Passing the
checks does not turn the address protocol into a physical law.

**Next gate:** derive, from an interaction/readout already available in the
native source, whether the cross-role current of R18.8 is retained and how
interacting targets share an operational event address. The interface must
retain NC's returning term, independent integer sheets, nonaffine clock term
and matched Smriti error control. Then establish an
actual controlled large-scale limit and its dimensions before assigning a
physical metric. The exact source front and the candidate symbol must remain
separate until that analysis chooses the relevant propagation notion. After
the common propagation/clock structure is selected, a dimensionless coupling
target such as alpha requires a source-derived normalization of phase charge
relative to the gauge response. No target physical constant is fitted here.
