# R20: native record writing, readout and returning interaction

Research owner: Monty Dabas. Development: 1 October 2026.

R19 derived the memory needed by current transport and calculated the effect
of cutting it. R20 constructs an actual reversible record-writing arrow from
the same cut-role source. A fresh slot at every event realizes R19's repeated
cut on its stated input class. Reusing one slot instead preserves every local
energy and current of the unrecorded single-origin propagation. This is an
exact interaction/readout mechanism, with no supplied coupling strength,
random reset, measurement axiom or classical evolution equation.

The distinction between **constructing an available arrow** and **selecting a
physical law** remains essential. Slot allocation, readout target and initial
preparation are explicit constructions here. Their consequences are proved;
their unique physical selection is not asserted. The certification rule does
not permit changing a definition's label to make that missing result disappear.

## Source, carrier and notation

Use unchanged R16 C1, C3–C8 and C13, R18.1–R18.4 and R19.1–R19.9. On the native
marked roles write e_0=e_k, e_1=e_c, H e_a=(-1)^a e_a and K e_a=e_(1-a).
Thus H²=K²=I, HK=-KH. These maps were constructed on free cut roles before
their matrices or represented iota. Define the already derived arrows

\[
\Pi_a=(I+(-1)^a H)/2,\quad F=H+K,\quad C=F/\sqrt2,
\quad L_a=\Pi_a F,\quad s_a=(-1)^a.
\]

The inherited address update is

\[
(U\psi)(x)=2^{-1/2}\sum_{a=0}^1L_a\psi(x-s_a).
\tag{20.0}
\]

R18 proves its pairing preservation and inverse. All fields have finite
support at each finite event. Coefficients are the native signed/refined
counts and, where necessary, their R16 C8 completion. The checker's rational
arrays encode these derived operations; an ordinary complex field is not a
new premise.

A record slot is another copy of the same **marked-role ledger**. Multiple
slots are constructed as the free ledger on finite tuples of those marks.
For system and record labels, form the free coefficients A_(x,a,m), with

\[
B_{\rm joint}(A,D)=\sum_{x,a,m}A_{x,a,m}D_{x,a,m}.
\tag{20.1}
\]

This repeats C6's matching-mark construction on tuples. The notation
v tensor r abbreviates the coefficient array v_a r_m; expansion gives
B(v tensor r,w tensor t)=B(v,w)B(r,t). No independently postulated Hilbert
space, tensor rule for nature, partial trace or Born probability is used.
The enlarged finite target is a construction, not a claim about the entire
generalized UGD carrier. Source words and independent sheet labels are not
identified merely because their evaluated coefficients agree.

## Written results

### R20.1 — The source constructs a reversible, nontrivial record gate

On a system role and one record role define

\[
\boxed{W=\Pi_0\mathbin{\mathrm{tensor}} I+
              \Pi_1\mathbin{\mathrm{tensor}} K.}
\tag{20.2}
\]

It leaves (0,b) fixed and sends (1,b) to (1,1-b). Consequently W²=I,
W^dagger=W, and W preserves (20.1). This is the unique linear extension of
that specified marked-role rule, not the only interaction allowed by the
source. Its action does not factor into separate arrows on the two ledgers.

**Proof.** The rule permutes the four ordered marks by one exchange. Applying
it twice restores every mark, and equality of marks is preserved. Alternatively,
Pi_a Pi_b=delta_ab Pi_a and K²=I give both identities by multiplication.
For the native prepared input C e_0 tensor e_0, the output is
(e_0 tensor e_0+e_1 tensor e_1)/sqrt(2). If it were v tensor r, its four
coefficients would obey A_00 A_11=A_01 A_10, whereas here they are 1/2 and 0.
Thus the output does not factor. Also W(C tensor I) and (C tensor I)W differ
on e_0 tensor e_0. The native preparation and writing operations genuinely
fail to commute. The map writes a role mark; it is not a copying map for an
arbitrary coefficient vector or its unknown full history.

### R20.2 — The retained cross coefficient is the record's native current

Define the record-unresolved pair readout by matching the record marks:

\[
G^\sharp_{ab}(x,y)=\sum_m A_{x,a,m}A_{y,b,m}.
\tag{20.3}
\]

This is a specified finite response target. Its local energy is
rho(x)=sum_a G^sharp_aa(x,x), and j(x)=G^sharp_01(x,x)+G^sharp_10(x,x).
For a nonzero record preparation mu=u e_0+v e_1, use its normalized value
mu/sqrt(E), E=u²+v². Writing one role mark multiplies the system pair block
G_ab(x,y) by 1 for a=b and by

\[
\boxed{\kappa=\frac{B(\mu,K\mu)}{B(\mu,\mu)}
              =\frac{2uv}{u^2+v^2}}
\tag{20.4}
\]

for a!=b. In particular -1<=kappa<=1. The source preparations e_0,
C e_0 and H C e_0 give 0, +1 and -1, respectively. Kappa is a calculated
record response, not a noise strength introduced as an extra parameter.

**Proof.** After W the record attached to role a is K^a mu/sqrt(E).
Substitution into (20.3) yields G_ab B(K^a mu,K^b mu)/E. The pairing-preserving
involution K gives 1 for equal roles and (20.4) otherwise. Positivity of
(u-v)² and (u+v)² gives the two bounds; E>0 follows from C6. Direct evaluation
gives the three displayed preparations, with the common sqrt(2) cancelling.
For several fresh slots, expanding the finite matching sum factors the
overlap into the product of their individual overlaps. This is a product of
derived coefficients, not a probabilistic independence postulate. A negative
overlap retains a reversed sign; it is not a negative probability.

### R20.3 — Fresh native records realize the repeated address cut on its exact domain

After each U step, apply W to a fresh e_0 slot; retain all earlier slots but
do not couple them again. Define the role-block cut
\(\mathcal J G=\sum_a\Pi_a G\Pi_a\), acting separately on each (x,y) block.
The exact record-unresolved update on any finite pair array is

\[
\boxed{G^\sharp_{n+1}=\mathcal J\mathcal T G^\sharp_n,}
\tag{20.5}
\]

where T is R19's derived pair propagation. For address-diagonal inputs
P G=G, it equals P T G and is again address- and role-diagonal. Therefore
the single-origin seed realizes R19's repeated-cut sequence (P T)^n G_0.
No physical destruction or random reset of the old record is required.

**Proof.** In each product contributing to (20.3), the fresh slot gives
B(K^a e_0,K^b e_0)=delta_ab. Summing the unchanged old record labels produces
the old G^sharp, proving (20.5), also for correlated old records. For an
address-diagonal input, the summand L_a G(x-s_a,y-s_b)L_b^dagger can survive
only if x-s_a=y-s_b. The role cut additionally requires a=b, whence x=y.
Conversely, x=y in T of such an input forces a=b. Thus J T G=P T G, with
only a single output role per surviving branch. Induction proves the claim.
For general off-address input, same-role pairs at different addresses survive
J. For example, the pair field built from e_0 at -1 and +1 has a nonzero
e_0 e_0^dagger block between 0 and 2 after (20.5). Hence J T is not P T
on the full preparation module. This scope condition is part of the theorem.

### R20.4 — Full fresh records retain source words and derive the count response

For psi_0(0)=e_0, let d=(a_1,...,a_n) record the successive output roles,
x(d)=sum_i s_(a_i), and epsilon(d)=+1 or -1 be the product of the native
F transition signs. Fresh writing gives exactly

\[
\Psi_n^{\rm fresh}=2^{-n/2}\sum_{d\in\{0,1\}^n}
 \epsilon(d)e_{x(d),a_n}\mathbin{\mathrm{tensor}}e_d\quad(n>=1).
\tag{20.6}
\]

R18's bijection between output directions and the original H/K words retains
all 2^n source tags. Let N_a(n,x) count the histories ending at (x,a). The
pair readout is diagonal on system ports, with diagonal N_a(n,x)/2^n.
Thus j_n(x)=0 and rho_n(x)=binom(n,(n+x)/2)/2^n, using zero outside the
native reachable addresses or at a nonintegral binomial index. More precisely,

\[
N_0(n,x)=\binom{n-1}{(n+x)/2-1},\qquad
N_1(n,x)=\binom{n-1}{(n+x)/2}.
\tag{20.7}
\]

**Proof.** In the marked-role basis every F transition has coefficient +1
except 1 to 1, which has -1. U appends one direction with coefficient
1/sqrt(2), and W appends precisely that output mark to the fresh slot.
Induction gives (20.6). Distinct d have zero matching pairing, so (20.3)
removes every cross-history term and squares each sign to one. To finish
with role 0 at x, fix the last positive step and count the remaining positive
steps among n-1 slots; this gives N_0. Fixing a last negative step gives N_1.
Their sum counts all words with the prescribed number of positive steps.
This proves (20.7) and the stated responses by finite counting. No binomial
probability distribution is assumed. Equal evaluated operator returns do
not equate different record words.

### R20.5 — Reusing one record slot preserves the complete local coherent response

Start from any native role v at the single address 0 and one normalized
record mu. Use the same slot in W after every U. At each reachable address
put m_n(x)=(n-x)/2, the number of negative steps. The joint field is

\[
\boxed{\Psi_n^{\rm reuse}(x)=\psi_n(x)\mathbin{\mathrm{tensor}}
                         K^{m_n(x)}\mu,\qquad\psi_n=U^n\psi_0.}
\tag{20.8}
\]

Consequently G^sharp(x,y)=G_n(x,y) when m_n(x)-m_n(y) is even and
kappa G_n(x,y) when it is odd, with kappa from (20.4). In particular,
every local rho_n and j_n is exactly the unrecorded R18 response at all n.
For mu=e_0, the surviving address pairs have x-y divisible by four.

**Proof.** A positive branch leaves m unchanged and the record unexchanged.
A negative branch increases m by one and applies K. Thus all paths arriving
at a fixed x have the same record K^m mu, independent of their order and
last role. Factor it out of their signed sum; induction gives (20.8), including
any cancellations of that sum. Matching the two record factors proves the
pair formula. At x=y their pairing is one, proving the local equality.
This statement uses a single initial address; a superposition of distinct
initial addresses need not have a common m_n(x), and is not silently included.
The shared slot is a returning record, not a fresh slot with the same name.

### R20.6 — The first return distinguishes record allocation and defeats local closure

For the e_0 seed and e_0 blanks, fresh and reused protocols have the same
record-unresolved pair field at event one. At event two their energies still
agree, but at x=0 the reused current is 1/2 and the fresh current is zero.
At event three their energies, ordered by x=-3,-1,1,3, are

\[
\rho^{\rm reuse}_3=(1,1,5,1)/8,\qquad
\rho^{\rm fresh}_3=(1,3,3,1)/8.
\tag{20.9}
\]

Even with a fixed subsequent gate on a specified slot, G^sharp alone is not
a closed state for every joint preparation.

**Proof.** Event one is
\((e_{1,0}e_0+e_{-1,1}e_1)/\sqrt2\). Its readout is address-diagonal in
both protocols. At event two the two paths returning to x=0 both give
record e_1 when the slot is reused, so their local role cross term survives.
R19.4, or multiplying F twice, gives j(0)=1/2. Fresh records 01 and 10 are
distinct, giving zero. R20.4 and R20.5 then give (20.9).

For the stronger fixed-gate witness append an untouched second slot c.
Prepare Omega=(e_(1,0)e_0^b+e_(-1,1)e_1^b)e_0^c/sqrt(2), or
Xi=(e_(1,0)e_0^b e_0^c+e_(-1,1)e_1^b e_1^c)/sqrt(2).
Their entire G^sharp fields agree, since b already distinguishes the two
branches. Apply exactly the same next U and W on b, leaving c untouched.
At x=0, b becomes e_1 on both returning paths in both experiments. In Omega,
c also agrees, giving current 1/2; in Xi, c still distinguishes the paths,
giving zero. Thus equal present system pair readouts do not determine the
next current even under this fixed next interaction. The missing information
is the joint record, not an externally assigned random force.

### R20.7 — Address records realize equality cuts and identify a harmless cut

For a finite address set S choose a code f:S -> {0,1}^b, built from predicates
on the existing integer address. On b blank slots exchange slot i by K exactly
when f_i(x)=1, leaving the system role unchanged. Extend the predicates to
other addresses by any stated values when a larger domain is needed. The
resulting arrow D_f is a pairing-preserving involution and its readout is

\[
\boxed{G^\sharp(x,y)=\mathbf1_{f(x)=f(y)}G(x,y).}
\tag{20.10}
\]

An injective code realizes R19's address cut P on fields supported in S. A
one-slot parity code realizes the even-separation cut E of R19.8; applying
this cut at any finite set of events leaves all future local rho and j
unchanged, for every finite initial pair field. A final address record is
not generally equivalent to fresh role records after every event.

**Proof.** Each controlled exchange is a permutation and distinct slots
commute. On blank slots the record becomes e_f(x); equality of these labels
gives (20.10) by coefficient matching. An injective code gives delta_xy.
For parity, equality is x-y even. R19.8 proves E T=T E and P E=P, so any
sequence of such cuts reduces to E T^n for local readouts and changes none
of them. Alternatively, each pair separation changes only by 0 or +/-2,
so the discarded sector never reaches x=y. Finally, recording the address
only at event two keeps j(0)=1/2 inside its local block, whereas R20.4 has
j(0)=0 after fresh recording at both events. Thus timing and recorded target
cannot be omitted from an interaction specification.

### R20.8 — Record-only reversible processing cannot erase a distinction

Every pairing-preserving arrow O acting only on the record leaves the entire
G^sharp unchanged. Coherence restoration therefore requires a coupled action
or a different readout, rather than merely renaming the record labels.
For the reused protocol, an address-controlled K^m_n(x) at fixed n uncomputes
its record exactly. For n fresh slots, let V_n denote the complete ordered
writing evolution on those slots. The available coupled decoder

\[
\boxed{D_n=(U^n\mathbin{\mathrm{tensor}} I)V_n^{-1}}
\tag{20.11}
\]

maps V_n(psi_0 tensor e_(0...0)) to U^n psi_0 tensor e_(0...0).

**Proof.** Regard each system port's record coefficients as a vector r_p.
The unresolved pair entry is B(r_p,r_q). After O it is B(O r_p,O r_q),
which equals the original entry by its stated pairing preservation. This
proves the first claim without a measurement theorem. In (20.8), the
address-controlled exchange squares K^m to I and so restores the blank;
it acts jointly on address and record, outside the record-only restriction.
For (20.11), each W has its derived inverse W and U has the inverse supplied
by R18. Reverse their order to construct V_n^-1, then apply U n times to
the system. Cancellation proves (20.11) for every admissible initial system
field with those blank slots. The protocol reverses and replays the dynamics;
it is not free disposal of recorded information or an instantaneous physical
eraser. No thermodynamic cost or availability in an experiment is asserted.

### R20.9 — Full history storage and one-event readout have different exact costs

In the fresh-slot construction at depth n, all 2^n history labels are
mutually orthogonal. Retaining their full signed ledger by an injective
linear record map requires rank at least 2^n; n two-role slots attain it.
This is not a bound on arbitrary nonlinear encodings of individual words.
By contrast, for n>=1
the single-seed final unresolved pair field of R20.4 has a joint coefficient
realization with **exactly 2n record directions**, and no smaller rank can
realize it. This second minimum concerns one final full pair target, not
causal online storage or retained history identity.

**Proof.** Pairing a putative linear relation among distinct history marks
with each individual mark forces its coefficient to vanish. Thus they are
linearly independent. A ledger of b two-role slots has 2^b basis marks by
successive concatenation, proving the first bound and its attainment.

At depth n there are n+1 reachable addresses. Both roles have positive
counts at each of the n-1 interior addresses; at each endpoint exactly one
role has positive count. Formula (20.7) proves this also at n=1. Hence exactly
2n system ports have positive diagonal entries, and all off-port entries
vanish. In any joint realization, the corresponding record vectors r_p have
B(r_p,r_q)=0 for p!=q and B(r_p,r_p)>0. Pairing a linear relation with r_q
again forces each coefficient to zero. They require at least 2n independent
record directions. Conversely, construct 2n distinct record marks t_p and
take r_p=sqrt(N_p/2^n)t_p. These square roots exist by R16 C8; matching marks
gives precisely the required pair field and proves attainment. This is a
finite construction proof, not an imported purification or rank theorem.
It does not preserve all 2^n histories and does not specify how to compress
a continuing experiment. A minimum rank is not a derived entropy, time unit,
metric scale or physical constant.

## Evidence and premise-labelled lineage

The [derivation ledger](../04-operator-evolution/R20_DERIVATION_LEDGER.json)
lists the nine written proofs and their exact dependencies. Only native
sources, explicit native constructions and derived results occur on proof
paths. The origin checker rejects classical/admitted/comparison premises,
cycles, unpinned sources and attempts to promote definitions to physical
selection. It checks declared metadata and bound proof text, not the semantic
truth of arbitrary prose. This is not proof-assistant formalization.

The [adapter](../04-operator-evolution/native_record_interaction.cjs) uses the
unchanged RKF arithmetic and symbolic replayer. It compares literal tagged
H/K histories, full joint response propagation and independently evolved
pair fields; checks fresh and reused allocations, fixed-gate nonclosure,
coarse address cuts, inverse protocols and exact minimum-rank witnesses;
and rejects false alternatives. The [wrapper](../04-operator-evolution/verify_r20.py)
replays the frozen R19/R18/R17/R16 chain into temporary outputs and checks
all 166 pre-existing non-navigation files byte for byte. Written all-depth
proofs, finite exact regressions, source provenance, formal verification and
physical validation remain separate evidence categories.

| Pinned source | Status in this derivation |
| --- | --- |
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native construction and written proof: free marks, refinements, H/K, matching pairing, completion and distinct histories. The finite marked-role carrier is explicitly scoped. |
| [R18.1–R18.4](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Native derivation within the constructed address and aggregation target; physical target selection not supplied. |
| [R19.1–R19.9](https://github.com/Parveen117/extra-ideas/blob/f46316e4407a65bca1e07f43baf32a419048f67e/02-relational-response/CURRENT_MEMORY_RETENTION_R19.md) | Native pair update, return, count cut and invisible parity sector; inherited and replayed. |
| [RKF canonical algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Native word equality/replay, used for scoped identities; no claim that every upstream certificate passes. |
| [Publications U27–U29](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/NATIVE_INTERACTION_AND_LAWFUL_CUTS.md) | Comparison only: supplied real carrier, graph, independently addressed controls and nonzero Cayley parameters. Those declared inputs are not source-selected premises here. Its proof certificate is not rerun by R20. |
| [Publications U24–U26](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/uncut-cut-measurement/NATIVE_GRADING_AND_SEAM_MEMORY.md) | Comparison only: supplied graph/port carrier and readout; warns that local native grading need not descend through an assembly. R20 makes no common global EMK-grading claim for its new tuple carrier. |
| [Publications NC-1–NC-7](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-ugd-propagation/THEOREM.md) | Comparison only: native sheet/return distinction retained; the chosen product arrow and clock are declared inputs, not physical selections imported here. |

In ordinary coordinates (20.2) has the familiar controlled-exchange form,
and (20.3) is a Gram contraction. These comparisons acknowledge established
mathematical forms; they are not imported premises or claims of historical
priority. The native derivations above are the certified content.

## Next derivation boundary

The record-writing arrow and its cut/return mechanism are now explicit.
Same-source admissible constructions give distinguishable responses depending
on which slot is available again. Therefore H/K identities and pairing
preservation alone do not select a memory-allocation law. This is an explicit
separation within that stated class, not a universal impossibility theorem
about the full Recognition framework.

Next derive the admissible continuation and identification of record slots
from retained cut-history relations, including which joint records can meet
again. That is the missing input to a source-selected operational interaction
and common event/metric target. Do not insert a classical environment, time
scale, locality rule or coupling to choose between the constructed protocols.
No value of physical c or alpha follows from this packet.
