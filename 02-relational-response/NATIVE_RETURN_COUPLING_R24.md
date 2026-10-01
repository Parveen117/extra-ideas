# R24: native first return, fixed response coefficients and record coupling

Research owner: Monty Dabas. Development: 1 October 2026.

R24 takes the next concrete step toward a source-fixed dimensionless
interaction. It eliminates the entire first-return tail of R18's native
transport, without supplying a damping coefficient, probability law,
continuum equation or classical return theorem. The coherent return map has
the exact positive scale

\[
g=\sqrt2-1.
\]

Keeping arrival marks instead gives a different, exactly specified energy
coefficient \(\eta\). It determines a return-controlled realization of R20's
record gate. Both coefficients have native derivations and certified rational
error bounds. Neither is identified with the electromagnetic fine-structure
constant. In particular, normalizing a returned record cancels its overall
scale: substituting \(g\) for R20's normalized record overlap is incorrect.

The transport, origin-return cut and record interface below are explicit
native constructions. Their consequences are proved, but nature's unique
choice of these constructions is not. A second available native memory
protocol gives a different coupling, proving that this selection cannot be
silently inferred from the present source. This is a scoped mathematical
development, not an unconditional physical-constant claim.

## Native source and declared target

Use R16 C1, C3–C8, C13; R18.1–R18.4; and R20.1–R20.4, R20.7. On the native
marked roles, \(H^2=K^2=I\), \(HK=-KH\), \(R=KH\), \(R^2=-I\), and
\(R^\dagger=-R\). Set

\[
F=H+K,\quad \Pi_a=(I+(-1)^aH)/2,\quad L_a=\Pi_aF,
\qquad (U\psi)(x)=2^{-1/2}\sum_{a=0}^1L_a\psi(x-(-1)^a).
\tag{24.0}
\]

These are the already derived arrows, not an ordinary complex starting
field. Work on the real signed/refined coefficient sector of the native
finite ledger and its C8 scalar completion. This sector does not replace
generalized UGD or identify distinct full source words. Exact check arrays
are a chart for the derived operations.

Embed a role at address zero by \(E\); let \(P=EE^\dagger\) be equality with
that address and \(Q=I-P\). At each step, keep the returned branch with a
fresh arrival mark and continue only the live \(Q\) branch. Matching marks
make different arrivals orthogonal. This construction uses R20's address
record rule and a declared live/returned routing; it is not a physical
measurement postulate. Within one arrival branch, paths still add with their
native signs. Full source words remain available as provenance, not as an
extra orthogonal record imposed on that coherent sum.

Define

\[
S_n=(QU)^nE,\quad A_n=E^\dagger U(QU)^{n-1}E\ (n\ge1),
\quad V_n=E^\dagger U^nE,\quad S_0=E,\quad V_0=I.
\tag{24.1}
\]

The symbol \(z\) below is a formal event-count marker. Evaluating its
convergent series at one does not assign a physical frequency or add a
decay parameter. Coherently summing arrival coefficients and retaining
orthogonal arrival marks are two explicitly different targets.

## Written results

### R24.1 — Stopping at first return retains the source norm and closes renewal

For every finite \(N\), the tuple arrow
\(J_Nv=(A_1v,\ldots,A_Nv,S_Nv)\) preserves the native matching pairing:

\[
\sum_{n=1}^N A_n^\dagger A_n+S_N^\dagger S_N=I.
\tag{24.2}
\]

No branch is destroyed. The unrestricted return coefficients satisfy

\[
V_n=\sum_{k=1}^nV_{n-k}A_k\qquad(n\ge1).
\tag{24.3}
\]

**Proof.** Split \(US_{n-1}\) into its disjoint native address cuts
\(EA_n\) and \(S_n\). Since \(U\) preserves pairing and the two cuts have
zero mutual pairing,
\(S_{n-1}^\dagger S_{n-1}=A_n^\dagger A_n+S_n^\dagger S_n\).
Finite telescoping proves (24.2). Fresh arrival marks realize precisely
this direct-sum pairing. For (24.3), partition the finite source words that
end at zero by their first positive return time \(k\). Their first portion
is \(A_k\); their remaining unrestricted portion is \(V_{n-k}\), applied
after it. Every returning word has exactly one such split, with its signs
and order retained. Expanding both sides therefore gives identical native
word sums. This is an exact memory elimination, not a fitted recurrence.

### R24.2 — Native excursion signs determine one algebraic return series

Put \(B=I+R\). Odd first returns vanish. At event \(2m\),

\[
A_{2m}=2^{-m}t_m B,
\qquad T(w)=\sum_{m\ge1}t_m w^m,
\qquad T^2-(1+2w)T+w=0,\quad T(0)=0.
\tag{24.4}
\]

Consequently the formal first-return map is

\[
\mathcal A(z)=\sum_{n\ge1}A_nz^n
 =\frac{1+z^2-\sqrt{1+z^4}}2\,B,
\tag{24.5}
\]

where the square-root series is the unique branch with constant coefficient
one. The equation and branch are derived; they are not a supplied quadratic
constitutive law.

**Proof.** Every direction changes address by one, so returning words have
equal counts of the two directions. A right first excursion begins with
role 0, ends with role 1 and stays positive in between. Its transition sign
is \((-1)^{n_{11}}\), because only the native transition \(1\to1\) has
negative sign. Its role matrix is a signed multiple of
\(e_1(e_0+e_1)^\dagger\). The reflected left excursion has the same
internal sign: in a balanced word beginning with 0 and ending with 1 the
numbers of 00 and 11 adjacencies agree. Its role matrix is the same signed
multiple of \(e_0(e_0-e_1)^\dagger\). Their sum is
\(\left(\begin{smallmatrix}1&-1\\1&1\end{smallmatrix}\right)=I+KH=B\).

Let \(D(w)\) count all nonnegative balanced direction words, including the
empty word, with these signed transition weights and one \(w\) per pair
of steps. Unique splitting at returns gives \(D=1+T+T^2+\cdots=1/(1-T)\):
between two nonempty excursions the seam is 10, of positive sign. A
primitive excursion is either 01, of weight \(w\), or 0 followed by a
nonempty balanced word followed by 1. In the latter case the final 11 seam
adds a minus sign, giving \(T=w-w(D-1)=2w-wD\).
Eliminating \(D\) proves (24.4). Each coefficient of positive degree has
coefficient -1 in the term linear in \(T\); the zero constant term
therefore determines all coefficients uniquely. Completing the square
gives (24.5) after \(w=z^2/2\). All manipulations are formal coefficient
identities with finitely many terms at each degree.

### R24.3 — Finite word counts give all coefficients and exact cancellations

Let \(c_0=1\) and \(c_{k+1}=\sum_{j=0}^kc_jc_{k-j}\). Then

\[
c_k=\frac1{k+1}\binom{2k}{k},\qquad
A_2=B/2,\qquad
A_{4(k+1)}=\frac{(-1)^{k+1}c_k}{4^{k+1}}B\quad(k\ge0),
\tag{24.6}
\]

and all other first-return coefficients are zero. Thus possible six-event
first-return words exist, while their coherent aggregate is exactly zero.
Define

\[
a_k=\binom{2k}{k}/4^k,\quad b_k=a_k^2.
\]

Their exact recurrences imply

\[
\frac{a_{k+1}}{a_k}=\frac{2k+1}{2k+2},\qquad
0<b_k\le\frac1{k+1},\qquad
\frac{c_k}{4^k}=\frac{a_k}{k+1}\le\frac1{(k+1)^{3/2}}.
\tag{24.7}
\]

**Proof.** Count nonnegative balanced words of length \(2k\), ignoring
signs for this count alone. Splitting at the partner of the first opening
gives the displayed convolution. Among all \(\binom{2k}{k}\) balanced
words, reflect the prefix through the first address -1 for every word that
goes below zero. This is a reversible correspondence with words having
\(k+1\) positive and \(k-1\) negative steps. Subtracting their count
\(\binom{2k}{k-1}\) gives \(c_k\). These are derived finite word counts;
the familiar name Catalan is only a label.

Writing \(T=w+Y\) in (24.4) gives \(Y^2-Y-w^2=0\). The count series
\(C(s)=\sum c_ks^k\) satisfies \(C=1+sC^2\), so
\(T=w-w^2C(-w^2)\) satisfies the same equation and constant-term
condition. Uniqueness proves (24.6). For (24.7), cancel finite factorial
counts to get the ratio. Induct from \(b_0=1\), using
\(((2k+1)/(2k+2))^2\le(k+1)/(k+2)\); after clearing positive
denominators its difference has numerator \(3k+2>0\). Taking the C8
positive square root gives the last bound.

### R24.4 — The completed coherent return has the fixed scale sqrt(2) minus one

The series \(\sum_n A_n\) is absolutely Cauchy in the native scalar
completion. Its value is

\[
\mathcal A(1)=f_*B=gJ,\quad
f_*=1-1/\sqrt2,\quad J=(I+R)/\sqrt2,\quad
g=\sqrt2-1,
\tag{24.8}
\]

with \(J^\dagger J=I\). The positive scale \(g\) is the unique positive
root of \(g^2+2g-1=0\); its squared scale is \(3-2\sqrt2\).
Neither a free depth cutoff nor a fitted strength remains in this completed
coherent target.

For the scalar partial sum
\(f_M=\frac12-\frac14\sum_{k=0}^M(-1)^kc_k/4^k\), the error has the sign
of the next term and

\[
|f_*-f_M|\le\frac{c_{M+1}}{4^{M+2}}.
\tag{24.9}
\]

**Proof.** By (24.7) the absolute tail beyond \(M\) is bounded by
\(\frac14\sum_{k>M}(k+1)^{-3/2}\le1/(2\sqrt{M+1})\).
For the last inequality sum the elementary telescoping bounds
\((k+1)^{-3/2}\le2(k^{-1/2}-(k+1)^{-1/2})\), for \(k\ge1\).
Thus completion C8 applies. Absolute convergence also permits the
coefficient products in (24.4) at \(w=1/2\), by bounding omitted products
by bounded partial sums times their absolute tails. Hence
\(f_*^2-2f_*+1/2=0\). The positive terms \(c_k/4^k\) strictly decrease
to zero; grouping successive pairs bounds the alternating sum between its
consecutive partial sums and proves (24.9). In particular
\(1/4\le f_*\le3/8\), which selects the smaller quadratic root.
Finally \((I-R)(I+R)=2I\) gives \(J^\dagger J=I\) and (24.8).

The coherent sum is a constructed amplitude target. Folding orthogonal
arrival marks into one coefficient is not asserted to be a reversible,
norm-preserving erasure operation. The next result retains those marks.

### R24.5 — Retained arrival memory fixes a different coefficient with a rigorous tail

At horizon \(4(M+1)\), the total returned energy of any input \(v\) is
\(\eta_M B(v,v)\), where

\[
\eta_M=\frac12+\frac18\sum_{k=0}^M\frac{c_k^2}{16^k}
 =\left(2M+\frac12+\frac1{8(M+1)^2}\right)b_M.
\tag{24.10}
\]

The native completed scalar \(\eta\) is independent of the input role and
of a common translation of the origin-return target. It obeys

\[
0<\eta-\eta_M\le\frac1{16(M+1)^2},\qquad
\frac58<\eta<\frac{11}{16}<1.
\tag{24.11}
\]

The live branch retains energy \((1-\eta_M)B(v,v)\); its limiting response
is positive. No limit of the moving live vector itself is needed or claimed.

**Proof.** Every nonzero \(A_n\) is a scalar multiple of \(B\), and
\(B^\dagger B=2I\). Squaring (24.6), then matching the distinct arrival
marks, proves the sum in (24.10). For its finite evaluation put
\(d_k=16k+4+(k+1)^{-2}\). Direct multiplication by positive denominators
proves
\[
d_kb_k-d_{k-1}b_{k-1}=b_k/(k+1)^2\quad(k\ge1),\qquad d_0b_0=5.
\]
Thus \(\sum_{k=0}^M c_k^2/16^k=d_Mb_M-4\), proving the second formula.
From (24.7), every omitted squared term is at most \((k+1)^{-3}\).
Sum
\((k+1)^{-3}\le\frac12(k^{-2}-(k+1)^{-2})\)
for \(k>M\) to get (24.11). All omitted terms are positive. At \(M=0\)
the lower sum is \(5/8\) and the bound is \(1/16\); the tail estimate
is strict here, giving the strict upper bound. C8 completes this explicitly
Cauchy sequence. Equation (24.2) gives the live energy; translating every
address leaves all direction words and their signs unchanged.

The series, recurrence and rational enclosures define \(\eta\) without an
imported integral, special-function identity or numerical constant. A
possible identification of this series with a separately derived period
formula is not a premise or a claim of this packet.

### R24.6 — Normalizing a returned record cancels its return scale

For every nonzero return branch, and also for the coherently completed
return, the normalized output record is \(Jv/\sqrt{B(v,v)}\), up to an
overall sign. Its R20 cross-retention coefficient is

\[
\kappa_{\rm returned}
 =\frac{B(v,Hv)}{B(v,v)}=:h(v),\qquad -1\le h(v)\le1.
\tag{24.12}
\]

The same coefficient follows when arrival marks are retained and the whole
returned branch is normalized. In particular, the native inputs
\(e_0,e_1,Ce_0\), with \(C=(H+K)/\sqrt2\), give \(+1,-1,0\).
The scale \(g\) and the returned-energy fraction \(\eta\) cannot replace
this normalized overlap.

**Proof.** Direct expansion of the native identities gives
\(B^\dagger K B=2H\) and \(B^\dagger B=2I\). For a scalar multiple of
\(Bv\), division of its current by its energy cancels the scalar square
and gives (24.12). For several retained arrivals, both numerator and
denominator have the same sum of scalar squares, so the result is unchanged.
For \(v=ue_0+te_1\),
the quotient is \((u^2-t^2)/(u^2+t^2)\); nonnegativity of the squares
proves the bounds. Evaluating the three stated source preparations proves
the values. No independent interaction-strength parameter is introduced.

### R24.7 — A return-controlled native gate has an exact completed coupling law

Prepare the record source as a nonzero native role \(v\) at the origin and
retain the full tuple \(J_Nv\). On an independent system role apply R20's
\(W=\Pi_0\otimes I+\Pi_1\otimes K\) only to returned branches; leave the
live branch unchanged. This is a declared controlled composition of native
arrows, preserving the joint pairing. Normalize the whole record source,
not each returned branch separately.

At horizon \(N=4(M+1)\), matching all record marks preserves diagonal
system pair blocks and multiplies cross blocks by

\[
\boxed{\kappa_M(v)=1-\eta_M+\eta_Mh(v).}
\tag{24.13}
\]

The completed response and its finite error are

\[
\kappa_*(v)=1-\eta(1-h(v)),\qquad
|\kappa_* -\kappa_M|\le\frac{1-h(v)}{16(M+1)^2}.
\tag{24.14}
\]

The canonical native preparations give

| Record source | Completed system cross retention |
| --- | --- |
| \(e_0\) | \(1\) |
| \(e_1\) | \(1-2\eta\) |
| \(Ce_0\) | \(1-\eta\) |

Thus this interface fixes its strength once its native source preparation
is specified. In particular, it supplies a nontrivial dimensionless
coefficient \(1-\eta\) for the already available balanced source, with no
empirical constant input or tunable decay rate.

**Proof.** On returned records the relative action of the two system roles
is \(K\); on live records it is \(I\). Distinct arrival/live marks make
their cross pairings vanish. By (24.2) and R24.5, the live branch contributes
\(1-\eta_M\). The returned branches contribute
\(\sum A_n^\dagger K A_n/B(v,v)=\eta_M h(v)\), by R24.6.
Equal system roles instead give the full norm one. This proves (24.13)
on arbitrary finite system pair fields independent of this prepared record,
including off-address blocks. The finite controlled record operator is a
direct sum of \(K\)'s and \(I\), hence is a pairing-preserving involution;
the system-controlled extension is exactly the native tuple construction
of R20.1. Apply (24.11) to get (24.14) and the listed values. These are
deterministic matching sums over retained branches, not a probabilistic
mixture or an assumed stochastic channel.

### R24.8 — Full word memory changes the limiting interaction and blocks physical promotion

Use the same transport and first-origin-return routing, but also write a
fresh R20 slot after every direction step. Keep the complete direction word
as an orthogonal record, including on stopped branches. At horizon \(2m\),
the returned-energy fraction is

\[
\vartheta_m=\sum_{j=1}^m\frac{2c_{j-1}}{4^j}
 =1-a_m\longrightarrow1,
\tag{24.15}
\]

for every normalized initial role. Its returned role current is zero.
Applying the same return-controlled gate to the returned role while keeping
its word slots therefore gives

\[
\kappa^{\rm words}_{2m}=a_m\longrightarrow0,
\qquad 0<a_m\le1/\sqrt{m+1}.
\tag{24.16}
\]

This differs from (24.14), including \(1-\eta>0\) for the native balanced
source. Equal arrows before the additional memory writes do not make the
two joint protocols identical. Their first returned-energy difference is
at event six: arrival-only \(5/8\), full-word \(11/16\).

**Proof.** For each initial departure side there are \(c_{j-1}\)
first-return direction words at event \(2j\), by removing the first and
last steps and applying the unsigned word count in R24.3. Each individual
word has the same magnitude factor \(2^{-(2j-1)/2}\) times the appropriate
component of \(Cv\). Fresh word marks remove every pairing between
distinct words. Summing the two squared departure components gives
\(B(Cv,Cv)=B(v,v)\), proving the summand \(2c_{j-1}/4^j\).
Every word has one final role and distinct words have distinct records;
there is no within-record cross-role entry, hence the current is zero.
The ratio in (24.7) gives
\(a_{j-1}-a_j=2c_{j-1}/4^j\); telescoping proves (24.15). The bound in
(24.7) gives its limit and (24.16). The unchanged live branch contributes
\(1-\vartheta_m\), and every returned branch contributes zero current,
which proves the gate formula. Directly the first three return fractions
are \(1/2,1/8,1/16\), whereas R24.3 has a zero six-event coherent
coefficient. Earlier energy totals agree and the stated difference follows.

The two constructions are explicit witnesses that availability of native
cuts, memory and continuation alone does not select one physical coupling
for these interfaces. This does not rule out a later native theorem
selecting a physical interface or preparation. It prevents the narrower
results here from being relabelled as that missing theorem.

## Evidence, lineage and next physical boundary

The [derivation ledger](../04-operator-evolution/R24_DERIVATION_LEDGER.json)
binds these eight written proof sections to their exact source and check
groups. The [native adapter](../04-operator-evolution/native_return_coupling.cjs)
uses the unchanged RKF arithmetic and symbolic replayer. Independent literal
H/K enumeration, unrestricted and stopped propagation, word-count grammar,
renewal, branch-norm balance, full joint record gates and rational tail
enclosures check the stated formulas. The
[wrapper](../04-operator-evolution/verify_r24.py) replays the frozen R23
chain and preserves all 194 prior non-navigation files.

The origin gate rejects classical/admitted/open/comparison proof premises,
unpinned citations, cycles and a promotion of target constructions to
physical selection. It checks declared dependency metadata and bound proof
text. Neither this gate nor a finite exact PASS is automatic semantic
proof-assistant verification. The all-depth and completion claims have the
written proofs above. No upstream whole-engine PASS is newly claimed.

| Source at exact commit | Premise status here |
| --- | --- |
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native marks, signed/refined counts, H/K, matching pairing and scalar completion. Equality, finite counting and induction remain the stated metamathematical infrastructure. Native analytic comparison citations in R16 are not needed for the new series proofs. |
| [R18.1–R18.4](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Native transport on a constructed address target; physical choice of that transport/target is not supplied. |
| [R20.1–R20.4, R20.7](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | Native record gate, matching readout, fresh-word records and address cuts. Preparations and allocation remain explicit constructions. |
| [RKF canonical algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Native finite equality and proof replay; unchanged active engine. This is not a claim of complete formal or physical certification of that repository. |
| [Publications AS-1–AS-5](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-alpha-selection/THEOREM.md) | **Comparison only, mixed premises:** native R2 response with supplied cell weights; AS-3 admits a positive quadratic physical adapter. Neither that imported adapter nor an electromagnetic normalization is a proof premise here. Research-branch commit, not the public main snapshot. |
| [Publications NI-1–NI-4](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-return-identification/THEOREM.md) | **Comparison only:** given positive cell weights and specified probes identify a parameter. Identification is not native physical parameter selection. Research-branch commit. |

R24 fixes return-response coefficients of the inherited source construction
and shows exactly how one enters an actual native record interaction. It
also proves why amplitude scale, normalized overlap and retained return
energy must be distinguished. No spectral observer theorem is rebuilt.
The next unresolved selection is which joint source, retained memory and
continuation interface a physical system realizes. A physical charged
source, its interaction channel and action normalization must then be
derived before any of these dimensionless quantities can be called alpha.
