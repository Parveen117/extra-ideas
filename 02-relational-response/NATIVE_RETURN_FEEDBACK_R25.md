# R25: endogenous return feedback, native memory closure and boundary selection

Research owner: Monty Dabas. Development: 1 October 2026.

R24 calculated a first-return source and a gate applied to its returned
records. R25 closes that feedback: whenever the source returns to its own
origin, apply the same gate and continue the source, retaining its arrival
marks. Gate times now come from the native path itself. No separate timing
sequence, damping parameter, noise law or replacement record is supplied.

The resulting response has a four-component native memory, an exact
eight-return cycle, and one uniquely completed feedback operator. Every
admissible boundary continuation placed after sufficiently many returns
has vanishing influence. A quantitative bound also proves convergence in
the original event count, rather than only at selected return times.

This selects the response of an explicit native feedback construction.
It does not select that construction as nature's unique interaction.
Another causal, pairing-preserving native construction, with full word
memory, gives a different limit. The two constructions test exactly which
selection claim has and has not been established. No value is fitted to
alpha or identified as an electromagnetic coupling.

## Source and the finite feedback construction

Inherit R16's native marked roles, signed/refined coefficients, matching
pairing and C8 completion; R18's transport; R20's record gate; and R24's
proved first-return coefficients. This packet works on their real signed
two-role response sector, preserving full source-word provenance. It does
not replace generalized UGD by an ordinary complex primitive.

Use the source identities

\[
H^2=K^2=I,\quad HK=-KH,\quad R=KH,\quad R^2=-I,
\quad J=(I+R)/\sqrt2,\quad C=(H+K)/\sqrt2.
\tag{25.0}
\]

The source role moves under the unchanged \(U\). A separate two-role
system controls the available gate
\(W=\Pi_0^{s}\otimes I+\Pi_1^{s}\otimes K\), where the superscript
distinguishes system cuts from source roles. After each positive event:

1. Write the equality predicate “source address is its initial origin”
   into a fresh native record slot.
2. On the equality branch apply \(W\), then continue both branches under
   the same source transport.

At any finite horizon this is a composition of the already constructed
record and controlled arrows. All flags, source coordinates and roles are
retained when matching the final environment. There is no declaration at
a finite event that a source will never return. The initial event is not
counted as a return. The initial source role is any nonzero \(v\).

Writing these arrival predicates and using this local controlled arrow
are target/interface constructions. Their unique physical selection is
not a premise. Event count remains source-word depth, not physical time.

For the first-return maps of R24 write
\(A_n=q_nJ\) and \(w_n=q_n^2\). A zero return coefficient gives \(w_n=0\).
These are derived squared coefficients, not postulated probabilities:

\[
w_2=\tfrac12,\quad
w_{4(k+1)}=\frac{c_k^2}{8\,16^k},\quad
c_k=\frac1{k+1}\binom{2k}{k},\quad
s_N=1-\sum_{n=1}^Nw_n,
\quad \eta=\sum_{n\ge1}w_n.
\tag{25.1}
\]

All other \(w_n\) vanish, \(s_0=1\), and R24 proves
\(5/8<\eta<11/16<1\), \(S_N^\dagger S_N=s_NI\), and a rigorous
tail for the count series. Source norms below are the matching norm already
derived in R16, not a supplied Hilbert or stochastic structure.

## Written results

### R25.1 — Native returns generate the feedback pulse and its clock

On a nonzero first-return branch, factoring out the common signed return
coefficient leaves the pairing-preserving pulse

\[
\mathcal T=W(I\otimes J)
 =\Pi_0^s\otimes J+\Pi_1^s\otimes C,
\quad
\mathcal T^4=-H_s\otimes I,\quad \mathcal T^8=I.
\tag{25.2}
\]

Four completed returns restore the source role and leave a system sign
contrast, independently of its prepared source role. Eight returns restore
the full pulse. Their timing is the path's return history; the four or eight
returns need not occupy a fixed number of transport events. The total
matching norm of all histories completing at least \(k\) returns is
\(\eta^k\) times the initial norm.

**Proof.** R24 gives \(A_n=q_nJ\). The system-0 branch leaves its
returned source unchanged; system-1 applies \(K\). Since \(KR=H\),
\(KJ=(K+H)/\sqrt2=C\), proving the pulse. Direct multiplication gives
\(J^\dagger J=C^\dagger C=I\), \(J^2=R\), \(J^4=-I\), and
\(C^2=I\). The system cuts are disjoint, so powers separate as in
(25.2). The source-0 block \(J\) has no smaller positive period than
eight: its even powers are \(I,R,-I,-R\), and odd powers have nonzero
coefficients of both \(I\) and \(R\). These native maps are distinct.

For a fixed list of return gaps \((n_1,\ldots,n_k)\), the branch map is
\(\prod_iq_{n_i}\) times \(\mathcal T^k\). Its norm multiplier is
\(\prod_iw_{n_i}\), because the pulse preserves pairing. Distinct gap
lists have distinct arrival marks. Finite matching sums factor, and the
positive scalar completion of those sums gives \(\eta^k\). This is a
derived word-ledger product, not an independence assumption. Four-return
phase restoration is not promoted to an unconditional four-event gate;
its returned branch carries the derived norm fraction \(\eta^4\).

### R25.2 — The exact feedback memory has four native components

For a source response matrix define the derived continuation arrow

\[
\mathscr L(X)=J^\dagger XC
 =\tfrac12(I-R)X(H+K).
\tag{25.3}
\]

It obeys

\[
I\ \mapsto\ H\ \mapsto\ -R\ \mapsto\ -K\ \mapsto\ -I,
\qquad \mathscr L^4=-\mathrm{id},\quad \mathscr L^8=\mathrm{id}.
\tag{25.4}
\]

The linear continuation closure generated by \(I\) is exactly
\(\operatorname{span}\{I,H,R,K\}\), of dimension four. Although
\(B(v,Rv)=0\) for every real source role, the \(R\) component cannot be
dropped before continuing: \(\mathscr L(R)=K\).

**Proof.** Expand (25.3) with \(R=KH\) and the two source involutions
to get the four displayed images. More generally
\(\mathscr L^k(X)=(J^\dagger)^kXC^k\), so the fourth and eighth
powers follow from R25.1. If
\(aI+bH+cR+dK=0\), the two diagonal native role entries force
\(a=b=0\), and the two off-diagonal entries force \(c=d=0\).
Thus the first four iterates of \(I\) are independent and already span
the full real two-role matrix target. Every continuation-closed linear
subspace containing \(I\) must contain them. No generic minimum-observer
theorem is reconstructed here; this is the explicit closure of this pulse.

Since \(R^\dagger=-R\), the real scalar \(B(v,Rv)\) is its own
negative and vanishes. Nevertheless \(\mathscr L(R)=K\), whose pairing
with \(Ce_0\) is nonzero. A presently invisible component can therefore
feed a future response. Replacing the full return arrow by a scalar
retention at each step loses a source-required component.

### R25.3 — The full finite-event evolution has an exact retained-memory recurrence

Let \(D_N\) be the source cross-response operator obtained by matching all
environment coordinates between the two system branches after \(N\)
events. Then

\[
D_0=I,\qquad
\boxed{D_N=s_NI+\sum_{n=1}^Nw_n\mathscr L(D_{N-n}).}
\tag{25.5}
\]

The diagonal system pair blocks are preserved. For an independent prepared
source \(v\), a cross block is multiplied by
\(\kappa_N(v)=B(v,D_Nv)/B(v,v)\). Equivalently,

\[
D_N=\sum_{k\ge0}p_{N,k}\mathscr L^k(I),\quad
p_{N,k}=\sum_{n_1+\cdots+n_k\le N}
 \left(\prod_iw_{n_i}\right)s_{N-\sum_i n_i},\quad
\sum_kp_{N,k}=1.
\tag{25.6}
\]

The empty-list term is \(p_{N,0}=s_N\); all sums at finite \(N\) are
finite. These coefficients arise from the actual joint record evolution.

**Proof.** Split native histories at their first positive return. Histories
with no return contribute \(S_N^\dagger S_N=s_NI\), because both system
branches undergo the same source evolution there. A first return at \(n\)
has source maps \(q_nJ\) and \(q_nC\). Its remaining cross response is
\(D_{N-n}\), yielding
\(w_nJ^\dagger D_{N-n}C\). Arrival flags prevent cross pairing between
different first-return times. Summing proves (25.5) with the correct
operator order. Repeatedly apply the same split to obtain (25.6).

For each fixed arrival list the final live source norm is
\(\prod_iw_{n_i}s_{N-\sum n_i}\) in either system branch, independent
of the prepared source role. Every finite step consists of pairing-
preserving arrows and retains all record flags. Its total norm is one,
so these coefficients sum to one. Thus no stochastic channel or renewed
independent record state was assumed. In particular
\(D_2=(I+H)/2\) and
\(D_4=3I/8+3H/8-R/4\); the antisymmetric memory is already present.

### R25.4 — Completing all returns fixes one exact feedback response

The finite-event response converges, in every native bilinear entry, to

\[
\boxed{D_*=
 \frac{1-\eta}{1+\eta^4}
 (I+\eta H-\eta^2R-\eta^3K).}
\tag{25.7}
\]

Writing
\(h(v)=B(v,Hv)/B(v,v)\),
\(j(v)=B(v,Kv)/B(v,v)\), its scalar cross retention is

\[
\boxed{\kappa_*(v)=
 \frac{(1-\eta)(1+\eta h(v)-\eta^3j(v))}{1+\eta^4}.}
\tag{25.8}
\]

No new constant is supplied: \(\eta\) is precisely the R24 count-series
coefficient. The coherent amplitude scale \(\sqrt2-1\) is not substituted
for this arrival-retained energy coefficient.

**Proof.** Complete the nonnegative scalar ledger of finite arrival lists
by assigning a terminal remainder coefficient \(1-\eta\) after the last
return. A list with \(k\) gaps has weight
\((1-\eta)\prod_iw_{n_i}\). Finite sums and their C8 limits give total
weight \((1-\eta)\sum_{k\ge0}\eta^k=1\); the omitted return-count tail
is \(\eta^{k+1}\). This construction concerns a completed response ledger,
not a finite-time declaration of escape or an assumed path probability.

Its marginal before the next return at age \(r\) has weight
\((1-\eta)+\sum_{n>r}w_n=s_r\). Consequently grouping its lists by
their prefix before event \(N\) gives exactly (25.6). For a fixed finite
list, that prefix eventually contains the whole list. To pass to the
completed sum, first bound the large return-count tail by a geometric
tail, then the long-gap tails of the finitely many remaining positions
by R24's Cauchy bound. Each matrix coefficient of
\(\mathscr L^k(I)=(J^\dagger)^kC^k\) is bounded by one from native
pairing preservation. Those tails therefore tend to zero uniformly.
This proves convergence to
\((1-\eta)\sum_{k\ge0}\eta^k\mathscr L^k(I)\); R25.6 gives an
explicit event-horizon bound.

Group successive sets of four terms and use \(\mathscr L^4=-\mathrm{id}\).
The convergent scalar series in \(-\eta^4\) gives
\[
(\mathrm{id}-\eta\mathscr L)^{-1}I
 =\frac{I+\eta H-\eta^2R-\eta^3K}{1+\eta^4}.
\]
This inverse is also checked directly by multiplication, so (25.7)
follows without a resolvent or convergence theorem imported from elsewhere.
R25.2 removes the antisymmetric term only at the final real quadratic
readout, proving (25.8).

### R25.5 — All admissible far-boundary continuations give the same limit

Suppose a terminal continuation after \(m\) completed returns is made from
native pairing-preserving arrows \(V_0,V_1\) on the two system branches.
Express these arrows on the two source input roles, including the common
normalized retained-prefix record embedding before any processing of those
records. Its relative response is \(G=V_0^\dagger V_1\). After completing the
intervening gap sums, its response is

\[
D^{[m]}(G)=(1-\eta)\sum_{k=0}^{m-1}\eta^k\mathscr L^k(I)
           +\eta^m\mathscr L^m(G).
\tag{25.9}
\]

For any two such terminal choices, and arbitrary source roles \(u,v\),

\[
|B(u,[D^{[m]}(G_1)-D^{[m]}(G_2)]v)|
 \le2\eta^m\|u\|\|v\|.
\tag{25.10}
\]

The same upper bound holds for the difference between \(D^{[m]}(G)\)
and \(D_*\). Furthermore \(D_*\) is the only real two-role solution of

\[
D=(1-\eta)I+\eta\mathscr L(D).
\tag{25.11}
\]

This removes terminal-boundary freedom in this feedback construction.
It does not establish that every native memory construction has the same
response.

**Proof.** Split lists into those ending before the \(m\)-th return and
those reaching it. R25.1 gives their respective weights in (25.9).
On the latter lists, all gap-dependent scalar coefficients and retained
marks form a common prefix ledger, while the two source role maps are
\(J^m\) and \(C^m\). Its norm is \(\eta^m\). Factoring that norm
defines the embedding used in \(V_0,V_1\); the same bound follows first
on finite prefix ledgers and then by their proved scalar norm tails.
Thus processing retained marks is included, not silently replaced by a
source-only terminal gate. The terminal response may vary with \(m\).
The native sum-of-squares inequality from R16 gives
\(|B(u,Gv)|=|B(V_0u,V_1v)|\le\|u\|\|v\|\).
Left and right multiplication by powers of \(J\) and \(C\) preserve this
bound. The triangle inequality proves (25.10). The completed tail has the
same bilinear bound, because it is a normalized positive sum of such
relative maps; hence comparison with \(D_*\) has the same bound.

Equation (25.7) solves (25.11) by direct multiplication. The difference
\(X\) of two solutions obeys \(X=\eta\mathscr L(X)\).
Applying this equation four times yields \(X=-\eta^4X\), whence
\((1+\eta^4)X=0\). Its positive scalar coefficient is invertible, so
\(X=0\). Uniqueness is derived algebraically, not imposed as a fixed-point
axiom. The boundary estimate further proves which finite native cutoffs
actually approach that unique value.

### R25.6 — The event-horizon error is bounded without an external decay law

For every integer \(N\ge0\) and source roles \(u,v\),

\[
|B(u,(D_N-D_*)v)|
\le\min\!\left\{2,
 \frac{8(1+4\eta+\eta^2)}{(1-\eta)^3(N+1)^2}\right\}\|u\|\|v\|
\le\min\!\left\{2,\frac{138368}{125(N+1)^2}\right\}\|u\|\|v\|.
\tag{25.12}
\]

Thus the completed response is selected by finite forward evolution as
well as by return-count truncation. Return count and event depth have
different derived error laws; neither is converted to a physical clock.

**Proof.** Let \(\epsilon_L=\sum_{n>L}w_n\). R24 proves, for
\(L\ge4\),
\(\epsilon_L\le1/(16\lfloor L/4\rfloor^2)\).
Since \(L+1\le8\lfloor L/4\rfloor\), this is at most
\(4/(L+1)^2\). For \(L=0,1,2,3\), use
\(\eta<11/16\) and \(w_2=1/2\) to check the same upper bound.

In the completed ledger of R25.4, a list of \(k\ge1\) return gaps can
have a return after event \(N\) only if some gap exceeds
\(\lfloor N/k\rfloor\). Sum the weight of that possibility over its
\(k\) positions, allowing overcounting. Its weight is at most
\[
(1-\eta)\,k\,\eta^{k-1}\epsilon_{\lfloor N/k\rfloor}
\le \frac{4(1-\eta)k^3\eta^{k-1}}{(N+1)^2}.
\]
Here \(\lfloor N/k\rfloor+1\ge(N+1)/k\). Outside these lists the
finite prefix and the full list give exactly the same iterate of
\(\mathscr L\); on each remaining list their bilinear difference is at
most \(2\|u\|\|v\|\). Summing proves the first nontrivial bound,
using
\[
\sum_{k\ge1}k^3\eta^{k-1}
 =\frac{1+4\eta+\eta^2}{(1-\eta)^4}.
\]
For completeness, this identity follows by multiplying the coefficient
series by \((1-\eta)^4\): the fourth finite difference of \(k^3\)
vanishes, and the first three coefficients are \(1,4,1\).
The resulting polynomial times geometric tails vanish; expand
\((1/\eta)^k=(1+(1/\eta-1))^k\) and retain a fixed binomial term of
degree at least five to bound them by a summable inverse-square tail.
Thus neither differentiation of an imported function nor a random-walk
tail theorem is required. Each full response has bilinear bound one,
giving the alternative bound two. Finally insert \(\eta<11/16\) in
the positive increasing numerator and reciprocal denominator to get
\(8(1+4(11/16)+(11/16)^2)/(1-11/16)^3=138368/125\).

### R25.7 — The completed response has a sharp source-preparation range

Put \(q=(1-\eta)/(1+\eta^4)\). For every nonzero real source preparation,

\[
q\bigl(1-\eta\sqrt{1+\eta^4}\bigr)
\le\kappa_*(v)\le
q\bigl(1+\eta\sqrt{1+\eta^4}\bigr),
\tag{25.13}
\]

and both extrema are attained by native completed preparations. The lower
bound is strictly positive; every value is less than one. Three already
available preparations give

| Source role | Completed feedback retention |
| --- | --- |
| \(e_0\) | \((1-\eta^2)/(1+\eta^4)\), approximately \(0.51081\) |
| \(e_1\) | \((1-\eta)^2/(1+\eta^4)\), approximately \(0.11342\) |
| \(Ce_0\) | \((1-\eta)(1-\eta^3)/(1+\eta^4)\), approximately \(0.23158\) |

These numbers are evaluated with rational enclosures in the certificate.
The response operator is fixed by the source construction; its scalar
readout still depends on the prepared role. In particular the first two
values differ by \(2q\eta>0\).

**Proof.** For \(v=ue_0+te_1\), expand
\(h=(u^2-t^2)/(u^2+t^2)\) and \(j=2ut/(u^2+t^2)\) to obtain
\(h^2+j^2=1\). The native two-component sum-of-squares identity gives
\(|h-\eta^2j|\le\sqrt{1+\eta^4}\), proving (25.13).
Equality occurs at
\((h,j)=\pm(1,-\eta^2)/\sqrt{1+\eta^4}\).
Construct the corresponding normalized coefficients by
\(u=\sqrt{(1+h)/2}\), \(t=j/(2u)\); both target \(h\)'s are greater
than -1. C8 supplies the square roots, and direct substitution verifies
the preparation. Since
\(\eta^2+\eta^6<(11/16)^2+(11/16)^6<1\), the lower bound is positive.
Each return-count contribution to a normalized quadratic response is at
most one, and the fourth contribution is -1 with positive weight
\((1-\eta)\eta^4\). Hence
\(\kappa_*\le1-2(1-\eta)\eta^4<1\).
Substitute the three native preparations in (25.8) to get the table and
their positive separation. No source-preparation independence is claimed.

### R25.8 — Full word memory gives a different causal continuation law

Keep the same return-controlled gate, source and origin test, but also
write a fresh direction record after every transport event, as in R24.8.
The within-excursion cross-continuation map becomes

\[
\mathscr L_{\rm word}(X)
 =X_{10}\frac{I+K}{2}+X_{01}\frac{I-K}{2},
\tag{25.14}
\]

so

\[
I\mapsto0,\quad H\mapsto0,\quad R\mapsto K\mapsto I\mapsto0,
\qquad \mathscr L_{\rm word}^3=0.
\tag{25.15}
\]

Its full finite-event feedback response is exactly

\[
D_N^{\rm word}=a_{\lfloor N/2\rfloor}I,
\quad a_m=\binom{2m}{m}/4^m,
\qquad D_N^{\rm word}\longrightarrow0.
\tag{25.16}
\]

The original arrival-only protocol instead has the strictly positive
scalar range (25.13). Both protocols are causal, preserve the full joint
pairing, trigger at the same native origin returns and use the same
available system/source gate. Those consistency properties alone therefore
do not force one memory interface or one physical interaction strength.

**Proof.** A fully recorded first-return direction word starts in direction
\(s\) and ends in role \(1-s\). Apart from its signed scalar magnitude,
its source map is \(e_{1-s}e_s^\dagger C\). With system-1's terminal
exchange, its cross response is
\[
(XK)_{1-s,1-s}\,C^\dagger\Pi_sC.
\]
The two reflected word families have the same unsigned count. Factoring
out their common returned-energy coefficient
\(\widetilde w_{2m}=2c_{m-1}/4^m\), and using
\(C^\dagger\Pi_sC=(I+(-1)^sK)/2\), yields (25.14).
Evaluation on the four native basis matrices proves (25.15), including
that the third power vanishes while the second need not.

The same first-return split as in R25.3 gives
\[
D_N^{\rm word}=\widetilde s_N I+
 \sum_{n=1}^N\widetilde w_n\mathscr L_{\rm word}(D_{N-n}^{\rm word}).
\]
R24.8 proves \(\widetilde s_N=a_{\lfloor N/2\rfloor}\) and
\(a_m\le1/\sqrt{m+1}\). Starting from \(I\), induction and
\(\mathscr L_{\rm word}(I)=0\) prove (25.16), with its explicit native
tail. Unlike the arrival-only ledger, this full-word ledger has completed
first-return fraction one; the escape-based completion argument is not
silently reused there. Its independent finite proof establishes the limit.

Writing either set of records is a forward native equality-controlled
exchange, and each controlled source/system gate preserves pairing. No
future return is consulted to decide an earlier operation. Thus both
examples meet the stated consistency tests. Their different responses
exclude promotion of those tests into unique physical selection within
these constructions. A further source theorem could still select a physical
interface; no universal impossibility claim is made.

## Certification, citations and the remaining selection

The [ledger](../04-operator-evolution/R25_DERIVATION_LEDGER.json) binds all
eight written proofs to exact native sources and scoped checks. The
[adapter](../04-operator-evolution/native_return_feedback.cjs) independently
propagates the full finite event/record fields and compares them with the
return-memory recurrence. It checks native cycles, tail bounds, terminal
continuations, completed response identities, rational enclosures and the
alternative full-word dynamics. The
[wrapper](../04-operator-evolution/verify_r25.py) replays the unchanged R24
and earlier certificate chain and preserves all 201 prior non-navigation
files.

Native definitions, mathematical deductions, finite checks, formal proof
verification and physical selection remain different statuses. The
dependency gate rejects admitted/classical/open/comparison proof premises,
unpinned citations and physical promotion of interface definitions. It
checks the declared graph and bound proof text, not arbitrary proof
semantics. This packet is not proof-assistant formalization and does not
claim a new whole-engine PASS.

| Pinned source | Premise status in R25 |
| --- | --- |
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native marked roles, H/K, signed/refined counts, matching pairing, completion and retained word identities. Equality, finite counting and induction are explicit metamathematical infrastructure. |
| [R18.1–R18.4](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Native transport on a constructed address target; no physical clock or unique physical transport selection. |
| [R20.1–R20.4, R20.7](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | Native tuple pairing, record gate, matching readout and available fresh/address records. Allocation and preparation remain constructions. |
| [R24.1–R24.8](https://github.com/Parveen117/extra-ideas/blob/d660bbf026e1de818fcb3dd61d46346a4ad618ae/02-relational-response/NATIVE_RETURN_COUPLING_R24.md) | Native signed first-return maps, norm balance, count coefficients and proved completion tails. Coherent amplitude, arrival energy and normalized overlap are distinct. No classical return theorem or electromagnetic identification imported. |
| [RKF canonical algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Unchanged native finite identity/replay engine; scoped replay only. |
| [Publications AS-1–AS-5](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-alpha-selection/THEOREM.md) | **Comparison only, mixed premises:** native R2 responses with supplied cell weights; AS-3 admits a positive quadratic physical adapter. That imported adapter and electromagnetic normalization are excluded from these proof paths. Research-branch commit. |

R25 selects the all-return response independently of an arbitrary terminal
continuation and derives its finite-event approach. It also narrows the
remaining physical question: causality and pairing preservation do not
choose whether the relevant source retains only arrival distinctions or
the full direction word. Selecting that joint source/memory interface, and
then deriving its physical source and interaction normalization, remains
necessary before claiming a physical alpha. The scalar readout's ordinary
preparation dependence must not be confused with a free coupling parameter.
