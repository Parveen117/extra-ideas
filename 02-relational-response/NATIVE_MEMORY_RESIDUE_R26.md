# R26: returning memory, native commutator residue and a propagation gap

Research owner: Monty Dabas. Development: 1 October 2026.

R26 makes memory an active, repeatedly used part of the source. It first
removes two large classes of memory choices that leave the declared response
unchanged. It then derives an observable commutator square for the remaining
choices. Reusing one native H/K record yields an exact propagation gap, a
new signed first-return law and a feedback channel ratio fixed at three.
No classical curvature, noise law, wave equation or physical mass term is
inserted. The needed scalar completions and their error bounds are derived.

The gap below is specifically a lower bound for the native two-event propagation
defect. It is not an identification with particle mass, a Hamiltonian mass
gap, the Yang–Mills problem or a physical unit. The direction-controlled
record interface is an explicit construction from the existing source;
its selection as a unique physical interaction remains open.

## Source and controls

Use R16 C1, C3–C8 and C13, R18.1–R18.4, R20's tuple pairing and controlled
records, and R24–R25's return decomposition. Keep the active RKF engine
unchanged. Work on the real signed finite source/record ledger, except for
the explicitly labelled native scalar-turn argument in R26.6 and R26.9.
Its iota is the already derived cut turn, not an ordinary complex primitive.

On the moving source role put
\(R_s=K_sH_s\), \(F_s=H_s+K_s\), \(C_s=F_s/\sqrt2\),
\(\Pi_a=(I+(-1)^aH_s)/2\), \(L_a=\Pi_aF_s\).
The address shift is \((S\psi)(x)=\psi(x-1)\). The unchanged transport
is \(U=(S\Pi_0+S^{-1}\Pi_1)C_s\).

A record has its own finite native role target. For pairing-preserving arrows
\(A,B\), define the available direction-controlled transport

\[
V_{A,B}=\bigl(S\Pi_0\otimes A+S^{-1}\Pi_1\otimes B\bigr)
        (C_s\otimes I).
\tag{26.0}
\]

Thus A follows an outgoing positive direction and B follows an outgoing
negative direction. These controls are not the original H/K letter tags:
either original letter can lead to either direction depending on the
incoming role. Full source words remain distinct as provenance. R26.9
separately examines copying the original letters.

All field norms and pairings are finite sums over native matching marks.
The shift permutes addresses. The two source cuts are disjoint; hence the
first factor of (26.0) and \(C_s\) preserve the joint pairing, and
\(V_{A,B}\) has their reversed inverse. No new tensor or measurement
postulate is required.

## Written results

### R26.1 — Flat returning records do not change the declared source response

Two classes are exactly removable from local source readouts for a source
initially at one address with an independent prepared record.

First, attach a pairing-preserving arrow \(A_x\) to the positive crossing
of edge \(x\to x+1\), and its inverse to the reverse crossing. Define
\(G_0=I\), \(G_{x+1}=A_xG_x\), and
\((\mathcal D\Psi)(x)=(I\otimes G_x)\Psi(x)\). The resulting transport
is

\[
V_{\rm reciprocal}=\mathcal D(U\otimes I)\mathcal D^{-1}.
\tag{26.1}
\]

Second, for homogeneous commuting \(A,B\), every path at event n and
address x has the same record arrow
\(G_{n,x}=A^{(n+x)/2}B^{(n-x)/2}\).
The joint field is the unrecorded source field tensored with
\(G_{n,x}\mu\). Every unresolved local energy and current is unchanged.

Both statements also hold with R25's origin-return flags and its separate
system-controlled K gate on the moving source role: every finite-event
system cross response is unchanged. Record labels themselves are retained.

**Proof.** On a positive crossing, (26.1) has record arrow
\(G_{x+1}G_x^{-1}=A_x\); on a negative crossing it has \(A_x^{-1}\).
The role coin commutes with record arrows, proving the equality. It supplies
the common factor \(G_x\) at a fixed endpoint. Commutation in the second
case allows every record word with the same direction counts to be reordered
into \(G_{n,x}\). Each address/event determines those counts, so the same
factor leaves the entire source path sum intact. Pairing preservation of G
then removes it from local matching readouts.

Origin flags depend only on address, and the separate gate acts only on
the moving source role and its system controller. Neither changes these
factorizations. Within each retained arrival history and final address,
the two controller branches acquire the same record arrow, which cancels
from their cross pairing. These are exact source-readout equivalences,
not an assertion that the record histories have been erased or that a
multi-address preparation with unrelated record factors has the same scope.

### R26.2 — A two-step return measures the native commutator square

Start the moving source in \(e_0\) at zero and a nonzero record \(\mu\).
Define

\[
\Delta_\mu(A,B)=\frac{\|(AB-BA)\mu\|^2}{\|\mu\|^2},\quad
0\le\Delta_\mu\le4.
\]

The two-event return energy and current are exactly

\[
\rho_2(0)=\tfrac12,\qquad
\boxed{j_2(0)=\tfrac12-\tfrac14\Delta_\mu(A,B).}
\tag{26.2}
\]

Zero defect on every record preparation is equivalent to \(AB=BA\).
Maximal defect four on every preparation is equivalent to \(AB=-BA\).
For a fixed preparation, forward-order independence instead requires zero
defect on its whole native word orbit, not merely on the initial vector.

The path residue cannot be removed by a pairing-preserving change of record
frame at each event/address. Both paths have the same endpoints, so their
relative arrow is only conjugated at the starting frame.

**Proof.** The returned field is
\((e_0\otimes AB\mu+e_1\otimes BA\mu)/(2\|\mu\|)\).
Each record arrow preserves norm, giving the energy. Its current is half
the normalized real pairing of \(AB\mu\) and \(BA\mu\).
Expanding the squared difference gives (26.2). Expanding the squared sum
gives \(\|(AB+BA)\mu\|^2/\|\mu\|^2=4-\Delta_\mu\).
Positivity proves the bounds and both all-preparation equivalences.

If the commutator vanishes on every word applied to \(\mu\), any adjacent
AB pair inside a longer word can be exchanged for BA: its right suffix
prepares exactly one of those orbit vectors. Repeated swaps prove order
independence. Conversely, the two extended words ABw and BAw test each
orbit vector. A single initial null response is not this complete condition.
Under frame changes, each full path arrow becomes
\(G_{\rm end}OG_{\rm start}^{-1}\). End frames cancel in the relative
arrow, and the difference norm is unchanged on the transformed preparation.

For native \(A=H_m,B=K_m\), the defect is four for every \(\mu\).
Thus \(j_2(0)=-1/2\), and the next source energies at
\((-3,-1,1,3)\) are \((1,5,1,1)/8\). The latter follows by expanding
the source coin and shifts, or from the inherited native energy/current
balance. The record phase changes an actual source readout.

### R26.3 — One reused H/K record changes path signs without a free record state

Take \(A=H_m,B=K_m\), with the inherited native relations, and write
\(V=V_{H_m,K_m}\). For a one-origin product preparation, at event n and
address x the joint field factors as

\[
\Psi_n(x)=\phi_n(x)\otimes H_m^{a}K_m^{b}\mu,
\quad a=(n+x)/2,\quad b=(n-x)/2.
\tag{26.3}
\]

Every local source response is therefore independent of the normalized
record preparation. The source coefficients obey the exact sign stencil

\[
\phi_{n+1}(x)=2^{-1/2}\left[
L_0\phi_n(x-1)+(-1)^{(n+x+1)/2}L_1\phi_n(x+1)\right]
\tag{26.4}
\]

on the native reachable parity. The two orderings of a return have relative
record arrow \(-I\); a common endpoint frame cannot remove this sign.

**Proof.** Normal-order every record word into \(H_m^aK_m^b\).
Each exchange of unlike adjacent letters contributes the derived sign
-1. On a positive step, left multiplication by H requires no new exchange.
On a negative step, \(K_mH_m^a=(-1)^aH_m^aK_m\), giving (26.4).
All paths at the fixed endpoint have the same unsigned counts, so this
common record factor can be taken out of their signed source sum.
Its norm and matching self-pairing are unchanged, proving record-state
independence of the local source response. Nevertheless HK and KH differ
by -1, and their relative arrow is -I. This is an order residue in source
amplitudes, not loss of norm into an unobserved random force.

### R26.4 — Anticommuting memory generates an exact native propagation gap

More generally let the two record arrows be pairing-preserving involutions
with \(AB=-BA\). Their transport obeys the operator identity

\[
\boxed{V^4-\tfrac12(S^2+S^{-2})V^2+I=0.}
\tag{26.5}
\]

For every finitely supported joint field,

\[
\boxed{\|(V^2-I)\Psi\|^2
 =\|\Psi\|^2+\tfrac12\|(S^2-I)\Psi\|^2,}
\tag{26.6}
\]

and consequently the ratio to \(\|\Psi\|^2\) lies between one and three.
The two bounds are sharp as an infimum and supremum. In particular there
is no nonzero finite-energy approximate two-event stationary sequence with
normalized defect tending to zero. No physical energy or mass unit has
been assigned to this defect.

**Proof.** Put \(X=SA\), \(Y=S^{-1}B\), \(a=S^2\), \(b=S^{-2}\),
and \(Z=AB\). The unnormalized source block is
\(M=\left(\begin{smallmatrix}X&X\\Y&-Y\end{smallmatrix}\right)\).
Here \(XY=-YX=Z\), \(Z^2=-I\), and \(ab=I\). Direct multiplication gives
\[
M^2=\begin{pmatrix}a+Z&a-Z\\-b-Z&b-Z\end{pmatrix},\qquad
M^4-(a+b)M^2+4I=0.
\]
Dividing by the derived powers of \(\sqrt2\) proves (26.5). This is a
direct native identity, not an imported wave or characteristic-polynomial
theorem. Pairing preservation supplies \(V^\dagger=V^{-1}\). Multiply
(26.5) by \(V^{-2}\), then expand the left-hand norm in (26.6).
It is the pairing of \(2I-(S^2+S^{-2})/2\), which is precisely its right
side. The native triangle inequality gives the upper bound three.
If finite-energy fields are included, they mean the completion of finite
fields under this derived matching norm. Both V and S preserve the norm,
so they carry Cauchy sequences and equivalent sequences to the same kind
of sequences. Their extensions, and the displayed norm identity, therefore
follow by taking native Cauchy limits; no independent infinite carrier
or Hilbert-space axiom is supplied.

For sharpness, take the same nonzero source/record vector at the M even
addresses \(2,4,\ldots,2M\) and zero elsewhere. Its shift-difference
energy is twice the one-vector energy, giving ratio \(1+1/M\).
Alternating the signs on those addresses gives ratio \(3-1/M\).
By comparison, for the unrecorded source prepared in role \(e_0\), the
same constant packet has \(\|(U^2-I)\psi\|^2/\|\psi\|^2=1/M\),
as direct expansion of U shows. The positive term in (26.6) is thus a
consequence of the specified anticommuting memory, not a supplied mass term.

### R26.5 — Signed first excursions determine the new return grammar

Return to the explicit H/K record. Let \(\mathcal A_n\) be the first-origin-
return map on source and record, with no earlier return, and put
\(B_s=I+R_s\), \(D_s=K_s-H_s\). Then

\[
\mathcal A_{2m}=\frac{t_m}{2^m}
 \begin{cases}B_s\otimes R_m^m,&m\text{ even},\\
 D_s\otimes R_m^m,&m\text{ odd},\end{cases}
\quad \mathcal A_{2m+1}=0.
\tag{26.7}
\]

Its signed word series \(T(w)=\sum_{m\ge1}t_mw^m\), with zero constant
term, is uniquely determined by

\[
T(w)=2w-\frac{w}{1-T(-w)},\qquad
(1+2w)T^2-(1+2w+4w^2)T+w(1+2w)=0.
\tag{26.8}
\]

If \(Q(x)=\sum_{k\ge1}q_kx^k\), then

\[
Q(x)=\frac{1-\sqrt{1-x+x^2}}{2(1-x)},\quad
(1-x)Q^2-Q+x/4=0,\quad
T(w)=w+(1-2w)Q(4w^2).
\tag{26.9}
\]

Thus \(t_1=1\), \(t_{2k}=a_k=4^kq_k\), and
\(t_{2k+1}=-2a_k\) for \(k\ge1\). The first values of \(a_k\) are
\(1,1,-2,-11,-14,58,316,397\).

**Proof.** A balanced record word with m positive and m negative steps is
a signed multiple of \(R_m^m\). In a right primitive excursion its final
source row is \(e_1(e_0+e_1)^\dagger\). Reflection interchanges the two
record letters and contributes \((-1)^{m^2}=(-1)^m\); the native source
internal transition sign is unchanged. The reflected source row is
\(e_0(e_0-e_1)^\dagger\). Adding the two sides gives (26.7).

Let \(D(w)\) count all nonnegative balanced words with these source signs
and record coefficient relative to \(R_m^m\). Splitting at returns gives
\(D=1/(1-T)\). The shortest primitive word contributes w. Enclosing a
nonempty balanced word of half-length r adds the source minus sign from
its final repeated negative direction, and the record factor
\(K_mR_m^rH_m=(-1)^rR_m^{r+1}\).
Hence \(T=w-w(D(-w)-1)\), proving the first equation. Substitute its
version at -w and clear the unit-constant denominators to get the quadratic.
Its coefficient linear in the current unknown coefficient of T is -1,
so coefficients are unique. Squaring or substituting (26.9) verifies the
same equation and initial term, establishing the parity relation. All
operations at a formal degree involve finite native word sums.

### R26.6 — The return series has native completion bounds and a unitary coherent value

The count coefficients satisfy

\[
|q_k|\le\frac8{(k+1)^{3/2}}\quad(k\ge1).
\tag{26.10}
\]

The returned-energy coefficients are
\(w_2=1/2\) and
\(w_{4k}=w_{4k+2}=2q_k^2\) for \(k\ge1\).
All other positive-time coefficients vanish. Their completed sum \(p\)
is independent of source and record preparation and satisfies

\[
p_M=\tfrac12+4\sum_{k=1}^M q_k^2,\quad
0\le p-p_M\le\frac{128}{(M+1)^2},\qquad
\tfrac34<p<\tfrac45.
\tag{26.11}
\]

The coherently summed first-return map, a different target, is

\[
\widehat{\mathcal A}_*
 =\frac{1-\sqrt3}{4}\,B_s\otimes I
  +\frac{1+\sqrt3}{4}\,D_s\otimes R_m,
\qquad
\widehat{\mathcal A}_*^\dagger\widehat{\mathcal A}_*=I.
\tag{26.12}
\]

Its unitary value does not imply that all source norm has returned:
retained arrival energy is p, strictly below one. Coherent folding of
arrival marks is not declared to be a norm-preserving erasure of those marks.

**Proof.** The following coefficient bound is derived entirely on the
native scalar completion. Let
\(\lambda=(1+\sqrt3\,\iota)/2\). Its native norm and that of
\(1-\lambda\) are one, and \(\lambda+\lambda^\dagger=1\).
The already counted balanced words give
\[
\sqrt{1-y}=1-\sum_{n\ge1}\beta_n y^n,\quad
\beta_n=\frac{c_{n-1}}{2\,4^{n-1}},\quad
0<\beta_n\le\frac1{2n^{3/2}},\quad \sum_{n\ge1}\beta_n=1.
\]
The square identity follows from the Catalan count recurrence of R24;
the bound follows from its central-count bound. Absolute completion at
y=1 then gives the last equality by the square identity and nonnegative
partial sums. The beta coefficients decrease.

Put \(b_0=1,b_n=-\beta_n\). Multiplying the two absolutely convergent
series with arguments \(\lambda x,\lambda^\dagger x\) gives
\(S(x)=\sqrt{1-x+x^2}=\sum s_nx^n\). At real x, this product is a
native scalar norm, fixing the nonnegative root; in particular S(1)=1.
Thus \(q_n=\tfrac12\sum_{j>n}s_j\).
The finite geometric identity gives
\(|\sum_{j=r}^t\lambda^j|\le2\). Finite summation by parts and
monotonic beta therefore give
\(|\sum_{j\ge r}\beta_j\lambda^j|\le2\beta_r\).
No classical analytic or oscillatory-sum theorem is assumed here.

For \(h=\lfloor n/2\rfloor\), split the double coefficient tail
\(\sum_{r+s>n}b_rb_s\lambda^r(\lambda^\dagger)^s\) into
\(r\le h\), \(s\le h\), and \(r,s>h\). These sets are disjoint;
the third is a full tail rectangle. Since \(\sum|b_r|=2\), the resulting
bound is
\[
|q_n|\le4\beta_{n-h+1}+2\beta_{h+1}^2
 \le\frac{5\sqrt2}{(n+1)^{3/2}}
 <\frac8{(n+1)^{3/2}}.
\]
This proves (26.10) and absolute completion of Q at both 1 and -1.

Each source matrix in (26.7) has squared norm matrix \(2I\); record
powers preserve pairing. Squaring its coefficients and using
\(t_{2k+1}=-2t_{2k}\) proves the displayed equal return weights.
Sum (26.10) squared and use the native telescoping bound
\(\sum_{k>M}(k+1)^{-3}\le1/(2(M+1)^2)\) to obtain (26.11)'s tail.
The first three nonzero return weights sum to 3/4, and later ones are
positive. For a reproducible strict upper bound, native coefficient
comparison in (26.9) gives
\[
a_0=0,\ a_1=a_2=1,\qquad
2k a_k=4(4k-3)a_{k-1}-16(4k-9)a_{k-2}+64(2k-6)a_{k-3}\quad(k\ge3).
\tag{26.13}
\]
One may obtain it by comparing coefficients of
\(2(1-x+x^2)S'=(-1+2x)S\); the prime is formal index multiplication.
The certificate provides the exact integer
\(Z_{128}=\sum_{k=1}^{128}a_k^2 16^{128-k}\) and verifies
\(2000Z_{128}<141\,16^{128}\). This is the finite rational inequality
\(p_{128}<391/500\). Together with \(128/129^2<9/1000\), it proves
\(p<791/1000<4/5\). This finite arithmetic witness supplements, rather
than replaces, the all-depth tail proof.

Finally, the even and odd source/record sums in (26.7) are
\(E=Q(-1)=(1-\sqrt3)/4\) and \(O=1/2-E\). The native scalar norm
construction above selects \(S(-1)=\sqrt3\). Absolute convergence
justifies these substitutions and proves (26.12). Its cross terms cancel
because \(B_s^\dagger D_s=D_s^\dagger B_s=2K_s\) and
\(R_m^\dagger=-R_m\); its diagonal terms give \(2(E^2+O^2)=1\).

### R26.7 — Return parity fixes one feedback coefficient exactly

Now retain origin-arrival flags and apply R25's separate system-controlled
K on the moving source after each return, while keeping the same reused
H/K memory. Define \(J=B_s/\sqrt2\), \(D=D_s/\sqrt2\).
Even half-length returns have source pulse branches \(J,C_s\); odd
half-length returns have branches \(D,J^\dagger\).

Their source cross-continuations are

\[
\mathscr L_e(X)=J^\dagger XC_s,\qquad
\mathscr L_o(X)=DXJ^\dagger.
\tag{26.14}
\]

Put \(p_e=\sum_{k\ge1}w_{4k}\),
\(p_o=w_2+\sum_{k\ge1}w_{4k+2}\). Then

\[
\boxed{p_e-p_o=-\tfrac12,\quad p_e+p_o=p.}
\tag{26.15}
\]

The completed one-return continuation \(\mathscr M=p_e\mathscr L_e+
p_o\mathscr L_o\) acts as

\[
I\mapsto-\tfrac12H_s,\quad H_s\mapsto\tfrac12R_s,\quad
R_s\mapsto pK_s,\quad K_s\mapsto pI,
\qquad \mathscr M^4=-\tfrac{p^2}{4}\mathrm{id}.
\tag{26.16}
\]

The finite-event source response has the exact recurrence

\[
D_0=I,\qquad D_N=s_NI+
 \sum_{2m\le N}w_{2m}\mathscr L_{m\bmod2}(D_{N-2m}),
\quad s_N=1-\sum_{n\le N}w_n.
\tag{26.17}
\]

On the joint source/record preparation this response is \(D_N\otimes I\).
Thus the reused record's prepared state does not supply an adjustable
coupling. Return parity, rather than just return count, is required before
the completed sum is taken.

**Proof.** Factor each first-return map using (26.7). On system-1's branch,
\(K_sJ=C_s\) and \(K_sD=J^\dagger\). Its record factor \(R_m^m\)
is common to both branches. Direct native multiplication gives
\[
\mathscr L_e:(I,H_s,R_s,K_s)\mapsto(H_s,-R_s,K_s,I),\qquad
\mathscr L_o:(I,H_s,R_s,K_s)\mapsto(-H_s,R_s,K_s,I).
\]
Equality of the paired return weights in R26.6 leaves exactly the unpaired
weight \(w_2=1/2\), proving (26.15) and then (26.16).
Splitting at the first return, with all arrival flags retained, proves
(26.17) exactly as a native matching sum. The live first-return complement
has squared norm matrix \(s_NI\) by finite telescoping of the source
isometry. Induction shows that every cross response has record factor I:
the common record power cancels on both sides of that factor. No scalar
reset, independence law or chosen coupling coefficient is used.

### R26.8 — Completed feedback yields a fixed channel ratio and certified limits

The finite response (26.17) converges to

\[
\boxed{D_*=
 \frac{1-p}{1+p^2/4}
 \left(I-\tfrac12H_s-\tfrac14R_s-\tfrac p4K_s\right).}
\tag{26.18}
\]

This is the unique solution of \(D=(1-p)I+\mathscr M(D)\).
For any nonzero real source role v, put
\(h=B(v,H_sv)/\|v\|^2\), \(j=B(v,K_sv)/\|v\|^2\),
and \(q=(1-p)/(1+p^2/4)\). Its completed cross retention is
\(q(1-h/2-pj/4)\), independently of the auxiliary record state.
In particular,

\[
\kappa_*(e_0)=q/2,\quad \kappa_*(e_1)=3q/2,\quad
\boxed{\kappa_*(e_1)/\kappa_*(e_0)=3},\quad
\kappa_*(C_se_0)=q(1-p/4).
\tag{26.19}
\]

The ratio is exact and contains no fitted return-tail value. These are
native response channels, not identified particle charges.

Any two pairing-preserving terminal continuations after k completed
returns change a normalized source cross response by at most \(2p^k\).
For event-horizon control, if \(N\ge k(4M+2)\), then

\[
|\kappa_N(v)-\kappa_*(v)|
\le2\left[p^{k+1}+\frac{128}{(M+1)^2(1-p)}\right].
\tag{26.20}
\]

**Proof.** Every normalized return pulse preserves pairing. The norm
weight of a gap list is the product of its derived \(w_n\)'s, regardless
of pulse parity. Summing lists of length k gives \(p^k\). The completed
ledger has terminal remainder \(1-p\), so its response is
\((1-p)\sum_{k\ge0}\mathscr M^k(I)\). The bilinear tail is bounded
by the scalar geometric tail because each ordered pulse product has norm
one. Equation (26.16) sums that series explicitly and gives (26.18).
If two fixed-point solutions differ by X, applying \(\mathscr M\) four
times gives \(X=-p^2X/4\); positivity forces X=0.

The skew R term has zero final real quadratic pairing. Substitute the
three prepared source roles to get (26.19). Since \(p<1\), q is positive
and the ratio is well-defined. Its common normalization cancels exactly.
For arbitrary terminal actions on retained histories, use the joint norm
weight \(p^k\) and native pairing inequality on each controller branch;
their difference is at most \(2p^k\). This bound does not require a
terminal action to be the same on different gap lists.

To compare finite event N with the completed ledger, lists with at least
\(k+1\) returns have weight \(p^{k+1}\). Among shorter lists, a return
after N requires a gap longer than \(4M+2\). Summing the weight of a
first such gap over positions is at most
\((p-p_M)\sum_{j\ge0}p^j\). Elsewhere the finite and completed ordered
pulse lists agree. Their bilinear difference is at most twice that weight,
which proves (26.20). Thus actual event evolution selects the same completed
response; algebraic solvability alone is not substituted for convergence.

### R26.9 — Copying source letters has its own derived normalization constraint

Direction-controlled memory above must not be confused with a simultaneous
copy of the original H/K letter. Let
\(\widehat H=H_s\otimes H_m\), \(\widehat K=K_s\otimes K_m\).
These two doubled letters commute. The unflagged coherent lift
\((\widehat H+\widehat K)/\sqrt2\) is not pairing-preserving and no
scalar rescaling can repair its nonzero kernel.

More generally use normalized native flags \(f_H,f_K\), with
\(\zeta=B(f_H,f_K)\), and construct
\[
\mathcal Fv=\bigl(\widehat Hv\otimes f_H+
                  \widehat Kv\otimes f_K\bigr)/\sqrt2.
\]
Its exact normalization rule is

\[
\boxed{\mathcal F^\dagger\mathcal F
 =I+\operatorname{Re}_{\rm cut}(\zeta)R_s\otimes R_m,\qquad
\mathcal F\text{ preserves pairing iff }
\operatorname{Re}_{\rm cut}(\zeta)=0.}
\tag{26.21}
\]

For real flags this forces orthogonality. On the native turn-valued scalar
completion, a derived iota phase also satisfies it: \(f_K=\iota f_H\).
Those flags occupy one scalar direction and do not distinguish the two
letters through their matching energies. Orthogonal flags do distinguish
them. Thus normalization derives a constraint, not a unique physical record
choice, and a real-sector conclusion is not promoted to all UGD targets.

**Proof.** Two source sign reversals cancel in the doubled tensor:
\(\widehat H\widehat K=\widehat K\widehat H=R_s\otimes R_m\).
Each doubled letter is a pairing-preserving involution. Expand the full
flag pairing; the two cross terms contribute
\((\zeta+\zeta^\dagger)R_s\otimes R_m/2\), proving (26.21).
The last matrix is an invertible involution, so its coefficient must
vanish. For identical real flags it is one. The vector
\(e_{00}-e_{11}\) is then in the lift's kernel, since both doubled letters
send it to opposite values. A nonzero kernel survives every scalar rescaling.
Orthogonal native role flags give zero overlap. Alternatively the already
derived \(\iota^\dagger=-\iota\), \(\iota^2=-1\) give unit norm and
zero cut-real overlap for an iota-shifted flag. These distinct repairs are
derived by native expansion; no copying, measurement or classical complex
axiom is invoked.

## Certification and scope

The [ledger](../04-operator-evolution/R26_DERIVATION_LEDGER.json) binds the
nine written proofs to native sources and exact checks. The
[adapter](../04-operator-evolution/native_memory_residue.cjs) uses the
unchanged RKF arithmetic/replayer. It independently checks full source/record
propagation, literal H/K histories, flat-record equivalence, the commutator
response, gap identities, first-return grammar, integer/tail witnesses,
feedback and flag normalization. The
[wrapper](../04-operator-evolution/verify_r26.py) replays R25 and its frozen
earlier chain and preserves all 208 prior non-navigation files.

The [verification](../04-operator-evolution/R26_VERIFICATION.json) records
12 exact check groups, 16 canonical native identity replays, 13 rejected
false alternatives and 9 invalid dependency-graph mutations. Checks include
full first-return operators through 64 events, 4,080 literal H/K prefixes,
4,096 coefficient bounds and 1,224 terminal bilinear bounds. These finite
checks accompany the written all-depth proofs. The outward enclosures are:

| Completed native quantity | Lower bound | Upper bound |
| --- | --- | --- |
| Retained first-return energy p | 0.782004429700 | 0.782012055372 |
| Feedback retention for e0 | 0.094540136129 | 0.094543687856 |
| Feedback retention for e1 | 0.283620408389 | 0.283631063567 |
| Feedback retention for Ce0 | 0.152114509174 | 0.152120584360 |

The exact e1/e0 ratio is three before any numerical enclosure. For a
concrete finite-event witness, (26.20) with k=64 and M=4096 gives an error
below 0.000070193186 once N is at least 1,048,704. This evaluates the proved
tail bound; it is not direct enumeration of a million-event history tree.

The provenance gate audits declared dependencies and bound proof text.
Finite checks are not a proof assistant, and neither is the graph gate.
All-depth/completion claims have the written proofs above. No upstream
whole-engine PASS or new physical validation is claimed.

| Source at exact commit | Premise status here |
| --- | --- |
| [R16 C1, C3–C8, C13](https://github.com/Parveen117/extra-ideas/blob/4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff/02-relational-response/emk-topology-foundation/emk_topology_foundation.tex) | Native roles, signs, matching pairing, derived iota and scalar completion. Equality, finite counting and induction remain explicit metamathematical infrastructure. |
| [R18.1–R18.4](https://github.com/Parveen117/extra-ideas/blob/7afb4f745fc4ac64b9eaa4fe2939b55ed22d3699/02-relational-response/CUT_TRANSPORT_METRIC_R18.md) | Native transport on a constructed address target; no physical metric or clock is selected. |
| [R20.1–R20.4, R20.7](https://github.com/Parveen117/extra-ideas/blob/c6d1810114129b6aa74addd05cceb9083de27dd5/02-relational-response/NATIVE_RECORD_INTERACTION_R20.md) | Native tuple construction, controlled records and matching readout; allocation choices remain constructions. |
| [R24.1–R24.8](https://github.com/Parveen117/extra-ideas/blob/d660bbf026e1de818fcb3dd61d46346a4ad618ae/02-relational-response/NATIVE_RETURN_COUPLING_R24.md) | Native first-return splitting, finite counts and scalar tails. New H/K-memory coefficients are derived here rather than equated with R24's eta. |
| [R25.1–R25.8](https://github.com/Parveen117/extra-ideas/blob/d4f1ae2481a72a04c465f1c17b37d154a6122eff/02-relational-response/NATIVE_RETURN_FEEDBACK_R25.md) | Native retained-arrival feedback, bilinear completion and boundary control. R26 retains parity-dependent pulses and derives their new weights explicitly. |
| [RKF canonical algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md) | Unchanged native equality/replay engine; scoped certificate use. |
| [Publications AS-1–AS-5](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-alpha-selection/THEOREM.md) | **Comparison only, mixed premises:** supplied positive cell weights; AS-3 admits a positive quadratic physical adapter. That imported adapter and electromagnetic normalization are excluded from R26 proof paths. Research-branch commit. |

The new selection data are concrete: flat memory classes are response-
equivalent, a two-step commutator square detects nonflat order, and the
native anticommuting choice fixes a propagation gap and channel ratio.
Physical selection of that choice and its mass/charge interpretation remain
open. The nine results do not alter any earlier RH or Yang–Mills boundary.
