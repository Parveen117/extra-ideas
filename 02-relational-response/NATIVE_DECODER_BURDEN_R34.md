# R34: native decoder burden and a balanced observer reduction

Research owner: Monty Dabas. Development: 2 October 2026.

The user's beta question concerns the burden of a design. A repository alone
does not have one mathematical beta: a target, observer and output pairing
must be named. The exact `1/4` example located in canonical RKF T43 section
11 uses a four-layer observer equal to the identity and the target
\(L=e_4^\dagger/2\). It is that example's squared decoder norm, not a
repository-wide RH constant or a fraction of unfinished proofs. T43's
abstract chapter admits Hilbert spaces and a bounded generator; it is a
comparison citation, excluded from the native proofs below. No claim about
the unresolved RH endpoint is changed.

R34 derives this burden directly from the native matching pairing and
R32--R33's source. For the existing eight-event increment, the sharp
full-state burden is **1**, and a unit local-role target has burden **1/2**.
A balanced bank of eight- and sixteen-event increments reduces these to
**1/2** and **0.36163703165179...**, respectively. The balance is uniquely
optimal for the worst-state burden within the declared two-increment family.
It is not selected by fitting a desired physical constant.

The output's total squared branch weight stays one. The longer branch costs
more source events. A separate cost calculation and the exact inherited-error
map below prevent gain, duplication or postprocessing from being presented
as free extra information. Physical noise, clocks and electromagnetic alpha
are not supplied.

## Native data and the design contract

Use the native real signed source/record sector and matching completion in
R16, R20, R32 and R33. Let \(B\) be R32's pairing-preserving eight-event
continuation, \(F=B-I\), and

\[
G=F^\dagger F=Q_xQ_y,\quad Q_i=I+D_i^\dagger D_i/4,
\quad D_i=T_i^2-I,\quad I\le G\le4I,
\quad B+B^\dagger=2I-G.
\tag{34.0}
\]

Unless a finite cycle is explicitly named, the addresses are the unwrapped
count plane with the native Cauchy completion. Set \(S_i=T_i^2\),
\(E_i=S_i+I\), \(C=G-I\), and \(D=4I-G\); the capital D without an
axis index is different from \(D_i\). Direct native difference expansion
gives

\[
C=\tfrac14D_x^\dagger D_x+\tfrac14D_y^\dagger D_y
 +\tfrac1{16}(D_xD_y)^\dagger(D_xD_y),
\]
\[
D=\tfrac12E_x^\dagger E_x+\tfrac14E_y^\dagger E_y
 +\tfrac1{16}(D_xE_y)^\dagger(D_xE_y).
\tag{34.1}
\]

The scalar-shift factors commute. Their product \(CD\) is a sum of the
nine pairwise product squares in (34.1). These identities will prove all
interval bounds without an imported spectral or minimax theorem.

## Written results

### R34.1 — Native pairing derives the target burden and its minimum decoder

For a specified observer \(A\) and scalar target \(L\), define

\[
\beta_A(L)=\sup_{A\psi\ne0}
 \frac{|L\psi|^2}{\|A\psi\|^2},
\tag{34.2}
\]

with value infinity when L does not vanish on the observer kernel. Suppose
\(M=A^\dagger A\ge mI\), \(m>0\), and its native inverse exists.
For a target written with its native coefficient field,
\(L\psi=\langle l,\psi\rangle\), then

\[
\boxed{\beta_A(L)=\langle l,M^{-1}l\rangle,\qquad
 c_L=AM^{-1}l.}
\tag{34.3}
\]

The observer-output vector \(c_L\) decodes L and has the least squared
norm among all such decoders. For the full source target define separately

\[
\beta_A^{\rm full}=\sup_{\psi\ne0}
 \frac{\|\psi\|^2}{\|A\psi\|^2}.
\tag{34.4}
\]

It equals the reciprocal of the sharp lower squared-norm bound of A.

**Proof.** Native dagger matching gives
\(\langle AM^{-1}l,A\psi\rangle=\langle l,\psi\rangle\), and
\(\|AM^{-1}l\|^2=\langle l,M^{-1}l\rangle\). The native
sum-of-squares inequality supplies the upper bound in (34.2). Equality is
attained at \(\psi=M^{-1}l\) for nonzero L. If \(c'\) is another
decoder, \(A^\dagger(c'-c_L)=0\), so its difference is orthogonal to
\(c_L\) and \(\|c'\|^2=\|c_L\|^2+\|c'-c_L\|^2\).
Finite matching proves every identity; bounded completion passes them to
the already constructed native limits. Formula (34.4) is the reciprocal
ratio itself; approaching the sharp lower bound proves equality. No
classical Hilbert or Riesz representation theorem is a premise: l is the
specified native target field, not an arbitrary unrepresented functional.

### R34.2 — The current increment has burden one globally and one half locally

For the existing observer \(A_0=F\),

\[
\boxed{\beta_F^{\rm full}=1.}
\tag{34.5}
\]

For any address p and unit four-role vector v, let
\(l=E_pv\) be its single-cut target. Then

\[
\boxed{\beta_F(L)=\tfrac12,\qquad
 c_L=FG^{-1}E_pv.}
\tag{34.6}
\]

Thus R33's \(\rho^2\), its spread fraction \(1/9\), and this local
burden \(1/2\) are three different defined quantities.

**Proof.** R32 proves the sharp increment bounds one and four. More
explicitly, normalized constant packets on M by M even-address squares
have Gram ratio \((1+1/(2M))^2\) tending to one. This proves (34.5).
R32's inverse kernel has diagonal \(I/2\), so (34.3) gives (34.6),
including attainment by its localized completed field.

These values are for the unwrapped target. On a component cycle of length
L, write \(Q_L=(3/2)I-(S+S^\dagger)/4\). The local burden is instead
\(((Q_L^{-1})_{00})^2\). For example L=2 gives \(9/16\), not
\(1/2\). The full-state value one is still attained by the uniform
cyclic field. Boundary choice is part of the observation contract.

### R34.3 — Balanced native increments uniquely minimize the two-horizon worst burden

Let \(F_2=B^2-I=(B+I)F\). Split the source into two retained native tags
with squared weights \(1-s,s\), where \(0\le s\le1\), and declare
the observer

\[
A_s\psi=(\sqrt{1-s}\,F\psi,\sqrt{s}\,F_2\psi).
\tag{34.7}
\]

The sum of squared branch weights is one. Native positive cut roots realize
the coefficients. At the selected balance they are the existing equal-tag
normalization \(1/\sqrt2\). Then

\[
F_2^\dagger F_2=G(4I-G),\qquad
M_s=A_s^\dagger A_s=G[(1+3s)I-sG],
\tag{34.8}
\]

and its sharp lower bound is

\[
m(s)=\min\{1+2s,\,4-4s\},\quad
\beta_{A_s}^{\rm full}=\frac1{m(s)}\ (s<1).
\tag{34.9}
\]

For \(s=1\) the full-state burden is infinite. The unique best balance
is \(s=1/2\), which gives

\[
\boxed{M_* =\tfrac12G(5I-G)
 =2I+\tfrac12(G-I)(4I-G),\qquad
 \beta_{A_*}^{\rm full}=\tfrac12.}
\tag{34.10}
\]

**Proof.** B commutes with its inverse and with G. Expand
\((B^2-I)^\dagger(B^2-I)=4I-(B+B^\dagger)^2\), then insert (34.0)
to obtain (34.8). Its exact chord-plus-residue identity is

\[
M_s=\frac{4I-G}{3}(1+2s)+\frac{G-I}{3}(4-4s)
 +s(G-I)(4I-G).
\tag{34.11}
\]

Equations (34.1) make all three factors nonnegative and the first two
weights sum to identity. Hence \(M_s\ge m(s)I\).

For sharpness, use constant and alternating \((-1)^{a+b}\) packets on
the component-count square. Away from its boundary G acts as I or 4I.
The norm of \((G-gI)\psi_M\) tends to zero for g=1 or 4: for M at
least two its squared norm is bounded by \(512/M\), since at most 8M
boundary sites remain and every residual coefficient is at most \(8/M\).
A finite polynomial in G therefore has the corresponding endpoint limit.
This proves the two sharp values and the infinite burden at s=1. On even
component cycles the two endpoint fields are exact.

For s below one half the first endpoint is the minimum and increases with
s; above one half the second is the minimum and decreases. They meet only
at one half, with common value two. This proves unique optimality within
this family and objective without an imported optimizer. It is not a claim
of optimality over every native observer or physical measurement.

### R34.4 — Native path counts certify the improved local burden

For the same unit point target as R34.2,

\[
\beta_*^{\rm point}=\frac15+\frac2{25}
 \sum_{n\ge0}\frac{m_n^2}{5^n},\qquad
m_n=4^{-n}\sum_{k=0}^{\lfloor n/2\rfloor}
 \binom n{2k}\binom{2k}k6^{n-2k}.
\tag{34.12}
\]

After terms zero through N the missing burden obeys

\[
0\le\beta_*^{\rm point}-\beta_{*,N}^{\rm point}
 \le\frac25(4/5)^{N+1}.
\tag{34.13}
\]

The N=160 exact native count certificate gives

\[
\boxed{0.3616370316517928
 <\beta_*^{\rm point}<0.3616370316517930<\tfrac12.}
\tag{34.14}
\]

This is a reduction greater than 27.67 percent under the declared matching
output norm. It holds for every unit role at every address. The balanced
bank was optimized for the full-state burden; no claim that it minimizes
this different point-target objective is made.

**Proof.** Direct multiplication gives

\[
M_*^{-1}=\frac25\{G^{-1}+(5I-G)^{-1}\},\qquad
(5I-G)^{-1}=\frac15\sum_{n\ge0}(G/5)^n.
\tag{34.15}
\]

The second inverse exists by finite geometric cancellation and the native
norm bound \(\|G/5\|\le4/5\). Its norm remainder is at most
\((4/5)^{N+1}\). The diagonal of \(G^{-1}\) is one half. Since
\(G=Q_xQ_y\), the origin coefficient of \(G^n\) is the square of
the origin coefficient of \(Q^n\). In the component shift,
\(Q=(6I-S-S^{-1})/4\). A returning length-n word has k positive and k
negative shifts. Choose their 2k positions and their k positive positions;
the remaining positions contribute \(6^{n-2k}\). The 2k minus signs
cancel. This counts exactly the positive finite sum in (34.12), without
Fourier integration, a classical probability measure or spectral inputs.
Every omitted diagonal term is nonnegative; the operator tail supplies
(34.13). Substitution of finite integer counts and rational comparisons
proves (34.14). Translation and the scalar four-role identity prove the
stated preparation invariance.

### R34.5 — Gain and event-cost audits separate a real design change from rescaling

The exact squared-norm ranges of the baseline and balanced banks are

\[
\|F\psi\|^2/\|\psi\|^2\in[1,4],\qquad
\|A_*\psi\|^2/\|\psi\|^2\in[2,25/8],
\tag{34.16}
\]

with all endpoints sharp. Their ratio of upper to lower response, which
cancels any common observer gain, improves from **4** to **25/16**.
Thus the design change is more than scalar amplification.

Nevertheless costs depend on their declared meaning. The branch-weighted
event count in (34.7) is \(8(1+s)\). If output is normalized by that count
relative to the baseline eight events, the full-state burden becomes

\[
\widehat\beta_s^{\rm full}=\frac{1+s}{m(s)}.
\tag{34.17}
\]

Its unique minimum is again at one half, with value **3/4** rather than one.
The similarly cost-normalized local burden is
\(3\beta_*^{\rm point}/2>1/2\), so its improvement does not survive
this cost convention. If the cost is maximum event horizon instead, the
balanced bank needs sixteen events and its corresponding normalized
full-state burden is \(2(1/2)=1\). No universal resource improvement is
asserted.

**Proof.** The lower bound is R34.3. Complete one scalar operator square:

\[
M_*=\frac{25}{8}I-\frac12(G-\tfrac52I)^2.
\tag{34.18}
\]

This gives the upper bound without a spectral theorem. It is sharp by a
finite native count pattern. On a component six-cycle, the values
\((2,1,-1,-2,-1,1)\) obey \((S+S^{-1})u=u\), hence
\(Q_xu=5u/4\). An alternating second-direction pattern has Q_y value
two. Their product is an exact G value \(5/2\) and therefore an M_*
value \(25/8\). Truncating larger repeated rectangles loses only boundary
sites, proving the same sharp supremum in the unwrapped completion. The
baseline upper bound four has the alternating packet witness from R32.

Under a scalar gain a, both endpoints multiply by \(a^2\), leaving
their ratio unchanged. Dividing the bank by \(\sqrt{1+s}\) proves
(34.17). For \(s\le1/2\), the ratio \((1+s)/(1+2s)\) strictly
decreases by direct cross multiplication; for \(s\ge1/2\),
\((1+s)/(4-4s)\) strictly increases. This proves the optimum. The
point-target and maximum-horizon conclusions follow by the same declared
normalization and (34.14). The branch average is a native weighted count,
not a derived physical clock or universal experimental cost model.

### R34.6 — The improved bank preserves the information and exposes inherited errors

The balanced output factors exactly through the old signed increment:

\[
A_*=JF,\quad Jd=\frac1{\sqrt2}(d,(B+I)d),\quad
J^\dagger J=\frac12(5I-G).
\tag{34.19}
\]

Therefore it has the same recognition kernel as F. If it is generated only
by postprocessing an already measured \(d+e\), its error is Je, with the
same shared e in both branches. Equip its range with the transported norm
\(\|Jd\|_{\rm tr}=\|d\|\); every burden is then exactly the old one.
The reduction in R34.3--R34.4 refers to the declared native product matching
norm on the two-tag bank, not to fictitious independent errors introduced
by duplicating existing data.

More generally,

\[
\beta_{aA}(L)=\beta_A(L)/a^2,\qquad
\beta_A(bL)=b^2\beta_A(L).
\tag{34.20}
\]

Unnormalized m-fold duplicate readout divides raw beta by m; the normalized
bank \((A,\ldots,A)/\sqrt m\) leaves it unchanged. A numerical comparison
with the T43 quarter therefore requires matching both target and observer
normalizations. In that example changing the target from \(e_4/2\) to
\(e_4\) changes beta from one quarter to one.

**Proof.** The first factorization follows from \(F_2=(B+I)F\).
Expand \((B+I)^\dagger(B+I)=4I-G\) to prove the Gram identity for J.
Its first component proves injectivity. Substitution of Je and the defined
transported norm proves the inherited-error and unchanged-burden claims
without any stochastic noise premise. Each identity in (34.20) follows by
substitution in the ratio (34.2). The duplicate bank has Gram
\(mA^\dagger A\), or \(A^\dagger A\) after normalization. Finally the
T43 example has identity Gram and target norm squared one quarter; this is
an arithmetic explanation of its cited value, not an imported theorem.

The present coefficient decrease is a certified mathematical design result
for a declared two-horizon observation. It neither decreases the number of
unproved physical selection obligations nor derives physical alpha. The
physical observation, its error carrier and the cost contract still need
their own native identification.

## Certification and source status

The [ledger](../04-operator-evolution/R34_DERIVATION_LEDGER.json) binds the
six proofs, their native source chain and scoped exact computations. The
[application](../04-operator-evolution/native_decoder_burden.cjs) uses the
unchanged canonical native arithmetic and word replayer. The
[verification record](../04-operator-evolution/R34_VERIFICATION.json)
includes complete Laurent identities, cyclic native inverses, explicit
sharp patterns, exact finite count tails and false-normalization controls.
It replays the frozen R33 chain and preserves all 264 earlier non-navigation
files. These checks and proof-source hashes are not a formal proof-assistant
verification of all prose or physical validation of the observation.

* **Native derived sources:** R16 C1, C3--C8, C13 at extra-ideas commit
  `4cc46d138c6e0610865ac3d92056cf9e5f7eb7ff`, and R20 tuple matching at
  `c6d1810114129b6aa74addd05cceb9083de27dd5`. The signed/refined sector
  and its completion remain explicit; no classical Hilbert primitive.
* **Native derived continuation with constructed interface:**
  [R32.4, R32.6--R32.7](NATIVE_RETAINED_LOOP_GAP_R32.md), commit
  `1105a43209b144bd34a81bfbb99acc7876948e95`: G, its sharp bounds and
  inverse follow from the retained H/K event construction. Its physical
  interpretation is not promoted to a premise.
* **Native cut-response source:** [R33.1--R33.2, R33.5--R33.7](NATIVE_INTERACTION_MEMORY_R33.md),
  commit `5ddce4544de7111a80c93a8ade7cf2399c7c90f3`: native inverse
  elimination, source-scale/readout distinctions and completed kernel.
* **Native implementation:** Recognition-Kernel-Framework commit
  `3cc5a33b05c16d59c90994ddda69dedc0d392424`, unchanged canonical
  `operator_foundation/core/native_operator.cjs`, `core/paninian_operator.cjs`
  and `core/workbench.cjs`; exact files retained in the source pins.
* **Comparison only, admitted classical representation, excluded from proof:**
  the same RKF commit,
  [`theorum/43_bilateral_jet_flow_recognition_capstone_theorem.md`](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/43_bilateral_jet_flow_recognition_capstone_theorem.md),
  sections 2, 8 and 11. The abstract statement begins with Hilbert spaces
  and a bounded generator; its quarter is a specified finite target example.
  R34 reconstructs the needed burden and decoder identities directly from
  native matching. Neither this citation nor its local finite PASS closes RH.
