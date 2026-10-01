# R14 — Native source, phase, event content, memory and recognition geometry

Author: Monty Dabas. Development: 1 October 2026.

## Result and exact scope

R14 moves R13's foundation back to **oriented cut arithmetic**. Dagger derives the radial/complementary split. Multiplication by a native source scalar derives the transport, including its radial gain. Native polar/logarithm results derive its phase. Refining the resulting paths derives event content. Eliminating the complementary coordinate derives memory. A finite source-refinement ledger derives the unresolved-force covariance. Finally, the native logarithm derives a local contrast geometry directly from that same source.

The central formulas, for a nonzero native scalar `z=a+iota b`, are

\[
d=N_\Sigma(z)=a^2+b^2,\quad
q=\frac{b^2}{d},\quad r=\frac{a^2}{d},\quad
K_\ell=-b^2a^\ell,
\]
\[
\operatorname{Cov}_\nu(\xi_n,\xi_m)=b^2a^{n+m}S_\nu,
\qquad
g_z(\delta z,\delta z)
=\frac{4(a\,\delta b-b\,\delta a)^2}{(a^2+b^2)^2}.
\]

These are **mathematical derivations inside a declared native source class**. There is no independently chosen noise variance, angle, stopping probability, Hilbert metric or Gaussian likelihood in these formulas. The actual source scalar and initial refinement ledger remain source data. Repackaging those data as one tuple would not prove that the primitive framework uniquely selects them. Section 10 proves precisely why that stronger conclusion does not follow from the algebra presently cited.

The native event object and its content are constructed here. A physical detector selecting one actual outcome with those frequencies is a further statement, not a synonym for content additivity. This distinction keeps the requested single-source programme testable.

## 1. Starting object and provenance

Use signed cut counts, their rational refinements, and the oriented two-component cut construction. A quarter-turn acts on a pair by

\[
R(a,b)=(-b,a),\qquad R^2(a,b)=(-a,-b).
\]

Writing the operator as `iota` gives `iota^2=-1`. Composing `aI+bR` and `cI+eR` gives

\[
(a+\iota b)(c+\iota e)=(ac-be)+\iota(ae+bc).
\]

Reversing turn orientation defines dagger, so

\[
(a+\iota b)^\dagger=a-\iota b,
\qquad N_\Sigma(z)=z^\dagger z=a^2+b^2.
\]

This is the existing F00 cut-pair construction, not a new claim of priority. The decision to construct an **oriented two-component cut carrier** is explicit; pure logic alone does not single it out. Its order/completion is the established input to F00-E. Finite statements below use its exact rational subfield; analytic statements use the completed radial field.

References travel with their input contracts: [F00/F00-E — native cut-pair arithmetic and finite audit; completed oriented field is the Euler theorem's input][F00], [F00-G — native logarithm; complete Archimedean ordered radial field][LOG], and [RH T01 — native path/recognition completion; typed paths and declared gauges][RH01]. A full origin/evidence map is in [the source ledger](DERIVATION_SOURCE_LEDGER_R14.md).

## 2. The split and positive geometry are derived from dagger

Regard the cut field as a two-dimensional module **over its radial subfield**. Set

\[
Jv=v^\dagger,\qquad P=\tfrac12(I+J),\qquad Q=\tfrac12(I-J).
\]

Then `P(a+iota b)=a`, `Q(a+iota b)=iota b`. Since `J^2=I`, expansion gives

\[
P^2=P,\quad Q^2=Q,\quad PQ=QP=0,\quad P+Q=I.
\]

These maps are radial-linear. **They are not cut-complex-linear:** `P(iota)=0` while `iota P(1)=iota`. The two real components here are not two independent complex channels. This avoids a type error in transporting the result to the complex-linear N03 catalogue.

The radial pairing is the derived polarization

\[
B(u,v)=\operatorname{rad}(u^\dagger v)
=\tfrac12\bigl(N(u+v)-N(u)-N(v)\bigr).
\]

Multiplying the cut pairs proves bilinearity, symmetry, `B(v,v)=N(v)`, and positivity. It also proves `B(Pu,Qv)=0` and

\[
N(v)=N(Pv)+N(Qv).
\tag{R14.1}
\]

Thus the finite positive pairing is earned from the cut square. No ordinary inner-product space is postulated. Positivity of this scalar pairing does not establish a Lorentzian spacetime metric.

The cut is canonical **relative to the already oriented scalar carrier**. R14 has not proved that every native path carrier or physical observer is this radial cut.

## 3. Source multiplication derives transport, scale and theta

For nonzero `z=a+iota b`, define `T_z v=zv`. Its columns in the derived basis `(1,iota)` give

\[
T_z=\begin{pmatrix}a&-b\\b&a\end{pmatrix},\qquad
N(T_zv)=dN(v),\quad d=a^2+b^2>0.
\tag{R14.2}
\]

This derives the transport and retains its gain; it does **not** assume `a^2+b^2=1`. The positive square root is obtained within the native completion: start with `[0,d+1]`, bisect, and retain the interval bracketing `t^2=d`. Its lengths tend to zero; continuity of multiplication identifies the common limit, and positivity gives uniqueness. Thus `rho=sqrt(d)` is earned and `u=z/rho` has unit norm.

[F00-J — native fundamental period; factorial completion and native order arguments][PERIOD] and [F00-K §§1–5 — native polar/logarithm theorem; same completed cut field][POLAR] give

\[
z=\rho\operatorname{Exp}_\Sigma(\iota\theta),\qquad
\boxed{\theta=\operatorname{Arg}_\Sigma(z)\in(-\pi_\Sigma,\pi_\Sigma].}
\tag{R14.3}
\]

Theta is therefore a coordinate **derived from the supplied native transition**, not an extra angle supplied alongside it. Its infinitesimal variation away from the principal branch seam is

\[
\boxed{d\theta=\frac{a\,db-b\,da}{a^2+b^2}.}
\tag{R14.4}
\]

**Proof.** Differentiate `z=rho Exp_Sigma(iota theta)` using F00-E/G's native rules. The turn part of `z^{-1}dz` is `dtheta`. Direct cut multiplication gives `(a db-b da)/d`. The radial part is `d rho/rho`. No ordinary `atan2` or trigonometric library enters. ∎

The native phase metric can already be constructed here, before defining any statistical contrast. Differentiate `u=z/rho` and use the pairing derived in §2:

\[
\boxed{h_z(\delta z,\delta z)=N(\delta u)
=\frac{(\delta a)^2+(\delta b)^2}{d}
-\frac{(a\delta a+b\delta b)^2}{d^2}
=\frac{(a\delta b-b\delta a)^2}{d^2}=d\theta^2.}
\]

Every equality is expansion of the native cut square. This is the metric induced by the source's own normalized amplitude map and native pairing. Its radial kernel follows from scale invariance. Section 8 proves that the independently constructed record-logarithm Hessian is **four times this same metric** on its positive-support domain; it is not needed to define `h`.

An explicit local angle construction is useful for verification. On the Cayley subchart `|h|<1`, put

\[
\Theta(h)=2\sum_{k\ge0}\frac{(-1)^kh^{2k+1}}{2k+1},\quad
u(h)=\frac{1+\iota h}{1-\iota h}.
\]

Geometric tails justify the derivative `Theta'=2/(1+h^2)`. Direct multiplication gives `u'=2 iota u/(1+h^2)`. Hence `Exp_Sigma(-iota Theta)u` has zero derivative and value one at zero, proving `u=Exp_Sigma(iota Theta)`. The native alternating-tail bound encloses theta using adjacent rational partial sums. The certificate checks those enclosures against independently evaluated finite Euler polynomials with explicit factorial tails. It does not claim to certify all of F00-J/K.

For a word of transitions, principal arguments alone omit an integer: the sum of their arguments differs from the principal argument of the product by `2 pi_Sigma k`. The source word retains this wrap ledger. Endpoint phase cannot reconstruct an erased word: four quarter-turns and an empty word both end at one. The phase ledger and the complementary state coordinate remain distinct kinds of Smriti.

## 4. One coefficient controls incompatibility, leakage and return

On the radial module of §2,

\[
[P,T_z]=\begin{pmatrix}0&-b\\-b&0\end{pmatrix},
\qquad [P,T_z]^2=b^2I.
\]

Moreover,

\[
dP-(PT_zP)^\dagger(PT_zP)=(QT_zP)^\dagger(QT_zP)=b^2P,
\tag{R14.5}
\]
\[
PT_z^2P-(PT_zP)^2=PT_zQT_zP=-b^2P.
\tag{R14.6}
\]

Here the operator dagger is derived from `B`; the coordinate transpose is justified by that pairing. Each identity follows by multiplication. R14.6 is the [native N06 returning-memory identity — algebraic compression; no positive pairing needed for the expansion][NATIVE]. R14.5 is also the finite radial realization of [RH F06 §7 — native compression square; complete positive nondegenerate pairing is a prerequisite][F06].

After dividing by source scale `d`, the common leakage strength is `q=b^2/d`. No noise parameter or geometric curvature tensor has been introduced. This cut/transport incompatibility is not automatically a Riemann tensor.

## 5. Memory and unresolved force without a stochastic assumption

Full coherent evolution reads

\[
x_{n+1}=ax_n-by_n,\qquad y_{n+1}=bx_n+ay_n.
\]

Induction gives `y_n=a^n y_0+b sum_(j<n) a^(n-1-j)x_j`. Substitution proves

\[
\boxed{x_{n+1}=ax_n- b^2\sum_{j<n}a^{n-1-j}x_j-ba^ny_0.}
\tag{R14.7}
\]

Thus `K_l=-b^2 a^l` and `f_n=-ba^n y_0` are derived from the same source multiplication. This holds for every finite step, including `|a|>=1`; no infinite summability claim is needed. The memory kernel is deterministic. Calling `f_n` stochastic requires a source ensemble or unresolved preparation.

Eliminating again gives

\[
\boxed{x_{n+2}=2ax_{n+1}-d x_n.}
\tag{R14.8}
\]

When `b!=0`, two radial records recover `y_0=(ax_0-x_1)/b`. Their observation rows are `(1,0)` and `(a,-b)`, with determinant `-b`; exactly one further radial scalar channel is necessary and sufficient. When `b=0`, the complementary state never returns to this observer. This is N03 on a radial module; the numerical row engine is used on radial coefficients, without silently doubling the field dimension.

The corresponding time-varying/multichannel identity remains the R13 ordered block elimination. A constant scalar source is the particular class proved and executed here, not a claim that every native evolution is scalar or stationary.

## 6. Deriving noise covariance from refinement, without imposing balance

Take an actual finite source ledger of `M` equally refined cells. On the observer fibre with common visible value `x_0`, cell `j` retains hidden value `y_j`. Equal elementary cell content is the native counting construction; repeated values may be compressed into integer multiplicities. [RH T01-E1/E3 — native refinement-count content; declared cell grid, finite additivity][COUNT] supplies this source, without Haar or a Gaussian measure.

Define finite averages by counting:

\[
\mu=\frac1M\sum_jy_j,\qquad
S_\nu=\frac1M\sum_j(y_j-\mu)^2
=\frac1{2M^2}\sum_{j,k}(y_j-y_k)^2\ge0.
\tag{R14.9}
\]

The last equality follows by expansion. Refining every cell equally or reordering cells preserves both expressions. The centered unresolved force is `xi_n(j)=-ba^n(y_j-mu)`. Consequently

\[
\boxed{\operatorname{Cov}_\nu(\xi_n,\xi_m)=b^2a^{n+m}S_\nu,\qquad
\operatorname{Cov}_\nu(\xi_l,\xi_0)=-S_\nu K_l.}
\tag{R14.10}
\]

**Proof.** Substitute the centered force and factor the fixed source coefficients out of the finite sum. ∎

No free covariance is installed. Unlike R13, neither equal sign weights nor unit initial norm is required. R13's balanced pair is recovered when the ledger actually has equal counts at `+y` and `-y`, and its formula `S=1-x_0^2` additionally uses that family's unit norm.

For example, three cells at `y=2` and one at `y=-2` derive `mu=1` and `S=3`; imposing zero mean would give the wrong source. A single retained hidden value derives `S=0`, although the return kernel can remain nonzero. Under exact two-record recovery with `b!=0`, the refined observer fibre has one hidden value, so this unresolved covariance vanishes. Noise is observer-relative unresolved returning residue in this class; this does not prove that all physical noise has that origin.

The initial cell assignments are not deduced from the transition scalar. Their selection is a remaining source question, made explicit rather than hidden in an assumed covariance.

## 7. Events as native record refinements; content as a theorem

An event record is a tagged native alternative `P` or `Q` after a transition. Retain both in a direct sum of labelled path sectors:

\[
\mathscr V_zv=(PT_zv,QT_zv).
\]

By R14.1–2 its total energy is `d N(v)`. Repeating the record construction `n` times gives path amplitudes

\[
A_hv=P_{h_n}T_z\cdots P_{h_1}T_zv,
\quad h\in\{P,Q\}^n,
\]
\[
\sum_hN(A_hv)=d^nN(v),\qquad
\sum_h A_hv=T_z^nv.
\tag{R14.11}
\]

**Proof.** For the first identity, each parent branch splits into two energies whose sum is its energy multiplied by `d`; induction gives the result. For the second, sum the final projectors to `I` and repeat backwards. ∎

Define the normalized recognition-energy content of a set of records by

\[
\nu_z(E)=\frac{\sum_{h\in E}N(A_hv)}{d^nN(v)}.
\tag{R14.12}
\]

Normalization, positivity, finite additivity and consistency under further record refinement now follow from R14.11. R13's separate assumptions about additive **norm** weights are therefore replaced by an explicit record-energy construction. The choice to use recognition energy is inherited from the source's `E_Sigma`, not a theorem that no other native content can be defined. Initial **count content** in §6 and record **energy content** here are different derived functionals; they are not silently identified.

Starting on the radial unit, the first `Q` after `k` continuations has raw amplitude `iota b a^k`. Dividing its energy by `d^(k+1)` gives

\[
\boxed{\nu_z(K=k)=qr^k,\qquad r=\frac{a^2}{d},\quad q=\frac{b^2}{d}.}
\tag{R14.13}
\]

The finite exit ledger is `sum_(k=0)^m q r^k+r^(m+1)=1`. It is exact even at the endpoints: `b=0` leaves tail one; `a=0` gives immediate exit. For `b!=0`, native geometric completion sends the tail to zero. This supplies the distribution on first-exit counts; it does not assert a general sigma-additivity theorem for every completed event algebra.

The record operation is constructed in the native path algebra. It is an available protocol, not a proof that all observations record every intermediate cut. Coherent summation of branch amplitudes and retention of their separate energy records give different readouts. For `z=3/5+4 iota/5`, two coherent steps have visible content `49/625`; the all-`P` record has `81/625`.

## 8. A local recognition metric derived using the existing native logarithm

For strictly positive normalized record contents `p_i`, define the logarithmic contrast

\[
\mathcal D(p\Vert\widetilde p)
=\sum_i p_i\operatorname{Log}_\Sigma(p_i/\widetilde p_i).
\]

This is a defined comparison functional built using F00-G. The selected functional is explicit; no claim of uniqueness among all possible contrasts is made. Put `tilde p=p+epsilon v+(epsilon^2/2)w+...`, with sums of `v,w` zero. Native log differentiation gives

\[
\mathcal D'(0)=-\sum_i v_i=0,\qquad
\boxed{\mathcal D''(0)=\sum_i\frac{v_i^2}{p_i}.}
\tag{R14.14}
\]

The second derivative term `-sum_i w_i` vanishes by normalization. Positivity and the null directions follow directly from this sum of native squares. The formula agrees with Fisher's form in a later statistical chart; no Fisher theorem is used to derive it.

Now insert the **derived** binary contents `(r,q)=(a^2/d,b^2/d)`. Direct differentiation and simplification give, for `ab!=0`,

\[
\boxed{g_z=\frac{dr^2}{r}+\frac{dq^2}{q}
=4\frac{(a\,db-b\,da)^2}{d^2}=4\,d\theta^2.}
\tag{R14.15}
\]

The raw source plane has a radial null direction: scaling `z` changes neither normalized record content nor this phase metric. On a local phase quotient it is positive. The squared readout still loses discrete sign/orientation information globally; a locally nondegenerate metric does not undo that blindness.

At `a=0` or `b=0`, the positive-support log-Hessian formula is not applicable. The source expression `4(a db-b da)^2/d^2` has an amplitude-geometric extension, but calling it the regular event Hessian at changing support would be incorrect. At `z=0`, normalized content and phase are undefined and are retained as a separate stratum.

Equation R14.15 derives phase recognition geometry. R13's two-dimensional seam metric also uses its sheet gain, register construction and tangent chart; R14 does not silently promote that more specific construction to a uniquely selected spacetime.

## 9. Exact worked source

Take the native transition `z=3/5+4 iota/5`, and the preparation ledger `{2,2,2,-2}` at one fixed visible value. Then

| Output | Derived value |
|---|---|
| Native source scale `d` | `1` |
| Principal phase | `2 sum_(k>=0) (-1)^k (1/2)^(2k+1)/(2k+1)` |
| Complementary event content `q` | `16/25` |
| Continuation content `r` | `9/25` |
| Memory kernel `K_l` | `-(16/25)(3/5)^l` |
| Hidden mean and variance from cell counts | `mu=1`, `S=3` |
| Unresolved-force covariance | `(48/25)(3/5)^(n+m)` |
| Event content at `k` | `(16/25)(9/25)^k` |
| Local angular contrast metric | `4 dtheta^2` |

The worked source is a witness, not a fitted physical constant.

## 10. What single-source closure can and cannot mean here

There are two different demands:

1. **Dependency closure:** construct the output from native operations and already earned theorems, with no independent angle/noise/stopping/metric parameter slipped in. Sections 2–8 establish this for the declared source class and chosen record/contrast constructions.
2. **Source selection:** derive the actual transition, preparation, observation policy and contrast selection uniquely from the primitive cut rules, without specifying a source instance. That is not established by the existing sources or by R14.

**Native non-uniqueness witness.** In the same oriented field, both `z1=3+4 iota` and `z2=5+12 iota` satisfy every scalar law used above, but give `q1=16/25` and `q2=144/169`. Even for fixed `z`, source ledgers `{+1,-1}` and `{+1}` give hidden variances one and zero. Both the coherent path and the fully tagged record tree are defined in the same native algebra and give different two-step visible contents. Therefore the cited algebra cannot imply a unique numeric event/noise law independently of transition, preparation and protocol. ∎

This witness identifies a missing **selection law**, rather than an excuse to import one. A future source-selection theorem must exclude or relate these alternatives by a rule derived from an explicitly identified native primitive. Declaring the alternatives a single source tuple, calling symmetry “automatic”, or renaming normalized norm content “physical probability” does not supply that theorem.

The starting oriented cut construction and its metatheory remain acknowledged. “Everything derived from no premises” is not claimed. The research target is one transparent native dependency chain with a proved selection mechanism, not an empty hypothesis list.

## 11. Verification and preservation

[Implementation](../04-operator-evolution/native_source_foundation.cjs), [tests](../04-operator-evolution/test_native_source_foundation.cjs), [source pins](../04-operator-evolution/R14_SOURCE_PINS.json) and [finite certificate](../04-operator-evolution/R14_VERIFICATION.json) are separate from this written proof. They use the pinned **canonical RKF engine unchanged**, through an adapter; no engine code is copied here.

```bash
node 04-operator-evolution/verify_r14.cjs \
  --rkf-root /path/to/Recognition-Kernel-Framework \
  --output /tmp/R14_VERIFICATION.json
```

The certificate exercises exact raw and unit-source histories, native observer completion, arbitrary finite preparation ledgers, entire finite record trees and prefix consistency, local geometry, native logarithm/phase enclosures, and mathematical mutation controls. A PASS is finite executable evidence; it is neither a universal proof, a Lean formalization, nor a physical experiment. R1–R13 files and source pins are preserved.

[F00]: https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/certificates/foundation/F00_F00E_EULER_V0_1/README.md
[LOG]: https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00G_NATIVE_LOGARITHM_AND_POWERS.md
[PERIOD]: https://github.com/Parveen117/RH-Framework/blob/c588ded973a395b5fce37c82f670616160979a2e/theorems/foundation/F00J_NATIVE_PI_AND_FUNDAMENTAL_CIRCULAR_PERIOD.md
[POLAR]: https://github.com/Parveen117/RH-Framework/blob/c588ded973a395b5fce37c82f670616160979a2e/theorems/foundation/F00K_NATIVE_POLAR_DECOMPOSITION_AND_COMPLEX_LOG.md
[RH01]: https://github.com/Parveen117/RH-Framework/blob/fd104f46cc8c99c05f0e05b10f18c9332e8219de/theorems/T01_NATIVE_SUMMABILITY_AND_RECOGNITION_COMPLETION.md
[COUNT]: https://github.com/Parveen117/RH-Framework/blob/fd104f46cc8c99c05f0e05b10f18c9332e8219de/theorems/T01_E_NATIVE_REFINEMENT_COUNT_MEASURE.md
[NATIVE]: https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md
[F06]: https://github.com/Parveen117/RH-Framework/blob/40c4d8730e6015f018be2ec0783d6922c182074e/theorems/foundation/F06_NATIVE_ADJOINTABLE_OPERATOR_ALGEBRA.md
