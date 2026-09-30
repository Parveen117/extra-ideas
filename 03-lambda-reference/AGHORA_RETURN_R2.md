# R2: Aghora return law and the effect of a projection

Research program: Monty Dabas / Extra Ideas. Development date: 30 September 2026.

**Result.** The source relations `A^2=I` and `AGA=-G` support a specified reversible return `R_t=A exp(tG)` with `R_t^2=I`. Its eigenvalues are confined to `+1` and `-1`, independently of the flow parameter. Inserting a bond projection can produce a variable effective response. The exact defect is

$$B-(BR_tB)^2=BR_t(I-B)R_tB.$$

This progresses from R1's arbitrary scalar cycle gains to a constrained operator return built from relations actually stated in `4ways.tex`. The flow specialization and return protocol are explicit modelling choices. The result does not select a unique sign, a universal physical constant, or the manuscript's still-unspecified six-stage closing operation.

## 1. Source relations and choices

Write `A` for the manuscript's Aghora operator, `G` for its distortion operator, and `B` for its bond projection.

| Item | Status in R2 |
| --- | --- |
| `A^2=I`, `AGA=-G` | Given in the Aghora Operator definition of [4ways.tex](../4ways.tex) |
| `B^2=B` | Algebraic content of the source's Bond Operator projection |
| A finite-dimensional real or complex vector space `V` | A stated setting that avoids unbounded-operator domain questions |
| A fixed generator `G`, with `U_t=exp(tG)` | A specialization of the source's general exponential-flow idea; the source does not identify its `D_v` with `G` |
| `R_t=A U_t` as the return protocol | A construction studied here; it is not established as the unique `YA -> OM` operation |
| Orthogonality, adjoints, and norms | Additional structure used only in the explicitly marked metric subsection |

The flow parameter `t` need not represent physical time. Keep `A`, `G`, and `t` the same when applying the same return twice. None of the algebraic results requires a Hilbert-space primitive.

## 2. The constructed return closes after two applications

**Theorem 1 (reversal and double return).** Under the assumptions above,

$$A U_t A=U_{-t}=U_t^{-1},\qquad R_t^2=I.$$

Moreover,

$$R_t=U_{-t/2}\,A\,U_{t/2},$$

so `R_t` is similar to `A`. For unequal parameters,

$$R_{t_2}R_{t_1}=U_{t_1-t_2}.$$

**Proof.** From `AGA=-G`, conjugation of each power gives `AG^kA=(-1)^kG^k`. The finite-matrix exponential series converges, so conjugating the series gives `A exp(tG) A=exp(-tG)`. Therefore `R_t^2=A U_t A U_t=U_{-t}U_t=I`. Applying the same reversal relation with the half-parameter proves the displayed similarity. Finally, `A U_{t_2}A U_{t_1}=U_{-t_2}U_{t_1}=U_{t_1-t_2}` because these exponentials have the same fixed generator. QED.

The unequal-parameter identity prevents overinterpreting the result: two arbitrary returns do not necessarily cancel. Nor does `AGA=-G` license an unrelated generator `D_v`. For example, choosing `U=2I` gives `(AU)^2=4I`.

## 3. Parameter-independent signs and a necessary two-mode structure

**Corollary 2 (return signs).** Every eigenvalue `r` of `R_t` satisfies

$$r^2=1.$$

The projectors

$$P_+(t)=\tfrac12(I+R_t),\qquad P_-(t)=\tfrac12(I-R_t)$$

satisfy `P_++P_-=I`, `P_+P_-=0`, `P_+^2=P_+`, `P_-^2=P_-`, and `R_t P_+=P_+`, `R_t P_-=-P_-`. The multiplicities of the signs equal those for `A` and do not depend on `t`.

**Proof.** Use `R_t^2=I`. The polynomial `(x-1)(x+1)` has distinct roots over the real or complex scalars, yielding the direct sum of its two eigenspaces and the stated projector formulas. Similarity to `A` fixes their dimensions. QED.

These are algebraic projectors; they need not be orthogonal before a compatible inner product is supplied.

**Proposition 3 (one scalar mode is insufficient).** If `G` is nonzero, `A` and therefore `R_t` have both sign sectors, so `dim V >= 2`. Relative to the `A` eigenspaces, `G` has the form

$$G=\begin{pmatrix}0&D\\E&0\end{pmatrix}.$$

**Proof.** The relation `AG=-GA` sends the positive eigenspace into the negative one and vice versa, proving the block form. If only one eigenspace existed, `A` would be `I` or `-I`, making `AGA=G=-G` and hence `G=0`. QED.

This is why the R2 return is an operator acting on modes, extending R1's scalar transport model. The exact sign set is constrained by the source involution, while the generator scale remains free: `G -> cG` preserves both source relations.

## 4. What would select a single sign?

**Proposition 4 (bond selection).** Let `B` be a nonzero idempotent and `epsilon` be `+1` or `-1`. Every vector in the retained space `ran B` has return multiplier `epsilon` if and only if

$$R_tB=\epsilon B.$$

Equivalently, `ran B` lies in the `epsilon` eigenspace. In that case `BR_tB=epsilon B`.

**Proof.** Every retained vector is `Bx`, so the statement about its full return is precisely `R_tBx=epsilon Bx` for every `x`. Premultiplying by `B` gives the compression identity. QED.

The choices `B=P_+(t)` and `B=P_-(t)` realize both signs. Commutation `[B,R_t]=0` alone does not choose one: `B=I` retains both sectors whenever `G` is nonzero. In higher rank, a commuting projection may also retain parts of both sectors.

Thus `-1` is an available and exact return multiplier, but the source has not supplied a rule selecting the negative bond sector. Merely naming that sector would add a selection assumption. A compression equation `BR_tB=epsilon B` also need not imply the stronger full-return equation in an arbitrary algebraic, nonorthogonal setting.

## 5. Exact defect caused by an intervening projection

Fix `t` and abbreviate `R=R_t`. Set `C=BRB`. The operator `C` acts within `ran B`, where `B` is the identity.

**Theorem 5 (projected double-return defect).** For any idempotent `B` and involution `R`,

$$\boxed{B-C^2=BR(I-B)RB.}$$

**Proof.** Write `I=B+(I-B)` between the two factors of `R` in `BR^2B`. Since `R^2=I` and `B^2=B`,

$$B=BRBRB+BR(I-B)RB=C^2+BR(I-B)RB.$$

Rearrange. QED.

The right-hand side follows the component that leaves the retained subspace under one return and re-enters under the next. Full evolution includes this contribution. Inserting `B` after each return removes the intermediate component. In the purely algebraic setting the defect has no asserted sign, norm interpretation, or entropy interpretation.

## 6. Optional metric interpretation

Now additionally choose an inner product and assume

$$A^*=A,\qquad G^*=-G,\qquad B^*=B.$$

Then `U_t` is unitary and `R_t` is a self-adjoint involution: `R_t^*=U_{-t}A=A U_t=R_t`. With `Q=I-B` and `L=QR_tB`, Theorem 5 becomes

$$B-(BR_tB)^2=L^*L\ \succeq\ 0.$$

Consequently the compression on `ran B` is a self-adjoint contraction with spectrum in `[-1,1]`. For `x` in `ran B`,

$$\|x\|^2-\|BR_tx\|^2=\|(I-B)R_tx\|^2.$$

For a normalized retained eigenvector with compressed eigenvalue `mu`, the retained squared norm is `mu^2` and the omitted squared norm is `1-mu^2`. This is a norm balance, not by itself a probability postulate or a definition of thermodynamic entropy.

**Proof.** Since `R_t` and `B` are self-adjoint, `L^*L=BR_tQ^2R_tB=BR_tQR_tB`. The defect identity gives positivity. Orthogonal decomposition of the unit-norm-preserving output `R_tx` gives the norm balance. QED.

In this metric setting, zero defect is equivalent to `L=0`, and self-adjointness then makes the retained subspace reducing, so `[B,R_t]=0`. The orthogonality of `B` alone is insufficient when `R_t` is not self-adjoint.

## 7. A fully exact two-mode example

Choose

$$A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
G=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
B=\begin{pmatrix}1&0\\0&0\end{pmatrix}.$$

These matrices satisfy the source relations and the optional Euclidean adjoint conditions. For exact rational evaluation, use the reversible Cayley family

$$U(s)=(I+sG)(I-sG)^{-1}.$$

It obeys `A U(s) A=U(s)^{-1}`, so `R(s)=A U(s)` is an involution. This family is not being identified with `exp(sG)` at the same parameter. Here it is the same rotation subgroup with a different parameterization.

Direct multiplication gives

$$R(s)=\frac{1}{1+s^2}
\begin{pmatrix}1-s^2&2s\\2s&s^2-1\end{pmatrix}.$$

The full return eigenvalues stay `+1,-1`, while the coefficient seen after applying `B` is

$$\mu_B(s)=\frac{1-s^2}{1+s^2}.$$

For the illustrative choice `s=1/3`,

$$R=\begin{pmatrix}4/5&3/5\\3/5&-4/5\end{pmatrix},\qquad
BRB=\tfrac45B,\qquad
B-(BRB)^2=\tfrac9{25}B.$$

The input `(1,0)` returns as `(4/5,3/5)` before the projection. The projection retains `(4/5,0)`. The omitted squared norm is `9/25`; two full returns still give the original input. Applying the projection after each return instead gives `(16/25,0)` after two passes.

Thus this model explicitly supports fixed full-return signs alongside a variable projected coefficient. The chosen `1/3` and resulting `4/5` are worked-example inputs and consequences, not derived physical constants.

## 8. The algebraic projected coefficient can remain arbitrary

**Proposition 6.** For `A=diag(1,-1)` and any real or complex scalar `mu`, define

$$B_\mu=\frac12
\begin{pmatrix}1+\mu&1-\mu\\1+\mu&1-\mu\end{pmatrix}.$$

Then `B_mu^2=B_mu`, `rank B_mu=1`, and

$$B_\mu A B_\mu=\mu B_\mu.$$

**Proof.** Write `u=(1,1)^T` and `v^T=((1+mu)/2,(1-mu)/2)`. Then `v^Tu=1`, so `B_mu=uv^T` is a rank-one idempotent. Since `v^TAu=mu`, the compression identity follows. QED.

These projections are generally oblique, so there is no conflict with the metric bound in Section 6. For a nonzero-parameter exponential return, conjugate `A` and `B_mu` together by `U_{-t/2}`; the same compressed scalar follows for `R_t`. Thus arbitrary algebraic compression values are compatible with the source reversal relations.

Even an orthogonal `B` does not suffice without a self-adjoint return. Take `G=[[0,1],[0,0]]`, `t=4`, and `B=(1/2)[[1,1],[1,1]]`. Then `R_t=[[1,4],[0,-1]]`, `R_t^2=I`, but `BR_tB=2B` and the defect is `-3B`.

## 9. Reference changes and the four routes

For any invertible coordinate change `S`, transform the entire model:

$$A'=SAS^{-1},\quad G'=SGS^{-1},\quad B'=SBS^{-1}.$$

Then `R_t'=SR_tS^{-1}`, `C'=SCS^{-1}`, and the defect transforms in the same way. The spectra of the full return and its retained-space compression are invariant. Changing only `R_t` while freezing the numerical matrix of `B` is generally a change in the model, not a consistent change of reference. Norm claims also require carrying the inner product along, or restricting to unitary coordinate changes.

| Way | R2 realization |
| --- | --- |
| 1: intrinsic modes | The direct sum of positive and negative Aghora sectors; a nonzero odd generator connects them |
| 2: response | The action on a chosen retained space, including `mu_B` and the omitted-channel term |
| 3: invariant values | The parameter-independent full-return sign set, with a separate bond-selection obligation |
| 4: operators | The reversible flow, `R_t=A exp(tG)`, and the exact projection-defect identity |

## 10. What has been selected, and what has not

The source involution now constrains this return family by the polynomial `r^2-1=0`. This removes R1's unrestricted multiplier freedom **within this specific family**. It does not constrain every possible closed process in the framework.

For a nonzero generator, both signs occur in the full space. A native rule for the bond or seam subspace would be needed to select one sign for retained states. A non-invariant bond can instead produce continuously varying effective responses, determined by the actual operator and projection data.

The scale of `G`, the flow parameter, the choice of return protocol, and the bond geometry remain undetermined by the two source equations. For example, replacing the rotation generator in Section 7 by `2G` at the same Cayley parameter `1/3` preserves the full-return signs but changes the compressed coefficient from `4/5` to `5/13`.

The next concrete obligation is to derive the bond projection and its relation to this return from a typed native seam rule. Either prove `R_tB=epsilon B` for a distinguished sector, or derive a particular noncommuting `B` and compute its projected response. The manuscript's names and its current untyped unified composition do not yet provide that selection law.

## 11. Reproduction and mathematical provenance

Run from the repository root with Python 3.11 or 3.12:

```bash
python3.12 -B 04-operator-evolution/verify_r2.py
```

The [implementation](../04-operator-evolution/aghora_return.py) uses exact rational matrices. The [focused tests](../04-operator-evolution/test_aghora_return.py) include nilpotent exponentials, a rational reversible family, sign projectors, projection defects, simultaneous reference changes, and counterexamples to missing hypotheses. The [verification record](../04-operator-evolution/R2_VERIFICATION.json) records file hashes and the scope of the check. The proofs above establish the general finite-dimensional statements; the finite examples do not replace those proofs.

Reversibility by an involution and factorization into two involutions are established mathematical constructions. See Michael Baake, *A brief guide to reversing and extended symmetries of dynamical systems*, Section 2, [arXiv:1803.06263](https://arxiv.org/html/1803.06263v2). R2 derives the relevant formulas directly from the manuscript's stated algebra and connects them to its bond projection. It does not claim a new discovery of involutions, operator compression, or fundamental physical constants.
