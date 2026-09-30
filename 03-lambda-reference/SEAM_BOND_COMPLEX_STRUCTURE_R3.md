# R3: A typed seam, its bond, and a linear complex structure

Research owner: **Monty Dabas**. Development: 30 September 2026.

## Result and scope

An explicit rule now fixes the bond that R2 left open: retain a seam subspace and remove its Aghora image, provided they are complementary. This rule gives a unique projection `B` and the exact identity

$$
K_\Gamma=[A,B],\qquad K_\Gamma^2=-I.
$$

In the real two-mode realization below, requiring Aghora to preserve an inner product and the bond to be orthogonal determines that inner product up to positive scale. An Aghora-odd generator preserving it must then have the form `G=omega K_Gamma`. Its rate `omega` remains free. The compressed exponential return has coefficient `sin(omega t)`; both return signs remain possible.

**The complement rule and metric requirements are added hypotheses.** The source's seam ratio alone does not imply them. This is a conditional construction with proofs and exact examples, not a selection of a universal `lambda_*`, a physical metric, or the complete six-stage closure.

## 1. Source data and the choices made here

The original [4ways.tex](../4ways.tex) supplies the seam ratio at lines 232–245, Aghora relations at lines 417–425, and a bond projection at lines 427–432. Its seam operator returns a scalar ratio; it is not yet a specified linear endomorphism that can be inserted in an operator product on states.

| Item | Status in R3 |
| --- | --- |
| `Gamma=lambda_p/lambda_v`; preservation under admissible transformations | Source statement; the admissible group still needs specification |
| `A^2=I`, `AGA=-G`, a bond projection | Source relations used in R2 |
| A real two-dimensional vector space, with one mode in each Aghora sector | Chosen minimal realization; broader state spaces are not ruled out |
| State coordinates `(x,y)` and the linear constraint `x-Gamma*y=0` | Added realization of the ratio idea; not a proved identification with the source's physical quantities |
| `range(B)=L`, `kernel(B)=A L`, with `V=L direct-sum A L` | Added bond-selection rule |
| A positive inner product making `A` an isometry and `B` self-adjoint | Additional metric-admissibility requirements |
| `G` skew-adjoint in this metric | Additional dynamical requirement |
| `U_t=exp(tG)` and `R_t=A U_t` | The particular return protocol from R2 |

Coordinates are real numbers after a choice of units. This construction neither supplies those units nor equates the coordinate ratio to a measured coupling without a further model. We write `mu` for the compressed response; it is not silently identified with `lambda_*`.

The new symbol `K_Gamma` denotes a complex-structure operator. It is distinct from the source's positive release generator named `K`.

## 2. Type the seam before composing it

Choose an Aghora eigenbasis, so

$$
A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
q_\Gamma:V\longrightarrow\mathbb R,\quad q_\Gamma(x,y)=x-\Gamma y.
$$

For a finite real nonzero `Gamma`, define

$$
L_\Gamma=\ker q_\Gamma=\operatorname{span}(u),\quad
u=(\Gamma,1)^T,\quad v=Au=(\Gamma,-1)^T.
$$

The constraint is defined even at the zero state, where a coordinate ratio would be undefined. The scalar map `q_Gamma` and a state map `B:V -> V` have different types. A bond associated with the seam should satisfy `q_Gamma B=0` and act as the identity on `L_Gamma`.

Those two conditions alone leave a choice: every rank-one map `B=u ell` with `ell(u)=1` is such a projection. In this basis `ell=(a,1-a Gamma)` for arbitrary real `a`. For example,

$$
\widetilde B=\begin{pmatrix}0&\Gamma\\0&1\end{pmatrix}
$$

retains the same seam but removes a coordinate axis. The source ratio alone therefore does not determine the bond.

## 3. Fix the removed space by an explicit rule

**Theorem 1.** Let `A^2=I` on a real vector space and suppose `V=L direct-sum A L`. There is a unique projection with range `L` and kernel `A L`. It obeys

$$
ABA=I-B.
$$

**Proof.** Every state has a unique decomposition `z=l+A m` with `l,m in L`. Set `Bz=l`. This proves existence, idempotence, range, kernel, and uniqueness. Also `Az=A l+m`, hence `ABA z=A m=(I-B)z`. No inner product is required. In finite dimension, this decomposition requires equally dimensional complementary spaces and therefore even total dimension.

For the two-mode seam, `det[u,v]=-2 Gamma`, so transversality is precisely `Gamma != 0` in this finite-slope chart. Solving `ell(u)=1`, `ell(v)=0` gives

$$
\boxed{B_\Gamma=\frac12
\begin{pmatrix}1&\Gamma\\\Gamma^{-1}&1\end{pmatrix}.}
$$

At `Gamma=0`, the seam is an Aghora eigenline and equals its own image. The infinite-slope eigenline fails for the same reason. These are excluded cases, not removable singularities of this selection rule.

The rule has a precise symmetry motivation: Aghora exchanges the retained and removed spaces. It remains a choice until the underlying model independently justifies this exchange property.

## 4. A complex structure follows without a metric

**Theorem 2.** Under Theorem 1, put `J=2B-I` and `K_Gamma=[A,B]`. Then

$$
J^2=I,\quad AJA=-J,\quad K_\Gamma=AJ,\quad
K_\Gamma^2=-I,\quad AK_\Gamma A=-K_\Gamma.
$$

**Proof.** Idempotence gives `J^2=I`, and the exchange law gives `AJA=2ABA-I=-J`. Also `BA=A(I-B)`, so `AB-BA=A(2B-I)=AJ`. Squaring yields `AJAJ=-J^2=-I`; conjugating by `A` gives the last relation.

In two dimensions,

$$
K_\Gamma=\begin{pmatrix}0&\Gamma\\-\Gamma^{-1}&0\end{pmatrix},
\qquad K_\Gamma u=v,\quad K_\Gamma v=-u.
$$

Thus `a+ib -> aI+bK_Gamma` represents the complex numbers as a real operator algebra:

$$
(aI+bK_\Gamma)(cI+dK_\Gamma)
=(ac-bd)I+(ad+bc)K_\Gamma.
$$

The representation is faithful because a real scalar multiple of `I` cannot square to `-I`. Every nonzero represented number has inverse `(aI-bK_Gamma)/(a^2+b^2)`.

This is the standard notion of a **linear complex structure**, here constructed from the specified pair `(A,B)`. It is not a claim to discover the complex numbers or to prove a new physical law. Reversing the commutator order or retaining the opposite seam changes `K_Gamma` to `-K_Gamma`; an absolute orientation has not been selected.

## 5. A relative metric, under stated compatibility conditions

**Theorem 3 (two real modes).** The symmetric positive matrices `H` satisfying

$$
A^T H A=H,\qquad B_\Gamma^T H=H B_\Gamma
$$

are exactly the positive multiples of

$$
\boxed{H_\Gamma=\begin{pmatrix}\Gamma^{-2}&0\\0&1\end{pmatrix}.}
$$

**Proof.** Write `H=[[h11,h12],[h12,h22]]`. The first equation forces `h12=0`; the second gives `h22=Gamma^2 h11`. Positivity leaves one positive multiplier. Conversely every such matrix satisfies both equations.

For `S=[u,v]`, direct substitution gives `S^T H_Gamma S=2I`. The seam and its Aghora image are therefore orthogonal and equally weighted. With the metric adjoint `M^dagger=H_Gamma^{-1} M^T H_Gamma`,

$$
A^\dagger=A,\quad B_\Gamma^\dagger=B_\Gamma,\quad
K_\Gamma^\dagger=-K_\Gamma,\quad
K_\Gamma^T H_\Gamma K_\Gamma=H_\Gamma.
$$

These statements derive the metric's shape **after imposing compatibility**, not from the source ratio alone. This is a positive inner product on a two-mode state space. It is not a spacetime metric and has no determined absolute length, time, or energy normalization. The uniqueness-up-to-scale statement is specific to two real modes; Theorem 1 does not by itself provide that uniqueness in higher dimension.

**Theorem 4.** A real two-by-two generator satisfying both `AGA=-G` and `G^dagger=-G` is precisely

$$
\boxed{G=\omega K_\Gamma,\qquad\omega\in\mathbb R.}
$$

**Proof.** The first equation gives `G=[[0,p],[q,0]]`. The second requires `q=-p/Gamma^2`. Set `omega=p/Gamma`. Conversely every `omega K_Gamma` satisfies both conditions. The zero generator is allowed.

The seam fixes the generator's relative matrix entries under these requirements, but not the rate `omega` or its sign.

## 6. The resulting return and the sign-selection condition

Set `theta=omega t`. Since `K_Gamma^2=-I`, its power series gives

$$
U_t=e^{tG}=\cos\theta\,I+\sin\theta\,K_\Gamma,
\qquad R_t=A U_t.
$$

In the seam basis `S=[u,v]`,

$$
S^{-1}AS=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
S^{-1}B_\Gamma S=\begin{pmatrix}1&0\\0&0\end{pmatrix},\quad
S^{-1}K_\Gamma S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
$$

$$
S^{-1}R_tS=\begin{pmatrix}\sin\theta&\cos\theta\\
\cos\theta&-\sin\theta\end{pmatrix}.
$$

Consequently,

$$
\boxed{B_\Gamma R_t B_\Gamma=\mu B_\Gamma,\quad\mu=\sin(\omega t),}
$$

$$
B_\Gamma-(B_\Gamma R_t B_\Gamma)^2=\cos^2(\omega t)B_\Gamma.
$$

For a seam input `u`, `R_t u=sin(theta)u+cos(theta)v`. Because the two directions are orthogonal with equal norm, the retained and removed squared-norm fractions are `sin^2(theta)` and `cos^2(theta)`. This is a norm decomposition, with no additional probability or entropy interpretation asserted.

The full return is an involution with both eigenvalues `+1,-1`. The chosen seam is a full return eigenspace exactly when `cos(theta)=0`; its eigenvalue is then `sin(theta)=+1` or `-1`. Thus requiring no leakage selects a discrete set of phases **after the protocol and rate are specified**, and still permits both signs. It does not calibrate `omega`, determine a preferred phase, or select a universal `lambda_*`.

## 7. Exact rational version and a reproducible example

The implementation uses the Cayley family

$$
U(s)=(I+sG)(I-sG)^{-1},\qquad z=\omega s.
$$

Since `G^2=-omega^2 I`, it is defined for every real `s`. It equals

$$
\frac{1-z^2}{1+z^2}I+\frac{2z}{1+z^2}K_\Gamma.
$$

This is not `exp(sG)`. As a rotation it corresponds to the angle `theta=2 arctan(z)` modulo a full turn. The rational return has

$$
\boxed{\mu(z)=\frac{2z}{1+z^2},\qquad
d(z)=\frac{(1-z^2)^2}{(1+z^2)^2}.}
$$

Here `BRB=mu B` and `B-(BRB)^2=d B`. The return in the seam basis is `[[mu,c],[c,-mu]]` with `c=(1-z^2)/(1+z^2)`. Zero leakage occurs at `z=+1` or `-1`, giving the corresponding sign on the retained seam. These numerical Cayley parameter values depend on the parametrization; the intrinsic criterion is preservation of the seam by the return.

Choose the illustrative inputs `Gamma=2`, `omega=3`, `s=1/9`. Then

$$
B=\begin{pmatrix}1/2&1\\1/4&1/2\end{pmatrix},\quad
K_\Gamma=\begin{pmatrix}0&2\\-1/2&0\end{pmatrix},\quad
H=\begin{pmatrix}1/4&0\\0&1\end{pmatrix},
$$

$$
G=\begin{pmatrix}0&6\\-3/2&0\end{pmatrix},\quad
R=\begin{pmatrix}4/5&6/5\\3/10&-4/5\end{pmatrix},\quad
S^{-1}RS=\begin{pmatrix}3/5&4/5\\4/5&-3/5\end{pmatrix}.
$$

Therefore `BRB=(3/5)B` and the defect is `(16/25)B`. For `u=(2,1)`, the retained squared-norm fraction is `9/25`, and the removed fraction is `16/25`. All inputs are chosen; none of these fractions is proposed as a universal constant.

## 8. Why the extra dynamical condition matters

Even after selecting this bond, the source's oddness condition alone does not make a zero scalar defect equivalent to no leakage. Set

$$
\Gamma=1,\quad B=\tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix},\quad
G=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad t=2.
$$

Here `G^2=0`, so the exact exponential is `I+tG`, and

$$
R=Ae^{2G}=\begin{pmatrix}1&2\\0&-1\end{pmatrix},\quad
R^2=I,\quad BRB=B,\quad B-(BRB)^2=0.
$$

Nevertheless `RB != B`. In the seam basis `R=[[1,0],[2,-1]]`: a seam state acquires a nonzero removed component which the later projection hides. This generator is Aghora-odd but is not skew-adjoint for `H=I`. At `t=4`, the same example gives `BRB=2B` and defect `-3B`, showing why positivity also needs the metric conditions.

For compatible dynamics instead, R2's identity becomes

$$
B-(BRB)^2=\mathcal L^\dagger\mathcal L,\qquad
\mathcal L=(I-B)RB.
$$

Positive definiteness then makes zero defect equivalent to zero leakage. This distinction is included in the tests.

## 9. Reference changes and what is actually invariant

For any invertible coordinate change `T`, transform the entire data together:

$$
A'=TAT^{-1},\quad L'=TL,\quad B'=TBT^{-1},\quad
K'=TK_\Gamma T^{-1},\quad H'=T^{-T}H_\Gamma T^{-1}.
$$

The subspace construction, complex-structure equation, return spectrum, and compressed coefficient are preserved. For the diagonal change `T=diag(a,b)` with nonzero real `a,b`, the canonical chart has

$$
\Gamma'=\frac ab\Gamma,\qquad
B_{\Gamma'}=TB_\Gamma T^{-1},\quad
K_{\Gamma'}=TK_\Gamma T^{-1},\quad
H_{\Gamma'}=b^2 T^{-T}H_\Gamma T^{-1}.
$$

The last multiplier is precisely the allowed overall metric scale. If the generator is transformed too, `omega=p/Gamma` is unchanged by these diagonal reference changes. Its value still has to be specified by a dynamical law and time convention.

This calculation does not prove the source's scalar `Gamma` invariant under every basis change. A coordinate ratio is covariant under this broader group. If admissible transformations must preserve its numerical value, that is a smaller group to be defined by the model. Under compatible dynamics the compressed response is independent of `Gamma`, depending instead on the phase or Cayley product `z`; varying the coordinate seam alone cannot select its value.

## 10. Connection to the four routes and the next open question

| Route | Concrete R3 object | Remaining limitation |
| --- | --- | --- |
| Way-1 | Real states `a u+b v` and complex operator action `aI+bK_Gamma` | This operator product is not the source's componentwise coordinate product |
| Way-2 | Typed constraint `q_Gamma`, response `mu`, and leakage | Physical interpretation and admissible transformations remain to be supplied |
| Way-3 | Reference-covariant operators and a return-sign criterion | No absolute scale, phase, rate, or universal `lambda_*` is selected |
| Way-4 | Unique conditional bond, complex structure, relative metric, compatible flow | The exchange rule and metric compatibility need independent motivation |

R3 removes arbitrary projection freedom **inside the stated seam-exchange model**. The next decisive question is whether the source's remaining operations force that exchange rule and its metric-compatible dynamics, or provide counterexamples. A competing complement or non-skew odd generator must not be excluded merely because it spoils a desired constant. The full source closure still requires correctly typed maps and an independently specified composition.

## Reproduction and mathematical provenance

Run from the repository root with Python 3.11 or 3.12:

```bash
python3.12 -B 04-operator-evolution/verify_r3.py
```

The [implementation](../04-operator-evolution/seam_bond.py) reuses R2's unchanged exact matrix arithmetic. The [18 focused tests](../04-operator-evolution/test_seam_bond.py) cover the construction, changes of reference, rational response, both signs, degeneracies, and the incompatible-generator counterexamples. [R3_VERIFICATION.json](../04-operator-evolution/R3_VERIFICATION.json) records the result and hashes of the proof, code, tests, and imported implementation. Finite examples check the implementation; the general identities follow from the proofs above.

Direct-sum projections and real linear complex structures are established linear algebra. For an explicit primary-source definition of an orthogonal complex structure, see Kirill Krasnov, *Geometry of Spin(10) Symmetry Breaking*, [arXiv:2209.05088v2](https://arxiv.org/pdf/2209.05088v2), Section 2.4, Definition 3. Only the terminology and defining equations are invoked here; none of that paper's particle-physics model is used. The contribution of this development is the conditional seam-to-bond construction and its audited consequences within this repository, not a novelty claim for those general mathematical notions.
