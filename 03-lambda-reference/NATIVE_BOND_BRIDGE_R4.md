# R4: The native source of R3's bond, and a closure test that normalization cannot replace

Monty Dabas · 30 September 2026.

## What advances, and what is inherited

R3 chose a seam and imposed its Aghora image as the removed space. The existing RKF paired-depth return supplies a particular bond without imposing that complement separately. At its exact-cut response, the existing grading exchanges the retained and removed spaces and its commutator with the bond recovers the existing EMK quarter-turn.

This closes a dependency for **that source family and aperture**. It is not a universal deduction from the four-way manuscript alone. The source's exact-cut synthesis, its quarter-turn, its scalar return recurrence, and its remaining free coupling are earlier results. They are consumed here, not announced again.

The additional operational result links the **raw** bond and complex-structure defects, transports the source's finite-aperture enclosures to them, and exposes two false closure tests: normalizing a commutator, and treating an accidentally closed finite aperture as a closed completed response.

See [CROSS_REPO_LINEAGE.md](../CROSS_REPO_LINEAGE.md) for the duplication assessment and [R4_SOURCE_PINS.json](../04-operator-evolution/R4_SOURCE_PINS.json) for exact source revisions and hashes.

## 1. Native input and the notation bridge

Use the canonical RKF paired-depth family at commit `3cc5a33b05c16d59c90994ddda69dedc0d392424`, specifically [R2, RD1–RD2, RI1–RI2 and SY1](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/research/recognition_return/r2/THEOREM.md). Its internal algebra has

$$
R^2=-I,\quad K^2=I,\quad KR=-RK,\quad L=KR,\quad L^2=I,\quad L^\dagger=L.
$$

`R` here is the source's EMK quarter-turn, not R2's return operator. The existing `K` grading is documented in [T55](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/55_generalized_euler_emk_dock_theorem.md). It is different from the radial/turn grading of the scalar dagger.

The positive paired-depth source has inverse corners

$$
F=I+xL,\qquad f_a(x)=\frac{a}{1+ax},\quad a>0.
$$

For the positive two-cell period `a,b,a,b,...`, the completed response is the unique positive solution of

$$
bx^2+x-a=0.
$$

Source SY1 already proves that its normalized response `F/2` is idempotent exactly when `x=1`, equivalently `a=b+1`. Its convergence is by finite-aperture inverse completion. It is not an assertion that a Neumann series or a full infinite-operator inverse exists.

We instantiate the four-way variables as follows:

| Four-way object | Source object | Status |
| --- | --- | --- |
| Aghora involution `A` | Existing grading `K` | Explicit bridge between the models |
| Bond `B` | Normalized source response `F/2` at an exact cut | Source SY1's already defined normalization |
| Retained seam | `W=range(B)` | Supplied by the response, rather than an independently chosen line |
| Complex-structure operator `[A,B]` | Source quarter-turn `R` at exact closure | Consequence proved below |
| Flow generator `G` | `omega R` in the real two-mode compatible class | Shape inherited from R3; rate remains unspecified |

The choice of native source family and the identification `A=K` are explicit model connections. Nothing here proves that every possible four-way Aghora or bond must be this source realization.

## 2. The source bond forces R3's complementary exchange

**Bridge theorem.** On the source's exact-cut family, put

$$
B=\frac{I+L}{2}=\frac{I+KR}{2},\qquad A=K.
$$

Then

$$
\boxed{B^2=B,\qquad ABA=I-B,\qquad [A,B]=R.}
$$

**Proof.** The first equation follows from `L^2=I`. Since `KLK=-L`, conjugation gives the second. Finally `KL=R` and `LK=-R`, so

$$
[K,B]=\tfrac12(KL-LK)=R.
$$

In any faithful nonzero finite realization, an idempotent gives `V=range(B) direct-sum kernel(B)`. The exchange identity sends `range(B)` onto `kernel(B)`, so

$$
\ker B=A\operatorname{range}B.
$$

Thus the complement that R3 imposed is now a consequence of the **selected native exact-cut response**. Within scalar multiples of the nonzero `F`, the only nonzero idempotent is `F/2`: because `F^2=2F`, the equation `(cF)^2=cF` forces `c=1/2`.

For every positive `b`, the period `(b+1,b)` gives this same bond. Positivity of the source response chooses `x=+1` rather than `-1` when idempotence is required. This selects the retained `L` sector inside the declared source family. It does not select a sign for every subsequent dynamical return.

## 3. Exact connection to R3's coordinates and metric

In the established real two-coordinate calibration, use

$$
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
K=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
L=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
B=\begin{pmatrix}1&0\\0&0\end{pmatrix}.
$$

These are existing EMK matrices, recorded in [native algebra N01](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md). For finite real nonzero `Gamma`, define

$$
T_\Gamma=\begin{pmatrix}\Gamma&\Gamma\\1&-1\end{pmatrix}.
$$

Direct multiplication yields

$$
T_\Gamma K T_\Gamma^{-1}=A_{\rm R3},\quad
T_\Gamma B T_\Gamma^{-1}=B_\Gamma,\quad
T_\Gamma R T_\Gamma^{-1}=K_\Gamma,
$$

where the three objects on the right are exactly R3's involution, bond, and commutator. In particular, the source's retained coordinate line becomes `span((Gamma,1))`, while its removed coordinate line becomes `span((Gamma,-1))`.

On the real radial subcarrier, the supplied componentwise native norm has the usual sum-of-squares calibration. Transporting it gives

$$
H'=T_\Gamma^{-T}T_\Gamma^{-1}
=\tfrac12\operatorname{diag}(\Gamma^{-2},1).
$$

R3's `H_Gamma` is `2H'`, the allowed overall metric rescaling. This identifies the source of its relative weights in this representation. It is not a new derivation of a spacetime metric or an absolute norm scale.

The matrix `R` realizes the familiar quarter-turn on this real pair. Do not identify it with the central scalar `iota I` throughout the full native matrix algebra: `R` anticommutes with `K`, whereas a central scalar commutes with every native-linear matrix. The previously proved scalar-to-real-pair representation is a representation on a specified module, not permission to erase that distinction.

All nonzero `Gamma` in this displayed family are related by these coordinate maps when the whole represented data are transported. A physical restriction on admissible reference changes would need an additional statement; the formulas alone do not make `Gamma` a selected coupling.

## 4. A raw response carries linked closure defects

For a general nonnegative source response, retain its actual amplitude:

$$
B_x=\frac{I+xL}{2},\qquad J_x=[K,B_x].
$$

Call `B_x` a **bond candidate** until idempotence is established. Algebra gives

$$
KB_xK=I-B_x,\quad J_x=xR,
$$

$$
E_x:=B_x^2-B_x=\frac{x^2-1}{4}I,
\qquad C_x:=J_x^2+I=(1-x^2)I.
$$

Hence the exact shared-defect relation is

$$
\boxed{C_x=-4E_x.}
$$

The source already gives the formula for `E_x`; the added check connects it to the **unrenormalized commutator** used in R3. For positive `x`, bond closure and raw quarter-turn closure hold together exactly at `x=1`. Complementary exchange alone holds at every `x`, so it does not establish idempotence.

For a completed two-cell response, define the source mismatch `delta=a-b-1`. Its existing quadratic equation gives the useful calibration

$$
\delta=(x-1)[b(x+1)+1],\qquad
\frac{x^2-1}{4}
=\frac{\delta(x+1)}{4[b(x+1)+1]}.
$$

The denominator is positive. A positive/negative cell mismatch therefore gives the same sign of bond defect, and zero mismatch is exactly closure. This is a corollary of source SY1, not another independently selected constant.

## 5. Why normalization can manufacture a passing test

For every `x>0`, closed or otherwise,

$$
\widehat J_x=J_x/x=R,\qquad \widehat J_x^2=-I.
$$

Likewise sharpening the response gives

$$
\widehat B_x=\tfrac12\left[I+\frac{2B_x-I}{x}\right]
=\frac{I+L}{2}.
$$

Both operations are algebraically legitimate. They remove precisely the response amplitude needed to test the original boundary response's closure. They must not be reported as proof that the original `B_x` was idempotent. At `x=0`, the commutator is zero and this division is undefined.

The existing native period `(30,20)` provides a sharp exact control. Source SY1/SY3 gives `x=6/5`, since `20(6/5)^2+6/5=30`. Then

$$
B_x^2-B_x=\frac{11}{100}I,\qquad
J_x^2+I=-\frac{11}{25}I,
$$

yet `widehat J_x^2=-I` and `widehat B_x^2=widehat B_x`. This control uses an actual native return; its failure is not an arbitrary matrix invented outside the source family.

## 6. Finite aperture, completed response, and certified error

The unchanged RKF solver encloses every nonnegative-tail response for a prefix in `[ell,u]`, with `0<=ell<=u`. Monotonicity of squaring on that interval gives

$$
\frac{\ell^2-1}{4}\le e_x\le\frac{u^2-1}{4},\qquad E_x=e_xI,
$$

$$
1-u^2\le c_x\le1-\ell^2,\qquad C_x=c_xI.
$$

The distance to the exact source bond also obeys

$$
\boxed{M(B_x-B_1)=\frac{|x-1|}{2}
\le\frac{\max\{|\ell-1|,|u-1|\}}{2}.}
$$

Here `M` is the native normal-form coefficient mass with `M(L)=1`; the same number is the induced operator norm in the stated real calibration. It is not a sum of all matrix entries or an unspecified physical error measure.

If the interval excludes one, it rules out exact closure for every admitted tail response. An interval containing one generally proves only a bound. Exact completed closure for `(b+1,b)` follows from the source's infinite-profile theorem, not from rounding a finite output.

There is an additional false-positive control: the uniform profile `a_j=1`. Its first zero-tail finite response is exactly `x=1`, hence its finite bond is idempotent. At the three-cell prefix, the native solver encloses all nonnegative continuations in

$$
[\ell,u]=[1/2,2/3],
$$

excluding one. The bounded positive infinite profile has a unique response by source RI2; it satisfies `x^2+x-1=0`, so it is not an exact cut. Thus an exactly closed finite aperture need not certify the completed source.

Source nesting makes these transported budgets shrink as appropriate prefixes are refined. No finite prefix alone proves convergence of an unknown infinite tail. The known positive periodic profile is what licenses the completed-response conclusion here.

## 7. The selection achieved and the freedoms that remain

R4 obtains R3's specific complementary bond from an existing native response under an explicit source bridge. It does not establish that the paired source family, its aperture, or exact-cut target is uniquely selected by all native laws.

The complete period `(b+1,b)` still contains a free `b>0`. The earlier [Publications AS-1](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-alpha-selection/THEOREM.md) already shows that these identical bonds have distinct probe slopes `1/(1+2b)`. R4 does not count this non-selection result again.

Likewise the dynamical return `T_t=K exp(t omega R)` has

$$
BT_tB=\sin(\omega t)B.
$$

Both signs remain available at different phases, even though the source fixed the retained `L` sector. The rate `omega` is not identified with `b` or its probe slope without another dynamical theorem. No universal `lambda_*`, numerical action scale, or electromagnetic coupling has been selected.

The next substantive target is therefore a source law selecting a profile, aperture, or calibrated dynamics—not another derivation of the existing quarter-turn or compression identity.

## Reproduction

Use Python 3.11 or 3.12 and Node.js with its built-in modules. Provide a checkout of the pinned RKF revision; no upstream implementation is copied into this repository.

```bash
git clone https://github.com/Parveen117/Recognition-Kernel-Framework ../rkf-r4
git -C ../rkf-r4 checkout --detach 3cc5a33b05c16d59c90994ddda69dedc0d392424
python3.12 -B 04-operator-evolution/verify_r4.py --rkf-root ../rkf-r4
```

The verifier checks four upstream runtime hashes before executing the unchanged native solver and Pāṇinian proof replayer. It runs 12 focused integration tests, including 16 native polynomial proof replays, the actual paired Schur witness, R3 coordinate/metric correspondence, finite error bounds, and the false-positive controls. The tests exercise exact finite instances; the all-parameter conclusions follow from the written arguments and pinned source theorems.
