# R11: spectral blindness, native curvature and the observer that closes the path

Research owner: **Monty Dabas**. Development and certificate date: **1 October 2026**.

R11 connects the two existing spectral papers in Publications to [R8](CURVATURE_OBSERVATION_R8.md), [R9](CURVATURE_BALANCE_R9.md) and [R10](NATIVE_CURVATURE_DESCENT_R10.md). A whole native curvature family has a flat exact Riemann quotient, constant characteristic data, and finite order loops with identity spectrum. A complete cross-sector marker bank recovers its hidden curvature with an exact minimum scalar count.

The same marker closes the hidden path in the spectral paper's cubic order channel and produces R10's visible excursion when admitted as an actual feedback transport. Reading a response and introducing a transport interaction are separately declared protocols.

## 1. Existing spectral foundation

**When Spectra Forget Order: Spectral Blindness, Observer-Marked Separation, and Stable Reconstruction of Ordered Matrix Products** (57 pages) proves two-factor characteristic blindness in Theorem 5.1, inserts markers before compression in Section 7, and gives model-relative reconstruction, stability and the endpoint/history boundary in Sections 15–19.

**Stable Recovery Beyond Spectral Blindness: Finite Sampling, Minimal Observer Lifts, and Completion-Stable Reconstruction** (32 pages) gives the minimum arbitrary linear scalar lift count in Theorem 7.2, the fixed-catalogue criterion in Theorem 4.4, and staged/target blindness in Sections 5–6.

These are existing results. R11 applies and certifies them on an explicit native curvature family. The degree-five and degree-six necklace calibrations are credited but not rerun here. No private MP source is copied or executed. The [pins](../04-operator-evolution/R11_SOURCE_PINS.json) bind both PDF bytes, native source bytes, public commits and the observer-facing reading scope.

## 2. A family that is spectrally flat and natively curved

Admit two visible tangent modes and h hidden modes, with finite h at least one. For any rational h-by-2 matrix L, define

\[
P_h=\begin{pmatrix}0_2&0\\0&I_h\end{pmatrix},
\qquad N_L=\begin{pmatrix}0_2&0\\L&0_h\end{pmatrix}.
\]

These are a hidden-sector projector and a visible-to-hidden recording map:

\[
P_h^2=P_h,\quad P_hN_L=N_L,\quad N_LP_h=0,\quad N_L^2=0.
\tag{1}
\]

**Theorem R11.1.** The admitted constant two-direction connection \(A_u=P_h,\ A_v=N_L\) has

\[
\boxed{F_{uv}=[P_h,N_L]=N_L.}
\tag{2}
\]

It is nonzero precisely when L is nonzero, but

\[
\boxed{\det(I-qF_{uv})=1\quad\text{for every }q.}
\tag{3}
\]

All its diagonal entries and all positive trace moments vanish.

**Proof.** Equation (1) gives (2). The characteristic matrix is block triangular with identity diagonal blocks. The strict lower shape has zero trace, and every higher power of N vanishes. This proves (3) and the trace claims.

The constant surjective observer \(C=(I_2\;0)\) satisfies \(CA_u=CA_v=0\), including all differentiated intertwinings. The complete R10 gate therefore certifies the flat two-mode Levi-Civita quotient:

\[
CF_{uv}=0=R_{uv}(g)C.
\tag{4}
\]

The family has **2h independent native curvature coordinates** beyond this flat exact classical quotient. Curvature follows from ordered composition; it is not inferred merely from an entry's position.

### Nonlinear spectral data and lawful rank statements

The characteristic map is nonlinear on general matrices. R11 does not assign it a linear kernel. On this declared linear family its entire output is exactly constant. Subtracting the L equal to zero reference gives the zero map on the 2h-dimensional coordinate space. The classical curvature readout is also zero there. This centered restriction is the lawful baseline for the source's linear minimum-lift theorem.

## 3. Nonidentity order loops with identity spectrum

Admit invertible rational moves \(X=I+aP_h,\ Y=I+bN_L\), where a differs from minus one. Their inverses are

\[
X^{-1}=I-\frac{a}{1+a}P_h,\qquad Y^{-1}=I-bN_L.
\]

**Theorem R11.2.**

\[
\boxed{T=XYX^{-1}Y^{-1}=I+abN_L.}
\tag{5}
\]

For nonzero ab and L, the full carrier does not return, although

\[
\boxed{\det(I-qT)=(1-q)^{2+h},\qquad CT=C.}
\tag{6}
\]

**Proof.** Multiply using (1). The return is block triangular with identity diagonal blocks. Its full difference from identity is ab times N. This also gives \(XY-YX=abN_L\), while the unequal XY and YX endpoints obey the source's two-factor spectral blindness.

R7's existing labelled moves are composed chronologically as inverse Y, inverse X, Y, X. Their final carrier is (5). Sheet winding is zero for this particular loop; independent sheet memory remains a separate type.

The finite protocol and the local connection share declared generators. A physical parameter-to-clock or global connection-integration law is not assumed.

## 4. The exact observer repair

Let the insertion marker \(M_{rc}=E_{c,\,2+r}\), for hidden row r and visible column c, return hidden input to the selected visible output. Then

\[
\boxed{\operatorname{Tr}(M_{rc}N_L)=L_{rc}.}
\tag{7}
\]

This uses the spectral paper's insertion convention trace(M D). Under its Hilbert–Schmidt convention the analysis vector is the insertion marker's adjoint; transposing the insertion itself would make this response blind.

**Theorem R11.3.**

\[
\boxed{\det(I-qM_{rc}F_{uv})=1-qL_{rc},}
\tag{8}
\]
\[
\boxed{\det(I-qM_{rc}T)=1-qabL_{rc}.}
\tag{9}
\]

**Proof.** Each insertion marker is rank one, so each product has rank at most one and determinant one minus q times its trace. Equation (7) gives (8); trace(M) equals zero and (5) gives (9).

At calibrated nonzero q, (8) recovers L. Equation (9) also requires nonzero ab. The protocol gain cannot be silently omitted.

### Exact minimum and dimensional scope

For the complete curvature target, the source target-relative theorem gives

\[
\boxed{r_{\min}=\operatorname{rank}(\Pi|_{\ker A})=2h.}
\tag{10}
\]

The target is identity on L and the centered baseline observation is zero. The 2h markers in (7) attain the lower bound with coordinate response matrix \(I_{2h}\). Omitting one leaves the corresponding matrix unit of L as an explicit nonzero blind witness.

For two hidden modes, **four independent linear scalar channels are necessary and sufficient**. This differs from R8's two-channel K-odd target because the target and carrier differ. It is not a universal instrument count: one device could return several independent scalar coefficients. Nor does it reconstruct arbitrary curvature matrices of size four.

No fixed finite scalar count can recover the entire family for every h; eventually 2h exceeds that count.

### A forbidden-feedback catalogue stays blind

If every admissible marker has no hidden-to-visible block,

\[
M=\begin{pmatrix}D&0\\B&H\end{pmatrix},
\qquad
MN_L=\begin{pmatrix}0&0\\HL&0\end{pmatrix}.
\]

All marked products remain nilpotent and their determinants remain one for every q. More parameter samples from this catalogue cannot recover any L coordinate. The missing ingredient is the reverse-sector leg.

Smaller targets have different budgets. For h equal to two, the target \(L_{11}+L_{22}\) needs one synthesized sum response, but an entry-only catalogue needs two responses. This is the source's catalogue overhead applied to curvature.

The real symmetric markers \(H_{rc}=M_{rc}+M_{rc}^{T}\) have the same response and give (8) on this curvature family because the added strict lower map annihilates N. They provide a Hermitian option; positive quantum effects or hardware implementation are not assumed. Equation (9) uses the rank-one insertion bank with its stated baseline.

## 5. Cubic spectral excitation and the R10 excursion

The existing spectral Theorem 6.1 states

\[
[xyz](\mathcal L_{CBA}-\mathcal L_{CAB})
=\frac{q}{(1-q)^2}\operatorname{Tr}(A[B,C]).
\tag{11}
\]

Choose \(A=P_h,\ B=N_L,\ C=M_{rc}\). By trace cyclicity,

\[
\operatorname{Tr}(P_h[N_L,M_{rc}])
=\operatorname{Tr}(M_{rc}[P_h,N_L])=L_{rc}.
\]

Thus

\[
\boxed{[xyz](\mathcal L_{MNP}-\mathcal L_{MPN})
=\frac{q}{(1-q)^2}L_{rc}.}
\tag{12}
\]

The pair is spectrally blind; a third controlled reverse leg creates a connected spectral response. Nonfeedback markers leave that channel zero.

R11 independently evaluates this multilinear coefficient using the exact formal trace-log expansion. W, W squared and W cubed are retained; higher powers cannot contribute to xyz. Truncation by \(x^2=y^2=z^2=0\) preserves this coefficient. No numerical exponential or global Magnus convergence is needed; q differs from one for the normalized expression.

Now write a strictly upper marker with block U. If admitted as an actual transport generator, its bracket with N is

\[
[M,N_L]=\operatorname{diag}(UL,-LU).
\]

The R10 decomposition for this pair and a flat visible metric gives

\[
R=0,\quad\mathcal D=0,\quad\mathcal E=UL,
\qquad
\boxed{\operatorname{Tr}\mathcal E=\operatorname{Tr}(MN_L).}
\tag{13}
\]

The spectral response and the curvature excursion use the same closed composition: the hidden record travels back through the observer's reverse leg. Operating that leg can alter the transport model; merely reading its declared response does not assert such an intervention.

## 6. Conditional balance, orientation and response geometry

The cut derived from P is

\[
J=I-2P_h,\qquad J^2=I,\qquad JN_LJ=-N_L.
\]

Within R9's declared branch protocol,

\[
\Phi_\theta(N_L)=(1-2\theta)N_L.
\tag{14}
\]

At theta equal to one half, the averaged readout vanishes throughout this family while the raw connection can remain curved. Keeping the realized branch permits exact inverse recovery by J. Dropping it at this endpoint destroys all 2h coordinates in that readout. This is conditional balance, not a universal statement about observation.

An active cut reverses fixed-marker responses. Simultaneously transporting carrier and marker frame preserves their trace and determinant responses. Orientation reversal and passive covariance remain distinct.

For \(\lambda=1-2\theta\), the complete response bank has

\[
B_\theta=\lambda I_{2h},\qquad
G_{\mathrm{response}}=B_\theta^TB_\theta=\lambda^2I_{2h}.
\tag{15}
\]

Nonzero lambda permits exact recovery, with growing noise amplification near zero. At zero the response rank vanishes. This Gram matrix belongs to declared response coordinates; it is not identified with a spacetime metric or quantum Fisher information.

The source blindness ledger avoids double counting. In the chain complete coordinate response, balanced erasure, classical readout, spectral summary, target blindness is \(0,2h,2h,2h\), with increments \(2h,0,0\).

## 7. Stable and reviewable curvature decisions

With known family gain kappa and supplied deterministic noise bound,

\[
d_{rc}=1-q\kappa L_{rc}+\eta_{rc},
\qquad|\eta_{rc}|\le\varepsilon,\qquad q\kappa\ne0.
\]

Decode \(\widehat L_{rc}=(1-d_{rc})/(q\kappa)\). Each entry error is bounded by epsilon divided by the absolute gain, and

\[
\boxed{\|\widehat L-L\|_F^2
\le\frac{2h\varepsilon^2}{|q\kappa|^2}.}
\tag{16}
\]

Extremal error in every channel attains this bound.

An entry interval excluding zero certifies nonzero curvature within the family. Exact zero responses with zero error certify zero curvature within the family. Zero responses with a nonzero error interval remain unresolved because small curvature may fit that bound.

The model restriction matters: a hidden RK block outside this family is invisible to this bank. The operator-family extractor rejects it. Four entry channels are not a universal full-matrix flatness certificate.

## 8. Endpoint completeness preserves the history boundary

Since N squared is zero, \((I+N)(I-N)=I=I\,I\). Every endpoint-only marked determinant agrees for these different histories. More endpoint markers cannot recover a history already erased at the endpoint map.

R7's independent sheet register can differ even at identity carrier endpoint. It requires retained sheet data or its explicit ledger. R11 does not infer it from spectral data.

## 9. Certification, lineage and reproduction

| Evidence | Coverage |
| --- | --- |
| Written proofs | Finite h curvature, flat tangent quotient, exact finite return, cross-sector repair, minimum scalar count and stable decisions |
| Native replay | Eight equalities and one distinctness witness, including cut grading and the rational four-move return |
| Exact rational tests | 31 tests; 81 lower-coupling cases; 405 ordinary and 1,620 marked curvature determinants |
| Boundary tests | All subbanks, forbidden-feedback catalogue, model refusal, orientation, noise and endpoint-history cases |
| Negative controls | Six mathematical mutations and two native alterations rejected |
| Preserved evidence | R1–R10 and master-review records and their recorded hashes |
| Further validation | No new Lean proof, private MP suite execution, external review or physical experiment claimed |

The new certified bridge is the shared explicit family across **native curvature, exact Riemann quotient, finite return, spectral blindness, minimal marked repair and closed-leg excitation**. Nilpotent blocks, rank-nullity, trace pairing, finite frame stability and the source spectral theorems keep their established lineage.

The canonical engine remains in RKF's operator_foundation. No engine implementation or spectral PDF body is copied into this repository.

Run with Python 3.11 or 3.12, Node and separate source checkouts at the pinned commits:

    python3.12 -B 04-operator-evolution/verify_r11.py \
      --publications-root ../Publications \
      --rkf-root ../Recognition-Kernel-Framework

See the [implementation](../04-operator-evolution/spectral_curvature_observer.py), [tests](../04-operator-evolution/test_spectral_curvature_observer.py), [source pins](../04-operator-evolution/R11_SOURCE_PINS.json), [native packet](../04-operator-evolution/R11_NATIVE_CERTIFICATE.json) and [verification record](../04-operator-evolution/R11_VERIFICATION.json).

The next physical selection question is which cross-sector catalogue a specified native response sector actually permits, and whether its couplings are measurements, transport interventions, or both. That choice determines which hidden curvature can be certified and how feedback modifies visible geometry.

## Source links and reading scope

- [When Spectra Forget Order](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/SPECTRAL%201.pdf): characteristic blindness, marked traces, cubic response, model-relative banks, stability and history scope.
- [Stable Recovery Beyond Spectral Blindness](https://github.com/Parveen117/Publications/blob/e1dc4e3773f063f14e56cf30222c8e8a504519cd/SPECTRAL_2.pdf): sampling, minimum target lifts, fixed-catalogue criteria and staged blindness.
- [Native algebra](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/01_NATIVE_ALGEBRA.md): lawful composition, observer kernels and the unchanged canonical engine.
- R8–R10 provide local curvature, conditional grading and complete tangent descent. Their proof and code bytes are unchanged.

