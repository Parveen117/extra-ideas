# Physics interface for the native compact gauge completion

Owner: Monty Dabas. Updated: 3 October 2026. Conditional application of Publications
[NCG-1–NCG-8][ncg], not a new primitive-only R47 certificate.

The useful new ingredient is a compatible internal gauge sector and a source
that actually transforms under it. Previously, SM's massive source carried
one native U(1) phase while the Yang–Mills benchmark supplied SU(2) separately.
NCG derives the permitted group inside the existing module/action contract.
The physical existence and selection of that contract are not inferred from
its algebraic consistency.

## Native construction and what it fixes

Take two independent EMK factors with R_i²=-1, K_i²=1 and K_iR_i=-R_iK_i.
Set C=R₁ and

\[
t_1=R_2/2,\quad t_2=R_1K_2/2,\quad
t_3=R_1R_2K_2/2.
\]

They are dagger-skew, commute with the massive multiplier C and satisfy
[t_a,t_b]=ε_abc t_c. The full norm- and multiplier-preserving internal
group is U(2), with this SU(2) sector and a remaining central phase. This uses
four copies of the existing real spin module, hence sixteen real coefficient
components. It is minimal only within that fixed massive-copy class.

For A_μ=A_μ^a t_a, the curvature is the native order defect:

\[
F^a_{\mu\nu}=\partial_\mu A^a_\nu-\partial_\nu A^a_\mu
                +\epsilon_{bca}A^b_\mu A^c_\nu.
\]

In particular A_x=t₁,A_y=t₂ gives F_xy=t₃ even though each scalar coefficient
is constant and its one-form locally exact. The mixed term is essential.
General noncommuting direction frames retain EMK-C1's additional -c_μν^ρ A_ρ.
Prior memory and observer compression retain NT's separate correction terms;
neither gauge covariance nor observation is a universal flattening operation.

## Consistent source, gauge equation and preserved constraint

On the admitted SM geometry, retain its massive commuting-field action and
replace the phase connection by A acting on the enlarged internal module.
With H=-C⊗Γ⁰ and γ^a=1₄⊗Γ^a, NCG-7 proves

\[
j_a^\mu=-\frac Z2\psi^T H\gamma^\mu(t_a\otimes1_4)\psi.
\]

Admit the positive internal pairing B(X,Y)=-4sc(XY) and the declared quadratic
gauge law S_A=-(4g_*²)⁻¹∫w B(F_μν,F^μν). NCG-6 then gives

\[
D_\mu(wF^{\mu\nu})=g_*^2 wj^\nu,
\qquad D_\nu(wj^\nu)=0\quad\text{on the matter equation}.
\]

The covariant divergence includes color precession; replacing it by an
ordinary divergence fails the independent exact control. The temporal Gauss
equation propagates when the spatial gauge equations and matter equation
hold. This supplies a classical consistency condition needed before a quantum
constraint construction. It does not already prove a quantum constraint
algebra, anomaly freedom or a continuum measure.

Twelve independent stable thermo coefficient channels can implement every
su(2) gauge variation, including around F=0. This uses CP's explicit Darboux
coefficient extraction; independently discarding every dχ term is not a lawful
non-Abelian gauge transformation. The coframe channels remain independent.

## How this joins the existing gravity sector

The internal generators commute with every spin Clifford operator. They also
commute with C, so they commute with the positive-Λ Cartan translation
coefficients C⊗Γ^a used in the corresponding LR enlargement. Thus both the
massive multiplier and that admitted geometric sector can coexist with the
new internal gauge action. This is a coefficient-level compatibility result,
not a quantum unification theorem or a selected positive cosmological constant.

Within the already admitted MG/LR coframe action, the new gauge action has
stress tensor

\[
T^A_{\mu\nu}=\frac1{g_*^2}
\left[B(F_{\mu\rho},F_\nu{}^\rho)
-\frac14q_{\mu\nu}B(F_{\rho\sigma},F^{\rho\sigma})\right].
\]

This is the conditional MG variation with the scalar product replaced by the
constant invariant B. In the orthonormal color basis it is a sum of its three
already-derived Maxwell metric variations, holding each F fixed during metric
variation. That proves the displayed source; the non-Abelian A variation is
the distinct NCG-6 theorem. Full independent coframe variations therefore give
the corresponding Einstein–gauge equation in that action class.

With spin matter, the independent spin connection still has the SM algebraic
torsion equation; its source is the enlarged module's total axial spin density.
The gauge kinetic term has no independent spin-connection source because F is
defined by exterior/gauge differentiation, not a torsionful affine curl.
Eliminating torsion must retain the resulting contact interaction. SM's special
eight-component Fierz formulas are not asserted for the new sixteen-component
field. The classical Hamiltonian's positive-quantum-energy problem is also not
removed merely by adding internal copies.

All of this is a **conditional action adapter**. Native group and curvature
identities have been derived; the base continuum, action/readout choice,
coframe sector, local coupling rule and numerical constants have not become
primitive consequences. In particular the positive gauge norm readout and LR's
oriented gravitational loop readout are distinct process choices.

## Verification and the next quantum target

[The pinned NCG certificate][cert] replays the native algebra and independently
checks the centralizer, invariant pairing, affine connection jets, variable
gauge transformations, action/source variations and off-shell current identity.
It reports 1,357 native checks, 1,073 independent rational checks, 15 word
replays and 23 rejected alternatives. The metric/torsion paragraph above is a
written conditional application of the cited prior variations; it is not a
fresh whole-MG/LR/SM rerun or a new physical certificate.

The physical gauge/state/action dictionary remains a separate quantum target.
Existing YM-12 supplies square-sourced transfer positivity.
[YM-50][ym50] now adds a proved compact reference/observable/transfer bridge
using NCG's quaternion sector and specified normalized Phi_Sigma record
readouts. Its coefficient recognition completion and symmetric native-turn
heat law represent the former SU(2) reference objects. This does not
identify the chain's interacting law with the full NCG quantum-field action.
The preceding YM-42 result advances the three-rail silence step and proves
why spatial contraction alone cannot establish the time gap. YM-43 now
supplies the separate arbitrary-source time-transfer bound on the full
heat-kernel chain, uniformly in finite spatial width, at three certified
coarse parameter cells. Its positive reference functional was admitted at
that stage; YM-50 now supplies the above scoped native construction.
YM-44 now constructs the full interacting time limit at every fixed finite
width, with an explicit operator-norm error on the declared linear trajectory.
YM-45 supplies the width/fine-time uniform gap on that full heat chain for
abs(theta)<1/1680 and kappa=theta a, and transports its positive rate to
every fixed-width time limit. YM-46 then constructs the positive infinite-volume
local-observable state and normalized time map at each fixed a, including the
actual square-sourced ordering, its gap and regulated reflection positivity.
YM-47 now joins those local vacuum history limits at arbitrary relative
cutoff rates, with both iterated limits equal. The positive history
functional retains the correlation gap and reflection positivity.
YM-48 now constructs a strongly continuous time action and self-adjoint
generator on the reflected history completion, retaining the all-source gap.
Its local coefficient row embedding and nonzero source do not yet identify
the full physical observable sector. YM-49 now constructs bounded time-zero
coefficient observable action, its product/dagger laws and ordered readouts.
The positive defect C(2t)-C(t)^2=L(t)^dagger L(t) characterizes closure
of the coefficient-history sector on its row. RKF's exact memory equation
has a gap-controlled remainder on the retained complement.
The actual chain's row defect remains unevaluated; finite closing/nonclosing
observers are controls, not a verdict on it. The nonzero time/observable
commutator is not automatically NCG gauge curvature.
The remaining distinct targets include this row-closure calculation,
physical gauge/state selection, general UGD integration and NCG quantum measure. The
[mission map](../PHYSICS_INGREDIENTS.md) links the commit-pinned proof/certificate
and tracks quantum matter, gravity constraints, the original Wilson E4D-C
problem and 4D continuum at their actual scopes. The small bridge window
remains a declared trajectory. YM-50's native record protocol and isotropic
clock are also specified choices, kept explicit in its origin ledger.
Raw endpoint records are not recognition-Cauchy despite convergence of
their scalar readouts. A compact representation theorem does not select
the physical state or construct quantum gravity.
Current YM verification uses Python 3.12 only.

[YM-51's tool audit][ym51] now isolates a remaining selection step.
Symmetric native protocols can preserve the quaternion brackets,
reference functional and all linear decay while quadratic decay varies.
The protocol second moment C is additional data. Covariance transforms
C with the frame; only a separately required invariance of the fixed
law forces C=cI. Six quadratic channels can recover general symmetric C.
The remaining scale is a heat-clock choice and changes the relative
interaction theta/c. This does not establish local gauge Ward identities,
physical spatial isotropy or a gap for the alternative interacting laws.
The mission map therefore tracks this selection gate separately from
the old fixed-protocol gap and the actual row-closure question.

[YM-52's energy and observer bridge][ym52] now derives the free protocol's
mean/record-noise and bounded-density entropy balances. Its weighted frame
brackets construct cof(C), with an exact centered rate
min(tr(C)/4,a+b) for ordered tensor eigenvalues a<=b<=c. Rank-two active
protocols can relax despite lacking one direct direction.
These raw order brackets remain distinct from the connection curvature
defined above.

A sign-insensitive observer has rate a+b. Maximizing this rate at fixed
trace uniquely selects isotropy and gives a sharp stability estimate;
maximizing the full-observer rate leaves many anisotropic maximizers.
The classical spectral formula is credited, and the native record
construction and observer criterion are explicit. The physical criterion,
energy/time calibration and transfer of these rates through interacting
chains are separate remaining steps. The full-observer objective cannot
be promoted into a unique physical isotropy principle.

[ncg]: https://github.com/Parveen117/Publications/blob/7b34bc657883272f58dc74e9e511ae5c568266ff/papers/native-compact-gauge/THEOREM.md
[cert]: https://github.com/Parveen117/Publications/blob/7b34bc657883272f58dc74e9e511ae5c568266ff/papers/native-compact-gauge/CERTIFICATE.json
[ym50]: https://github.com/Parveen117/Publications/blob/ef4a0b53295456b13c905e31359212905dd5de3d/papers/yang-mills-certified-benchmark/YM50_NATIVE_REFERENCE_HEAT_BRIDGE.md
[ym51]: https://github.com/Parveen117/Publications/blob/a0dee6e72f97a1cd03e652f853eacf763e7daf06/papers/yang-mills-certified-benchmark/YM51_TOOL_DEPENDENCY_AND_HEAT_SELECTION.md
[ym52]: https://github.com/Parveen117/Publications/blob/fa87843c10220ee344d490db87ac15d376c52506/papers/yang-mills-certified-benchmark/YM52_ENERGY_NOISE_AND_OBSERVER_GAP.md
