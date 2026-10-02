# Physics interface for the native compact gauge completion

Owner: Monty Dabas. 2 October 2026. Conditional application of Publications
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

The gauge/observable dictionary remains a separate quantum target. Existing
YM-12 already supplies square-sourced transfer positivity; a common SU(2)
algebra does not yet identify its measure and transfer with this NCG module.
The preceding YM-42 result advances the three-rail silence step and proves
why spatial contraction alone cannot establish the time gap. YM-43 now
supplies the separate arbitrary-source time-transfer bound on the full
heat-kernel chain, uniformly in finite spatial width, at three certified
coarse parameter cells. Its positive reference functional remains admitted.
YM-44 now constructs the full interacting time limit at every fixed finite
width, with an explicit operator-norm error on the declared linear trajectory.
The remaining distinct targets are the volume-uniform fine-step gap and the
native measure/NCG dictionary. The [mission map](../PHYSICS_INGREDIENTS.md)
links the commit-pinned proof/certificate and tracks quantum matter, gravity
constraints, E4D-C and the interacting continuum at their actual scopes.
Current YM verification uses Python 3.12 only.

[ncg]: https://github.com/Parveen117/Publications/blob/7b34bc657883272f58dc74e9e511ae5c568266ff/papers/native-compact-gauge/THEOREM.md
[cert]: https://github.com/Parveen117/Publications/blob/7b34bc657883272f58dc74e9e511ae5c568266ff/papers/native-compact-gauge/CERTIFICATE.json
