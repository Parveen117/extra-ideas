# Signed compass arrows: handoff to the Yang–Mills spatial comparison

9 October 2026. Research clarification after YC6 (`0b3b9d2`), not a new
spectral certificate. The owner clarified that observation exchanges T and V
**and reverses their arrow signs**: the third quadrant goes to the first,
while the second and fourth stay in their quadrants. The following is an
exact chart realizing that pattern. It does not establish that a physical
measurement implements this chart.

## 1. The two observation maps have different arrow actions

For two typed, normalized arrow components let

\[
\Pi_2=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
J=-\Pi_2=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}.
\]

YC2 recorded the unsigned swap. The owner's current pattern is realized by
the signed swap \(a'=Ja\), that is \((x,y)\mapsto(-y,-x)\).

| Before | Unsigned swap | Signed swap |
|---|---|---|
| I: (+,+) | I | III |
| II: (−,+) | IV | II |
| III: (−,−) | III | I |
| IV: (+,−) | II | IV |

The axes are excluded from this quadrant table. An entire quadrant is mapped
as stated; an individual arrow need not stay at the same angle. The fixed
line of J is x=−y, whereas the fixed line of the unsigned swap is x=y.
J²=I and det J=−1: this is an orientation-reversing reflection, not a singular
map or a physical reversal of temperature or volume. Labels carry units;
the orthogonal matrix statement uses normalized components.

No assignment of a particular thermodynamic potential to the owner's drawn
quadrants is inferred without that drawing. The four identities below are
indexed by their potentials, independently of their position on a sketch.

## 2. All four Maxwell identities, before and after

Use the thermodynamic convention dU=T dS−P dV, commuting smooth derivatives,
and nonsingular local charts for each held-fixed derivative. Denote the
signed chart by

\[
(t,v,s,p)=(-V,-T,S,P).
\]

Here t and v are new oriented chart labels, not new thermodynamic temperature
and volume. A derivative and its held-fixed variable must both be carried
through this map. The chain rule gives:

| Potential | Before | After |
|---|---|---|
| U | \((\partial T/\partial V)_S=-(\partial P/\partial S)_V\) | \((\partial v/\partial t)_s=-(\partial p/\partial s)_t\) |
| F | \((\partial S/\partial V)_T=(\partial P/\partial T)_V\) | \((\partial s/\partial t)_v=(\partial p/\partial v)_t\) |
| H | \((\partial T/\partial P)_S=(\partial V/\partial S)_P\) | \((\partial v/\partial p)_s=(\partial t/\partial s)_p\) |
| G | \((\partial S/\partial P)_T=-(\partial V/\partial T)_P\) | \((\partial s/\partial p)_v=-(\partial t/\partial v)_p\) |

For example the F relation acquires one minus on each side, so those signs
cancel. Thus arrow locations change, but the equality signs of the four
transported identities do not acquire an extra physical defect. These are
thermodynamic Maxwell relations, not the dynamical Maxwell field equations.
The classical input is David Tong, *Statistical Physics*, section 4.4.2,
equations (4.176)–(4.179):
https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S4.html

For noncommuting cuts retain UP4/UP5's actual defect. A passive chart change
cannot discard \([D_S,D_V]\Phi\). Exchanging two operations reverses their
commutator sign; it does not set it to zero.

## 3. The central ST scalar travels with its conjugate pairing

On YC1's source-response carrier, the convention is instead

\[
d\Psi=T\,dS+P\,dV,\quad
(T,V,S,P)=(\Psi_\kappa,h,\kappa,\Psi_h),\quad
\Psi_{\rm mid}=\Psi-\tfrac12(ST+VP).
\]

The plus sign here is the source convention, already distinguished from
thermodynamic pressure in CA2. Carrying the same signed chart gives

\[
d\Psi'=-v\,ds-p\,dt,\qquad
ST+VP=-sv-tp,\qquad
\boxed{\Psi'_{\rm mid}=\Psi'+\tfrac12(sv+tp)=\Psi_{\rm mid}.}
\]

The final equality compares the corresponding point of the two charts.
For example \(\Psi'(s,t)=\Psi(s,-t)\). On YC1's h=0 line, t=0,
ST=−sv and the same quartic-leading scalar survives. Replacing ST by st
while keeping the old potential would change the paired observable rather
than transport it. No new quartic coefficient or gap follows from this map.

## 4. Mixed response: sign memory versus magnitude

If input and output response frames are both transformed by the same J,
then a general response block transforms as

\[
L=\begin{pmatrix}A&B_1\\B_2&C\end{pmatrix},\qquad
L'=JLJ^T=\begin{pmatrix}C&B_2\\B_1&A\end{pmatrix}.
\]

Consequently

\[
\chi'=\frac{\det L'}{L'_{11}L'_{22}}=\chi,
\qquad \bar B'=\bar B,
\qquad w'=-w,\quad w=(B_2-B_1)/2.
\]

The oriented turn reading changes sign; the determinant ratio and magnitude
of the lost channel do not. If only one frame changes, use
\(L'=J_{\rm out}LJ_{\rm in}^{-1}\); the congruence formula above must not
be applied to that different operation. This abstract response-frame
statement is not an identification of J with an operation on YM link fields.

## 5. Exact spectral handoff: the full hidden return is covariant

For a self-adjoint YM operator on a declared sector and orthogonal retained
cut P with Q=I−P, write

\[
H=\begin{pmatrix}A&r^*\\r&D\end{pmatrix},\quad
\Sigma(z)=r^*(D-z)^{-1}r,\quad S(z)=A-z-\Sigma(z),
\]

at a real z below the full hidden floor. All compressions, domains and the
coupling in r are included; this r is the full off-diagonal block, not YC6's
coupling-stripped source. Let a change of representation be a unitary
\(\mathcal U=U_P\oplus U_Q\), with domains carried along. Then

\[
r'=U_QrU_P^*,\quad D'=U_QDU_Q^*,\quad
\boxed{\Sigma'(z)=U_P\Sigma(z)U_P^*},\quad
S'(z)=U_PS(z)U_P^*.
\]

Proof: substitute the first two formulas and
\((D'-z)^{-1}=U_Q(D-z)^{-1}U_Q^*\). Thus the inertia used in YC6 and the
full spectrum are preserved. A general unitary also works if the retained
cut is transported to P'=UPU*; the block notation above uses the corresponding
retained/hidden frames.

In particular, replacing the entire source r by −r leaves its return
unchanged. Flipping only one hidden component without transforming D can
change cross terms and is not the same operator. The executable witness
retains a nondiagonal D to detect exactly this mistake. The signed compass
is therefore useful for tracking orientation and complete cross-source
memory; it cannot erase an interaction merely by changing its arrow.

## 6. Continuing Yang–Mills

The next spectral target remains YC6's **spatial two-cell source/return
comparison**, with shared-link gauge constraints, the actual vacuum pairing,
all interface cross sources and the full hidden remainder retained. Use
the signed rule above when constructing the corresponding response frames.
To improve a lower bound one must still estimate that actual return and
control the physical scale as cells are added. Neither an orientation
convention nor a static centre scalar alone provides the estimate.

Status: exact chart, chain-rule and unitary-covariance identities; **no new
numerical gap, spatial gluing bound or continuum result**. The physical
implementation of the signed observation rule remains an identification
to establish. No predecessor packet is edited.

Reproduce the focused algebra checks from the repository root:

```text
python physics/check_signed_compass_handoff.py
```

Read alongside UP1 (centre degree filter), UP4/UP5 (noncommuting defects),
CA1/CA2 (typed responses), YC1 (centre scalar), YC2 section 4 (unsigned
observation map) and YC6 (actual vacuum gap and source return), at base
revision `0b3b9d2`.
