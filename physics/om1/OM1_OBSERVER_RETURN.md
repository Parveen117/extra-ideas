# OM1 — a computed observer with a certified full hidden floor and return

Independent continuation of CM1 → SG1, incorporating GC1 at
`b332a71894a87484015342197d519fe2f41fccdb`, 9 October 2026.

**The observer-return equation now runs on the actual free physical core.**
A computed Gaussian-polynomial reading has a proved floor on its entire
orthogonal complement. One or two trial directions in that hidden space give
a return enclosure with a complete residual bound. The resulting ground
windows improve GC1's:

| Turns | GC1 ground window | OM1 ground window | Full hidden floor for the OM1 cut |
|---|---|---|---|
| 3 | [5.1865, 5.1868] | **[5.186691205, 5.186742055]** | **5.52076320** |
| 4 | [7.9942, 8.0029] | **[8.002732159, 8.002896171]** | **8.00516987** |

The retained cut is the span of **one** computed reading ψ: P=|ψ⟩⟨ψ|.
Its polynomial has 67 coefficients for d3 and 102 for d4; coefficient count
is not retained rank. The missing source Qhψ consequently has rank one.
CM1's rank 26 for its entire 67-reading cut remains a different statement.
The full hidden dynamics is infinite dimensional in both constructions.

## 1. What GC1 adds and what OM1 uses

GC1 and SG1 independently reach the same gauge-singlet comparison threshold:
E1≥5.5210 (three turns), E1≥8.0060 (four). GC1 also uses the exact trial
variance to give Temple ground bounds, improves the displayed four-turn gap
lower to 0.0031, and proves monotonicity of
ρd=E0(d)/(d(d−1)^(1/3)). Its 42 checks and 8 tests were replayed here.

OM1 uses SG1's pinned, full gauge-sector E1 lower, GC1's trial-finding
procedure, and CM1's exact Gaussian action/metric. It does not repeat the
core-gap existence stage. Its new obligation is the true hidden floor for
an **actual computed reading**, followed by a source-resolved hidden solve.

Precision note for the shared record: GC1's prose rounds its angular
expression to 3.2451, but its executable certificate uses the valid lower
P0=3.245. The expression (9/4)·3^(1/3) is slightly below 3.2451, so the
rounded display must not be used as a lower bound. GC1's main results use
3.245 and are unaffected. OM1 uses SG1's recorded E1 thresholds directly.

## 2. Carrier and retained reading

The carrier and realization are exactly SG1's free bosonic core:

\[
h=-\tfrac12\Delta_C+e_2(CC^T),\qquad C\in\mathbb R^{3\times d},
\]

with its Friedrichs form on the SO(3) row-gauge-invariant space. All column
direction sectors remain in the physical target. Denote the simple ground
by ψ0, with E0 its energy, and use the SG1 bound E1≥z.

The computed reading is ψ=p(e1,e2,e3) exp(−w e1/2)/√N, where N is its
Gaussian norm. It is smooth, rapidly decaying, gauge invariant and in the
domain of every power of h used below. For d3 use degree 10 and w=8/5;
for d4 use degree 12 and w=9/5. The frozen coefficients in OM1_TRIALS.json
were proposed by GC1's six-step floating inverse iteration, then stored as
exact rational numbers. No floating solve is part of certificate replay.

Write

\[
\eta=\langle\psi,h\psi\rangle,\quad
v=\|(h-\eta)\psi\|^2,\quad
P=|\psi\rangle\langle\psi|,\quad Q=I-P,\quad
D=QhQ,\quad r=Qh\psi.
\]

D denotes the self-adjoint compression associated with the restricted closed
form. Since ψ is in Dom(h), the off-diagonal source is a bounded rank-one
map and the usual block/resolvent formulas apply for energies below D.
All pairings are in the declared Gaussian C-space metric. The approximate
reading ψ is never identified with ψ0.

## 3. A full hidden floor from a second-level count

Assume η<z≤E1. The nonnegative spectral operator (h−E0)(h−z) gives

\[
E_0\ge L:=\eta-\frac{v}{z-\eta}.
\tag{1}
\]

This is the Temple step used in GC1. Let p0=|⟨ψ0,ψ⟩|². The spectral
decomposition also gives η≥E0 p0+z(1−p0), hence

\[
p_0\ge\frac{z-\eta}{z-E_0}
\ge\frac{(z-\eta)^2}{(z-\eta)^2+v}.
\tag{2}
\]

For every unit f orthogonal to ψ, |⟨ψ0,f⟩|²≤1−p0. Consequently

\[
\begin{aligned}
\langle f,hf\rangle
&\ge E_0+(z-E_0)p_0\\
&\ge E_0+z-\eta\\
&\ge \delta:=z-\frac{v}{z-\eta}.
\end{aligned}
\qquad\boxed{D\ge\delta I.}
\tag{3}
\]

This is a bound for the **whole physical complement**, not a finite hidden
matrix. It does not assume that ψ is the exact vacuum. In particular, using
z itself as D's floor would be false: for h=diag(1,3,5) and
ψ=(2,1,0)/√5, E1=3 but D has eigenvalues 2.6 and 5. Formula (3) gives the
correct 2.6 in that example.

The useful condition δ>η is precisely v<(z−η)². GC1's degree-10 d4 trial
does not meet this condition. Raising that reading to degree 12 reduces its
variance enough to cross it. The exact data, with outward display rounding,
are:

| d | η, approximately | Variance upper | Vacuum overlap lower | Recorded δ lower |
|---|---:|---:|---:|---:|
| 3 | 5.186743366131 | 0.000079148825 | 0.999292091 | 5.52076320 |
| 4 | 8.002896205804 | 0.000002576549 | 0.788981948 | 8.00516987 |

The recorded δ is rounded down from (3) to eight decimal places. Formula
(3) also applies to any larger retained polynomial cut containing ψ, because
its hidden space is a subspace of ψ⊥: degrees at least 10 for d3 and 12 for
d4 inherit this floor. This does not supply CM1's much larger conditional
floors at test energies 6 and 9, or a new excited-rate lower near its upper.

## 4. The lost part returns to the retained axis

For E<δ, elimination of the entire hidden space gives the scalar equation

\[
F(E)=\eta-E-\Sigma(E),\qquad
\Sigma(E)=\langle r,(D-E)^{-1}r\rangle.
\tag{4}
\]

The unique zero of F below δ is E0. For E<η, the framework's normalized
retained/returned comparison is

\[
\mathcal M(E)=\frac{\Sigma(E)}{\eta-E},\qquad
F(E)=0\ \Longleftrightarrow\ \mathcal M(E)=1.
\]

This is the concrete seam for this declared observer: the retained reading
equals the hidden return. Its source, metric and energy dependence are all
specified. The TVSP→VTSP chart interpretation remains as recorded in SG1;
the present result supplies a core operator realization, not a thermodynamic
intertwiner inferred from a coordinate swap.

Choose a finite hidden trial Y for (D−E)Y=r, and keep the entire residual
R=r−(D−E)Y. Completing the square, as in Publications YM75 and the transfer
map, gives

\[
\Sigma(E)=T_Y+\langle R,(D-E)^{-1}R\rangle,
\quad T_Y=2\operatorname{Re}\langle r,Y\rangle-\langle Y,(D-E)Y\rangle.
\]

Thus

\[
T_Y\le\Sigma(E)\le T_Y+\frac{\|R\|^2}{\delta-E}.
\tag{5}
\]

The independent bare upper v/(δ−E) is also valid; the code takes the smaller
of the two upper bounds. It never drops a positive residual term.

For a driven equation with retained source f and hidden source g, the right
side is f−⟨r,(D−E)^−1g⟩. The same Y gives the certified approximation
f−⟨Y,g⟩, with error at most ||R||·||g||/(δ−E). Hence the source correction
can travel with the operator correction. Its numerical value requires an
actual specified hidden input g and its norm; none is invented here.

## 5. Exact source moments and two hidden directions

Let μk=⟨r,D^k r⟩. The trials are Y in span{r}, then in span{r,Dr}.
For size m=1 or 2, set

\[
S_{ij}=\mu_{i+j},\quad H_{ij}=\mu_{i+j+1},\quad
H^{(2)}_{ij}=\mu_{i+j+2},\quad b_i=\mu_i,\quad c_i=\mu_{i+1}.
\]

Solve (H−ES)a=b by exact positive LDL elimination. Then TY=a^T b and

\[
\|R\|^2=\mu_0-2a^T(c-Eb)
       +a^T(H^{(2)}-2EH+E^2S)a.
\tag{6}
\]

Every moment through μ4 is computed exactly from ⟨ψ,h^kψ⟩ through k=6.
Only p,hp,h²p,h³p are needed: symmetry supplies the other pairings. Projecting
after each h action removes the actual ψ component, so Dr is not silently
replaced by hr. All means use CR1's Gaussian recursion. The code separately
checks full-action symmetry and retained-source orthogonality.

Tests reconstruct the complete residual by direct projected polynomial
action in C-space and compare its Gaussian norm with (6). Nonorthogonal
changes of the hidden trial coordinates preserve all return bounds when
S,H,H²,b,c are transported together. Substituting the identity for S fails
the negative control.

## 6. Ground windows from certified Schur signs

Let Σ− and Σ+ be the bounds in (5). At each proposed lower endpoint a,
the exact test η−a−Σ+(a)>0 certifies E0>a. At each upper endpoint b,
η−b−Σ−(b)<0 certifies E0<b. The endpoints are checked separately; the
search does not assume that its residual-error estimate is monotone.

| d | Hidden trial directions | Ground lower | Ground upper |
|---|---:|---:|---:|
| 3 | 1 | 5.186649374 | 5.186742329 |
| 3 | 2 | 5.186691205 | 5.186742055 |
| 4 | 1 | 8.002547441 | 8.002896178 |
| 4 | 2 | 8.002732159 | 8.002896171 |

At the recorded test energies 5.186743 and 8.002896, the two-direction return
uppers are 0.000052168444 and 0.000175864667; the corresponding bare uppers
are 0.000236958198 and 0.001133111747. These are enclosures at those energies,
not assertions that the upper bounds equal the actual memory.

Combining the new ground windows with SG1's excitation lower and CM1's
excitation upper gives

| d | Gap lower | Gap upper |
|---|---:|---:|
| 3 | 0.334257945 | 2.392008795 |
| 4 | 0.003103829 | 2.971067841 |

The main development is the actual hidden-floor/return contract and tighter
ground enclosure. These gap windows remain broad because the excitation
lower has not been sharpened beyond the pair comparison.

## 7. What is discharged and what remains

| Obligation | Status |
|---|---|
| Hidden floor behind a computed approximate vacuum | PROVED by (1)–(3), on the full physical complement |
| Correct source and metric for the retained observer | EXACT; frozen rational polynomial and full Gaussian pairings |
| Hidden return, with all residual contributions | ENCLOSED by exact moments/solves and the true floor |
| Source dressing in a driven equation | Identity and norm-error bound; numerical source g remains an input |
| Computed reading equals the true vacuum | Not assumed; only the proved overlap lower is used |
| CM1's 67-dimensional cut has source rank 26 | Preserved; OM1 chooses a different retained cut |
| Sharp first excitation or exact full hidden inverse | OPEN |
| Nonconstant modes, compact-core matching, volume/cutoff uniformity | OPEN |
| Continuum four-dimensional Yang–Mills gap | OPEN |

The core scale is still g^(2/3)/L. Equation (5) now has a real core floor
instead of an unproved placeholder, but transferring it to nonconstant
field modes requires their actual operator, source and uniform estimates.

## Sources and replay

GC1 `b332a71894a87484015342197d519fe2f41fccdb` supplied the reviewed Temple
step, trial-finding method and shared progress. Its note and implementation
hashes are embedded in OM1_TRIALS.json. The GC1 physics branch was read and
tested in a detached checkout; its code is not copied into this packet.

CM1 is pinned at `de14f48c19cb5f84baac3e3d13e7de158317c738`; SG1 at
`9157f54275a19914f53fcdee37634a3523305c7c`. The existing
[transfer map](../RH_YM_TRANSFER_MAP.md), section 4B, identifies Publications
YM75 at `e75e453b2fdb082570aa24e25a10488ec5d87b41` as the source of the
residual-return identity. OM1 supplies its actual core inputs. These spectral,
Temple and Schur tools are standard mathematical interfaces, not claims of
new priority for their abstract identities.

OM1_RESULT.json stores full rational moments and SHA-256 pins of all local
proof inputs, trials, implementation, note and tests. Replay is Python 3.12,
stdlib only:

```sh
python3.12 -B physics/om1/om1_observer_memory.py --check
python3.12 -B -m unittest discover -s physics/om1 -p 'test_*.py'
```

`--write` intentionally regenerates the result; `--check` rejects a stale
result or source pin. Frozen trial coefficients are certificate data: they
are not regenerated by the verifier.
