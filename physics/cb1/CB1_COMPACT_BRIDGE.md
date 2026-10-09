# CB1 — compact return at every coupling, and a correlated weak-side lift

9 October 2026. Independent continuation of OL1 at 10d80e9, using the
three-turn reading frozen in OM1 at 1c8a886. Python 3.12, standard library.

**The one-site vacuum-sector gap now has an explicit positive lower bound at
every finite coupling.** The correlated compact trial also improves the useful
weak-side bounds: their starting point falls from 10⁷ to **2·10⁶**, the bound
at and beyond 10⁷ rises, and the asymptotic lower coefficient more than doubles.

| Target, in OL1's H units | Certified result |
|---|---|
| One-site vacuum sector, every finite θ≥0 | Δvac ≥ (32/27) exp(−3θ) > 0 |
| All one-site gauge sectors together, every finite θ≥0 | Δgauge ≥ (54/125) exp(−6θ) > 0 |
| Vacuum sector, every θ≥2·10⁶ | Δvac ≥ 0.11 θ^(1/3) − 15/2 > 0 |
| Vacuum sector, every θ≥10⁷ | Δvac ≥ 0.36 θ^(1/3) − 15/2 > 0 |
| Vacuum sector, θ→∞ | liminf Δvac/θ^(1/3) > 0.6685 |

The exponential bridge is intentionally reported as such: it can be extremely
small in the former undecided interval. It proves positivity there, not a
useful order-one margin. OL1's stronger small-coupling certificate remains
available on [0,15/4], with least recorded margin 0.3098. The weak asymptotic
statement is a **liminf**; no existence of a limiting spectral ratio is assumed.

For a finite SU(2) lattice with L links and P plaquettes, the same positive
return gives Δ ≥ 2(3/5)^L exp(−2θP) in the corresponding H units. This is an
explicit finite-volume theorem whose bound decays with size. It supplies no
volume-uniform or continuum mass gap.

## 1. Carrier, target and normalization

Keep OL1's exact operator on normalized Haar measure on (S³)³:

\[
H_\theta=-\sum_{i=1}^3\Delta_i+V_\theta,
\qquad V_\theta=2\theta\sum_{i<j}|u_i\times u_j|^2,
\quad U_i=(a_i,u_i)\in S^3.
\tag{1}
\]

Use its Friedrichs realization. In the first and third through fifth rows
above, the target is invariant under simultaneous conjugation and each
antipodal centre flip U_i↦−U_i. The second row includes all centre parities,
while still imposing gauge invariance. A centre-odd state is orthogonal to
the unique positive vacuum, so the second row bounds its energy **above the
vacuum**. It does not compute the spacing between two levels inside that odd
sector, nor give a uniform lower bound as θ→∞.

The compact carrier has discrete spectrum; its strictly positive heat kernel
gives a unique positive ground state. Since the symmetries commute with H and
preserve positivity, that ground is gauge invariant and even under all three
flips. We establish kernel estimates before gauge restriction, avoiding any
assumed regular coordinate system on the gauge quotient.

The YM-line conversion is unchanged:

\[
\theta=4\theta_{\rm YM},\qquad A=H/4-3\theta_{\rm YM},
\qquad\Delta(A)=\Delta(H)/4.
\]

The improved weak-side starting point is θYM=500,000. The global bounds become
(8/27)exp(−12θYM) in the vacuum sector and (27/250)exp(−24θYM) for all gauge
sectors together. The stronger OL1/CB1 weak bounds apply only to the vacuum
sector; they must not be transferred to the centre-odd floors.

## 2. A retained vacuum ray and a bounded positive return

Here is the spectral interface, proved directly for continuous kernels. Let
T=exp(−tH), with a continuous kernel satisfying 0<m≤K(x,y)≤M, m<M, on a compact
connected carrier. Let ψ0>0 be its normalized ground eigenfunction and λ0
its largest eigenvalue. Its actual-vacuum transform is the Markov kernel

\[
P(x,dy)=\frac{K(x,y)\psi_0(y)}{\lambda_0\psi_0(x)}\,d\mu(y).
\tag{2}
\]

This is the retained vacuum frame. The full nonconstant return acts on
functions modulo constants, measured by oscillation sup f−inf f. We do not
replace ψ0 by a computed trial. For two rows x,x′, their likelihood ratio has
maximum divided by minimum at most (M/m)²: the unknown ground factors cancel
from this ratio of ratios.

For two probability measures with likelihood ratio in [a,b], a≤1≤b,
their total variation is at most (b−1)(1−a)/(b−a). With b/a≤r², this is at
most (r−1)/(r+1). One can enlarge b to ar² and verify the last step by

\[
\frac{r-1}{r+1}-\frac{(r^2a-1)(1-a)}{a(r^2-1)}
=\frac{(ra-1)^2}{a(r^2-1)}\ge0.
\]

Consequently the entire nonconstant return has oscillation contraction
at most q=(1−β)/(1+β), β=m/M. For any excited eigenfunction, φ/ψ0 is bounded,
nonconstant and has eigenvalue exp(−t(E−E0)) under P. Thus

\[
\Delta\ge\frac1t\log\frac{1+\beta}{1-\beta}
\ge\frac{2\beta}{t}>0.
\tag{3}
\]

The second inequality follows by integrating 2/(1−s²)≥2 from 0 to β.
This is the Hopf/Birkhoff contraction mechanism, expressed here through an
elementary row-overlap proof so that its infinite-dimensional use is explicit.
It is a different typed return from OM1's energy-dependent Schur memory:
no equality between their ranks, norms or kernels is asserted.

## 3. The actual compact kernel and the all-coupling bridge

For one unit S³, the free heat kernel relative to Haar probability is

\[
p_t(x,y)=\sum_{n=0}^\infty (n+1)e^{-n(n+2)t}U_n(x\cdot y),
\qquad |U_n(s)|\le n+1.
\tag{4}
\]

Here U_n is the Chebyshev polynomial of the second kind, equivalently the
SU(2) character. The eigenvalues are n(n+2), with multiplicity (n+1)².
This normalization agrees with the elliptic kernel in Baudoin–Bonnefont,
equation (3.2) and Lemma 3.3; their main subelliptic operator is not used.

On the centre quotient, average p_t(x,y) and p_t(x,−y): only even n remain,
still with constant term one. Exact positive Taylor sums and a geometric
bound for the entire omitted series give

\[
\begin{array}{lll}
t=1, & |p_t-1|\le0.202172332505<1/4,&\text{all parities};\\
t=1/2,& |p_t^+-1|\le0.16499435716<1/5,&\text{centre even}.
\end{array}
\tag{5}
\]

The certificate uses exp(−x)≤1/(Σ_{j=0}^{40}x^j/j!). After four character
terms, the omitted-term ratio is bounded by
((n+s+1)/(n+1))² exp(−ts(2n+s+2)), with s=1 or 2. This bound decreases
through the tail and is strictly below one. No finite heat sum is treated as
the full kernel, and the potential does not need a harmonic cutoff.

For (1), 0≤Vθ≤6θ: |u_i×u_j|²≤1, and the upper endpoint is achieved by
three orthogonal unit vectors. Positivity and the bounded-potential product
formula give

\[
e^{-6\theta t}p_t^{\otimes3}(x,y)
\le K_\theta(t;x,y)\le p_t^{\otimes3}(x,y).
\tag{6}
\]

Indeed each positive multiplication factor exp(−tV/n) lies between
exp(−6θt/n) and one; compare the product kernels and pass to the semigroup
limit. The same argument works on the centre quotient because V is even.
The compact smooth kernel makes these initially operator/a.e. bounds
pointwise. This is also the usual bounded-potential Feynman–Kac domination.

Equations (5)–(6) give β≥(3/5)³exp(−6θ) on the full carrier at t=1,
and β≥(2/3)³exp(−3θ) on the centre quotient at t=1/2. Applying (3) and
restricting to the gauge target proves the first two rows of the result table.
Adding a constant to V changes λ0 but cancels from P and from the gap.

More generally, on any specified finite SU(2) link lattice use
H=−Σ_{e=1}^L Δ_e+θΣ_{p=1}^P(1−Wp). Since −1≤Wp≤1, the same comparison
has β≥(3/5)^Lexp(−2θP). The positive ground is gauge invariant, so (3)
proves the finite-lattice bound stated above. Independent link flips need
not be symmetries on this larger lattice; the centre-even improvement is
used **only at one site**. At cubic volume N, L=P=3N, this particular reserve
decays exponentially with N. It identifies the cost of the global comparison,
not a closing of the actual thermodynamic gap.

## 4. Transporting the correlated free reading to the compact links

For the quantitative weak improvement we use the free h=−ΔC/2+X(C) on
C∈R^(3×3), with X=e2(CCᵀ), and the frozen normalized reading
ψ(C)=p(e1,e2,e3)exp(−w e1/2)/√N, w=8/5, degree 10, 67 coefficients.
These are the exact d3 coefficients frozen by OM1, proposed originally by
GC1's inverse iteration. The verifier repeats their Gaussian moments without
rerunning a floating search. They give

\[
\begin{split}
\eta=\langle\psi,h\psi\rangle&<5.186743367,\\
T=\tfrac12\int|\nabla\psi|^2&\simeq3.457825883684,\\
T_1=\tfrac12\int|C|^2|\nabla\psi|^2&\simeq13.442409536225,\\
E_{\rm w}=\int |C|^2\bigl(\tfrac12|\nabla\psi|^2+X|\psi|^2\bigr)
&\simeq20.561089049263.
\end{split}
\tag{7}
\]

The displayed latter three decimals are approximate; all computations use the
exact fractions. For the weighted kinetic integral the Gaussian action gives
T1=⟨ψ,|C|²(−Δ/2)ψ⟩+9/2. A separate gradient-polynomial test verifies the
9/2 term. Dropping it understates the compact error.

Let χ(C) be one for |C|≤R, linear down to zero for R<|C|<R+b, and zero
outside. It is an admissible form-domain cutoff with |∇χ|≤1/b. Put
ε=θ^(−1/6) and, near each of the eight central configurations, define

\[
U_i=(\pm\sqrt{1-|u_i|^2},u_i),\qquad u_i=\epsilon c_i,
\qquad\Phi_\theta(U)=\chi(C)\psi(C).
\tag{8}
\]

Assume ε²(R+b)²<1. The reading vanishes near every equator and extends by
zero; all eight patches are included. It is invariant under gauge rotations
and independent column sign flips, hence lies in OL1's vacuum sector.
This is a concrete map of the trial, source carrier and pairing, not an
identification of the full flat and compact spectra.

On each link the inverse metric is I−uuᵀ, and the measure is
du/√(1−|u|²). After rescaling, the product density is
Jε(C)=∏i(1−ε²|c_i|²)^(−1/2). The common chart normalization and eight-patch
factor cancel in the Rayleigh quotient. The exact form quotient divided by
2θ^(1/3) is

\[
\frac{\int J_\epsilon\left[
 \tfrac12\sum_i\bigl(|\nabla_i(\chi\psi)|^2
       -\epsilon^2|c_i\cdot\nabla_i(\chi\psi)|^2\bigr)
       +X|\chi\psi|^2\right]dC}
 {\int J_\epsilon|\chi\psi|^2dC}.
\tag{9}
\]

The negative inverse-metric correction can be discarded for an **upper**
bound. The density and cutoff terms cannot be discarded. For A=R+b,
s=ε²A²<1, put J*=(1−s)^(−1/2), c*=J*³/2. Since
∏i(1−ε²|c_i|²)≥1−ε²|C|²,

\[
1\le J_\epsilon\le J_*,
\qquad J_\epsilon-1\le c_*\epsilon^2|C|^2.
\]

Let τ=M_m/R^(2m), where M_m=∫|C|^(2m)|ψ|² is an exact Gaussian moment.
Markov's bound gives tail mass outside R at most τ. For any a>0, Young's
inequality applied to ∇(χψ) bounds (9) by

\[
G(\theta;R,b,m,a)=
\frac{\eta+aT+c_*\epsilon^2(E_{\rm w}+aT_1)
 +(1+1/a)J_*\tau/(2b^2)}{1-\tau}.
\qquad E_0(H_\theta)\le2G\theta^{1/3}.
\tag{10}
\]

Every omitted contribution is now displayed: outside mass, cutoff derivative,
weighted energy, metric and normalization. This bound decreases with θ when
the other parameters are fixed. Its input ψ need not be positive or equal
to the true vacuum.

## 5. Certified weak windows and the improved asymptotic coefficient

OL1's comparison proof works for any 0<α<π/2 with
4√θ α³σ≥9³, σ≤sin α/α. Its two class-level sign certificates and angular
floor then give, on the whole vacuum sector,

\[
E_1(H_\theta)\ge C_\alpha\theta^{1/3}-15/2,
\quad C_\alpha=(2\cdot2.3381+4.0879)(2\sigma^2)^{1/3}.
\tag{11}
\]

CB1 replays both GC1 sign certificates and the angular-floor inequality.
The change from OL1's α=0.4 is certified through the same reach condition,
not by extrapolating a numerical eigenvalue. All sines and square/cube roots
are bounded outward with rational or integer arithmetic.

| θ start | α | R | b | m | Young a | Cα lower | 2G upper | Available gap coefficient lower |
|---|---|---|---|---|---|---:|---:|---:|
| 2·10⁶ | 0.513 | 4.75 | 0.25 | 18 | 0.000442 | 10.721070381114 | 10.607127221291 | 0.113943159824 |
| 10⁷ | 0.39 | 5 | 0.25 | 20 | 1/6000 | 10.856091718607 | 10.492954743968 | 0.363136974639 |

The cutoff-tail uppers are respectively 7.7113647543·10⁻⁸ and
1.0825438242·10⁻⁸. Rounding the coefficient down to 0.11 / 0.36 gives the
uniform weak inequalities in the first table. At their starting points the
gap lower bounds are respectively **6.359131548843** and **70.059648841147**.
The reach condition and all compact errors improve with θ, so these are
interval theorems through infinity, not isolated spot checks.

For the asymptotic claim, first take θ→∞ for a fixed cutoff and Young a.
Then increase R (for example with b=R) and send a→0. Gaussian-polynomial
tails and all needed moments are finite, so (10) implies
limsup E0(Hθ)/θ^(1/3)≤2η. Independently take α→0 after the θ→∞ limit
in (11), giving liminf E1(Hθ)/θ^(1/3)≥8.7641·2^(1/3). Therefore

\[
\liminf_{\theta\to\infty}\frac{\Delta_{\rm vac}(H_\theta)}{\theta^{1/3}}
\ge 8.7641\,2^{1/3}-2\eta
>\boxed{0.6685}.
\tag{12}
\]

The exact outward evaluation is at least 0.668587341114. This improves OL1's
product-reading coefficient 0.3271 by retaining correlations already present
in the flat core. It is not a new claim that the θ^(1/3) scaling itself was
unknown. Nor does (12) prove full norm-resolvent convergence or identify the
exact first excitation.

## 6. Claim ledger and next obligation

| Obligation | Status after CB1 |
|---|---|
| Gap somewhere in the former one-site undecided interval | Explicit positive bound at every finite coupling; often extremely small |
| Centre-odd energies separated from the vacuum | Explicit positive lower at each finite θ; no weak-coupling uniformity |
| Correlated flat reading used on actual compact links | Explicit form-domain map, Haar measure, inverse metric and full error bound |
| Useful weak-side vacuum estimate | Starts at 2·10⁶; strengthened again from 10⁷; liminf coefficient >0.6685 |
| Arbitrary fixed finite lattice | Explicit positive bound decaying exponentially with links/plaquettes |
| Comparable quantitative margins throughout intermediate coupling | OPEN; the heat-return bridge alone is too small |
| Uniform interaction/return control as volume increases | OPEN; replacing the global potential oscillation cost is necessary |
| Continuum Yang–Mills mass gap | OPEN |

The classical ingredients are positive elliptic semigroups, Hopf/Birkhoff
contraction, Rayleigh–Ritz and a localized trial with exact moments. CB1's
contribution within this program is their explicit carrier/sector adapter,
constants, verified error ledger and new coupling coverage. Finite-volume
positivity by itself is a standard compact-operator fact, not a solution of
the Yang–Mills problem. The quantitative route toward that problem must
replace the extensive factor exp(−2θP)(3/5)^L with a controlled local
interaction or source-return bound that survives volume growth.

## Sources, analytic inputs and replay

- OL1: 10d80e926c96e1a6b9edbdff7ac6baa9b9c133bc on
  physics-ph1-ph2-response-holonomy; its 19 checks and
  7 tests were replayed. Earlier packets are not rewritten by CB1.
- OM1: 1c8a8860fb70d177619ebe7fdcb3bb437ac9ae08, specifically its d3
  frozen reading. CB1_TRIAL.json records the original file's SHA-256. The
  larger OM1 code/packet is not required to replay this derived certificate.
- CR1/CM1: exact Gaussian means and the full invariant-polynomial action.
  GC1: the two half-line sign certificates and angular floor used by OL1.
- [Baudoin–Bonnefont, arXiv:0802.3320](https://arxiv.org/abs/0802.3320),
  equation (3.2), Lemma 3.3: the elliptic S³ kernel and normalization.
- [W. Han–G. Han, arXiv:1906.04875](https://arxiv.org/abs/1906.04875):
  finite-matrix Hopf inequality and Birkhoff spectral-ratio context. Section 2
  above supplies the continuous-kernel argument actually used here.

The analytic imports are compact elliptic spectral theory/regularity,
positivity and bounded-potential semigroup product convergence, harmonic
completeness on S³, and the form min–max principle. The written proof uses
these named interfaces; finite arithmetic checks do not purport to prove
them. CB1_RESULT.json pins the local proof inputs, note, code, trial and tests.

~~~sh
python3.12 -B physics/cb1/cb1_compact_bridge.py --check
python3.12 -B -m unittest discover -s physics/cb1 -p 'test_*.py'
~~~

The --write option intentionally regenerates the packet; --check refuses a
stale result or input hash. The all-coupling exponential factors remain symbolic,
so large θ cannot silently underflow to a numerical zero in the certificate.
