# YC6 — vacuum gap, the full return step, and the physical scale gate

9 October 2026. Continuation of YC5 at `f1e6787`, preserving its frozen packet.
The research target remains Yang–Mills. The owner's proposed dimensionless
information flow motivates this stage; its useful scale question is retained
while its proposed identifications are checked separately.

**New finite-site result:** on the gauge-invariant, all-centre-even sector,
the actual compact operator has

\[
\boxed{\Delta_{\rm vac}(\theta)\ge1\quad(0\le\theta\le14).}
\tag{1}
\]

This is an internal vacuum-sector gap, not YC5's gap across all centre sectors.
The same complete retained vacuum cut is used. A second result is an exact
coupling-step identity and alternating response tower for the **full hidden
return**. Neither result is called a spatial renormalization map or a 4D theorem.

## 1. Exact carrier and a normalization correction

Keep

\[
\mathcal H_{\rm vac}=L^2((S^3)^3,d\mu)^{\operatorname{Ad}SU(2),\,Z_1,Z_2,Z_3},
\quad H_\theta=H_0+\theta V,
\quad H_0=-\sum_i\Delta_{S^3,i},
\quad V=2\sum_{i<j}|u_i\times u_j|^2,\quad0\le V\le6.
\]

The three centre flips are \(Z_i:U_i\mapsto-U_i\). Write
\(E_0\le E_1\le\cdots\) for this sector's eigenvalues, counted with
multiplicity, and \(\Delta_{\rm vac}=E_1-E_0\). No rotational or spatial
permutation restriction is added inside the gauge-singlet vacuum sector.

OL1's exact earlier normalization is

\[
A_{\theta_{\rm YM}}=\tfrac14H_{4\theta_{\rm YM}}-3\theta_{\rm YM}I,
\quad \Delta(A)=\tfrac14\Delta(H).
\tag{2}
\]

**Correction to YC5 section 1:** its shorthand \(H_\theta=4A_{\theta/4}\)
omitted the additive scalar \(3\theta I\). Its gap conversion and every
certificate sign in \(H\) units remain correct. The frozen YC5 packet is
not edited here. In (2), (1) reads
\(\Delta_{\rm vac}(A)\ge1/4\) for \(0\le\theta_{\rm YM}\le7/2\).
Equation (2) still supplies no physical length or energy calibration.

## 2. Vacuum certificate without dropping other vacuum directions

YC4's complete gauge-singlet free cut \(P\), with energy below 24 in this
sector, has rank 13. Its entire complement obeys \(QH_\theta Q\ge24Q\).
All four scalar-reflection blocks \(000,011,101,110\) are counted. YC5's
full-source resolution and its actual hidden-interaction refinement are
reused without changing a polynomial, matrix, norm or hidden floor.

For \(z<24\), let \(L_\theta(z)\preceq S_\theta(z)\preceq U_\theta(z)\)
be YC5's two-sided Schur forms. Exact negative inertia of \(L\) at most
\(j\) proves \(E_j\ge z\). At least \(j+1\) nonpositive directions of
\(U\) proves \(E_j\le z\). The positive hidden block supplies the
full-space congruence and reconstruction, so the lower is not a Ritz value.

For the larger excited upper thresholds, the tangent Schur upper need not
be useful below 24. We then use the ordinary retained Ritz form: at least
two nonpositive directions of \(PH_\theta P-zP\) prove \(E_1\le z\),
even when \(z\ge24\). Such an upper does not use a hidden inverse. Each
endpoint records which of these two valid upper arguments was used.

There are 113 rational nodes \(\theta_n=n/8\), \(0\le n\le112\).
At every node the replay proves \(l_j(\theta_n)\le E_j(\theta_n)\le
u_j(\theta_n)\), \(j=0,1\). Since \(V\ge0\), each \(E_j\) is
nondecreasing. A lower comparison estimate itself need not be nondecreasing.
Consequently use the already proved envelope

\[
\widehat l_1(\theta_n)=\max_{k\le n}l_1(\theta_k).
\]

On every full cell \([a,b]\),

\[
\Delta_{\rm vac}(\theta)\ge\widehat l_1(a)-u_0(b).
\tag{3}
\]

The minimum across all 112 cells is exactly \(2113/2000=1.0565>1\),
proving (1). This uses monotonicity of the ordered energies, not of their
difference, and never carries a later lower bound backward.

| \(\theta\) | Exact lower on vacuum gap | Exact upper |
|---:|---:|---:|
| 0 | 8 | 8 |
| 2 | 7.0791 | 7.0897 |
| 4 | 6.4900 | 6.6288 |
| 6 | 5.9228 | 6.5259 |
| 8 | 5.1212 | 6.7520 |
| 10 | 3.9868 | 7.3412 |
| 12 | 2.5826 | 8.3877 |
| 14 | 1.0565 | 11.3959 |

These terminating decimals are rationals. The loss of sharpness at large
\(\theta\) is not a proof that the actual vacuum gap tends to zero.
For example the early decrease from 8 to below 7.0897 is certified, while
CM2 proves eventual \(\theta^{1/3}\) growth. Global monotonicity of the
vacuum gap would therefore be a wrong assumption.

## 3. A full return step and its response tower

This theorem applies to any fixed finite retained cut commuting with \(H_0\),
on a declared sector. Define

\[
D_0=QH_0Q\ge\Lambda Q,\quad W=QVQ,\quad0\le W\le6Q,
\quad r=QVP,
\quad R_\theta(z)=(D_0+\theta W-z)^{-1},\quad z<\Lambda,
\]

where \(r\) maps the retained space to the entire hidden space. Matrix
representations use the same Haar pairing as the retained Gram form \(G\).
The return is \(\Sigma_\theta(z)=r^*R_\theta(z)r\), and the retained
form is \(S_\theta=T_\theta-\theta^2\Sigma_\theta\).

**Theorem YC6-S.** At fixed \(z<\Lambda\), for every \(h\ge0\),

\[
\Sigma_{\theta+h}(z)=\Sigma_\theta(z)
-h\,r^*R_{\theta+h}(z)W R_\theta(z)r.
\tag{4}
\]

For the proposed fourfold coupling step,

\[
\boxed{S_{4\theta}(z)=S_\theta(z)+3\theta M
-15\theta^2\Sigma_\theta(z)
+48\theta^3r^*R_{4\theta}(z)W R_\theta(z)r.}
\tag{5}
\]

Proof: multiply out the second resolvent identity
\(R_{\theta+h}-R_\theta=-hR_{\theta+h}WR_\theta\), then compress
with the complete source \(r\). Substitute in
\(S_{4\theta}=T_\theta+3\theta M-16\theta^2\Sigma_{4\theta}\).
The mixed resolvent product in (4) is self-adjoint after the displayed
compression, by the identity itself. No commutation of \(D_0\) and \(W\)
is assumed. In particular \(\Sigma_{\theta+h}\preceq\Sigma_\theta\).

The differentiated tower is

\[
(-1)^n\partial_\theta^n\Sigma_\theta(z)
=n!\,r^*R_\theta(z)[W R_\theta(z)]^n r\succeq0.
\tag{6}
\]

To see both positivity and a full remainder, put
\(A=R_\theta^{1/2}WR_\theta^{1/2}\ge0\) and
\(B=R_\theta^{1/2}r\). Then

\[
\Sigma_{\theta+h}=B^*(I+hA)^{-1}B
=\sum_{j=0}^n(-h)^j B^*A^jB
 +(-h)^{n+1}B^*A^{n+1}(I+hA)^{-1}B.
\tag{7}
\]

This is a finite exact identity for every \(h\ge0\); convergence of an
infinite geometric series is not required. Thus odd partial sums give lower
return bounds and even partial sums give upper return bounds, with the
entire omitted interaction retained in the final term. The magnitude of
that term is bounded in form by
\((6h/(\Lambda-z))^{n+1}\Sigma_\theta\). The bound may be coarse at a
large step, but its sign and domain are exact.

At \(\theta=0\), the first two moments in (7) are YC5's resolved
\(R(z)\) and interaction \(C(z)\). This is a typed response tower for the
actual hidden kinetic return. Identifying it with a particular native
\(\lambda\)-tower requires a further intertwiner; a shared word is not one.

Equations (4)–(7) retain hidden dynamics. They are **not a scalar closed
recurrence for the gap**: the new hidden inverse remains in (4). The packet
also encloses \(\Delta_{\rm vac}(4\theta)/\Delta_{\rm vac}(\theta)\)
at seven base couplings \(\theta=1/2,1,\ldots,7/2\), by dividing positive
certified gap intervals. These are finite coupling comparisons on the same
carrier, not spatial coarse-graining data.

From CM2, with \(d=2(e_1-e_0)>1.0112\),

\[
\Delta_{\rm vac}(\theta)/\theta^{1/3}\longrightarrow d>0
\quad\Longrightarrow\quad
\frac{\Delta_{\rm vac}(4\theta)}{\Delta_{\rm vac}(\theta)}
\longrightarrow4^{1/3}.
\tag{8}
\]

This is an asymptotic consequence of CM2, not an exact finite-coupling
power law and not a newly proved beta function.

## 4. The honest information-frame gap

Let \(\psi_\theta>0\) be the normalized actual ground of \(H_\theta\).
It is gauge invariant and centre even. Set
\(d\nu_\theta=\psi_\theta^2d\mu\). The exact ground-state transform gives

\[
\langle\psi_\theta f,(H_\theta-E_0)\psi_\theta f\rangle_\mu
=\int\sum_i|\nabla_i f|^2d\nu_\theta.
\tag{9}
\]

For smooth functions this follows by integration by parts, expanding the
gradient and using the ground equation; extension follows by form closure.
Multiplication by \(\psi_\theta\) is unitary from \(L^2(\nu_\theta)\)
to \(L^2(\mu)\) and preserves the declared gauge/centre restrictions.
Consequently

\[
\boxed{\Delta_{\rm vac}(\theta)=
\inf_{\substack{f\in\mathcal H_{\rm vac}\cap H^1(\nu_\theta)\\
\int f\,d\nu_\theta=0,\ f\ne0}}
\frac{\int\sum_i|\nabla_i f|^2d\nu_\theta}
{\int|f|^2d\nu_\theta}.}
\tag{10}
\]

Thus (1) is also a Poincare inequality for the **actual vacuum probability
measure**, with its kinetic metric and sector restriction. This makes the
information/probability interpretation precise. It does not compute that
unknown measure explicitly, identify a static response Hessian with the
kinetic form, or derive a gap from a dimensionless shape ratio alone.

A sharp control: a two-state reversible process with stationary probabilities
\((1/5,4/5)\), conductance \(2/5\), has gap \(5/2\). Multiplying every
rate by \(N^{-2}\) leaves those probabilities and their static information
unchanged, while the **dimensionless** gap becomes \(5/(2N^2)\to0\).
The generator's clock is extra data. Logarithmic coordinates or removal of
physical units do not remove dependence on dimensionless size \(N=L/a\).

## 5. What the proposed continuum argument does and does not establish

The following distinctions prevent a wrong no-go theorem.

**Normalized centre splitting.** CM2 proves only
\(\Delta_{\rm all}/\theta^{1/3}\to0\). The positive functions
\(\theta^{-1},1,\theta^{1/6}\) all satisfy this limit, with three different
absolute behaviors. CM2 does not select among them. YC5's widening error
intervals likewise do not determine an absolute tunnelling rate. Moving the
present computational target to the vacuum sector is useful; calling the
centre splitting already proved to die is unsupported.

**Homogeneous core.** TC1 proves, under a unitary dilation,

\[
-\tfrac{g^2}{2}\Delta+g^{-2}X\simeq g^{2/3}h.
\]

Refining bounds on \(h\)'s fixed pure numbers cannot change this exponent.
But the stronger claim “no logarithm can occur in a core calculation” is
false even in this repository: TC1-C5 already proves a logarithmic cutoff
term for a four-turn **integral**. That is different from a beta function
or the homogeneous spectral gap. Their carriers and targets must stay typed.

**Physical units and two different limits.** To state a conditional cutoff
claim without guessing a coupling convention, explicitly choose

\[
\mathcal E_{a,g}=\frac{\gamma g^2}{a}H_{\kappa/g^4}+c(a,g)I,
\qquad \gamma,\kappa>0\text{ fixed}.
\tag{11}
\]

This is a declared family, not an inferred calibration of the YM-line
parameter in (2). For example \(\gamma=1/2,\kappa=1\) matches TC1's local
quartic coefficients with cell unit \(a\). CM2 implies for (11)

\[
\operatorname{gap}_{\rm vac}(\mathcal E_{a,g})
\sim\gamma\kappa^{1/3}d\,\frac{g^{2/3}}a.
\tag{12}
\]

The physical gap in lattice units is therefore
\(\delta=a\operatorname{gap}(\mathcal E)=\gamma g^2\Delta_{\rm vac}(H)\),
not \(\Delta_{\rm vac}(H)\) itself. Both are dimensionless; they differ
by a coupling-dependent clock normalization.

If a proposed ultraviolet trajectory additionally has
\(g(a)^2\sim[2b_0\log(1/(a\Lambda_{\rm RG}))]^{-1}\), \(b_0>0\),
then (12) **diverges** as
\(1/[a\log(1/a)^{1/3}]\). It does not vanish.
More generally any fixed finite power \(g^p/a\) diverges on this trajectory:
write \(t=\log(1/(a\Lambda_{\rm RG}))\), obtaining a constant times
\(e^t t^{-p/2}\); its logarithm tends to infinity.

This rules out a finite positive continuum mass from the isolated family
(11) along that declared trajectory, using only homogeneous power-law
refinements. It does **not** rule out contributions from constant modes
inside the correctly matched field theory. A one-site lattice has physical
box \(L=a\); at \(a\to0\) its whole box shrinks. At fixed physical
\(L\), a continuum lattice instead has \(N=L/a\to\infty\). The
small-box effective core uses \(g_R(1/L)\), not an unqualified bare
\(g_0(a)\). Substituting the latter and obtaining \(g_0(a)^{2/3}/L\to0\)
does not establish a statement about the correctly matched effective theory.

For primary-source context, van Baal,
[QCD in a Finite Volume](https://arxiv.org/abs/hep-ph/0008206), section 3,
equations (9)–(10), explicitly distinguishes the overall \(1/L\) from
the renormalized \(g(L)\) dependence; its introduction explains the role
of the eliminated high-energy modes. We use its one-loop running as a named
conditional input, not as a nonperturbative mass-gap proof.

**The gap's scaling is not the coupling's logarithm.** If an actual physical
mass tends to \(m_*>0\), then its dimensionless gap satisfies
\(\delta(a)=a m_*+o(a)\). Under \(a\mapsto a/b\),
\(\delta(a/b)/\delta(a)\to1/b\): a canonical **power** law. Logarithmic
running belongs to the coupling. A lower bound \(K(a)\le\delta(a)\)
can vanish faster and still be true; to prove a positive limiting physical
lower one needs \(\liminf K(a)/a>0\), plus the field/limit construction.

At one-loop order let \(x=1/g^2\), so a fixed refinement gives
\(x'=x+2b_0\log b\). If \(\theta=\kappa x^2\), then

\[
\sqrt{\theta'}=\sqrt\theta+2b_0\sqrt\kappa\log b,
\qquad\theta'\ne4\theta\text{ in general}.
\tag{13}
\]

The exact one-loop solution has \(a'\Lambda_{\rm RG}=(a\Lambda_{\rm RG})^2\)
when \(\theta'=4\theta\); the required refinement factor depends on the
starting scale. A fixed coupling multiplication is not automatically a fixed
spatial blocking step. Actual step scaling must name a renormalized observable,
scale ratio and continuum extrapolation: see Luescher,
[Step scaling and the Yang–Mills gradient flow](https://arxiv.org/abs/1404.5930),
sections 1, 4 and 6. This is background for the missing interface, not evidence
that YC6 has constructed it.

No identification of “flat native \(\lambda\)” with a frozen RG flow is
assumed here. A chart's derivative tower, a kinetic response tower and a
spatial renormalization trajectory are distinct until their maps are proved.
The quoted YM98 \(\rho_N-2\beta\) inequality was not available in the
checked-out source set or located by the repository search in this turn;
it is not silently reinterpreted as (10).

## 6. The next field-theory target, without moving the present goalposts

For a many-cell lattice the relevant next statement would bound the full
source generated by eliminating one spatial scale, in the actual vacuum
pairing, with controlled dependence on \(N\), \(a\) and the chosen
coupling trajectory. In the language of (9)–(10), it needs a proved
comparison of the kinetic forms **and** vacuum measures under a declared
coarse-graining map. A positive ratio or a coordinate rotation alone gives
neither comparison. The required physical lower has order \(a\) in lattice
units on the correct trajectory; a nontrivial limiting field remains required.

YC6 supplies a finite-site vacuum certificate and a full-source coupling
identity that can be inputs to that task. It supplies no volume-uniform
reserve, no lattice beta function, and no 4D continuum construction. The
next concrete calculation is a spatial two-cell source/return comparison;
the quantitative vacuum bounds are a local input, not its conclusion.

## 7. Reproduce

```bash
PYTHONDONTWRITEBYTECODE=1 python physics/yc6/yc6_vacuum_step.py --check
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s physics/yc6 -p 'test_*.py'
```

Replay uses Python 3.12 standard-library rationals and YC5's unchanged exact
Haar engine. `explore_yc6.py` uses NumPy/SciPy only to propose rational
thresholds; it is not imported by the certificate. Tests include a genuinely
noncommuting resolvent step, the full alternating remainder at a large step,
an explicit ground-state/Dirichlet transform, constant information with
vanishing dimensionless gap, and controls on the limit inferences. The
operator and asymptotic statements have written proofs above; finite checks
are not formal proof-assistant verification or simulations of a 4D field.
