# DR2 — every radial level, every turn count, and the core gap scale

9 October 2026. Continuation of [DR1](../dr1/DR1_CORE_BY_ROWS.md), reviewed at
`28776c9` with index fix `035c943503e8298e822c2d39228148f9ae97328c`.
DR1's 20 checks and 8 tests were replayed before this development.
The publication parent is `8665cc810805abc249746ac15c71d9d335773879`, which
integrates the OM1/CB1 handoff `ad44e6a` into the physics branch. None of those
frozen packets is changed.

**Result.** On the free bosonic SO(3) gauge-singlet core, for every integer
number of columns (d\ge3), write (w=[2(d-1)]^{1/3}). Then

\[
\boxed{\frac{w}{64}\le E_1(d)-E_0(d)
\le w\left(\frac{11}{4}-\frac{15}{32(d-1)}\right).}                 \tag{1}
\]

There is a stronger explicit lower bound for (d\ge9), equation (16).
In particular the gap has order (d^{1/3}), with

\[
\liminf_{d\to\infty}\frac{E_1-E_0}{d^{1/3}}
\ge 2^{1/3}\frac{69}{1042}>0.0834,
\qquad
\limsup_{d\to\infty}\frac{E_1-E_0}{d^{1/3}}
\le 2^{1/3}\frac{11}{4}<3.4648.                                  \tag{2}
\]

This does not identify the limit or the lowest excited sector. DR1's stronger
small-d numerical bounds remain in force. The new coverage is all (d),
including an upper bound on the physical excitation for all (d\ge3).

## 1. Carrier and retained source

The carrier, measure, pairing and source remain

\[
\mathscr H_d=L^2(\mathbb R^{3\times d},dC),\quad
\mathscr H_d^{\rm phys}=\{f:f(GC)=f(C),\ G\in SO(3)\},\quad
h_d=-\tfrac12\Delta_C+X(C),\quad X=e_2(CC^T).
\]

All operators use their Friedrichs realizations. (E_0,E_1) are the first two
levels on this physical carrier, counting multiplicity. Column rotation
invariance is not imposed on the carrier. The lower bound covers all its
column direction sectors.

DR1 supplies the full-space form comparison

\[
h_d\ge K_d=\sum_{a=1}^3 k_a,\qquad
k_a=-\tfrac14\Delta_{r_a}+\tfrac{d-1}{2}|r_a|.                    \tag{3}
\]

The linearly confining (K_d) has compact resolvent, hence so does (h_d)
by the compact form embedding. Positivity improvement of the connected
elliptic Schrödinger semigroup gives a unique strictly positive ground state.
Uniqueness makes it invariant under row gauge and all column orthogonal
transformations; it is therefore the physical ground. These are the standard
analytic interfaces already used in TC1, CM1 and GC1.

The framework cut used here is the normalized Gaussian reading (P=|g\rangle
\langle g|). Its lost source is the actual (Q\bar h g=rg), computed below
in the same Lebesgue (L^2) metric. One multiple of that source is returned to
the trial. No change of metric, quotient or identification with an exact
vacuum is made. This is a two-direction variational use of the source, not an
evaluation of the full hidden inverse from OM1.

## 2. A closed floor for every radial level

After the row radial/angular decomposition in DR1, put (m=d+2l\ge3) and

\[
A_m=-\frac{d^2}{ds^2}+\frac{(m-1)(m-3)}{4s^2}+s,
\qquad u_d^3=(d-1)^2/16.
\]

The row levels are (u_d\nu_j(m)), with (j=0,1,\ldots) the radial index.
The origin uses the Friedrichs/regular condition. The extension of DR1 is

\[
\boxed{\nu_j(m)\ge \frac34(2m-1+6j)^{2/3}.}                     \tag{4}
\]

Here is a nodal form proof, rather than a positive-ground-state argument
applied incorrectly to an excited trial. Let

\[
p=(m-1)/2,\quad \kappa=\sqrt2/3,\quad
\alpha=(4p-2)/3,\quad x=2\kappa s^{3/2},\qquad
\phi_j=s^p e^{-\kappa s^{3/2}}L_j^{(\alpha)}(x).
\]

The generalized Laguerre equation is

\[
xq''+(\alpha+1-x)q'+jq=0.
\]

Differentiating and substituting it gives the exact differential identity

\[
A_m\phi_j=
\left[\frac{3\kappa}{4}(4p+1+6j)s^{-1/2}+\frac{s}{2}\right]\phi_j. \tag{5}
\]

It holds across the zeros as an identity of smooth functions; the apparent
poles from dividing by (L_j^{(\alpha)}) cancel. Set

\[
a=\tfrac{3\kappa}{4}(4p+1+6j).
\]

The minimum of (a/\sqrt s+s/2) is (3a^{2/3}/2), whose cube is

\[
\frac{27}{8}a^2=\frac{27}{64}(2m-1+6j)^2.                        \tag{6}
\]

Since (\alpha>-1), the Laguerre polynomial has exactly (j) simple positive
zeros. Require a form-domain function to vanish at their corresponding (s)
values. These are (j) independent interior trace conditions, making a
codimension-(j) form subspace. On every resulting interval, (\phi_j) has
constant sign. The ground-state representation, initially on smooth functions
supported away from the interval endpoints, says

\[
\langle f,A_m f\rangle=
\int\phi_j^2\left|(f/\phi_j)'\right|^2ds+
\int\frac{A_m\phi_j}{\phi_j}|f|^2ds
\ge \frac34(2m-1+6j)^{2/3}\|f\|^2.
\]

Closure extends the inequality to the restricted Friedrichs form. Thus at most
(j) eigenvalues lie strictly below this floor. Form min–max proves (4).
The Laguerre ODE and positive simple zeros are classical inputs; see
[DLMF 18.8, Table 18.8.1](https://dlmf.nist.gov/18.8) and
[DLMF 18.16(iv), especially 18.16.9](https://dlmf.nist.gov/18.16#iv).
No priority claim about this one-dimensional bound is needed.

For an exact algebra check, the polynomial can be normalized to (c_0=1):

\[
c_{k+1}=-\frac{j-k}{(k+1)(\alpha+k+1)}c_k.
\]

The recurrence proves the ODE coefficient by coefficient for general (j).
The executable checks representative degrees through 12 in addition to this
written identity; the finite checks alone are not the all-(j) proof.

## 3. Which count transfers to the gauge carrier

The physical carrier is contained in the subspace unchanged by all three
half turns that flip two rows. In this larger subspace the three row parities
are either all even or all odd. The row-separable (K_d) preserves this
parity subspace. It need not preserve the full row-gauge-singlet carrier.
We use inclusion followed by form min–max, not a restriction of (K_d) to a
claimed invariant gauge carrier.

The only possible level below the first parity excitation is the product of
three (l=0,j=0) row grounds. Consequently

\[
E_1(h_d|_{\rm phys})\ge u_d\min\left\{
2\nu_0(d)+\min[\nu_1(d),\nu_0(d+4)],\quad3\nu_0(d+2)\right\}.     \tag{7}
\]

Indeed an even-row excitation is radial or has (l\ge2); an all-odd product
has (l\ge1) in each row. This is exactly DR1's count. Let
(L(m)=\frac34(2m-1)^{2/3}). Equation (4) supplies
(\nu_1(d)\ge L(d+3)), and (L(d+4)>L(d+3)). We now have the closed floor

\[
E_1\ge u_d\min\{2L(d)+L(d+3),\;3L(d+2)\}.                      \tag{8}
\]

No numerical radial shooting is required for this all-(d) statement.

## 4. Return one actual source to the Gaussian

Write (n=d-1), (w=(2n)^{1/3}), (t=\sqrt w C). The unitary dilation is

\[
(U_w f)(t)=w^{-3d/4}f(t/\sqrt w),\qquad
U_w h_dU_w^{-1}=w\bar h_d,\qquad
\bar h_d=-\tfrac12\Delta_t+\frac{e_2(tt^T)}{2n}.                 \tag{9}
\]

For (g(t)=\pi^{-3d/4}e^{-e_1/2}), put

\[
\eta=\langle g,\bar h_dg\rangle=\frac{9d}{8},\qquad
r=\frac{(\bar h_d-\eta)g}{g}
=\frac{3d}{8}-\frac{e_1}{2}+\frac{e_2}{2n}.
\]

The exact source identities are

\[
\langle r\rangle=0,\qquad
v=\langle r^2\rangle=\frac{9d}{32n},\qquad
b=\langle rg,(\bar h_d-\eta)rg\rangle
=\frac{3d(25d-13)}{64n^2}.                                    \tag{10}
\]

All brackets on polynomials use (g^2dt). These formulas follow by Gaussian
integration, or by CR1's exact invariant moment recursion
(\langle P\rangle=\langle L_dP\rangle/(2\deg P)) for homogeneous invariant
monomials at width one. The verifier performs that recursion over
(\mathbb Q[d]); it does not interpolate a formula from sampled dimensions.
Clearing denominators with (N=nr) gives the polynomial identities

\[
\langle N\rangle=0,\quad
\langle N^2\rangle=\frac{9dn}{32},\quad
\langle N,(n\bar h_d-n\eta)N\rangle=\frac{3dn(25d-13)}{64}.
\]

For any real (a), the trial ((1-ar)g) has norm squared (1+a^2v) and
Rayleigh value

\[
\eta-\frac{2av-a^2b}{1+a^2v}.
\]

Use the single rational choice (a=1/4). It yields

\[
\boxed{\frac{E_0(d)}w\le\frac{9d}{8}-G(d),\qquad
G(d)=\frac{3d(23d-35)}{2(d-1)(521d-512)}>0.}                     \tag{11}
\]

The trial is gauge invariant and also invariant under column rotations.
It need not be everywhere positive to supply a ground upper bound. Its
source is exact; its return coefficient is a chosen variational coefficient,
not the exact hidden response. The plain Gaussian misses the positive
constant reserve (G(d)\to69/1042).

## 5. A positive count reserve for every dimension

For every (x\ge0), not merely asymptotically,

\[
(1+x)^{2/3}\ge1+\frac23x-\frac19x^2,                           \tag{12}
\]

because the third derivative is positive. The radial expression in (8),
divided by (w), is

\[
\frac{3n}{8}\left[2\left(1+\frac1{2n}\right)^{2/3}
+\left(1+\frac7{2n}\right)^{2/3}\right]
\ge\frac{9d}{8}-\frac{17}{32n}.                                \tag{13}
\]

The odd expression gives

\[
\frac{9n}{8}\left(1+\frac5{2n}\right)^{2/3}
\ge\frac{9d}{8}+\frac34-\frac{25}{32n}>\frac{9d}{8}
\quad(n\ge2).                                                 \tag{14}
\]

Thus (E_1/w\ge9d/8-17/(32n)). With (11),

\[
\frac{E_1-E_0}{w}\ge C(d)=G(d)-\frac{17}{32(d-1)}
=\frac{1104d^2-10537d+8704}{32(d-1)(521d-512)}.                 \tag{15}
\]

At (d=9+x), the numerator is (1104x^2+9335x+3295>0). In particular

\[
\boxed{E_1-E_0\ge [2(d-1)]^{1/3}C(d)>0\quad(d\ge9).}          \tag{16}
\]

For (d\ge11), (C(d)>1/64): after clearing its positive denominator, the
needed polynomial is

\[
1687d^2-20041d+16896
=1687(d-11)^2+17073(d-11)+572>0.                               \tag{17}
\]

For the remaining (d=3,\ldots,10), the verifier reads DR1's frozen rational
one-row floors and replays their exact sign certificates, including the full
ODE tails. Using (7) and the explicit trial (11), each normalized gap is
strictly larger than (1/64). The smallest of those eight certified margins
is greater than (0.1758), at (d=3). This completes the lower bound in (1)
without a dimension cutoff conjecture or an unverified large-(d) fit.

## 6. An excited trial with exact vacuum orthogonality

Let (t_i\in\mathbb R^3) denote columns and take

\[
q(t)=|t_1|^2-|t_2|^2,\qquad \chi=qg.
\]

This trial is row-gauge invariant. The exchange of columns 1 and 2 makes it
odd, while the unique positive ground is even. Hence its overlap with the
true vacuum is exactly zero. No approximate observer is used to deflate it.

The independent column radii (s_i=|t_i|^2) under (g^2dt) have moments

\[
\langle s_i\rangle=\frac32,\quad
\langle s_i^2\rangle=\frac{15}4,\quad
\langle s_i^3\rangle=\frac{105}8.
\]

Thus (\langle q^2\rangle=3). Since (q) is harmonic and homogeneous of
degree two,

\[
\frac{\langle qg,-\frac12\Delta(qg)\rangle}{\|qg\|^2}
=\frac{3d}{4}+1.
\]

Conditioned on radii, the angular mean of a column cross square is
(\frac23s_is_j). The radial means multiplied by (q^2) are
(45/4) for the pair ((1,2)) and for any pair with exactly one index 1 or 2,
and (27/4) for a pair using neither. There are (2d-3) pairs of the first
kind and ((d-2)(d-3)/2) of the second. Consequently

\[
\frac{\langle Xq^2\rangle}{\langle q^2\rangle}
=\frac{3d^2+5d-12}{4},\qquad
\boxed{\frac{E_1(d)}w\le\frac{9d}{8}+2-\frac1{2(d-1)}.}        \tag{18}
\]

The ground lower inherited from DR1 and (12) is

\[
\frac{E_0(d)}w\ge\frac{9n}{8}\left(1+\frac1{2n}\right)^{2/3}
\ge\frac{9d}{8}-\frac34-\frac1{32n}.                           \tag{19}
\]

Subtracting (19) from (18) gives the upper bound in (1). Together with (16)
it gives (2). The trial supplies an upper bound; it does not identify the
lowest excitation as belonging to this turning sector.

## 7. Numbers, reproducibility and claim boundary

The universal formulas prioritize coverage and transparency. They do not
replace the sharper finite ladders of CM1, OM1 or DR1. For example:

| d | Retained DR1 gap lower | DR2 simple-source gap lower | DR2 gap upper |
|---:|---:|---:|---:|
| 3 | 0.3342 | 0.27910 | 3.99331 |
| 4 | 0.4425 | 0.40023 | 4.71316 |
| 10 | 0.8097 | 0.78379 | 7.07055 |
| 11 | — | 0.04287 | 7.33742 |
| 20 | — | 0.13369 | 9.16249 |
| 100 | — | 0.35645 | 16.00072 |
| 1000 | — | 0.82771 | 34.63037 |

All lower decimals round down and upper decimals up. The apparent drop after
d=10 comes from switching proof estimates; it is not a claim that the actual
gap drops. The smaller universal (w/64) lower bound holds throughout.

Run from the repository root:

```sh
python physics/dr2/dr2_all_turn_gap.py --check
python -m unittest discover -s physics/dr2 -p 'test_*.py'
```

`DR2_RESULT.json` binds this note, code, tests and imported inputs by SHA-256.
The packet has 22 exact/outward checks and 8 unit tests.
The executable uses only Python stdlib and exact fractions. Symbolic identities
in (d), shifted positive coefficients, and outward integer cube roots
support the written proof. The small-d sign replays use recorded floors, not
a new floating search. Tests also compute the turning trial directly in
Cartesian coordinates, independently of the invariant moment generator.

| Item | Status |
|---|---|
| Every radial index (j\ge0), (m\ge3) | Closed nodal form bound (4); classical Laguerre/min–max inputs named |
| Every integer (d\ge3), physical free core gap | Explicit two-sided bounds (1), improved lower (16) |
| Large-d gap scale | Order (d^{1/3}); exact coefficient and convergence not established |
| Small-d sharper floors | DR1/OM1/CM1 retained unchanged |
| Actual Gaussian lost source | Computed exactly; one explicit trial correction, no full hidden inverse |
| Lowest excited sector | Not identified |
| Nonconstant modes, volume/cutoff uniformity, 4D continuum | Not obtained here |

Here (d) counts free-core columns/turn directions. It is not the number of
lattice cells or the spatial volume. Restoring the TC1 torus units still
gives (g^{2/3}(E_1-E_0)/L); fixing (d,g) does not remove its (1/L)
factor. CB1's compact-link statements remain separate. A useful next physical
target remains a local interaction/source-return estimate whose reserve
survives cell growth and the required nonconstant-mode matching.
