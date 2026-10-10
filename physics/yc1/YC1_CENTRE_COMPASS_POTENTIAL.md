# YC1 — the Yang–Mills centre record has an exact TVSP response potential

9 October 2026. Exact SU(2) Haar moments and formal series; the global bound
uses the proved RW1 winding inequality. The executable certificate is
`yc1_centre_compass_potential.py`; `test_yc1.py` contains the independent
focused tests.

The result is an exact symmetry theorem on a declared Yang–Mills carrier. One
SU(2) face supplies a gauge-invariant centre coordinate

\[
c(U)=\frac12\operatorname {Tr}U\in[-1,1],\qquad
c(gUg^{-1})=c(U),\qquad c(-U)=-c(U).
\]

Choose the centre-odd and centre-even observables \(c,c^2\), with sources
\(\kappa,h\). Their log-partition potential is

\[
\boxed{\Psi(\kappa,h)=
\log\int_{SU(2)}e^{\kappa c(U)+h c(U)^2}\,dU.}                 \tag{1}
\]

The literal CA2 substitution gives the pre-observation compass

\[
\boxed{(T,V,S,P)=
(\Psi_\kappa,h,\kappa,\Psi_h).}                              \tag{2}
\]

After observation the first two places exchange, with the sources and metric
carried through the chart:

\[
(T,V,S,P)\longmapsto(V,T,S,P)
=(h,\Psi_\kappa,\kappa,\Psi_h).                              \tag{3}
\]

Thus the owner's centre \(ST\) scalar is not an analogy on this carrier. It is

\[
\boxed{ST=\kappa\Psi_\kappa.}                               \tag{4}
\]

At \(h=0\), CT1 identifies it exactly as

\[
\Psi(\kappa,0)=\log\frac{2I_1(\kappa)}{\kappa},\qquad
\Psi_\kappa=\frac{I_2(\kappa)}{I_1(\kappa)}=\tanh k_c,
\qquad ST=\kappa\tanh k_c.                                  \tag{5}
\]

## 1. Carrier, symmetry and metric

| Required item | YC1 declaration |
|---|---|
| carrier | one SU(2) face with normalized Haar measure |
| symmetry | centre involution \(J:U\mapsto-U\), hence \(c\mapsto-c\) |
| target | response of the selected observables \((c,c^2)\) |
| source | \((\kappa,h)\) |
| metric | \(H=\nabla^2\Psi=\operatorname{Cov}(c,c^2)\) |
| intertwiner | \(\Psi(-\kappa,h)=\Psi(\kappa,h)\), induced by \(J\) |
| hypothesis | finite real \(\kappa,h\); normalized Haar measure |

The Hessian is positive definite. For every nonzero \((a,b)\), the polynomial
\(ac+bc^2\) is nonconstant on the full interval supporting the Haar law, so
its variance is positive.

The centre acts on the two selected channels by

\[
R=\begin{pmatrix}-1&0\\0&1\end{pmatrix}.
\]

At the fixed seam \(\kappa=0\), invariance gives \(R^THR=H\). Therefore

\[
\operatorname{Cov}_{0,h}(c,c^2)=0,
\qquad
\boxed{\chi_{\rm centre}(0,h)
=\frac{\det H}{H_{11}H_{22}}=1}                              \tag{6}
\]

for every finite \(h\). This is the precise sense in which the chosen
centre-symmetric compass has no lost mixed response. It does not say that one
face contains the full information of the Yang–Mills field.

At the Haar origin, the Catalan moments

\[
\langle c^{2n}\rangle=\frac{\operatorname{Cat}_n}{4^n}
\]

give

\[
H(0,0)=
\begin{pmatrix}
\operatorname{Var}(c)&\operatorname{Cov}(c,c^2)\\
\operatorname{Cov}(c,c^2)&\operatorname{Var}(c^2)
\end{pmatrix}
=\begin{pmatrix}1/4&0\\0&1/16\end{pmatrix}.                 \tag{7}
\]

Moving away from the seam turns on precisely the normalized lost channel:

\[
1-\chi_{\rm centre}
=\operatorname{Corr}(c,c^2)^2
=\frac{\kappa^2}{4}-\frac{\kappa^4}{16}
+\frac{41\kappa^6}{2560}+O(\kappa^8).                       \tag{8}
\]

## 2. The potential at the centre of the TVSP diagram

UP1 assigns a two-pair potential the diagram-centre reading

\[
\boxed{\Psi_{\rm mid}(\kappa,h)
=\Psi-\frac12\left(\kappa\Psi_\kappa+h\Psi_h\right).}       \tag{9}
\]

On the physical centre-source line \(h=0\), equations (4) and (5) give

\[
\boxed{\Psi_{\rm mid}(\kappa)
=\log\frac{2I_1(\kappa)}{\kappa}
-\frac{\kappa}{2}\frac{I_2(\kappa)}{I_1(\kappa)}
=\Psi-\frac12ST.}                                           \tag{10}
\]

This is the requested scalar function at the centre of the TVSP diagram.
It is gauge invariant and centre even.

Its weak-source series exposes what the centre does:

\[
\begin{aligned}
\Psi(\kappa)&=\frac{\kappa^2}{8}
-\frac{\kappa^4}{384}
+\frac{\kappa^6}{9216}
-\frac{\kappa^8}{184320}+O(\kappa^{10}),\\
ST&=\frac{\kappa^2}{4}
-\frac{\kappa^4}{96}
+\frac{\kappa^6}{1536}
-\frac{\kappa^8}{23040}+O(\kappa^{10}),\\
\boxed{\Psi_{\rm mid}(\kappa)}
&=\boxed{\frac{\kappa^4}{384}
-\frac{\kappa^6}{4608}
+\frac{\kappa^8}{61440}+O(\kappa^{10}).}                    \tag{11}
\end{aligned}
\]

The quadratic, Gaussian response is exactly silent at the diagram centre.
The first retained term has degree four. This is UP1's degree filter applied
to an actual SU(2) centre potential rather than a free formal polynomial.

## 3. A global theorem for the centre scalar

Put \(r=I_2/I_1=\tanh k_c\). Bessel recurrence and differentiation give

\[
r'=1-r^2-\frac{3r}{\kappa}.                                  \tag{12}
\]

Hence

\[
\Psi_{\rm mid}'
=\frac12(r-\kappa r')
=\frac12\left(4r-\kappa(1-r^2)\right).                      \tag{13}
\]

RW1 proves, for every \(\kappa>0\),

\[
\frac{\kappa}{2}<\sinh(2k_c)
=\frac{2r}{1-r^2}<\frac{2\kappa}{3}.
\]

Equivalently,

\[
3r<\kappa(1-r^2)<4r.
\]

Substitution in (13) proves the global bound

\[
\boxed{0<\Psi_{\rm mid}'(\kappa)<\frac12r(\kappa)
\quad(\kappa>0).}                                           \tag{14}
\]

Both \(\Psi_{\rm mid}\) and \(\Psi\) vanish at zero, while \(\Psi'=r\).
Integration therefore gives two equivalent scalar sandwiches:

\[
\boxed{0<\Psi_{\rm mid}<\frac12\Psi},
\qquad
\boxed{\Psi<ST<2\Psi}
\quad(\kappa>0).                                             \tag{15}
\]

So the centre scalar is positive at every nonzero coupling, grows
monotonically, and is quantitatively controlled by the uncut potential. This
statement uses no fitted constant.

## 4. What transfers to the Yang–Mills core

The exact transfer obtained here is limited but useful:

| YC1 face result | TC1/CM2 object | What is preserved |
|---|---|---|
| centre cut removes degree two | TC1 potential \(X=e_2(CC^T)\) has degree four | first visible degree |
| \(c\mapsto-c\) with even retained scalar | compact link centre flips; folded core is centre invariant | parity of the retained target |
| \(\chi=1\) on the centre seam | CM2 has the same leading core in all eight fixed centre sectors | symmetry-compatible leading blindness |

The middle row is a structural handoff, not an equality of operators. YC1's
\(c\) is a one-face trace coordinate. TC1's \(C\) is a flat constant-mode
matrix, and its quartic is \(e_2(CC^T)\). CM2's \(\theta\) is also not
identified with the one-face source \(\kappa\). No map in this stage sends
the coefficient \(1/384\) to the normalization of the TC1 core.

What YC1 adds to the Yang–Mills route is a declared positive two-source
response Hessian, which CA1 had left as the next requirement. It also explains
why a centre reading can be blind to the Gaussian layer while retaining a
quartic signal. It does not improve CM2's spectral gap bounds, determine the
absolute splitting of its eight centre-sector bottoms, or produce a
volume-uniform continuum gap.

The next measurable obligation is an intertwiner between this local
source-response observable and the compact multi-link Hamiltonian: introduce
the face sources in the one-site compact form, differentiate its actual ground
or heat-kernel potential, and bound the shared-link covariance error. Only
after that step can the YC1 scalar be used as a certified sector-splitting or
gap input instead of a local response theorem.

## Sources used unchanged

- CA2 for the literal two-source TVSP substitution and response ratio.
- UP1 for the potential's diagram-centre degree filter.
- CT1 for the SU(2) Haar moments, \(I_2/I_1=\tanh k_c\), and centre record.
- RW1 for the all-coupling two-sided bound on \(\sinh(2k_c)\).
- TC1 for the quartic free-core form and scaling boundary.
- CM2 for the fixed-centre-sector compact/core theorem and its unsolved
  all-sector splitting rate.

The result JSON pins these six source files by SHA-256.
