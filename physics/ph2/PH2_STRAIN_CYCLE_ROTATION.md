# PH2 — The return angle is the rotation left by a cycle of pure strains

Monty Dabas. 7 October 2026. Python 3.12. Continues PH1.

PH1 left one gate open: a physical process whose outcome is the return
angle Θ. This stage identifies the class of such processes and gives the
numbers for one of them.

## 1. Statement

Let H(t) be a closed path of positive responses and let F(t) be driven by
pure strains only,

```text
dF = D F ,      D = ½ F⁻ᵀ dH F⁻¹   (symmetric),      F₀ᵀF₀ = H₀ .
```

**T1.** FᵀF = H(t) along the path, and F⁻¹F₀ is the propagator of the
native response transport dw = −½H⁻¹dH·w (NT-3, RMG1).

**T2.** When H returns, F(T) = R·F₀ with R a rotation, and R is the
inverse of the return of RMG1: its angle is −Θ(γ).

**T3 (finite chain).** For samples H₀, H₁, …, H_N = H₀ put
M_k = (F_{k−1}⁻ᵀ H_k F_{k−1}⁻¹)^{1/2}, F_k = M_k F_{k−1}. Every M_k is
symmetric positive, M_N⋯M₁ is exactly a rotation R_N, and R_N → R. With
F₀ upper triangular the elements M_k do not change when the units of the
state variables are changed.

*Proof.* d(FᵀF) = Fᵀ(D + Dᵀ)F = dH. d(F⁻¹) = −F⁻¹D = −½F⁻¹F⁻ᵀdH·F⁻¹ = −½H⁻¹dH·F⁻¹,
which is the transport equation. At return F(T)ᵀF(T) = F₀ᵀF₀, so F(T)F₀⁻¹
is orthogonal; the transport return in the H-orthonormal frame is
F₀F(T)⁻¹ = R⁻¹. For T3, F_kᵀF_k = F_{k−1}ᵀM_k²F_{k−1} = H_k; under H → DHD
with D diagonal, an upper triangular F₀ becomes F₀D and M_k² is unchanged. ∎

So no element of the chain rotates anything, and the chain as a whole
does. The rotation is the return angle.

## 2. Where this is physics

A symmetric positive 2×2 matrix acting in sequence is:

- a partial polarizer (linear diattenuator) acting on the polarization of
  light — the chain is N partial polarizers at different azimuths, and R
  is a rotation of the plane of polarization with no rotator present;
- a pure stretch of a material element — the chain is a cycle of
  vorticity-free strains, and R is the net rotation of the element;
- a squeezing step of a mode pair, with R the residual phase-space rotation.

The classical name of R is the rotation left by composing non-collinear
boosts. That identification is stated from general knowledge; the
literature was not checked in this stage.

## 3. Numbers: the CO₂ cycle of PH1 as a chain

```text
elements N    rotation          weakest element (weak/strong amplitude)
     6        +25.70°           0.172
    12        +41.45°           0.287
    24        +49.40°           0.492
    48        +52.07°           0.694
    96        +52.81°           0.831
   384        +53.04°           0.955
  8000        +53.058°          0.998        (PH1: Θ = −0.926029 rad = −53.058°)
```

- The product is orthogonal to 10⁻¹⁵ at every N.
- Changing the entropy and volume units by (37, 0.004) changes no element
  (10⁻¹⁵) and no angle.
- Each finite N is its own exact prediction: twelve elements must give
  41.45°, not 53.06°.

**Twelve-element recipe** (strong-axis azimuth, weak/strong amplitude ratio):

```text
 1  108.22°  0.5216      5   23.89°  0.2872       9  167.09°  0.9367
 2   97.34°  0.7034      6   65.72°  0.3161      10   22.88°  0.6275
 3  163.60°  0.8980      7  106.87°  0.3479      11    9.21°  0.4982
 4    2.63°  0.5204      8  124.21°  0.5732      12  168.74°  0.4697
```

**Exact witness** (rational arithmetic): diag(2, ½) followed by diag(3, ⅓)
turned by the (3,4,5) angle is a strain times a rotation with
tan θ = 288/541; the same two strains aligned leave tan θ = 0.

## 4. Claim boundary

```text
RETURN ANGLE = ROTATION LEFT BY A PURE-STRAIN CYCLE              PROVED (T1–T3), checked to 1e-15
CO₂ CYCLE AS A 6…8000 ELEMENT CHAIN, UNIT-BLIND                  COMPUTED
A BENCH MEASUREMENT OF THE CHAIN                                 NOT DONE; feasibility not assessed
THE FLUID ITSELF PERFORMING THE CHAIN                            NOT IDENTIFIED
ANY CLASSICAL EQUATION MODIFIED                                  NO
```

The bench would measure the rotation of a programmed path; the fluid enters
through the program (c_v, pressure coefficient, speed of sound along the
cycle). A process in which the fluid's own responses compose in sequence —
so that the fluid performs the chain — is the remaining gate.

## 5. Reproduce

```text
pip install CoolProp numpy
python ph2_strain_cycle_rotation.py
python -m unittest test_ph2_exact
```
