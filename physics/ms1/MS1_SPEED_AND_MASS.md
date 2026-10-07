# MS1 — More speed, less mass: the state law, the flow law and the combination law

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues CL2.

Owner's statement for this stage: if the speed is higher the mass is
less; as the speed reduces, mass begins to combine.

Sources read before building: **PR1-T2** (speed² + (curvature/2)² = 1 for
a coin), **PR4-T1/T3** (currents, σ² = n²·4p(1−p)), **CL1-T2/T4**
(mass² = clock curvature; the coin that varies), **QC3-T3** (the
character law), **QC4-T1** (memory = variance).

## 1. The state law

For a state ψ = (ψ₁, ψ₂) with density n, current j, proper density σ and
share p of the first light-like reading:

**T1.**

```text
speed v = j/n = 2p − 1 ,          memory = 4p(1−p) = 1 − v² ,          n = σ / √memory .
```

The speed of a state is the imbalance between its two light-like
readings; its memory is what is left: **speed² + memory = 1**. A pure
reading moves at the full speed and carries no memory, no proper
density, no source (PR4). The balanced state is at rest and all memory.
This is PR1's law for a coin, now for a state.

## 2. The flow law

A static flow through a region of flip rate g(x) obeys ∂_xψ = g·Kψ
(from the sector law, since A·R = −K). Let G = ∫g dx be the accumulated
flip.

**T2.** The current j is constant along the flow, and for a light-like
inflow (a pure reading):

```text
density = cosh 2G ,    proper density = sinh 2G ,    speed = sech 2G ,    memory = tanh² 2G .
```

In general σ = j·√(1 − v²)/v. As the flow slows, proper density — the
source τ = gσ of PR4 — piles up; as it speeds up, it thins out.

```text
e^G        1       3/2       2        3         5
speed      1       72/97     8/17     9/41      25/313
σ          0       65/72     15/8     40/9      312/25
```

A pure reading entering a region with flip is turned into a shared
reading; the more flip it has crossed, the slower it is and the more
mass it carries.

## 3. The combination law

Composing coin turns adds rest turns. With half-turn parts
(E_i, O_i) = (cos θ_i/2, sin θ_i/2) and M_i = 2·O_i the clock-curvature
mass of CL1:

**T3.**

```text
speed          c₁₂ = c₁ c₂ − s₁ s₂ < min(c₁, c₂)
mass           M₁₂ = M₁ E₂ + M₂ E₁
mass defect    M₁ + M₂ − M₁₂ = M₁ (1 − E₂) + M₂ (1 − E₁) > 0
ratio t = O/E  t₁₂ = (t₁ + t₂) / (1 − t₁ t₂)
```

Combining always lowers the speed, and the combined mass is less than
the sum: each mass is reduced by the record defect of the other's half
turn. The ratio t adds by the circular twin of the velocity law
(u₁ + u₂)/(1 + u₁u₂) of PR1: velocities saturate at the speed limit,
mass ratios do not.

**T4 (N equal coins).**

```text
M_N = M₁ · χ_N(θ/2) ,      χ_N = sin(Nθ/2)/sin(θ/2) ≤ N ,      binding = M₁ (N − χ_N) ,      speed = cos Nθ .
```

The character of QC3 once more: N combined masses weigh less than N
masses by N − χ_N. For a coin of speed 1519/1681 ≈ 0.904:

```text
coins     1         2         3         4
speed     0.904     0.633     0.241     −0.198
mass      0.439     0.857     1.232     1.548        (N × one mass: 0.439, 0.878, 1.317, 1.756)
```

Four of them combined no longer propagate forward: the speed has passed
through zero.

## 4. What this says about gravity

CL2 asked whether histories through a region where the coin varies show
something like gravity. They do not, in this model. A region of larger
rest turn slows a flow and loads it with memory (T2); the speed and the
rest turn move in opposite directions there (PR1-T2), whereas PR3 showed
that a change of the unit of time moves them together. A varying coin is
a profile of mass, not a change of clock. Refused, with T2 as the exact
witness.

## 5. Certificate

T1 at four rational states (pure, balanced and two mixed). T2 for five
accumulated flips with a light-like inflow — generator, semigroup law,
(cosh 2G, 1, sinh 2G), monotone slowing — and two general inflows. T3 on
ten pairs of rational coins, each composite also run through PR1's walk
identity. T4 for six coins. Six tests; a circular generator in the flow
and a split half-turn in the tower are rejected.

## 6. Claim boundary

```text
STATE: SPEED² + MEMORY = 1 ; PURE READING LIGHT-LIKE, NO SOURCE          PROVED
STATIC FLOW: j CONSTANT ; SPEED = sech 2G FOR A LIGHT-LIKE INFLOW        PROVED
COMBINATION: SPEED DROPS ; MASS DEFECT ; CIRCULAR ADDITION LAW           PROVED
N COINS: M_N = M₁ χ_N ≤ N M₁ ; SPEED cos Nθ                               PROVED
A VARYING COIN ACTS LIKE GRAVITY                                         REFUSED
WHICH COINS COMBINE IN NATURE; ANY PHYSICAL BINDING ENERGY               NO
TIME-DEPENDENT FLOWS; SCATTERING OFF A REGION                            NOT BUILT
```

## 7. Reproduce

```text
python ms1_speed_and_mass.py
python -m unittest test_ms1_exact
```
