# GR1 — Gravity as memory: the static dictionary, and what G/c² is in it

Monty Dabas. 7 October 2026. Python 3.12. Exact rational arithmetic for
the theorems; ordinary floats only for the illustrative table.

Owner's thesis for this stage: gravity is the memory of lost information
that cannot be recovered; an observer does not read change as such —
change appears only because of lost information; otherwise everything
tends to the uncut potential, to recovery. Perhaps a relation between c
and G comes out.

Sources read before building: **MS1-T1/T2** (state law; static flow),
**PR4-T1/T3** (dual residue τ; σ² = n²·4p(1−p)), **PR3-T2** (a clock
factor multiplies the whole generator), **SY1-S6/S7** (record variance;
null state tensor), and the information-invariance ledger's **D3′**
(emergent gravity withdrawn as a derivation; retained only as a
one-function-sector programme) and **Theorem F** (no length from c and G
alone).

## 1. The local law: change is flip × √memory

Static flow, light-like inflow, accumulated flip G (MS1-T2). Put
N := v = j/n and m := memory = 1 − N².

**T1.** On any flip profile g(x):

```text
∂_x n = 2 g σ ,      ∂_x σ = 2 g n ,      ∂_x j = 0 ,        hence        ∂_x N = − 2 g N √m .
```

The rate at which the reading changes from place to place is the flip
rate times N times the root of the memory. **With no memory nothing
changes, whatever the flip rate** (certified: g = 7, m = 0 ⇒ ∂N = 0).
This is the thesis in exact form for this model: change is visible only
through lost information.

## 2. The dictionary

Read N as the rate of a static clock against a far one, and take as the
**named physical input** the Newtonian law in three space dimensions,

```text
m(r) = r_s / r ,        r_s = 2 G_N M / c² .
```

**T2.** Then, exactly:

```text
N = √(1 − r_s/r)                        the static clock factor
n = 1/N                                 blueshift of the inflow
p = (1 + N)/2                           share of the incoming reading;  memory = 4p(1−p) = r_s/r
acceleration / c² = ∂_r N = r_s / (2 r² N)
flip profile       g(r) = β / (4r(1 − β²)) ,  β² = r_s/r        (half the gradient of the fall rapidity)
horizon            m = 1  ⟺  N = 0  ⟺  p = ½ : the balanced share, all memory
```

The clock factor, the blueshift and the static acceleration are the
known ones for a point mass. The horizon is the balanced share of MS1 —
the state that is at rest and all memory — and it is not reached at any
finite accumulated flip (certified).

Illustration (standard constants):

```text
                         memory r_s/r        share p             acceleration
Earth surface            1.39·10⁻⁹           1 − 3.5·10⁻¹⁰       9.82 m/s²
Sun surface              4.25·10⁻⁶           1 − 1.1·10⁻⁶        274 m/s²
neutron star (1.4, 12km) 0.345               0.905               1.6·10¹² m/s²
```

## 3. c and G

In the dictionary the gravitational potential is Φ = −½c²·m, and

```text
2 G_N / c²  =  m · r / M  =  1.485·10⁻²⁷ m/kg .
```

G/c² is the length over which a unit mass produces total memory: the
conversion between mass and lost information at a distance.

**A relation between c and G alone is refused**, for an exact reason:
c and G do not combine into a pure number or a length (Theorem F of the
information-invariance ledger). Any relation needs a third quantity — a
mass, a length, or a unit of action. With a unit of action the third
quantity is the mass at which a sector's memory at its own clock
distance reaches one; that is the known Planck mass, and nothing here
derives it.

## 4. What is derived and what is put in

```text
∂N = −2 g N √m ; NO MEMORY ⇒ NO CHANGE                                  PROVED (model)
N² + m = 1 ; SHARE p = (1+N)/2 ; HORIZON = BALANCED SHARE               PROVED (model)
STATIC CLOCK FACTOR, BLUESHIFT, ACCELERATION OF A POINT MASS            REPRODUCED given m = r_s/r
m = r_s/r (THE 1/r LAW) AND THE CONSTANT 2G/c²                           INPUT — not derived
A RELATION BETWEEN c AND G                                              REFUSED (dimensional)
MOTION OF BODIES, LIGHT BENDING, WAVES, FIELD EQUATIONS                 NOT TOUCHED
MORE THAN A ONE-FUNCTION STATIC SECTOR                                  NO — as D3′ of the ledger requires
```

The stage does not contradict the ledger's withdrawal of emergent
gravity: it is static, carries one function, and takes the Newtonian law
as input. What it adds is the reading of that one function: the clock
factor is the imbalance of the two light-like readings, and the
potential is half c² times their memory.

## 5. Next gate

Why memory would fall as 1/r. In three space dimensions that is the
statement that the gradient of memory has no sources outside matter. The
2×2 block holds one or two space dimensions; the question needs three.

## 6. Reproduce

```text
python gr1_memory_dictionary.py
python -m unittest test_gr1_exact
```
