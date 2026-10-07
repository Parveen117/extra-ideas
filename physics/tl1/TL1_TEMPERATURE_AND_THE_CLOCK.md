# TL1 — Heat and the clock factor: the unit of fluctuation is carried by the clock

Monty Dabas. 8 October 2026. Python 3.12, standard library only. Exact
rational arithmetic; floats only in the illustration.

The owner's statement: energy is the bridge between the thermal side
and gravity. This stage makes that bridge exact with pieces already
certified.

Sources read before building: **R43.3** (exchange keeps the summed
readout for every preparation iff the unit factors are equal; eq. 43.7)
and its re-certification **DC1-K1**; **R43.7** (the common scale stays
free); **GR1-T2 / CL1-T2** (a reading at clock factor N counts g·N per
unit of static time); **PR3** (accelerated frame: the law is the
inertial one with ∂_t → (1/ρ)∂_η); **QC2** (fluctuation unit κ; weights
exp(−energy/κ)); **HB1**, **SC1**.

## Results

**L1 — a turn carried through static time.** A turn of local rate w_A at
clock factor N_A arrives at clock factor N_B with local rate

```text
w_B = w_A · N_A / N_B .
```

**L2 — units follow the clock.** Two ledgers at different places that
exchange through static time keep N_A E_A + N_B E_B. With local unit
factors u_A, u_B, R43.3's identity holds with b = u·N, so the summed
readout is kept for every preparation exactly when

```text
u_A N_A = u_B N_B .
```

DC1-K1 (one unit for frames that exchange) is the case of one place.
Between places the units are not equal; their products with the clock
factor are.

**L3 — no net flow.** Two fluctuation ensembles exchanging quanta of
static size e have zero net flow exactly when their weights per static
quantum agree, that is when

```text
κ_A N_A = κ_B N_B ,
```

and otherwise the flow runs toward the smaller κN.

**L4 — the law.** In equilibrium in a static field

```text
κ(r) · N(r) = constant .
```

Weak field: d ln κ / dh = g/c². Accelerated frame (PR3): κ·ρ = constant.

**L5 — heat has weight.** A ledger of local energy U at clock factor N
has static energy N·U; its pull is U r_s/(2r²N): that of a mass U/c².

## Numbers

```text
Earth's surface      1.09·10⁻¹⁶ per metre        (0.03 nK per km at 300 K)
Sun, surface to far  2.12·10⁻⁶
weight of one joule  1.09·10⁻¹⁶ N at the Earth's surface
```

## Reading

Temperature is not the same at every height in equilibrium; κ·N is.
What equalises between places that exchange is the unit of fluctuation
*as carried by the clock*. Thermal energy enters gravity as energy
(L5), and gravity enters heat through the clock factor (L4): the bridge
is energy, in both directions.

This is the known equilibrium law of temperature in a static field. The
route here is the framework's: R43.3 for units, GR1 for the clock, QC2
for the weights. The effect at the Earth's surface is far below present
thermometry; nothing new is predicted.

## Claim boundary

```text
TRANSPORT OF A TURN w_B = w_A N_A / N_B                             PROVED from GR1/CL1
UNITS: u_A N_A = u_B N_B ⇔ SUMMED READOUT KEPT                       PROVED (R43.3 with b = uN)
NO NET FLOW ⇔ κ_A N_A = κ_B N_B ; DIRECTION OF FLOW                  PROVED for the two-ledger exchange with geometric weights
κ N = CONSTANT IN EQUILIBRIUM                                       FOLLOWS from the above
WEIGHT OF A LOCAL ENERGY                                            PROVED from N² = 1 − r_s/r
RESPONSE CURVATURE OF MATTER AS A SOURCE OF GRAVITY                 NOT CLAIMED — bodies of different materials fall alike
A PREDICTION BEYOND THE KNOWN LAW                                   NONE
```

## Reproduce

```text
python tl1_temperature_and_the_clock.py
python -m unittest test_tl1_exact
```
