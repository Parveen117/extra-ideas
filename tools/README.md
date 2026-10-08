# Tools

Mathematical tools kept from two uploaded drafts of "new mathematical objects" (23 objects), 9 October 2026.
An object is kept only in a form that is proved and checked here; what is wrong, undefined, or already
certified elsewhere is not rewritten. The [source audit](SOURCE_AUDIT.md) lists every object.

Evidence here is exact (integers, rationals, cut-complex rationals; directed rational enclosures where a
transcendental number enters). Python 3.12, standard library only. It is outside the frozen R1–R46 / MP1–MP2
register and its workflows. One workflow, `tools.yml`, replays the four certificates on Python 3.12.

| Tool | From the drafts' | Statement |
| --- | --- | --- |
| [RW1 — returned winding count](rw1/RW1_RETURNED_WINDING_COUNT.md) | winding equation dw/dh = K − w²/h (behind its "emotional potential") | The equation is ⟨n(n+ν)⟩ = y for the returned count of two strands. ½·mean < spread < mean in every sector. For the face ladder: 4c+3 < κ(1−r²)/r < 4c+4 at every coupling; for CT1's centre record: κ/2 < sinh 2k_c < 2κ/3. |
| [JR1 — jet reading of the change operator](jr1/JR1_JET_READING_OF_THE_CHANGE_OPERATOR.md) | "λ-residue" | On tails with poles the reading obeys Rd(C f) = λ·Rd(f). A simple pole carries no λ; a pole of order m needs m readings. On the draft's stated domain it is zero. |
| [PM1 — prime-turn series](pm1/PM1_PRIME_TURN_SERIES.md) | "prime lattice", "prime-modular form", "seam zeta function" | The draft's series are not kept (one diverges, the other has no seam law). The series built on PT1's prime turns has whole coefficients, an exact product over primes, and a seam law Θ(1/t) = t^(4k+1)Θ(t) whose phase +1 is fixed by the four units; the two sides agree to forty digits at the tested points. |
| [CY1 — cyclic form of three readings](cy1/CY1_THREE_READING_CYCLIC_FORM.md) | "cyclic tensor" | T − Tᵗ = d ln(I/K) ∧ d ln(J/K), unchanged by a common factor of the readings; the symmetric part changes, except for one special factor. |

What each tool is for:

- **RW1** gives two-sided bounds on the face-ladder ratio at every coupling (theorum/75-T1 is one-sided) and
  turns CT1-C3's numerical check into a proved inequality. Bears on: theorum/75, CT1, MG1.
- **PM1** is the next case after LAM-2: a series whose data per prime are exact turns and whose crossing phase
  follows from the four units. Bears on: PT1, the LAM line.
- **JR1** is an exact reading of a pole on which the change operator of GE1 acts as a number. No stage uses it
  yet.
- **CY1** separates the reference-free part of three readings from the part that depends on the reference. No
  stage uses it yet.

## Claim boundary

```text
RW1, JR1, CY1                         PROVED (written proofs ; exact checks)
PM1 P1, P3                            PROVED
PM1 P2, P4                            PROVED to n = 1500 ; classical in general
PM1 P5 (seam law)                     AGREEMENT TO 40 DIGITS CERTIFIED at the stated points ; equality classical, pinned
N1 / N2 / N3 , K₀ / L₀ , RH , 4D      NOT TOUCHED ; OPEN
```

Exact finite checks support the written proofs; the few float comparisons in the tests are cross-checks and
are marked as such. This is not formal or independent expert certification. The bounds of RW1 and the series
of PM1 are known mathematics (credited in each stage); the stages say what is the line's own.

## Reproduce

```text
python tools/run_all_tests.py          # every tool's tests
```
