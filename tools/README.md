# Tools

For the next application of DS1, see [RH → Yang–Mills](../physics/RH_YM_TRANSFER_MAP.md): the exact diagonal map, applicable theorems, remaining core-gap hypotheses and pinned sources.

Mathematical tools from two uploaded drafts of "new mathematical objects" (23 objects), 9 October 2026.
The drafts are read by the owner's rule — the seam is the diagonal, the 45° line — and matched against the
one law of the framework. [DS1](ds1/DS1_ONE_DIAGONAL.md) is that match. An object has a tool of its own only
in a form that is proved and checked here and not already certified elsewhere. The
[source audit](SOURCE_AUDIT.md) says what each of the 23 objects is in the law.

Evidence here is exact (integers, rationals, cut-complex rationals; directed rational enclosures where a
transcendental number enters). Python 3.12, standard library only. It is outside the frozen R1–R46 / MP1–MP2
register and its workflows. One workflow, `tools.yml`, replays the certificates on Python 3.12.

| Tool | From the drafts' | Statement |
| --- | --- | --- |
| [DS1 — one diagonal](ds1/DS1_ONE_DIAGONAL.md) | seam, cut J(z) = i·z̄, seam measure e^(−w²), and the owner's statement on the 45° line | F = R − D is the pairing of a record with its mirror. Seen = lost is the 45° line. Positivity is a count: n₋(R − D) = n₊(B − 1) = zeros of det(1 − zB) inside z = 1, at most the rank of the memory. A mirror form is non-negative iff every point is on the line; the critical-line mirror and the 45° mirror are one, in the chart z = (1 − ι)(s − ½). |
| [RW1 — returned winding count](rw1/RW1_RETURNED_WINDING_COUNT.md) | winding equation dw/dh = K − w²/h (behind its "emotional potential") | The equation is ⟨n(n+ν)⟩ = y for the returned count of two strands. ½·mean < spread < mean in every sector. For the face ladder: 4c+3 < κ(1−r²)/r < 4c+4 at every coupling; for CT1's centre record: κ/2 < sinh 2k_c < 2κ/3. |
| [JR1 — jet reading of the change operator](jr1/JR1_JET_READING_OF_THE_CHANGE_OPERATOR.md) | "λ-residue" | On tails with poles the reading obeys Rd(C f) = λ·Rd(f). A simple pole carries no λ; a pole of order m needs m readings. On the draft's stated domain it is zero. |
| [PM1 — prime-turn series](pm1/PM1_PRIME_TURN_SERIES.md) | "prime lattice", "prime-modular form", "seam zeta function" | The draft's series are not kept (one diverges, the other has no seam law). The series built on PT1's prime turns has whole coefficients, an exact product over primes, and a seam law Θ(1/t) = t^(4k+1)Θ(t) whose phase +1 is fixed by the four units; the two sides agree to forty digits at the tested points. |
| [CY1 — cyclic form of three readings](cy1/CY1_THREE_READING_CYCLIC_FORM.md) | "cyclic tensor" | T − Tᵗ = d ln(I/K) ∧ d ln(J/K), unchanged by a common factor of the readings; the symmetric part changes, except for one special factor. |

What each tool is for:

- **DS1** is the common structure: the same law, line and count in the framework's block, the physics line
  and the Riemann line, with the table that matches them row by row. Bears on: ONE_LAW, theorum/02, 24, 28.
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
DS1 D1–D7                             PROVED (finite ; exact checks) ; the match: each row certified where named
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
