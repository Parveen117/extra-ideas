# OR1 — One rule for records: least record memory among neighbours gives the whole numbers and the motion; the field of gravity is a different kind of law

Monty Dabas. 8 October 2026. Python 3.12. Exact rational arithmetic on the rational block; sympy for identities
with symbols.

The line carries three premises that no stage derives:

```text
QC1 / LC1     a record must be silent: the return a whole number of turns              "premise"
MO1           free readings take the largest count in a varying field                  "assumed"
MC1           the cost between neighbours is quadratic, and memory is what is spread   "assumed"
```

The question put to this stage: are they one rule — conditions on one quantity — or not?

Sources read before building: Publications **GE2-T2** (record lift: mean S = Σ p_s U_s, memory M = I − S†S =
Σ p_s‖U_s f − S f‖², for independent records without feedback; T5: a deterministic record has M = 0), **QC4-T1**
(curvature of a mean of flat readings = −variance; share law), **QC4-T4** (Wilson density = record defect of a
loop arrow; field strength² = its memory), **QC3**, **QC1-T3**, **LC1**, **WQ1-W2** (recordable ⇔ area = 2πκ·n),
**CL1-T1** (clock curvature 2 − ζ − ζ̄), **CL2-T1/T2** (count of a leg = g·τ; straight history largest, in one
frame), **MO1** (circles, E, L; claim boundary), **PT1** (exact turns; P4: no exact turn repeats), **OA1-T3**,
**DG1**, **RC1**, **TH1**, **NC1** (curvature-square law refused by Mercury), **SE1**, RKF **theorum/28** §4, §9
(a claim about all K needs an explicit bound).

Carrier: the rational block z = a + bR, det z = a² + b² (PT1). On turns GE2-T2 reads M = 1 − det S.

## 1. The quantity

**A1, A2.** For turns u_s with weights p_s: M = 1 − det S = Σ p_s det(u_s − S). For two readings,
M = p(1 − p)·det(u_a − u_b): the share law. Checked exactly on PT1's turns.

## 2. Neighbours in repetition: closure

Let a record be shared by K successive returns of one turn u (which return it was is not kept).

**B1 — a turn that does not close builds no record.** For every K, exactly,

```text
det( 1 + u + … + u^(K−1) ) · det( u − 1 )  =  det( u^K − 1 )  ≤  4 .
```

The sum never grows. The seen part of the record is det S_K ≤ 4 / (K²·det(u − 1)) for every K — an explicit
bound, no limit taken — and the memory goes to 1. For PT1's prime turns, which never return (P4), this is checked
for K up to 300; for u₅ the sum stays below det = 5 forever.

**B2 — a turn that closes builds one.** If u^q = 1 the content-q record has sum K, S = 1, memory 0; the lower
contents sum to zero after q returns (QC1-T3).

**B3 — the cost of one repetition.** Between two successive returns the memory is ¼·det(u − 1) = sin²(Θ/2):
a quarter of CL1's clock curvature.

**B4 — least.** That cost is zero exactly at a whole number of turns, and those are its minima.

So QC1's premise is the statement: *a record shared by repetitions keeps exactly what has no memory between
them.* The owner's words — rational: closure; irrational: bounded — are B2 and B1.

**B5 — the rungs.** On RC1's rung p/q the in-out phase closes for contents that are multiples of q and for no
other: on 1/3, content 3 is the lowest content the record over circuits keeps (TH1-H4).

## 3. Neighbours among histories: motion

A reading with clock g on a history of count τ carries the turn Exp(g·τ·R) (CL2-T1). Let neighbouring histories
τ(ε) = τ₀ + a·ε + b·ε² share a record.

**C1, C2.** The memory between the neighbours ±ε is g²a²ε² to leading order; among three neighbours (2/3)·g²a²ε².
It vanishes to this order exactly when a = 0.

So MO1's assumption is the same statement: *a record shared by neighbouring histories keeps the one whose count
is stationary* — in any field, with no reference to one frame.

## 4. One function on MO1's circles

With T the far period of a circle and τ its count per circuit:

```text
dτ/dT = E                 the conserved energy of MO1 is the slope of the count               D1
W = E·T − τ = 2π·L        the transform of the count is the area of the round pair            D2
dW/dE = T
```

A record over circuits at fixed far energy turns by g·W per circuit, so by §2 it is kept when g·L is whole: WQ1's
whole area. The count is stationary for the motion (§3) and its transform is whole for the record (§2): two
conditions on one function.

**D3 — the carried direction against the clock.** OA1's turn per circuit is 1 − τ/T = W/T + (1 − E): the area
share plus the binding share; to first order 3x/4 = x/2 + x/4. The two closures coincide where the binding
vanishes: x = 1/2, DG1's diagonal.

## 5. Where the rule stops: the fields

```text
light        its law is of the memory kind: field strength² is the memory of the loop arrow      QC4-T4 (cited)
gravity      its law is not: on flat slices it is e₂(stretch) = 0, a ratio variance : mean² = d − 1    SE1
             and the least-memory (curvature-square) law was tried and refused by Mercury              NC1
```

MC1's cost gives the right exponent because, for a static memory, the time–time contracted curvature is minus
half its variance operator (CV1-T5); it is not the law's general form.

## What this settles

```text
before                                          after
QC1's premise, MO1's assumption                 one rule: a shared record keeps what has least memory (GE2-T2)
                                                    among repetitions  →  closure, whole numbers
                                                    among histories    →  stationary count
light's field law                               the same quantity on loop arrows (QC4-T4)
gravity's field law                             a different kind: a ratio fixed by equivalence (TP1, SE1)
```

Three premises become two: the record rule, and equivalence of local frames for the frame field.

## What is put in

- GE2-T2 with its hypotheses (independent records, no feedback) and the reading of repetitions and of
  neighbouring histories as such records.
- The rule itself — that what is recorded is what a shared record keeps — is the remaining premise. It is not
  derived here.

## What is not shown

- §3 gives a stationary count; that it is the largest is CL2's, in one frame.
- §2 is about a record folded on one period. A reading that is free to analyse the rates separately is not such
  a record, so RC1's rule on the X-ray rates (QP1) stays a tested rule and is not derived by B1–B5.
- The statements of §2–§4 are the known pair "stationary action, whole action" (general knowledge). The line's
  part: both are GE2's memory, with exact bounds on the rational block, and the place where gravity's field law
  leaves the pattern is named.
- Why the frame field obeys a ratio and the phase field a least memory is open.

## Claim boundary

```text
GE2-T2 ON TURNS: M = 1 − det S = VARIANCE ; SHARE LAW                               PROVED (exact)
NON-CLOSING TURN: SUM BOUNDED FOR EVERY K, det S_K ≤ 4/(K² det(u−1))                PROVED (identity; checked to K = 300 on PT1's turns)
CLOSING CONTENT: SUM = K, MEMORY 0 ; COST OF ONE REPETITION = sin²(Θ/2)             PROVED
MEMORY AMONG NEIGHBOURING HISTORIES VANISHES AT LEADING ORDER ⇔ COUNT STATIONARY    PROVED (series)
ON MO1's CIRCLES: dτ/dT = E , E·T − τ = 2πL , dW/dE = T                             PROVED (symbolic)
TURN OF A CARRIED DIRECTION = AREA SHARE + BINDING SHARE ; EQUAL CLOSURES AT r = 2r_s   PROVED
QC1's PREMISE AND MO1's ASSUMPTION ARE ONE RULE                                     SHOWN, given GE2-T2's record model
THE RULE ITSELF                                                                     PREMISE — not derived
GRAVITY'S FIELD LAW AS LEAST MEMORY                                                 REFUSED (SE1: a ratio law; NC1: Mercury)
```

## Reproduce

```text
python or1_one_rule_for_records.py
python -m unittest test_or1
```

## Later note (MM1, 8 October)

§5 calls gravity's field law a different kind. MM1 refines this: it is not a least variance of turns, but on flat
slices it is the same quadratic form e₂ as GE2's and IN1's memories, set to zero between the cuts, and TP1's
coefficients follow there from two no-memory conditions. Nothing above is changed.
