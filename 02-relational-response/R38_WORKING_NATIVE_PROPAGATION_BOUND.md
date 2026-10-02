# R38 working target: native retained-information propagation upper bound

Research owner: Monty Dabas. Development started: 2 October 2026 (India).

Status: **WORKING / NOT CERTIFIED**.

This file starts the next certification cycle from the frozen R37 chain. It
does not promote a conjecture to a result. The purpose is to state the
strongest target that is compatible with the already certified negative
controls before a proof, certificate or README claim is written.

## Frozen facts that constrain the target

R35 proves that the exact microscopic support front is faster than the
sharp phase coefficient: a nonzero corner survives at component speed
sqrt(2). Therefore no theorem of the form

```
all nonzero support lies inside |x| <= n/sqrt(2)
```

can be certified for the existing source.

R36 derives, for its native leading-flow class,

[
\rho=\langle u,u\rangle,qquad j_i=\langle u,A_i u\rangle,
qquad j_x^2+j_y^2\le \rho^2/2,
]

and a native local conservation identity. R36 explicitly does **not** call
this the exact one-block current of the lattice source.

R37 constructs exact finite-support source packets with reliable
positive-threshold arrival. Their retained-information speeds approach
(1/\sqrt2) from below. R37 explicitly leaves a universal all-packet
ballistic upper theorem open.

These three facts must all survive R38.

## No-classical-premise rule

The R38 proof path may use only pinned native source roles, matching,
completion, finite count sums, source Laurent coefficients, native phase
evaluation already earned in R16/R35, the R36 observer, and exact R37
packet/detector definitions.

The proof path must reject as premises:

- Fourier completeness or a Fourier representation theorem;
- stationary phase, dispersive PDE or continuum finite-speed theorems;
- quantum-walk limit theorems or spectral theorems;
- probability/stochastic transport laws;
- Lorentz invariance, a physical light cone, rods or clocks;
- a physical energy condition;
- any imported theorem that supplies the desired propagation bound.

A comparison citation may be recorded later, but it cannot be a proof
dependency.

## Candidate admissibility definition

A **native reliable localized flight family** is a sequence of exact source
fields (u^{(L)}_n=Z^n u^{(L)}_0), finite at preparation, together with
native address cuts (D^{(L)}_n), such that:

1. (\|u^{(L)}_0\|>0) and source matching norm is the normalization;
2. the preparation diameter (w_L) and observation time (N_L) satisfy
   (w_L/N_L\to0);
3. for one fixed native rational threshold (0<\theta<1), independent of
   L, the detector retains at least theta of the matching norm squared at
   the declared arrival;
4. the detector diameter is (o(N_L));
5. its displacement divided by (N_L) has a native completed limit v.

The definition is a target construction, not a physical detector selection.
It intentionally excludes zero-threshold first-nonzero-support arrival,
because the certified R35 fast corner would otherwise be a counterexample.

Before certification this definition must be red-teamed for dependence on
detector shape, threshold, cancellation, multi-packet splitting and
preparation-dependent recentering.

## R38 target theorem

The desired statement is

[
\boxed{\mathscr L(v)\le1}
]

in the R36 response-dual count length, equivalently

[
\boxed{|v|\le1/\sqrt2}
]

in component Euclidean count units, for every admissible reliable localized
flight family.

Together with R37.6 this would imply sharpness:

[
\sup_{\text{admissible reliable localized flights}}|v|=1/\sqrt2.
]

This is **not yet proved**.

## Native proof route to test

The next algebraic task is not a moving-domain PDE argument. It is to
derive a finite-source displacement/moment identity directly from the
Laurent source.

For a finite field f define native address moments by finite count sums,

[
M_i(f)=\sum_x x_i\,\langle f(x),f(x)\rangle.
]

For one exact source block compute, coefficient by coefficient,

[
M_i(Zf)-M_i(f).
]

The candidate route is to rewrite this increment as a native pairing with
a source-derived displacement operator (V_i), plus any exact finite
boundary/interference residue forced by the Laurent coefficients. The
following questions must be answered by exact source algebra rather than
assumed:

1. Does the displacement pair (V_x,V_y) obey a sharp joint matching bound?
2. Is its long-flight average controlled by the already certified R35
   phase-gradient bound, or is an additional native residue present?
3. Can a fixed positive retained-information detector moving with velocity
   v force the normalized displacement moment to approach v despite
   arbitrary small fast tails?
4. Do split packets or oppositely signed role components defeat a first
   moment argument? If yes, replace the moment by a native truncated
   transport functional and certify that instead.

Only after these questions close may an R38 theorem be written.

## Mandatory negative controls

Any eventual R38 certificate must reject at least the following false
alternatives:

1. exact support front is (1/\sqrt2);
2. R36 leading current is silently the exact lattice current;
3. phase-slope bound alone proves an all-packet theorem;
4. zero-threshold arrival obeys the proposed bound;
5. exponentially small fast tails may be set to zero;
6. detector threshold may depend on L so as to chase the fast corner;
7. detector width of order N counts as localized;
8. physical c is identified by the native bound;
9. SI rods/clocks are fixed by source counts;
10. electromagnetic field selection has been proved.

R35's explicit fast-corner witness and R37's reliable packets must be
replayed unchanged.

## Certification gate

R38 remains WORKING until all of the following exist and pass:

- written native proof with explicit scope;
- exact source application reconstructing the new identities;
- adversarial packet and detector tests;
- derivation ledger with zero imported classical premises on proof paths;
- graph-mutation rejection;
- pinned source hashes;
- frozen R37 replay with old certificates unchanged;
- altered-source and altered-certificate rejection;
- explicit physical-selection flags remaining false.

No PASS status, README promotion or physical-constant claim is permitted
before that gate.
