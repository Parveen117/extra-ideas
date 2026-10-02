# R46 source reuse and physical-interface audit

This pass read the relevant sources below. It is not a full audit of every
repository or of the private device programme. No private source body is
included in this application.

| Source at pinned commit | Consumed result or decision |
|---|---|
| [RKF R2](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/research/recognition_return/r2/THEOREM.md) | Paired native update, positive finite-aperture inverses, all-tail enclosure and periodic response. The existing exact solver is imported unchanged. |
| [Publications NI](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-return-identification/THEOREM.md) | Existing two-reading identification, its symmetry/readout hypotheses and uncertainty boundary. |
| [Publications CR](https://github.com/Parveen117/Publications/blob/f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56/papers/native-critical-response/THEOREM.md) | Existing third polynomial and distinction between return and inverse-return pole. These are credited rather than reintroduced as R46 discoveries. |
| [Materials response](https://github.com/Parveen117/EMK-material-response-public/blob/5dfd7f2cbdcd7fe7559ae2334c495dea947fa175/docs/ASSUMPTIONS_AND_LIMITATIONS.md) | Current public descriptors and transport coefficients are screening constructions. The inspected v1 implementation defines the reciprocity asymmetry using its curvature score; that is not an independent experimental test of the native return law. |
| [Raman T02](https://github.com/Parveen117/Thermodynamics-Reproducibility/blob/118596980f3dfd72329a3966c867867bf51fb99a/results/T02_RESULT.md) | Existing thermal-cycle data establish the reported hysteresis/nonclosure in that analysis. They do not supply this paper's paired-component probe or port admittance. No new Raman fit or validation is claimed. |
| [T01 measurement contract](https://github.com/Parveen117/Thermodynamics-Reproducibility/blob/118596980f3dfd72329a3966c867867bf51fb99a/protocols/T01_MEASUREMENT_CONTRACT.md) | Reuse the discipline of separately declared observables, units, controls and uncertainty; do not manufacture an empirical identity from one shared constructed response. |

The source pins include the actual consumed public files. Historical source
statements keep their own scope. R46 adds an ideal electrical target matching
the native recurrence, not a claim that existing materials or Raman data have
already instantiated it.

## Established physical comparison

[OpenStax, University Physics 2, section 10.2](https://openstax.org/books/university-physics-volume-2/pages/10-2-resistors-in-series-and-parallel)
supplies the ordinary ideal-resistor series/parallel interpretation and power
balance. [MIT 6.091, lecture 4](https://www.ocw.mit.edu/courses/6-091-hands-on-introduction-to-electrical-engineering-lab-skills-january-iap-2008/1ebc7a3c1f59bab488ae4a8b41c31296_lec4a.pdf)
contains an established repeated resistor-ladder example. These references are
comparison/implementation lineage, not additional native primitive axioms.
The particular R2 control mapping is proved directly in R46; the general
resistor-ladder mechanism is not claimed as new physics.

## Premise labels

* **Existing native theorem:** R2's paired products, update and limit.
* **Constructed readout target:** positive finite path weights and balanced flows.
* **Derived in this target:** elimination, native-return equality, positivity,
  control mapping, finite error and target faithfulness.
* **Physical interface assumption:** ideal ohmic components in a steady-DC
  configuration with known units and no unmodelled loading.
* **Constructed fixture:** component values and the previous NI/CR rational example.
* **Not supplied:** measured third reading, selected natural source, full EMK
  hardware dynamics or a fundamental constant.
