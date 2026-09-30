# R3: Complex action on two real seam coordinates

R3 chooses a real two-mode space with Aghora `A=diag(1,-1)`, a seam vector `u=(Gamma,1)`, and its image `v=Au=(Gamma,-1)`, for finite real nonzero `Gamma`. These two vectors form a basis.

With the **added** rule to retain `u` and remove `v`, the unique projection `B` satisfies `ABA=I-B`. Consequently `K_Gamma=[A,B]` has `K_Gamma u=v`, `K_Gamma v=-u`, and `K_Gamma^2=-I`.

Write a state as `a u+b v` and assign it the complex coordinate `w=a+ib`. The operator actions are then:

| Real operator | Complex-coordinate action |
| --- | --- |
| `K_Gamma` | `w -> i w` |
| `B` | `w -> Re(w)` |
| `2B-I` | `w -> conjugate(w)` |
| `A` | `w -> i conjugate(w)` |
| `exp(t omega K_Gamma)` | `w -> exp(i omega t) w` |
| `A exp(t omega K_Gamma)` | `w -> i exp(-i omega t) conjugate(w)` |

The first four actions follow from the seam construction. Choosing the exponential dynamics with real rate `omega` requires the extra generator condition explained in the [R3 proof](../03-lambda-reference/SEAM_BOND_COMPLEX_STRUCTURE_R3.md).

The real operator algebra `span{I,K_Gamma}` is isomorphic to the complex numbers under **operator composition**. This does not make the source's componentwise product algebra `R x R` a field: that coordinate product still has zero divisors. The two multiplication laws are different, and no general four-way equivalence follows by renaming them.

If Aghora must be an isometry and the bond orthogonal, the compatible positive metric is proportional to `diag(Gamma^-2,1)`. In the seam coordinates its squared norm is proportional to `a^2+b^2`. That shape is derived under the stated requirements; its overall scale and the dynamical rate remain free.

Replacing the retained seam by its Aghora image reverses `K_Gamma`. Thus this construction supplies a complex action relative to a chosen retained seam, without selecting an absolute complex orientation, a physical metric, or a universal constant.
