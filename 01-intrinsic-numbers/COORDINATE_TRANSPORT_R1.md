# R1: A coordinate realization of the response network

Choose a scalar field `F` and a labelled direct sum `W=F^n`, one coordinate line for each vertex of the finite network. This gives a precise Way-1 representation for the [R1 return model](../03-lambda-reference/RETURN_INVARIANTS_R1.md).

The projector `P_i` selects the line at vertex `i`. A directed edge `i->j` acts through `T_ji=a_ij E_ji`. Thus a path acts on its starting coordinate and deposits the transported value in its ending coordinate. A closed path acts as `Lambda_C P_i` on the full direct sum, and as multiplication by `Lambda_C` on its starting line.

The reverse edge is inverse on the appropriate source and target lines; its product with the forward edge gives `P_i` on `W`, not the identity on every coordinate. Keeping these types explicit repairs a common ambiguity in an untyped product of operators.

Independent coordinate rescalings give a diagonal map `D_g`; edge operators transform by `D_g T_ji D_g^{-1}`. A labelled tree-normalized network and its return coefficients describe the reference-equivalence class. Adding the tree transport factors reconstructs the original coordinate representative.

This construction uses a specified field and finite vector spaces. It does not depend on the source's incomplete valued-field definition or prove a new field construction. The original coordinate algebra remains available; the added labelled transport data supplies the response and composition information that a bare coordinate vector would not contain.
