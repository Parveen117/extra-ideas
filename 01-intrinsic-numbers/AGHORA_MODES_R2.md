# R2: The Aghora relations require two sectors for nonzero evolution

Over real or complex scalars, `A^2=I` splits the state space into the positive and negative eigenspaces of `A`. The relation `AG=-GA` forces the generator to send each sector into the other. In a basis respecting this split,

`A=diag(I_+,-I_-)`, and `G=[[0,D],[E,0]]`.

If only one sector exists, `A=I` or `A=-I`, and the source relation forces `G=0`. Thus a nonzero generator requires at least two modes. A single scalar transport coefficient cannot encode this internal structure.

For the constructed return `R_t=A exp(tG)`, the sector projectors are `(I+R_t)/2` and `(I-R_t)/2`. They are algebraic projectors and are related to the original Aghora projectors by the same similarity that relates `R_t` to `A`.

Choosing one projector as a bond would select its return sign. Both choices are compatible with the source's algebra; deciding which subspace the native seam retains remains an additional problem. A general bond can mix the two sectors and produce a projected coefficient different from either full-return sign.

See [AGHORA_RETURN_R2.md](../03-lambda-reference/AGHORA_RETURN_R2.md) for the proofs, exact examples, and the conditions needed for a norm interpretation.
