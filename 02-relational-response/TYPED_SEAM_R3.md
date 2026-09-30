# R3: A seam constraint and a projected response

The source's seam ratio is a scalar relation. Before inserting it into a state-operator product, specify its type. In R3's chosen real two-mode realization:

$$
q_\Gamma:\mathbb R^2\to\mathbb R,\qquad q_\Gamma(x,y)=x-\Gamma y,
\qquad L_\Gamma=\ker q_\Gamma.
$$

This linear constraint is an added realization of the ratio idea. Its coordinates have not been identified with the source's physical quantities. It is also defined at the zero state, where `x/y` is not.

The constraint fixes a retained line but leaves many projections onto it. Add the explicit rule that the removed line is its Aghora image, using `A=diag(1,-1)` and finite real `Gamma != 0`. The unique bond is

$$
B_\Gamma=\frac12\begin{pmatrix}1&\Gamma\\\Gamma^{-1}&1\end{pmatrix},
\qquad q_\Gamma B_\Gamma=0,\qquad AB_\Gamma A=I-B_\Gamma.
$$

Now `B_Gamma:V -> V` is a composable state operator, while `q_Gamma:V -> R` tests the relation. They are different maps; no scalar limit is silently treated as an endomorphism.

Under the additional metric and generator conditions in the [R3 proof](../03-lambda-reference/SEAM_BOND_COMPLEX_STRUCTURE_R3.md), `G=omega[A,B_Gamma]` and the exponential return `R_t=A exp(tG)` satisfies

$$
B_\Gamma R_tB_\Gamma=\sin(\omega t)B_\Gamma,
\qquad B_\Gamma-(B_\Gamma R_tB_\Gamma)^2=\cos^2(\omega t)B_\Gamma.
$$

This is an observable compressed response for a specified rate and time. Both pure return signs are possible when there is no leakage. The coefficient is not a universal fixed reference. Squared-norm fractions follow only from the stated compatibility conditions; a counterexample in the proof has zero scalar defect but nonzero leakage when these conditions fail.

For the exact Cayley family, the response is `2z/(1+z^2)` with `z=omega*s`. The chosen example `Gamma=2, omega=3, s=1/9` gives coefficient `3/5` and removed squared-norm fraction `16/25`. The source's full composition and the law fixing `omega` remain open.
