# Weight-function-on-WeylGroup
Quantum Bruhat graph arose in the context of enumerative geometry; in particular, it was first introduced in [BFP'99] to describe the multiplicative structure of the quantum cohomology ring of the complex flag variety, in particular the Chevalley-Monk rule for multiplying by a divisor class. In recent years, it has proved to be useful in the study many Lie theoretic and arithmetic geometric problems. To name a few instances, it encodes the covering relation in extended affine Weyl group $\widetilde{W}$- as first shown by [M21] and later in [S21+]; it also features in the formula for Demazure product of elements in $\widetilde{W}$, cf. [HN21+]. It has also been utilized to give explicit description of admissible sets in $\widetilde{W}$ and generic Newton point in Iwahori double cosets in loop group of a reductive algebraic group, cf. [HY21], [S21+]. 

The quantum Bruhat graph is a directed graph that is associated to a root system $\Phi$. Its vertices are labelled by the elements of the Weyl group generated by $\Phi$ and the upward edges encode the usual Bruhat graph relations while the downward edges occur between two elements if they differ by multiplying a quantum reflection that decreases the length as much as possible; the downward edges are then labelled by the coroot associated with the reflection, called weight of the edge. This gives rise to a function on $W \times W$: at $(x,y)$ it outputs the sum of weights of edges appearing in a(ny) shortest path from $x$ to $y$. This function is the main novel feature that appears in the contexts mentioned above, and we implement the computation of this function by utilizing a greedy algorithm.

References:
[BFP'99]: F. Brenti, S. Fomin, and A. Postnikov, Mixed Bruhat operators and Yang– Baxter equations for Weyl groups, Int. Math. Res. Not. 8 (1999), 419–441.

[M21]: E. Mili´cevi´c, Maximal Newton points and the quantum Bruhat graph, Michigan Math. J. Advance Publication (2021), 1-52.

[S21+]: A. Sadhukhan, Affine Deligne-Lusztig variety and quantum Bruhat graph, https://arxiv.org/pdf/2110.02172.pdf.

[HN21+]: X. He, S. Nie, Demazure product of the affine Weyl groups, https://arxiv.org/pdf/2112.06376.pdf.

[HY21]:  X. He and Q. Yu, Dimension formula for the affine Deligne–Lusztig variety X(µ, b), Mathematische Annalen 379 (2021), 1747–1765.

