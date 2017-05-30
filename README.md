
<!-- README.md is generated from README.Rmd. Please edit that file -->
mdpeer
======

**Graph-Constrained Regression with Enhanced Regularization Parameters**

Performs graph-constrained regularization in which regularization parameters are selected with the use of a known fact of equivalence between penalized regression and Linear Mixed Model solutions. Provides implementation of three regression methods where graph-constraints among coefficients are accounted for.

1.  `riPEERc` (ridgified Partially Empirical Eigenvectors for Regression with constant) method utilizes additional Ridge term to handle the non-invertibility of a graph Laplacian matrix.

2.  `vrPEER` (variable reducted PEER) method performs variable-reduction procedure to handle the non-invertibility of a graph Laplacian matrix.

3.  `riPEER` (ridgified Partially Empirical Eigenvectors for Regression) method employs a penalty term being a linear combination of graph-originated and ridge-originated penalty terms, whose two regularization parameters are ML estimators from corresponding Linear Mixed Model solution.

Notably, in `riPEER` method a graph-originated penalty term allows imposing similarity between coefficients based on graph information given whereas additional ridge-originated penalty term facilitates parameters estimation: it reduces computational issues arising from singularity in a graph- originated penalty matrix and yields plausible results in situations when graph information is not informative or when it is unclear whether connectivities represented by a graph reflect similarities among corresponding coefficients.
