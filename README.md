# Stencils2D
Stencils for 2D differential operators. Mathematica file to supplement the article *Stencils with isotropic discretisation error for differential operators*, Michael Patra and Mikko Karttunen, [Numerical Methods for Partial Differential Equations 22, 936-953 (2005)]( https://doi.org/10.1002/num.20129). 

In numerical solutions of physical/mathematical problems, one often has to discretize the Laplacian and bi-Laplacian operators. There are, however, several ways of choosing them even in two dimensions For example, stencil can be isotropic or anisotropic. The picture shows error propagation induced by a particular stencil. It is important to choose a stencil with some care since error can propagate anisotropically and lead to creation of structures and patterns that are not due to the underlying physics but rather due to bad numerical solution of the problem. Such issues are important, for example in simulations of Turing patterns.

[Error propagation](error-propagation.png)
