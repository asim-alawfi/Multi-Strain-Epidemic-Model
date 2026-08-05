# Multi-Strain-Epidemic-Model

## Equilibria
The script [compute_equilibria_grobner.m](compute_equilibria_grobner.m) 
uses the symbolic toolbox to find **all** equilibria with a certain symmetry.
The variables $z$ and $w$ can be eliminated and the $y$ equations can be 
transformed into a system of polynomial equations. Systems of polynomial
equations can be transformed into a root-finding problem for a 
single-variable polynomial using the function `gbasis` 
(finding a Groebner basis).

The script [plot_equilibria.m](plot_equilibria.m) then plots a 
bifurcation diagram.