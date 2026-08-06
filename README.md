# Multi-Strain-Epidemic-Model

Model equations from 

[B14] K Blyuss (2014) *Analysis of symmetries in models of multi-strain
infections*. J. Math. Biol. (2014) 69:1431–1459,
DOI: [10.1007/s00285-013-0740-0](https://doi.org/10.1007/s00285-013-0740-0).

In contrast to [B14] we rescale the variable $y$ by thr parameter $e$, so now study
$$
\dot y_i  =  ry_i(1-\gamma w_i-(1-\gamma)z_i)-y_i
$$
$$
\dot z_i  =  e (ry_i(1-z_i)-z_i)
$$
$$
\dot w_i  =  e ((1-w_i)r[Ay]_i-w_i),
$$
where $i=1,\ldots,4$, $e$ is small, and the coupling matrix $A$, given by

    A =
     1     1     0     1
     1     1     1     0
     0     1     1     1
     1     0     1     1


