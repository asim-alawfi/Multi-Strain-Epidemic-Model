# Multi-Strain-Epidemic-Model - periodic orbit computations

## Right-hand side generation
The script [gen_dengue_fever_sym.m](gen_dengue_fever_sym.m) creates the 
right-hand side for computations (needs to be run first) in file 
[dengue_fever_sym.m](dengue_fever_sym.m).

## Main script
Run [POs_dengue.m](POs_dengue.m) to perform periodic orbit computations. 
It automatically loads the DDE-Biftool path.

**Warning:** This script uses a modified small parameter $e=0.019$ (instead of $e=0.0019$).
The variable $y$ is also scaled by $1/e$.

## Dependence on steady-state computations
At the beginning the file loads [../stst/stst_sym_e_0_019.mat](../stst/stst_sym_e_0_019.mat),
so one needs to run [../stst/compute_equilibria_grobner.m](../stst/compute_equilibria_grobner.m)
first.

1. The script detects the Hopf bifurcation at $r_H$ on the $(13)(42)$-symmetric 
   equilibrium branch, and branches off from there. The periodic orbits are 
   unstable and end in a connecting orbit to the $(1234)$-symmetric equilibrium.

2. For $r>r_H$ there is no stable equilibrium. We run a simulation that 
   settles in a periodic orbit that is $(13)(24)$ symmetric with half-period
   phase shift.

3. From there we initialise a branch of periodic orbits, which is possible 
   to compute for the larger parameter $e$. At the low-$r$ end the parameter
   also ends in a connecting orbit to the $(1234)$-symmetric equilibrium. 
   At the large-$r$ end the periodic orbits look similar to those for smaller
   parameter $e$.

The resulting bifurcation diagram is [bifdiagram_r.html](bifdiagram_r.html).