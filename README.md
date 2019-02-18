# multi-start-4D
Creates molecular movies of 4-atom molecules from Coulomb explosion imaging experiments.

This was an extension of the triatomic molecular geometry reconstruction to 4-atom molecules. It also uses MATLAB's MultiStart global optimization solver to locate multiple local solutions in hopes of finding the global solution. I tested it on acetylene (C2H2) and it was able to reconstruct some geometries but convergence was difficult. In hindsight this may have been because I used a logarithmic objective function (which comes with a singularity). A logarithmic objective function worked for triatomic molecules (3-dimensional parameter space) but might work badly in higher dimensional spaces (i.e. for 4+ atom molecules).
