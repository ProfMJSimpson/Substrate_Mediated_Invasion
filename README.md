# Substrate_Mediated_Invasion
Julia and Matlab codes to simulated all problems in El-Hachem, McCue and Simpson (2021)

2DSolver.jl reproduces the simulation results in Figure 1 by solving (1)--(2) on a square domain.  The code plots both u(x,y,t) and s(x,y,t) 
as in Figure 1, and also produces plots of u(x,L/2,t) and s(x,L/2,t) to show the cross-section through the evolving density profiles.  Use lines 36--40 
to alter the parameter values.  Default values reproduce results in Figure 1.

Subsstrate1D_PDE_solve.m is a MATLAB code that solves (3)--(5). Parameters are defined on lines 20-21.

SlowManifold_ODE_solver.m is a MATLAB code that solves (39)--(40). Parameters are defined on lines 15-17.
