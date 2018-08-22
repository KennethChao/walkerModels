# SLIP (Spring Loaded Inverted Pendulum) Fixed-Point Finder
SLIP Fixed-Point Finder formulates the process of finding fixed-point at the Poincare section as a 
unconstrained nonlinear optimization to minimize the cost associated with the periodic conditions.

## Main funtionalities
- __SLIP_fixedPointFinderMain.m__: Main script of fixed point finder where the most of the parameters for fixed-point finder are speicifed here. The search result will stored in the 'result' folder.
- __showAllFixedPointsSLIP.m__: The script to show how to use helper functions to plot fixed points of SLIP model.
- Note: All the kinematic and dynamic variables (except variables related to angle) used in SLIP Fixed-Point Finder are dimensionless, i.e. the legnth is scaled by the leg length 'l0', and the mass is scaled by body mass. For more information please refer to:

## Dependencies
It requires 'Optimization Toolbox' and 'Parallel Computing Toolbox' (optional).

## In the 'result/Comparison' folder
- the fixed-point data for SLIP and SLIPPER model comparison.

## Others details
