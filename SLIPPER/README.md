# SLIPPER - Fixed Point Finder for SLIP with PEndulum Runner (SLIPPER)
SLIPPER - Fixed Point Finder formulate the process of finding fixed-point at the Poincare section as a 
unconstrained nonlinear optimization to minimize the cost associated with the periodic conditions and control.

## Main Funtionalities
- __SLIPPER_fixedPointFinderMain-__ Main script of fixed point finder where the most of the parameters for fixed-point finder are speicifed here. The search result will stored in the 'result' folder.
- __showFixedPointResultSLIPPER-__ The script to plot fixed points for showing the SLIP and SLIPPER comparison.
- __showPhasePortraitSLIPPER-__ The script to plot phase portrait of the set of stable fixed points.
- __showAnimationSLIPPER-__ The script to show animation of the selected stable fixed point.

## Dependencies
It requires 'Optimization Toolbox' and 'Parallel Computing Toolbox'(optional).

### Others
## In 'dynamics' folder
- dynamicsGEN.m is the main function to generate the kinematic/dynamic function calculated using symbolic tool.
- runUNittests.m contains the test cases for functions of the kinematics and the dynamics which are used to generate the kinematic/dynamic function of SLIPPER model.

## Others details
