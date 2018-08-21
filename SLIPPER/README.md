# SLIPPER (SLIP with PEndulum Runner) Fixed-Point Finder
SLIPPER Fixed-Point Finder formulates the process of finding fixed-point at the Poincare section as a 
unconstrained nonlinear optimization to minimize the cost associated with the periodic conditions and control.

## Main Funtionalities
- __SLIPPER_fixedPointFinderMain.m__: Main script of fixed point finder where the most of the parameters for fixed-point finder are speicifed here. The search result will stored in the 'result' folder.
- __showFixedPointResultSLIPPER.m__: The script to plot fixed points for showing the SLIP and SLIPPER comparison.
- __showPhasePortraitSLIPPER.m__: The script to plot phase portrait of the set of stable fixed points.
- __showAnimationSLIPPER.m__: The script to show animation of the selected stable fixed point, and generate the animated gif file.

## Dependencies
It requires 'Optimization Toolbox' and 'Parallel Computing Toolbox' (optional).

## In the 'dynamics' folder
- dynamicsGEN.m is the main function to generate the kinematic/dynamic function calculated using symbolic tool.
- runUnittests.m contains the test cases for functions of the kinematic and the dynamic equation derivations, which are used to generate the kinematic/dynamic function of SLIPPER model.

## In the 'result/SLIPPER_PControl' folder
- the fixed-point data for SLIP and SLIPPER model comparison.

## Others details
