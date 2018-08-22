# SLIPPER (SLIP with PEndulum Runner) Fixed-Point Finder
SLIPPER Fixed-Point Finder formulates the process of finding fixed-point at the Poincare section as a 
unconstrained nonlinear optimization to minimize the cost associated with the periodic conditions and control.

## Main funtionalities
- __SLIPPER_fixedPointFinderMain.m__: Main script of fixed point finder where the most of the parameters are speicifed here. The search result will stored in the 'result' folder.
- __showFixedPointResultSLIPPER.m__: The script to show how to plot fixed points for showing the SLIP and SLIPPER comparison using helper functions.
- __showPhasePortraitSLIPPER.m__: The script to show how to  plot phase portrait of the set of stable fixed points using helper functions.
- __showAnimationSLIPPER.m__: The script to show animation of the selected stable fixed point, and generate the animated gif file.
- Note: All the kinematic and dynamic variables (except variables related to angle) used in SLIPPER Fixed-Point Finder are dimensionless, i.e. the legnth is scaled by the unstretched leg length 'l_0', and the mass is scaled by body mass. For more information please refer to:
Reference:
   Shen Z, Seipel J. A Piecewise-Linear Approximation of the Canonical 
   Spring-Loaded Inverted Pendulum Model of Legged Locomotion. ASME. J. 
   Comput. Nonlinear Dynam., 2016.

## Dependencies
It requires 'Optimization Toolbox' and 'Parallel Computing Toolbox' (optional).

## In the 'dynamics' folder
- dynamicsGEN.m is the main function to generate the kinematic/dynamic function calculated using symbolic tool.
- runUnittests.m contains the test cases for functions of the kinematic and the dynamic equation derivations, which are used to generate the kinematic/dynamic function of SLIPPER model.
- Note: for the kinematics and dynamics used in SLIPPER Fixed-Point Finder, the legnth is scaled by the leg length 'l0', and the mass is scaled by body mass. (The event function 'eventFncTouchDownSLIPPER' is also using the same scaling). To change the mass scaling, both 'dynamicsGEN.m' and 'eventFncTouchDownSLIPPER.m' need to be modified.

## In the 'result/SLIPPER_PControl' folder
The folder conatins the fixed-point data for SLIP and SLIPPER model comparison.

## Others details