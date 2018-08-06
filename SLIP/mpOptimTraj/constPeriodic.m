function c = constPeriodic(x,dx,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% State in Cartesian space at the start of the stance phase
xStance0 = [x(1,1),...
      dx(1,1),...
      x(2,1),...
      dx(2,1)];
  
    zS0 = xStance0( 1) * sin(xStance0( 3));
    zdS0 = xStance0( 2) * sin(xStance0( 3)) + xStance0( 1) * xStance0( 4) * cos(xStance0( 3));
    xS0 = -xStance0( 1) * cos(xStance0( 3)); %#ok<NASGU>
    xdS0 = -xStance0( 2) * cos(xStance0( 3)) + xStance0( 1) * xStance0( 4) * sin(xStance0( 3));

    velVecS = [xSd0, -ySd0];
    deltaOld = atan2(velVecS(2), velVecS(1)); 

%% State in Cartesian space at the end of the stance phase
xStanceEnd = [x(1,parms.phase(1).KnotNumber),...
      dx(1,parms.phase(1).KnotNumber),...
      x(2,parms.phase(1).KnotNumber),...
      dx(2,parms.phase(1).KnotNumber)];

    zF0 = xStanceEnd( 1) * sin(xStanceEnd( 3));
    zdF0 = xStanceEnd( 2) * sin(xStanceEnd( 3)) + xStanceEnd( 1) * xStanceEnd( 4) * cos(xStanceEnd( 3));
    xF0 = -xStanceEnd( 1) * cos(xStanceEnd( 3)); 
    xdF0 = -xStanceEnd( 2) * cos(xStanceEnd( 3)) + xStanceEnd( 1) * xStanceEnd( 4) * sin(xStanceEnd( 3));

%% State in Cartesian space at the end of the flight phase
xFlightEnd = [x(1,parms.totalKnotNumber),...
      dx(1,parms.totalKnotNumber),...
      x(2,parms.totalKnotNumber),...
      dx(2,parms.totalKnotNumber)];    
    
  
    velVecF = [xFlightEnd(2), -xFlightEnd(4)];
    deltaNew = atan2(velVecF(2), velVecF(1));   
%%

%transition Condition    
cTransition = [xStanceEnd(1)-1;
               xF0 - x(1,parms.phase(1).KnotNumber);...
               xdF0 - dx(1,parms.phase(1).KnotNumber);...
               zF0 - x(2,parms.phase(1).KnotNumber);...
               zdF0 - dx(2,parms.phase(1).KnotNumber)];
           
%boundary Condition           
cBoundary = [xStance0( 1)-1; %l0=1
             xStance0( 2)-parms.beta; %theta0=beta
             xdS0^2 + zdS0^2-1; %v0 = 1
             xFlightEnd(2)^2 + xFlightEnd(4)^2 -1; % vF = 1
             deltaNew - deltaOld; % deltaNew = deltaOld
               ]
    
c = [cTransition;cBoundary];
end

