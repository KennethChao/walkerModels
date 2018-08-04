function [c, ceq] = periodicGait(xF,x0,parms)
%[c, ceq] = periodicGait(zBefore,zAfter,p)
%
% Puts a periodic constraint on the gait, including the heel-strike and
% foot switch.
%
% INPUTS:
%   zBefore = [4,1] = state before heel-strike
%   zAfter = [4,1] = state after heel-strike
%   p = parameter struct:
%       .m1 = hip mass
%       .m2 = foot mass
%       .g = gravitational acceleration
%       .l = leg length
%
% OUTPUTS:
%   c = [];
%   ceq = [4, 1] = defect constraint on periodic step map
% 
% NOTES:
%   
%   states:
%       1 = q1 = first link angle
%       2 = q2 = second link angle
%       3 = dq1 = first link angular rate
%       4 = dq2 = second link angular rate
%
%   angles: measured from negative j axis with positive convention
%

[velVec,deltaNew, deltaOld,costFly] = dymFlightDimensionless(x0,xF,1,parms);


c = [];
% ceq = [x0(1)-1;
%         xF(1)-1;
%        x0(3)-parms.beta;
%     norm(velVec)-1;
%       deltaNew-deltaOld];

ceq = [xF(1)-1;
       norm(velVec)-1;
      deltaNew-deltaOld];
  
  
end