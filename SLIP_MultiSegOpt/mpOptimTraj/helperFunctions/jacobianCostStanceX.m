function jacobianCostStanceX = jacobianCostStanceX(l,theta,ld,thetad,h)
%JACOBIANCOSTSTANCEX
%    JACOBIANCOSTSTANCEX = JACOBIANCOSTSTANCEX(L,THETA,LD,THETAD,H)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    03-Sep-2018 19:41:32

jacobianCostStanceX = [h.*l.*thetad.^2.*2.0,0.0,h.*ld.*2.0,h.*l.^2.*thetad.*2.0,0.0,0.0];
