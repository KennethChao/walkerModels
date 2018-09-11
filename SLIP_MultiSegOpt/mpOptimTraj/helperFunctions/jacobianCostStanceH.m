function jacobianCostStanceH = jacobianCostStanceH(l,theta,ld,thetad,h,g)
%JACOBIANCOSTSTANCEH
%    JACOBIANCOSTSTANCEH = JACOBIANCOSTSTANCEH(L,THETA,LD,THETAD,H,G)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Sep-2018 22:35:00

t2 = sin(ld);
t4 = cos(ld);
t3 = t2.*theta+l.*t4.*thetad;
jacobianCostStanceH = g.*(t4.*theta-l.*t2.*thetad)+t3.^2+l.^2.*t2.^2;
