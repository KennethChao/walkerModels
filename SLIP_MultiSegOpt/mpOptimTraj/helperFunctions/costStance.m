function costStance = costStance(l,theta,ld,thetad,h)
%COSTSTANCE
%    COSTSTANCE = COSTSTANCE(L,THETA,LD,THETAD,H)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    03-Sep-2018 19:41:32

costStance = h.*(ld.^2+l.^2.*thetad.^2);
