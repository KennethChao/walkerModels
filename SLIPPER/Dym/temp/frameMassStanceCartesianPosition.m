function mfPosStance = frameMassStanceCartesianPosition(l,dl,theta,dtheta,phi,dphi,rc,mf)
%FRAMEMASSSTANCECARTESIANPOSITION
%    MFPOSSTANCE = FRAMEMASSSTANCECARTESIANPOSITION(L,DL,THETA,DTHETA,PHI,DPHI,RC,MF)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Aug-2018 10:17:01

mfPosStance = [-l.*cos(theta);l.*sin(theta)];
