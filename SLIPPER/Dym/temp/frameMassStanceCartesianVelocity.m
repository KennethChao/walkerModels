function mfVelStance = frameMassStanceCartesianVelocity(l,dl,theta,dtheta,phi,dphi,rc,mf)
%FRAMEMASSSTANCECARTESIANVELOCITY
%    MFVELSTANCE = FRAMEMASSSTANCECARTESIANVELOCITY(L,DL,THETA,DTHETA,PHI,DPHI,RC,MF)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Aug-2018 10:17:01

t2 = sin(theta);
t3 = cos(theta);
mfVelStance = [-dl.*t3+dtheta.*l.*t2;dl.*t2+dtheta.*l.*t3];
