function mbMotion = bodyMassCartesianMotionStance(l,theta,phi,ld,thetad,phid,rc)
%BODYMASSCARTESIANMOTIONSTANCE
%    MBMOTION = BODYMASSCARTESIANMOTIONSTANCE(L,THETA,PHI,LD,THETAD,PHID,RC)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Aug-2018 12:35:47

t2 = cos(theta);
t3 = cos(phi);
t4 = sin(theta);
t5 = sin(phi);
mbMotion = [-l.*t2-rc.*t5;l.*t4-rc.*t3;-ld.*t2+l.*t4.*thetad-phid.*rc.*t3;ld.*t4+l.*t2.*thetad+phid.*rc.*t5];
