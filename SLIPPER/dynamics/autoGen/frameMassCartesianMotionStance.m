function mfMotion = frameMassCartesianMotionStance(l,theta,phi,ld,thetad,phid,rc)
%FRAMEMASSCARTESIANMOTIONSTANCE
%    MFMOTION = FRAMEMASSCARTESIANMOTIONSTANCE(L,THETA,PHI,LD,THETAD,PHID,RC)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Aug-2018 09:47:18

t2 = cos(theta);
t3 = sin(theta);
mfMotion = [-l.*t2;l.*t3;-ld.*t2+l.*t3.*thetad;ld.*t3+l.*t2.*thetad];