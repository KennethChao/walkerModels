function jacobianStance = jacobianControlledStanceDym(l,theta,ld,thetad,u,g,k)
%JACOBIANCONTROLLEDSTANCEDYM
%    JACOBIANSTANCE = JACOBIANCONTROLLEDSTANCEDYM(L,THETA,LD,THETAD,U,G,K)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    05-Sep-2018 19:47:22

t2 = cos(theta);
t3 = t2.*u;
t4 = sin(theta);
t5 = 1.0./l;
jacobianStance = reshape([k-thetad.^2,-1.0./l.^2.*(t3+ld.*thetad.*2.0),t3,-t4.*t5.*u,0.0,t5.*thetad.*2.0,l.*thetad.*-2.0,ld.*t5.*2.0,1.0,0.0,0.0,1.0,t4,t2.*t5],[2,7]);
