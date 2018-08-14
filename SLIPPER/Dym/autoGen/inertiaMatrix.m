function Mmat = inertiaMatrix(l,theta,phi,ld,thetad,phid,g,k,mf,rc)
%INERTIAMATRIX
%    MMAT = INERTIAMATRIX(L,THETA,PHI,LD,THETAD,PHID,G,K,MF,RC)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Aug-2018 18:43:32

t2 = mf+1.0;
t3 = phi-theta;
t4 = cos(t3);
t5 = rc.*t4;
t6 = sin(t3);
t7 = l.*rc.*t6;
Mmat = reshape([t2,0.0,t5,0.0,l.^2.*t2,t7,t5,t7,rc.^2],[3,3]);
