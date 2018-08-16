function bvec = nonInertiaTerms(l,theta,phi,ld,thetad,phid,g,k,mf,rc,I)
%NONINERTIATERMS
%    BVEC = NONINERTIATERMS(L,THETA,PHI,LD,THETAD,PHID,G,K,MF,RC,I)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Aug-2018 09:47:20

t2 = thetad.^2;
t3 = sin(theta);
t4 = cos(theta);
t5 = phid.^2;
t6 = phi-theta;
t7 = sin(t6);
t8 = cos(t6);
bvec = [k-k.*l-g.*t3+l.*t2-g.*mf.*t3+l.*mf.*t2+rc.*t5.*t7;-l.*(g.*t4+ld.*thetad.*2.0+g.*mf.*t4+ld.*mf.*thetad.*2.0+rc.*t5.*t8);-rc.*(g.*sin(phi)-l.*t2.*t8+ld.*t7.*thetad.*2.0)];
