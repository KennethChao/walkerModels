function costStance = costStance(l,theta,ld,thetad,h,g,k)
%COSTSTANCE
%    COSTSTANCE = COSTSTANCE(L,THETA,LD,THETAD,H,G,K)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Sep-2018 11:40:03

t3 = cos(theta);
t4 = sin(theta);
t2 = ld.*t4+l.*t3.*thetad;
t5 = ld.*t3-l.*t4.*thetad;
t6 = l-1.0;
costStance = -h.*(k.*t6.^2.*(1.0./2.0)-t2.^2.*(1.0./2.0)-t5.^2.*(1.0./2.0)+g.*l.*t4);
