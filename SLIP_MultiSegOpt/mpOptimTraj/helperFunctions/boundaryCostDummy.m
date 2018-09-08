function c = boundaryCostDummy(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%BOUNDARYCOSTDUMMY
%    C = BOUNDARYCOSTDUMMY(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    08-Sep-2018 15:20:41

t2 = lEnd-1.0;
t3 = l0d.^2+l0.^2.*theta0d.^2-1.0;
t5 = sin(theta0);
t6 = sin(thetaEnd);
t7 = cos(theta0);
t8 = cos(thetaEnd);
t4 = l0d.*t7-lEndd.*t8.*5.0-l0.*t5.*theta0d+lEnd.*t6.*thetaEndd.*5.0;
t9 = l0d.*t5+lEndd.*t6.*5.0+l0.*t7.*theta0d+lEnd.*t8.*thetaEndd.*5.0;
c = [t2.^2,t3.^2,t4.^2,t9.^2];
