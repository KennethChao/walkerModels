function c = boundaryCostDummy(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%BOUNDARYCOSTDUMMY
%    C = BOUNDARYCOSTDUMMY(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    09-Sep-2018 13:35:15

t2 = l0d.^2+l0.^2.*theta0d.^2-1.0;
t4 = sin(theta0);
t5 = sin(thetaEnd);
t6 = cos(theta0);
t7 = cos(thetaEnd);
t3 = l0d.*t6-lEndd.*t7.*5.0-l0.*t4.*theta0d+lEnd.*t5.*thetaEndd.*5.0;
t8 = l0d.*t4+lEndd.*t5.*5.0+l0.*t6.*theta0d+lEnd.*t7.*thetaEndd.*5.0;
c = [t2.^2,t3.^2,t8.^2];