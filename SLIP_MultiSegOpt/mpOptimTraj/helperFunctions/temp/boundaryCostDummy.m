function c = boundaryCostDummy(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%BOUNDARYCOSTDUMMY
%    C = BOUNDARYCOSTDUMMY(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Sep-2018 11:17:33

t3 = sin(theta0);
t4 = sin(thetaEnd);
t5 = cos(theta0);
t6 = cos(thetaEnd);
t2 = l0d.*t5-lEndd.*t6.*5.0-l0.*t3.*theta0d+lEnd.*t4.*thetaEndd.*5.0;
t7 = l0d.*t3+lEndd.*t4.*5.0+l0.*t5.*theta0d+lEnd.*t6.*thetaEndd.*5.0;
c = [t2.^2,t7.^2];