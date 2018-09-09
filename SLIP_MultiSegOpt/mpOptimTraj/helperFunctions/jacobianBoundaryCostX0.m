function jacobianStance0 = jacobianBoundaryCostX0(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%JACOBIANBOUNDARYCOSTX0
%    JACOBIANSTANCE0 = JACOBIANBOUNDARYCOSTX0(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    09-Sep-2018 13:35:14

t2 = theta0d.^2;
t3 = l0d.^2;
t4 = l0.^2;
t5 = t2.*t4;
t6 = t3+t5-1.0;
t7 = sin(theta0);
t8 = cos(theta0);
t9 = l0d.*t8;
t10 = cos(thetaEnd);
t11 = sin(thetaEnd);
t12 = lEnd.*t11.*thetaEndd;
t14 = lEndd.*t10;
t15 = l0.*t7.*theta0d;
t13 = t9+t12-t14-t15;
t16 = l0d.*t7;
t17 = l0.*t8.*theta0d;
t18 = lEndd.^2;
t19 = t11.^2;
t20 = t18.*t19;
t21 = sin(beta);
t22 = g.*t21.*2.0;
t23 = lEnd.^2;
t24 = thetaEndd.^2;
t25 = t10.^2;
t26 = t23.*t24.*t25;
t27 = lEnd.*lEndd.*t10.*t11.*thetaEndd.*2.0;
t32 = g.*lEnd.*t11.*2.0;
t28 = t20+t22+t26+t27-t32;
t29 = sqrt(t28);
t30 = conj(t29);
t31 = t16+t17-t30;
jacobianStance0 = reshape([l0.*t2.*t6.*4.0,t7.*t13.*theta0d.*-2.0,t8.*t31.*theta0d.*2.0,0.0,t13.*(t16+t17).*-2.0,t31.*(t9-t15).*2.0,l0d.*t6.*4.0,t8.*t13.*2.0,t7.*t31.*2.0,t4.*t6.*theta0d.*4.0,l0.*t7.*t13.*-2.0,l0.*t8.*t31.*2.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
