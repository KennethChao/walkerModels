function expr =  kinematicsGEN(write2Function,massScaling)
%KINEMATICS Summary of this function goes here
%   Detailed explanation goes here
if nargin==0  
    write2Function=true;
    massScaling = 'bodyMass';  
end
%%
syms l ld theta thetad phi phid
syms ldd thetadd phidd
syms rc
%% Forward Kinematics

% motion of mass of frame
zf = l*sin(theta);
xf = -l*cos(theta);

xVec = [l;theta;phi];
xVecd = [ld;thetad;phid];
xVecdd = [ldd;thetadd;phidd];

zfd = timeDiff(zf,xVec,xVecd,xVecdd);
xfd = timeDiff(xf,xVec,xVecd,xVecdd);

% motion of body mass

% zb = zf -rc*sin(phi);
% xb = xf +rc*cos(phi) ;

zb = zf -rc*cos(phi);
xb = xf -rc*sin(phi) ;

zbd = timeDiff(zb,xVec,xVecd,xVecdd);
xbd = timeDiff(xb,xVec,xVecd,xVecdd);

% 
expr.zf = zf;
expr.xf = xf;
expr.zfd = zfd;
expr.xfd = xfd;
expr.zb = zb;
expr.xb = xb;
expr.zbd = zbd;
expr.xbd = xbd;

%%
if write2Function

syms l ld theta thetad phi phid

mfPosStance = [xf;zf];
mfVelStance = [xfd;zfd];
mfMotion = [mfPosStance;mfVelStance];

mbPosStance = [xb;zb];
mbVelStance = [xbd;zbd];
mbMotion = [mbPosStance;mbVelStance];


%%
cd ./autoGen
variableVector = [l, theta, phi, ld, thetad, phid, rc];
matlabFunction(mfMotion,'File','frameMassCartesianMotionStance','Vars',variableVector);

matlabFunction(mbMotion,'File','bodyMassCartesianMotionStance','Vars',variableVector);
cd ../

%
if strcmp(massScaling,'bodyMass') || strcmp(massScaling,'totalMass')
pointMass2COM_KinematicsGen(true, massScaling);
com2pointMass_KinematicsGen(true, massScaling);
end

end