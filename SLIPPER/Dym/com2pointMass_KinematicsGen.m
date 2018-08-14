function com2pointMass_KinematicsGen(normalizedMass, massScaling)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin==0
    normalizedMass=true;
    massScaling = 'bodyMass';
end
% COM motion to Point mass motion
% input
syms xc zc dxc dzc phi dphi
% output
% syms xf zf

% parameters
syms mf mb rc

[posFrame,posBody,velFrame,velBody] = pointMassMotionSynthesis(xc,zc, phi, dxc,dzc, dphi,mf,mb, rc);

if normalizedMass
posFrame = massNomralization(posFrame,massScaling);
posBody = massNomralization(posBody,massScaling);
velFrame = massNomralization(velFrame,massScaling);
velBody = massNomralization(velBody,massScaling);

end

%
pointMassMotionVector = [posFrame;posBody;velFrame;velBody];

if strcmp(massScaling, 'bodyMass')
variableVector = [xc zc phi dxc dzc dphi rc mf];
else
variableVector = [xc zc phi dxc dzc dphi rc mf mb];    
end

cd ./autoGen
matlabFunction(pointMassMotionVector,'File','motionCOM2PointMass','Vars',variableVector);
cd ../
end