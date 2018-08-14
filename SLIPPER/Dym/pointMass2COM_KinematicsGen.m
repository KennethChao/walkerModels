function pointMass2COM_KinematicsGen(massNormalization, massScaling)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin==0
    massNormalization=true;
    massScaling = 'bodyMass';
end

% Point mass motion to COM motion
syms xf zf xb zb xfd zfd xbd zbd

syms mf mb rc

syms xc zc

xc = comMotionSynthesis(xf,xb,mf,mb);
zc = comMotionSynthesis(zf,zb,mf,mb);

xcd = comMotionSynthesis(xfd,xbd,mf,mb);
zcd = comMotionSynthesis(zfd,zbd,mf,mb);

if massNormalization
xc = massNomralization(xc,massScaling);
zc = massNomralization(zc,massScaling);

xcd = massNomralization(xcd,massScaling);
zcd = massNomralization(zcd,massScaling);
end

%
comMotionVector = [xc;zc;xcd;zcd];
if strcmp(massScaling, 'bodyMass')
variableVector = [xf zf xb zb xfd zfd xbd zbd rc mf];
else
variableVector = [xf zf xb zb xfd zfd xbd zbd rc mf mb];    
end

cd ./autoGen
matlabFunction(comMotionVector,'File','motionPointMass2COM','Vars',variableVector);
cd ../

end

