function [posFrame,posBody,velFrame,velBody] = pointMassMotionSynthesis(xc,zc, phi, xcd,zcd, phid,massFrame,massBody, rc)
%POINTMASSMOTIONSYNTHESIS Summary of this function goes here
%   Detailed explanation goes here
lengthC2f = rc*massBody/(massFrame + massBody);
lengthC2b = rc*massFrame/(massFrame + massBody);

xf = xc + lengthC2f*sin(phi);
zf = zc + lengthC2f*cos(phi);

xb = xc - lengthC2b*sin(phi);
zb = zc - lengthC2b*cos(phi);

vecC2f = [ xf-xc, 0,zf-zc ];
vecC2b = [ xb-xc, 0,zb-zc ];

omega = [0,phid,0];

velCOM = [xcd,0,zcd];

velF = velCOM + cross(omega,vecC2f);
velB = velCOM + cross(omega,vecC2b);


posFrame = [xf;zf];
posBody = [xb;zb];

velFrame = [velF(1);velF(3)];
velBody = [velB(1);velB(3)];
end

