function kinematicData = extractKinematicDatafromResult(result, parms)
%EXTRACTKINEMATICDATAFROMRESULT Extract kinematic data from the result
%   The original results contains states in different cooridnates, this
%   function is used to generate vector with consistent coordinate as a
%   helper function for 'plotPhasePortraitSLIPPER', or 'showAnimationSLIPPER'
%

% Get original state variable and time data
x = result.x;
t = result.t;
x2 = result.x2;
t2 = result.t2;
rc = parms.rc;
mf = parms.mf;
beta = parms.beta;

%% Extract kinematic data

% Kinematic data in stance phase
l = x(:, 1)';
ld = x(:, 2)';
theta = x(:, 3)';
thetad = x(:, 4)';
phi = x(:, 5)';
phid = x(:, 6)';

mfMotion = frameMassCartesianMotionStance(l, theta, phi, ld, thetad, phid, rc);
xf = mfMotion(1, :);
zf = mfMotion(2, :);
xfd = mfMotion(3, :);
zfd = mfMotion(4, :);

mbMotion = bodyMassCartesianMotionStance(l, theta, phi, ld, thetad, phid, rc);
xb = mbMotion(1, :);
zb = mbMotion(2, :);
xbd = mbMotion(3, :);
zbd = mbMotion(4, :);

comMotion = motionPointMass2COM(xf, zf, xb, zb, xfd, zfd, xbd, zbd, rc, mf);

% Kinematic data in flight phase
zc = x2(:, 1)';
zcd = x2(:, 2)';
phi = x2(:, 3)';
phid = x2(:, 4)';
xcd = ones(size(zc)) * comMotion(3, end);
xc = comMotion(1, end) + t2' .* xcd;
pointMassMotion = motionCOM2PointMass(xc, zc, phi, xcd, zcd, phid, rc, mf);

%% Combine Stance and Flight kinematic data

% frame mass motion
xfVecStance = xf;
xfVecFlight = pointMassMotion(1, 2:end);

zfVecStance = zf;
zfVecFlight = pointMassMotion(2, 2:end);

xfVec = [xfVecStance, xfVecFlight];
xfdVec = [xfd, pointMassMotion(5, 2:end)];
zfVec = [zfVecStance, zfVecFlight];
zfdVec = [zfd, pointMassMotion(6, 2:end)];

% body mass motion
xbVecStance = xb;
xbVecFlight = pointMassMotion(3, 2:end);

zbVecStance = zb;
zbVecFlight = pointMassMotion(4, 2:end);

xbVec = [xbVecStance, xbVecFlight];
xbdVec = [xbd, pointMassMotion(7, 2:end)];

zbVec = [zbVecStance, zbVecFlight];
zbdVec = [zbd, pointMassMotion(8, 2:end)];

% COM motion
xcVec = [comMotion(1, :), xc(2:end)];
xcdVec = [comMotion(3, :), xcd(2:end)];
zcVec = [comMotion(2, :), zc(2:end)];
zcdVec = [comMotion(4, :), zcd(2:end)];

% Pendulum motion
phiVec = [x(:, 5); x2(2:end, 3)]';
phidVec = [x(:, 6); x2(2:end, 4)]';

% Foot position
zfootVec = [zeros(size(zfVecStance)), zfVecFlight - 1 * sin(beta)];
xfootVec = [zeros(size(xfVecStance)), xfVecFlight + 1 * cos(beta)];

% Simulation time
tVec = [t; t2(2:end) + t(end)]';

%% Store the kienmatic data into the struct

kinematicData.tVec = tVec;

% Pendulum motion
kinematicData.phiVec = phiVec;
kinematicData.phidVec = phidVec;

% Frame mass motion
kinematicData.zfVec = zfVec;
kinematicData.zfdVec = zfdVec;
kinematicData.xfVec = xfVec;
kinematicData.xfdVec = xfdVec;

% Body mass motion
kinematicData.zbVec = zbVec;
kinematicData.zbdVec = zbdVec;
kinematicData.xbVec = xbVec;
kinematicData.xbdVec = xbdVec;

% COM motion
kinematicData.zcVec = zcVec;
kinematicData.zcdVec = zcdVec;
kinematicData.xcVec = xcVec;
kinematicData.xcdVec = xcdVec;

% Foot position
kinematicData.zfootVec = zfootVec;
kinematicData.xfootVec = xfootVec;

end
