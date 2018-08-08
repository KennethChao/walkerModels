function [meshgridK,meshgridDelta] = createMeshGrid(optParms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x = linspace(optParms.kMin, optParms.kMax, optParms.sampledNumberK);
y = linspace(optParms.deltaMin, optParms.deltaMax, optParms.sampledNumberDelta);
[meshgridK, meshgridDelta] = meshgrid(x, y);
end

