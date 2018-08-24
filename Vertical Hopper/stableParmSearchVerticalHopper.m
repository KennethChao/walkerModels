clc;
clear;
close all;

%% Control setting and parameter range



% 'discreteP' 'discretePD' 'continuousP' 'continuousPD'
controlType = 'continuousPD';


kpMax = 1;
kdMax = 1;
sampledNum = 200;


%% Setting for the choice of third parameter

% 'hoppingPeriod' 'dutyFactor'
searchVariable = 'hoppingPeriod';
if strcmp(searchVariable, 'hoppingPeriod')
    aMin = 0.1;
    aMax = 2;
    dutyFactor = 0.4;
elseif strcmp(searchVariable, 'dutyFactor')
    T = 0.25;
    aMin =0.2*T;
    aMax = 0.8*T;    
end

sampledNum2 = 20;

%% Model parameter and Map Initialization
m = 1;
I = 0.5;

kpSampled = 0:kpMax / (sampledNum - 1):kpMax;
kdSampled = 0:kdMax / (sampledNum - 1):kdMax;
aSampled = aMin:(aMax - aMin) / (sampledNum2 - 1):aMax;

if aMin == aMax
    [X, Y] = meshgrid(kpSampled, kdSampled);
    Z = aMin;
    sampledNum2 = 1;
else
    [X, Y, Z] = meshgrid(kpSampled, kdSampled, aSampled);
end
Map = X;

%% Build the map

for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        for k = 1:sampledNum2
            
            if aMin == aMax
                kp = X(i, j);
                kd = Y(i, j);
                if strcmp(searchVariable, 'hoppingPeriod')
                    T = Z;
                    a = T * dutyFactor;
                elseif strcmp(searchVariable, 'dutyFactor')
                    a = Z;
                end
                
                
                f = T / a * m * 9.81;
                K = 1 / 2 * f / I * kp;
                C = 1 / 2 * f / I * kd;
                
                A = PoincareMapExpression(T, a, K, C, controlType);
                try
                    qq = abs(eig(A));
                    if qq(1) < 1 && qq(2) < 1
                        Map(i, j) = max(qq);
                    else
                        Map(i, j) = nan;
                    end
                catch
                    Map(i, j) = nan;
                end
            else
                kp = X(i, j, k);
                kd = Y(i, j, k);
                
                if strcmp(searchVariable, 'hoppingPeriod')
                    T = Z(i, j, k);
                    a = T * dutyFactor;
                elseif strcmp(searchVariable, 'dutyFactor')
                    a = Z(i, j, k);
                end
                
                f = T / a * m * 9.81;
                K = 1 / 2 * f / I * kp;
                C = 1 / 2 * f / I * kd;
                
                A = PoincareMapExpression(T, a, K, C, controlType);
                try
                    qq = abs(eig(A));
                    if qq(1) <= 1 && qq(2) <= 1
                        Map(i, j, k) = max(qq);
                    else
                        Map(i, j, k) = nan;
                    end
                catch
                    Map(i, j, k) = nan;
                end
            end
        end
    end
end

%% Plot result

if aMin == aMax
    h = surf(X, Y, Map);
    set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');
    
    set(h, 'DefaultTextFontName', 'Times', 'DefaultTextFontSize', 18, ...
        'DefaultAxesFontName', 'Times', 'DefaultAxesFontSize', 18, ...
        'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)
    caxis([0 1])
    colorbar
    str = sprintf('$T = %.3f$, $a = %.3f$, $f = %0.3f$', T, a, f);
    title(str, 'Interpreter', 'latex', 'FontName', 'Times', 'FontWeight', 'Normal')
    
    xlabel('$kp$', 'Interpreter', 'latex')
    ylabel('$kd$', 'Interpreter', 'latex')
    view(2)
    axis([0, kpMax, 0, kdMax])
else
    h = slice(X, Y, Z, Map, kpSampled, kdSampled, aSampled);
    
    set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');
    alpha(0.4);
    axis([0, kpMax, 0, kdMax, aMin, aMax])
    set(h, 'DefaultTextFontName', 'Times', 'DefaultTextFontSize', 18, ...
        'DefaultAxesFontName', 'Times', 'DefaultAxesFontSize', 18, ...
        'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)
    colorbar
    xlabel('$kp$', 'Interpreter', 'latex')
    ylabel('$kd$', 'Interpreter', 'latex')
    if strcmp(searchVariable,'hoppingPeriod')
        zlabel('$T$', 'Interpreter', 'latex')
    elseif strcmp(searchVariable, 'dutyFactor')
        zlabel('$\alpha$', 'Interpreter', 'latex')
    end
    
    az = -45;
    el = 69;
    view(az, el);
    view(2)
end


function A = PoincareMapExpression(T, a, K, C, type)

switch type
    case 'discreteP'
        A = [1 - a^2 * K, T; ...
            -2 * a * K, 1];
    case 'discretePD'
        A = [1 - a^2 * K, T - a^2 * C; ...
            -2 * a * K, 1 - 2 * a * C];
    case 'continuousP'
        A = [exp(2^(1 / 2)*(-K)^(1 / 2)*a) / 2 + exp(-2^(1 / 2)*(-K)^(1 / 2)*a) / 2, (exp(2^(1 / 2)*(-K)^(1 / 2)*a) / 2 + exp(-2^(1 / 2)*(-K)^(1 / 2)*a) / 2) * (T - a) + (2^(1 / 2) * exp(2^(1 / 2)*(-K)^(1 / 2)*a) - 2^(1 / 2) * exp(-2^(1 / 2)*(-K)^(1 / 2)*a)) / (4 * (-K)^(1 / 2)); ...
            (2^(1 / 2) * (-K)^(1 / 2) * exp(2^(1 / 2)*(-K)^(1 / 2)*a)) / 2 - (2^(1 / 2) * (-K)^(1 / 2) * exp(-2^(1 / 2)*(-K)^(1 / 2)*a)) / 2, exp(2^(1 / 2)*(-K)^(1 / 2)*a) / 2 + exp(-2^(1 / 2)*(-K)^(1 / 2)*a) / 2 + ((2^(1 / 2) * (-K)^(1 / 2) * exp(2^(1 / 2)*(-K)^(1 / 2)*a)) / 2 - (2^(1 / 2) * (-K)^(1 / 2) * exp(-2^(1 / 2)*(-K)^(1 / 2)*a)) / 2) * (T - a)];
    case 'continuousPD'
        A = [(exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) - C * exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) + exp(-a*(C - (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) + C * exp(-a*(C - (C^2 - 2 * K)^(1 / 2)))) / (2 * (C^2 - 2 * K)^(1 / 2)), (exp(-a*(C - (C^2 - 2 * K)^(1 / 2))) - exp(-a*(C + (C^2 - 2 * K)^(1 / 2)))) / (2 * (C^2 - 2 * K)^(1 / 2)) + ((T - a) * (exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) - C * exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) + exp(-a*(C - (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) + C * exp(-a*(C - (C^2 - 2 * K)^(1 / 2))))) / (2 * (C^2 - 2 * K)^(1 / 2)); ...
            (K * exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) - K * exp(-a*(C - (C^2 - 2 * K)^(1 / 2)))) / (C^2 - 2 * K)^(1 / 2), (exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) + C * exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) + exp(-a*(C - (C^2 - 2 * K)^(1 / 2))) * (C^2 - 2 * K)^(1 / 2) - C * exp(-a*(C - (C^2 - 2 * K)^(1 / 2)))) / (2 * (C^2 - 2 * K)^(1 / 2)) + ((T - a) * (K * exp(-a*(C + (C^2 - 2 * K)^(1 / 2))) - K * exp(-a*(C - (C^2 - 2 * K)^(1 / 2))))) / (C^2 - 2 * K)^(1 / 2)];
end

end