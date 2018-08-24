% The script to analyze the numerical Poincare map using exported 
% SCS simulation data and plot the stable region for 2 DOF vertical hopper.
close all;
clc;
clear;


% Change cd path for plotting result with different parameters or controls

cd Data/continuousPD/Normal
% cd Data/continuousPD/Period0.5
% cd Data/continuousPD/Period1.0
% cd Data/discretePD/Period1.0
listing = dir(pwd);

kpMax = 5;
kdMax = 5;
sampledNum = 20;


kpSampled = 0:kpMax/(sampledNum-1):kpMax;
kdSampled = 0:kdMax/(sampledNum-1):kdMax;


[X,Y] = meshgrid(kpSampled,kdSampled);

Map = nan(size(X));
for i = 1:length(listing)
    [filepath,name,ext] = fileparts(listing(i).name);
    if strcmp(ext,'.mat')
        load(listing(i).name)

        q_pitch = root.Runner.q_pitch-pi;
        qd_pitch = root.Runner.qd_pitch;
        fz = root.Runner.contactPoint0_fZ;
        diffFz = diff(fz);
        positiveIndex = find(diffFz>0)+1;
%%
        q_pitchPoincareSection = q_pitch(positiveIndex);
        q_pitchPoincareSectionStepN = q_pitchPoincareSection(1:(end-1));
        q_pitchPoincareSectionStepNplus1 = q_pitchPoincareSection(2:end);

        qd_pitchPoincareSection = qd_pitch(positiveIndex);
        qd_pitchPoincareSectionStepN = qd_pitchPoincareSection(1:(end-1));
        qd_pitchPoincareSectionStepNplus1 = qd_pitchPoincareSection(2:end);

        xStepN = [q_pitchPoincareSectionStepN;qd_pitchPoincareSectionStepN];
        xStepNplus1 = [q_pitchPoincareSectionStepNplus1;qd_pitchPoincareSectionStepNplus1];

        A = xStepN'\xStepNplus1';

        [ev, D]=eig(A);
        x = diag(D);
        clc;
        absx = abs(x)

%%
%         unshiftedX = [q_pitch(1:(end-1));qd_pitch(1:(end-1))];
%         shiftedX = [q_pitch(2:end);qd_pitch(2:end)];
%         
%         [evals,modes,Atilde] = tdmd(unshiftedX,shiftedX,10);
%         evals
%%
        str = name;
        C = strsplit(str,'_');
        kpIndex = str2double(C(3))+1;
        kdIndex = str2double(C(5))+1;
        
        
        if absx(1)<1 && absx(2)<1
        Map(kdIndex,kpIndex) = max(absx);            
        else
        Map(kdIndex,kpIndex) = nan;
        
        end
%         if rank(A) ==1
%             Map(kdIndex,kpIndex) = nan;
%         end
%         figure()
% 
%         plot(q_pitchPoincareSection,qd_pitchPoincareSection,'o')
    end
    
end
figure()
    h=surf(X,Y,Map);
    set(h,'EdgeColor','none','FaceColor','interp');%,
    
    set(h,'DefaultTextFontName','Times','DefaultTextFontSize',18,...
       'DefaultAxesFontName','Times','DefaultAxesFontSize',18,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)
    caxis([0 1])
    colorbar
    %str = sprintf('$T = %.3f$, $a = %.3f$, $f = %0.3f$', T,a,f);
    %title(str,'Interpreter','latex','FontName','Times','FontWeight','Normal')    
    
    xlabel('$kp$','Interpreter','latex')
    ylabel('$kd$','Interpreter','latex')
    view(2)
    axis([0 kpMax 0 kdMax])
cd ../../../