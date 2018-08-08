close all
h = figure()

set(h,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

if strcmp(    storedQuantity,'delta')
    axis([kMinFig, kMaxFig, deltaMinFig, deltaMaxFig])
elseif( strcmp(    storedQuantity,'lambda'))
    axis([kMinFig, kMaxFig, -1.5, 2.5])
end

hold on
% for k = 1:searchingVarLength
%     for i = 1:sampledNumber1
%         plot(stablePhiStack(1,:,i,k)-pi/2,stablePhiStack(2,:,i,k),'-o')
%     end
%     axis([-0.15 0 0  0.05])
% end
        cbh=colorbar();
        titleString = 'max(abs($\lambda$))';

ylabel(cbh, titleString,'Interpreter','latex')
        if strcmp(searchingVar,'g')
msg = sprintf('g = %0.2f',gVec(1));
        elseif strcmp(searchingVar,'beta')
msg = sprintf('beta = %0.2^of',betaVec(1));
        end
t = text(kMaxPlot/2,deltaMinPlot/3,msg);
xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta*$','Interpreter','latex')

% Parameters for animated gif
if strcmp(searchingVar,'g')
filename = 'SLIPPER_FixedPointAnimation_VaryingG.gif'; % Specify the output file name
else
filename = 'SLIPPER_FixedPointAnimation_VaryingBeta.gif'; % Specify the output file name    
end
DT = 1;

for k = 1:searchingVarLength
   C = contourf(X, Y, stableDataStack(:, :, k),'LineStyle',':')
%    C = surf(X, Y, stableDataStack(:, :, k))
%     contourf(X, Y, stableDataStack(:, :, k))
colormap parula

        if strcmp(searchingVar,'g')
msg = sprintf('g = %0.2f',gVec(k));
        elseif strcmp(searchingVar,'beta')
msg = sprintf('beta = %0.2f^o',betaVec(k)/pi*180);
        end
set(t,'String',msg)
caxis([round(min(min(min(stableDataStack))),1) 1])      

      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

      if k == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DT); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DT); 
      end    
pause(1.0)
end
%
    axis([kMinPlot kMaxPlot deltaMinPlot deltaMaxPlot])
