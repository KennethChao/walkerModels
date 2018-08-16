close all
h = figure()

set(h,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

% if strcmp(    storedQuantity,'delta')
%     axis([kMinFig, kMaxFig, deltaMinFig, deltaMaxFig])
% elseif( strcmp(    storedQuantity,'lambda'))
%     axis([kMinFig, kMaxFig, -1.5, 2.5])
% end
c = parula(optParms.searchingVarLength);
hold on
% axis([-0.15 0 -5e-3  5e-3])
        cbh=colorbar();
        titleString = 'max(abs($\lambda$))';
%         if strcmp(searchingVar,'g')
%         tempLabel = get(cbh,'YTickLabel');
% %         set(cbh,'YTickLabel',linspace(gVec(1),gVec(end), length(tempLabel)));
%         titleString = '$\tilde g$';
%         elseif strcmp(searchingVar,'beta')
%         tempLabel = get(cbh,'YTickLabel');
%         set(cbh,'YTickLabel',linspace(betaVec(1),betaVec(end), length(tempLabel))/pi*180);
%         titleString = '$\beta(^o)$';
%         end
if strcmp(optParms.searchingVar,'g')
filename = 'SLIPPER_PendulumFixedPointAnimation_VaryingG.gif'; % Specify the output file name
else
filename = 'SLIPPER_PendulumFixedPointAnimation_VaryingBeta.gif'; % Specify the output file name    
end
DT = 1;
        yLimit = max(max(max(abs(result.stablePhi(2,:,:,:)))))+1e-3;

%         axis([min(min(min(stablePhiStack(1,:,:,:))))-pi/2, 0, -4e-3, 4e-3])
%         caxis([min(min(min(stableDataStack))), 1])

        if strcmp(optParms.searchingVar,'g')
msg = sprintf('g = %0.2f',gVec(1));
        elseif strcmp(optParms.searchingVar,'beta')
msg = sprintf('beta = %0.2^of',betaVec(1));
        else
            msg ='';
        end
t = text(min(min(min(result.stablePhi(1,:,:,:))))+0.1,yLimit-1e-3,msg);
ylabel(cbh, titleString,'Interpreter','latex')
xlabel('$\phi$','Interpreter','latex')
ylabel('$\dot\phi$','Interpreter','latex')
% titleString = sprintf('   $$\\beta = %.1f^o$$',betaVec/pi*180);
% % titleString = '$$\beta$$'
% title(titleString,'Interpreter','latex');
% title('$$Q \geq \frac{I_h H}{I_h H+I_z C}, b_1 \geq b_2$$','interpreter','latex')
for k = optParms.searchingVarLength:-1:1
    resultBuf = nan(size(result.meshgridK));
    for i = 1:optParms.sampledNumberK
        for j = 1:optParms.sampledNumberDelta
            resultBuf(j,i) = result.stableData(j,i, k).maxAbsEigenValue;
%               resultBuf(j,i) = result.stableData(j,i, k).netWork;
%               resultBuf(j,i) = result.stableData(j,i, k).dutyFactor;
%               resultBuf(j,i) = result.stableData(j,i, k).fval;
        end
    end
    
    
        C = contourf(squeeze(result.stablePhi(1,:,:,k)), squeeze(result.stablePhi(2,:,:,k)), resultBuf,'LineStyle',':')

%       axis([0, 1.5*1e-2, 0.05, 0.06])
      caxis([0.0, 1])
      
% msg = sprintf('$$\\tilde g = %0.2f,  \\beta = %0.2f^o$$', optParms.g(k), optParms.beta(k)/pi*180);
%     textHandle = text(0.0052, 0.0522, msg, 'Interpreter', 'latex');
%     msg2 = sprintf('$$\\tilde m_f = %0.2f,  \\tilde r_c = %0.2f$$',optParms.mf, optParms.rc);  
% set(t,'String',msg)      
%           textHandle = text(0.0052, 0.051, msg2, 'Interpreter', 'latex');
%         grid on
             frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      
    


end

