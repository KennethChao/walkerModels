% close all
% figure()
hold on
% for k = 1:searchingVarLength
%     for i = 1:sampledNumber1
%         plot(stablePhiStack(1,:,i,k)-pi/2,stablePhiStack(2,:,i,k),'-o')
%     end
%     axis([-0.15 0 0  0.05])
% end

for k = 1:searchingVarLength
    surf(X, Y, unstableDataStack(:, :, k))
end
%
%     axis([kMinPlot kMaxPlot deltaMinPlot deltaMaxPlot])
