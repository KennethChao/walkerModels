close all
figure()
hold on
for k = 1:length(betaVec)
%     figure()
%     hold on
for i = 1:20
%     plot3(xStack(1,:,i),xStack(2,:,i),xStack(3,:,i),'-o')
    plot(xStack(1,:,i,k)-pi/2,xStack(2,:,i,k),'-o')
end
%     hold off
    axis([-0.15 0 0  0.25])
end