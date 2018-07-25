clc;
clear B;


B = nan(2,size(solution,2));
for i = 1:size(solution,2)
buf =solution(:,i);
buf(~isnan(buf))
uniqeBuf = uniquetol(buf(~isnan(buf)),1e-2);
if ~isempty(uniqeBuf)
B(1:length(uniqeBuf),i) =flipud(uniqeBuf);
end



end

combinedDelta = [B(1,:) fliplr(B(2,:))];
combinedKvec = [kVec fliplr(kVec)];
nanIndex = find(isnan(combinedDelta));
combinedDelta(nanIndex) = [];
combinedKvec(nanIndex) = [];


plot(combinedKvec,combinedDelta)