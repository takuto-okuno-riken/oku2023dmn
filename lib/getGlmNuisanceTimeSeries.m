function Xn = getGlmNuisanceTimeSeries(V, csfV, wmV, gsV)
    Xn = nan(size(V,4),4);
    for i=1:size(V,4)
        Vi = V(:,:,:,i);
        V1 = Vi .* single(csfV); V1(V1<=0) = nan;
        Xn(i,1) = nanmean(V1(:)); % csf mean
        V1 = Vi .* single(wmV); V1(V1<=0) = nan;
        Xn(i,2) = nanmean(V1(:)); % wm mean
        V1 = Vi .* single(gsV); V1(V1<=0) = nan;
        Xn(i,3) = nanmean(V1(:)); % global signal
        Xn(i,4) = nanmean(Vi(:)); % global mean
    end
    Xn = Xn - mean(Xn,1);
    Xn = Xn / (std(Xn(:,3),1)*4); % normalize around [-1 1] range
end
