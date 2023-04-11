%%
% Calculate GLM Contrast Image with prewhitening.
% based on K.J.Friston et al. (2000), M.W.Woolrich et al. (2001), K.J.Worsley (2001)
% returns cells of T-value matrix (Ts)
% input:
%  Cs           cells of contrast vectors (contrasts (predictor size) x 1)
%  B            predictor variables (node x predictor variables)
%  RSS          Residual Sum of Squares (node x 1)
%  X2is         Vector or single value of inv(X' * X) for contrast
%  tRs          Vector or single value of trace(R) for contrast

function [Ts] = calcGlmContrastImage(Cs, B, RSS, X2is, tRs)
    % GLM contrast image
    % T = (c' * B) / sqrt(c' * X' * inv(V) * X * c * (RSS / trace(R)))
    Ts = {};
    roiNum = size(RSS,1);
    T2 = zeros(roiNum,1,'single');

    for j=1:length(Cs)
        c = Cs{j};
%        for i=1:roiNum
        parfor i=1:roiNum
            if size(X2is,1) == roiNum
                X2i = squeeze(X2is(i,:,:));
                d = sqrt(c' * X2i * c);
            else
                d = sqrt(c' * X2is * c);
            end
            if size(tRs,1) == roiNum
                se2 = sqrt(RSS(i) / tRs(i));
            else
                se2 = sqrt(RSS(i) / tRs);
            end
            T2(i) = (c' * B(i,:)') / (d * se2);
        end
        Ts{j} = T2;
    end
end
