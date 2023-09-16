%%
% Calculate General Linear Model with prewhitening (Tukey-Taper, Auto-correlation estimation Version)
% based on M.W.Woolrich et al. (2001)
%   Y = X * B + e
%   e ~ N(0, s^2 * V)
%   S * Y = S * X * B + r
%   S = inv(K)
%   V = K * K'
%   V is auto-corerlation matrix estimated by Tukey-Taper
%   r ~ N(0, s^2 * SVS') therefore N(0, s^2 * I)
% This function is used for single session.
% returns predictor variables (B), Residual Sum of Squares (RSS), degree of freedom (df),
%   inv(X' * X) for contrast (X2is), trace(R) for contrast (tRs), full of Residuals (R)
% input:
%  Y         ROI or voxel time series (time series x node)
%  X         design matrix (time series x predictor variables)
%  tuM       Tukey-Taper window size (default: sqrt(time series length))

function [B, RSS, df, X2is, tRs, R] = calcGlmTukey(Y, X, tuM)
    if nargin < 3, tuM = floor(sqrt(size(X, 1))); end

    disp(['process GLM with Tukey-Taper(' num2str(tuM) ') estimation ...']);
    [~, ~, perm, RiQ, dR2i] = regressPrepare(X);

    roiNum = size(Y,2);
    xsz = size(X,2);
    X2is = nan(roiNum,xsz,xsz,'single');
    tRs = nan(roiNum,1,'single');
    B = nan(roiNum,xsz,'single');
    RSS = nan(roiNum,1,'single');
    CR = cell(roiNum,1);
    isOutR = (nargout > 5);
    df = size(X,1) - size(X,2);

    % make Tukey window
    tuWin = zeros(1, tuM);
    for k=1:tuM-1, tuWin(1,k) = 0.5 * (1 + cos(pi*k/tuM)); end

    tc = tic; % check running time 
%    for i=1:roiNum
    parfor i=1:roiNum
        % 1st step OLS regression
        [~, r] = regressLinear(Y(:,i), X, [], [], perm, RiQ, dR2i);

        % if residuals are all zero, probably original signal is
        % just zero. ignore this voxel
        if sum(r==0) == size(r, 1)
            B(i,:) = 0;
            RSS(i) = 0;
            if isOutR, CR{i} = r; end
            continue;
        end

        % calc AR coefficients by Tukey-Taper of frequency domain
%        C = xcov(r',tuM,'unbiased'); % this does not affect 
%        Rxx = C(tuM+1:end) / var(r', 1);
        C = xcorr(r',tuM,'normalized');
        Rxx = C(tuM+1:end);
        Pxx = zeros(1,size(r,1));
        Pxx(1:tuM) = Rxx(1:tuM) .* tuWin(1,1:tuM);
        V1 = toeplitz(Pxx(1:end));
        % chol type
%        Vi1 = inv(V1); % this order does not affect 
%        Ki1 = chol(Vi1,'lower');
        K1 = chol(V1,'lower');
        Ki1 = inv(K1);
        % yule-walker type
%        br = inv(V1) * [Pxx(2:end) 0]';
%        Ki1 = tril(toeplitz([1 -br(1:end-1)']));

        % second time regression
        Ya = Ki1 * Y(:,i);
        Xt = Ki1 * X;
        [b, r] = regressLinear(Ya, Xt);

%        C = xcov(r',tuM,'unbiased'); % this does not affect 
%        Rxx = C(tuM+1:end) / var(r',1);
        C = xcorr(r',tuM,'normalized');
        Rxx = C(tuM+1:end);
        Pxx = zeros(1,size(r,1));
        Pxx(1:tuM) = Rxx(1:tuM) .* tuWin(1,1:tuM);
        V2 = toeplitz(Pxx(1:end));
%{
        K2 = chol(V2,'lower');
        Ki2 = inv(K2);
%       br = Vi * [Pxx(2:end) 0]';
%       Ki2 = tril(toeplitz([1 -br(1:end-1)']));
%}
        B(i,:) = b;
        RSS(i) = r' * r;
        if isOutR, CR{i} = r; end
%{
        C = xcov(r',tuM,'unbiased');
        Rxx = C(tuM+1:end) / var(r',1);
        Pxx(1:tuM) = Rxx(1:tuM) .* tuWin(2,1:tuM);
        V3 = toeplitz(Pxx(1:end));
%}
        % used for contrast
%{
        % get K-inverse by AR yule-walker
        Arc = aryule(R2(i,:),256);
        br = Arc(2:end);
        Ki2 = tril(toeplitz([1 -br zeros(1,size(Xall,1)-256-1)])); % Ki = I - A
        Xt2 = Ki2 * [Xall, ones(size(Xall,1),1)];
        Xt2is(i,:,:) = inv(Xt2' * Xt2);
%}
        X2is(i,:,:) = inv(X' * (V2 \ X));
        IR = eye(size(Xt,1)) - Xt * ((Xt'*Xt) \ Xt');
        tRs(i) = trace(IR);
%{
        % this does not change Tmax and tcnt, but twice slow
        Vi2 = inv(V2);
        Ki2 = chol(Vi2,'lower');
        Xb = Ki2 * [Xall, ones(size(Xall,1),1)];
        Xt2i = inv(Xb'*Xb);
        Xt2is(i,:,:) = Xt2i;
%            tIR = trace(eye(size(Xb,1)) - Xb * Xt2i * Xb');
        tIR = size(Xb,1) - sum((Xb * Xt2i)' .* Xb','all');
        tR2s(i) = tIR;
%}
    end
    if isOutR
        R = nan(roiNum,size(Y,1),'single');
        for i=1:roiNum
            R(i,:) = CR{i};
        end
    end
    t = toc(tc);
    disp(['done t=' num2str(t) 'sec'])
end
