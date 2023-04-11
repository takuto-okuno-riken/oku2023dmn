%%
% Calculate 2nd-Level Mixed-Effects General Linear Model with prewhitening (Tukey-Taper, Auto-correlation estimation Version)
% based on M.W.Woolrich et al. (2001), K.J.Friston et al. (2002),
% C.F.Beckmann (2003), J.A.Mumford (2007).
%   Y = X * B1 + e
%   B1 = X2 * B2 + e2
% This function is used for multi session.
% returns predictor variables (B), Residual Sum of Squares (RSS), degree of freedom (df),
%   inv(X' * X) for contrast (X2is), trace(R) for contrast (tRs)
% input:
%  Y         cells of ROI or voxel time series (time series x node)
%  X         cells of design matrix (time series x predictor variables)
%  tuM       Tukey-Taper window size (default: sqrt(time series length))
%  contLen   contrast length (to ignore nuisance) (option)

function [B, RSS, df, X2is, tRs] = calcGlmMixedTukey(CY, CX, tuM, contLen)
    if nargin < 4, contLen = size(CX{1},2); end
    if nargin < 3, tuM = floor(sqrt(size(X, 1))); end

    disp(['process Mixed-Effects GLM with Tukey-Taper(' num2str(tuM) ') estimation ...']);
    tc = tic;

    B1 = [];
    X2 = [];
    for j=1:length(CY)
        % 1st-level estimation
        B2 = calcGlmTukey(CY{j}, CX{j}, tuM);

        % 2nd-level Y vector
        B2 = B2(:,1:contLen);
        B1 = [B1; B2'];

        % 2nd-level design matrix
        X2 = [X2; eye(size(B2,2))];
    end

    % calc 2nd-level estimation
    [B, RSS, df, X2is, tRs] = calcGlmTukey(B1, X2, tuM);
    t = toc(tc);
    disp(['done t=' num2str(t) 'sec'])
end
