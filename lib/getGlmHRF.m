%%
% get GLM Canonical hemodynamic response function
% returns time range (t), HRF time-series (hrf)
% input:
%  dt                 time resolution (sec) (default:0.045)
%  responseDelay      delay of response (gamma a)(sec) (default:6)
%  underShootDelay    delay of undershoot (gamma a)(sec) (default:16)
%  kernelSec          kernel time length (sec) (default: 32)
%  underShootRatio    ratio of underShoot (default: 0.167)
%  hrfScale           HRF scale (gamma b) (default: 0.9)

function [t, hrf] = getGlmHRF(dt, responseDelay, underShootDelay, kernelSec, underShootRatio, hrfScale)
    if nargin < 6, hrfScale = 0.9; end
    if nargin < 5, underShootRatio = 0.167; end
    if nargin < 4, kernelSec = 32; end
    if nargin < 3, underShootDelay = 16; end
    if nargin < 2, responseDelay = 6; end
    if nargin < 1, dt = 0.045; end
    
    t = 0:dt:kernelSec;
    hrf = gampdf(t,responseDelay,hrfScale) - gampdf(t,underShootDelay,hrfScale) * underShootRatio;
    hrf = hrf'/sum(hrf);
end
