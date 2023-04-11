%%
% get GLM HRF (hemodynamic response function) design matrix
% input:
%  onsets       cells of task set start time
%  durations    cells of task set duration
%  frames       fMRI time frames
%  TR           fMRI TR
%  res          sampling resolution of HRF
%  sp           sampling starting point (in resolution)
%  hrf          Canonical hemodynamic response function (optional)

function [X, U] = getGlmHRFDesignMatrix(onsets, durations, frames, TR, res, sp, hrf)
    if nargin < 7, hrf = []; end

    % get Canonical hemodynamic response function
    if isempty(hrf)
        dt = TR / res;
        [t, hrf] = getGlmHRF(dt);
%       figure; plot(t,hrf); % example plot
    end

    taskNum = length(onsets);
    X = zeros(frames*res, taskNum);
    U = zeros(frames*res, taskNum);
    for k=1:taskNum
        onset = onsets{k};
        duration = durations{k};

        for i=1:length(onset)
            t1 = ceil(onset(i) / TR * res);
            t2 = ceil(t1 + duration(i) / TR * res);
            U(t1:t2,k) = 1;
        end

        % get design matrix
        C = conv(U(:,k), hrf);
        X(:,k) = C(1:frames*res);
    end
    
    % final sampling
    X = X(sp:res:end,:);
    U = U(sp:res:end,:);
end

