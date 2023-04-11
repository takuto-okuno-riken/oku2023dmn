% GLM with cube Atlas (generated from resampling) ROIs
% by marmoset audio-task subjects
% single session. with Tukey 8 window
function marmoAudGLMindividual
    res = 16; % HRF sampling resolution
    sp = 8;   % HRF sampling starting point
    TR = 3;
    
    atlasSize = 1;
    smooth = 's34';
    prefix = [smooth 'wa'];

    tfmribase = {'data'};
    sbjs = {
        'M3_1', 'M3_2', 'M3_3', ...
        'M4_1', 'M4_2', 'M4_3', ...
        };

    path = 'results/glm/';
    if ~exist(path, 'dir')
        mkdir(path)
    end

    hpfTh = 1 / 128; % highpass filter threthold (Hz)
    tuM = 8; % for file check

    % get marmoset HRF
    dt = TR / res;
    [t, hrf] = getGlmHRF(dt, 4.5, 4.5 * 16/6, 24, 0.167, 0.9);
    [t2, hrf2] = getGlmHRF(dt); % human's HRF;
    figure; plot(t,hrf); hold on; plot(t2,hrf2); hold off;
    title('Canonical hemodynamic response function for Marmoset');

    frames = 205;
    for k=1:20
        ons(k) = 15 + (k-1)*30;
        dur(k) = 15;
    end
    onsets{1} = ons;
    durations{1} = dur;
    [Xorg, U] = getGlmHRFDesignMatrix(onsets, durations, frames, TR, res, sp, hrf);
%    figure; imagesc([Xorg ones(size(Xorg,1),1)]); colorbar;

    % load background nii
    tempNii = 'data\sp2_avg_mri_exvivo_t2wi_v1.0.0Audio.nii.gz';
    tempinfo = niftiinfo(tempNii);
    tempV = niftiread(tempinfo);

    % generate atlas of cube clusters
    cubename = ['marmoAuCube' num2str(atlasSize)'];
    atlas = ['data/' cubename 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);

    % for Nuisance Signal Regression
    csfF = 'data\csfAudio.nii.gz';
    csfinfo = niftiinfo(csfF);
    csfV = niftiread(csfinfo);
    csfV = single(csfV) / 255; % to [0 1] range
    wmF = 'data\whiteAudio.nii.gz';
    wminfo = niftiinfo(wmF);
    wmV = niftiread(wminfo);
    wmV = single(wmV) / 255; % to [0 1] range
    gsV = tempV;
    gsV(gsV>=1) = 1;
    gsV(gsV<1) = 0;

    for i=1:length(sbjs)
        % check session processed or not
        betaBmat = [path cubename prefix sbjs{i} 'C-Tukey' num2str(tuM) '.mat'];
        if exist(betaBmat,'file')
            disp(['file found : ' betaBmat]);
            continue;
        end
        
        % read task-fMRI volume
        tfmri = [tfmribase{1} '/' prefix sbjs{i} '.nii.gz'];
        if ~exist(tfmri,'file')
            disp(['file not found. please put NIfTI data. (skipped) : ' tfmri]);
            continue;
        end
        disp(['loading : ' tfmri]);
        info = niftiinfo(tfmri);
        V = single(niftiread(info));

        % ROI time-series from task-fMRI
        disp('apply mask atlas...');
        if atlasSize == 1
            aIdx = find(atlasV(:) > 0);
            A = reshape(V,[],size(V,4));
            Z = A(aIdx,:);
        else
            Z = getRoiTSFromNifti4D(V, atlasV, 'mean');
        end
        Z = Z - nanmean(Z,2);
%                Z = convert2SigmoidSignal(Z); % both work, convert or not convert
        
        frames = size(Z,2);
        TR = info.PixelDimensions(4);
%{
        % sample plot of voxels
        figure;
        for k=1:10
            hold on; plot(squeeze(Z(15+k,:))); hold off;
        end
%}

        % get design matrix
        X = Xorg;
        
        % high pass filter as preprocessing step (M.W.Woolrich, 2001) type.
        disp(['apply highpass filter (' num2str(hpfTh) ' Hz) : tfMRI and design matrix...']);
        Z = highpass(Z',hpfTh,1/TR);
        Z = Z';
        % get Nuisance time-series (CSF, WM, Global Signal, Global Mean)
        Xn = getGlmNuisanceTimeSeries(V, csfV, wmV, gsV);
        figure; imagesc([X Xn ones(size(X,1),1)], [-0.4, 1.2]); colorbar;
        % K*W*Y = K*W*B*X + K*W*e case. (SPM type. K.J.Friston, 2000)
        % actually, FLS (FEAT) do this as default option
        X = highpass(X,hpfTh,1/TR);
        X = [X Xn]; % Nuisance should be raw

        % check Tukey range
        checkTukeyRange(Z, X, path, cubename, atlasSize, prefix, atlasV, tempV, sbjs{i});
    end
end


function checkTukeyRange(Zall, Xall, path, cubename, atlasSize, prefix,  atlasV, tempV, subject)
    % contrast image params
    contnames = {'audio'};
    contrasts = {[1 0 0 0 0 0]'};
    Pth = 0.05; % pvalue threshold

    tuMrange = 8:8;

    % ---------------------------------------------------------------------
    % AR estimation with Tukey-Taper of frequency domain
    for tuM = tuMrange
        betaBmat = [path cubename prefix subject 'C-Tukey' num2str(tuM) '.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            Xt = [Xall, ones(size(Xall,1),1)];
            [B2, RSS, df, X2is, tRs] = calcGlmTukey(Zall', Xt, tuM);

            % output beta matrix
            save(betaBmat,'B2','RSS','X2is','tRs','df','-v7.3');
        end
    
        % get contrast image
        plotGlmContrastImage(contnames, contrasts, B2, RSS, X2is, tRs, df, Pth, atlasV, tempV, (atlasSize==1), ...
            ['GLM6 ' cubename prefix subject 'CTukey' num2str(tuM)]);
    end
end
