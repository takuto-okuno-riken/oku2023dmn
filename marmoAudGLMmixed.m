% GLM with cube Atlas (generated from resampling) ROIs
% by marmoset audio-task subjects
% multi-subject. mixed-effects with Tukey-Taper.
% need to exec marmoAudGLMmixed.m first.
function marmoAudGLMmixed
    atlasSize = 1;
    smooth = 's34';
    filter = 'hf';
    prefix = [smooth 'wa'];

    tfmribase = {'data'};
    sbjs = {
        'M3_1', 'M3_2', 'M3_3', ...
        'M4_1', 'M4_2', 'M4_3', ...
        };
    path = 'results/glm/';

    % template nii for save data
    tempNii = 'data\s34waM3_1.nii.gz';
    % load background nii
    backNii = 'data\sp2_avg_mri_exvivo_t2wi_v1.0.0Audio.nii.gz';
    backinfo = niftiinfo(backNii);
    backV = niftiread(backinfo);

    % load atlas of cube clusters
    cubename = ['marmoAuCube' num2str(atlasSize)'];
    atlas = ['data/' cubename 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);

    % contrast image params
    contnames = {'audio'};
    contrasts = {[1 0]'};
    Pth = 0.001; % pvalue threshold
    rangePlus = [nan 15];
    rangeMinus = [nan 15];

    tuM = 8; % might be best Tukey-window
    isRtoL = true;  % this is SPM12 output

    % calc 2nd-level estimation
    B1 = [];
    X2 = [];
    FWHMs = [];
    for i=1:length(sbjs)
        betaBmat = [path cubename prefix sbjs{i} 'C-Tukey' num2str(tuM) '.mat'];
        if ~exist(betaBmat,'file')
            disp(['file not found. please calc individual sessions first : ' betaBmat])
            continue;
        end

        % load beta volumes
        f = load(betaBmat);
        % 2nd-level Y vector
        B2 = f.B2(:,[1,6]); % include design and intercept (we need more than 8 length for tukey taper)
        B1 = [B1; B2'];
        FWHMs = [FWHMs; f.FWHM];

        % 2nd-level design matrix
        X2 = [X2; eye(size(B2,2))];
    end
    B1(isnan(B1)) = 0; % there might be nan
    FWHMs = mean(FWHMs,1); % let's take the mean of FWHM.

    for tuM = 8:8
        betaBmat = [path cubename prefix 'D-Tukey' num2str(tuM) 'full.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            % calc 2nd-level estimation
            [B, RSS, df, X2is, tRs, R] = calcGlmTukey(B1, X2, tuM);

            [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, atlasV);

            % output beta matrix
            save(betaBmat,'B','RSS','X2is','tRs','recel','FWHM','df','-v7.3');
        end

        % GLM contrast images
        Ts = calcGlmContrastImage(contrasts, B, RSS, X2is, tRs);

        % GLM contrast image
        thParam = {df, Pth};
        clParam = {69, FWHMs}; % clustering parameter for GLM contrast
        [Tth, Vts, Vfs, Tmaxs, Tcnts] = plotGlmContrastImage(contnames, Ts, thParam, clParam, atlasV, (atlasSize==1), isRtoL, backV, ...
            ['GLM6marmoAudD ' '2nd-mix-Tukey' num2str(tuM) 'full'], rangePlus, rangeMinus, [], [], []);

        % save T-value NIfTI volume
        saveContrastNii(tempNii,contnames,Vts,path,[cubename prefix 'D_2nd-mix-Tukey' num2str(tuM) 'th' 'full']);
    end
end

%%
function saveContrastNii(tfmri, contnames, V2s, path, outname)
    info = niftiinfo(tfmri);
    info.ImageSize = info.ImageSize(1:3);
    info.PixelDimensions = info.PixelDimensions(1:3);
    info.raw.dim(1) = 3;
    info.raw.dim(5) = 1;
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    for j=1:length(contnames)
        fname = [path outname '_' contnames{j} '.nii'];
        niftiwrite(V2s{j},fname,info,'Compressed',true);
    end
end

