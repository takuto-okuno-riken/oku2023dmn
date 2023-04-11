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

    % load background nii
    tempNii = 'data\sp2_avg_mri_exvivo_t2wi_v1.0.0Audio.nii.gz';
    tempinfo = niftiinfo(tempNii);
    tempV = niftiread(tempinfo);

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

    % calc 2nd-level estimation
    B1 = [];
    X2 = [];
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

        % 2nd-level design matrix
        X2 = [X2; eye(size(B2,2))];
    end
    B1(isnan(B1)) = 0; % there might be nan
    if isempty(B1), return; end

    for tuM = 8:8
        betaBmat = [path cubename prefix 'D-Tukey' num2str(tuM) 'full.mat'];
        if exist(betaBmat,'file')
            % load beta volumes
            load(betaBmat);
        else
            % calc 2nd-level estimation
            [B, RSS, df, X2is, tRs] = calcGlmTukey(B1, X2, tuM);

            % output beta matrix
            save(betaBmat,'B','RSS','X2is','tRs','df','-v7.3');
        end

        % GLM contrast image
        [Ts, Tth, Vts, Vfs, Tmaxs, Tcnts, mrvs] = plotGlmContrastImage(contnames, contrasts, B, RSS, X2is, tRs, df, Pth, atlasV, tempV, ...
            (atlasSize==1), ['GLM6marmoAudD ' '2nd-mix-Tukey' num2str(tuM) 'full'], 'none', rangePlus, rangeMinus, [], [], []);

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

