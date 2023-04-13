% Seed Correlation (mixed-effects) with sub-cortical Atlas voxels
% by 4 marmosets of resting-state fMRI

function marmoAwakeFingerprint
    compName = 'Aw';
    compNum = 30;
    stat = 'Awake';
    compnames = {'FPN','DMN'};
    compids = [8, 10];

    % atlas of cube clusters
    checkSeedFingerprint(1, compName, compNum, stat, compnames, compids);
end

function checkSeedFingerprint(atlasSize, compName, compNum, stat, compnames, compids)
    path = 'results/marmo/';
    cmap = turbo;
%    cmap = hot;

    % load background nii
    tempNii = 'data\sp2_avg_mri_exvivo_t2wi_v1.0.0AnethfMRI.nii.gz';
    tempinfo = niftiinfo(tempNii);
    tempV = niftiread(tempinfo);
    % load subcortical atlas 
    sub = 'data\sp2_label_512_v1.0.0SubAnethfMRI.nii';
    subinfo = niftiinfo([sub '.gz']);
    subV = niftiread(subinfo);
    subIdx = find(subV(:)>0);

    sessionName = ['seedFingerprintMarmo' stat];

    % load component time-series (Zall) (calculation is omitted demo version)
    zallmat = [path sessionName 'Zall.mat'];
    load(zallmat);

    % contrasts
    contrasts = cell(length(compnames),1);
    for i=1:length(compnames)
        contrasts{i} = nan(compNum,1);
        contrasts{i}(compids(i)) = 1;
    end
    Pth = 0.05;

    zallmat = [path sessionName '.mat'];
    if exist(zallmat,'file')
        % load beta volumes
        load(zallmat);
    else
        [B2, RSS2, T2, df] = calcSeedCorrMixed(CY, CS);

        % output beta matrix
        save(zallmat,'df','B2','RSS2','T2','-v7.3');
    end

    % show histogram
    T3 = T2(:);
    figure; histogram(T3,'EdgeColor','none'); title([sessionName ' T-values']);
    
    % get contrast image
    [Vsp, ~, ~, ~, tm, tc, mr] = plotSeedCorrImage(compnames, contrasts, df, T2, Pth, subV, tempV, (atlasSize==1), [sessionName 'mix-Corr'], 'none', [NaN 20], [NaN 20], cmap, [9:31]);

    % save T-value NIfTI volume
%    saveSeedGlmContrastNii(tempinfo, contnames, Vsp, [path sessionName 'mixCorr']);
    
    % finger print check
    fingerPrintCheck(compnames, contrasts, T2, B2, subV, path, sessionName)
end

function fingerPrintCheck(contnames, contrasts, T2, B2, atlasV, path, sessionName)
    ROIs = {'CAU','PUT','HIPPO','AMY','SC','IC','LGN','ANT','LD','MD','VA','VL','VP','PUL'};
    % prepare
    aIdx = find(atlasV(:) > 0);
    V = single(atlasV);
    V(:) = nan;
    amax = max(atlasV(:));
    fgIdx = {};
    fgmaxsz = 0;
    for j=1:amax
        fgIdx{j} = find(atlasV(:) == j);
        if fgmaxsz < length(fgIdx{j}), fgmaxsz = length(fgIdx{j}); end
    end

    % get fingerprint
    Ts={};
    for i=1:length(contrasts)
        T = nan(fgmaxsz,amax);
        tg = contrasts{i};
        Vt = V;
        T2c = T2 .* tg';
        T2p = max(T2c,[],2);
        Vt(aIdx) = T2p;
        for j=1:amax
            T(1:length(fgIdx{j}),j) = Vt(fgIdx{j});
        end
        figure; boxplot(T,ROIs); title(['marmo awake rest finger print (T-value): ' contnames{i}])
        ylim([-5 20]);
        Ts{i} = T;
    end
end
