%%
% Plot GLM Contrast Image with prewhitening.
% based on K.J.Friston et al. (2000), M.W.Woolrich et al. (2001),
% K.J.Worsley (2001), A.Kawaguchi (2017). 
% returns threshold value for T-value (Tth),
%   cells of range thresholded T-value 3D voxels (Vts), cells of full T-value 3D voxels (Vfs)
% input:
%  contnames    cells of contrast names
%  Ts           cells of T-value matrix
%  thParam      threshold params for T-value matrix {df, Pth, corrMeth}
%            df:       degree of freedom
%            Pth:      P-value threshold for T-value matrix
%            corrMeth: family-wise error rate (FWER) correction method for P-value ([default] 'none','bonf','sidak','holm-bonf','holm-sidak')
%  clParam      clustering threshold params {extK, rFWHM, cdt, conn, corrMeth}
%            extK:     extent threshold (voxels)
%            rFWHM:    resel FWHM (should be estimated by estimateSmoothFWHM.m)
%            cdt:      cluster defining threshold (default: 3.3)
%            conn:     adjacent voxel connection number (6, [default] 18, 26)
%            corrMeth: family-wise error rate (FWER) correction method for clustering ('bonf', 'sidak', [default] 'poisson')
%  atlasV       3D atlas volume (or {0, 1} mask volume)
%  isMask       if atlasV is mask volume, then true. otherwise, atlasV has regional info, and also used as mask.
%  isRtoL       X axis is right to left (default: true)
%  backV        3D background volume for plotNifti3DAxes (default: [])
%  sessionName  session name for title (optional)
%  rangePlus    plus T-value range ([min max]) (optional)
%  rangeMinus   minus T-value range ([min max]) (optional)
%  cmap         color map for 3D plot (default: turbo)
%  zidx         show z-index slices (default: [])
%  flatXY       show functional flat map (default: [])

function [Tth, Vts, Vfs, Tmaxs, Tcnts] = plotGlmContrastImage(contnames, Ts, thParam, clParam, atlasV, isMask, isRtoL, backV, sessionName, rangePlus, rangeMinus, cmap, zidx, flatXY)
    if nargin < 14, flatXY = []; end
    if nargin < 13, zidx = []; end
    if nargin < 12, cmap = []; end
    if nargin < 11, rangeMinus = []; end
    if nargin < 10, rangePlus = []; end
    if nargin < 9, sessionName = ''; end
    if nargin < 8, backV = []; end
    if nargin < 7, isRtoL = true; end

    % init threshold parameters
    df = thParam{1}; Pth = thParam{2}; corrMeth = 'none';
    if length(thParam)>=3, corrMeth = thParam{3}; end
    extK = 0; rFW = []; conn = 18; cdt = 4; clCorrMeth = 'poisson';
    if length(clParam)>=1, extK = clParam{1}; end
    if length(clParam)>=2, rFW = clParam{2}; end
    if length(clParam)>=3, cdt = clParam{3}; end
    if length(clParam)>=4, conn = clParam{4}; end
    if length(clParam)>=5, clCorrMeth = clParam{5}; end
    resel = prod(rFW);

    if isempty(rangePlus), rangePlus = [nan nan]; end
    if isempty(rangeMinus), rangeMinus = [nan nan]; end
    if isempty(cmap), cmap = turbo; end
    orgRp = rangePlus; orgRm = rangeMinus;

    % T-value threshold (Bonferroni correction / Šidák correction)
    if strcmp(corrMeth,'bonf') || strcmp(corrMeth,'sidak')
        m = size(Ts{1},1);
        if strcmp(corrMeth,'bonf')
            BPth = Pth / m;
        else
            % Šidák correction - almost same as bonferroni
            BPth = 1 - power(1 - Pth, 1/m);
        end
        Tth = abs(tinv(BPth,df));
        disp(['Height threshold: P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), T-value=' num2str(Tth)])
    elseif ~strcmp(corrMeth,'holm-bonf') && ~strcmp(corrMeth,'holm-sidak')
        Tth = abs(tinv(Pth,df));
        disp(['Height threshold: P-value=' num2str(Pth) ', T-value=' num2str(Tth)])
    end

    Tmaxs = [];
    Tcnts = [];

    for j=1:length(contnames)
        T2 = Ts{j};

        % transform to Z
%        s = std(T2,1);
%        m = mean(T2);
%        Z = (T2 - m)/ s;
%        figure; histogram(Z);
%        Z(Z>-Zth & Z<Zth) = nan;
%        T2 = Z;
%        T2(T2>-Tth & T2<Tth) = nan;

        if isMask == 1
            aIdx = find(atlasV(:) > 0);
            V2 = single(atlasV);
            V2(:) = nan;
            V2(aIdx) = T2;
        else
            V2 = getNifti4DFromRoiTS(T2, atlasV);
        end

        % Holm–Bonferroni method
        if strcmp(corrMeth,'holm-bonf') || strcmp(corrMeth,'holm-sidak')
            m = size(T2,1);
            T2s = sort(T2(:),'descend');
            for k=1:length(T2s)
                if strcmp(corrMeth,'holm-bonf')
                    BPth = Pth / (m + 1 - k);
                else
                    BPth = 1 - power(1 - Pth, 1/(m + 1 - k));
                end
                Tth = abs(tinv(BPth,df));
                if T2s(k) <= Tth
                    break;
                end
            end
            disp(['Height threshold: P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), T-value=' num2str(Tth)])
        end
        V2p = V2; V2m = -V2;
        V2p(V2p<Tth) = 0;   % thresholded T-value
        V2m(V2m<Tth) = 0;   % thresholded T-value

        % clustering threshold
        if extK > 0
            b = 0.5334942;
            filter = images.internal.getBinaryConnectivityMatrix(conn);
            cl = {};
            peaks = {};
            for ii=1:2
                if ii==1, BW = V2p; else BW = V2m; end
                BW(isnan(BW)) = 0;
                BWconv = ( convn(single(BW),filter,'same') >= cdt ) & logical(BW); % find clusters of voxels with more than CDT neighbors
                L = bwlabeln(BWconv,conn);
                clFull = max(L(:));
                for k=1:clFull
                    idx = find(L==k);
                    s = length(idx);
                    if s <= extK
                        if ii==1, V2p(idx) = 0; else V2m(idx) = 0; end
                        continue;
                    end
                    cl{end+1} = s;
                    if ii==1
                        [pV, I] = max(V2p(idx));
                    else
                        [pV, I] = max(V2m(idx)); pV = -pV;
                    end
                    pI = idx(I);
                    [i1,i2,i3] = ind2sub(size(V2p), pI);
                    peaks{end+1} = [i1, i2, i3, pV];
                end
            end
            clExp = length(cl); %expected cluster number
            p = exp(-b * power(extK/resel * cdt*(cdt*cdt-1),2/3)); % uncorrected
            if strcmp(clCorrMeth,'poisson'), Pcr = min(1, 1-poisscdf(0,(clExp + eps)*p)); % Poisson clumping heuristic (used in spm_P_RF.m of SPM12)
            elseif strcmp(clCorrMeth,'sidak'), Pcr = min(1, 1-(1-p)^clExp); % Šidák correction
            else, Pcr = min(1, p * clExp); end % Bonferroni correction
            disp(['Extent threshold: k=' num2str(extK) ' voxels, P-value=' num2str(p) ' (' num2str(Pcr) ')'])
            disp(['Expected cluster list (num=' num2str(length(cl)) '), degree of freedom=' num2str(df) ', FWHM={' num2str(rFW(1)) ' ' num2str(rFW(2)) ' ' num2str(rFW(3)) '} voxels']);
            for k=1:length(cl)
                p = exp(-b * power(cl{k}/resel * cdt*(cdt*cdt-1),2/3)); % uncorrected
                if strcmp(clCorrMeth,'poisson'), Pcr = min(1, 1-poisscdf(0,(clExp + eps)*p)); % Poisson clumping heuristic
                elseif strcmp(clCorrMeth,'sidak'), Pcr = min(1, 1-(1-p)^clExp); % Šidák correction
                else, Pcr = min(1, p * clExp); end % Bonferroni correction
                pk = peaks{k};
                disp([num2str(k) ') k=' num2str(cl{k}) ' voxels, P uncorr=' num2str(p) ', P fwe-corr=' num2str(Pcr) ', peak (' num2str(pk(1)) ',' num2str(pk(2)) ',' num2str(pk(3)) ')=' num2str(pk(4))]);
            end
        end
        V2p(V2p==0) = nan; % set nan to visibility
        V2m(V2m==0) = nan; % set nan to visibility

        % set auto range
        if isnan(orgRp(1))
            rangePlus(1) = Tth;
        end
        if isnan(orgRp(2))
            m = max(V2p(~isinf(V2p)));
            if m>Tth, rangePlus(2) = m; else, rangePlus(2) = Tth + 1; end
        end
        if isnan(orgRm(1))
            rangeMinus(1) = Tth;
        end
        if isnan(orgRm(2))
            m = max(V2m(~isinf(V2m)));
            if m>Tth, rangeMinus(2) = m; else, rangeMinus(2) = Tth + 1; end
        end

        % show figure
        figure; plotNifti3DAxes(V2p,'max',rangePlus,cmap,backV,gray,isRtoL);
        sgtitle(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}],'Color','white');
        figure; plotNifti3DAxes(V2m,'max',rangeMinus,cmap,backV,gray,isRtoL);
        sgtitle(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}],'Color','white');

        if ~isempty(zidx)
            figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}],'Color','white');

            figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}],'Color','white');
        end

        if ~isempty(flatXY)
            T2p = T2; T2m = -T2;
            T2p(T2p<Tth) = nan; T2m(T2m<Tth) = nan;   % thresholded T-value
            figure; plotNifti3Dflatmap(T2p, atlasV, isMask, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}]);

            figure; plotNifti3Dflatmap(T2m, atlasV, isMask, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}]);
        end

        % show Tmax
        tmax = max(V2(:));
        tcnt = length(find(V2>=Tth));
        disp(['Tmax of ' sessionName ' : ' contnames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt)])

        % output 
        V2p(isnan(V2p)) = 0;
        V2m(isnan(V2m)) = 0;
        V2t = V2p - V2m;
        V2t(V2t>rangePlus(2)) = rangePlus(2);
        V2t(V2t<-rangeMinus(2)) = -rangeMinus(2);
        V2t(-rangeMinus(1)<V2t & V2t < rangePlus(1)) = 0;
        Vts{j} = V2t;
        Vfs{j} = V2;
        Tmaxs(j) = tmax;
        Tcnts(j) = tcnt;
    end
end
