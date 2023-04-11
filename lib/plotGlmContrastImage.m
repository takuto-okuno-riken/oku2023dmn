%%
% Plot GLM Contrast Image with prewhitening.
% based on K.J.Friston et al. (2000), M.W.Woolrich et al. (2001), K.J.Worsley (2001)
% returns cells of T-value matrix (Ts), threshold value for T-value (Tth),
%   cells of range thresholded T-value 3D voxels (Vts), cells of full T-value 3D voxels (Vfs)
% returns 
% input:
%  contnames    cells of contrast names
%  Cs           cells of contrast vectors (contrasts (predictor size) x 1)
%  B            predictor variables (node x predictor variables)
%  RSS          Residual Sum of Squares (node x 1)
%  X2is         Vector or single value of inv(X' * X) for contrast
%  tRs          Vector or single value of trace(R) for contrast
%  df           degree of freedom (df)
%  Pth          P-value threshold for T-value matrix
%  maskV        3D mask volume
%  backV        3D background volume 
%  isFullVoxel  full voxel atlas or not
%  sessionName  session name for title (optional)
%  corrMeth     family-wise error rate (FWER) correction method ('none','bonf','sidak','holm-bonf','holm-sidak')
%  rangePlus    plus T-value range ([min max]) (optional)
%  rangeMinus   minus T-value range ([min max]) (optional)
%  cmap         color map for 3D plot (default: turbo)
%  zidx         show z-index slices (default: [])
%  flatXY       show functional flat map (default: [])

function [Ts, Tth, Vts, Vfs, Tmaxs, Tcnts, mrvs] = plotGlmContrastImage(contnames, Cs, B, RSS, X2is, tRs, df, Pth, maskV, backV, isFullVoxel, sessionName, corrMeth, rangePlus, rangeMinus, cmap, zidx, flatXY)
    if nargin < 18, flatXY = []; end
    if nargin < 17, zidx = []; end
    if nargin < 16, cmap = []; end
    if nargin < 15, rangeMinus = []; end
    if nargin < 14, rangePlus = []; end
    if nargin < 13, corrMeth = 'none'; end
    if nargin < 12, sessionName = ''; end

    if isempty(rangePlus), rangePlus = [nan nan]; end
    if isempty(rangeMinus), rangeMinus = [nan nan]; end
    if isempty(cmap), cmap = turbo; end

    % GLM contrast image
    Ts = calcGlmContrastImage(Cs, B, RSS, X2is, tRs);

    % T-value threshold (Bonferroni correction / Šidák correction)
    if strcmp(corrMeth,'bonf') || strcmp(corrMeth,'sidak')
        m = size(RSS,1);
        if strcmp(corrMeth,'bonf')
            BPth = Pth / m;
        else
            % Šidák correction - almost same as bonferroni
            BPth = 1 - power(1 - Pth, 1/m);
        end
        Tth = abs(tinv(BPth,df));
        disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), T-value=' num2str(Tth)])
    elseif ~strcmp(corrMeth,'holm-bonf') && ~strcmp(corrMeth,'holm-sidak')
        Tth = abs(tinv(Pth,df));
        disp(['P-value=' num2str(Pth) ', T-value=' num2str(Tth)])
    end

    Tmaxs = [];
    Tcnts = [];
    mrvs = [];

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

        if isFullVoxel == 1
            aIdx = find(maskV(:) > 0);
            V2 = single(maskV);
            V2(:) = nan;
            V2(aIdx) = T2;
        else
            V2 = getNifti4DFromRoiTS(T2, maskV);
        end

        % Holm–Bonferroni method
        if strcmp(corrMeth,'holm-bonf') || strcmp(corrMeth,'holm-sidak')
            m = size(RSS,1);
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
            disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), T-value=' num2str(Tth)])
        end
        V2p = V2; V2m = -V2;
        V2p(V2p<Tth) = nan;   % thresholded T-value
        V2m(V2m<Tth) = nan;   % thresholded T-value

        % set auto range
        if isnan(rangePlus(1))
            rangePlus(1) = Tth;
        end
        if isnan(rangePlus(2))
            m = max(V2p(~isinf(V2p)));
            if m>Tth, rangePlus(2) = m; else, rangePlus(2) = Tth + 1; end
        end
        if isnan(rangeMinus(1))
            rangeMinus(1) = Tth;
        end
        if isnan(rangeMinus(2))
            m = max(V2m(~isinf(V2m)));
            if m>Tth, rangeMinus(2) = m; else, rangeMinus(2) = Tth + 1; end
        end

        % show figure
        figure; plotNifti3DAxes(V2p,'max',rangePlus,cmap,backV);
        sgtitle(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}],'Color','white');
        figure; plotNifti3DAxes(V2m,'max',rangeMinus,cmap,backV);
        sgtitle(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}],'Color','white');

        if ~isempty(zidx)
            figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,0.15);
            sgtitle(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}],'Color','white');

            figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,0.15);
            sgtitle(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}],'Color','white');
        end

        if ~isempty(flatXY)
            T2p = T2; T2m = -T2;
            T2p(T2p<Tth) = nan;   % thresholded T-value
            T2m(T2m<Tth) = nan;   % thresholded T-value
            figure; plotNifti3Dflatmap(T2p, maskV, isFullVoxel, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['GLM contrast (plus) of ' sessionName ' : ' contnames{j}]);

            figure; plotNifti3Dflatmap(T2m, maskV, isFullVoxel, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['GLM contrast (minus) of ' sessionName ' : ' contnames{j}]);
        end

        % show Tmax
        tmax = max(V2(:));
        tcnt = length(find(V2>=Tth));
        mrv = nanmean(RSS);
        disp(['Tmax of ' sessionName ' : ' contnames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt) ', mrv=' num2str(mrv)])

        % output 
        V2(isnan(V2)) = 0;
        V2t = V2;
        V2t(V2t>rangePlus(2)) = rangePlus(2);
        V2t(V2t<-rangeMinus(2)) = -rangeMinus(2);
        V2t(-rangeMinus(1)<V2t & V2t < rangePlus(1)) = 0;
        Vts{j} = V2t;
        Vfs{j} = V2;
        Tmaxs(j) = tmax;
        Tcnts(j) = tcnt;
        mrvs(j) = mrv;
    end
end
