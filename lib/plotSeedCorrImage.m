%%
% Plot multi seed analysis results
% returns
%   cells of thresholded T or Z value 3D voxels (seeds are source, plus side) (CVsp), cells of thresholded T or Z value 3D voxels (seeds are source, minus side) (CVsm)
%   cells of thresholded T or Z value 3D voxels (seeds are target, plus side) (CVtp), cells of thresholded T or Z value 3D voxels (seeds are target, minus side) (CVtm)
% input:
%  tgNames      cells of target region names
%  Tgs          cells of target region vectors (node x 1)
%  T2           T-value or Z-value matrix (node x node)
%  df           degree of freedom (df)
%  Pth          P-value threshold for T or Z value matrix
%  maskV        3D mask volume
%  backV        3D background volume 
%  isFullVoxel  full voxel atlas or not
%  isRtoL       X axis is right to left (default: true)
%  sessionName  session name for title (optional)
%  corrMeth     family-wise error rate (FWER) correction method ('none','bonf','sidak','holm-bonf','holm-sidak')
%  rangePlus    plus T-value range ([min max]) (optional)
%  rangeMinus   minus T-value range ([min max]) (optional)
%  cmap         color map for 3D plot (default: hot)
%  zidx         show z-index slices (default: [])
%  flatXY       show functional flat map (default: [])

function [CVsp, CVsm, CVtp, CVtm, Tmaxs, Tcnts, mrvs] = plotSeedCorrImage(tgNames, Tgs, df, T2, Pth, maskV, backV, isFullVoxel, isRtoL, sessionName, corrMeth, rangePlus, rangeMinus, cmap, zidx, flatXY)
    if nargin < 16, flatXY = []; end
    if nargin < 15, zidx = []; end
    if nargin < 14, cmap = hot; end
    if nargin < 13, rangeMinus = [3 9]; end
    if nargin < 12, rangePlus = [3 10]; end
    if nargin < 11, corrMeth = 'none'; end
    if nargin < 10, sessionName = ''; end
    if nargin < 9, isRtoL = true; end

    if isinf(df), Tstr='Z'; else Tstr='T'; end % T-value or Z-value

    Tmaxs = [];
    Tcnts = [];
    mrvs = [];
    CVsp = {}; CVsm = {};
    CVtp = {}; CVtm = {};

    for j=1:length(tgNames)
        tg = Tgs{j};

        % 'target' is source contrast 
        if size(T2,2) == size(tg,1)
            T2c = T2 .* tg';
        else
            T2c = T2 .* [tg; 0]';
        end
        T2p = max(T2c,[],2);
        T2m = max(-T2c,[],2);

        % T-value threshold (Bonferroni correction / Šidák correction)
        if strcmp(corrMeth,'bonf') || strcmp(corrMeth,'sidak')
            % only masked ROI is used
            if size(T2,1) == size(T2,2)
                m = size(T2,1) * sum(tg~=0);
            else
                m = size(T2,1); % basically, column ROIs are independent.
            end
            if strcmp(corrMeth,'bonf')
                BPth = Pth / m;
            else
                % Šidák correction - almost same as bonferroni
                BPth = 1 - power(1 - Pth, 1/m);
            end
            if isinf(df)
                Tth = abs(norminv(BPth));
            else
                Tth = abs(tinv(BPth,df));
            end
            disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), ' Tstr '-value=' num2str(Tth)])
        % Holm–Bonferroni method
        elseif strcmp(corrMeth,'holm-bonf') || strcmp(corrMeth,'holm-sidak')
            % check symmetry or not
            tg2 = tg;
            tg2(tg2==0) = nan;
            if size(T2,2) == size(tg,1)
                T2c2 = T2 .* tg2';
            else
                T2c2 = T2 .* [tg2; 0]';
            end
            T2c2(isnan(T2c2)) = []; % remove NaN
            T2s = sort(T2c2(:),'descend'); % full matrix T-test
            for k=1:length(T2s)
                if strcmp(corrMeth,'holm-bonf')
                    BPth = Pth / (length(T2s) + 1 - k);
                else
                    BPth = 1 - power(1 - Pth, 1/(length(T2s) + 1 - k));
                end
                if isinf(df)
                    Tth = abs(norminv(BPth));
                else
                    Tth = abs(tinv(BPth,df));
                end
                if T2s(k) <= Tth
                    break;
                end
            end
            disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), ' Tstr '-value=' num2str(Tth)])
        else
            if isinf(df)
                Tth = abs(norminv(Pth));
            else
                Tth = abs(tinv(Pth,df));
            end
            disp(['P-value=' num2str(Pth) ', ' Tstr '-value=' num2str(Tth)])
        end
        % set auto range
        if isnan(rangePlus(1))
            rangePlus(1) = Tth;
        end
        if isnan(rangeMinus(1))
            rangeMinus(1) = Tth;
        end

        T2p(T2p<Tth) = nan;   % thresholded T-value
        T2m(T2m<Tth) = nan;   % thresholded T-value
        if isFullVoxel == 1
            aIdx = find(maskV(:) > 0);
            Vp = single(maskV);
            Vp(:) = nan;
            Vm = Vp;
            Vp(aIdx) = T2p;
            Vm(aIdx) = T2m;
        else
            Vp = getNifti4DFromRoiTS(T2p, maskV);
            Vm = getNifti4DFromRoiTS(T2m, maskV);
        end

        % show figure
        figure; plotNifti3DAxes(Vp,'max',rangePlus,cmap,backV,gray,isRtoL);
        sgtitle(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
        figure; plotNifti3DAxes(Vm,'max',rangeMinus,cmap,backV,gray,isRtoL);
        sgtitle(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');

        if ~isempty(zidx)
            figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');

            figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
        end

        if ~isempty(flatXY)
            figure; plotNifti3Dflatmap(T2p, maskV, isFullVoxel, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}]);

            figure; plotNifti3Dflatmap(T2m, maskV, isFullVoxel, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}]);
        end
        
        % show Tmax
        tmax = max(Vp(:));
        tcnt = length(find(Vp>=Tth));
        mrv = NaN;
        disp([Tstr 'max (source) of ' sessionName ' : ' tgNames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt) ', mrv=' num2str(mrv)])

        Vp(Vp > rangePlus(2)) = rangePlus(2); Vp = Vp - rangePlus(1); Vp(Vp<0) = 0;
        Vm(Vm > rangeMinus(2)) = rangeMinus(2); Vm = Vm - rangeMinus(1); Vm(Vm<0) = 0;
        CVsp{j} = Vp;
        CVsm{j} = Vm;

        % 'target' is target voxels
        if size(T2,1) == size(tg,1)
            T2c = T2 .* tg;
            T2p = max(T2c,[],1);
            T2m = max(-T2c,[],1);
            T2p = T2p(1:end-1)';   % remove intercept
            T2m = T2m(1:end-1)';   % remove intercept
    
            T2p(T2p<Tth) = nan;   % thresholded T-value
            T2m(T2m<Tth) = nan;   % thresholded T-value
            if isFullVoxel == 1
                aIdx = find(maskV(:) > 0);
                Vp = single(maskV);
                Vp(:) = nan;
                Vm = Vp;
                Vp(aIdx) = T2p;
                Vm(aIdx) = T2m;
            else
                Vp = getNifti4DFromRoiTS(T2p, maskV);
                Vm = getNifti4DFromRoiTS(T2m, maskV);
            end
    
            % show figure
            figure; plotNifti3DAxes(Vp,'max',rangePlus,cmap,backV,gray,isRtoL);
            sgtitle(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
            figure; plotNifti3DAxes(Vm,'max',rangeMinus,cmap,backV,gray,isRtoL);
            sgtitle(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
    
            if ~isempty(zidx)
                figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,isRtoL,0.15);
                sgtitle(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
    
                figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,isRtoL,0.15);
                sgtitle(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
            end

            if ~isempty(flatXY)
                figure; plotNifti3Dflatmap(T2p, maskV, isFullVoxel, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
                title(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}]);
    
                figure; plotNifti3Dflatmap(T2m, maskV, isFullVoxel, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
                title(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}]);
            end

            % show Tmax
            tmax = max(Vp(:));
            tcnt = length(find(Vp>=Tth));
            mrv = NaN;
            disp([Tstr 'max (target) of ' sessionName ' : ' tgNames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt) ', mrv=' num2str(mrv)])

            Vp(Vp > rangePlus(2)) = rangePlus(2); Vp = Vp - rangePlus(1); Vp(Vp<0) = 0;
            Vm(Vm > rangeMinus(2)) = rangeMinus(2); Vm = Vm - rangeMinus(1); Vm(Vm<0) = 0;
            CVtp{j} = Vp;
            CVtm{j} = Vm;
        end

        % output 
        Tmaxs(j) = tmax;
        Tcnts(j) = tcnt;
        mrvs(j) = mrv;
    end
end
