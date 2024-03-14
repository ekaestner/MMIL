load([ dta_dir '/' 'performance.mat' ])
load([ dta_dir '/' 'saliency.mat' ])

sal_out_dir = [ out_dir '/' 'Saliency' '/']; ejk_chk_dir(sal_out_dir);

%% Side by Side plots
slc = 1:2:size(sal_dta.x3d_model_original_data.cn_vs_ep,2);

% iR_one = 38;
% iR_two = 42;
% nft_atl_lbl([ iR_one iR_two ],:)
% [~, cor_ind_one ] = max(sum(squeeze(sum(nft_atl==iR_one,1)),2));
% [~, cor_ind_two ] = max(sum(squeeze(sum(nft_atl==iR_two,1)),2));

for cor_ind_one = slc   
    
    thr_dim_hld = rot90(squeeze(mean(sal_dta.x3d_model_original_data.cn_vs_ep(:,cor_ind_one,:,:),4)));
    two_dim_hld = rot90(squeeze(mean(sal_dta.x2d_model_original_data.cn_vs_ep(:,cor_ind_one,:,:),4)));
    
    imAlpha=ones(size(thr_dim_hld));
    imAlpha(isnan(thr_dim_hld))=0;
    
    figure('Visible','off');
    subplot(2,2,1)
    img_one = imagesc(thr_dim_hld,'AlphaData',imAlpha,[-1 4]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    set(gcf, 'InvertHardcopy', 'off')
    subplot(2,2,2)
    img_two = imagesc(two_dim_hld,'AlphaData',imAlpha,[-1 4]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    subplot(2,2,3)
    img_thr = imagesc(thr_dim_hld-two_dim_hld,'AlphaData',imAlpha,[-1 1]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    tightfig;
    print( [ sal_out_dir '/' '3d_vs_2d' '/' 'slice' '_' num2str(cor_ind_one) '.png'] ,'-dpng')
    close all
end

%% Create mega ROIs
% cell2csv([ atl_dir '/' atl_nme '.csv' ],nft_atl_lbl);

% type
mga_typ = 'lateralized'; % bilateral lateralized

ejk_chk_dir([ sal_out_dir '/' 'mega_atlas' '_' mga_typ '/']);

% MegaROIs
mga_roi_csv = mmil_readtext([ atl_dir '/' atl_nme '_mega.csv' ]);
if strcmpi(mga_typ,'bilateral')    
    mga_roi_nme = unique(mga_roi_csv(:,4));
elseif strcmpi(mga_typ,'lateralized')
    mga_roi_nme = unique(strcat(mga_roi_csv(:,4),'_',mga_roi_csv(:,5)));
    mga_roi_csv(:,4) = strcat(mga_roi_csv(:,4),'_',mga_roi_csv(:,5));
end
mga_roi_num = 1001:1000+numel(mga_roi_nme);

mga_roi_atl = nan(size(nft_atl));
for iNR = 1:numel(mga_roi_nme)    
    roi_ind = find(strcmpi(mga_roi_csv(:,4),mga_roi_nme{iNR}));
    for iR = 1:numel(roi_ind)
        mga_roi_atl(nft_atl==mga_roi_csv{roi_ind(iR),1}) = mga_roi_num(iNR);
    end
end

% Plot
cor_slc = 10:5:125;
for iCS = 1:numel(cor_slc)
    figure('Visible','off','Position',[0 0 1920 1080]);
    subplot(1,2,1)
    imagesc(rot90(squeeze(nft_atl(:,cor_slc(iCS),:))))
    subplot(1,2,2)
    imagesc(rot90(squeeze(mga_roi_atl(:,cor_slc(iCS),:))))
    colormap('jet')
    print( [ sal_out_dir '/' 'mega_atlas' '_' mga_typ '/' 'atlas' '_' num2str(cor_slc(iCS)) '.png'] ,'-dpng')
    close all
end

%% Collate ROI data
for iC = 1:numel(cat_nii_nme)
    
    roi_med = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
    roi_men = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
    roi_max = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
        
    for iR = 1:size(roi_med,1)
        
        roi_med{iR,1} = nft_atl_lbl{iR,3};
        roi_men{iR,1} = nft_atl_lbl{iR,3};
        roi_max{iR,1} = nft_atl_lbl{iR,3};
        
        roi_med_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_iqr_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_men_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_std_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_max_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        
        for iSA = 1:numel(sal_grp)
            for iS = 1:size(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iC}),4)
                sal_hld = squeeze(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iC})(:,:,:,iS));
                sal_hld = sal_hld(nft_atl==str2num(nft_atl_lbl{iR,1}));
                
                roi_med_hld(iS,iSA) = nanmedian(sal_hld);
                roi_iqr_hld(iS,iSA) = iqr(sal_hld);                
                roi_men_hld(iS,iSA) = nanmean(sal_hld);
                roi_std_hld(iS,iSA) = nanstd(sal_hld);                
                roi_max_hld(iS,iSA) = nanmax(sal_hld);
                
            end
            
            % ToDo: Save median & mean for 
            roi_med{iR,iSA+1} = [ num2str(roundsd(nanmedian(roi_med_hld(iS,iSA)),2))];% ' (' num2str(roundsd(nanmean(roi_iqr_hld(iS,iSA)),2)) ')' ];
            roi_men{iR,iSA+1} = [ num2str(roundsd(nanmean(roi_men_hld(iS,iSA)),2))];%   ' (' num2str(roundsd(nanmean(roi_std_hld(iS,iSA)),2)) ')' ];
            roi_max{iR,iSA+1} = [ num2str(roundsd(nanmean(roi_max_hld(iS,iSA)),2)) ];
        end
        
        [ ~, roi_med{iR,end-1} ] = ttest2(roi_med_hld(:,1),roi_med_hld(:,3));
        roi_med{iR,end}          = signrank(roi_med_hld(:,1),roi_med_hld(:,3));
        
        [ ~, roi_men{iR,end-1} ] = ttest2(roi_men_hld(:,1),roi_men_hld(:,3));
        roi_men{iR,end}          = signrank(roi_men_hld(:,1),roi_men_hld(:,3));
        
        [ ~, roi_max{iR,end-1} ] = ttest2(roi_max_hld(:,1),roi_max_hld(:,3));
        roi_max{iR,end}          = signrank(roi_max_hld(:,1),roi_max_hld(:,3));       
        
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).median = roi_med_hld;
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).mean   = roi_men_hld;
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).max    = roi_max_hld;
        
    end

    cell2csv( [ sal_out_dir '/' 'saliency' '_' cat_nii_nme{iC} '_median.csv'], [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_med ])
    cell2csv( [ sal_out_dir '/' 'saliency' '_' cat_nii_nme{iC} '_mean.csv'],   [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_men ])
    cell2csv( [ sal_out_dir '/' 'saliency' '_' cat_nii_nme{iC} '_max.csv'],    [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_max ])
end

%% Plots
mes_int = { 'median' 'mean' 'max' };
xdt_val = [ 1 3 5];
xdt_adj = [-.3 0 .3];

ttl_hld = strcat(cat_nii_nme,' VS ');
ttl_hld = cat(2,ttl_hld{:});
ttl_hld = strrep(ttl_hld,'VS','VS ');
ttl_hld = ttl_hld(1:end-4);

for iM = 1:numel(mes_int)
    plt_out_dir = [ sal_out_dir '/' mes_int{iM} '/']; ejk_chk_dir(plt_out_dir);
    
    for iR = 1:size(roi_med,1)
        
        pcfg = [];
        
        cnt = 1;
        for iC = 1:numel(cat_nii_nme)
            for iS = 1:numel(sal_grp)
                pcfg.ydt{cnt} = plt_hld.(roi_med{iR,1}).(cat_nii_nme{iC}).(mes_int{iM})(:,iS);
                pcfg.xdt{cnt} = xdt_val(iC) + xdt_adj(iS);
                
                pcfg.fce_col{cnt}     = mdl_int_col{strcmpi(mdl_int_nme,sal_nme{iS})};
                pcfg.box_plt_col{cnt} = mdl_int_col{strcmpi(mdl_int_nme,sal_nme{iS})};
                
                pcfg.xlb{cnt} = sal_grp{iS};
                
                cnt = cnt + 1;
            end
        end
                
        pcfg.edg_col = repmat({[0 0 0]},1,numel(pcfg.xdt));
        pcfg.box_plt = ones(1,numel(pcfg.xdt));
        
        pcfg.xlm = [ 0.5 max(cat(1,pcfg.xdt{:}))+0.5 ];
        pcfg.ylb = {[ roi_med{iR} ' : ' upper(mes_int{iM}) ]};
        
        pcfg.ttl = ttl_hld;
        
        pcfg.mkr_sze = repmat(40,1,numel(pcfg.xdt));
        pcfg.aph_val = 0.80;
        pcfg.jtr_wdt = 0.15;
        pcfg.box_wdt = 0.20;
        
        pcfg.out_dir = plt_out_dir;
        pcfg.out_nme = [ 'saliency' '_' roi_med{iR,1} '_' mes_int{iM} ];
        
        try ejk_scatter(pcfg); catch; end
        
    end
end

%% Mega-ROIs
nft_atl_mga_lbl = [ num2cell(mga_roi_num') cell(numel(mga_roi_nme),1) mga_roi_nme];

for iC = 1:numel(cat_nii_nme)
    
    roi_med = cell(size(nft_atl_mga_lbl,1),numel(sal_grp)+3);
    roi_men = cell(size(nft_atl_mga_lbl,1),numel(sal_grp)+3);
    roi_max = cell(size(nft_atl_mga_lbl,1),numel(sal_grp)+3);
        
    for iR = 1:size(roi_med,1)
        
        roi_med{iR,1} = nft_atl_mga_lbl{iR,3};
        roi_men{iR,1} = nft_atl_mga_lbl{iR,3};
        roi_max{iR,1} = nft_atl_mga_lbl{iR,3};
        
        roi_med_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_iqr_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_men_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_std_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        roi_max_hld = nan(size(sal_dta.(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
        
        for iSA = 1:numel(sal_grp)
            for iS = 1:size(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iC}),4)
                sal_hld = squeeze(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iC})(:,:,:,iS));
                sal_hld = sal_hld(mga_roi_atl==nft_atl_mga_lbl{iR,1});
                
                roi_med_hld(iS,iSA) = nanmedian(sal_hld);
                roi_iqr_hld(iS,iSA) = iqr(sal_hld);                
                roi_men_hld(iS,iSA) = nanmean(sal_hld);
                roi_std_hld(iS,iSA) = nanstd(sal_hld);                
                roi_max_hld(iS,iSA) = nanmax(sal_hld);
                
            end
            
            % ToDo: Save median & mean for 
            roi_med{iR,iSA+1} = [ num2str(roundsd(nanmedian(roi_med_hld(iS,iSA)),2))];% ' (' num2str(roundsd(nanmean(roi_iqr_hld(iS,iSA)),2)) ')' ];
            roi_men{iR,iSA+1} = [ num2str(roundsd(nanmean(roi_men_hld(iS,iSA)),2))];%   ' (' num2str(roundsd(nanmean(roi_std_hld(iS,iSA)),2)) ')' ];
            roi_max{iR,iSA+1} = [ num2str(roundsd(nanmean(roi_max_hld(iS,iSA)),2)) ];
        end
        
        [ ~, roi_med{iR,end-1} ] = ttest2(roi_med_hld(:,1),roi_med_hld(:,3));
        roi_med{iR,end}          = signrank(roi_med_hld(:,1),roi_med_hld(:,3));
        
        [ ~, roi_men{iR,end-1} ] = ttest2(roi_men_hld(:,1),roi_men_hld(:,3));
        roi_men{iR,end}          = signrank(roi_men_hld(:,1),roi_men_hld(:,3));
        
        [ ~, roi_max{iR,end-1} ] = ttest2(roi_max_hld(:,1),roi_max_hld(:,3));
        roi_max{iR,end}          = signrank(roi_max_hld(:,1),roi_max_hld(:,3));       
        
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).median = roi_med_hld;
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).mean   = roi_men_hld;
        plt_hld.(roi_med{iR}).(cat_nii_nme{iC}).max    = roi_max_hld;
        
    end

    cell2csv( [ sal_out_dir '/' 'mega' '_' mga_typ '_saliency' '_' cat_nii_nme{iC} '_median.csv'], [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_med ])
    cell2csv( [ sal_out_dir '/' 'mega' '_' mga_typ '_saliency' '_' cat_nii_nme{iC} '_mean.csv'],   [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_men ])
    cell2csv( [ sal_out_dir '/' 'mega' '_' mga_typ '_saliency' '_' cat_nii_nme{iC} '_max.csv'],    [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_max ])
end

%% Mega-ROIs Plots
mes_int = { 'median' 'mean' 'max' };
xdt_val = [ 1 3 5];
xdt_adj = [-.3 0 .3];

ttl_hld = strcat(cat_nii_nme,' VS ');
ttl_hld = cat(2,ttl_hld{:});
ttl_hld = strrep(ttl_hld,'VS','VS ');
ttl_hld = ttl_hld(1:end-4);

for iM = 1:numel(mes_int)
    plt_out_dir = [ sal_out_dir '/' 'mega' '_' mga_typ '_' mes_int{iM} '/']; ejk_chk_dir(plt_out_dir);
    
    for iR = 1:size(roi_med,1)
        
        pcfg = [];
        
        cnt = 1;
        for iC = 1:numel(cat_nii_nme)
            for iS = 1:numel(sal_grp)
                pcfg.ydt{cnt} = plt_hld.(roi_med{iR,1}).(cat_nii_nme{iC}).(mes_int{iM})(:,iS);
                pcfg.xdt{cnt} = xdt_val(iC) + xdt_adj(iS);
                
                pcfg.fce_col{cnt}     = mdl_int_col{strcmpi(mdl_int_nme,sal_nme{iS})};
                pcfg.box_plt_col{cnt} = mdl_int_col{strcmpi(mdl_int_nme,sal_nme{iS})};
                
                pcfg.xlb{cnt} = sal_grp{iS};
                
                cnt = cnt + 1;
            end
        end
                
        pcfg.edg_col = repmat({[0 0 0]},1,numel(pcfg.xdt));
        pcfg.box_plt = ones(1,numel(pcfg.xdt));
        
        pcfg.xlm = [ 0.5 max(cat(1,pcfg.xdt{:}))+0.5 ];
        pcfg.ylb = {[ roi_med{iR} ' : ' upper(mes_int{iM}) ]};
        
        pcfg.ttl = ttl_hld;
        
        pcfg.mkr_sze = repmat(40,1,numel(pcfg.xdt));
        pcfg.aph_val = 0.80;
        pcfg.jtr_wdt = 0.15;
        pcfg.box_wdt = 0.20;
        
        pcfg.out_dir = plt_out_dir;
        pcfg.out_nme = [ 'saliency' '_' roi_med{iR,1} '_' mes_int{iM} ];
        
        ejk_scatter(pcfg)
        
    end
end

