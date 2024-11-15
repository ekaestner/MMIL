load([ dta_dir '/' 'performance.mat' ])
load([ dta_dir '/' 'saliency.mat' ])

sal_out_dir = [ out_dir '/' 'Saliency_v2' '/']; ejk_chk_dir(sal_out_dir);

%% Normalize
% All data
for iN = 1:numel(nor_dta_nme)
    
    % Put together data structure
    for iZ = 1:numel(sal_grp)
        for iCT = 1:numel(cat_nii_nme)
            
            sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}) = sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT});
            
            for iF = 1:size(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
               
                % Values to play with
                sal_val_hld = squeeze(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
                
                
                
                % Re-assign
                sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF)   = sal_val_hld;
                sal_dta_nor.(nor_dta_nme{iN}).remove.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld;
                
            end
        end
    end
end

% Original Data
sal_dta_nor.org_dta = sal_dta;

for iZ = 1:numel(sal_grp)
    for iCT = 1:numel(cat_nii_nme)
        for iF = 1:size(sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
            
            % Values to play with
            sal_val_hld = squeeze(sal_dta_nor.org_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
            
            % Remove either middle or low values
            if low_adj
                if cat_typ(iCT)==0
                    spl_pct     = (100-low_pct)/2;
                    low_zsc_pct = 0+spl_pct;
                    hgh_zsc_pct = 100-spl_pct;
                    
                    low_pct_hld = prctile(sal_val_hld(:),low_zsc_pct);
                    hgh_pct_hld = prctile(sal_val_hld(:),hgh_zsc_pct);
                    
                    sal_val_hld( (sal_val_hld>low_pct_hld) & (sal_val_hld<hgh_pct_hld) ) = NaN;
                    
                elseif cat_typ(iCT)==1
                    
                    pct_hld = prctile(sal_val_hld(:),low_pct);
                    sal_val_hld(sal_val_hld<pct_hld) = NaN;
                    
                end
            end
            
            % Re-assign
            sal_dta_nor.org_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld;
            
        end
    end
end

% z-score
sal_dta_nor.zsc_dta = sal_dta;

for iZ = 1:numel(sal_grp)
    for iCT = 1:numel(cat_nii_nme)
        for iF = 1:size(sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
            
        % Values to play with
        sal_val_hld = squeeze(sal_dta_nor.zsc_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
    
        % Z-score
        sal_val_men = nanmean(sal_val_hld(:));
        sal_val_std = nanstd(sal_val_hld(:));
        sal_val_hld = (sal_val_hld-sal_val_men) ./ sal_val_std;
        
        % Remove either middle or low values
        if low_adj
            if cat_typ(iCT)==0
                spl_pct     = (100-low_pct)/2;
                low_zsc_pct = 0+spl_pct;
                hgh_zsc_pct = 100-spl_pct;
                
                low_pct_hld = prctile(sal_val_hld(:),low_zsc_pct);
                hgh_pct_hld = prctile(sal_val_hld(:),hgh_zsc_pct);
                
                sal_val_hld( (sal_val_hld>low_pct_hld) & (sal_val_hld<hgh_pct_hld) ) = NaN;
                
            elseif cat_typ(iCT)==1
                
                pct_hld = prctile(sal_val_hld(:),low_pct);
                sal_val_hld(sal_val_hld<pct_hld) = NaN;
                
            end
        end
        
        % Re-assign
        sal_dta_nor.zsc_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld;
        
        end                   
    end
end


% min-max
sal_dta_nor.min_max_dta = sal_dta;

for iZ = 1:numel(sal_grp)
    for iCT = 1:numel(cat_nii_nme)
        for iF = 1:size(sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
            
            % Values to play with
            sal_val_hld = squeeze(sal_dta_nor.min_max_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
            
            % Minmax
            sal_val_max = max(sal_val_hld(:));
            sal_val_min = min(sal_val_hld(:));
            sal_val_hld = (sal_val_hld-sal_val_min) ./ (sal_val_max - sal_val_min);
            
            % Remove either middle or low values
            if low_adj
            if cat_typ(iCT)==0
                spl_pct     = (100-low_pct)/2;
                low_zsc_pct = 0+spl_pct;
                hgh_zsc_pct = 100-spl_pct;
                
                low_pct_hld = prctile(sal_val_hld(:),low_zsc_pct);
                hgh_pct_hld = prctile(sal_val_hld(:),hgh_zsc_pct);
                
                sal_val_hld( (sal_val_hld>low_pct_hld) & (sal_val_hld<hgh_pct_hld) ) = NaN;
                
            elseif cat_typ(iCT)==1
                
                pct_hld = prctile(sal_val_hld(:),low_pct);
                sal_val_hld(sal_val_hld<pct_hld) = NaN;
                
            end
            end
            
            % Re-assign
            sal_dta_nor.min_max_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld;
            
        end
    end
end

% Make histograms
for iZ = 1:numel(sal_grp)
    for iCT = 1:numel(cat_nii_nme)
        
        figure('Visible','off')
        subplot(numel(nor_dta_nme),1,1); hold on;
        title(mmil_spec_char([sal_grp{iZ} ' ' cat_nii_nme{iCT}],{'_'},{' '}))
        for iNO = 1:numel(nor_dta_nme)
            subplot(numel(nor_dta_nme),1,iNO); hold on;
            for iF = 1:size(sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
                sal_val_hld = squeeze(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
                histogram(sal_val_hld(:),200);
            end
            ylabel(mmil_spec_char(nor_dta_nme{iNO},{'_'},{' '}))
        end
        tightfig();        
        print([ sal_out_dir '/' 'histograms' '/' sal_grp{iZ} '_' cat_nii_nme{iCT} '.png'],'-dpng')
        
    end
end

%% Side by Side plots
out_nme = 'Reihaneh_test';
cmp_nme = { 'x3D' 'x3D_1st_method' };
iNO     = [ 3     3 ];
typ_nme = 'ep';
ejk_chk_dir([ sal_out_dir '/'  'side_by_side' '/' out_nme '/' nor_dta_nme{iNO(1)} '_' nor_dta_nme{iNO(2)} '/' ]);

slc = 1:2:size(sal_dta_nor.(nor_dta_nme{1}).(cmp_nme{1}).(typ_nme),2);

ejk_chk_dir([ sal_out_dir '/' 'side_by_side' '/' out_nme '/' nor_dta_nme{iNO(1)}]);
ylm_low_one = prctile(sal_dta_nor.(nor_dta_nme{iNO(1)}).(cmp_nme{1}).(typ_nme)(:),10);
ylm_hgh_one = prctile(sal_dta_nor.(nor_dta_nme{iNO(1)}).(cmp_nme{1}).(typ_nme)(:),90);

ejk_chk_dir([ sal_out_dir '/' 'side_by_side' '/' out_nme '/' nor_dta_nme{iNO(2)}]);
ylm_low_two = prctile(sal_dta_nor.(nor_dta_nme{iNO(2)}).(cmp_nme{2}).(typ_nme)(:),10);
ylm_hgh_two = prctile(sal_dta_nor.(nor_dta_nme{iNO(2)}).(cmp_nme{2}).(typ_nme)(:),90);

sub_hld = mean(sal_dta_nor.(nor_dta_nme{iNO(1)}).(cmp_nme{1}).(typ_nme),4) - mean(sal_dta_nor.(nor_dta_nme{iNO(2)}).(cmp_nme{2}).(typ_nme),4);
ylm_dff_low = prctile(sub_hld(:),10);
ylm_dff_hgh = prctile(sub_hld(:),90);

for cor_ind_one = slc
    
    thr_dim_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_dta_nme{iNO(1)}).(cmp_nme{1}).(typ_nme)(:,cor_ind_one,:,:),4)));
    two_dim_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_dta_nme{iNO(2)}).(cmp_nme{2}).(typ_nme)(:,cor_ind_one,:,:),4)));
    
    imAlpha_one=ones(size(thr_dim_hld));
    imAlpha_one(isnan(thr_dim_hld))=0;
    
    imAlpha_two=ones(size(two_dim_hld));
    imAlpha_two(isnan(two_dim_hld))=0;
    
    imAlpha_dff=zeros(size(thr_dim_hld));
    imAlpha_dff( ~isnan(thr_dim_hld) & ~isnan(two_dim_hld) )=1;
    
    figure('Visible','off');
    subplot(2,2,1)
    img_one = imagesc(thr_dim_hld,'AlphaData',imAlpha_one,[ ylm_low_one ylm_hgh_one ]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    set(gcf, 'InvertHardcopy', 'off')
    subplot(2,2,2)
    img_two = imagesc(two_dim_hld,'AlphaData',imAlpha_two,[ ylm_low_two ylm_hgh_two ]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    subplot(2,2,3)
    img_thr = imagesc(thr_dim_hld-two_dim_hld,'AlphaData',imAlpha_dff,[ ylm_dff_low ylm_dff_hgh ]);
    set(gca,'color',rgb('bluish grey')+0.40);
    colorbar
    axis off;
    %         tightfig;
    print( [ sal_out_dir '/'  'side_by_side' '/' out_nme '/' nor_dta_nme{iNO(1)} '_' nor_dta_nme{iNO(2)} '/' 'slice' '_' num2str(cor_ind_one) '.png'] ,'-dpng')
    close all
end

%% Collate ROI data
for iNO = 1:numel(nor_dta_nme)
    for iC = 1:numel(cat_nii_nme)
        
        ejk_chk_dir([ sal_out_dir '/' 'median' '/' nor_dta_nme{iNO} '/' ]);
        ejk_chk_dir([ sal_out_dir '/' 'mean'   '/' nor_dta_nme{iNO} '/' ]);
        ejk_chk_dir([ sal_out_dir '/' 'max'    '/' nor_dta_nme{iNO} '/' ]);
        
        roi_med = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
        roi_men = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
        roi_max = cell(size(nft_atl_lbl,1),numel(sal_grp)+3);
        
        for iR = 1:size(roi_med,1)
            
            roi_med{iR,1} = nft_atl_lbl{iR,3};
            roi_men{iR,1} = nft_atl_lbl{iR,3};
            roi_max{iR,1} = nft_atl_lbl{iR,3};
            
            roi_med_hld = nan(size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
            roi_iqr_hld = nan(size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
            roi_men_hld = nan(size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
            roi_std_hld = nan(size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
            roi_max_hld = nan(size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{1}]).(cat_nii_nme{iC}),4),numel(sal_grp));
            
            for iSA = 1:numel(sal_grp)
                for iS = 1:size(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{iSA}]).(cat_nii_nme{iC}),4)
                    sal_hld = squeeze(sal_dta_nor.(nor_dta_nme{iNO}).(['x' sal_grp{iSA}]).(cat_nii_nme{iC})(:,:,:,iS));
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
            
%             [ ~, roi_med{iR,end-1} ] = ttest2(roi_med_hld(:,1),roi_med_hld(:,3));
%             roi_med{iR,end}          = signrank(roi_med_hld(:,1),roi_med_hld(:,3));
%             
%             [ ~, roi_men{iR,end-1} ] = ttest2(roi_men_hld(:,1),roi_men_hld(:,3));
%             roi_men{iR,end}          = signrank(roi_men_hld(:,1),roi_men_hld(:,3));
%             
%             [ ~, roi_max{iR,end-1} ] = ttest2(roi_max_hld(:,1),roi_max_hld(:,3));
%             roi_max{iR,end}          = signrank(roi_max_hld(:,1),roi_max_hld(:,3));
            
            plt_hld.(nor_dta_nme{iNO}).(roi_med{iR}).(cat_nii_nme{iC}).median = roi_med_hld;
            plt_hld.(nor_dta_nme{iNO}).(roi_med{iR}).(cat_nii_nme{iC}).mean   = roi_men_hld;
            plt_hld.(nor_dta_nme{iNO}).(roi_med{iR}).(cat_nii_nme{iC}).max    = roi_max_hld;
            
        end
        
        cell2csv( [ sal_out_dir '/' 'median' '/' nor_dta_nme{iNO} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_median.csv'], [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_med ])
        cell2csv( [ sal_out_dir '/' 'mean'   '/' nor_dta_nme{iNO} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_mean.csv'],   [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_men ])
        cell2csv( [ sal_out_dir '/' 'max'    '/' nor_dta_nme{iNO} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_max.csv'],    [ 'ROI' sal_grp 'ttest p' 'wilcoxon p' ; roi_max ])
    end
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
    for iNO = 1:numel(nor_dta_nme)
        plt_out_dir = [ sal_out_dir '/' mes_int{iM} '/' nor_dta_nme{iNO} '/']; ejk_chk_dir(plt_out_dir);
        
        for iR = 1:size(roi_med,1)
            
            pcfg = [];
            
            cnt = 1;
            for iC = 1:numel(cat_nii_nme)
                for iS = 1:numel(sal_grp)
                    pcfg.ydt{cnt} = plt_hld.(nor_dta_nme{iNO}).(roi_med{iR,1}).(cat_nii_nme{iC}).(mes_int{iM})(:,iS);
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
            pcfg.out_nme = [ 'saliency' '_' roi_med{iR,1} '_' mes_int{iM} '_' nor_dta_nme{iNO} ];
            
            try ejk_scatter(pcfg); catch; end
            
        end
    end
end

