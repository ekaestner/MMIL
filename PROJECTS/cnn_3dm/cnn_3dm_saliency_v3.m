load([ dta_dir '/' 'performance.mat' ])
load([ dta_dir '/' 'saliency.mat' ])

sal_out_dir = [ out_dir '/' 'Saliency_v2' '/']; ejk_chk_dir(sal_out_dir);

cat_nii_nme_int = find(strcmpi(cat_nii_nme,'ep'));

%% Normalize
% All data
for iN = 1:numel(nor_dta_nme)
    
    % Put together data structure
    for iZ = 1:numel(sal_grp)
        for iCT = cat_nii_nme_int %1:numel(cat_nii_nme)
            
            sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}) = sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT});
            
            if size(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4) > 0
            for iF = 1:size(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
               
                % Values to play with
                sal_val_hld = squeeze(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
                
                % Normalize
                if strcmpi(nor_dta_nme{iN}, 'zsc_dta')
                    sal_val_men = nanmean(sal_val_hld(:));
                    sal_val_std = nanstd(sal_val_hld(:));
                    sal_val_hld = (sal_val_hld-sal_val_men) ./ sal_val_std;
                elseif strcmpi(nor_dta_nme{iN}, 'min_max_dta')
                    sal_val_max = max(sal_val_hld(:));
                    sal_val_min = min(sal_val_hld(:));
                    sal_val_hld = (sal_val_hld-sal_val_min) ./ (sal_val_max - sal_val_min);
                end
                    
                % Remove values
                sal_val_hld_rmv     = sal_val_hld;
                sal_val_hld_rmv_mid = sal_val_hld;
                
                % Remove middle
                spl_pct     = (100-low_pct)/2;
                low_zsc_pct = 0+spl_pct;
                hgh_zsc_pct = 100-spl_pct;
                
                low_pct_hld = prctile(sal_val_hld(:),low_zsc_pct);
                hgh_pct_hld = prctile(sal_val_hld(:),hgh_zsc_pct);
                
                sal_val_hld_rmv_mid( (sal_val_hld>low_pct_hld) & (sal_val_hld<hgh_pct_hld) ) = NaN;
                
                % Remove bottom
                pct_hld = prctile(sal_val_hld(:),low_pct);
                sal_val_hld_rmv(sal_val_hld<pct_hld) = NaN;
                                
                % Re-assign
                sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF)   = sal_val_hld;
                sal_dta_nor.(nor_dta_nme{iN}).remove.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld_rmv;
                sal_dta_nor.(nor_dta_nme{iN}).remove_mid.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF) = sal_val_hld_rmv_mid;
                
                clear sal_val_hld sal_val_hld_rmv sal_val_hld_rmv_mid
            
            end
            elseif size(sal_dta_nor.(nor_dta_nme{iN}).keep.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4) == 0
                sal_dta_nor.(nor_dta_nme{iN}).remove.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})     = sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT});
                sal_dta_nor.(nor_dta_nme{iN}).remove_mid.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}) = sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT});
            end
        end
    end
end

%% Make histograms
ejk_chk_dir([ sal_out_dir '/' 'histograms' '/' ])
for iK = 1:numel(kep_dta_nme)
    for iZ = 1:numel(sal_grp)
        for iCT = cat_nii_nme_int%1:numel(cat_nii_nme)
            
            figure('Visible','off')
            subplot(numel(nor_dta_nme),1,1); hold on;
            title(mmil_spec_char([sal_grp{iZ} ' ' cat_nii_nme{iCT}],{'_'},{' '}))
            for iNO = 1:numel(nor_dta_nme)
                subplot(numel(nor_dta_nme),1,iNO); hold on;
                for iF = 1:size(sal_dta.(['x' sal_grp{iZ}]).(cat_nii_nme{iCT}),4)
                    sal_val_hld = squeeze(sal_dta_nor.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(['x' sal_grp{iZ}]).(cat_nii_nme{iCT})(:,:,:,iF));
                    histogram(sal_val_hld(:),200);
                end
                ylabel(mmil_spec_char(nor_dta_nme{iNO},{'_'},{' '}))
            end
            tightfig();
            print([ sal_out_dir '/' 'histograms' '/' sal_grp{iZ} '_' cat_nii_nme{iCT} '_' kep_dta_nme{iK} '.png'],'-dpng')
            
        end
    end
end

%% Side by Side plots
out_nme = { '3D_vs_2D_ep_zsc_keep_set_lim_hot_focus'    };% '3D_vs_2D_dff_zsc' };
cmp_nme = { {'x3D'      'x2D'}    };% {'x3D'         'x2D'} };
typ_nme = { {'ep'       'ep' }    };% {'cn-ep'       'cn-ep' } };
nor_nme = { {'zsc_dta'  'zsc_dta'}};% {'min_max_dta' 'org_dta'} };
kep_nme = { {'keep'     'keep'}   };% {'keep'        'keep'} };
col_nme = { {'copper'   'copper' } };
grp_lim = { {[0 1.5] [0 1.5]} };
dff_lim = { [-2 2] };
hot_col = hot;
col_col = cool; col_col =  [ col_col(1:127,:) ; repmat(col_col(128,:),50,1) ; col_col(129:end,:) ];
for iP = 1:numel(out_nme)
    
    ejk_chk_dir([ sal_out_dir '/'  'side_by_side' '/' out_nme{iP} '/' ]);
    
    txt_out(1,1:2) = { cmp_nme{iP}{1} cmp_nme{iP}{2}};
    txt_out(2,1:2) = { typ_nme{iP}{1} typ_nme{iP}{2}};
    txt_out(3,1:2) = { nor_nme{iP}{1} nor_nme{iP}{2}};
    txt_out(4,1:2) = { kep_nme{iP}{1} kep_nme{iP}{2}};
    cell2csv([ sal_out_dir '/'  'side_by_side' '/'  out_nme{iP} '/' '00_values_.csv'],txt_out)
    
    slc = 1:2:size(sal_dta_nor.(nor_nme{iP}{1}).(kep_nme{iP}{1}).(cmp_nme{iP}{1}).(typ_nme{iP}{1}),2);
    
    if isempty(grp_lim{iP}{1})
        ylm_low_one = prctile(sal_dta_nor.(nor_nme{iP}{1}).(kep_nme{iP}{1}).(cmp_nme{iP}{1}).(typ_nme{iP}{1})(:),10);
        ylm_hgh_one = prctile(sal_dta_nor.(nor_nme{iP}{1}).(kep_nme{iP}{1}).(cmp_nme{iP}{1}).(typ_nme{iP}{1})(:),90);
    else
        ylm_low_one = grp_lim{iP}{1}(1);
        ylm_hgh_one = grp_lim{iP}{1}(2);
    end
    
    if isempty(grp_lim{iP})
        ylm_low_two = prctile(sal_dta_nor.(nor_nme{iP}{2}).(kep_nme{iP}{2}).(cmp_nme{iP}{2}).(typ_nme{iP}{2})(:),10);
        ylm_hgh_two = prctile(sal_dta_nor.(nor_nme{iP}{2}).(kep_nme{iP}{2}).(cmp_nme{iP}{2}).(typ_nme{iP}{2})(:),90);
    else
        ylm_low_two = grp_lim{iP}{1}(1);
        ylm_hgh_two = grp_lim{iP}{1}(2);
    end
    
    if isempty(dff_lim{iP})
        sub_hld = mean(sal_dta_nor.(nor_nme{iP}{1}).(kep_nme{iP}{1}).(cmp_nme{iP}{1}).(typ_nme{iP}{1}),4) - mean(sal_dta_nor.(nor_nme{iP}{2}).(kep_nme{iP}{2}).(cmp_nme{iP}{2}).(typ_nme{iP}{2}),4);
        ylm_dff_low = prctile(sub_hld(:),10);
        ylm_dff_hgh = prctile(sub_hld(:),90);
    else
        ylm_dff_low = dff_lim{iP}(1);
        ylm_dff_hgh = dff_lim{iP}(2);
    end
    
    for cor_ind_one = slc
        
        one_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{iP}{1}).(kep_nme{iP}{1}).(cmp_nme{iP}{1}).(typ_nme{iP}{1})(11:end-9,cor_ind_one,12:end-12,:),4)));
        two_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{iP}{2}).(kep_nme{iP}{2}).(cmp_nme{iP}{2}).(typ_nme{iP}{2})(11:end-9,cor_ind_one,12:end-12,:),4)));
                
        imAlpha_one=ones(size(one_hld));
        imAlpha_one(isnan(one_hld))=0;
        
        imAlpha_two=ones(size(two_hld));
        imAlpha_two(isnan(two_hld))=0;
        
        imAlpha_dff=zeros(size(one_hld));
        imAlpha_dff( ~isnan(one_hld) & ~isnan(two_hld) )=1;
        
        figure('Visible','off');
        
        ax(1) = subplot(3,3,1);
        img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ ylm_low_one ylm_hgh_one ]);
        colormap(ax(1),hot_col(1:end-45,:));%col_nme{iP}{1}
        set(gca,'color',rgb('bluish grey')+0.40);
        axis off;
        set(gcf, 'InvertHardcopy', 'off')
        ax(2) = subplot(3,3,2);
        img_two = imagesc(two_hld,'AlphaData',imAlpha_two,[ ylm_low_two ylm_hgh_two ]);
        colormap(ax(2),hot_col(1:end-45,:)); %col_nme{iP}{2}
        set(gca,'color',rgb('bluish grey')+0.40);
        axis off;
        ax(3) = subplot(3,3,3);
        img_thr = imagesc(one_hld-two_hld,'AlphaData',imAlpha_dff,[ ylm_dff_low ylm_dff_hgh ]);
        colormap(ax(3),col_col)
        set(gca,'color',rgb('bluish grey')+0.40);
        axis off;
        
        ax(4) = subplot(3,3,4);
        img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ ylm_low_one ylm_hgh_one ]);
        colormap(ax(4),hot_col(1:end-45,:));%col_nme{iP}{1}
        set(gca,'color',rgb('bluish grey')+0.40);
        colorbar
        axis off;
        set(gcf, 'InvertHardcopy', 'off')
        ax(5) = subplot(3,3,5);
        img_two = imagesc(two_hld,'AlphaData',imAlpha_two,[ ylm_low_two ylm_hgh_two ]);
        colormap(ax(5),hot_col(1:end-45,:)); %col_nme{iP}{2}
        set(gca,'color',rgb('bluish grey')+0.40);
        colorbar
        axis off;
        ax(6) = subplot(3,3,6);
        img_thr = imagesc(one_hld-two_hld,'AlphaData',imAlpha_dff,[ ylm_dff_low ylm_dff_hgh ]);
        colormap(ax(6),col_col)
        set(gca,'color',rgb('bluish grey')+0.40);
        colorbar
        axis off;
        
        
        tightfig();
        set(gcf','Position',[0 0 1080 1080]);
        print( [ sal_out_dir '/'  'side_by_side' '/'  out_nme{iP} '/' 'slice' '_' num2str(cor_ind_one) '_'  '.png'] ,'-dpng')
        close all
        
    end
end

%% Axial plots
cmp_nme = {'x3D'       'x2D' };% {'x3D'         'x2D'} };
typ_nme = {'ep'        'ep'  };% {'cn-ep'       'cn-ep' } };
nor_nme = { 'zsc_dta'  'zsc_dta' };% {'min_max_dta' 'org_dta'} };
kep_nme = { 'keep'     'keep'    };% {'keep'        'keep'} };
col_nme = { 'hot'   'hot'  };
grp_lim = { [0 1.0] [0 1.0] };
hot_col = hot;

for iP = 1:numel(cmp_nme)

    ejk_chk_dir([ sal_out_dir '/'  'axial' '/' ]);
    
    slc = 18:6:size(sal_dta_nor.(nor_nme{iP}).(kep_nme{iP}).(cmp_nme{iP}).(typ_nme{iP}),3)-18;
    dim = ceil(sqrt(numel(slc)));
    
    figure('Visible','off');
    
    for iSB = 1:numel(slc)
        
        one_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{iP}).(kep_nme{iP}).(cmp_nme{iP}).(typ_nme{iP})(11:end-9,11:end-8,slc(iSB),:),4)));
        
        imAlpha_one=ones(size(one_hld));
        imAlpha_one(isnan(one_hld))=0;
        
        ax(iSB) = subplot(dim,dim,iSB);
        img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ grp_lim{iP}(1) grp_lim{iP}(2) ]);
        colormap(ax(iSB),hot_col(1:end-45,:)); %col_nme{iP}
        set(gca,'color',rgb('white'));
        axis off;
        set(gcf, 'InvertHardcopy', 'off')
        
    end
    
    tightfig();
    set(gcf','Position',[0 0 1080 1080]);
    print( [ sal_out_dir '/'  'axial' '/' cmp_nme{iP} '_' typ_nme{iP} '_' nor_nme{iP} '_' kep_nme{iP}  '_' col_nme{iP} '.png'] ,'-dpng')
    close all
    
end

%% Collate ROI data
cmp_nme = { {'3D'         '2D'} };
cmp_num = numel(cmp_nme);
stt_nme = [ sal_grp cell(1,numel(cmp_nme)*2) ];

for iNO = 1:numel(nor_dta_nme)
    for iK = 1:numel(kep_dta_nme)
        
        ejk_chk_dir([ sal_out_dir '/' 'ROIs' '/' 'median' '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' ]);
        ejk_chk_dir([ sal_out_dir '/' 'ROIs' '/' 'mean'   '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' ]);
        ejk_chk_dir([ sal_out_dir '/' 'ROIs' '/' 'max'    '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' ]);
        
        for iC = cat_nii_nme_int%:numel(cat_nii_nme)
            
            roi_med = cell(size(nft_atl_lbl,1),1+numel(sal_grp)+2*cmp_num);
            roi_men = cell(size(nft_atl_lbl,1),1+numel(sal_grp)+2*cmp_num);
            roi_max = cell(size(nft_atl_lbl,1),1+numel(sal_grp)+2*cmp_num);
            
            for iR = 1:size(roi_med,1)
                
                roi_med{iR,1} = nft_atl_lbl{iR,3};
                roi_men{iR,1} = nft_atl_lbl{iR,3};
                roi_max{iR,1} = nft_atl_lbl{iR,3};                

                sal_dta_tmp = sal_dta_nor.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(['x' sal_grp{1}]).(cat_nii_nme{iC});
                
                roi_med_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
                roi_iqr_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
                roi_men_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
                roi_std_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
                roi_max_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
                
                for iSA = 1:numel(sal_grp)
                    for iS = 1:size(sal_dta_nor.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(['x' sal_grp{iSA}]).(cat_nii_nme{iC}),4)
                        sal_hld = squeeze(sal_dta_nor.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(['x' sal_grp{iSA}]).(cat_nii_nme{iC})(:,:,:,iS));
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
                
                % Tests
                for iT = 1:numel(cmp_nme)
                    
                    ind = 0+((iT-1)*2):1+((iT-1)*2);
                    
                    col_one = strcmpi(sal_grp,cmp_nme{iT}{1});
                    col_two = strcmpi(sal_grp,cmp_nme{iT}{2});
                    
                    stt_nme(1,end-ind(1)) = {[ cmp_nme{iT}{1} '_VS_' cmp_nme{iT}{2} '_ttest' ]};
                    stt_nme(1,end-ind(2)) = {[ cmp_nme{iT}{1} '_VS_' cmp_nme{iT}{2} '_signrank']};
                    
                    [ ~, roi_med{iR,end-ind(1)} ] = ttest2(roi_med_hld(:,col_one),roi_med_hld(:,col_two));
                    roi_med{iR,end-ind(2)}        = signrank(roi_med_hld(:,col_one),roi_med_hld(:,col_two));
                    
                    [ ~, roi_men{iR,end-ind(1)} ] = ttest2(roi_men_hld(:,col_one),roi_men_hld(:,col_two));
                    roi_men{iR,end-ind(2)}        = signrank(roi_men_hld(:,col_one),roi_men_hld(:,col_two));
                    
                    [ ~, roi_max{iR,end-ind(1)} ] = ttest2(roi_max_hld(:,col_one),roi_max_hld(:,col_two));
                    roi_max{iR,end-ind(2)}        = signrank(roi_max_hld(:,col_one),roi_max_hld(:,col_two));
                end
                
                % Hold
                plt_hld.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(roi_med{iR}).(cat_nii_nme{iC}).median = roi_med_hld;
                plt_hld.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(roi_med{iR}).(cat_nii_nme{iC}).mean   = roi_men_hld;
                plt_hld.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(roi_med{iR}).(cat_nii_nme{iC}).max    = roi_max_hld;
                
            end
            
            cell2csv( [ sal_out_dir '/' 'ROIs' '/' 'median' '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_median.csv'], [ 'ROI' stt_nme ; roi_med ])
            cell2csv( [ sal_out_dir '/' 'ROIs' '/' 'mean' '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_mean.csv'],     [ 'ROI' stt_nme ; roi_men ])
            cell2csv( [ sal_out_dir '/' 'ROIs' '/' 'max' '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' '00_saliency' '_' cat_nii_nme{iC} '_' nor_dta_nme{iNO} '_max.csv'],       [ 'ROI' stt_nme ; roi_max ])
        end
    end
end

%% Plots
mes_int = { 'median' };

xdt_val = [ 1 1 1];
xdt_adj = [-.3 0 .3];

ttl_hld = cat_nii_nme{cat_nii_nme_int};
% ttl_hld = strcat(cat_nii_nme,' VS ');
% ttl_hld = cat(2,ttl_hld{:});
% ttl_hld = strrep(ttl_hld,'VS','VS ');
% ttl_hld = ttl_hld(1:end-4);

for iM = 1:numel(mes_int)
    for iNO = 2%1:numel(nor_dta_nme)
        for iK = 2%1:numel(kep_dta_nme)
            plt_out_dir = [ sal_out_dir '/' 'ROIs' '/' 'median' '/' nor_dta_nme{iNO} '_' kep_dta_nme{iK} '/' ]; ejk_chk_dir(plt_out_dir);
            sprintf('Working on %s',plt_out_dir)
            
            for iR = 1:size(roi_med,1)
                
                pcfg = [];
                
                cnt = 1;
                for iC = cat_nii_nme_int%1:numel(cat_nii_nme)
                    for iS = 1:numel(sal_grp)
                        pcfg.ydt{cnt} = plt_hld.(nor_dta_nme{iNO}).(kep_dta_nme{iK}).(roi_med{iR,1}).(cat_nii_nme{iC}).(mes_int{iM})(:,iS);
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
                
                pcfg.out_dir = [ plt_out_dir '/' cat_nii_nme{cat_nii_nme_int} '/' ];
                pcfg.out_nme = [ 'saliency' '_' cat_nii_nme{cat_nii_nme_int} '_' roi_med{iR,1} '_' mes_int{iM} '_' nor_dta_nme{iNO} ];
                
                try ejk_scatter(pcfg); catch; end
                
            end
        end
    end
end

%% Correlations
% Load & Create Data
fcfg = [];
fcfg.dta_loc = [ sal_out_dir '/' 'ROIs' '/' 'median' '/' 'zsc_dta' '_' 'remove' '/' '00_saliency' '_' 'ep' '_' 'zsc_dta' '_median.csv'];
[ dta_csv, roi_csv, col_csv ] = ejk_dta_frm(fcfg);
    dta_csv = cell2mat(dta_csv);
    
col_csv(1,end+1) = {'3D_vs_2D'};
dta_csv(:,end+1) = dta_csv(:,strcmpi(col_csv,'3D')) - dta_csv(:,strcmpi(col_csv,'2D'));

% Stats
[rho_bet,pvl_bet] = corr(dta_csv(:,strcmpi(col_csv,'2D')),dta_csv(:,strcmpi(col_csv,'3D')));
[rho_dff,pvl_dff] = corr(dta_csv(:,strcmpi(col_csv,'3D_vs_2D')),dta_csv(:,strcmpi(col_csv,'3D')));

cell2csv([ sal_out_dir '/' 'ROIs' '/' 'median' '/' 'zsc_dta' '_' 'remove' '/' 'stats.csv'],[ {'2D BY 3D'} rho_bet pvl_bet ; {'3D-2D BY 3D'} rho_dff pvl_dff])

% Scatterplot
pcfg = [];

pcfg.ydt = {dta_csv(:,strcmpi(col_csv,'2D'))};
pcfg.xdt = {dta_csv(:,strcmpi(col_csv,'3D'))};

pcfg.fce_col = {rgb('dark grey')};

pcfg.edg_col = repmat({[0 0 0]},1,numel(pcfg.xdt));

pcfg.xlb = {'3D Saliency (z)'};
pcfg.ylb = {'2D Saliency (z)'};

pcfg.mkr_sze = repmat(40,1,numel(pcfg.xdt));
pcfg.aph_val = 0.80;
pcfg.trd_lne = 1;

pcfg.xlm = [-0.5 6.5];
pcfg.ylm = [-0.5 6.5];

pcfg.out_dir = [ sal_out_dir '/' 'ROIs' '/' 'median' '/' 'zsc_dta' '_' 'remove' '/' ];
pcfg.out_nme = [ 'scatterplot_importance_between.png' ];

ejk_scatter(pcfg);

% Scatterplot
pcfg = [];

pcfg.ydt = {dta_csv(:,strcmpi(col_csv,'3D_vs_2D'))};
pcfg.xdt = {dta_csv(:,strcmpi(col_csv,'3D'))};

pcfg.fce_col = {rgb('dark grey')};

pcfg.edg_col = repmat({[0 0 0]},1,numel(pcfg.xdt));

pcfg.xlb = {'3D Saliency (z)'};
pcfg.ylb = {'3D - 2D (z)'};

pcfg.xlm = [-0.5 6.5];

pcfg.mkr_sze = repmat(40,1,numel(pcfg.xdt));
pcfg.aph_val = 0.80;
pcfg.trd_lne = 1;

pcfg.out_dir = [ sal_out_dir '/' 'ROIs' '/' 'median' '/' 'zsc_dta' '_' 'remove' '/' ];
pcfg.out_nme = [ 'scatterplot_importance_subtraction.png' ];

ejk_scatter(pcfg);






