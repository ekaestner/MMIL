% Get dummy size
nft_fle = dir([ nii_dir '/' sal_grp{1} '/' '*.nii']); nft_fle = {nft_fle(1).name};
dum_dta = niftiread([ nii_dir '/' sal_grp{1} '/' nft_fle{1}]);

%% Load Atlas
% Reslice Atlas
flg = [];
flg.prefix = 'atlas_to_map_reslice_';
flg.interp = 0;
flg.which  = 1;
flg.mean   = 0;
spm_reslice( [ {[ nii_dir '/' sal_grp{1} '/' nft_fle{1}]}; {[ atl_dir '/' atl_nme '.nii' ]} ], flg);

% Load Atlas
nft_atl = niftiread([ atl_dir '/' flg.prefix atl_nme '.nii' ]);

nft_atl_lbl = mmil_readtext([ atl_dir '/' atl_nme '.txt' ]);
    nft_atl_lbl = regexpi(nft_atl_lbl,'\|','split');
    nft_atl_lbl = cat(1,nft_atl_lbl{:});
    
%% Collate ROI data
nor_nme = { 'zsc' } ; % 'all'
kep_nme = { 'remove' }; % 'keep' 

cmp_nme = { sal_grp };
cmp_num = numel(cmp_nme);
stt_nme = [ sal_grp cell(1,numel(cmp_nme)*2) ];

one_num = size(sal_dta_nor.(nor_nme{1}).(kep_nme{1}).(['x' sal_grp{1}]),4);
two_num = size(sal_dta_nor.(nor_nme{1}).(kep_nme{1}).(['x' sal_grp{2}]),4);
sbj_num = one_num + two_num;

for iNO = 1:numel(nor_nme)
    for iK = 1:numel(kep_nme)
        
        ejk_chk_dir([ out_fld '/' 'ROIs' '/' 'median' '/' nor_nme{iNO} '_' kep_nme{iK} '/' ]);
        
        roi_med = cell(size(nft_atl_lbl,1),1+numel(sal_grp)+2*cmp_num);
        roi_sbj = cell( sbj_num+1, size(nft_atl_lbl,1)+2);
                
        for iR = 1:size(roi_med,1)
        
            roi_sbj_int = 2;
            
            roi_med{iR,1}   = nft_atl_lbl{iR,3};
            roi_sbj{1,iR+2} = nft_atl_lbl{iR,3};
            
            sal_dta_tmp = sal_dta_nor.(nor_nme{iNO}).(kep_nme{iK}).(['x' sal_grp{1}]);
            
            roi_med_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
            roi_iqr_hld = nan(size(sal_dta_tmp,4),numel(sal_grp));
            
            for iSA = 1:numel(sal_grp)
                for iS = 1:size(sal_dta_nor.(nor_nme{iNO}).(kep_nme{iK}).(['x' sal_grp{iSA}]),4)
                    
                    sal_hld = squeeze(sal_dta_nor.(nor_nme{iNO}).(kep_nme{iK}).(['x' sal_grp{iSA}])(:,:,:,iS));
                    sal_hld = sal_hld(nft_atl==str2num(nft_atl_lbl{iR,1}));
                    
                    roi_med_hld(iS,iSA) = nanmedian(sal_hld);
                    roi_iqr_hld(iS,iSA) = iqr(sal_hld);

                    roi_sbj{roi_sbj_int,2}    = sal_grp{iSA};
                    roi_sbj{roi_sbj_int,iR+2} = roi_med_hld(iS,iSA);
                    roi_sbj_int = roi_sbj_int + 1;
                    
                end

                roi_med{iR,iSA+1} = [ num2str(roundsd(nanmedian(roi_med_hld(iS,iSA)),2))];% ' (' num2str(roundsd(nanmean(roi_iqr_hld(iS,iSA)),2)) ')' ];

            end
            
            % Tests
            for iT = 1:numel(cmp_nme)
                
                ind = 0+((iT-1)*2):1+((iT-1)*2);
                
                col_one = strcmpi(sal_grp,cmp_nme{iT}{1});
                col_two = strcmpi(sal_grp,cmp_nme{iT}{2});
                
                stt_nme(1,end-ind(1)) = {[ cmp_nme{iT}{1} '_VS_' cmp_nme{iT}{2} '_ttest' ]};
                stt_nme(1,end-ind(2)) = {[ cmp_nme{iT}{1} '_VS_' cmp_nme{iT}{2} '_signrank']};
                
                [ ~, roi_med{iR,end-ind(1)} ] = ttest2(roi_med_hld(1:one_num,col_one),roi_med_hld(1:two_num,col_two));
                roi_med{iR,end-ind(2)}        = ranksum(roi_med_hld(1:one_num,col_one),roi_med_hld(1:two_num,col_two));
                
            end
            
            % Hold
            plt_hld.(nor_nme{iNO}).(kep_nme{iK}).(roi_med{iR}).median = roi_med_hld;
            
        end
        
        cell2csv( [ out_fld '/' 'ROIs' '/' 'median' '/' nor_nme{iNO} '_' kep_nme{iK}  '/' '00_median_saliency' '_' nor_nme{iNO} '_' kep_nme{iK} '_median.csv'],            [ 'ROI' stt_nme ; roi_med ]);
        cell2csv( [ out_fld '/' 'ROIs' '/' 'median' '/' nor_nme{iNO} '_' kep_nme{iK}  '/' '00_median_saliency' '_' nor_nme{iNO} '_' kep_nme{iK} '_median_individual.csv'], roi_sbj );
    
    end
end

%% Preliminary plots
mes_int = { 'median' };
mdl_int_col = { rgb('dark red') rgb('green') };

for iM = 1:numel(mes_int)
    for iNO = 1%1:numel(nor_nme)
        for iK = 1%1:numel(kep_nme)
            
            plt_out_dir = [ out_fld '/' 'plots' '/' nor_nme{iNO} '_' kep_nme{iK} '/' ]; ejk_chk_dir(plt_out_dir);
            sprintf('Working on %s',plt_out_dir)
            
            for iR = 1:size(roi_med,1)
                
                pcfg = [];
                
                cnt = 1;
                for iS = 1:numel(sal_grp)
                    pcfg.ydt{cnt} = plt_hld.(nor_nme{iNO}).(kep_nme{iK}).(roi_med{iR,1}).(mes_int{iM})(:,iS);
                        pcfg.ydt{cnt}(pcfg.ydt{cnt}==0) = [];
                    pcfg.xdt{cnt} = cnt;
                    
                    pcfg.fce_col{cnt}     = mdl_int_col{iS};
                    pcfg.box_plt_col{cnt} = mdl_int_col{iS};
                    
                    pcfg.xlb{cnt} = sal_grp{iS};
                    
                    cnt = cnt + 1;
                end
                
                pcfg.edg_col = repmat({[0 0 0]},1,numel(pcfg.xdt));
                pcfg.box_plt = ones(1,numel(pcfg.xdt));
                
                pcfg.xlm = [ 0.5 max(cat(1,pcfg.xdt{:}))+0.5 ];
                pcfg.ylb = {[ roi_med{iR} ' : ' upper(mes_int{iM}) ]};
                                
                pcfg.mkr_sze = repmat(40,1,numel(pcfg.xdt));
                pcfg.aph_val = 0.80;
                pcfg.jtr_wdt = 0.15;
                pcfg.box_wdt = 0.20;
                
                pcfg.out_dir = plt_out_dir ;
                pcfg.out_nme = [ 'saliency' '_' roi_med{iR,1} '_' mes_int{iM} '_' nor_nme{iNO} ];
                
                try ejk_scatter(pcfg); catch; end
                
            end
        end
    end
end

