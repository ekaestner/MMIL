
%% Get dummy size
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
    
%% Load Data
for iSA = 1:numel(sal_grp)
    nft_fle = dir([ nii_dir '/' sal_grp{iSA} '/' '*.nii']); nft_fle = {nft_fle(:).name};
    
    for iNI = 1:numel(cat_nii_nme)   
        fle_nme.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI}) = nft_fle(string_find(nft_fle,cat_nii{iNI}));
        
        sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI}) = nan( [ size(dum_dta) numel( fle_nme.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})) ] );
        
        for iF = 1:numel(fle_nme.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI}))           
            sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})(:,:,:,iF) = niftiread([ nii_dir '/' sal_grp{iSA} '/' fle_nme.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI}){iF} ]);
                        
            % Remove 0's
            sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})==0) = NaN;
            
            % Values to play with
            sal_val_hld = squeeze(sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})(:,:,:,iF));
            
            % Remove lower values
            low_pct = 75;
            pct_hld = prctile(sal_val_hld(:),low_pct);
            sal_val_hld(sal_val_hld<pct_hld) = NaN;
            
            % Minmax            
            sal_val_max = max(sal_val_hld(:));
            sal_val_min = min(sal_val_hld(:));
            sal_val_hld = (sal_val_hld-sal_val_min) ./ (sal_val_max - sal_val_min);
           
            % Z-score
%             sal_val_men = nanmean(sal_val_hld(:));
%             sal_val_std = nanstd(sal_val_hld(:));
%             sal_val_hld = (sal_val_hld-sal_val_men) ./ sal_val_std;
            
            % Re-assign
            sal_dta.(['x' sal_grp{iSA}]).(cat_nii_nme{iNI})(:,:,:,iF) = sal_val_hld;
            
        end
    end
end

%% Save
save( [ dta_dir '/' 'saliency.mat' ], 'sal_dta', 'sal_grp', 'sal_nme', 'cat_nii_nme', 'cat_col', 'nft_atl', 'nft_atl_lbl')

