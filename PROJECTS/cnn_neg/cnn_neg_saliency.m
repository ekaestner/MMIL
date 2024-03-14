
%% Load
% Get dummy size
nft_fle = dir([ nii_dir '/' sal_grp{1} '/' '*.nii']); nft_fle = {nft_fle(1).name};
dum_dta = niftiread([ nii_dir '/' sal_grp{1} '/' nft_fle{1}]);
clear nft_fle

% Load
for iSA = 1:numel(sal_grp)
    
    nft_fle.(sal_grp{iSA}) = dir([ nii_dir '/' sal_grp{iSA} '/' '*.nii']); nft_fle.(sal_grp{iSA}) = {nft_fle.(sal_grp{iSA})(:).name};    
    fle_nme.(['x' sal_grp{iSA}]) = nft_fle.(sal_grp{iSA});
    
    sal_dta.(['x' sal_grp{iSA}]) = nan( [ size(dum_dta) numel( fle_nme.(['x' sal_grp{iSA}])) ] );
    
    for iF = 1:numel(fle_nme.(['x' sal_grp{iSA}]))
        sal_dta.(['x' sal_grp{iSA}])(:,:,:,iF) = niftiread([ nii_dir '/' sal_grp{iSA} '/' fle_nme.(['x' sal_grp{iSA}]){iF} ]);
        
        % Remove 0's
        sal_dta.(['x' sal_grp{iSA}])(sal_dta.(['x' sal_grp{iSA}])==0) = NaN;
        
    end
end

%% Normalize
low_pct = 75;

% Put together data structure
for iZ = 1:numel(sal_grp)
    
    sal_dta_nor.all.keep.(['x' sal_grp{iZ}]) = sal_dta.(['x' sal_grp{iZ}]);
    
    for iF = 1:size(sal_dta_nor.all.keep.(['x' sal_grp{iZ}]),4)
        
        % Values to play with
        sal_val_hld = squeeze(sal_dta_nor.all.keep.(['x' sal_grp{iZ}])(:,:,:,iF));
        
        % Normalize
        sal_val_men = nanmean(sal_val_hld(:));
        sal_val_std = nanstd(sal_val_hld(:));
        sal_val_hld = (sal_val_hld-sal_val_men) ./ sal_val_std;
        
        % Remove values
        sal_val_hld_rmv     = sal_val_hld;
        sal_val_hld_rmv_mid = sal_val_hld;
        
        % Remove bottom
        pct_hld = prctile(sal_val_hld(:),low_pct);
        sal_val_hld_rmv(sal_val_hld<pct_hld) = NaN;
        
        % Re-assign
        sal_dta_nor.zsc.keep.(['x' sal_grp{iZ}])(:,:,:,iF)       = sal_val_hld;
        sal_dta_nor.zsc.remove.(['x' sal_grp{iZ}])(:,:,:,iF)     = sal_val_hld_rmv;
        
        clear sal_val_hld sal_val_hld_rmv
        
    end
end
