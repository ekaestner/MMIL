clear; clc; 

dta_dir = '/space/mcdonald-syn01/1/projects/ekaestner/cnn_fig/new_data/_cat12_t1w_data_repo/';

sbj_nme = { 'enigma_conglom_sub-USZENIGMAZURPAT24_ses-1_cat12_mwp1.nii' ; ...
            'enigma_conglom_sub-LangoneNY0255_ses-01_cat12_mwp1.nii' ; ...
            'ECP_dataset_bids_sub-ECPEC2100_ses-1A_cat12_mwp1.nii' ; ... ...
            'CAPES_sub-NYUC0007_ses-01_cat12_mwp1.nii' };

for iS = 1:numel(sbj_nme)
    subplot(4,4,iS)
    sbj_one = niftiread([dta_dir '/' sbj_nme{iS}]);
    imagesc(rot90(squeeze(sbj_one(:,71,:))))
end