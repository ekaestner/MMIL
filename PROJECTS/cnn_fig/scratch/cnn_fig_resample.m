org_atl_dir = '/home/darienzo/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym';
new_atl_dir = '/space/mcdonald-syn01/1/projects/ekaestner/resample_test';
exp_sbj = 'enigma_conglom_sub-UCSDepd046_ses-1_cat12_mwp1.nii';
slc_num = 71;

nii_exp_sbj = niftiread([new_atl_dir '/' exp_sbj]);
nii_exp_sbj_slc_prp = 71 / size(nii_exp_sbj,2);

thr_dim_rsm_cmd{1} = '3dresample -rmode NN -dxyz ';
thr_dim_rsm_cmd{2} = '-input ';
thr_dim_rsm_cmd{3} = '-prefix ';

%%
lbl_nme = dir(org_atl_dir); lbl_nme = {lbl_nme(:).name}; lbl_nme = lbl_nme(string_find(lbl_nme,'.txt'));
for iL = 1:numel(lbl_nme)
    copyfile( [org_atl_dir '/' lbl_nme{iL}], [new_atl_dir '/' 'citations' '/' lbl_nme{iL}] );
end

%%
lbl_nme = dir(org_atl_dir); lbl_nme = {lbl_nme(:).name}; lbl_nme = lbl_nme(string_find(lbl_nme,'.csv'));
for iL = 1:numel(lbl_nme)
    copyfile( [org_atl_dir '/' lbl_nme{iL}], [new_atl_dir '/' 'labels' '/' lbl_nme{iL}] );
end

%%
atl_nme = dir(org_atl_dir); atl_nme = {atl_nme(:).name}; atl_nme = atl_nme(string_find(atl_nme,'.nii'));
atl_rpt = cell(numel(atl_nme),7);
for iA = 1:numel(atl_nme)
    %
    copyfile( [org_atl_dir '/' atl_nme{iA}], [new_atl_dir '/' 'orig' '/' atl_nme{iA}] );
    
    %
    nii_org = niftiread([new_atl_dir '/' 'orig' '/' atl_nme{iA}]);
    scl_256 = size(nii_org,1) / 256;
        scl_256_str = num2str(scl_256);
    scl_113 = size(nii_org,1) / 113;
        scl_113_str = num2str(scl_113);

    %
    cmd_256 = [ thr_dim_rsm_cmd{1} scl_256_str ' ' scl_256_str ' ' scl_256_str ' ' thr_dim_rsm_cmd{2} [new_atl_dir '/' 'orig' '/' atl_nme{iA}] ' ' thr_dim_rsm_cmd{3} [new_atl_dir '/' '256' '/' atl_nme{iA}] ];
    cmd_113 = [ thr_dim_rsm_cmd{1} scl_113_str ' ' scl_113_str ' ' scl_113_str ' ' thr_dim_rsm_cmd{2} [new_atl_dir '/' 'orig' '/' atl_nme{iA}] ' ' thr_dim_rsm_cmd{3} [new_atl_dir '/' '113' '/' atl_nme{iA}] ];

    %
    system(cmd_256)
    system(cmd_113)

    %
    nii_256 = niftiread([new_atl_dir '/' '256' '/' atl_nme{iA}]);
    nii_113 = niftiread([new_atl_dir '/' '113' '/' atl_nme{iA}]);

    %
    slc_num_org = round(size(nii_org,2) * nii_exp_sbj_slc_prp);
    slc_num_256 = round(size(nii_256,2) * nii_exp_sbj_slc_prp);
    slc_num_113 = round(size(nii_113,2) * nii_exp_sbj_slc_prp);

    %
    figure('Visible','off');
    subplot(2,2,1);
    imagesc(rot90(squeeze(nii_exp_sbj(:,slc_num,:))));
    subtitle('example subject');
    subplot(2,2,2);
    imagesc(rot90(squeeze(nii_org(:,slc_num_org,:))));
    subtitle([ mmil_spec_char(atl_nme{iA},{'_'},{' '}) ' ' 'original space']);
    subplot(2,2,3);
    imagesc(rot90(squeeze(nii_256(:,slc_num_256,:))));
    subtitle([ mmil_spec_char(atl_nme{iA},{'_'},{' '}) ' ' '256x256 space']);
    subplot(2,2,4);
    imagesc(rot90(squeeze(nii_113(:,slc_num_113,:))));
    subtitle([ mmil_spec_char(atl_nme{iA},{'_'},{' '}) ' ' '113/113 space']);
    print([new_atl_dir '/' 'images' '/' atl_nme{iA}],'-dpng');
    close all

    %
    atl_rpt{iA,1} = atl_nme{iA};
    atl_rpt{iA,2} = size(nii_org);
    atl_rpt{iA,3} = slc_num_org;
    atl_rpt{iA,4} = size(nii_256);
    atl_rpt{iA,5} = slc_num_256;
    atl_rpt{iA,6} = size(nii_113);
    atl_rpt{iA,7} = slc_num_113;

end
cell2csv([ new_atl_dir '/' 'report' '_' dte_str '.csv'], [ {'Atlas'} {'original size'} {['original slice ' num2str(slc_num)]} {'256x256 size'} {['256x256 slice ' num2str(slc_num)]} {'113x113 size'} {['113x113 slice ' num2str(slc_num)]}  ; atl_rpt])

%%
% prj_dir = '/home/darienzo/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/';
% cat_dir = '/space/mcdonald-syn01/1/projects/ekaestner/BIDS_dev/CAT12_investigation/BIDS_dev_cat12_v2/derivatives/CAT12.9';
% 
% nii_sch_100 = 'Schaefer2018_100Parcels_17Networks_order.nii';
% brn_msk = 'brainmask.nii';
% jul_brn = 'julichbrain.nii';
% aal_brn = 'aal3.nii';
% 
% nii_sch_100_dta = niftiread([ prj_dir '/' nii_sch_100 ]);
%     size(nii_sch_100_dta)
% brn_msk_dta = niftiread([ prj_dir '/' brn_msk ]);
%     size(brn_msk_dta)
% jul_brn_dta = niftiread([ prj_dir '/' jul_brn ]);
%     size(jul_brn_dta)
% aal_brn_dta = niftiread([ prj_dir '/' aal_brn ]);
%     size(jul_brn_dta)
% 
% subplot(2,2,1)
% imagesc(rot90(squeeze(nii_sch_100_dta(:,102,:)))) 
% subplot(2,2,2)
% imagesc(rot90(squeeze(brn_msk_dta(:,56,:))))
% 
% sze_chk_256 = niftiread([new_atl_dir '/' '256' '/' atl_nme{iA}]);
% size(sze_chk_256)
% 
% 3dresample -dxyz 1.0 1.0 0.9 -prefix 119.dset -input in+tlrc