clear; clc;

%%
load( [ '/home/ekaestne/PROJECTS/OUTPUT' '/' 'PostOperative/Naming' '/' 'groups.mat' ] );

cog_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'Cognitive'                '_' 'QC' '.csv'];;
cog_dta = mmil_readtext(cog_dta_nme);
    cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
    cog_dta_sbj = cog_dta(2:end,1);
    cog_dta     = cell2mat(cog_dta(2:end,2:end));

wmp_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'];
wmp_dta = mmil_readtext(wmp_dta_nme);
    wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
    wmp_dta_sbj = wmp_dta(2:end,1);
    wmp_dta_rcn = wmp_dta(2:end,2);
    wmp_dta     = wmp_dta(2:end,5:end);

fsr_dir     = '/home/mmilmcdRSI/data/fsurf/';
fsr_dir_hld = dir(fsr_dir);
    fsr_dir_hld = {fsr_dir_hld(:).name};

avg_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/';
avg_sbj = 'fsaverage';
    
%% Update ROI to include less dead space
% Load
[ lhs_prc_loc, lhs_prc_lbl, lhs_prc_ctb ]=fs_read_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus' '.aparc.annot' ]);
[ lhs_old_loc, ~, lhs_old_ctb ]          =fs_read_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.aparc.annot' ]);

[ rhs_prc_loc, rhs_prc_lbl, rhs_prc_ctb ]=fs_read_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus' '.aparc.annot' ]);
[ rhs_old_loc, ~, rhs_old_ctb ]          =fs_read_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'rh' '.aparc.annot' ]);

% Label lhs
prc_ded_reg = find(lhs_prc_loc==0);
org_ded_reg = find(lhs_old_loc==0);
prc_ded_reg = setxor(prc_ded_reg, org_ded_reg);

lhs_prc_loc( prc_ded_reg ) = 7;
lhs_prc_lbl{7} = 'deadzone';
    
lhs_prc_ctb.numEntries = 7;
lhs_prc_ctb.struct_names{7} = 'deadzone';
lhs_prc_ctb.table = [ lhs_prc_ctb.table ; lhs_old_ctb.table(20,:) ];

% Label rhs
prc_ded_reg = find(rhs_prc_loc==0);
org_ded_reg = find(rhs_old_loc==0);
prc_ded_reg = setxor(prc_ded_reg, org_ded_reg);

rhs_prc_loc( prc_ded_reg ) = 7;
rhs_prc_lbl{7} = 'deadzone';
    
rhs_prc_ctb.numEntries = 7;
rhs_prc_ctb.struct_names{7} = 'deadzone';
rhs_prc_ctb.table = [ rhs_prc_ctb.table ; rhs_old_ctb.table(20,:) ];

% Ctab nonsense
lhs_sve_ctb = zeros(size(lhs_prc_loc));
rhs_sve_ctb = zeros(size(rhs_prc_loc));
for iV = 1:size(lhs_sve_ctb,1)

    if lhs_prc_loc(iV)==0
        lhs_sve_ctb(iV)=0;
    else
        lhs_sve_ctb(iV) = lhs_prc_ctb.table( lhs_prc_loc(iV),1) + lhs_prc_ctb.table( lhs_prc_loc(iV),2)*2^8 + lhs_prc_ctb.table( lhs_prc_loc(iV),3)*2^16;
    end
    
    if rhs_prc_loc(iV)==0
        rhs_sve_ctb(iV)=0;
    else
        rhs_sve_ctb(iV) = rhs_prc_ctb.table( rhs_prc_loc(iV),1) + rhs_prc_ctb.table( rhs_prc_loc(iV),2)*2^8 + rhs_prc_ctb.table( rhs_prc_loc(iV),3)*2^16;
    end
    
end

% Save
fs_write_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus' '.aparc.annot' ], ...
                     lhs_prc_loc, lhs_sve_ctb, lhs_prc_ctb)
fs_write_annotation([ avg_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus' '.aparc.annot' ], ...
                      rhs_prc_loc, rhs_sve_ctb, rhs_prc_ctb)

%%
for iS = 1:numel(grp.tle_controls_pre_3T_allSurg_all)
    
    idv_sbj = fsr_dir_hld{string_find(fsr_dir_hld,wmp_dta_rcn{grp.tle_controls_pre_3T_allSurg_all(iS)})};
    
    unix(['cp' ' ' '-r' ' ' fsr_dir '/' idv_sbj ' ' avg_dir '/'])
    
    % Resample (fsaverage to individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_parms = [];
    tmp_parms.outdir      = [ avg_dir '/' idv_sbj '/' 'label' '/' ];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj        = idv_sbj;
    tmp_parms.subjdir     = [ avg_dir '/' ];
    tmp_parms.verbose     = 1;
    tmp_parms.forceflag   = 0;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot( [ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus' '.aparc.annot' ], args{:});
    
    % RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_parms = [];
    tmp_parms.outdir      = [ avg_dir '/' idv_sbj '/' 'label' '/' ];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj        = idv_sbj;
    tmp_parms.subjdir     = [ avg_dir '/' ];
    tmp_parms.verbose     = 1;
    tmp_parms.forceflag   = 0;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot( [ avg_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus' '.aparc.annot' ], args{:});
    
    % Create Volume (individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname_tmp = sprintf('%s/mri/%s+aseg.mgz', [ avg_dir '/' idv_sbj '/' ], ['aparc.new_fus']);
    
    cmd = [];
    cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
    cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj);
    cmd = sprintf('%s   --annot %s', cmd, ['new_fus' '.aparc']);
    cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['new_fus' '.aparc'] );
    cmd = sprintf('%s   --volmask',cmd);
    cmd = sprintf('%s   --o %s\n',cmd,fname_tmp);
    [s,r] = unix(cmd);

    % Create wmparc (individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname_out = sprintf('%s/%s.mgz',[ avg_dir '/' idv_sbj '/' 'mri' '/' ],['wmparc' '_' 'new_fus']);
    
    cmd = [];
    cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
    cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj );
    cmd = sprintf('%s   --annot %s', cmd, ['new_fus' '.aparc'] );
    cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['new_fus' '.aparc']);
    cmd = sprintf('%s   --wmparc-dmax 5', cmd);
    cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown', cmd);
    cmd = sprintf('%s   --volmask', cmd);
    cmd = sprintf('%s   --ctxseg %s+aseg.mgz', cmd, ['aparc.new_fus']);
    cmd = sprintf('%s   --o %s', cmd, fname_out);
    [s,r] = unix(cmd);
    
    %
%     fname_out = sprintf('%s/%s.mgz',[ avg_dir '/' idv_sbj '/' 'mri' '/' ],['wmparc' '_' 'new_fus' '_' 'bare']);
%     
%     cmd = [];
%     cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
%     cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj );
%     cmd = sprintf('%s   --annot %s', cmd, ['new_fus' '.aparc'] );
%     cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['new_fus' '.aparc']);
%     cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown ', cmd); %--ribbon
%     cmd = sprintf('%s   --ctxseg %s+aseg.mgz', cmd, ['aparc.new_fus']);
%     cmd = sprintf('%s   --o %s', cmd, fname_out);
%     [s,r] = unix(cmd);
%     
%     %
%     fname_out = sprintf('%s/%s.mgz',[ avg_dir '/' idv_sbj '/' 'mri' '/' ],['wmparc' '_' 'old_fus' '_' 'bare']);
%         
%     cmd = [];
%     cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
%     cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj );
%     cmd = sprintf('%s   --annot %s', cmd, ['aparc'] );
%     cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['aparc']);
%     cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown ', cmd); %--ribbon
%     cmd = sprintf('%s   --ctxseg %s+aseg.mgz', cmd, ['aparc']);
%     cmd = sprintf('%s   --o %s', cmd, fname_out);
%     [s,r] = unix(cmd);
    
    % %%%%%%%%%%%%%%%%%%
%     [vol_wmparc,M_wmparc] = fs_load_mgh(fname_out);
%     [vol_wmparc_res,M_res] = mmil_resample_vol(vol_wmparc,M_wmparc,...
%         'nvox_ref',size(parms.vol),'M_ref',parms.M,...
%         'interpm',0,'M_reg',inv(parms.M_reg));
%     fs_save_mgh(vol_wmparc_res,fname_wmparc_res,M_res);
    
    % Extract new wmparc (individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unix(['rm' ' ' '-r' ' ' avg_dir '/' idv_sbj '/'])
    
end

%% %%%%%%%%%%%%%%%%%%
dti_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';

sbj_nme = 'DTIPROC_epd005_fmri_BLOCK_130507_20130507.153624_1';

% org_seg_dta = ejk_load_vol([ dti_dir '/' sbj_nme '/' 'DTanalysis' '/' 'wmparc_resDTI.mgz' ]);

% [RegInfo,fname_reg,errcode] = DTI_MMIL_Load_RegInfo([ dti_dir '/' sbj_nme], ...
%                               'infix','corr');
load([ dti_dir '/' sbj_nme '/' 'DTI1_f0_corr_regT1_regT1.mat'])

[vol_wmparc,M_wmparc] = fs_load_mgh(['/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri' '/' 'wmparc_new_fus.mgz']);
[vol_wmparc_res,M_res] = mmil_resample_vol(vol_wmparc,M_wmparc,...
    'nvox_ref',RegInfo.volsz_T2,'M_ref',RegInfo.M_T2,...
    'interpm',0,'M_reg',inv(RegInfo.M_T1_to_T2));
% fs_save_mgh(vol_wmparc_res,fname_wmparc_res,M_res);

% MPR_res_info.mat
% M: 4x4 vox2ras matrix for vol
% [vol_wmparc_res,M_res] = mmil_resample_vol(vol_wmparc,M_wmparc,...
%     'nvox_ref',size(parms.vol),'M_ref',parms.M,...
%     'interpm',0,'M_reg',inv(parms.M_reg));

% Test
col_loc = mmil_readtext(['/home/ekaestne/PROJECTS/SCRIPTS/streamint' '/' 'wmparc' '.aparc' '' '.' 'annot']);

%% Check
idv_sbj = 'DTIPROC_epd005_fmri_BLOCK_130507_20130507.153624_1';

dti_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';

vol_dta_nme = [dti_dir '/' 'DTIPROC' '_' idv_sbj(9:end) '/' 'DTcalc' '/' 'DTI1_crev_corr_regT1_DT_' 'FA' '.mgz'];
vol_dta_sze = ejk_read_vol_sze(vol_dta_nme,1);
vol_dta = ejk_load_vol(vol_dta_nme);

seg_dta_nme = [dti_dir '/' 'DTIPROC' '_' idv_sbj(9:end) '/' 'DTanalysis' '/' 'wmparc_resDTI.mgz']; % - EJK - Fix to allow multiple label types
seg_dta_sze = ejk_read_vol_sze(seg_dta_nme,1);
seg_dta = ejk_load_vol(seg_dta_nme);

seg_dta_new_nme = ['/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri' '/' 'wmparc_new_fus.mgz']; % - EJK - Fix to allow multiple label types
seg_dta_new_sze = ejk_read_vol_sze(seg_dta_new_nme,1);
seg_dta_new = ejk_load_vol(seg_dta_new_nme);

min_val = 1e-06;
scl_fct = 1;

[ lhs_prc_loc, lhs_prc_lbl, lhs_prc_ctb ]=fs_read_annotation( [ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.aparc.annot' ] );
[ lhs_prc_lbl([ 6 13 16 25 32 ]) ...
col_loc([ 6 13 16 25 32 ],2) ]

% Check fsaverage annot %%%%%%%%%%%%
% lhs
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/fsaverage/label/';

[ old_prc_loc, old_prc_lbl, old_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'lh' '.aparc.annot' ] );

[ new_prc_loc, new_prc_lbl, new_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'lh' '.new_fus.aparc.annot' ] );

new_ind_use = [ 1 2  3  4  5 ]; % 3006 ];
org_ind_use = [ 6 13 16 25 32 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end

% rhs
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/fsaverage/label/';

[ old_prc_loc, old_prc_lbl, old_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'rh' '.aparc.annot' ] );

[ new_prc_loc, new_prc_lbl, new_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'rh' '.new_fus.aparc.annot' ] );

new_ind_use = [ 1 2  3  4  5 ]; % 3006 ];
org_ind_use = [ 6 13 16 25 32 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end

% Check individual annot %%%%%%%%%%%%
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/label/';

[ old_prc_loc, old_prc_lbl, old_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'lh' '.aparc.annot' ] );

[ new_prc_loc, new_prc_lbl, new_prc_ctb ] = ...
    fs_read_annotation( [ dir_hld '/' 'lh' '.new_fus.aparc.annot' ] );

new_ind_use = [ 2 3  4  5  6 ]; 
org_ind_use = [ 6 13 16 25 32 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end

% Check individual aparc+aseg %%%%%%%%%%%%
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri/';

old_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'aparc+aseg.mgz' ] );

new_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'aparc.new_fus+aseg.mgz' ] );

new_ind_use = [ 1001 1002 1003 1004 1005 ]; 
org_ind_use = [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end

% Check individual aparc+aseg - MRIres %%%%%%%%%%%%
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri/';

old_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc.mgz' ] );

new_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc_new_fus.mgz' ] );

new_ind_use = [ 3001 3002 3003 3004 3005 ]; % [ 1001 1002 1003 1004 1005 ]; 
org_ind_use = [ 3005 3012 3015 3024 3031 ]; % [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end
out_hld(:,1) ./ out_hld(:,2)

% Check individual aparc+aseg - MRIres bare call %%%%%%%%%%%%
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri/';

old_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc.mgz' ] );

new_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc_new_fus_bare.mgz' ] );

new_ind_use = [ 3001 3002 3003 3004 3005 ]; % [ 1001 1002 1003 1004 1005 ]; 
org_ind_use = [ 3005 3012 3015 3024 3031 ]; % [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end
out_hld(:,1) ./ out_hld(:,2)

new_ind_use = [ 2 11 63 251 255 ]; % [ 1001 1002 1003 1004 1005 ]; 
org_ind_use = [ 2 11 63 251 255 ]; % [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end
out_hld(:,1) ./ out_hld(:,2)

tt1 = tabulate(new_prc_loc(:));
tt2 = tabulate(old_prc_loc(:));

[tt1(1:15,1:2) tt2(1:15,2)]

[tt1(31:45,1:2) tt2(31:45,2)]

[tt1(ismember(tt1(:,1),[1001:1005 3001:3005]),1:2) tt2(ismember(tt2(:,1),[ 1005 1012 1015 1024 1031 3005 3012 3015 3024 3031 ]),[2 1])]

hld_new_ind = find(new_prc_loc(:)==1001);
hld_old_ind = find(old_prc_loc(:)==1005);
numel( intersect(hld_new_ind,hld_old_ind)) / numel(hld_old_ind)

hld_new_ind = find(new_prc_loc(:)==3001);
hld_old_ind = find(old_prc_loc(:)==3005);
numel( intersect(hld_new_ind,hld_old_ind)) / numel(hld_old_ind)

hld_new_ind = find(new_prc_loc(:)==3002);
hld_old_ind = find(old_prc_loc(:)==3012);
numel( intersect(hld_new_ind,hld_old_ind)) / numel(hld_old_ind)

hld_new_ind = find(new_prc_loc(:)==3005);
hld_old_ind = find(old_prc_loc(:)==3031);
numel( intersect(hld_new_ind,hld_old_ind)) / numel(hld_old_ind)

% Check individual aparc+aseg - MRIres orig bare call %%%%%%%%%%%%
dir_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/FSURF_epd005_fmri_BLOCK_130507_20130507.153624_1/mri/';

old_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc.mgz' ] );

new_prc_loc = ...
    ejk_load_vol( [ dir_hld '/' 'wmparc_old_fus_bare.mgz' ] );

new_ind_use = [ 3005 3012 3015 3024 3031 ]; % [ 1001 1002 1003 1004 1005 ]; 
org_ind_use = [ 3005 3012 3015 3024 3031 ]; % [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( new_prc_loc(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( old_prc_loc(:) == org_ind_use(iI) );
    
end
out_hld(:,1) ./ out_hld(:,2)

% Check Individual WMPARC - DTIres %%%%%%%%%%%%
new_ind_use = [ 3001 3002 3003 3004 3005 ]; % 3006 ];
org_ind_use = [ 3005 3012 3015 3024 3031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( vol_wmparc_res(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( seg_dta(:)        == org_ind_use(iI) );
    
    ovr_lap(iI) = numel(intersect( find(vol_wmparc_res(:) == new_ind_use(iI)), find(seg_dta(:)        == org_ind_use(iI)) )) / ...
        sum(seg_dta(:)        == org_ind_use(iI));
    
end
out_hld(:,1) ./ out_hld(:,2)



new_ind_use = [ 4001 4002 4003 4004 4005 ]; % 3006 ];
org_ind_use = [ 4005 4012 4015 4024 4031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( vol_wmparc_res(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( seg_dta(:)        == org_ind_use(iI) );
    
    ovr_lap(iI) = numel(intersect( find(vol_wmparc_res(:) == new_ind_use(iI)), find(seg_dta(:)        == org_ind_use(iI)) )) / ...
        sum(seg_dta(:)        == org_ind_use(iI));
    
end
out_hld(:,1) ./ out_hld(:,2)


new_ind_use = [ 1001 1002 1003 1004 1005 ]; % 3006 ];
org_ind_use = [ 1005 1012 1015 1024 1031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( vol_wmparc_res(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( seg_dta(:)        == org_ind_use(iI) );
   
    ovr_lap(iI) = numel(intersect( find(vol_wmparc_res(:) == new_ind_use(iI)), find(seg_dta(:)        == org_ind_use(iI)) )) / ...
        sum(seg_dta(:)        == org_ind_use(iI));
    
end
out_hld(:,1) ./ out_hld(:,2)


new_ind_use = [ 2001 2002 2003 2004 2005 ]; % 3006 ];
org_ind_use = [ 2005 2012 2015 2024 2031 ];
for iI = 1:numel(new_ind_use)
    
    out_hld(iI,1) = sum( vol_wmparc_res(:) == new_ind_use(iI) );
    out_hld(iI,2) = sum( seg_dta(:)        == org_ind_use(iI) );
    
    ovr_lap(iI) = numel(intersect( find(vol_wmparc_res(:) == new_ind_use(iI)), find(seg_dta(:)        == org_ind_use(iI)) )) / ...
        sum(seg_dta(:)        == org_ind_use(iI));
    
end
out_hld(:,1) ./ out_hld(:,2)

% %%%%%%%%%%%%
ind_use = [ 3001 3002 3003 3004 3005 3006 ];
for iI = 1:numel(ind_use)
    
    val_hld = vol_dta(find(seg_dta_new==ind_use(iI))) * scl_fct;
    num_val = size(val_hld,1);
    
    ind_vld = find(min(abs(val_hld),[],2)>=cfg.min_val & ~isnan(sum(val_hld,2)));
    ind_invalid = setdiff([1:num_val],ind_vld);
    ind_nan = isnan(sum(val_hld,2));
    
    val_hld(ind_nan) = 0;
    
    wmp_avg(1,iFC) = nanmean(val_hld(ind_vld,:),1);
    wmp_med(1,iFC) = nanmedian(val_hld(ind_vld,:),1);
    if numel(ind_vld)>1
        wmp_std(1,iFC) = nanstd(val_hld(ind_vld,:),1);
    end
    
end























