clear; clc;

%% Setup subjects
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

dti_dir     = '/home/mmilmcdRSI/data/proc_dti/';
    dti_dir_hld = dir(dti_dir);
    dti_dir_hld = {dti_dir_hld(:).name};
    
avg_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/';
avg_sbj = 'fsaverage';

%%
ejk_chk_dir( [ avg_dir '/' 'output' '/' ] );

[ ~, prc_lbl, ~ ]=fs_read_annotation( [ avg_dir '/' avg_sbj '/' 'label' '/' 'lh.new_fus_cor.aparc.annot' ] );
prc_lbl = [ strcat( 'lh_', prc_lbl )' strcat(  'rh_', prc_lbl )'  ];

for iS = 19:numel(grp.tle_post_3T_ATLonly_left)
    
    sbj_hld = string_find(fsr_dir_hld,wmp_dta_rcn{grp.tle_post_3T_ATLonly_left(iS)});
    if numel(sbj_hld)==2; sbj_hld = sbj_hld(end); end
    
    idv_sbj = fsr_dir_hld{sbj_hld};
    
    dti_sbj = dti_dir_hld{string_find(dti_dir_hld,wmp_dta_rcn{grp.tle_post_3T_ATLonly_left(iS)})};
    
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
    fname_out = fs_annot2annot( [ avg_dir '/' avg_sbj '/' 'label' '/' 'lh' '.' 'new_fus_cor' '.aparc.annot' ], args{:});
    
    % RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_parms = [];
    tmp_parms.outdir      = [ avg_dir '/' idv_sbj '/' 'label' '/' ];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj        = idv_sbj;
    tmp_parms.subjdir     = [ avg_dir '/' ];
    tmp_parms.verbose     = 1;
    tmp_parms.forceflag   = 0;
    args = mmil_parms2args(tmp_parms);
    fname_out = fs_annot2annot( [ avg_dir '/' avg_sbj '/' 'label' '/' 'rh' '.' 'new_fus_cor' '.aparc.annot' ], args{:});
    
    % Create Volume (individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname_tmp = sprintf('%s/mri/%s+aseg.mgz', [ avg_dir '/' idv_sbj '/' ], ['aparc.new_fus_cor']);
    
    cmd = [];
    cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
    cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj);
    cmd = sprintf('%s   --annot %s', cmd, ['new_fus_cor' '.aparc']);
    cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['new_fus_cor' '.aparc'] );
    cmd = sprintf('%s   --volmask',cmd);
    cmd = sprintf('%s   --o %s\n',cmd,fname_tmp);
    [s,r] = unix(cmd);

    % Create wmparc (individual subject) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname_out = sprintf('%s/%s.mgz',[ avg_dir '/' idv_sbj '/' 'mri' '/' ],['wmparc' '_' 'new_fus_cor']);
    
    cmd = [];
    cmd = sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir); % sprintf('%s setenv SUBJECTS_DIR %s\n',cmd , avg_dir) % sprintf('%s export SUBJECTS_DIR=''%s''\n',cmd , avg_dir)
    cmd = sprintf('%s mri_aparc2aseg --s %s', cmd, idv_sbj );
    cmd = sprintf('%s   --annot %s', cmd, ['new_fus_cor' '.aparc'] );
    cmd = sprintf('%s   --annot-table %s/label/%s.ctab', cmd, [ avg_dir '/' idv_sbj], ['new_fus_cor' '.aparc']);
    cmd = sprintf('%s   --wmparc-dmax 5', cmd);
    cmd = sprintf('%s   --labelwm --hypo-as-wm --rip-unknown', cmd);
    cmd = sprintf('%s   --volmask', cmd);
    cmd = sprintf('%s   --ctxseg %s+aseg.mgz', cmd, ['aparc.new_fus_cor']);
    cmd = sprintf('%s   --o %s', cmd, fname_out);
    [s,r] = unix(cmd);
    
    % Resample to DTI space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([ dti_dir '/' dti_sbj '/' 'DTI1_f0_corr_regT1_regT1.mat'])
    
    [vol_wmparc,M_wmparc] = fs_load_mgh([ avg_dir '/' idv_sbj '/' 'mri' '/' 'wmparc_new_fus_cor.mgz']);
    [vol_wmparc_res,M_res] = mmil_resample_vol(vol_wmparc,M_wmparc,...
        'nvox_ref',RegInfo.volsz_T2,'M_ref',RegInfo.M_T2,...
        'interpm',0,'M_reg',inv(RegInfo.M_T1_to_T2));
    fs_save_mgh(vol_wmparc_res,[ avg_dir '/' idv_sbj '/' 'mri' '/' 'wmparc_new_fus_cor_resDTI.mgz'],M_res);
    
    % Save relevant pieces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    lbl_dir = [ avg_dir '/' idv_sbj '/' 'label' '/' ];
    org_dir = [ avg_dir '/' idv_sbj '/' 'mri' '/' ];
    new_dir = [ avg_dir '/' 'output' '/' idv_sbj '/' ];
        ejk_chk_dir(new_dir);
    
    unix(['cp' ' ' '-r' ' ' lbl_dir '/' 'lh.new_fus_cor.aparc.annot'    ' ' new_dir '/' 'lh.new_fus_cor.aparc.annot' ]);
         unix(['cp' ' ' '-r' ' ' lbl_dir '/' 'lh.aparc.annot'    ' ' new_dir '/' 'lh.aparc.annot' ]);
    unix(['cp' ' ' '-r' ' ' lbl_dir '/' 'rh.new_fus_cor.aparc.annot'    ' ' new_dir '/' 'rh.new_fus_cor.aparc.annot' ]);
        unix(['cp' ' ' '-r' ' ' lbl_dir '/' 'rh.aparc.annot'    ' ' new_dir '/' 'rh.aparc.annot' ]);
    unix(['cp' ' ' '-r' ' ' org_dir '/' 'aparc.new_fus_cor+aseg.mgz'    ' ' new_dir '/' 'aparc.new_fus_cor+aseg.mgz']);
        unix(['cp' ' ' '-r' ' ' org_dir '/' 'aparc+aseg.mgz'    ' ' new_dir '/' 'aparc+aseg.mgz']);
    unix(['cp' ' ' '-r' ' ' org_dir '/' 'wmparc_new_fus_cor.mgz'        ' ' new_dir '/' 'wmparc_new_fus_cor.mgz' ]);
        unix(['cp' ' ' '-r' ' ' org_dir '/' 'wmparc.mgz'        ' ' new_dir '/' 'wmparc.mgz' ]);
    unix(['cp' ' ' '-r' ' ' org_dir '/' 'wmparc_new_fus_cor_resDTI.mgz' ' ' new_dir '/' 'wmparc_new_fus_cor_resDTI.mgz' ]);
    
    % Extract new wmparc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vol_dta = ejk_load_vol( [ dti_dir '/' dti_sbj '/' 'DTcalc' '/' 'DTI1_crev_corr_regT1_DT_FA.mgz' ] );
    
    seg_dta = ejk_load_vol( [ new_dir '/' 'wmparc_new_fus_cor_resDTI.mgz' ] );
      
    col_loc = unique( seg_dta(:) );
    col_loc = col_loc( col_loc>3000 & col_loc<5000 );
   
    wmp_hld{iS,1} = wmp_dta_rcn{grp.tle_post_3T_ATLonly_left(iS)};
    for iFC = 1:size(col_loc,1)
        if iFC ~= 8
            wmp_hld{iS, iFC+1 } = nanmean( nanmean(vol_dta(seg_dta==col_loc(iFC)),1) );
        elseif iFC == 8
            wmp_hld{iS, iFC+1 } = nanmean( nanmean( vol_dta(seg_dta==col_loc(6) | seg_dta==col_loc(8) ),1) );
        end
    end
    
    % Clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unix(['rm' ' ' '-r' ' ' avg_dir '/' idv_sbj '/'])
    
end

wmp_hld( cellfun(@isempty,wmp_hld)) = {NaN};

[wmp_hld(:,1) cog_dta_sbj(grp.tle_post_3T_ATLonly_left,1)]

cell2csv( '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/WMPARC_TRY/NewWmparcHold/new_roi.csv', [ 'SubjID' 'recon' prc_lbl ; cog_dta_sbj(grp.tle_post_3T_ATLonly_left,1) wmp_hld ] );



%%



%%

















