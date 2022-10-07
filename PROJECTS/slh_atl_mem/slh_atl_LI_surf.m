clear; clc;

slh_atl_mem_constants

out_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/surf/';

dta_loc = '/home/mmilmcd2/data/MCD_MRI/proc_dti/';
dta_fld = 'DTanalysis';

fsr_avg_dir = '/space/syn02/1/data/MMILDB/BACKUP/md7/1/pubsw/packages/freesurfer/RH4-x86_64-R711/subjects/';

lhs_str_reg_one = [ fsr_avg_dir '/' 'fsaverage'     '/'             'surf' '/' 'lh.fsaverage_sym.sphere.reg' ];
lhs_str_reg_two = [ fsr_avg_dir '/' 'fsaverage_sym' '/'             'surf' '/' 'lh.sphere.reg' ];
rhs_str_reg_one = [ fsr_avg_dir '/' 'fsaverage'     '/' 'xhemi' '/' 'surf' '/' 'lh.fsaverage_sym.sphere.reg' ];
rhs_str_reg_two = [ fsr_avg_dir '/' 'fsaverage_sym' '/'             'surf' '/' 'lh.sphere.reg' ];

smt_stp = { '16' '256' };
grp_nme = { 'ltle_slah' 'ltle_atl' };

%% Setup
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'])

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ ~, cln_dta_sbj, ~] = ejk_dta_frm( fcfg );

rcn_fle = mmil_readtext('/home/ekaestne/gitrep/MMIL/EXTERNAL/McD/mmilmcdRSI_freesurfer_recons.csv');

%% Flip L/R
sbj_ind = sort([ grp.surgery.pst_cog_dti.ltle_atl ; grp.surgery.pst_cog_dti.ltle_slah ]);

dti_prc_dir = dir(dta_loc); dti_prc_dir = {dti_prc_dir(:).name};

for iSM = 1:numel(smt_stp)
    for iS = 1:numel(sbj_ind)
        
        fsr_ind = find(strcmpi( rcn_fle(:,1), cln_dta_sbj{sbj_ind(iS)}));
        dti_prc_ind = string_find(dti_prc_dir,['DTIPROC' '_' rcn_fle{fsr_ind,3}]);
        if numel(dti_prc_ind)>1; error(''); end
        
        lhs_fle = [ dta_loc '/' dti_prc_dir{dti_prc_ind}  '/' dta_fld '/' 'FA_gwcsurf_wm-sphere-sm' smt_stp{iSM} '-lh.mgz' ];
        lhs_fle_out = [ out_dir '/' cln_dta_sbj{sbj_ind(iS)} '_' 'left_on_left_sm' smt_stp{iSM} '.mgz' ];
        rhs_fle = [ dta_loc '/' dti_prc_dir{dti_prc_ind}  '/' dta_fld '/' 'FA_gwcsurf_wm-sphere-sm' smt_stp{iSM} '-rh.mgz' ];
        rhs_fle_out = [ out_dir '/' cln_dta_sbj{sbj_ind(iS)} '_' 'right_on_left_sm' smt_stp{iSM} '.mgz' ];
        
        lhs_cmd = [ 'mris_apply_reg' ' ' ...
            '--src' ' ' lhs_fle ' ' ...
            '--trg' ' ' lhs_fle_out ' ' ...
            '--streg' ' ' lhs_str_reg_one ' ' lhs_str_reg_two ];
        unix(lhs_cmd)
        
        rhs_cmd = [ 'mris_apply_reg' ' ' ...
            '--src' ' ' rhs_fle ' ' ...
            '--trg' ' ' rhs_fle_out ' ' ...
            '--streg' ' ' rhs_str_reg_one ' ' rhs_str_reg_two ];
        unix(rhs_cmd)
        
    end
end

%% Load and Calculate LI
ref_dta = mmil_rowvec(fs_load_mgh([ out_dir '/' cln_dta_sbj{sbj_ind(1)} '_' 'left_on_left_sm' smt_stp{1} '.mgz' ]));

for iSM = 1:numel(smt_stp)
    
    lhs_gwc_srf = nan(numel(cln_dta_sbj),size(ref_dta,2));
    rhs_gwc_srf = nan(numel(cln_dta_sbj),size(ref_dta,2));
    lat_gwc_srf = nan(numel(cln_dta_sbj),size(ref_dta,2));
    
    for iS = 1:numel(sbj_ind)
        
        lhs_gwc_srf(sbj_ind(iS),:) = mmil_rowvec(fs_load_mgh([ out_dir '/' cln_dta_sbj{sbj_ind(iS)} '_' 'left_on_left_sm' smt_stp{iSM} '.mgz' ]));
        rhs_gwc_srf(sbj_ind(iS),:) = mmil_rowvec(fs_load_mgh([ out_dir '/' cln_dta_sbj{sbj_ind(iS)} '_' 'right_on_left_sm' smt_stp{iSM} '.mgz' ]));
        
        lat_gwc_srf(sbj_ind(iS),:) = ( lhs_gwc_srf(sbj_ind(iS),:) - rhs_gwc_srf(sbj_ind(iS),:) ) ./ ...
            ( lhs_gwc_srf(sbj_ind(iS),:) + rhs_gwc_srf(sbj_ind(iS),:) ) ;
        
    end
    
    srf_dta_sbj = cln_dta_sbj;
    srf_dta     = lhs_gwc_srf;
    save([ out_dir '/' 'left_on_left_gwc_sm' smt_stp{iSM} '.mat' ] ,'srf_dta','srf_dta_sbj');
    
    srf_dta     = rhs_gwc_srf;
    save([ out_dir '/' 'right_on_left_gwc_sm' smt_stp{iSM} '.mat' ],'srf_dta','srf_dta_sbj');
    
    srf_dta     = lat_gwc_srf;
    save([ out_dir '/' 'laterality_gwc_sm' smt_stp{iSM} '.mat' ]   ,'srf_dta','srf_dta_sbj');
end

%% Correlate LM2
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj,cog_dta_col] = ejk_dta_frm( fcfg );

% Surface Correlations
for iG = 1:numel(grp_nme)
    for iSM = 1:numel(smt_stp)
        
        smt_stp_use = str2double(smt_stp{iSM});
        pvl_chs = .05;
        pvl_cls = .05;
        
        fcfg = [];
        
        fcfg.smt_stp = smt_stp_use;
        fcfg.pvl_chs = pvl_chs;
        fcfg.pvl_cls = pvl_cls;
        
        fcfg.sbj_nme = cog_dta_sbj( grp.surgery.pst_cog_dti.(grp_nme{iG}), 1);
        
        fcfg.dta_lhs = [ out_dir '/' 'laterality_gwc_sm' smt_stp{iSM} '.mat']; %
        fcfg.dta_rhs = [ out_dir '/' 'laterality_gwc_sm' smt_stp{iSM} '.mat']; %
        
        fcfg.cor     = cell2mat(cog_dta( grp.surgery.pst_cog_dti.(grp_nme{iG}), strcmpi(cog_dta_col,'lm2_chg')) );
        fcfg.cor_nme = cog_dta_col( 1, strcmpi(cog_dta_col,'lm2_chg'));
        
        fcfg.grp     = repmat({grp_nme{iG}},numel(fcfg.sbj_nme),1);
        fcfg.grp_nme = {[ (grp_nme{iG}) '_grp']};
        fcfg.grp_cmp = {(grp_nme{iG})};
        
        fcfg.cov     = [];
        fcfg.cov_nme = [];
        
        fcfg.out_dir = out_dir;
        fcfg.out_pre = [ 'lm2_chg' '_' grp_nme{iG} '_' smt_stp{iSM}];
        
        ejk_surface_correlations_spearman( fcfg );
        
    end
end

%% Correlate LM1
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cog_LM1.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj,cog_dta_col] = ejk_dta_frm( fcfg );

% Surface Correlations
for iG = 1:numel(grp_nme)
    for iSM = 1:numel(smt_stp)
        
        smt_stp_use = str2double(smt_stp{iSM});
        pvl_chs = .05;
        pvl_cls = .05;
        
        fcfg = [];
        
        fcfg.smt_stp = smt_stp_use;
        fcfg.pvl_chs = pvl_chs;
        fcfg.pvl_cls = pvl_cls;
        
        fcfg.sbj_nme = cog_dta_sbj( grp.surgery.pst_cog_dti.(grp_nme{iG}), 1);
        
        fcfg.dta_lhs = [ out_dir '/' 'laterality_gwc_sm' smt_stp{iSM} '.mat']; %
        fcfg.dta_rhs = [ out_dir '/' 'laterality_gwc_sm' smt_stp{iSM} '.mat']; %
        
        fcfg.cor     = cell2mat(cog_dta( grp.surgery.pst_cog_dti.(grp_nme{iG}), strcmpi(cog_dta_col,'lm1_chg')) );
        fcfg.cor_nme = cog_dta_col( 1, strcmpi(cog_dta_col,'lm1_chg'));
        
        fcfg.grp     = repmat({grp_nme{iG}},numel(fcfg.sbj_nme),1);
        fcfg.grp_nme = {[ (grp_nme{iG}) '_grp']};
        fcfg.grp_cmp = {(grp_nme{iG})};
        
        fcfg.cov     = [];
        fcfg.cov_nme = [];
        
        fcfg.out_dir = out_dir;
        fcfg.out_pre = [ 'lm1_chg' '_' grp_nme{iG} '_' smt_stp{iSM}];
        
        ejk_surface_correlations_spearman( fcfg );
        
    end
end
