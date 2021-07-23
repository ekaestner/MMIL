load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' '/' ]; 

out_put = [ prj_dir '/' prj_nme '/' 'SurfaceCorrelation' '/'];
    ejk_chk_dir( out_put );

run_grp = fieldnames(grp);

hms_hld = { 'lhs' 'rhs' };

cmp_out = { 'TLE_Controls_pre_pre' ... 
            'LTLE_post_post' ...  
            'RTLE_post_post' };

cmp_nme = { 'tle_controls_pre_3T_allSurg_all' ... 
            'tle_post_3T_ATLonly_left'  ...  
            'tle_post_3T_ATLonly_right' };
 
smt_stp = 313;
pvl_chs = .05;
pvl_cls = .05;   

%% WMPARC
cln_dta = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ]);
    cln_dta_col = cln_dta(1,2:end);
    cln_dta_sbj = cln_dta(2:end,1);
    cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
    cog_dta_col = cog_dta(1,2:end);
    cog_dta_sbj = cog_dta(2:end,1);
    cog_dta     = cog_dta(2:end,2:end);

dta_hld = mmil_readtext( [prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC' '.csv']);
    dta_sbj = dta_hld(2:end,1);
    dta_lbl = dta_hld(1,5:end);
    dta_hld = cell2mat(dta_hld(2:end,5:end));
    
% Plot
for iG = 1:numel(cmp_nme)
    for iT = 1:2%numel(cog_col.(cmp_nme{iG}))
        
        % Surface Correlations
        fcfg = [];
        
        fcfg.smt_stp = smt_stp;
        fcfg.pvl_chs = pvl_chs;
        fcfg.pvl_cls = pvl_cls;
        
        fcfg.sbj_nme = cln_dta_sbj( grp.(cmp_nme{iG}), 1);
        
        fcfg.dta_lhs = [ dta_dir '/' 'surf_wmparc_fa_lhs_sm' num2str(smt_stp) '.mat']; %
        fcfg.dta_rhs = [ dta_dir '/' 'surf_wmparc_fa_rhs_sm' num2str(smt_stp) '.mat']; %
        
        fcfg.cor     = cell2mat(cog_dta( grp.(cmp_nme{iG}), cog_col.(cmp_nme{iG})(iT)));
        fcfg.cor_nme = cog_dta_col( cog_col.(cmp_nme{iG})(iT));
        
        fcfg.grp     = repmat({cmp_nme{iG}},numel(fcfg.sbj_nme),1);
        fcfg.grp_nme = cmp_out(iG);
        fcfg.grp_cmp = cmp_nme(iG);
        
        fcfg.cov     = [];
        fcfg.cov_nme = [];
        
        fcfg.out_dir = [ out_put '/' cmp_out{iG} '/']; ejk_chk_dir([ out_put '/' cmp_out{iG} '/'])            
        fcfg.out_pre = [ cog_dta_col{ cog_col.(cmp_nme{iG})(iT)} '_' 'novelpatients' ];
        
        ejk_surface_correlations( fcfg );
        
        
    end
end



