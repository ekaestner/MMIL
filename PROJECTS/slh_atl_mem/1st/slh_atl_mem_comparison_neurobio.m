out_put = [ prj_dir '/' prj_nme '/' 'InitialAnalysis'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% Group Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subset ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roi_sub_ind = 0;

roi_sub_mse = { 'cort_thick_ctx' 'fiber_FA' 'subcort_vol' 'wmparc_FA_wm' };
roi_sub_nme = { { 'entorhinal' 'parahippocampal' } ...
                { 'Unc' 'ILF' } ...
                { 'hippocampus' 'Thalamus_Proper' } ...
                { 'entorhinal' 'parahippocampal' } };
            
% Subset Neural Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_sub_ind = 0;

dta_sub_mse = { 'subcort_vol_238_LateralityIndex'              'subcort_vol_dev_LateralityIndex' ...
                'wmparc_FA_wm_aparc_annot_238_LateralityIndex' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex' ...
                'fiber_FA_238_LateralityIndex' 'fiber_FA_dev_LateralityIndex' };

% Group Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_grp = { { 'controls_pre_3T_allSurg_all'   'tle_pre_3T_allSurg_left'   'tle_pre_3T_allSurg_right'  } ...
            { 'controls_pre_allT_allSurg_all' 'tle_pre_allT_allSurg_left' 'tle_pre_allT_allSurg_right' }...
            { 'tle_pre_3T_allSurg_left'       'tle_pre_3T_allSurg_right' } ...
            { 'tle_pre_allT_allSurg_left'     'tle_pre_allT_allSurg_right' }...
            { 'tle_post_3T_ATLonly_left'      'tle_post_3T_ATLonly_right' } ...
            { 'tle_post_allT_ATLonly_left'    'tle_post_allT_ATLonly_right' } };
cmp_nme = { { 'HC_pre_3T'                     'LTLE_pre_3T'               'RTLE_pre_3T' } ...
            { 'HC_pre_allT'                   'LTLE_pre_allT'             'RTLE_pre_allT'} ...
            { 'LTLE_pre_3T'                   'RTLE_pre_3T' } ...
            { 'LTLE_pre_allT'                 'RTLE_pre_allT'  } ...
            { 'LTLE_post_3T'                  'RTLE_post_3T' } ...
            { 'LTLE_post_allT'                'RTLE_post_allT' } };
cmp_out = { 'TLE_Controls_pre_3T' ...
            'TLE_Controls_pre_allT' ...
            'TLE_Only_pre_3T' ...
            'TLE_Only_pre_allT' ...
            'TLE_Only_post_3T' ...
            'TLE_Only_post_allT' };     
       
ttb_use = [ 3 4 5 6 ];
anv_use = [ 1 2 ];    

% Correlation Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmp_out_cor = { 'TLE_Controls_pre_cog_3T'  ...
                ...'LTLE_Controls_pre_cog_3T' ...
                ...'RTLE_Controls_pre_cog_3T' ...
                'LTLE_post_cog_3T'         ...  
                'RTLE_post_cog_3T'         ...
                ...'TLE_post_cog_3T' ...
                'TLE_Controls_pre_cog_allT'  ...
                ...'LTLE_Controls_pre_cog_allT' ...
                ...'RTLE_Controls_pre_cog_allT' ...
                'LTLE_post_cog_allT'         ...  
                'RTLE_post_cog_allT'         ...
                };... 'TLE_post_cog_allT' };

cmp_nme_cor = { 'tle_controls_pre_3T_allSurg_all' ... 
                ...'tle_controls_pre_3T_allSurg_left'  ...
                ...'tle_controls_pre_3T_allSurg_right' ...
                'tle_post_3T_ATLonly_left'  ...  
                'tle_post_3T_ATLonly_right' ...
                ...'tle_post_3T_ATLonly_all' ...
                'tle_controls_pre_allT_allSurg_all' ... 
                ...'tle_controls_pre_allT_allSurg_left'  ...
                ...'tle_controls_pre_allT_allSurg_right' ...
                'tle_post_allT_ATLonly_left'  ...  
                'tle_post_allT_ATLonly_right' ...
                };...'tle_post_allT_ATLonly_all' };

% Correlation w/ Clinical Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
age_col = find(strcmpi(cln_dta_col,'AgeAtImaging'));
edu_col = find(strcmpi(cln_dta_col,'Educ'));
aos_col = find(strcmpi(cln_dta_col,'AgeOfSeizureOnset'));
asm_col = find(strcmpi(cln_dta_col,'NumAEDs'));
sze_col = find(strcmpi(cln_dta_col,'SeizureFreq'));

cln_col_one = repmat( {[ age_col edu_col aos_col asm_col sze_col ]}, numel(cmp_out_cor) ,1 );

% Correlation w/ Cognitive Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_tst_lm2_pre = find(strcmpi(cog_dta_col,'log_mem_nor_scr_two'));
cog_tst_vpa_pre = find(strcmpi(cog_dta_col,'vp2_nor_scr'));
cog_tst_lm2_pst = find(strcmpi(cog_dta_col,'log_mem_nor_scr_two_pst'));
cog_tst_vpa_pst = find(strcmpi(cog_dta_col,'vp2_nor_scr_pst'));
            
cog_col_one = { [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                ...[ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                ...[ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pst cog_tst_vpa_pst ] ...
                [ cog_tst_lm2_pst cog_tst_vpa_pst ] ...
                ...[ cog_tst_lm2_pst cog_tst_vpa_pst ] ...
                [ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                ...[ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                ...[ cog_tst_lm2_pre cog_tst_vpa_pre ] ...
                [ cog_tst_lm2_pst cog_tst_vpa_pst ] ...
                [ cog_tst_lm2_pst cog_tst_vpa_pst ] ...
                }; ...[ cog_tst_lm2_pst cog_tst_vpa_pst ] };

%% Load
inp_fle     = { mri_238_fle dti_238_fle mri_dev_fle dti_dev_fle };
inp_pre     = { mri_238_pre dti_238_pre mri_dev_pre dti_dev_pre };
inp_suf     = { mri_238_suf dti_238_suf mri_dev_suf dti_dev_suf };
inp_mse     = { mri_238_mse dti_238_mse mri_dev_mse dti_dev_mse };
inp_roi     = { mri_238_roi dti_238_roi mri_dev_roi dti_dev_roi };
inp_lat     = { mri_238_lat dti_238_lat mri_dev_lat dti_dev_lat };
inp_icv     = { mri_238_icv dti_238_icv mri_dev_icv dti_dev_icv };

lod_cnt = 1;
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                nme_hld = [inp_mse{iF}{iM} '_' inp_suf{iF}{iM}];
            else
                roi_fle_nme = [prj_dir '/' prj_nme '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                nme_hld = [inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM}];
            end
                        
            fcfg = [];
            fcfg.dta_loc = [roi_fle_nme(1:end-4) '_QC.csv'];
            fcfg.dta_col = 5;
            [ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );
            
            dta_hld{lod_cnt} = neu_dta;
            dta_sbj{lod_cnt} = neu_dta_sbj;
            dta_col{lod_cnt} = neu_dta_col;
            dta_nme{lod_cnt} = nme_hld;
            
            lod_cnt = lod_cnt + 1;
            
            if ~strcmpi(inp_icv{iF}{iM},'')
                fle_nme = [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '_QC.csv'];
                
                fcfg = [];
                fcfg.dta_loc = fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );
                
                dta_hld{lod_cnt} = neu_dta;
                dta_sbj{lod_cnt} = neu_dta_sbj;
                dta_col{lod_cnt} = neu_dta_col;
                dta_nme{lod_cnt} = [nme_hld '_' 'norm' '_' inp_icv{iF}{iM}];
                
                lod_cnt = lod_cnt + 1;
            end
            
            if inp_lat{iF}(iM)==1
                fle_nme = [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '_QC.csv'];
                
                fcfg = [];
                fcfg.dta_loc = fle_nme;
                fcfg.dta_col = 5;
                [ neu_dta, neu_dta_sbj, neu_dta_col] = ejk_dta_frm( fcfg );
                
                dta_hld{lod_cnt} = neu_dta;
                dta_sbj{lod_cnt} = neu_dta_sbj;
                dta_col{lod_cnt} = neu_dta_col;
                dta_nme{lod_cnt} = [nme_hld '_' 'LateralityIndex'];
                
                lod_cnt = lod_cnt + 1;
            end
            
        end
    end
end

if dta_sub_ind
    
    [ ~, dta_ind, ~ ] = intersect( dta_nme, dta_sub_mse);
    
    dta_hld = dta_hld(dta_ind);
    dta_sbj = dta_sbj(dta_ind);
    dta_col = dta_col(dta_ind);
    dta_nme = dta_nme(dta_ind);
    
end

%% t-tests independent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iD = 1:numel(dta_hld) 
                
    for iC = 1:numel(ttb_use)

        if roi_sub_ind==1
            dta_typ_use = { { } ...
                            { } ...
                            dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) ...
                            dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) ...
                            dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) ...
                            dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) };
            out_dir = [ out_put '/' 'Neurobio_compare' '/' 'ttest2_subset' '/' dta_nme{iD} '/' cmp_out{ttb_use(iC)} ];
        else
            dta_typ_use = { { } ...
                            { } ...
                            dta_col{iD} ...
                            dta_col{iD} ...
                            dta_col{iD} ...
                            dta_col{iD} };
            out_dir = [ out_put '/' 'Neurobio_compare' '/' 'ttest2' '/' dta_nme{iD} '/' cmp_out{ttb_use(iC)} ];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use{ttb_use(iC)} );
        
        fcfg = [];
        fcfg.grp     = grp;
        fcfg.grp_inc = cmp_grp(ttb_use(iC));
        fcfg.grp_nme = cmp_nme(ttb_use(iC));
        fcfg.dta = dta_hld{iD}(:,use_dta_col);
        fcfg.sbj = dta_sbj{iD};
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        % Test
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = dta_col{iD}(:,use_dta_col);
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = cmp_out(1);
        fcfg.out_dir = out_dir;
        ejk_ttest2_independent( fcfg );
        
    end
end

%% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iD = 1:numel(dta_hld)
                               
    for iC = 1:numel(anv_use)
        
        if roi_sub_ind==1
            dta_typ_use = { dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) ...
                            dta_col{iD}( string_find(dta_col{iD},roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))})) ...
                            { } ...
                            { } ...
                            { } ...
                            { } };
            out_dir = [ out_put '/' 'Neurobio_compare' '/' 'ANOVA_subset' '/' dta_nme{iD} '/' cmp_out{anv_use(iC)}];
        else
            dta_typ_use = { dta_col{iD} ...
                            dta_col{iD} ...
                            {} ...
                            {} ...
                            {} ...
                            {} };
            out_dir = [ out_put '/' 'Neurobio_compare' '/' 'ANOVA' '/' dta_nme{iD} '/' cmp_out{anv_use(iC)}];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use{anv_use(iC)} );
        
        fcfg = [];
        fcfg.grp     = grp;
        fcfg.grp_inc = cmp_grp(anv_use(iC));
        fcfg.grp_nme = cmp_nme(anv_use(iC));
        fcfg.dta = dta_hld{iD}(:,use_dta_col);
        fcfg.sbj = dta_sbj{iD};
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = dta_col{iD}(use_dta_col);
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = cmp_out(1);
        fcfg.out_dir = out_dir;
        ejk_1way_anova( fcfg )
        
    end
end

%% Correlation comparing Neurobio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iG = 1:numel(cmp_out_cor)
% 
%     for iD = 1:numel(dta_hld)/2
%         
%         % Get columns of interest
%         [ ~, ind_one, ind_two ] = intersect( dta_col{iD}, dta_col{iD+numel(dta_hld)/2});
%         
%         % Correlation
%         fcfg = [];
%         
%         fcfg.sbj_nme = dta_sbj{iD}( grp.(cmp_nme_cor{iG}), 1);
%         
%         fcfg.dta_one = cell2mat(dta_hld{iD}( grp.(cmp_nme_cor{iG}), ind_one));
%         fcfg.lbl_one = dta_col{iD}(ind_one);
%         
%         fcfg.cor_typ = 'spearman';
%         
%         fcfg.dta_two = cell2mat(dta_hld{iD+numel(dta_hld)/2}( grp.(cmp_nme_cor{iG}), ind_two));
%         fcfg.lbl_two = dta_col{iD+numel(dta_hld)/2}(ind_two);
%         
%         fcfg.pvl_cut = 0.05;
%         fcfg.pvl_lib = 0.20;
%         
%         fcfg.force_plot = 1;
%         
%         fcfg.out_dir = [ out_put '/' 'Neurobio_compare' '/' 'Correlation' '/' 'Compare' '/' cmp_out_cor{iG} '/' strrep(dta_nme{iD},'_238','') '/' ];
%         
%         ejk_cross_cor( fcfg );
%         
%     end
% end

%% Correlation of Neurobio w/ Clinical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out_cor)

    for iD = 1:numel(dta_hld)

        if roi_sub_ind==1
            dta_typ_use = roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))};
            dta_typ_use = dta_col{iD}( string_find(dta_col{iD},dta_typ_use) );
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation_subset' '/' 'Clinical' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        else
            dta_typ_use = dta_col{iD};
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation' '/' 'Clinical' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use );
        
        fcfg = [];
        
        fcfg.sbj_nme = dta_sbj{iD}( grp.(cmp_nme_cor{iG}), 1);
        
        fcfg.dta_one = cell2mat(dta_hld{iD}( grp.(cmp_nme_cor{iG}), :));
        fcfg.lbl_one = dta_col{iD}(use_dta_col);
        
        fcfg.cor_typ = 'pearson';
        
        fcfg.dta_two = cell2mat(cln_dta( grp.(cmp_nme_cor{iG}), cln_col_one{iG}));
        fcfg.lbl_two = cln_dta_col(cln_col_one{iG});
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.force_plot = 1;
        
        fcfg.out_dir = out_dir;
        
        ejk_cross_cor( fcfg );
        
    end
end

%% Correlation of Neurobio w/ Cognitive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out_cor)

    for iD = 1:numel(dta_hld)
        
        if roi_sub_ind==1
            dta_typ_use = roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))};
            dta_typ_use = dta_col{iD}( string_find(dta_col{iD},dta_typ_use) );
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation_subset' '/' 'Cognitive' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        else
            dta_typ_use = dta_col{iD};
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation' '/' 'Cognitive' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use );
        
        fcfg = [];
        
        fcfg.sbj_nme = dta_sbj{iD}( grp.(cmp_nme_cor{iG}), 1);
        
        fcfg.dta_two = cell2mat(dta_hld{iD}( grp.(cmp_nme_cor{iG}), use_dta_col));
        fcfg.lbl_two = dta_col{iD}(use_dta_col);
        
        fcfg.cor_typ = 'pearson';
                
        fcfg.dta_one = cell2mat(cog_dta( grp.(cmp_nme_cor{iG}), cog_col_one{iG}));
        fcfg.lbl_one = cog_dta_col(cog_col_one{iG});
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.force_plot = 1;
        
        fcfg.out_dir = out_dir;
        
        ejk_cross_cor( fcfg );
        
    end
end

%% Correlation of Neurobio w/ Clinical - SPEARMAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out_cor)

    for iD = 1:numel(dta_hld)

        if roi_sub_ind==1
            dta_typ_use = roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))};
            dta_typ_use = dta_col{iD}( string_find(dta_col{iD},dta_typ_use) );
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation_subset' '/' 'Clinical_spearman' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        else
            dta_typ_use = dta_col{iD};
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation' '/' 'Clinical_spearman' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use );
        
        fcfg = [];
        
        fcfg.sbj_nme = dta_sbj{iD}( grp.(cmp_nme_cor{iG}), 1);
        
        fcfg.dta_one = cell2mat(dta_hld{iD}( grp.(cmp_nme_cor{iG}), :));
        fcfg.lbl_one = dta_col{iD}(use_dta_col);
        
        fcfg.cor_typ = 'spearman';
        
        fcfg.dta_two = cell2mat(cln_dta( grp.(cmp_nme_cor{iG}), cln_col_one{iG}));
        fcfg.lbl_two = cln_dta_col(cln_col_one{iG});
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.force_plot = 1;
        
        fcfg.out_dir = out_dir;
        
        ejk_cross_cor( fcfg );
        
    end
end

%% Correlation of Neurobio w/ Cognitive - SPEARMAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(cmp_out_cor)

    for iD = 1:numel(dta_hld)
        
        if roi_sub_ind==1
            dta_typ_use = roi_sub_nme{ ~cellfun(@isempty,cellfun(@(x) strfind( dta_nme{iD}, x ), roi_sub_mse, 'uni', 0))};
            dta_typ_use = dta_col{iD}( string_find(dta_col{iD},dta_typ_use) );
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation_subset' '/' 'Cognitive_spearman' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        else
            dta_typ_use = dta_col{iD};
            out_dir     = [ out_put '/' 'Neurobio_compare' '/' 'Correlation' '/' 'Cognitive_spearman' '/' cmp_out_cor{iG} '/' dta_nme{iD} '/' ];
        end
        
        [~, use_dta_col ] = intersect( dta_col{iD}, dta_typ_use );
        
        fcfg = [];
        
        fcfg.sbj_nme = dta_sbj{iD}( grp.(cmp_nme_cor{iG}), 1);
        
        fcfg.dta_two = cell2mat(dta_hld{iD}( grp.(cmp_nme_cor{iG}), use_dta_col));
        fcfg.lbl_two = dta_col{iD}(use_dta_col);
        
        fcfg.cor_typ = 'spearman';
                
        fcfg.dta_one = cell2mat(cog_dta( grp.(cmp_nme_cor{iG}), cog_col_one{iG}));
        fcfg.lbl_one = cog_dta_col(cog_col_one{iG});
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.force_plot = 1;
        
        fcfg.out_dir = out_dir;
        
        ejk_cross_cor( fcfg );
        
    end
end














