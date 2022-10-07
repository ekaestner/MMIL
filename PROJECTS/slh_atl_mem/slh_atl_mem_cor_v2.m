out_dir = [ prj_dir '/' prj_nme '/' 'FinalAnalysis_restricted' '/'];

%% Final Data
% Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

dta_nme = { 'Clinical.csv'                                          2 ; ...
            'Cognitive.csv'                                         2 ; ...
            'subcort_vol_dev_LateralityIndex.csv'                   5 ; ...
            'fiber_FA_dev_LateralityIndex.csv'                      5 ; ...
            'gwcsurf_FA_wm_ctx_aparc_annot_dev_LateralityIndex.csv' 5 ; ...
            'subcort_vol_dev_norm_IntracranialVolume_ComBat.csv'    5 ; ...
            'fiber_FA_dev_ComBat.csv'                               5 ; ...
            'cort_thick_ctx_aparc_annot_dev_LateralityIndex_lbl.csv' 5 ; ...
            'gwcsurf_FA_wm_ctx_aparc_annot_dev_ComBat.csv'          5 };
dta_fld = { { 'Location' 'SurgeryType' 'SideOfSeizureFocus' 'AgeAtSurgery' 'Educ' 'AgeOfSeizureOnset' 'NumAEDs' 'Sex' 'Handedness' 'MTS' 'SeizureFreq' } ...
            { 'lm2_pre' 'lm2_chg' 'lm2_rci' } ...
            { 'Hippocampus' 'Thalamus_Proper' } ...
            { 'Unc' 'ILF' } ...
            { 'lateralorbitofrontal' 'fusiform' } ...
            { 'xLeft_Hippocampus' 'xRight_Hippocampus' } ...
            { 'xL_Unc' 'xR_Unc' 'xL_ILF' 'xR_ILF' } ...
            { 'lateralorbitofrontal_thick' 'fusiform_thick' } ...
            { 'lh_lateralorbitofrontal' 'rh_lateralorbitofrontal' 'lh_fusiform' 'rh_fusiform' } };

fnl_dta     = cell(numel(cln_dta_sbj),sum(cellfun(@numel,dta_fld)));
fnl_dta_sbj = cln_dta_sbj;
fnl_dta_col = cell(1,sum(cellfun(@numel,dta_fld)));
ind_col = 1;
for iD = 1:size(dta_nme,1)
    
    % Load
    fcfg = [];
    fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' dta_nme{iD,1} ];
    fcfg.dta_col = dta_nme{iD,2};
    [ add_dta, add_dta_sbj, add_dta_col] = ejk_dta_frm( fcfg );  
    
    % Add
    fnl_dta(:,ind_col:ind_col+numel(dta_fld{iD})-1)   = add_dta(:,ismember(add_dta_col,dta_fld{iD}));
    fnl_dta_col(ind_col:ind_col+numel(dta_fld{iD})-1) = add_dta_col(ismember(add_dta_col,dta_fld{iD}));
    
    ind_col = ind_col+numel(dta_fld{iD});
end

cell2csv([ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv' ], [ 'sbj_nme' fnl_dta_col ; fnl_dta_sbj fnl_dta ])

%% Final groupings
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

% Add L surgery only
grp.surgery_left.pst_cog_dti.ltle_slah = grp.surgery.pst_cog_dti.ltle_slah;
grp.surgery_left.pst_cog_dti.ltle_atl  = grp.surgery.pst_cog_dti.ltle_atl;

% Add TLE/HC by surgery
grp.diagnosis_surgery.pst_cog_dti.control = grp.diagnosis.pre_cog_dti.control;
grp.diagnosis_surgery.pst_cog_dti.ltle_slah = grp.surgery.pst_cog_dti.ltle_slah;
grp.diagnosis_surgery.pst_cog_dti.ltle_atl  = grp.surgery.pst_cog_dti.ltle_atl;
grp.diagnosis_surgery.pst_cog_dti.rtle_slah = grp.surgery.pst_cog_dti.rtle_slah;
grp.diagnosis_surgery.pst_cog_dti.rtle_atl  = grp.surgery.pst_cog_dti.rtle_atl;

% Add TLE/HC by laterality
grp.diagnosis_surgery_laterality.pst_cog_dti.control = grp.diagnosis.pre_cog_dti.control;
grp.diagnosis_surgery_laterality.pst_cog_dti.ltle    = [ grp.surgery.pst_cog_dti.ltle_slah ; grp.surgery.pst_cog_dti.ltle_atl ];
grp.diagnosis_surgery_laterality.pst_cog_dti.rtle    = [ grp.surgery.pst_cog_dti.rtle_slah ; grp.surgery.pst_cog_dti.rtle_atl ];

%% CHECK FOR NAN/NAN/0 RCI problem
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv' ];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] =  ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'additional_variables.csv'];
fcfg.dta_col = 2;
[ cln_nui_dta, cln_nui_dta_sbj, cln_nui_dta_col] = ejk_dta_frm( fcfg );

% Combine
fnl_dta     = [ fnl_dta     cln_nui_dta];
fnl_dta_col = [ fnl_dta_col cln_nui_dta_col ];

fnl_dta(cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_pre'))) & cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_chg'))),strcmpi(fnl_dta_col,'lm2_rci')) = {NaN};

%% Restrict grp
lm2_chg_col = strcmpi(fnl_dta_col,'lm2_chg');

grp_nme = fieldnames(grp);
for iG = 1:numel(grp_nme)
    grp_sub_nme = fieldnames(grp.(grp_nme{iG}));
    for iGS = 1:numel(grp_sub_nme)
        if ~isempty(strfind(grp_sub_nme{iGS},{'pst_cog'}))
            
            fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_sub_nme{iGS}));
            for iF = 1:numel(fld_nme)
                grp.(grp_nme{iG}).(grp_sub_nme{iGS}).(fld_nme{iF})( cellfun(@isnan,fnl_dta(grp.(grp_nme{iG}).(grp_sub_nme{iGS}).(fld_nme{iF}),lm2_chg_col)) ) = [];
            end
            
        end
    end
end

%% Clinical 
dta_inp{1} = [ fnl_dta_col ; fnl_dta ];

% build table
tbl_typ = { 'count'    'Location'                'HC/UCSD/UCSF/Emory' 1; ...
            'mean/std' 'AgeAtSurgery'            ''         1 ; ...
            'mean/std' 'Educ'                    ''         1 ; ...
            'count'    'Sex'                     'M/F'      1 ; ...
            'count'    'Handedness'              'R/L'      1 ; ...
            'mean/std' 'AgeOfSeizureOnset'       ''         1 ; ...
            'mean/std' 'NumAEDs'                 ''         1 ; ...
            'mean/std' 'sbj_sze_frq'             ''         1 ; ...
            'count'    'sbj_mts'                 'yes/no'   1 ; ...
            'count'    'SideOfSeizureFocus'      'R/L'      1 };

% stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_nme     = {   'surgery'      'surgery_left' };
grp_nme_sub = { { 'pst_cog_dti'} {  'pst_cog_dti' } };

cln_var = tbl_typ(:,2);
col_ind = find(ismember(fnl_dta_col,cln_var));

for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        fsh_var = cellfun(@ischar,fnl_dta(1,col_ind));
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        fcfg.grp_inc = {fld_nme};
        fcfg.grp_nme = {fld_nme};
        fcfg.dta = fnl_dta(:,col_ind(anv_var));
        fcfg.sbj = fnl_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        if numel(fld_nme)>2
            
            % ANOVA            
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_1way_anova( fcfg )
            
            % 2-way ANOVA            
            grp_typ_hld = regexp(grp_typ{1},'_','split');
            grp_typ{1}  = cat(1,grp_typ_hld{:});
            
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp_one     = grp_typ{1}(:,1);
            fcfg.grp_nme_one = {'side'};
            fcfg.grp_two     = grp_typ{1}(:,2);
            fcfg.grp_nme_two = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_2way_anova( fcfg )
            
        elseif numel(fld_nme)==2
            
            % t-test
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical_ttest2' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_ttest2_independent( fcfg );
            
        end
        
        % Fishers     
        [~, lvl_ind] = ismember(cln_var,fnl_dta_col(col_ind(fsh_var)));
        clear cln_lvl; for iLI = 1:numel(lvl_ind); if lvl_ind(iLI); cln_lvl(lvl_ind(iLI)) = tbl_typ(iLI,3); end; end
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        fcfg.grp_inc = {fld_nme};
        fcfg.grp_nme = {fld_nme};
        fcfg.dta = fnl_dta(:,col_ind(fsh_var));
        fcfg.sbj = fnl_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        grp_dta{1}(cellfun(@isempty,grp_dta{1})) = {'N/A'};
        
        fcfg = [];
        fcfg.sbj = grp_sbj{1};
        fcfg.grp_nme = grp_nme{iG};
        fcfg.dta_one = grp_dta{1};
        fcfg.lbl_one = fnl_dta_col(col_ind(fsh_var));
        fcfg.lvl     = cln_lvl;
        fcfg.dta_two = repmat(grp_typ{1},1,sum(fsh_var));
        fcfg.lbl_two = strcat( 'group_', fnl_dta_col(col_ind(fsh_var)));
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical_Fishers' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
        ejk_fisher_test( fcfg );
        
    end
end

% Build 3 Clinical tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_inp{1} = [ fnl_dta_col ; fnl_dta ];
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));        
        if numel(fld_nme)>2
            dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' grp_nme{iG} '/' 'output_table.csv']);
            dta_inp{4} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' 'side' '/' 'output_table.csv']);
        elseif numel(fld_nme)==2
            dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical_ttest2' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
        end
        dta_inp{3} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical_Fishers'   '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' grp_nme{iG} '/' 'output_table.csv']);
                            
        % Create Table
        fcfg = [];
        
        for iR = 1:size(tbl_typ)
            
            if strcmpi(tbl_typ{iR,1},'mean/std')
                for iN = 1:numel(fld_nme)
                    fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
                end
                if numel(fld_nme)>2
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                    fcfg.tbl{iR,iN+2} = ['copy,4,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_one'];
                    fcfg.tbl{iR,iN+3} = ['copy,4,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_two'];
                    fcfg.tbl{iR,iN+4} = ['copy,4,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_int'];
                elseif numel(fld_nme)==2
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                end
            elseif  strcmpi(tbl_typ{iR,1},'count')
                for iN = 1:numel(fld_nme)
                    fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
                end
                fcfg.tbl{iR,iN+1} = ['copy,3,' tbl_typ{iR,2} ',report'];
                if numel(fld_nme)>2
                    fcfg.tbl{iR,iN+2} = 'empty';
                    fcfg.tbl{iR,iN+3} = 'empty';
                    fcfg.tbl{iR,iN+4} = 'empty';
                end
            end
        end
        
        fcfg.dta = dta_inp;
        fcfg.grp = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        tbl_out  = ejk_create_table( fcfg );
        
        % Add-ons
        num_sbj{1} = 'N';
        tbl_lbl{1} = '';
        for iN = 1:numel(fld_nme)
            num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iN}));
            tbl_lbl{iN+1} = fld_nme{iN};
        end
        num_sbj{iN+2} = '-';
        tbl_lbl{iN+2} = 'Stats';
        if numel(fld_nme)>2
        num_sbj{iN+3} = '-';num_sbj{iN+4} = '-';num_sbj{iN+5} = '-';
        tbl_lbl{iN+3} = 'Side';
        tbl_lbl{iN+4} = 'Surgery';
        tbl_lbl{iN+5} = 'Side*Surgery';
        end      
        % Out
        tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
        
        cell2csv([ out_dir '/' 'tables' '/' 'Table1'  '/' 'Clinical' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
        
        clear tbl_lbl num_sbj
        
    end
end

%% Cognitive
% build table
tbl_typ = { 'mean/std' 'lm2_pre' '' 1; ...
            'mean/std' 'lm2_chg' '' 1 ; ...
            'mean/std' 'lm2_rci' '' 1};

% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_nme     = {   'surgery'      'surgery_left' };
grp_nme_sub = { { 'pst_cog_dti'} {  'pst_cog_dti' } };

cog_var = tbl_typ(:,2);
col_ind = find(ismember(fnl_dta_col,cog_var));

for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        % get numeric group data
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        fcfg.grp_inc = {fld_nme};
        fcfg.grp_nme = {fld_nme};
        fcfg.dta = fnl_dta(:,col_ind(anv_var));
        fcfg.sbj = fnl_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        if numel(fld_nme)>2
            
            % 1-way ANOVA
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_1way_anova( fcfg )
            
            % 2-way ANOVA
            grp_typ_hld = regexp(grp_typ{1},'_','split');
            grp_typ{1}  = cat(1,grp_typ_hld{:});
            
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp_one     = grp_typ{1}(:,1);
            fcfg.grp_nme_one = {'side'};
            fcfg.grp_two     = grp_typ{1}(:,2);
            fcfg.grp_nme_two = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_2way_anova( fcfg )
            
        elseif numel(fld_nme)==2
            
            % t-test
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = grp_nme(iG);
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_ttest2' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/'  ];
            ejk_ttest2_independent( fcfg );
            
        end
        
    end
end

% 1-sample ttest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        % get numeric group data
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));

        % ANOVA
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        
        for iGN = 1:numel(fld_nme)
            fcfg = [];
            fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
            fcfg.grp_inc = {fld_nme(iGN)};
            fcfg.grp_nme = {fld_nme(iGN)};
            fcfg.dta = fnl_dta(:,col_ind(anv_var));
            fcfg.sbj = fnl_dta_sbj;
            [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
            
            % Test
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(:,col_ind(anv_var));
            fcfg.men     = 0;
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_ttest1' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' fld_nme{iGN} '/' ];
            ejk_ttest1( fcfg );
        end
        
    end
end

% Build Cognitive tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})

        % Create comparison Table
        dta_inp{1} = [ fnl_dta_col ; fnl_dta ];
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        if numel(fld_nme)>2
            dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' grp_nme{iG} '/' 'output_table.csv']);
            dta_inp{3} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' 'side' '/' 'output_table.csv']);
        elseif numel(fld_nme)==2
            dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_ttest2' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
        end
                
        fcfg = [];
        
        for iR = 1:size(tbl_typ)
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
            end
                if numel(fld_nme)>2
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                    fcfg.tbl{iR,iN+2} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_one'];
                    fcfg.tbl{iR,iN+3} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_two'];
                    fcfg.tbl{iR,iN+4} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_int'];
                elseif numel(fld_nme)==2
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                end
        end
        
        fcfg.dta = dta_inp;
        fcfg.grp = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        tbl_out  = ejk_create_table( fcfg );
        
        % Add-ons
        num_sbj{1} = 'N';
        tbl_lbl{1} = '';
        for iN = 1:numel(fld_nme)
            num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iN}));
            tbl_lbl{iN+1} = fld_nme{iN};
        end
        num_sbj{iN+2} = '-';
        tbl_lbl{iN+2} = 'Stats';
        if numel(fld_nme)>2
            num_sbj{iN+3} = '-';num_sbj{iN+4} = '-';num_sbj{iN+5} = '-';
            tbl_lbl{iN+3} = 'Side';
            tbl_lbl{iN+4} = 'Surgery';
            tbl_lbl{iN+5} = 'Side*Surgery';
        end
        
        % Out
        tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
        
        cell2csv([ out_dir '/' 'tables' '/' 'Table2'  '/' 'Cognitive' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
        
        clear tbl_lbl num_sbj dta_inp
        
        % Create diff from 0 Table
        for iD = 1:numel(fld_nme)
            dta_inp{iD} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_ttest1' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' fld_nme{iD} '/' 'output_table.csv']);
        end
        fcfg = [];
        for iR = 1:size(tbl_typ)
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = ['copy,' num2str(iN) ',' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
            end
        end
        fcfg.dta = dta_inp;
        fcfg.grp = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        tbl_out  = ejk_create_table( fcfg );
        
        % Add-ons
        num_sbj{1} = 'N';
        tbl_lbl{1} = '';
        for iN = 1:numel(fld_nme)
            num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iN}));
            tbl_lbl{iN+1} = fld_nme{iN};
        end
        
        % Out
        tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
        
        cell2csv([ out_dir '/' 'tables' '/' 'Table2'  '/' 'Cognitive' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '_ttest1.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
        
        clear tbl_lbl num_sbj dta_inp
        
    end
end

%% Correlations
cor_typ = { 'spearman' }; % 'pearson' };

grp_nme     = { 'surgery' };
grp_nme_sub = { { 'pst_cog_dti' } };

cog_var = tbl_typ(:,2);
cog_ind = find(ismember(fnl_dta_col,cog_var));

cor_ind = 1:size(fnl_dta_col,2); cor_ind = setxor(cor_ind,cog_ind); cor_ind = cor_ind(cellfun(@isnumeric,fnl_dta(1,cor_ind))); cor_ind = [ find(strcmpi(fnl_dta_col,'lm2_pre')) cor_ind ];
cor_var = fnl_dta_col(cor_ind);

for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            
            fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
            
            for iFN = 1:numel(fld_nme)
                grp_ind = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iFN});
                
                fcfg = [];
                
                fcfg.sbj_nme = fnl_dta_sbj( grp_ind, 1);
                
                fcfg.dta_two = cell2mat(fnl_dta( grp_ind, cog_ind));
                fcfg.lbl_two = strcat('cog_',fnl_dta_col(cog_ind));
                
                fcfg.cor_typ = cor_typ{iCT};
                
                fcfg.dta_one = cell2mat(fnl_dta( grp_ind, cor_ind));
                fcfg.lbl_one = fnl_dta_col(cor_ind);
                
                fcfg.pvl_cut = 0.05;
                fcfg.pvl_lib = 0.10;
                
                fcfg.force_plot = 0;
                
                fcfg.out_dir = [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' ];
                
                ejk_cross_cor( fcfg );
            end
            
        end
    end
end

% Fisher's Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_lod = { 'ltle_atl' 'ltle_slah' };
for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            
            for iGL = 1:numel(grp_lod)
                grp_rvl{iGL} = mmil_readtext( [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' grp_lod{iGL}  '/' 'cross_correlation_rvalues.csv'] );
                    grp_rvl{iGL}( strcmpi(grp_rvl{iGL},'NA') ) = {0};
                grp_num{iGL} = mmil_readtext( [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' grp_lod{iGL}  '/' 'cross_correlation_n.csv'] );
            end
            
            fcfg = [];
            
            fcfg.col_lbl = grp_rvl{1}(1,2:end);
            fcfg.row_lbl = grp_rvl{1}(2:end,1);
            
            fcfg.rvl_one = grp_rvl{1}(2:end,2:end);
            fcfg.rvl_two = grp_rvl{2}(2:end,2:end);
            
            fcfg.num_one = grp_num{1}(2:end,2:end);
            fcfg.num_two = grp_num{2}(2:end,2:end);
            
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' 'FishersZ'  '/' ];
            fcfg.out_pre = 'ltle_surgery';
            
            ejk_fishersZ(fcfg);
            
        end
    end
end

% Make Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neu_var_nme = { 'Hippocampus' 'Thalamus.Proper' 'fusiform.thick' 'lateralorbitofrontal.thick' 'Unc' 'ILF' 'fusiform' 'lateralorbitofrontal' };

for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            for iCV = 1:numel(cog_var)
                
                fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
                
                ref_tbl = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{1}  '/' 'cross_correlation_rvalues.csv' ]);
                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+2 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                
                for iFN = 1:numel(fld_nme)
                    
                    bil_tbl_hld{1,iFN+1} = fld_nme{iFN};
                    
                    cor_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' 'cross_correlation_rvalues.csv' ]);
                    cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                    pvl_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' 'cross_correlation_pvalues.csv' ]);
                    pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                    num_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' 'cross_correlation_n.csv' ]);
                    num_hld = num_hld(2:end,strcmpi(num_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    num_hld(cellfun(@isstr,num_hld)) = {NaN};
                                                            
                    bil_tbl_hld(2:end,iFN+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                end
                
                fsh_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '/' 'FishersZ'  '/' 'ltle_surgery_fishersZ.csv' ]);
                bil_tbl_hld(:,iFN+2) = fsh_hld(:,strcmpi(fsh_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                
                cell2csv([ out_dir '/' 'tables' '/' 'Table3'  '/' 'Correlations' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '_' cog_var{iCV} '_' cor_typ{iCT} '.csv'],bil_tbl_hld);
                
                % Make Neuro-Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for iR = 1:numel(neu_var_nme); neu_loc_ind(iR) = find(strcmpi(bil_tbl_hld(:,1),neu_var_nme{iR})); end
                neu_bil_tbl_hld = bil_tbl_hld(neu_loc_ind,:);
                
                cell2csv([ out_dir '/' 'tables' '/' 'Table3'  '/' 'Correlations' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '_' cog_var{iCV} '_' cor_typ{iCT} '_neurobio.csv'],bil_tbl_hld);
                
            end
        end
    end
end

% Partial Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_col = { 'lm2_chg' 'lm2_rci' };
neu_col = { 'Thalamus.Proper' 'fusiform.thick' 'lateralorbitofrontal.thick' 'Unc' 'ILF' 'fusiform' 'lateralorbitofrontal' };
cov_col = { 'lm2_pre' 'AgeOfSeizureOnset' 'Hippocampus' };

cog_col = find(ismember(fnl_dta_col,cog_col));
neu_col = find(ismember(fnl_dta_col,neu_col));
cov_col = find(ismember(fnl_dta_col,cov_col));

% Run Correlations
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        for iFN = 1:numel(fld_nme)
            grp_ind = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iFN});
            
            fcfg = [];
            
            fcfg.sbj_nme = fnl_dta_sbj( grp_ind, 1);
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_one = cell2mat(fnl_dta( grp_ind, cog_col));
            fcfg.lbl_one = fnl_dta_col(cog_col);
            
            fcfg.dta_two = cell2mat(fnl_dta( grp_ind, neu_col));
            fcfg.lbl_two = fnl_dta_col(neu_col);
            
            fcfg.dta_cov = cell2mat(fnl_dta( grp_ind, cov_col));
            fcfg.lbl_cov = fnl_dta_col(cov_col);
            
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'partial_correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' fld_nme{iFN}  '/' ];
            
            ejk_partial_correlation( fcfg );
            
        end
        
    end
end







