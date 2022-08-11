out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'neuropsych' '/'];

% Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = [ cog_dta_col ; cog_dta ];

%% Groups
cmp_typ     = { 'QC_none' 'QC_visual_imaging' 'QC_post_imaging' 'QC_imaging_LowCog' };
grp_nme     = { 'grp.mat' 'grp_vis_qal.mat'   'grp_img_qal.mat' 'grp_img_cog_qal.mat' };
grp_run     = { [1 2]     [1 2]               [1 2]             [NaN 2]};

grp_use.nme     = { 'diagnosis'                'surgery' };
grp_use.sub_nme = { {'pre_cog' 'pre_cog_dti'} {'pst_cog' 'pst_cog_dti'} };

cog_col{1} = string_find(cog_dta_col,'_pre');
cog_col{2} = [ string_find(cog_dta_col,'_chg') string_find(cog_dta_col,'_pct')];

%% Statistics
for iCT = 1:numel(cmp_typ)
    
    load([ prj_dir '/' prj_nme '/' 'Data' '/' grp_nme{iCT}])
    
    for iGU = 1:numel(grp_use.nme)
        if any(ismember(iGU,grp_run{iCT}))
            for iGS = 1:numel(grp_use.sub_nme{iGU})
                
                % get group data
                fld_nme = fieldnames( grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}) );
                
                fcfg = [];
                fcfg.grp     = grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS});
                fcfg.grp_inc = {fld_nme};
                fcfg.grp_nme = {fld_nme};
                fcfg.dta = cog_dta(:,cog_col{grp_run{iCT}(iGU)});
                fcfg.sbj = cog_dta_sbj;
                [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
                
                % ANOVA
                fcfg = [];
                fcfg.sbj_nme = grp_sbj{1};
                fcfg.dta     = grp_dta{1};
                fcfg.dta_nme = cog_dta_col(:,cog_col{grp_run{iCT}(iGU)});
                fcfg.grp     = grp_typ{1};
                fcfg.grp_nme = grp_use.nme(iGU);
                fcfg.out_dir = [ out_dir '/' 'ANOVA' '/' grp_use.nme{iGU} '/' grp_use.sub_nme{iGU}{iGS} '_' cmp_typ{iCT} '/'];
                ejk_1way_anova( fcfg )
                
            end
        end
    end
end

%% Table
for iCT = 1:numel(cmp_typ)
    
    % Load Groups
    load([ prj_dir '/' prj_nme '/' 'Data' '/' grp_nme{iCT}])
    
    for iGU = 1:numel(grp_use.nme)
        if any(ismember(iGU,grp_run{iCT}))
        for iGS = 1:numel(grp_use.sub_nme{iGU})
                
            fld_nme = fieldnames( grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}) );
    
            dta_inp{2} = mmil_readtext([ out_dir '/' 'ANOVA' '/' grp_use.nme{iGU} '/' grp_use.sub_nme{iGU}{iGS} '_' cmp_typ{iCT} '/' grp_use.nme{iGU} '/' 'output_table.csv']);
            
            % Create Table
            fcfg = [];
            for iR = 1:numel(cog_col{grp_run{iCT}(iGU)})
                for iN = 1:numel(fld_nme)
                    fcfg.tbl{iR,iN} = [ 'mean/std' ',' '1' ',' fld_nme{iN} ',' cog_dta_col{cog_col{grp_run{iCT}(iGU)}(iR)}];
                end
                fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(cog_dta_col{cog_col{grp_run{iCT}(iGU)}(iR)},{'_'},{'.'}) ',report'];
            end
            
            fcfg.dta = dta_inp;
            fcfg.grp = grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS});
            tbl_out  = ejk_create_table( fcfg );
            
            % Add-ons
            num_sbj{1} = 'N';
            tbl_lbl{1} = '';
            for iN = 1:numel(fld_nme)
                num_sbj{iN+1} = numel(grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}).(fld_nme{iN}));
                tbl_lbl{iN+1} = fld_nme{iN};
            end
            num_sbj{iN+2} = '-';
            tbl_lbl{iN+2} = 'Stats';
            
            tbl_out = [ tbl_lbl; num_sbj; cog_dta_col(cog_col{grp_run{iCT}(iGU)})' tbl_out ];
    
            cell2csv([ out_dir '/' 'neuropsych_table' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '_' cmp_typ{iCT} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
        
        end
        end
    end    
    
end

