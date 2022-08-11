out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Clinical' '/'];

% Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

%% Groups
cmp_typ     = { 'QC_none' 'QC_visual_imaging' 'QC_post_imaging' 'QC_imaging_LowCog' };
grp_nme     = { 'grp.mat' 'grp_vis_qal.mat'   'grp_img_qal.mat' 'grp_img_cog_qal.mat' };
grp_run     = { [1 2]     [1 2]               [1 2]             [2]};

grp_use.nme     = { 'diagnosis'                'surgery' };
grp_use.sub_nme = { {'pre_cog' 'pre_cog_dti'} {'pst_cog' 'pst_cog_dti'} };

%% Specify table
tbl_typ = { 'count'    'Location'                'HC/UCSD/UCSF/Emory' 1; ...
            'mean/std' 'AgeAtSurgery'            ''         1 ; ...
            'mean/std' 'Educ'                    ''         1 ; ...
            'count'    'Sex'                     'M/F'      1 ; ...
            'count'    'Handedness'              'R/L'      1 ; ...
            'mean/std' 'AgeOfSeizureOnset'       ''         1 ; ...
            'mean/std' 'NumAEDs'                 ''         1 ; ...
            'count'    'SideOfSeizureFocus'      'R/L'      1 ; ...
            'mean/std' 'SeizureFreq'             ''         1 ; ...
            'count'    'SurgeryType'             'ATL/SLAH' 1 };

dta_inp{1} = [ cln_dta_col ; cln_dta ];

%% Statistics

%%
for iCT = 1:numel(cmp_typ)
    
    % Load Groups
    load([ prj_dir '/' prj_nme '/' 'Data' '/' grp_nme{iCT}])

    for iGU = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iGU})
    
            fld_nme = fieldnames( grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}) );
    
            % Create Table
            fcfg = [];
            
            for iR = 1:size(tbl_typ)
                if strcmpi(tbl_typ{iR,1},'mean/std')
                    for iN = 1:numel(fld_nme)
                        fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
                    end
                    fcfg.tbl{iR,iN+1} = 'empty';%['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                elseif  strcmpi(tbl_typ{iR,1},'count')
                    for iN = 1:numel(fld_nme)
                        fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
                    end
                    fcfg.tbl{iR,iN+1} = 'empty';% ['copy,3,' tbl_typ{iR,2} ',report'];
                end
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
            
            % Out
            tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
    
            cell2csv([ out_dir '/' 'clinical_table' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '_' cmp_typ{iCT} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
            
        end
    end
end

%%
% open emory_explore