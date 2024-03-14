out_dir = [ prj_dir '/' 'out' '/' 'Clinical' '/'];
stt_dir = [ out_dir '/' 'statistics' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(stt_dir);

%% Clinical
load([ dta_dir '/' 'group_initial.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%%
cmp_typ     = { 'QC_nuisance' };
grp_nme     = { 'group_quality_check.mat' };

grp_use.nme     = { 'srg_sde_ana'                                          'sde_ana' };
grp_use.sub_nme = { {'initial_QC' 'eng_one_QC' 'eng_two_QC' 'lng_lft_QC' } {'initial_QC' 'eng_one_QC' 'eng_two_QC' 'lng_lft_QC' }  };

%%
% unique(tot_dta(:,strcmpi(tot_col,'lng_lat_mth')))
tbl_typ = { ...'count'    'Location'                'HC/UCSD/UCSF/Emory' 1; ... %%Replace%%
            'mean/std' 'age'                     ''            1 ; ... 
            'mean/std' 'edu'                     ''            1 ; ... 
            'count'    'sex'                     'M/F'         1 ; ... 
            'count'    'hnd'                     'R/L'         1 ; ...
            'mean/std' 'age_sze_ons'             ''            1 ; ... 
            'mean/std' 'num_aed'                 ''            1 ; ... 
            'count'    'sde_sze_ons'             'R/L'         1 ; ...
            'mean/std' 'sze_frq'                 ''            1 ; ... 
            'count'    'mts'                     'yes/no'      1 ; ...
            'count'    'srg_typ'                 'ATL/SLAH'    1 ; ... 
            'count'    'eng_out'                 'I/II/III/IV' 1 ; ...
            'count'    'eng_out_one'             'I/II+' 1 ; ...
            'count'    'eng_out_two'             'I_II/III_IV' 1 ; ...
            'count'    'lng_lat'                 'L/Atypical'       1 ; ...
            }; ...'count'    'lng_lat_mth'             ''            1 };    %%Needs Work%%

dta_inp{1} = [ tot_col ; tot_dta ];

%% Statistics
anv_var = tbl_typ(strcmpi(tbl_typ(:,1),'mean/std'),2);
fsh_use = tbl_typ(strcmpi(tbl_typ(:,1),'count'),2);

for iCT = 1:numel(cmp_typ)
    
    % Load Groups
    load([ dta_dir '/' grp_nme{iCT}])
    
    for iGU = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iGU})
            
            fld_nme = fieldnames( grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}) );
            
            % ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, use_dta_col ] = intersect( tot_col, anv_var );
            
            fcfg = [];
            fcfg.grp     = grp.(grp_use.nme{iGU}).(grp_use.sub_nme {iGU}{iGS});
            fcfg.grp_inc = {fld_nme};
            fcfg.grp_nme = {fld_nme};
            fcfg.dta = tot_dta(:,use_dta_col);
            fcfg.sbj = tot_sbj;
            [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
                                  
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = tot_col(use_dta_col);
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = {[ grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ]};
            fcfg.out_dir = [ stt_dir '/' 'anova' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ];
            ejk_1way_anova( fcfg )
            
            % fisher test %%%%%%%%%%%%%%%%%%%%%            
            [~, use_dta_col ] = intersect( tot_col, fsh_use );
            
            fcfg = [];
            fcfg.grp     = grp.(grp_use.nme{iGU}).(grp_use.sub_nme {iGU}{iGS});
            fcfg.grp_inc = {fld_nme};
            fcfg.grp_nme = {fld_nme};
            fcfg.dta = tot_dta(:,use_dta_col);
            fcfg.sbj = tot_sbj;
            [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
            
            grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {''};
            
            [~, lvl_ind] = ismember(tbl_typ(:,2),tot_col(use_dta_col));
            clear cln_lvl; for iLI = 1:numel(lvl_ind); if lvl_ind(iLI); cln_lvl(lvl_ind(iLI)) = tbl_typ(iLI,3); end; end
            
            fcfg = [];
            
            fcfg.sbj = grp_sbj{1};
            
            fcfg.dta_one = grp_dta{1};
            fcfg.lbl_one = tot_col(use_dta_col);
            fcfg.lvl     = cln_lvl;
            
            fcfg.dta_two = repmat(grp_typ{1},1,numel(use_dta_col));
            fcfg.lbl_two = strcat( 'group_', tot_col(use_dta_col));
            
            fcfg.out_dir = [ stt_dir '/' ];
            fcfg.grp_nme = [ 'fisher' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ];
            
            ejk_fisher_test( fcfg );
            
        end
    end
end

%% Make Table
for iCT = 1:numel(cmp_typ)
    
    % Load Groups
    load([ dta_dir '/' grp_nme{iCT}])

    for iGU = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iGU})
    
            fld_nme = fieldnames( grp.(grp_use.nme{iGU}).(grp_use.sub_nme{iGU}{iGS}) );
    
            % Load stats
            dta_inp{2} = mmil_readtext([ stt_dir '/' 'anova' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '/' mmil_spec_char(grp_use.nme{iGU},{'_'},{'.'}) '.' mmil_spec_char(grp_use.sub_nme{iGU}{iGS},{'_'},{'.'}) '/' 'output_table.csv' ]);
            dta_inp{3} = mmil_readtext([stt_dir '/' 'fisher' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '/' 'output_table.csv' ]);
            
            % Create Table
            fcfg = [];
            
            for iR = 1:size(tbl_typ)
                if strcmpi(tbl_typ{iR,1},'mean/std')
                    for iN = 1:numel(fld_nme)
                        fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
                    end
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                elseif  strcmpi(tbl_typ{iR,1},'count')
                    for iN = 1:numel(fld_nme)
                        fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
                    end
                    fcfg.tbl{iR,iN+1} = ['copy,3,' tbl_typ{iR,2} ',report'];
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

%% Document Missing









