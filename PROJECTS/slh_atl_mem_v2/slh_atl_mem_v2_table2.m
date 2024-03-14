out_dir = [ prj_dir '/' 'out' '/' 'Cognitive' '/'];
stt_dir = [ out_dir '/' 'statistics' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(stt_dir);

%%
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
tbl_typ = { 'mean/std' 'lm2_pre' '' 1 ; ... 
            'mean/std' 'lm2_pst' '' 1 ; ... 
            'mean/std' 'lm2_chg' '' 1 };

dta_inp{1} = [ tot_col ; tot_dta ];

%% Statistics
anv_var = tbl_typ(strcmpi(tbl_typ(:,1),'mean/std'),2);

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
              
            % One way ANOVA
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = tot_col(use_dta_col);
            fcfg.grp     = grp_typ{1};
            fcfg.grp_nme = {[ grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ]};
            fcfg.out_dir = [ stt_dir '/' 'anova' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ];
            ejk_1way_anova( fcfg )
            
            % Two way ANOVA
            grp_typ_hld = regexp(grp_typ{1},'_','split');
            grp_typ{1}  = cat(1,grp_typ_hld{:});
            
            if size(grp_typ{1},2)==2
                fcfg = [];
                fcfg.sbj_nme = grp_sbj{1};
                fcfg.dta     = grp_dta{1};
                fcfg.dta_nme = tot_col(use_dta_col);
                fcfg.grp_one     = grp_typ{1}(:,1);
                fcfg.grp_nme_one = {'side'};
                fcfg.grp_two     = grp_typ{1}(:,2);
                fcfg.grp_nme_two = {[ grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} ]};
                fcfg.out_dir = [ stt_dir '/' 'anova2way' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS}  ];
                ejk_2way_anova( fcfg )
            end
            
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
            dta_inp{2} = mmil_readtext([ stt_dir '/' 'anova' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS}     '/' mmil_spec_char(grp_use.nme{iGU},{'_'},{'.'}) '.' mmil_spec_char(grp_use.sub_nme{iGU}{iGS},{'_'},{'.'}) '/' 'output_table.csv' ]);
            if strcmpi(grp_use.nme{iGU},'srg_sde_ana')
                dta_inp{3} = mmil_readtext([ stt_dir '/' 'anova2way' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '/' 'side' '/' 'output_table.csv']);
            end
            
            % Create Table
            fcfg = [];
            
            for iR = 1:size(tbl_typ)
                if strcmpi(tbl_typ{iR,1},'mean/std')
                    for iN = 1:numel(fld_nme)
                        fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
                    end
                    fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
                    if strcmpi(grp_use.nme{iGU},'srg_sde_ana')
                        fcfg.tbl{iR,iN+2} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_one'];
                        fcfg.tbl{iR,iN+3} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_two'];
                        fcfg.tbl{iR,iN+4} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_int'];
                    end
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
            if strcmpi(grp_use.nme{iGU},'srg_sde_ana')
                num_sbj{iN+3} = '-';num_sbj{iN+4} = '-';num_sbj{iN+5} = '-';
                tbl_lbl{iN+3} = 'Side';
                tbl_lbl{iN+4} = 'Surgery';
                tbl_lbl{iN+5} = 'Side*Surgery';
            end
            
            % Out
            tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
            
            cell2csv([ out_dir '/' 'cognitive_table' '_' grp_use.nme{iGU} '_' grp_use.sub_nme{iGU}{iGS} '_' cmp_typ{iCT} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
            
        end
    end
end

