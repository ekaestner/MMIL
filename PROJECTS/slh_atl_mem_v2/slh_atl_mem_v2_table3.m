out_dir = [ prj_dir '/' 'out' '/' 'Correlations' '/'];
stt_dir = [ out_dir '/' 'statistics' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(stt_dir);

%%
load([ dta_dir '/' 'group_quality_check.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%%
grp_use.nme     = { 'srg_sde_ana'                                          'sde_ana' };
grp_use.sub_nme = { {'initial_QC' 'eng_one_QC' 'eng_two_QC' 'lng_lft_QC'} {'initial_QC' 'eng_one_QC' 'eng_two_QC' 'lng_lft_QC'} };
grp_use.fsh_zsc = { { {'lft_atl' 'lft_slh'} {'lft_atl' 'rgh_atl'} {'lft_slh' 'rgh_slh'} } ...
                    { {'lft' 'rgh'} } };

cor_typ = { 'spearman' };
cog_var = { 'lm2_pre' 'lm2_chg' 'lm2_pst'};

%% Bivariate Correlations
cog_ind = find(ismember(tot_col,cog_var));

cor_ind = 1:size(tot_col,2); cor_ind = setxor(cor_ind,cog_ind); 
cor_ind = cor_ind(all(cellfun(@isnumeric,tot_dta(:,cor_ind)))); 
cor_ind = cor_ind(~all(cellfun(@isempty,tot_dta(:,cor_ind)))); 
cor_ind = [ find(strcmpi(tot_col,'lm2_pre')) cor_ind ];

cor_var = tot_col(cor_ind);

for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iG})
            
            fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
            
            for iFN = 1:numel(fld_nme)
                grp_ind = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iFN});
                
                fcfg = [];
                
                fcfg.sbj_nme = tot_sbj( grp_ind, 1);
                
                fcfg.dta_two = cell2mat(tot_dta( grp_ind, cog_ind));
                fcfg.lbl_two = strcat('cog_',tot_col(cog_ind));
                
                fcfg.cor_typ = cor_typ{iCT};
                
                fcfg.dta_one = cell2mat(tot_dta( grp_ind, cor_ind));
                fcfg.lbl_one = tot_col(cor_ind);
                
                fcfg.pvl_cut = 0.05;
                fcfg.pvl_lib = 0.10;
                
                fcfg.force_plot = 0;
                
                fcfg.out_dir = [ stt_dir '/' 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' ];
                
                ejk_cross_cor( fcfg );
            end
            
        end
    end
end

%% Robust correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_ind = find(ismember(tot_col,cog_var));

cor_ind = 1:size(tot_col,2); cor_ind = setxor(cor_ind,cog_ind);
cor_ind = cor_ind(all(cellfun(@isnumeric,tot_dta(:,cor_ind))));
cor_ind = cor_ind(~all(cellfun(@isempty,tot_dta(:,cor_ind))));
cor_ind = [ find(strcmpi(tot_col,'lm2_pre')) cor_ind ];

cor_var = tot_col(cor_ind);

for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iG})
            
            fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
            
            for iFN = 1:numel(fld_nme)
                
                grp_ind = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iFN});
                
                fcfg = [];
                
                fcfg.sbj_nme = tot_sbj( grp_ind, 1);
                
                fcfg.dta_two = cell2mat(tot_dta( grp_ind, cog_ind));
                fcfg.lbl_two = strcat('cog_',tot_col(cog_ind));
                
                fcfg.dta_one = cell2mat(tot_dta( grp_ind, cor_ind));
                fcfg.lbl_one = tot_col(cor_ind);
                
                fcfg.pvl_cut = 0.05;
                fcfg.pvl_lib = 0.10;
                
                fcfg.force_plot = 0;
                
                fcfg.out_dir = [ stt_dir '/' 'robustcorrelation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' ];
                
                ejk_cross_cor_robust( fcfg );
            end
            
        end
    end
end

%% Fisher's Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ejk_chk_dir([  stt_dir '/' 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' 'FishersZ'  '/' ]);

rob_typ = { '' 'robust' };

for iRB = 1%:numel(rob_typ)
    for iCT = 1:numel(cor_typ)
        for iG = 1:numel(grp_use.nme)
            for iGS = 1:numel(grp_use.sub_nme{iG})
                
                for iCM = 1:numel( grp_use.fsh_zsc{iG} )
                    
                    for iGL = 1:numel(grp_use.fsh_zsc{iG}{iCM})
                        grp_rvl{iGL} = mmil_readtext( [ stt_dir '/' rob_typ{iRB} 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' grp_use.fsh_zsc{iG}{iCM}{iGL}  '/' 'cross_correlation_rvalues.csv'] );
                        grp_rvl{iGL}( strcmpi(grp_rvl{iGL},'NA') ) = {0};
                        grp_num{iGL} = mmil_readtext( [ stt_dir '/' rob_typ{iRB} 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' grp_use.fsh_zsc{iG}{iCM}{iGL}  '/' 'cross_correlation_n.csv'] );
                    end
                    
                    
                    fcfg = [];
                    
                    fcfg.col_lbl = grp_rvl{1}(1,2:end);
                    fcfg.row_lbl = grp_rvl{1}(2:end,1);
                    
                    fcfg.rvl_one = grp_rvl{1}(2:end,2:end);
                    fcfg.rvl_two = grp_rvl{2}(2:end,2:end);
                    
                    fcfg.num_one = grp_num{1}(2:end,2:end);
                    fcfg.num_two = grp_num{2}(2:end,2:end);
                    
                    fcfg.out_dir = [  stt_dir '/' rob_typ{iRB} 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' 'FishersZ'  '/' ];
                    fcfg.out_pre = [ grp_use.fsh_zsc{iG}{iCM}{1} '_VS_' grp_use.fsh_zsc{iG}{iCM}{2}];
                    
                    ejk_fishersZ(fcfg);
                    
                    clear grp_rvl grp_num
                    
                end
            end
        end
    end
end

%% Partial Correlations
% Partial Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cog_col = { 'lm2_chg' };
neu_col = { 'Thalamus_subcort_vol_Laterality' 'fusiform_cort_thick_ctx_Laterality' 'lateralorbitofrontal_cort_thick_ctx_Laterality' 'Unc_fiber_FA_Laterality' 'ILF_fiber_FA_Laterality' 'IFO_fiber_FA_Laterality' 'tSLF_fiber_FA_Laterality' 'fusiform_gwcsurf_FA_wm_ctx_Laterality' 'lateralorbitofrontal_gwcsurf_FA_wm_ctx_Laterality' };
cov_col = { 'lm2_pre' 'age_sze_ons' 'Hippocampus_subcort_vol_Laterality' };

cog_col = find(ismember(tot_col,cog_col));
neu_col = find(ismember(tot_col,neu_col));
cov_col = find(ismember(tot_col,cov_col));

% Run Correlations
for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_use.nme)
        for iGS = 1:numel(grp_use.sub_nme{iG})
            
            fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
            
            for iFN = 1:numel(fld_nme)
                grp_ind = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iFN});
                
                fcfg = [];
                
                fcfg.sbj_nme = tot_sbj( grp_ind, 1);
                
                fcfg.cor_typ = 'spearman';
                
                fcfg.dta_one = cell2mat(tot_dta( grp_ind, cog_col));
                fcfg.lbl_one = tot_col(cog_col);
                
                fcfg.dta_two = cell2mat(tot_dta( grp_ind, neu_col));
                fcfg.lbl_two = tot_col(neu_col);
                
                fcfg.dta_cov = cell2mat(tot_dta( grp_ind, cov_col));
                fcfg.lbl_cov = tot_col(cov_col);
                
                fcfg.out_dir = [ stt_dir '/' 'partial_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} '/' fld_nme{iFN}  '/' ];
                
                ejk_partial_correlation( fcfg );
                
            end
            
        end
    end
end
    
%% Bivariate Table
rob_typ = { '' 'robust' };

for iRB = 1%:numel(rob_typ)
    for iCT = 1:numel(cor_typ)
        for iG = 1:numel(grp_use.nme)
            for iGS = 1:numel(grp_use.sub_nme{iG})
                
                col_add = 0;
                fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
                dta_tbl = cell(numel(cor_var), numel(cog_var)*numel(fld_nme)+numel(cog_var)-1 );
                col_tbl = cell(1,numel(cog_var)*numel(fld_nme)+numel(cog_var)-1);
                
                for iC = 1:numel(cog_var)
                    for iFN = 1:numel(fld_nme)
                        
                        if strcmpi(rob_typ{iRB},'');           fdr_nme = [ rob_typ{iRB} 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '_' cor_typ{iCT} ]; 
                        elseif strcmpi(rob_typ{iRB},'robust'); fdr_nme = [ rob_typ{iRB} 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} ]; end
                        
                        rvl_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_rvalues.csv' ] );
                        rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
                        pvl_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_pvalues.csv' ] );
                        pvl_hld(strcmpi(pvl_hld,'NA')) = {NaN};
                        num_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_n.csv' ] );
                        
                        col_tbl{iFN+col_add} = [ cog_var{iC} '_' fld_nme{iFN}];
                        for iV = 1:numel(cor_var)
                            
                            col_ind = strcmpi(rvl_hld(1,:),[ 'cog' '.' mmil_spec_char(cog_var{iC},{'_'},{'.'}) ]);
                            row_ind = strcmpi(rvl_hld(:,1),mmil_spec_char(cor_var{iV},{'_'},{'.'}) );
                            
                            dta_tbl{iV,iFN+col_add} = [ 'p = ' num2str(roundsd(pvl_hld{row_ind,col_ind},2)) '; r(' num2str(num_hld{row_ind,col_ind}-2) ') = ' num2str(roundsd(rvl_hld{row_ind,col_ind},2))];
                            
                        end
                    end
                    col_add = col_add + numel(fld_nme) + 1;
                    
                end
                
                cell2csv( [ stt_dir '/' fdr_nme '/' 'correlation_table.csv'], [ {''} col_tbl ; cor_var' dta_tbl ]);
                
            end
        end
    end
end

%% Organize
fld_tbl_nme = { 'sde_ana' 'srg_sde_ana'};
smp_tbl_nme = {'initial_QC' 'eng_one_QC' 'eng_two_QC' 'lng_lft_QC'};
cog_tbl_nme = { 'lm2_pre' 'lm2_chg' };
grp_tbl_nme = { {'lft' 'rgh'} {'lft_slh' 'lft_atl' 'rgh_slh' 'lft_slh'} };

cor_tbl_nme = { 'ILF_fiber_FA_Laterality' ; ...
                'Unc_fiber_FA_Laterality' ; ...
                '' ; ...
                'IFO_fiber_FA_Laterality' ; ...
                'tSLF_fiber_FA_Laterality' ; ...
                '' ; ...
                };

out_dir
for iFT = 1:numel(fld_tbl_nme)
    
    
    
    
    
end






