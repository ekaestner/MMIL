out_dir = [ prj_dir '/' 'out' '/' 'Correlations' '/'];
stt_dir = [ out_dir '/' 'statistics' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(stt_dir);

%%
cog_val         = 'lm2.chg';
cor_typ         = 'spearman';
grp_use_nme     = { 'srg_sde_ana' 'srg_sde_ana' };
grp_use_sub_nme = { 'initial_QC'  'initial_QC' };
grp_fld_nme     = { 'lft_slh'     'lft_atl' };
neu_bio_nme     = { 'Unc_fiber_FA_Laterality' 'ILF_fiber_FA_Laterality' 'IFO_fiber_FA_Laterality' 'tSLF_fiber_FA_Laterality' };
num_cov         = 3;

%%
tbl_itr = num_cov + 1 + 1;
tbl_out = cell( (tbl_itr*numel(neu_bio_nme)), numel(grp_use_nme)+1);

for iG = 1:numel(grp_use_nme)
    
    tbl_ind = 1;
    for iN = 1:numel(neu_bio_nme)
        
        pvl_hld = mmil_readtext( [ stt_dir '/' 'partial_' '_' grp_use_nme{iG} '_' grp_use_sub_nme{iG} '_' cor_typ '/' grp_fld_nme{iG} '/' cog_val '/' 'partial' '/' 'partial_correlation' '_' mmil_spec_char(neu_bio_nme{iN},{'_'},{'.'}) '_' 'pvalues.csv' ] );
        rvl_hld = mmil_readtext( [ stt_dir '/' 'partial_' '_' grp_use_nme{iG} '_' grp_use_sub_nme{iG} '_' cor_typ '/' grp_fld_nme{iG} '/' cog_val '/' 'partial' '/' 'partial_correlation' '_' mmil_spec_char(neu_bio_nme{iN},{'_'},{'.'}) '_' 'rvalues.csv' ] );
        
        for iC = 1:tbl_itr-1
            tbl_out{tbl_ind+iC-1,1}   = pvl_hld{iC+2,1};
            tbl_out{tbl_ind+iC-1,iG+1} = [ 'rs=' ' ' num2str(roundsd(rvl_hld{iC+2,2},2)) ' ; ' 'p=' ' ' num2str(roundsd(pvl_hld{iC+2,2},2)) ];
        end
        tbl_ind = tbl_ind + tbl_itr;
    end
    
end
 
