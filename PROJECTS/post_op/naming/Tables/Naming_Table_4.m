out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Tables/Parts';

cor_dir = 'SpecificCor';

hld_dir.cln_dir     = [ prj_dir '/' prj_nme '/' cor_dir '/' 'Clinical' '/' ];
    hld_dir.cln_nme_dir = '';
    hld_dir.cln_grp_dir = 'TLE_Controls_pre_pre';
hld_dir.fib_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/' 'DTI' '/' 'fiber_FA' '/' ];
    hld_dir.fib_nme_dir = 'Raw';
    hld_dir.fib_grp_dir = 'tle_controls_pre_3T_allSurg_all';
hld_dir.wmp_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/'  'DTI' '/' 'wmparc_FA_wm' '/' ];
    hld_dir.wmp_nme_dir = 'Raw';
    hld_dir.wmp_grp_dir = 'tle_controls_pre_3T_allSurg_all';
hld_dir.icv_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/'  'MRI' '/' 'subcort_vol_ICV_cor' '/' ];
    hld_dir.icv_nme_dir = 'Raw';
    hld_dir.icv_grp_dir = 'tle_controls_pre_3T_allSurg_all';
    
bnt_col = 2;
ant_col = 3;

%%
ovr_nme = { 'cln_dir'     'icv_dta_dir' 'fib_dta_dir' 'wmp_dta_dir' };
nme_dir = { 'cln_nme_dir' 'icv_nme_dir' 'fib_nme_dir' 'wmp_nme_dir' };
grp_dir = { 'cln_grp_dir' 'icv_grp_dir' 'fib_grp_dir' 'wmp_grp_dir' };

cor_nme = { { 'AgeAtSurgery' 'Educ' } ...
            { 'xLeft.Hippocampus' 'xRight.Hippocampus' } ...
            { 'xL.IFO' 'xR.IFO' 'xL.ILF' 'xR.ILF' } ...
            { 'xlh.fusiform' 'xrh.fusiform' } };

row_ind = 1;
for iO = 1:numel(ovr_nme)
    
    rvl_hld = mmil_readtext( [ hld_dir.(ovr_nme{iO}) '/' hld_dir.(nme_dir{iO}) '/' hld_dir.(grp_dir{iO}) '/' 'cross_correlation_rvalues.csv' ] );
    pvl_hld = mmil_readtext( [ hld_dir.(ovr_nme{iO}) '/' hld_dir.(nme_dir{iO}) '/' hld_dir.(grp_dir{iO}) '/' 'cross_correlation_pvalues.csv' ] );
    
    for iC = 1:numel(cor_nme{iO})
    
        cor_ind = strcmpi(rvl_hld(:,1), cor_nme{iO}{iC});
        
        rvl_tbl_out{row_ind,1} = cor_nme{iO}{iC};
        
        rvl_tbl_out{row_ind,2} = num2str(roundsd(rvl_hld{cor_ind,bnt_col},2));
        if rvl_hld{cor_ind,bnt_col}<0; rvl_tbl_out{row_ind,2}=['-' rvl_tbl_out{row_ind,2}(3:end)]; elseif rvl_hld{cor_ind,bnt_col}>=0; rvl_tbl_out{row_ind,2}=rvl_tbl_out{row_ind,2}(2:end); end
        if pvl_hld{cor_ind,bnt_col}<.01
            rvl_tbl_out{row_ind,2} = [ rvl_tbl_out{row_ind,2} '**' ];
        elseif pvl_hld{cor_ind,bnt_col}<.05
            rvl_tbl_out{row_ind,2} = [ rvl_tbl_out{row_ind,2} '*' ];
        end
        
        rvl_tbl_out{row_ind,3} = num2str(roundsd(rvl_hld{cor_ind,ant_col},2));
        if rvl_hld{cor_ind,ant_col}<0; rvl_tbl_out{row_ind,3}=['-' rvl_tbl_out{row_ind,3}(3:end)]; elseif rvl_hld{cor_ind,ant_col}>=0; rvl_tbl_out{row_ind,3}=rvl_tbl_out{row_ind,3}(3:end); end
        if pvl_hld{cor_ind,ant_col}<.01
            rvl_tbl_out{row_ind,3} = [ rvl_tbl_out{row_ind,3} '**' ];
        elseif pvl_hld{cor_ind,ant_col}<.05
            rvl_tbl_out{row_ind,3} = [ rvl_tbl_out{row_ind,3} '*' ];
        end
        
        row_ind = row_ind + 1;
        
    end
end

cell2csv([ out_dir '/' 'table3_preoperative.csv'], rvl_tbl_out, ';')
clear rvl_tbl_out








