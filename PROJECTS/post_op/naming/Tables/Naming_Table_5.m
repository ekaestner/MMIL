out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Tables/Parts';

cor_dir = 'SpecificCor';

hld_dir.cog_dir     = [ prj_dir '/' prj_nme '/' cor_dir '/' 'Cognitive' '/' ];
    hld_dir.cog_nme_dir = '';
    hld_dir.cog_grp_dir_lft = 'LTLE_pre_post';
    hld_dir.cog_grp_dir_rgh = 'RTLE_pre_post';
hld_dir.cln_dir     = [ prj_dir '/' prj_nme '/' cor_dir '/' 'Clinical' '/' ];
    hld_dir.cln_nme_dir = '';
    hld_dir.cln_grp_dir_lft = 'LTLE_post_post';
    hld_dir.cln_grp_dir_rgh = 'RTLE_post_post';
hld_dir.fib_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/' 'DTI' '/' 'fiber_FA' '/' ];
    hld_dir.fib_nme_dir = 'Raw';
    hld_dir.fib_grp_dir_lft = 'tle_post_3T_ATLonly_left';
    hld_dir.fib_grp_dir_rgh = 'tle_post_3T_ATLonly_right';
hld_dir.wmp_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/'  'DTI' '/' 'wmparc_FA_wm' '/' ];
    hld_dir.wmp_nme_dir = 'Raw';
    hld_dir.wmp_grp_dir_lft = 'tle_post_3T_ATLonly_left';
    hld_dir.wmp_grp_dir_rgh = 'tle_post_3T_ATLonly_right';
hld_dir.icv_dta_dir = [ prj_dir '/' prj_nme '/' cor_dir '/'  'MRI' '/' 'subcort_vol_ICV_cor' '/' ];
    hld_dir.icv_nme_dir = 'Raw';
    hld_dir.icv_grp_dir_lft = 'tle_post_3T_ATLonly_left';
    hld_dir.icv_grp_dir_rgh = 'tle_post_3T_ATLonly_right';
    
bnt_col = 2;
ant_col = 3;

%% LTLE
ovr_nme = { 'cog_dir'         'cln_dir'         'icv_dta_dir'     'fib_dta_dir'     'wmp_dta_dir' };
nme_dir = { 'cog_nme_dir'     'cln_nme_dir'     'icv_nme_dir'     'fib_nme_dir'     'wmp_nme_dir' };
grp_dir = { 'cog_grp_dir_lft' 'cln_grp_dir_lft' 'icv_grp_dir_lft' 'fib_grp_dir_lft' 'wmp_grp_dir_lft' };

cor_nme = { { 'bnt.raw.scr' 'ant.mem.raw.scr' }...
            { 'AgeAtSurgery' 'Educ' 'AgeOfSeizureOnset' } ...
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

cell2csv([ out_dir '/' 'table4_LTLE.csv'], rvl_tbl_out, ';')
clear rvl_tbl_out

%% RTLE
ovr_nme = { 'cog_dir'         'cln_dir'         'icv_dta_dir'     'fib_dta_dir'     'wmp_dta_dir' };
nme_dir = { 'cog_nme_dir'     'cln_nme_dir'     'icv_nme_dir'     'fib_nme_dir'     'wmp_nme_dir' };
grp_dir = { 'cog_grp_dir_rgh' 'cln_grp_dir_rgh' 'icv_grp_dir_rgh' 'fib_grp_dir_rgh' 'wmp_grp_dir_rgh' };

cor_nme = { { 'bnt.raw.scr' 'ant.mem.raw.scr' }...
            { 'AgeAtSurgery' 'Educ' 'AgeOfSeizureOnset' } ...
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

cell2csv([ out_dir '/' 'table4_RTLE.csv'], rvl_tbl_out, ';')
clear rvl_tbl_out