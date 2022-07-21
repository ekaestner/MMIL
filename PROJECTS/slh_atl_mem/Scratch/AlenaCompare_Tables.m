out_dir =  '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/InitialAnalysis/Tables';

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

cor_cln_cog_loc = '/InitialAnalysis/Clinical/Correlation/';
cor_cog_cog_loc = '/InitialAnalysis/Cognitive/Correlation/';
cor_neu_cog_loc = '/InitialAnalysis/Neurobio_compare/Correlation/Cognitive/';

%% Table - Post LTLE - Laterality  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl

% Data Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_loc = { [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'LTLE_post_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'LTLE_post_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'LTLE_post_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'LTLE_post_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'subcort_vol_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'subcort_vol_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'subcort_vol_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'subcort_vol_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'fiber_FA_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'fiber_FA_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'fiber_FA_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'fiber_FA_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex'  ] };

% Table Construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_nme = { 'Educ' 3; ...
            'Educ' 4; ...
            'AgeOfSeizureOnset' 3; ...
            'AgeOfSeizureOnset' 4; ...
            'AgeAtSurgery' 3; ...
            'AgeAtSurgery' 4; ...
            'hippocampus' 5 ; ...
            'hippocampus' 6 ; ...
            'hippocampus' 7 ; ...
            'hippocampus' 8 ; ...
            'ILF' 9 ; ...
            'ILF' 10 ; ...
            'ILF' 11 ; ...
            'ILF' 12 ; ...
            'Unc' 9 ; ...
            'Unc' 10 ; ...
            'Unc' 11 ; ...
            'Unc' 12 ; ...
            'parahippocampal' 13 ; ...
            'parahippocampal' 14 ; ...
            'parahippocampal' 15 ; ...
            'parahippocampal' 16 ; ... 
            'entorhinal' 13 ; ...
            'entorhinal' 14 ; ...
            'entorhinal' 15 ; ...
            'entorhinal' 16 };
        
col_lbl = { '' 'Post LM2' 'Post VPA2' };
row_lbl = { 'Education - allT' 2 ;...
            'Education - 3T' 2 ;...
            'AgeOfSeizureOnset - allT' 2 ;...
            'AgeOfSeizureOnset - 3T' 2 ;...
            'AgeAtSurgery - allT' 2 ;...
            'AgeAtSurgery - 3T' 2 ;...
            'Hippocampus - 238; allT' 2 ;...
            'Hippocampus - 238; 3T' 2 ;...
            'Hippocampus - DEV; allT' 2 ;...
            'Hippocampus - DEV; 3T' 2 ;...
            'ILF - 238; allT' 2 ;...
            'ILF - 238; 3T' 2 ;...
            'ILF - DEV; allT' 2 ;...
            'ILF - DEV; 3T' 2 ;...
            'Unc - 238; allT' 2 ;...
            'Unc - 238; 3T' 2 ;...
            'Unc - DEV; allT' 2 ;...
            'Unc - DEV; 3T' 2 ;...
            'parahippocampal - 238; allT' 2 ;...
            'parahippocampal - 238; 3T' 2 ;...
            'parahippocampal - DEV; allT' 2 ;...
            'parahippocampal - DEV; 3T' 2 ;...
            'entorhinal - 238; allT' 2 ;...
            'entorhinal - 238; 3T' 2 ;...
            'entorhinal - DEV; allT' 2 ;...
            'entorhinal - DEV; 3T' 2 };
        
% Out Table - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' 'log_mem_nor_scr_two_pst' ',' num2str(tbl_nme{iR,1})] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' 'vp2_nor_scr_pst'         ',' num2str(tbl_nme{iR,1})] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

% Out Table - OUTPUT (p-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_pvl = ejk_create_table( fcfg );

% Out Stat - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end; end; end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore' '/' 'Alena_postoperative_Laterality.csv'], out_tbl)

%% Table - Pre LTLE/RTLE/HC - Laterality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl

% Data Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_loc = { [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'TLE_Controls_pre_cln' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'subcort_vol_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'subcort_vol_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'subcort_vol_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'subcort_vol_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'fiber_FA_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'fiber_FA_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'fiber_FA_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'fiber_FA_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_238_LateralityIndex' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_238_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex'  ] };

% Table Construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_nme = { 'hippocampus' 2 ; ...
            'hippocampus' 3 ; ...
            'hippocampus' 4 ; ...
            'hippocampus' 5 ; ...
            'ILF' 6 ; ...
            'ILF' 7 ; ...
            'ILF' 8 ; ...
            'ILF' 9 ; ...
            'Unc' 6 ; ...
            'Unc' 7 ; ...
            'Unc' 8 ; ...
            'Unc' 9 ; ...
            'parahippocampal' 10 ; ...
            'parahippocampal' 11 ; ...
            'parahippocampal' 12 ; ...
            'parahippocampal' 13 ; ... 
            'entorhinal' 10 ; ...
            'entorhinal' 11 ; ...
            'entorhinal' 12 ; ...
            'entorhinal' 13 };
        
col_lbl = { '' 'Post LM2' 'Post VPA2' };
row_lbl = { 'Hippocampus - 238; allT' 2 ;...
            'Hippocampus - 238; 3T' 2 ;...
            'Hippocampus - DEV; allT' 2 ;...
            'Hippocampus - DEV; 3T' 2 ;...
            'ILF - 238; allT' 2 ;...
            'ILF - 238; 3T' 2 ;...
            'ILF - DEV; allT' 2 ;...
            'ILF - DEV; 3T' 2 ;...
            'Unc - 238; allT' 2 ;...
            'Unc - 238; 3T' 2 ;...
            'Unc - DEV; allT' 2 ;...
            'Unc - DEV; 3T' 2 ;...
            'parahippocampal - 238; allT' 2 ;...
            'parahippocampal - 238; 3T' 2 ;...
            'parahippocampal - DEV; allT' 2 ;...
            'parahippocampal - DEV; 3T' 2 ;...
            'entorhinal - 238; allT' 2 ;...
            'entorhinal - 238; 3T' 2 ;...
            'entorhinal - DEV; allT' 2 ;...
            'entorhinal - DEV; 3T' 2 };
        
% Out Table - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' 'log_mem_nor_scr_two' ',' num2str(tbl_nme{iR,1})] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' 'vp2_nor_scr'         ',' num2str(tbl_nme{iR,1})] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

% Out Table - OUTPUT (p-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_pvl = ejk_create_table( fcfg );

% Out Stat - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end; end; end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore' '/' 'Alena_preoperative_Laterality.csv'], out_tbl)

%% Table - Post LTLE - Raw  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl

% Data Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_loc = { [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'LTLE_post_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'LTLE_post_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'LTLE_post_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'LTLE_post_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'subcort_vol_238_norm_IntracranialVolume' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'subcort_vol_238_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'subcort_vol_dev_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'subcort_vol_dev_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'fiber_FA_238' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'fiber_FA_238'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'fiber_FA_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'fiber_FA_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_238' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_238'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'LTLE_post_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_dev'  ] };

% Table Construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_nme = { 'Educ' 3; ...
            'Educ' 4; ...
            'AgeOfSeizureOnset' 3; ...
            'AgeOfSeizureOnset' 4; ...
            'AgeAtSurgery' 3; ...
            'AgeAtSurgery' 4; ...
            'Left.hippocampus' 5 ; ...
            'Left.hippocampus' 6 ; ...
            'Left.hippocampus' 7 ; ...
            'Left.hippocampus' 8 ; ...
            'Right.hippocampus' 5 ; ...
            'Right.hippocampus' 6 ; ...
            'Right.hippocampus' 7 ; ...
            'Right.hippocampus' 8 ; ...
            'L.ILF' 9 ; ...
            'L.ILF' 10 ; ...
            'L.ILF' 11 ; ...
            'L.ILF' 12 ; ...
            'R.ILF' 9 ; ...
            'R.ILF' 10 ; ...
            'R.ILF' 11 ; ...
            'R.ILF' 12 ; ...
            'L.Unc' 9 ; ...
            'L.Unc' 10 ; ...
            'L.Unc' 11 ; ...
            'L.Unc' 12 ; ...
            'R.Unc' 9 ; ...
            'R.Unc' 10 ; ...
            'R.Unc' 11 ; ...
            'R.Unc' 12 ; ...
            'lh.parahippocampal' 13 ; ...
            'lh.parahippocampal' 14 ; ...
            'lh.parahippocampal' 15 ; ...
            'lh.parahippocampal' 16 ; ... 
            'rh.parahippocampal' 13 ; ...
            'rh.parahippocampal' 14 ; ...
            'rh.parahippocampal' 15 ; ...
            'rh.parahippocampal' 16 ; ... 
            'lh.entorhinal' 13 ; ...
            'lh.entorhinal' 14 ; ...
            'lh.entorhinal' 15 ; ...
            'lh.entorhinal' 16 ; ...
            'rh.entorhinal' 13 ; ...
            'rh.entorhinal' 14 ; ...
            'rh.entorhinal' 15 ; ...
            'rh.entorhinal' 16};
        
col_lbl = { '' 'Post LM2' 'Post VPA2' };
row_lbl = { 'Education - allT' 2 ;...
            'Education - 3T' 2 ;...
            'AgeOfSeizureOnset - allT' 2 ;...
            'AgeOfSeizureOnset - 3T' 2 ;...
            'AgeAtSurgery - allT' 2 ;...
            'AgeAtSurgery - 3T' 2 ;...
            'Left Hippocampus - 238; allT' 2 ;...
            'Left Hippocampus - 238; 3T' 2 ;...
            'Left Hippocampus - DEV; allT' 2 ;...
            'Left Hippocampus - DEV; 3T' 2 ;...
            'Right Hippocampus - 238; allT' 2 ;...
            'Right Hippocampus - 238; 3T' 2 ;...
            'Right Hippocampus - DEV; allT' 2 ;...
            'Right Hippocampus - DEV; 3T' 2 ;...
            'Left ILF - 238; allT' 2 ;...
            'Left ILF - 238; 3T' 2 ;...
            'Left ILF - DEV; allT' 2 ;...
            'Left ILF - DEV; 3T' 2 ;...
            'Right ILF - 238; allT' 2 ;...
            'Right ILF - 238; 3T' 2 ;...
            'Right ILF - DEV; allT' 2 ;...
            'Right ILF - DEV; 3T' 2 ;...
            'Left Unc - 238; allT' 2 ;...
            'Left Unc - 238; 3T' 2 ;...
            'Left Unc - DEV; allT' 2 ;...
            'Left Unc - DEV; 3T' 2 ;...
            'Right Unc - 238; allT' 2 ;...
            'Right Unc - 238; 3T' 2 ;...
            'Right Unc - DEV; allT' 2 ;...
            'Right Unc - DEV; 3T' 2 ; ...
            'Left parahippocampal - 238; allT' 2 ;...
            'Left parahippocampal - 238; 3T' 2 ;...
            'Left parahippocampal - DEV; allT' 2 ;...
            'Left parahippocampal - DEV; 3T' 2 ;...
             'Right parahippocampal - 238; allT' 2 ;...
            'Right parahippocampal - 238; 3T' 2 ;...
            'Right parahippocampal - DEV; allT' 2 ;...
            'Right parahippocampal - DEV; 3T' 2 ;...
            'Left entorhinal - 238; allT' 2 ;...
            'Left entorhinal - 238; 3T' 2 ;...
            'Left entorhinal - DEV; allT' 2 ;...
            'Left entorhinal - DEV; 3T' 2 ; ...
            'Right entorhinal - 238; allT' 2 ;...
            'Right entorhinal - 238; 3T' 2 ;...
            'Right entorhinal - DEV; allT' 2 ; ...
            'Right entorhinal - DEV; 3T' 2 ; };
        
% Out Table - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' 'log_mem_nor_scr_two_pst' ',' num2str(tbl_nme{iR,1})] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' 'vp2_nor_scr_pst'         ',' num2str(tbl_nme{iR,1})] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

% Out Table - OUTPUT (p-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_pvl = ejk_create_table( fcfg );

% Out Stat - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end; end; end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore' '/' 'Alena_postoperative_Raw.csv'], out_tbl)

%% Table - Pre LTLE/RTLE/HC - RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl

% Data Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_loc = { [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cog_cog_loc '/' 'TLE_Controls_pre_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_cln_cog_loc '/' 'TLE_Controls_pre_cog_3T' '/'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'subcort_vol_238_norm_IntracranialVolume' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'subcort_vol_238_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'subcort_vol_dev_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'subcort_vol_dev_norm_IntracranialVolume'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'fiber_FA_238' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'fiber_FA_238'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'fiber_FA_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'fiber_FA_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_238' ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_238'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_allT' '/' 'wmparc_FA_wm_aparc_annot_dev'  ] ...
            [ prj_dir '/' prj_nme '/' cor_neu_cog_loc '/' 'TLE_Controls_pre_cog_3T'   '/' 'wmparc_FA_wm_aparc_annot_dev'  ] };

% Table Construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_nme = { 'Educ' 3; ...
            'Educ' 4; ...
            'AgeOfSeizureOnset' 3; ...
            'AgeOfSeizureOnset' 4; ...
            'AgeAtSurgery' 3; ...
            'AgeAtSurgery' 4; ...
            'Left.hippocampus' 5 ; ...
            'Left.hippocampus' 6 ; ...
            'Left.hippocampus' 7 ; ...
            'Left.hippocampus' 8 ; ...
            'Right.hippocampus' 5 ; ...
            'Right.hippocampus' 6 ; ...
            'Right.hippocampus' 7 ; ...
            'Right.hippocampus' 8 ; ...
            'L.ILF' 9 ; ...
            'L.ILF' 10 ; ...
            'L.ILF' 11 ; ...
            'L.ILF' 12 ; ...
            'R.ILF' 9 ; ...
            'R.ILF' 10 ; ...
            'R.ILF' 11 ; ...
            'R.ILF' 12 ; ...
            'L.Unc' 9 ; ...
            'L.Unc' 10 ; ...
            'L.Unc' 11 ; ...
            'L.Unc' 12 ; ...
            'R.Unc' 9 ; ...
            'R.Unc' 10 ; ...
            'R.Unc' 11 ; ...
            'R.Unc' 12 ; ...
            'lh.parahippocampal' 13 ; ...
            'lh.parahippocampal' 14 ; ...
            'lh.parahippocampal' 15 ; ...
            'lh.parahippocampal' 16 ; ... 
            'rh.parahippocampal' 13 ; ...
            'rh.parahippocampal' 14 ; ...
            'rh.parahippocampal' 15 ; ...
            'rh.parahippocampal' 16 ; ... 
            'lh.entorhinal' 13 ; ...
            'lh.entorhinal' 14 ; ...
            'lh.entorhinal' 15 ; ...
            'lh.entorhinal' 16 ; ...
            'rh.entorhinal' 13 ; ...
            'rh.entorhinal' 14 ; ...
            'rh.entorhinal' 15 ; ...
            'rh.entorhinal' 16};
        
col_lbl = { '' 'Post LM2' 'Post VPA2' };
row_lbl = { 'Education - allT' 2 ;...
            'Education - 3T' 2 ;...
            'AgeOfSeizureOnset - allT' 2 ;...
            'AgeOfSeizureOnset - 3T' 2 ;...
            'AgeAtSurgery - allT' 2 ;...
            'AgeAtSurgery - 3T' 2 ;...
            'Left Hippocampus - 238; allT' 2 ;...
            'Left Hippocampus - 238; 3T' 2 ;...
            'Left Hippocampus - DEV; allT' 2 ;...
            'Left Hippocampus - DEV; 3T' 2 ;...
            'Right Hippocampus - 238; allT' 2 ;...
            'Right Hippocampus - 238; 3T' 2 ;...
            'Right Hippocampus - DEV; allT' 2 ;...
            'Right Hippocampus - DEV; 3T' 2 ;...
            'Left ILF - 238; allT' 2 ;...
            'Left ILF - 238; 3T' 2 ;...
            'Left ILF - DEV; allT' 2 ;...
            'Left ILF - DEV; 3T' 2 ;...
            'Right ILF - 238; allT' 2 ;...
            'Right ILF - 238; 3T' 2 ;...
            'Right ILF - DEV; allT' 2 ;...
            'Right ILF - DEV; 3T' 2 ;...
            'Left Unc - 238; allT' 2 ;...
            'Left Unc - 238; 3T' 2 ;...
            'Left Unc - DEV; allT' 2 ;...
            'Left Unc - DEV; 3T' 2 ;...
            'Right Unc - 238; allT' 2 ;...
            'Right Unc - 238; 3T' 2 ;...
            'Right Unc - DEV; allT' 2 ;...
            'Right Unc - DEV; 3T' 2 ; ...
            'Left parahippocampal - 238; allT' 2 ;...
            'Left parahippocampal - 238; 3T' 2 ;...
            'Left parahippocampal - DEV; allT' 2 ;...
            'Left parahippocampal - DEV; 3T' 2 ;...
             'Right parahippocampal - 238; allT' 2 ;...
            'Right parahippocampal - 238; 3T' 2 ;...
            'Right parahippocampal - DEV; allT' 2 ;...
            'Right parahippocampal - DEV; 3T' 2 ;...
            'Left entorhinal - 238; allT' 2 ;...
            'Left entorhinal - 238; 3T' 2 ;...
            'Left entorhinal - DEV; allT' 2 ;...
            'Left entorhinal - DEV; 3T' 2 ; ...
            'Right entorhinal - 238; allT' 2 ;...
            'Right entorhinal - 238; 3T' 2 ;...
            'Right entorhinal - DEV; allT' 2 ; ...
            'Right entorhinal - DEV; 3T' 2 ; };
        
% Out Table - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' 'log_mem_nor_scr_two' ',' num2str(tbl_nme{iR,1})] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' 'vp2_nor_scr'         ',' num2str(tbl_nme{iR,1})] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

% Out Table - OUTPUT (p-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_pvl = ejk_create_table( fcfg );

% Out Stat - OUTPUT (r-value) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end; end; end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore' '/' 'Alena_preoperative_Raw.csv'], out_tbl)

