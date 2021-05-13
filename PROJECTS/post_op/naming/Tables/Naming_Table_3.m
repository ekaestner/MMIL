out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Tables/Parts';

cor_dir = 'SpecificCor';

cog_dir = [ prj_dir '/' prj_nme '/' cor_dir '/' 'Cognitive' '/' ];

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];
cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];

pre_opr_dir         = 'TLE_Controls_pre_pre';
pst_opr_lft_tle_dir = 'LTLE_post_post';
pst_opr_rgh_tle_dir = 'RTLE_post_post';

pre_anv_dir = 'Pre_Omnibus_anova';
    pre_anv_cmp = 'HC.TLE';
pst_tts_dir = 'Post_TLE_ttest';
    pst_tts_cmp = 'TLE.post';

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

cog_dta = mmil_readtext(cog_dta_nme);
cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cell2mat(cog_dta(2:end,2:end));

%% Pre-operative scores
anv_tbl = mmil_readtext( [ prj_dir '/' prj_nme '/' cor_dir '/' 'Cognitive' '/' pre_anv_dir '/' pre_anv_cmp '/' 'output_table.csv'] );

tbl_anv_out{1,1} = '';
tbl_anv_out{1,2} = 'HC';
tbl_anv_out{1,3} = 'L-TLE';
tbl_anv_out{1,4} = 'R-TLE';
tbl_anv_out{1,5} = 'ANOVA';
tbl_anv_out{1,6} = 'L-TLE / HC';
tbl_anv_out{1,7} = 'R-TLE / HC';
tbl_anv_out{1,8} = 'L-TLE / R-TLE';

% Preoperative BNT Scores (L- vs R-TLE)
tbl_anv_out{2,1} = 'Pre BNT';
tbl_anv_out{2,2} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_controls_pre_3T_allSurg_all,1)),3)) ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_controls_pre_3T_allSurg_all,1)),2)) ')' ];
tbl_anv_out{2,3} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_left,1)),3))        ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_left,1)),2)) ')' ];
tbl_anv_out{2,4} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_right,1)),3))       ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_right,1)),2)) ')' ];
tbl_anv_out{2,5} = [anv_tbl{2,6} ',' anv_tbl{2,7} ', ' anv_tbl{2,8} ];
tbl_anv_out{2,6} = ['p = ' num2str(roundsd(anv_tbl{2,16},2)) ];
tbl_anv_out{2,7} = ['p = ' num2str(roundsd(anv_tbl{2,17},2)) ];
tbl_anv_out{2,8} = ['p = ' num2str(roundsd(anv_tbl{2,18},2)) ];

% Preoperative ANT Scores (L- vs R-TLE)
tbl_anv_out{3,1} = 'Pre ANT';
tbl_anv_out{3,2} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_controls_pre_3T_allSurg_all,2)),3)) ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_controls_pre_3T_allSurg_all,2)),2)) ')' ];
tbl_anv_out{3,3} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_left,2)),3))        ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_left,2)),2)) ')' ];
tbl_anv_out{3,4} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_right,2)),3))       ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_right,2)),2)) ')' ];
tbl_anv_out{3,5} = [anv_tbl{3,6} ',' anv_tbl{3,7} ', ' anv_tbl{3,8} ];
tbl_anv_out{3,6} = ['p = ' num2str(roundsd(anv_tbl{3,16},2)) ];
tbl_anv_out{3,7} = ['p = ' num2str(roundsd(anv_tbl{3,17},2)) ];
tbl_anv_out{3,8} = ['p = ' num2str(roundsd(anv_tbl{3,18},2)) ];

cell2csv([ out_dir '/' 'table1_preoperative.csv'], tbl_anv_out, ';')
clear tbl_anv_out

%% Post-operative scores
tts_tbl = mmil_readtext( [ prj_dir '/' prj_nme '/' cor_dir '/' 'Cognitive' '/' pst_tts_dir '/' pst_tts_cmp '/' 'output_table.csv'] );

tbl_tts_out{1,1} = '';
tbl_tts_out{1,2} = 'L-TLE';
tbl_tts_out{1,3} = 'R-TLE';
tbl_tts_out{1,4} = 't-test';

% Post-operative BNT Scores (L- vs R-TLE)
tbl_tts_out{2,1} = 'Post BNT';
tbl_tts_out{2,2} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_left,4)),3))  ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_left,4)),2)) ')' ];
tbl_tts_out{2,3} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_right,4)),3)) ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_right,4)),2)) ')' ];
tbl_tts_out{2,4} = [ tts_tbl{2,7} ', ' tts_tbl{2,8} ];

[ ~, pvl ] = ttest(cog_dta(grp.tle_post_3T_ATLonly_left,4),0);
    if pvl<.05; tbl_tts_out{2,2} = [tbl_tts_out{2,2} '*']; end
    if pvl<.01; tbl_tts_out{2,2} = [tbl_tts_out{2,2} '*']; end
[ ~, pvl ] = ttest(cog_dta(grp.tle_post_3T_ATLonly_right,4),0);
    if pvl<.05; tbl_tts_out{2,3} = [tbl_tts_out{2,3} '*']; end
    if pvl<.01; tbl_tts_out{2,3} = [tbl_tts_out{2,3} '*']; end

tbl_tts_out{3,1} = 'BNT Decliners';
tbl_tts_out{3,2} = [ num2str( sum(cog_dta(grp.tle_post_3T_ATLonly_left,4)  <=-1.5)) ' ('...
    num2str(roundsd(sum(cog_dta(grp.tle_post_3T_ATLonly_left,4)  <=-1.5) / sum(~isnan(cog_dta(grp.tle_post_3T_ATLonly_left,4)))*100,2)) '%)'];
tbl_tts_out{3,3} = [ num2str( sum(cog_dta(grp.tle_post_3T_ATLonly_right,4)  <=-1.5)) ' ('...
    num2str(roundsd(sum(cog_dta(grp.tle_post_3T_ATLonly_right,4)  <=-1.5) / sum(~isnan(cog_dta(grp.tle_post_3T_ATLonly_right,4)))*100,2)) '%)'];


% Post-operative ANT Scores (L- vs R-TLE)
tbl_tts_out{4,1} = 'Post ANT';
tbl_tts_out{4,2} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_left,5)),3))  ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_left,5)),2)) ')' ];
tbl_tts_out{4,3} = [ num2str(roundsd(nanmean(cog_dta(grp.tle_post_3T_ATLonly_right,5)),3)) ' (' num2str(roundsd(nanstd(cog_dta(grp.tle_post_3T_ATLonly_right,5)),2)) ')' ];
tbl_tts_out{4,4} = [ tts_tbl{3,7} ', ' tts_tbl{3,8} ];

[ ~, pvl ] = ttest(cog_dta(grp.tle_post_3T_ATLonly_left,5),0);
    if pvl<.05; tbl_tts_out{3,2} = [tbl_tts_out{3,2} '*']; end
    if pvl<.01; tbl_tts_out{3,2} = [tbl_tts_out{3,2} '*']; end
[ ~, pvl ] = ttest(cog_dta(grp.tle_post_3T_ATLonly_right,5),0);
    if pvl<.05; tbl_tts_out{3,3} = [tbl_tts_out{3,3} '*']; end
    if pvl<.01; tbl_tts_out{3,3} = [tbl_tts_out{3,3} '*']; end
    
tbl_tts_out{5,1} = 'ANT Decliners';
tbl_tts_out{5,2} = [ num2str( sum(cog_dta(grp.tle_post_3T_ATLonly_left,5)  <=-1.5)) ' ('...
    num2str(roundsd(sum(cog_dta(grp.tle_post_3T_ATLonly_left,5)  <=-1.5) / sum(~isnan(cog_dta(grp.tle_post_3T_ATLonly_left,5)))*100,2)) '%)'];
tbl_tts_out{5,3} = [ num2str( sum(cog_dta(grp.tle_post_3T_ATLonly_right,5)  <=-1.5)) ' ('...
    num2str(roundsd(sum(cog_dta(grp.tle_post_3T_ATLonly_right,5)  <=-1.5) / sum(~isnan(cog_dta(grp.tle_post_3T_ATLonly_right,5)))*100,2)) '%)'];

cell2csv([ out_dir '/' 'table2_postoperative.csv'], tbl_tts_out, ';')
clear tbl_tts_out

