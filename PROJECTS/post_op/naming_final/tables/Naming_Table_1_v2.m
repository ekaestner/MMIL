out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/Tables';

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv']);
dta_inp{2} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/ANOVA/Pre_Omnibus_anova/TLE.Controls.pre.anova/output_table.csv');
dta_inp{3} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/output_table.csv');

%%
fcfg = [];

fcfg.tbl(1,:) = { 'mean/std,1,controls_pre_3T_allSurg_all,AgeAtSurgery' ...
                 'mean/std,1,tle_pre_3T_allSurg_left,AgeAtSurgery' ...
                 'mean/std,1,tle_pre_3T_allSurg_right,AgeAtSurgery' ...
                 'copy,2,AgeAtSurgery,report'};

fcfg.tbl(2,:) = { 'mean/std,1,controls_pre_3T_allSurg_all,Educ' ...
                 'mean/std,1,tle_pre_3T_allSurg_left,Educ' ...
                 'mean/std,1,tle_pre_3T_allSurg_right,Educ' ...
                 'copy,2,Educ,report'};

fcfg.tbl(3,:) = { 'count,1,controls_pre_3T_allSurg_all,Sex,M/F' ...
                 'count,1,tle_pre_3T_allSurg_left,Sex,M/F' ...
                 'count,1,tle_pre_3T_allSurg_right,Sex,M/F' ...
                 'copy,3,Sex,report'};
             
fcfg.tbl(4,:) = { 'count,1,controls_pre_3T_allSurg_all,Handedness,R/L' ...
                 'count,1,tle_pre_3T_allSurg_left,Handedness,R/L' ...
                 'count,1,tle_pre_3T_allSurg_right,Handedness,R/L' ...
                 'copy,3,Handedness,report'};
             
fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%%
col_lbl = { '' 'HC' 'L-TLE' 'R-TLE' 'Test' };
row_lbl = { 'Age' 'Education' 'Sex (M/F)' 'Handedness (R/L)' }';

num_sbj = { 'N' numel(grp.('controls_pre_3T_allSurg_all')) numel(grp.('tle_pre_3T_allSurg_left')) numel(grp.('tle_pre_3T_allSurg_right')) ''};

out_tbl = [ col_lbl ; num_sbj ; row_lbl tbl_out ];

%%
cell2csv( [ out_dir '/' 'Table1.csv' ], out_tbl)







