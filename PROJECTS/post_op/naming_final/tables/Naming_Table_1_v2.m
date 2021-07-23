out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/tables/Table1';

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv']);
dta_inp{2} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/ANOVA/Pre_Omnibus_anova/TLE.Controls.pre.anova/output_table.csv');
dta_inp{3} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/output_table.csv');
dta_inp{4} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/ttest/Pre_TLE_ttest/tle.pre.3T.allSurg.left.right/output_table.csv');
dta_inp{5} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/Fisher/TLE_pre_fishers/output_table.csv');

%%
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'L') ,14) = {'Y'};
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'R') ,14) = {'Y'};
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'N/A') ,14) = {'N'};
dta_inp{1}( cellfun(@isempty,dta_inp{1}(:,14)) ,14) = {'N'};

%%
fcfg = [];

fcfg.tbl(1,:) = { 'mean/std,1,controls_pre_3T_allSurg_all,AgeAtImaging' ...
                 'mean/std,1,tle_pre_3T_allSurg_left,AgeAtImaging' ...
                 'mean/std,1,tle_pre_3T_allSurg_right,AgeAtImaging' ...
                 'copy,2,AgeAtImaging,report'};

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
             
fcfg.tbl(5,:) = { 'count,1,controls_pre_3T_allSurg_all,LanguageDominance,T/A' ...
                  'count,1,tle_pre_3T_allSurg_left,LanguageDominance,T/A' ...
                  'count,1,tle_pre_3T_allSurg_right,LanguageDominance,T/A' ...
                  'copy,3,LanguageDominance,report'};
              
fcfg.tbl(6,:) = { 'empty' ...
                  'mean/std,1,tle_pre_3T_allSurg_left,AgeOfSeizureOnset' ...
                  'mean/std,1,tle_pre_3T_allSurg_right,AgeOfSeizureOnset' ...
                  'copy,4,AgeOfSeizureOnset,report'};
             
fcfg.tbl(7,:) = { 'empty' ...
                  'mean/std,1,tle_pre_3T_allSurg_left,NumAEDs' ...
                  'mean/std,1,tle_pre_3T_allSurg_right,NumAEDs' ...
                  'copy,4,NumAEDs,report'};
             
fcfg.tbl(8,:) = { 'empty' ...
                  'count,1,tle_pre_3T_allSurg_left,MTS,Y/N' ...
                  'count,1,tle_pre_3T_allSurg_right,MTS,Y/N' ...
                  'copy,5,MTS,report'};
             
fcfg.tbl(9,:) = { 'empty' ...
                  'mean/std,1,tle_pre_3T_allSurg_left,SeizureFreq' ...
                  'mean/std,1,tle_pre_3T_allSurg_right,SeizureFreq' ...
                  'copy,4,SeizureFreq,report'};              
             
fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%%
col_lbl = { '' 'HC' 'L-TLE' 'R-TLE' 'Test' };
row_lbl = { 'Age' 'Education' 'Sex (M/F)' 'Handedness (R/L)' 'Dominance (T/A)'  'Age of Onset' '# current ASMs' 'MTS (Y/N)' 'Seizure Frequency' }';

num_sbj = { 'N' numel(grp.('controls_pre_3T_allSurg_all')) numel(grp.('tle_pre_3T_allSurg_left')) numel(grp.('tle_pre_3T_allSurg_right')) ''};

out_tbl = [ col_lbl ; num_sbj ; row_lbl tbl_out ];

%%
cell2csv( [ out_dir '/' 'Table1.csv' ], out_tbl)







