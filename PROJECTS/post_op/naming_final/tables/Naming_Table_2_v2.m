out_dir =  [ prj_dir '/' prj_nme '/' 'Tables' ];

clear out_tbl

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv']);
dta_inp{2} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/ttest/Post_TLE_ttest/TLE.post.ttest/output_table.csv');
dta_inp{3} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_post_fishers/output_table.csv');

%%
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'L') ,14) = {'Y'};
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'R') ,14) = {'Y'};
dta_inp{1}( strcmpi(dta_inp{1}(:,14),'N/A') ,14) = {'N'};
dta_inp{1}( cellfun(@isempty,dta_inp{1}(:,14)) ,14) = {'N'};

dta_inp{1}( strcmpi(dta_inp{1}(:,15),'II') ,15) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,15),'III') ,15) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,15),'IV') ,15) = {'II+'};

%%
fcfg = [];

fcfg.tbl(1,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,AgeAtSurgery' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,AgeAtSurgery' ...
                 'copy,2,AgeAtSurgery,report'};

fcfg.tbl(2,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,Educ' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,Educ' ...
                 'copy,2,Educ,report'};

fcfg.tbl(3,:) = { 'count,1,tle_post_3T_ATLonly_left,Sex,M/F' ...
                 'count,1,tle_post_3T_ATLonly_right,Sex,M/F' ...
                 'copy,3,Sex,report'};
             
fcfg.tbl(4,:) = { 'count,1,tle_post_3T_ATLonly_left,Handedness,R/L' ...
                 'count,1,tle_post_3T_ATLonly_right,Handedness,R/L' ...
                 'copy,3,Handedness,report'};
  
fcfg.tbl(5,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,AgeOfSeizureOnset' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,AgeOfSeizureOnset' ...
                 'copy,2,AgeOfSeizureOnset,report'};
             
fcfg.tbl(6,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,NumAEDs' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,NumAEDs' ...
                 'copy,2,NumAEDs,report'};
             
fcfg.tbl(7,:) = { 'count,1,tle_post_3T_ATLonly_left,MTS,Y/N' ...
                 'count,1,tle_post_3T_ATLonly_right,MTS,Y/N' ...
                 'copy,3,MTS,report'};
             
fcfg.tbl(8,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,SeizureFreq' ...
                 'mean/std,1,tle_post_3T_ATLonly_right,SeizureFreq' ...
                 'copy,2,SeizureFreq,report'};
 
fcfg.tbl(9,:) = { 'count,1,tle_post_3T_ATLonly_left,EngelOutcome,I/II+' ...
                 'count,1,tle_post_3T_ATLonly_right,EngelOutcome,I/II+' ...
                 'copy,3,EngelOutcome,report'};
             
fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

%%
col_lbl = { '' 'L-TLE' 'R-TLE' 'Test' };
row_lbl = { 'Age' 'Education' 'Sex (M/F)' 'Handedness (R/L)' 'Age of Onset' '# current ASMs' 'MTS (Y/N)' 'Seizure Frequency' 'Engel Outcome (I/II+)' }';

num_sbj = { 'N' numel(grp.('tle_post_3T_ATLonly_left')) numel(grp.('tle_post_3T_ATLonly_right')) ''};

out_tbl = [ col_lbl ; num_sbj ; row_lbl tbl_out ];

%%
cell2csv( [ out_dir '/' 'Table2.csv' ], out_tbl)













