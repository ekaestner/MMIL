
%% UCSF vs UCSD
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

numel(string_find( [cln_dta_sbj(grp.tle_post_3T_ATLonly_left);cln_dta_sbj(grp.tle_post_3T_ATLonly_right)], {'ucsf'}))
numel(string_find( [cln_dta_sbj(grp.tle_post_3T_ATLonly_left);cln_dta_sbj(grp.tle_post_3T_ATLonly_right)], {'epd'}))
numel(cln_dta_sbj(grp.tle_post_3T_ATLonly_left))
numel(cln_dta_sbj(grp.tle_post_3T_ATLonly_right))

%% Patients with both tests
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

(sum(( ~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,1)) & ...
       ~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,2)))) / numel(grp.tle_post_3T_ATLonly_all))*100

(sum(( ~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,3)) & ...
       ~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,4)))) / numel(grp.tle_post_3T_ATLonly_all))*100
   

setxor(cog_dta_sbj(~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,3))), cog_dta_sbj(~cellfun(@isnan,cog_dta(grp.tle_post_3T_ATLonly_all,4))))
[cog_dta_sbj(grp.tle_post_3T_ATLonly_left) cog_dta( grp.tle_post_3T_ATLonly_left, 3:4 )]
[cog_dta_sbj(grp.tle_post_3T_ATLonly_right) cog_dta( grp.tle_post_3T_ATLonly_right, 3:4 )]

%% Distance between imaging and surgery
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );


