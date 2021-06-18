
%% Subject #'s
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

sbj_num(1) = numel( grp.controls_pre_3T_allSurg_all ); % 
sbj_num(2) = numel( grp.tle_pre_3T_allSurg_left ); % 
sbj_num(3) = numel( grp.tle_pre_3T_allSurg_right ); % 
sbj_num(4) = numel( grp.tle_post_3T_ATLonly_left ); % 
sbj_num(5) = numel( grp.tle_post_3T_ATLonly_right ); % 

sbj_num