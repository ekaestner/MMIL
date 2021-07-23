load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/figures/SupplementaryFigure_2' ]; ejk_chk_dir( out_dir );

cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];

%% Load Data
cog_dta = mmil_readtext(cog_dta_nme);
cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

%% Pre-Operative - all sample
ydt{1} = cell2mat( cog_dta( grp.controls_pre_3T_allSurg_all, 1));  xdt{1} = 0.5; % L-TLE BNT
ydt{2} = cell2mat( cog_dta( grp.tle_pre_3T_allSurg_left, 1));      xdt{2} = 1; % L-TLE BNT
ydt{3} = cell2mat( cog_dta( grp.tle_pre_3T_allSurg_right, 1));     xdt{3} = 1.5; % R-TLE BNT

ydt{4} = cell2mat( cog_dta( grp.controls_pre_3T_allSurg_all, 2));  xdt{4} = 0.5; % L-TLE ANT
ydt{5} = cell2mat( cog_dta( grp.tle_pre_3T_allSurg_left, 2));      xdt{5} = 1.0; % L-TLE ANT
ydt{6} = cell2mat( cog_dta( grp.tle_pre_3T_allSurg_right, 2));     xdt{6} = 1.5; % R-TLE ANT

% BNT
fcfg = [];
        
fcfg.xdt     = xdt(1:3);
fcfg.ydt     = ydt(1:3);

fcfg.box_plt     = [ 1 1 1 ];
fcfg.box_plt_col = {  rgb('grey') rgb('royal purple') rgb('light magenta') };

fcfg.fce_col = { rgb('black') rgb('black') rgb('black') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlm = [0 2];

fcfg.xlb = { 'HC' 'LTLE' 'RTLE' };
fcfg.ylb = { 'pre-operative score' };

fcfg.ttl = ['BNT'];

fcfg.hln     = -1.5;
fcfg.hln_col = rgb('light red');

fcfg.out_dir = out_dir;
fcfg.out_nme = 'Pre_totalsample_Cognitive_bnt';

ejk_scatter(fcfg)

% ANT
fcfg = [];
        
fcfg.xdt     = xdt(4:6);
fcfg.ydt     = ydt(4:6);

fcfg.box_plt     = [ 1 1 1];
fcfg.box_plt_col = { rgb('grey') rgb('royal purple') rgb('light magenta') };

fcfg.fce_col = { rgb('black') rgb('black') rgb('black') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlm = [0 2];

fcfg.xlb = { 'HC' 'LTLE' 'RTLE' };
fcfg.ylb = { 'pre-operative score' };

fcfg.ttl = ['ANT'];

fcfg.hln     = -1.5;
fcfg.hln_col = rgb('light red');

fcfg.out_dir = out_dir;
fcfg.out_nme = 'Pre_totalsample_Cognitive_ant';

ejk_scatter(fcfg)

%% Post-Operative
clear ydt xdt

ydt{1} = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, 3));  xdt{1} = 0.7; % L-TLE BNT
ydt{2} = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_right, 3));  xdt{2} = 1.3; % R-TLE BNT
ydt{3} = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_left, 4)); xdt{3} = 0.7; % L-TLE ANT
ydt{4} = cell2mat( cog_dta( grp.tle_post_3T_ATLonly_right, 4)); xdt{4} = 1.3; % R-TLE ANT

%
% BNT
fcfg = [];
        
fcfg.xdt     = xdt(1:2);
fcfg.ydt     = ydt(1:2);

fcfg.box_plt     = [ 1 1 ];
fcfg.box_plt_col = {  rgb('royal purple') rgb('light magenta') };

fcfg.fce_col = { rgb('black') rgb('black') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlm = [0 2];

fcfg.xlb = { 'LTLE' 'RTLE' };
fcfg.ylb = { 'post-operative score' };

fcfg.ttl = ['BNT'];

fcfg.hln     = -1.5;
fcfg.hln_col = rgb('light red');

fcfg.out_dir = out_dir;
fcfg.out_nme = 'Post_Cognitive_bnt';

ejk_scatter(fcfg)

% ANT
fcfg = [];
        
fcfg.xdt     = xdt(3:4);
fcfg.ydt     = ydt(3:4);

fcfg.box_plt     = [  1 1];
fcfg.box_plt_col = { rgb('royal purple') rgb('light magenta') };

fcfg.fce_col = { rgb('black') rgb('black') };
fcfg.edg_col = { rgb('black') rgb('black') };

fcfg.xlm = [0 2];

fcfg.xlb = { 'LTLE' 'RTLE' };
fcfg.ylb = { 'post-operative score' };

fcfg.ttl = ['ANT'];

fcfg.hln     = -1.5;
fcfg.hln_col = rgb('light red');

fcfg.out_dir = out_dir;
fcfg.out_nme = 'Post_Cognitive_ant';

ejk_scatter(fcfg)




