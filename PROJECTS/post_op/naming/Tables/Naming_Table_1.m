
%%
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Tables/Parts';

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

cln_dta_nme = [ dta_dir '/' 'Clinical'                          '.csv'];

cln_dta = mmil_readtext(cln_dta_nme);
cln_dta_col = ejk_fix_column_names(cln_dta(1,2:end));
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

%%
tbl_out = cell( 9, 6);

grp_nme = { 'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right' };

tbl_out(1,:) = { '' 'HC' 'L-TLE' 'R-TLE' 'statistic' 'p' };

% N %%%%%%%%%%%%

tbl_out{2,1} = 'N';
tbl_out{2,2} = numel(grp.(grp_nme{1}));
tbl_out{2,3} = numel(grp.(grp_nme{2}));
tbl_out{2,4} = numel(grp.(grp_nme{3}));

% Age %%%%%%%%%%%%
row_num = 3;
col_num = 4;

tbl_out{row_num,1} = 'Age';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = nanmean( cell2mat(cln_dta( grp.(grp_nme{3}), col_num)) );

[ pvl, anv_hld ] = anova1( [ cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) ; cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) ; cell2mat(cln_dta( grp.(grp_nme{3}), col_num)) ], ...
                           [ repmat({'HC'},tbl_out{2,2},1) ; repmat({'L-TLE'},tbl_out{2,3},1) ; repmat({'R-TLE'},tbl_out{2,4},1) ] );
tbl_out{row_num,5} = [ 'F(' num2str(anv_hld{2,3}) ',' num2str(anv_hld{3,3}) ')=' num2str(roundsd(anv_hld{2,5},2))];
tbl_out{row_num,6} = num2str(roundsd(pvl,2));
    tbl_out{row_num,6} = tbl_out{row_num,6}(2:end);
                       
% Education %%%%%%%%%%%%
row_num = 4;
col_num = 5;

tbl_out{row_num,1} = 'Education';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = nanmean( cell2mat(cln_dta( grp.(grp_nme{3}), col_num)) );

[ pvl, anv_hld ] = anova1( [ cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) ; cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) ; cell2mat(cln_dta( grp.(grp_nme{3}), col_num)) ], ...
                           [ repmat({'HC'},tbl_out{2,2},1) ; repmat({'L-TLE'},tbl_out{2,3},1) ; repmat({'R-TLE'},tbl_out{2,4},1) ] );
tbl_out{row_num,5} = [ 'F(' num2str(anv_hld{2,3}) ',' num2str(anv_hld{3,3}) ')=' num2str(roundsd(anv_hld{2,5},2))];
tbl_out{row_num,6} = num2str(roundsd(pvl,2));
    tbl_out{row_num,6} = tbl_out{row_num,6}(2:end);
    
% Sex %%%%%%%%%%%%
row_num = 5;
col_num = 11;

tbl_out{row_num,1} = 'Sex (M/F)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'M'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'F'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'M'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'F'))) ];
tbl_out{row_num,4} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{3}), col_num),'M'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{3}), col_num),'F'))) ];

tbl_out{row_num,5} = '';
tbl_out{row_num,6} = '.74';
    
% Handedness %%%%%%%%%%%%
row_num = 6;
col_num = 12;

tbl_out{row_num,1} = 'Hand (R/L)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'R'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'L'))) ];
tbl_out{row_num,4} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{3}), col_num),'R'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{3}), col_num),'L'))) ];

tbl_out{row_num,5} = '';
tbl_out{row_num,6} = '.098';

% Race %%%%%%%%%%%%
row_num = 7;

tbl_out{row_num,1} = 'Race';

% Ethnicity %%%%%%%%%%%%
row_num = 8;

tbl_out{row_num,1} = 'Ethnicity';

% WTAR %%%%%%%%%%%%
row_num = 9;

tbl_out{row_num,1} = 'WTAR';

%%
cell2csv( [ out_dir '/' 'Table1_Preoperative.csv'], tbl_out)
