
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
tbl_out = cell( 13, 5);

grp_nme = { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

tbl_out(1,:) = { '' 'L-TLE' 'R-TLE' 'statistic' 'p' };

% N %%%%%%%%%%%%

tbl_out{2,1} = 'N';
tbl_out{2,2} = numel(grp.(grp_nme{1}));
tbl_out{2,3} = numel(grp.(grp_nme{2}));

% Age %%%%%%%%%%%%
row_num = 3;
col_num = 4;

tbl_out{row_num,1} = 'Age';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );

[ ~, pvl, ~, tst_hld ] = ttest2( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)), cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = [ 't(' num2str(tst_hld.df) ')=' num2str(roundsd(tst_hld.tstat,2))];
tbl_out{row_num,5} = num2str(roundsd(pvl,2));
    tbl_out{row_num,5} = tbl_out{row_num,5}(2:end);
                       
% Education %%%%%%%%%%%%
row_num = 4;
col_num = 5;

tbl_out{row_num,1} = 'Education';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );

[ ~, pvl, ~, tst_hld ] = ttest2( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)), cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = [ 't(' num2str(tst_hld.df) ')=' num2str(roundsd(tst_hld.tstat,2))];
tbl_out{row_num,5} = num2str(roundsd(pvl,2));
    tbl_out{row_num,5} = tbl_out{row_num,5}(2:end);
    
% Sex %%%%%%%%%%%%
row_num = 5;
col_num = 11;

tbl_out{row_num,1} = 'Sex (M/F)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'M'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'F'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'M'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'F'))) ];

tbl_out{row_num,4} = '';
tbl_out{row_num,5} = '1.0';
    
% Handedness %%%%%%%%%%%%
row_num = 6;
col_num = 12;

tbl_out{row_num,1} = 'Hand (R/L)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'R'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) '/' num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'L'))) ];

tbl_out{row_num,4} = '';
tbl_out{row_num,5} = '.34';

% Age of Onset %%%%%%%%%%%%
row_num = 7;
col_num = 6;

tbl_out{row_num,1} = 'Age of Onset';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );

[ ~, pvl, ~, tst_hld ] = ttest2( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)), cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = [ 't(' num2str(tst_hld.df) ')=' num2str(roundsd(tst_hld.tstat,2))];
tbl_out{row_num,5} = num2str(roundsd(pvl,2));
    tbl_out{row_num,5} = tbl_out{row_num,5}(2:end);

% ASMs %%%%%%%%%%%%
row_num = 8;
col_num = 7;

tbl_out{row_num,1} = '# ASMs';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );

[ ~, pvl, ~, tst_hld ] = ttest2( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)), cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = [ 't(' num2str(tst_hld.df) ')=' num2str(roundsd(tst_hld.tstat,2))];
tbl_out{row_num,5} = num2str(roundsd(pvl,2));
    tbl_out{row_num,5} = tbl_out{row_num,5}(2:end);

% MTS %%%%%%%%%%%%
row_num = 9;
col_num = 13;

tbl_out{row_num,1} = 'MTS (Y/N)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) '/' num2str(sum(~strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) '/' num2str(sum(~strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) ];

tbl_out{row_num,4} = '';
tbl_out{row_num,5} = '.34';

% Seizure Frequency %%%%%%%%%%%%
row_num = 10;
col_num = 8;

tbl_out{row_num,1} = 'Seizure Frequency';
tbl_out{row_num,2} = nanmean( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)) );
tbl_out{row_num,3} = nanmean( cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );

[ ~, pvl, ~, tst_hld ] = ttest2( cell2mat(cln_dta( grp.(grp_nme{1}), col_num)), cell2mat(cln_dta( grp.(grp_nme{2}), col_num)) );
tbl_out{row_num,4} = [ 't(' num2str(tst_hld.df) ')=' num2str(roundsd(tst_hld.tstat,2))];
tbl_out{row_num,5} = num2str(roundsd(pvl,2));
    tbl_out{row_num,5} = tbl_out{row_num,5}(2:end);

% Engel Outcome %%%%%%%%%%%%
row_num = 10;
col_num = 13;

tbl_out{row_num,1} = 'MTS (Y/N)';
tbl_out{row_num,2} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) '/' num2str(sum(~strcmpi(cln_dta( grp.(grp_nme{1}), col_num),'L'))) ];
tbl_out{row_num,3} = [ num2str(sum(strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) '/' num2str(sum(~strcmpi(cln_dta( grp.(grp_nme{2}), col_num),'R'))) ];

tbl_out{row_num,4} = '';
tbl_out{row_num,5} = '.34';

% Race %%%%%%%%%%%%
row_num = 11;

tbl_out{row_num,1} = 'Race';

% Ethnicity %%%%%%%%%%%%
row_num = 12;

tbl_out{row_num,1} = 'Ethnicity';

% WTAR %%%%%%%%%%%%
row_num = 13;

tbl_out{row_num,1} = 'WTAR';

%%
cell2csv( [ out_dir '/' 'Table1_Preoperative.csv'], tbl_out)
