load([ dta_dir '/' 'demographics.mat' ]); % 'cmb_dta', 'cmb_sbj', 'cmb_col'
    cmb_dta(:,strcmpi(cmb_col,'Race/Ethnicity')) = [];
    cmb_col(strcmpi(cmb_col,'Race/Ethnicity'))   = [];

all_sbj_ind = find(sum(~cellfun(@isempty,cmb_dta),2)==numel(cmb_col),1);

out_dir_dem = [ out_dir '/' 'Tables' '/' 'Demographics' '/'];

%% Create grp variable
grp.main.HC  = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'Group')),'HC'));
grp.main.TLE = find(strcmpi(cmb_dta(:,strcmpi(cmb_col,'Group')),'TLE'));

%% Stats
% ttest2 Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_dta_col = find(cellfun(@isnumeric, cmb_dta(all_sbj_ind,:)));

use_dta_col(sum(cellfun(@isnan, cmb_dta(grp.main.HC,use_dta_col)))==numel(grp.main.HC)) = [];

fcfg = [];
fcfg.grp     = grp.main;
fcfg.grp_inc = {{'HC' 'TLE'}};
fcfg.grp_nme = {{'HC' 'TLE'}};
fcfg.dta = cmb_dta(:,use_dta_col);
fcfg.sbj = cmb_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

fcfg = [];
fcfg.sbj_nme = grp_sbj{1};
fcfg.dta     = grp_dta{1};
fcfg.dta_nme = cmb_col(:,use_dta_col);
fcfg.grp     = grp_typ{1};
fcfg.grp_nme = {'ttest'};
fcfg.out_dir = [ out_dir_dem '/' 'stats' '/' ];

fcfg.dta_nme = mmil_spec_char(fcfg.dta_nme,{' '},{'_'});

ejk_ttest2_independent( fcfg );

% fisher's Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_dta_col = ~cellfun(@isnumeric, cmb_dta(all_sbj_ind,:));

fcfg = [];
fcfg.grp     = grp.main;
fcfg.grp_inc = {{'HC' 'TLE'}};
fcfg.grp_nme = {{'HC' 'TLE'}};
fcfg.dta = cmb_dta(:,use_dta_col);
fcfg.sbj = cmb_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};

fcfg = [];
fcfg.sbj = grp_sbj{1};
fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cmb_col(:,use_dta_col);
fcfg.dta_two = repmat(grp_typ{1},1,sum(use_dta_col));
fcfg.lbl_two = strcat( 'group_', cmb_col(:,use_dta_col));
fcfg.grp_nme = 'fisher';
fcfg.out_dir = [ out_dir_dem '/' 'stats' '/' ];

fcfg.lbl_one = mmil_spec_char(fcfg.lbl_one,{' '},{'_'});
fcfg.lbl_two = mmil_spec_char(fcfg.lbl_two,{' '},{'_'});

ejk_fisher_test( fcfg );

%% Table design
grp_typ = {'HC' 'TLE'};
col_hld = cmb_col';
for iR = 1:numel(col_hld)
    if isnumeric(cmb_dta{all_sbj_ind,iR})
        col_hld{iR,2} = 'num';
    elseif ~isnumeric(cmb_dta{all_sbj_ind,iR})
        col_hld{iR,2} = 'str';
    end
end

%% Make Table
tst_tbl = mmil_readtext([ out_dir_dem '/' 'stats' '/' 'ttest'  '/' 'output_table.csv']);
fsh_tbl = mmil_readtext([ out_dir_dem '/' 'stats' '/' 'fisher' '/' 'output_table.csv']);

% Setup table
clear tbl_dsg
for iR = 1:numel(col_hld(:,1))
    for iC = 1:numel(grp_typ)
        if strcmpi(col_hld{iR,end},'num')
            tbl_dsg{iR,iC} = [ 'mean/std' ','  '1' ',' grp_typ{iC} ',' col_hld{iR,1} ];
            row_lbl{iR} = col_hld{iR,1};
        elseif strcmpi(col_hld{iR,end},'str')
            cmb_dta(cellfun(@isnumeric,cmb_dta(:,strcmpi(cmb_col,col_hld{iR,1}))),strcmpi(cmb_col,col_hld{iR,1})) = {''};            
            hld_cat = unique(cmb_dta(:,strcmpi(cmb_col,col_hld{iR,1})));
            hld_cat(cellfun(@isempty,hld_cat)) = [];
            hld_cat = strcat(hld_cat,'/');
            hld_cat = cat(2,hld_cat{:});
            hld_cat = hld_cat(1:end-1);
            tbl_dsg{iR,iC} = [ 'count'    ','  '1' ',' grp_typ{iC} ',' col_hld{iR,1} ','  hld_cat];
            row_lbl{iR} = [ col_hld{iR,1} ' (' hld_cat ')'];
        end
    end
    if strcmpi(col_hld{iR,end},'num')
        tbl_dsg{iR,iC+1} = ['copy' ',' '2' ',' col_hld{iR,1} ',' 'report'];
    elseif strcmpi(col_hld{iR,end},'str')
        tbl_dsg{iR,iC+1} = ['copy' ',' '3' ',' col_hld{iR,1} ',' 'report' ];
    end
end

fcfg = [];
fcfg.tbl = tbl_dsg;
fcfg.dta = {[ cmb_col ; cmb_dta ] tst_tbl fsh_tbl };
fcfg.grp = grp.main;
tbl_out = ejk_create_table( fcfg );

%% Save out
cell2csv([ out_dir_dem '/' 'demographic_table.csv'],[ {''} grp_typ {'Stats'} ;  col_hld(:,1) tbl_out ]);









