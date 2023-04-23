
%% Load
fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_TLE.csv' ];
[ tle_dem, tle_dem_sbj, tle_dem_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_HC.csv' ];
[ con_dem, con_dem_sbj, con_dem_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'predictions.csv' ];
[ prd_dta, prd_dta_sbj, prd_dta_col ] = ejk_dta_frm(fcfg);

%% Demographics: Make sure tables match, add duration
% Calculate Duration
sze_ons_col = strcmpi(tle_dem_col,'Age at Onset');
age_col     = strcmpi(tle_dem_col,'Age');
dur_col = numel(tle_dem_col)+1;

tle_dem_col = [ tle_dem_col 'Duration' ];
for iS = 1:numel(tle_dem_sbj)
    if ~isempty(tle_dem{iS,sze_ons_col}) && ~isempty(tle_dem{iS,age_col})
        tle_dem{iS,dur_col} = tle_dem{iS,age_col} - tle_dem{iS,sze_ons_col};
    else
        tle_dem{iS,dur_col} = NaN;     
    end
end

% Combine EPD/HC
con_dem_col = [ con_dem_col 'Group' ];
con_dem     = [ con_dem repmat({'HC'},numel(con_dem_sbj),1) ];

tle_dem_col = [ tle_dem_col 'Group' ];
tle_dem     = [ tle_dem repmat({'TLE'},numel(tle_dem_sbj),1) ];

cmb_col = unique([tle_dem_col con_dem_col]);
cmb_sbj = [ tle_dem_sbj ; con_dem_sbj ];
cmb_dta = cell( numel(cmb_sbj), numel(cmb_col) );
for iC = 1:numel(cmb_col)
    con_col_ind = strcmpi(con_dem_col,cmb_col{iC});
    tle_col_ind = strcmpi(tle_dem_col,cmb_col{iC});
    if sum(con_col_ind)==0
        cmb_dta(:,iC) = [ tle_dem(:,tle_col_ind) ; num2cell(nan(numel(con_dem_sbj),1)) ];
    elseif sum(tle_col_ind)==0
        error('line 63')
    else
        cmb_dta(:,iC) = [ tle_dem(:,tle_col_ind) ; con_dem(:,con_col_ind) ];
    end
end

% Check tables match
emp_cnt = 0;
cmb_dta_use = cell( numel(prd_dta_sbj),  numel(cmb_col));
for iS = 1:numel(prd_dta_sbj)
    if ~sum(strcmpi(cmb_sbj,prd_dta_sbj{iS}))
        emp_cnt = emp_cnt + 1;
    else
        cmb_dta_use(iS,:) = cmb_dta( strcmpi(cmb_sbj,prd_dta_sbj{iS}),:);
    end
end
cmb_dta = cmb_dta_use;
cmb_sbj = prd_dta_sbj;

%% Recode Values
rcd_nme    = { 'Gender' 'Handedness' 'Side of Epilepsy' };
rcd_val{1} = { 0 'Male' ; 1 'Female' };
rcd_val{2} = { 0 'Left' ; 1 'Right'  };
rcd_val{3} = { 1 'LTLE' ; 2 'RTLE'  };

fcfg = [];

fcfg.dta     = cmb_dta;
fcfg.dta_col = cmb_col;

fcfg.rcd_nme = rcd_nme;
fcfg.rcd_val = rcd_val;

fcfg.swt_val = 1;

cmb_dta_rcd = ejk_recode(fcfg);

cmb_dta = cmb_dta_rcd;

%% Prediction: Calculate Performances
prd_scr = nan(numel(prd_dta_sbj),1);
for iS = 1:numel(prd_dta_sbj)
    prd_scr(iS) = (sum(cell2mat(prd_dta(iS,:))) / numel(cell2mat(prd_dta(iS,:)))) * 100;
end

%% Save
save([ dta_dir '/' 'demographics.mat' ],'cmb_dta', 'cmb_sbj', 'cmb_col');

save([ dta_dir '/' 'prediction.mat' ],'prd_scr', 'prd_dta_sbj');

