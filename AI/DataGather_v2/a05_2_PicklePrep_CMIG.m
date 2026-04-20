ejk_chk_dir([ new_dta_dir '/' '05_PicklePreps' '/'])

load([ new_dta_dir '/' '04_labels' '/' typ_nme{1} '/' typ_nme{1} '_' 'group.mat'], 'grp')

fcfg = [];
fcfg.dta_loc = [ new_dta_dir '/' '03_Datasets' '/' '00_' typ_nme{1} '_' 'dataset' '_' dte_str '.csv'];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%% SETUP
inc_vec = { ...
    [ grp.con_var_col.HC ; ... % INCLUDE ALL HC &
      grp.con_var_col.EPD ] }; ... % & INCLUDE ALL EPD

%% SUBSELECT
sub_sel_ind = inc_vec{1};
for i = 2:numel(inc_vec)
    sub_sel_ind = intersect(sub_sel_ind, inc_vec{i});
end


%% SETUP LABELS
% TLE
tle_vec = { ...
    grp.con_var_col.EPD ... % EPD
    grp.loc_var_col.temporal ... % TLE
    [ grp.lat_var_col.left ; ... % LTLE &
      grp.lat_var_col.right ] ...% & RTLE
      sub_sel_ind }; ... 
tle_ind = tle_vec{1};
for i = 2:numel(tle_vec)
    tle_ind = intersect(tle_ind, tle_vec{i});
end

% HC
con_vec = { ...
    grp.con_var_col.HC ... % HC
      sub_sel_ind }; ... 
con_ind = con_vec{1};
for i = 2:numel(con_vec)
    con_ind = intersect(con_ind, con_vec{i});
end

% L-TLE
lft_tle_vec = { ...
    grp.con_var_col.EPD ... % EPD
    grp.loc_var_col.temporal ... % TLE
    grp.lat_var_col.left ... % LTLE 
      sub_sel_ind }; ... 
lft_tle_ind = lft_tle_vec{1};
for i = 2:numel(lft_tle_vec)
    lft_tle_ind = intersect(lft_tle_ind, lft_tle_vec{i});
end

% R-TLE
rgh_tle_vec = { ...
    grp.con_var_col.EPD ... % EPD
    grp.loc_var_col.temporal ... % TLE
    grp.lat_var_col.right ... % RTLE
      sub_sel_ind }; ... 
rgh_tle_ind = rgh_tle_vec{1};
for i = 2:numel(rgh_tle_vec)
    rgh_tle_ind = intersect(rgh_tle_ind, rgh_tle_vec{i});
end

%% ADD LABELS
tot_col(end+1) = { 'TLE_versus_HC' };
tot_dta(tle_ind,end+1) = { 'TLE' };
tot_dta(con_ind,end) = { 'HC' };

tot_col(end+1) = { 'Left_TLE_versus_Right_TLE' };
tot_dta(lft_tle_ind,end+1) = { 'Left_TLE' };
tot_dta(rgh_tle_ind,end) = { 'Right_TLE' };

tot_dta = tot_dta([tle_ind ; con_ind],:);
tot_sbj = tot_sbj([tle_ind ; con_ind],:);

%% SAVE
out_dta = [ {'sbj_nme'} tot_col ; tot_sbj tot_dta];
cell2csv([ new_dta_dir '/' '05_PicklePreps' '/' 'CMIG_structure_function.csv'], out_dta)


