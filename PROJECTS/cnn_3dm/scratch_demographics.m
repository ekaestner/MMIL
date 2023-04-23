clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm/';

dem_dir = [ prj_dir '/' 'Data/demographics' ]; % 2022_12_12

out_dir = [ prj_dir '/' 'Data/demographics/tables' ];

%% Load
fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_TLE.csv' ];
[ tle_dem, tle_dem_sbj, tle_dem_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_HC.csv' ];
[ con_dem, con_dem_sbj, con_dem_col ] = ejk_dta_frm(fcfg);

%% Calculate missing data
% TLE
ste_nme     = unique(tle_dem(:,strcmpi(tle_dem_col,'Site')));
tle_mss_out = cell(numel(ste_nme)+1,numel(tle_dem_col));

for iC = 1:numel(tle_dem_col)
    if iC ~= 1
        tle_mss_out{1,iC} = [ num2str(sum(~cellfun(@isempty,tle_dem(:,iC)))) ' (' num2str(round((sum(~cellfun(@isempty,tle_dem(:,iC))) / tle_mss_out{1,1})*100)) '%)' ];
    else
        tle_mss_out{1,iC} = sum(~cellfun(@isempty,tle_dem(:,iC)));
    end
    
    for iST = 1:numel(ste_nme)
        sbj_ind = strcmpi( tle_dem(:,strcmpi(tle_dem_col,'Site')), ste_nme{iST} );
        if iC ~= 1
            tle_mss_out{iST+1,iC} = [ num2str(sum(~cellfun(@isempty,tle_dem(sbj_ind,iC)))) ' (' num2str(round((sum(~cellfun(@isempty,tle_dem(sbj_ind,iC))) / tle_mss_out{iST+1,1})*100)) '%)' ];
        else
            tle_mss_out{iST+1,iC} = sum(~cellfun(@isempty,tle_dem(sbj_ind,iC)));
        end
    end
end

cell2csv( [out_dir '/' 'tle_missing.csv'], [ [tle_dem_col(1) 'N' tle_dem_col(2:end)] ;  [ 'Total' ; ste_nme] tle_mss_out ])

% HC
ste_nme     = unique(con_dem(:,strcmpi(con_dem_col,'Site')));
con_mss_out = cell(numel(ste_nme)+1,numel(con_dem_col));

for iC = 1:numel(con_dem_col)
    if iC ~= 1
        con_mss_out{1,iC} = [ num2str(sum(~cellfun(@isempty,con_dem(:,iC)))) ' (' num2str(round((sum(~cellfun(@isempty,con_dem(:,iC))) / con_mss_out{1,1})*100)) '%)' ];
    else
        con_mss_out{1,iC} = sum(~cellfun(@isempty,con_dem(:,iC)));
    end
    
    for iST = 1:numel(ste_nme)
        sbj_ind = strcmpi( con_dem(:,strcmpi(con_dem_col,'Site')), ste_nme{iST} );
        if iC ~= 1
            con_mss_out{iST+1,iC} = [ num2str(sum(~cellfun(@isempty,con_dem(sbj_ind,iC)))) ' (' num2str(round((sum(~cellfun(@isempty,con_dem(sbj_ind,iC))) / con_mss_out{iST+1,1})*100)) '%)' ];
        else
            con_mss_out{iST+1,iC} = sum(~cellfun(@isempty,con_dem(sbj_ind,iC)));
        end
    end
end

cell2csv( [out_dir '/' 'con_missing.csv'], [ [con_dem_col(1) 'N' con_dem_col(2:end)] ;  [ 'Total' ; ste_nme] con_mss_out ])

%% Recode
