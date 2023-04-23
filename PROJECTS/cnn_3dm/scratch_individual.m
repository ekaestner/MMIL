clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm/';

dem_dir = [ prj_dir '/' 'Data/demographics' ]; % 2022_12_12

prd_dir = [ prj_dir '/' 'Data/2023_02_17/' ];

out_dir = [ prj_dir '/' 'Data/2023_02_17/individual' ];

cor_nme = { 'Age' 'Age at Onset' 'Duration' 'Education Years' };
grp_nme = { 'Gender' 'Group' 'Handedness' 'Side of Epilepsy' 'Site' };

%% Load Data
fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_TLE.csv' ];
[ tle_dem, tle_dem_sbj, tle_dem_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dem_dir '/' 'Data_forAnees_2022_09_07_HC.csv' ];
[ con_dem, con_dem_sbj, con_dem_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ prd_dir '/' 'predictions.csv' ];
[ prd_dta, prd_dta_sbj, prd_dta_col ] = ejk_dta_frm(fcfg);

%% Prediction: Calculate Performances
prd_scr = nan(numel(prd_dta_sbj),1);
for iS = 1:numel(prd_dta_sbj)
    prd_scr(iS) = (sum(cell2mat(prd_dta(iS,:))) / numel(cell2mat(prd_dta(iS,:)))) * 100;
end

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

%% Correlations: Age, Education, Age of Onset, Duration
% Split by patient population
for iCR = 1:numel(cor_nme)
    
    % Collate
    cor_col = strcmpi(cmb_col,cor_nme{iCR});
    cor_row = ~cellfun( @isempty, cmb_dta(:,cor_col));
    
    cor_dta     = cell2mat(cmb_dta(cor_row,cor_col));
    cor_sbj     = cmb_sbj(cor_row,1);
    prd_scr_use = prd_scr(cor_row,1);
    
    [ rho, poh ] = corr(cor_dta,prd_scr_use,'Type','Spearman','rows','complete');
    [ rvl, pvl ] = corr(cor_dta,prd_scr_use,'Type','Pearson','rows','complete');
    
    % Plot
    fcfg = [];
    
    fcfg.xdt = { cor_dta };
    fcfg.ydt = { prd_scr_use };
    
    fcfg.fce_col     = { rgb('light orange') };
    fcfg.edg_col     = { [0 0 0]              };
    fcfg.box_plt_col = { rgb('dark orange')  };
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));
    fcfg.xlb = { cor_nme{iCR} };
    fcfg.ylb = {'Accuracy'};
    fcfg.ylm = [ 0 100 ];
    
    fcfg.ttl = [ cor_nme{iCR} ' ' 'pvl:' ' ' num2str(roundsd(pvl,2)) ' ' num2str(roundsd(poh,2)) ];

    fcfg.trd_lne = [1 1];
    
    fcfg.out_dir = out_dir;
    fcfg.out_nme = [ 'corr' '_' 'p' num2str(iCR) '_' mmil_spec_char(cor_nme{iCR},{' '},{'_'}) ];
    
    ejk_scatter(fcfg)
    
end

%% Group comparisons: Gender, Side of Epilepsy, Handedness, 
for iAN = 1:numel(grp_nme)
    
    % Collate
    anv_col = strcmpi(cmb_col,grp_nme{iAN});
    anv_row = ~strcmpi( cmb_dta(:,anv_col),'');
    
    anv_dta     = cmb_dta(anv_row,anv_col);
    cor_sbj     = cmb_sbj(anv_row,1);
    prd_scr_use = prd_scr(anv_row,1);
    
    [ pvl_hld, tst_hld, stt_hld ] = anova1(prd_scr_use, anv_dta,'off');
    
    % Plot
    lvl_nme = unique(anv_dta);
    
    fcfg = [];
    
    for iL = 1:numel(lvl_nme)
        fcfg.xdt{iL} = iL;
        fcfg.ydt{iL} = prd_scr_use(strcmpi(anv_dta,lvl_nme{iL}));
        fcfg.xlb{iL} = lvl_nme{iL};
    end
    
    fcfg.fce_col     = repmat( {[0 0 0]}, 1, numel(fcfg.xdt) );
    fcfg.edg_col     = repmat( {[0 0 0]}, 1, numel(fcfg.xdt) );
    fcfg.box_plt_col = repmat( {[0.6 0.6 0.6]}, 1, numel(fcfg.xdt) );
    
    fcfg.box_plt = ones(1,numel(fcfg.xdt));

    fcfg.ylb = {'Accuracy'};
    fcfg.ylm = [ 0 100 ];
    
    fcfg.mkr_sze = repmat(15,1,numel(fcfg.xdt));
    fcfg.aph_val = 0.45;
    
    fcfg.ttl = [ grp_nme{iAN} ' ' 'pvl:' ' ' num2str(roundsd(pvl_hld,2)) ];
    
    fcfg.out_dir = out_dir;
    fcfg.out_nme = [ 'anova' '_' 'p' num2str(iAN) '_' mmil_spec_char(grp_nme{iAN},{' '},{'_'}) ];
    
    ejk_scatter(fcfg)
    
end


%% Check sites
iAN = 5;
anv_col = strcmpi(cmb_col,grp_nme{iAN});
anv_row = ~strcmpi( cmb_dta(:,anv_col),'');

anv_dta     = cmb_dta(anv_row,anv_col);
cor_sbj     = cmb_sbj(anv_row,1);
prd_scr_use = prd_scr(anv_row,1);

[ pvl_hld, tst_hld, stt_hld ] = anova1(prd_scr_use, anv_dta,'off');
[c,m,h,gnames] = multcompare(stt_hld);




