clear; clc;

dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/performance_final';
dem_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/';

out_plt_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/scratch/individual';

%% Load Data
fcfg = [];
fcfg.dta_loc = [ dem_loc '/' 'demographic_table.csv'];
[ dem_dta, dem_dta_sbj, dem_dta_col ] = ejk_dta_frm(fcfg);

dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')) = strcat('rush',cellfun(@num2str,dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')),'uni',0));

cell2csv([ dem_loc '/' 'demographic_table_rushname.csv'], [ {'sbj_nme'} dem_dta_col ; dem_dta_sbj dem_dta ]);

%
lft_sbj_nme = mmil_readtext([ dta_loc '/' 'fn_left.csv' ]);
rgh_sbj_nme = mmil_readtext([ dta_loc '/' 'fn_right.csv' ]);

cnn_lft_sbj = mmil_readtext([ dta_loc '/' 'exp_6_scores_left_sbj.csv' ]);
    cnn_lft_sbj(cellfun(@ischar,cnn_lft_sbj)) = {NaN};
    cnn_lft_sbj = cell2mat(cnn_lft_sbj);
cnn_rgh_sbj = mmil_readtext([ dta_loc '/' 'exp_6_scores_right_sbj.csv' ]);
    cnn_rgh_sbj(cellfun(@ischar,cnn_rgh_sbj)) = {NaN};
    cnn_rgh_sbj = cell2mat(cnn_rgh_sbj);

ind_prf_sbj = [ lft_sbj_nme  ; rgh_sbj_nme ];
ind_prf     = [ cnn_lft_sbj ; cnn_rgh_sbj ];

% hippocampal volume
fcfg = [];
fcfg.dta_loc = [ dta_loc '/' 'hippocamapal_volume_ordered.csv'];
[ hip_dta, hip_dta_sbj, hip_dta_col ] = ejk_dta_frm(fcfg);
hip_dta = cell2mat(hip_dta(:,2:3));
hip_dta(:,3) = (hip_dta(:,1)-hip_dta(:,2)) ./ (hip_dta(:,1)+hip_dta(:,2));

% automatedQC
fcfg = [];
fcfg.dta_loc = [ dta_loc '/' 'automatedQC_ordered.csv'];
[ qal_dta, qal_dta_sbj, qal_dta_col ] = ejk_dta_frm(fcfg);
qal_dta = cell2mat(qal_dta);

%% Make grp
% Side
grp.side.L = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Side')),'L') );
grp.side.R = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Side')),'R') );

% Sex
grp.sex.T.M  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'M') );
grp.sex.T.F  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'F') );

grp.sex.L.M  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'M') ),grp.side.L);
grp.sex.L.F   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'F') ),grp.side.L);

grp.sex.R.M  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'M') ),grp.side.R);
grp.sex.R.F   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Sex')),'F') ),grp.side.R);

% Handedness
grp.hand.T.L  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'L') );
grp.hand.T.R  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'R') );

grp.hand.L.L  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'L') ),grp.side.L);
grp.hand.L.R   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'R') ),grp.side.L);

grp.hand.R.L  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'L') ),grp.side.R);
grp.hand.R.R   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Handedness')),'R') ),grp.side.R);

% MTS Status
grp.mts.T.yes  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'yes') );
grp.mts.T.no   = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'no') );

grp.mts.L.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'yes') ),grp.side.L);
grp.mts.L.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'no') ),grp.side.L);

grp.mts.R.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'yes') ),grp.side.R);
grp.mts.R.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'no') ),grp.side.R);

% Surgical Status
grp.surgery.T.yes  = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'yes') );
grp.surgery.T.no   = find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'no') );

grp.surgery.L.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'yes') ),grp.side.L);
grp.surgery.L.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'no') ),grp.side.L);

grp.surgery.R.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'yes') ),grp.side.R);
grp.surgery.R.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'no') ),grp.side.R);

% Surgical Status
eng_one = strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Engel')),'I');
eng_two = strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Engel')),'II') | ...
          strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Engel')),'III') | ...
          strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Engel')),'IV');

grp.engel.T.I        = find( eng_one );
grp.engel.T.IIplus   = find( eng_two );

grp.engel.L.I  = intersect( find( eng_one ),grp.side.L);
grp.engel.L.IIplus   = intersect( find( eng_two ),grp.side.L);

grp.engel.R.I  = intersect( find( eng_one ),grp.side.R);
grp.engel.R.IIplus   = intersect( find( eng_two ),grp.side.R);

% Site
grp.site.L.Bonn  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Bonn') ),grp.side.L);
grp.site.L.CCF   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'CCF') ),grp.side.L);
grp.site.L.Emory = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Emory') ),grp.side.L);
grp.site.L.MUSC  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'MUSC') ),grp.side.L);
grp.site.L.Rush  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush') ),grp.side.L);
grp.site.L.UCSD  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'UCSD') ),grp.side.L);
grp.site.L.UCSF  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'UCSF') ),grp.side.L);

grp.site.R.Bonn  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Bonn') ),grp.side.R);
grp.site.R.CCF   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'CCF') ),grp.side.R);
grp.site.R.Emory = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Emory') ),grp.side.R);
grp.site.R.MUSC  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'MUSC') ),grp.side.R);
grp.site.R.Rush  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush') ),grp.side.R);
grp.site.R.UCSD  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'UCSD') ),grp.side.R);
grp.site.R.UCSF  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'UCSF') ),grp.side.R);

%% Make scores
scr_hld = NaN(size(dem_dta_sbj,1),1);
scr_num = NaN(size(dem_dta_sbj,1),1);
for iS=1:size(dem_dta_sbj,1)
    sbj_ind = strcmpi(ind_prf_sbj(:,1),dem_dta_sbj{iS});
    if ~(sum(sbj_ind)==0)
        scr_hld(iS,1) = round((nansum(ind_prf(sbj_ind,:)) / sum(~isnan(ind_prf(sbj_ind,:))))*100);
        scr_num(iS,1) = sum(~isnan(ind_prf(sbj_ind,:)));
        out_ind_prf(iS,:) = ind_prf(sbj_ind,:);
    end    
end

cell2csv([out_plt_dir '/' 'individual_performance.csv'],[ dem_dta_sbj num2cell(out_ind_prf)])
cell2csv([out_plt_dir '/' 'individual_performance_summary.csv'],[ dem_dta_sbj num2cell(scr_hld)])


%% Plot 0: Legend
fcfg = [];

fcfg.xdt = { 1 2 };
fcfg.ydt = { 1 1 };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.xlb = { 'L' 'R' };
fcfg.xlm = [ 0.5 2.5 ];
fcfg.ylb = {''};
fcfg.ylm = [ 0.75 1.25 ];

fcfg.mkr_sze = [2000 2000];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p0_Legends';

ejk_scatter(fcfg)

%% Plot 0: Number of times
fcfg = [];

fcfg.xdt = { 1                   2 };
fcfg.ydt = { scr_num(grp.side.L) scr_num(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L' 'R' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 30 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p0_Number';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_num(grp.side.L), scr_num(grp.side.R));
    stt_tbl(1,:) = { 'Number of Times: L vs R' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };

%% Plot 1: L- & R-TLE scatter
fcfg = [];

fcfg.xdt = { 1                   2 };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L' 'R' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p1_Side';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2(scr_hld(grp.side.L),scr_hld(grp.side.R));
    stt_tbl(3,:) = { 'Correct: L vs R' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };

%% Plot 2: Scatter Score -by- Number
fcfg = [];

fcfg.xdt = { scr_num(grp.side.L) scr_num(grp.side.R) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Number Trials' };
fcfg.xlm = [ 0 30 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p2_Side_by_Number';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ]  = corrcoef(scr_num,scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(scr_num)), find(~isnan(scr_hld)) ));
    stt_tbl(5,:) = { 'Correct -BY- Number: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ]  = corrcoef(scr_num(grp.side.L),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(scr_num(grp.side.L))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(6,:) = { 'Correct -BY- Number: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(scr_num(grp.side.R),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(scr_num(grp.side.R))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(7,:) = { 'Correct -BY- Number: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
    
%% Plot 3: Scatter Score -by- AutomatedQC
fcfg = [];

fcfg.xdt = { qal_dta(grp.side.L) qal_dta(grp.side.R) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'AutomatedQC' };
fcfg.xlm = [ 70 100 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p3_Side_by_automatedQC';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(qal_dta,scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(qal_dta)), find(~isnan(scr_hld)) ));
    stt_tbl(9,:) = { 'Correct -BY- QC: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(qal_dta(grp.side.L),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(qal_dta(grp.side.L))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(10,:) = { 'Correct -BY- QC: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(qal_dta(grp.side.R),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(qal_dta(grp.side.R))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(11,:) = { 'Correct -BY- QC: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

%% Plot 4: Scatter Score -by- HippocampalVolume 
% Left
fcfg = [];

fcfg.xdt = { hip_dta(grp.side.L,1) hip_dta(grp.side.R,1) };
fcfg.ydt = { scr_hld(grp.side.L)   scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Left Hippocampal Volume' };
% fcfg.xlm = [ 70 100 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p4_Side_by_hippocampus_left';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(:,1),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(:,1))), find(~isnan(scr_hld)) ));
    stt_tbl(13,:) = { 'Correct -BY- L Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.L,1),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.L,1))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(14,:) = { 'Correct -BY- L Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.R,1),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.R,1))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(15,:) = { 'Correct -BY- L Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

% Right
fcfg = [];

fcfg.xdt = { hip_dta(grp.side.L,2) hip_dta(grp.side.R,2) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Right Hippocampal Volume' };
% fcfg.xlm = [ 70 100 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p4_Side_by_hippocampus_right';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(:,2),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(:,2))), find(~isnan(scr_hld)) ));
    stt_tbl(17,:) = { 'Correct -BY- R Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.L,2),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.L,2))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(18,:) = { 'Correct -BY- R Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.R,2),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.R,2))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(19,:) = { 'Correct -BY- R Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

% LI
fcfg = [];

fcfg.xdt = { hip_dta(grp.side.L,3) hip_dta(grp.side.R,3) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'LI Hippocampal Volume' };
fcfg.xlm = [ -1 1 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p4_Side_by_hippocampus_LI';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(:,3),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(:,3))), find(~isnan(scr_hld)) ));
    stt_tbl(21,:) = { 'Correct -BY- LI Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.L,3),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.L,2))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(22,:) = { 'Correct -BY- LI Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(hip_dta(grp.side.R,3),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(hip_dta(grp.side.R,2))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(23,:) = { 'Correct -BY- LI Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

%% Plot 5: BarScatter Score -by- Site
col_hld = distinguishable_colors(7);

% LEFT
fcfg = [];

fcfg.xdt = { 1                        2                       3                         4                        5                        6                         7 };
fcfg.ydt = { scr_hld(grp.site.L.Bonn) scr_hld(grp.site.L.CCF) scr_hld(grp.site.L.Emory) scr_hld(grp.site.L.MUSC) scr_hld(grp.site.L.Rush) scr_hld(grp.site.L.UCSD)  scr_hld(grp.site.L.UCSF) };

fcfg.fce_col     = { col_hld(1,:) col_hld(2,:) col_hld(3,:) col_hld(4,:) col_hld(5,:) col_hld(6,:) col_hld(7,:)};
fcfg.edg_col     = { [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]       };
fcfg.box_plt_col = { rgb('black') rgb('black') rgb('black') rgb('black') rgb('black') rgb('black') rgb('black')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = fieldnames(grp.site.L);
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p5_Site_L';

ejk_scatter(fcfg)

% RIGHT
fcfg = [];

fcfg.xdt = { 1                        2                       3                         4                        5                        6                         7 };
fcfg.ydt = { scr_hld(grp.site.R.Bonn) scr_hld(grp.site.R.CCF) scr_hld(grp.site.R.Emory) scr_hld(grp.site.R.MUSC) scr_hld(grp.site.R.Rush) scr_hld(grp.site.R.UCSD)  scr_hld(grp.site.R.UCSF) };

fcfg.fce_col     = { col_hld(1,:) col_hld(2,:) col_hld(3,:) col_hld(4,:) col_hld(5,:) col_hld(6,:) col_hld(7,:)};
fcfg.edg_col     = { [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]      [0 0 0]       };
fcfg.box_plt_col = { rgb('black') rgb('black') rgb('black') rgb('black') rgb('black') rgb('black') rgb('black')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = fieldnames(grp.site.R);
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p5_Site_R';

ejk_scatter(fcfg)


[ pvl_hld, tst_hld, stt_hld ] = anova1(scr_hld, dem_dta(:,strcmpi(dem_dta_col,'Site')),'off');
    stt_tbl(25,:) = { 'Correct: Site; all' ['F(' num2str(tst_hld{2,3}) ';' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[ pvl_hld, tst_hld, stt_hld ] = anova1(scr_hld(grp.side.L), dem_dta(grp.side.L,strcmpi(dem_dta_col,'Site')),'off');
    stt_tbl(26,:) = { 'Correct: Site; LTLE' ['F(' num2str(tst_hld{2,3}) ';' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[ pvl_hld, tst_hld, stt_hld ] = anova1(scr_hld(grp.side.R), dem_dta(grp.side.R,strcmpi(dem_dta_col,'Site')),'off');
    stt_tbl(27,:) = { 'Correct: Site; RTLE' ['F(' num2str(tst_hld{2,3}) ';' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
  
%% Plot 6: Scatter Score -by- Age     
fcfg = [];

fcfg.xdt = { cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Age'))) cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Age'))) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Age' };
fcfg.xlm = [ 0 80 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p6_Side_by_Age';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'Age'))),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'Age'))))), find(~isnan(scr_hld)) ));
    stt_tbl(29,:) = { 'Correct -BY- Age: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Age'))),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Age'))))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(30,:) = { 'Correct -BY- Age: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Age'))),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Age'))))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(31,:) = { 'Correct -BY- Age: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };  

%% Plot 7: BarScatter Score -by- Sex
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.sex.L.M) scr_hld(grp.sex.L.F) scr_hld(grp.sex.R.M) scr_hld(grp.sex.R.F) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('dark magenta') rgb('dark violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0]};
fcfg.box_plt_col = { rgb('black')         rgb('black')        rgb('black')        rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L Male' 'L Female' 'R Male' 'R Female' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p7_Sex';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.sex.T.M), scr_hld(grp.sex.T.F));
    stt_tbl(33,:) = { 'Score: TLE; M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.sex.L.M), scr_hld(grp.sex.L.F));
    stt_tbl(34,:) = { 'Score: LTLE; M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.sex.R.M), scr_hld(grp.sex.R.F));
    stt_tbl(35,:) = { 'Score: RTLE; M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };

%% Plot 8: BarScatter Score -by- Handedness
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.hand.L.L) scr_hld(grp.hand.L.R) scr_hld(grp.hand.R.L) scr_hld(grp.hand.R.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('dark magenta') rgb('dark violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0]};
fcfg.box_plt_col = { rgb('black')         rgb('black')        rgb('black')        rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L RightH' 'L LeftH' 'R RightH' 'R LeftH' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p8_Handedness';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.hand.T.L), scr_hld(grp.hand.T.R));
    stt_tbl(37,:) = { 'Score: TLE; Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.hand.L.L), scr_hld(grp.hand.L.R));
    stt_tbl(38,:) = { 'Score: LTLE; Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.hand.R.L), scr_hld(grp.hand.R.R));
    stt_tbl(39,:) = { 'Score: RTLE; Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
    
%% Plot 9: Scatter Score -by- Education     
fcfg = [];

fcfg.xdt = { cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Education'))) cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Education'))) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Education' };
fcfg.xlm = [ 0 25 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p9_Side_by_Education';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'Education'))),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'Education'))))), find(~isnan(scr_hld)) ));
    stt_tbl(41,:) = { 'Correct -BY- Education: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Education'))),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'Education'))))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(42,:) = { 'Correct -BY- Education: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Education'))),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'Education'))))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(43,:) = { 'Correct -BY- Education: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };      
    
%% Plot 10: Scatter Score -by- AgeOfOnset     
fcfg = [];

fcfg.xdt = { cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'AO'))) cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'AO'))) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'AgeOfOnset' };
fcfg.xlm = [ 0 70 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p10_Side_by_AgeOfOnset';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'AO'))),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'AO'))))), find(~isnan(scr_hld)) ));
    stt_tbl(45,:) = { 'Correct -BY- AgeOfOnset: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'AO'))),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'AO'))))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(46,:) = { 'Correct -BY- AgeOfOnset: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'AO'))),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'AO'))))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(47,:) = { 'Correct -BY- AgeOfOnset: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };     
    
%% Plot 11: Scatter Score -by- Duration     
fcfg = [];

fcfg.xdt = { cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'DURILL'))) cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'DURILL'))) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Duration' };
fcfg.xlm = [ 0 65 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p11_Side_by_Duration';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'DURILL'))),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'DURILL'))))), find(~isnan(scr_hld)) ));
    stt_tbl(45,:) = { 'Correct -BY- AgeOfOnset: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'DURILL'))),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'DURILL'))))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(46,:) = { 'Correct -BY- AgeOfOnset: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'DURILL'))),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'DURILL'))))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(47,:) = { 'Correct -BY- AgeOfOnset: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };    
    
%% Plot 12: BarScatter Score -by- MTSstatus
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.mts.L.yes) scr_hld(grp.mts.L.no) scr_hld(grp.mts.R.yes) scr_hld(grp.mts.R.no) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('dark magenta') rgb('dark violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0]};
fcfg.box_plt_col = { rgb('black')         rgb('black')        rgb('black')        rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L MTS' 'L Non-MTS' 'R MTS' 'R Non-MTS' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p12_MTSStatus';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.mts.T.yes), scr_hld(grp.mts.T.no));
    stt_tbl(49,:) = { 'Score: TLE; MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
    tot_eff_sze     = (nanmean(scr_hld(grp.mts.T.yes))-nanmean(scr_hld(grp.mts.T.no))) / nanstd(scr_hld([grp.mts.T.yes ; grp.mts.T.no]))
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.mts.L.yes), scr_hld(grp.mts.L.no));
    stt_tbl(50,:) = { 'Score: LTLE; MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.mts.R.yes), scr_hld(grp.mts.R.no));
    stt_tbl(51,:) = { 'Score: RTLE; MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };  
    
%% Plot 13: Scatter Score -by- ASMs     
fcfg = [];

fcfg.xdt = { cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'ASMs'))) cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'ASMs'))) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { rgb('dark magenta')  rgb('dark violet')  };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'ASMs' };
fcfg.xlm = [ 0 7 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p13_Side_by_ASMs';

ejk_scatter(fcfg)

[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'ASMs'))),scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(:,strcmpi(dem_dta_col,'ASMs'))))), find(~isnan(scr_hld)) ));
    stt_tbl(53,:) = { 'Correct -BY- ASMs: TLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'ASMs'))),scr_hld(grp.side.L),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.L,strcmpi(dem_dta_col,'ASMs'))))), find(~isnan(scr_hld(grp.side.L))) ));
    stt_tbl(54,:) = { 'Correct -BY- ASMs: LTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
[ rvl_hld, pvl_hld ] = corrcoef(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'ASMs'))),scr_hld(grp.side.R),'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(cell2mat(dem_dta(grp.side.R,strcmpi(dem_dta_col,'ASMs'))))), find(~isnan(scr_hld(grp.side.R))) ));
    stt_tbl(55,:) = { 'Correct -BY- ASMs: RTLE' ['r(' num2str(deg_fre-2) ')' '=' num2str(roundsd(rvl_hld(1,2),2)) ';' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };   
    
%% Plot 14: BarScatter Score -by- SurgicalStatus 
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.surgery.L.yes) scr_hld(grp.surgery.L.no) scr_hld(grp.surgery.R.yes) scr_hld(grp.surgery.R.no) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('dark magenta') rgb('dark violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0]};
fcfg.box_plt_col = { rgb('black')         rgb('black')        rgb('black')        rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L Surgery' 'L Non-Surgical' 'R Surgery' 'R Non-Surgical' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p14_SurgicalStatus';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.surgery.T.yes), scr_hld(grp.surgery.T.no));
    stt_tbl(57,:) = { 'Score: TLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.surgery.L.yes), scr_hld(grp.surgery.L.no));
    stt_tbl(58,:) = { 'Score: LTLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.surgery.R.yes), scr_hld(grp.surgery.R.no));
    stt_tbl(59,:) = { 'Score: RTLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };  

%% Plot 15: BarScatter Score -by- SurgicalOutcome 
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.engel.L.I) scr_hld(grp.engel.L.IIplus) scr_hld(grp.engel.R.I) scr_hld(grp.engel.R.IIplus) };

fcfg.fce_col     = { rgb('light magenta') rgb('light violet') rgb('dark magenta') rgb('dark violet') };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             [0 0 0]             [0 0 0]};
fcfg.box_plt_col = { rgb('black')         rgb('black')        rgb('black')        rgb('black') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L Engel I' 'L Engel II+' 'R Engel I' 'R Engel II+' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p15_SurgicalOutcome';

ejk_scatter(fcfg)

[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.engel.T.I), scr_hld(grp.engel.T.IIplus));
    stt_tbl(61,:) = { 'Score: TLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.engel.L.I), scr_hld(grp.engel.L.IIplus));
    stt_tbl(62,:) = { 'Score: LTLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] };
[~, pvl_hld, ~, stt_hld] = ttest2( scr_hld(grp.engel.R.I), scr_hld(grp.engel.R.IIplus));
    stt_tbl(63,:) = { 'Score: RTLE; Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ';' 'p=' num2str(roundsd(pvl_hld,2))] }; 
   
%% Resubmission 
fcfg = [];

fcfg.xdt = { 1 2 };
fcfg.ydt = { 1 1 };

fcfg.fce_col     = { [0/255 84/255 95/255] [163/255 2/255 52/255] };
fcfg.edg_col     = { [0 0 0]              [0 0 0]             };
fcfg.box_plt_col = { [0/255 84/255 95/255] [163/255 2/255 52/255] };

fcfg.xlb = { 'L' 'R' };
fcfg.xlm = [ 0.5 2.5 ];
fcfg.ylb = {''};
fcfg.ylm = [ 0.75 1.25 ];

fcfg.mkr_sze = [2000 2000];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p0_Legends_resubmission';

ejk_scatter(fcfg)

% p1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 1                   2 };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { [0/255 84/255 95/255] [163/255 2/255 52/255]  };
fcfg.edg_col     = { [0 0 0]               [0 0 0]             };
fcfg.box_plt_col = { [0/255 84/255 95/255] [163/255 2/255 52/255]   };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L' 'R' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');


fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p1_Side_resubmission';

ejk_scatter(fcfg)

% p4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { hip_dta(grp.side.L,3) hip_dta(grp.side.R,3) };
fcfg.ydt = { scr_hld(grp.side.L) scr_hld(grp.side.R) };

fcfg.fce_col     = { [0/255 84/255 95/255] [163/255 2/255 52/255]  };
fcfg.edg_col     = { [0 0 0]               [0 0 0]             };
fcfg.box_plt_col = { [0/255 84/255 95/255] [163/255 2/255 52/255]   };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'LI Hippocampal Volume' };
fcfg.xlm = [ -0.75 0.75 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.trd_lne = [1 1];

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p4_Side_by_hippocampus_LI_resubmission';

ejk_scatter(fcfg)

% p12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt = { 1                          2                         4                          5 };
fcfg.ydt = { scr_hld(grp.mts.L.yes) scr_hld(grp.mts.L.no) scr_hld(grp.mts.R.yes) scr_hld(grp.mts.R.no) };

fcfg.fce_col     = { [0/255 47/255 48/255] [86/255 152/255 163/255] [81/255 29/255 36/255] [206/255 128/255 128/255] };
fcfg.edg_col     = { [0 0 0]               [0 0 0]                  [0 0 0]                [0 0 0]};
fcfg.box_plt_col = { [0/255 47/255 48/255] [86/255 152/255 163/255] [81/255 29/255 36/255] [206/255 128/255 128/255] };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'L MTS' 'L Non-MTS' 'R MTS' 'R Non-MTS' };
fcfg.xlm = [ 0.5 7.5 ];
fcfg.ylb = {'Accuracy'};
fcfg.ylm = [ 0 100 ];

fcfg.mkr_sze = [20 20 20 20];
fcfg.aph_val = 0.45;

fcfg.hln = 0;
fcfg.hln_col = rgb('black');

fcfg.out_dir = out_plt_dir;
fcfg.out_nme = 'p12_MTSStatus_resubmission';

ejk_scatter(fcfg)    

% Resubmission numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_mts_chn = computeCohen_d(scr_hld(grp.mts.T.yes),scr_hld(grp.mts.T.no),'independent');
lft_mts_chn = computeCohen_d(scr_hld(grp.mts.L.yes),scr_hld(grp.mts.L.no),'independent');
rgh_mts_chn = computeCohen_d(scr_hld(grp.mts.R.yes),scr_hld(grp.mts.R.no),'independent');

%
rgh_zro_num = sum(scr_hld(grp.side.R)==0);
rgh_zro_pct = sum(scr_hld(grp.side.R)==0) / numel(scr_hld(grp.side.R));

%
men_tot_yes = nanmean( scr_hld(grp.mts.T.yes) );
med_tot_yes = nanmedian( scr_hld(grp.mts.T.yes) );
std_tot_yes = nanstd( scr_hld(grp.mts.T.yes) );

men_tot_no  = nanmean( scr_hld(grp.mts.T.no) );
med_tot_no  = nanmedian( scr_hld(grp.mts.T.no) );
std_tot_no = nanstd( scr_hld(grp.mts.T.no) );

men_lft_yes = nanmean( scr_hld(grp.mts.L.yes) );
med_lft_yes = nanmedian( scr_hld(grp.mts.L.yes) );
std_lft_yes = nanstd( scr_hld(grp.mts.L.yes) );

men_lft_no = nanmean( scr_hld(grp.mts.L.no) );
med_lft_no = nanmedian( scr_hld(grp.mts.L.no) );
std_lft_no = nanstd( scr_hld(grp.mts.L.no) );

men_rgh_yes = nanmean( scr_hld(grp.mts.R.yes) );
med_rgh_yes = nanmedian( scr_hld(grp.mts.R.yes) );
std_rgh_yes = nanstd( scr_hld(grp.mts.R.yes) );

men_rgh_no = nanmean( scr_hld(grp.mts.R.no) );
med_rgh_no = nanmedian( scr_hld(grp.mts.R.no) );
std_rgh_no = nanstd( scr_hld(grp.mts.R.no) );

% 
out_mts = [ {'rgh_zro_num' ' = ' num2str(rgh_zro_num)} ; ...
            {'rgh_zro_pct' ' = ' num2str(rgh_zro_pct)} ; ...
            {'men_tot_yes' ' = ' num2str(men_tot_yes)} ; ...
            {'med_tot_yes' ' = ' num2str(med_tot_yes)} ; ...
            {'std_tot_yes' ' = ' num2str(std_tot_yes)} ; ...
            {'men_tot_no ' ' = ' num2str(men_tot_no)} ; ...
            {'med_tot_no ' ' = ' num2str(med_tot_no)} ; ...
            {'std_tot_no ' ' = ' num2str(std_tot_no)} ; ...
            {'men_lft_yes' ' = ' num2str(men_lft_yes)} ; ...
            {'med_lft_yes' ' = ' num2str(med_lft_yes)} ; ...
            {'std_lft_yes' ' = ' num2str(std_lft_yes)} ; ...
            {'men_lft_no ' ' = ' num2str(men_lft_no)} ; ...
            {'med_lft_no ' ' = ' num2str(med_lft_no)} ; ...
            {'std_lft_no ' ' = ' num2str(std_lft_no)} ; ...
            {'men_rgh_yes' ' = ' num2str(men_rgh_yes)} ; ...
            {'med_rgh_yes' ' = ' num2str(med_rgh_yes)} ; ...
            {'std_rgh_yes' ' = ' num2str(std_rgh_yes)} ; ...
            {'men_rgh_no ' ' = ' num2str(men_rgh_no)} ; ...
            {'med_rgh_no ' ' = ' num2str(med_rgh_no)} ; ...
            {'std_rgh_no ' ' = ' num2str(std_rgh_no)} ; ...
            {'tot_mts_chn' ' = ' num2str(tot_mts_chn)} ; ...
            {'lft_mts_chn' ' = ' num2str(lft_mts_chn)} ; ...
            {'rgh_mts_chn' ' = ' num2str(rgh_mts_chn)} ];
cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/submission/2022_09_Neurology/Reviews/Figures/NewPieces' '/' 'out_mts.csv'],out_mts); 
 
%
cut_off = 50;
cnn_acc     = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_test.csv' ]))*100;
    cnn_acc = cnn_acc(:,cut_off);
log_acc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_scores_log.csv' ]))*100;
shf_acc = cell2mat(mmil_readtext([ dta_loc '/' 'exp_6_score_test_rand.csv' ]))*100;
    shf_acc = shf_acc(:,cut_off);
cell2csv( [ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/submission/2022_09_Neurology/Reviews/Figures/NewPieces' '/' 'scores.csv'],num2cell([ cnn_acc log_acc shf_acc ]) )
    
%% Make table of outputs
tot_tbl = { ''                           'TLE'         'L-TLE'        'R-TLE' ; ...
  'Correct -BY- Number'         stt_tbl{5,2}  stt_tbl{6,2}  stt_tbl{7,2} ; ...
  'Correct -BY- QC'             stt_tbl{9,2}  stt_tbl{10,2} stt_tbl{11,2} ; ...
  'Correct -BY- L-Hippocampus'  stt_tbl{13,2} stt_tbl{14,2} stt_tbl{15,2} ; ...
  'Correct -BY- R-Hippocampus'  stt_tbl{17,2} stt_tbl{18,2} stt_tbl{19,2} ; ...
  'Correct -BY- LI-Hippocampus' stt_tbl{21,2} stt_tbl{22,2} stt_tbl{23,2} ; ...
  'Site'                        stt_tbl{25,2} stt_tbl{26,2} stt_tbl{27,2} ; ...
  'Correct -BY- Age'            stt_tbl{29,2} stt_tbl{30,2} stt_tbl{31,2} ; ...
  'Sex'                         stt_tbl{33,2} stt_tbl{34,2} stt_tbl{35,2} ; ...
  'Handedness'                  stt_tbl{37,2} stt_tbl{38,2} stt_tbl{39,2} ; ...
  'Correct -BY- Education'      stt_tbl{41,2} stt_tbl{42,2} stt_tbl{43,2} ; ...
  'Correct -BY- AgeOfOnset'     stt_tbl{45,2} stt_tbl{46,2} stt_tbl{47,2} ; ...
  'MTS'                         stt_tbl{49,2} stt_tbl{50,2} stt_tbl{51,2} ; ...
  'Correct -BY- ASMs'           stt_tbl{53,2} stt_tbl{54,2} stt_tbl{55,2} ; ...
  'Surgery'                     stt_tbl{57,2} stt_tbl{58,2} stt_tbl{59,2} ; ...
  'Engel'                       stt_tbl{61,2} stt_tbl{62,2} stt_tbl{63,2} };


cell2csv([ out_plt_dir '/' 'stat_table.csv'],stt_tbl);  
cell2csv([ out_plt_dir '/' 'stat_table_organized.csv'],tot_tbl);  

%% Misc out
msc_out{1,1} = sprintf('Avereage number of groups: %.2f (%.2f; %i - %i)', nanmean(scr_num), nanstd(scr_num), min(scr_num), max(scr_num));
msc_out{2,1} = sprintf('%.2f above 90%',(sum(scr_hld>=90)/numel(scr_hld))*100 );
msc_out{3,1} = sprintf('%.2f above 80%',(sum(scr_hld>=80)/numel(scr_hld))*100 );
msc_out{4,1} = sprintf('%.2f below 20%',(sum(scr_hld<=20)/numel(scr_hld))*100 );
msc_out{5,1} = sprintf('%.2f below 40%',(sum(scr_hld<=40)/numel(scr_hld))*100 );

cell2csv([ out_plt_dir '/' 'stat_table_misc.csv'],stt_tbl); 

%% Stat table out
% stt_out = { 'pl0_stt_num' pl0_stt_num ; ...
%             '' '' ; ...
%             'pl1_stt_sde' pl1_stt_sde ; ...
%             '' '' ; ...
%             'pl2_stt_scr_cor_num_lft' pl2_stt_scr_cor_num_lft(1,2) ; ...
%             'pl2_stt_scr_cor_num_rgh' pl2_stt_scr_cor_num_rgh(1,2); ...
%             '' '' ; ...
%             'pl3_stt_scr_cor_qal_lft' pl3_stt_scr_cor_qal_lft(1,2); ...
%             'pl3_stt_scr_cor_qal_rgh' pl3_stt_scr_cor_qal_rgh(1,2); ...
%             '' '' ; ...
%             'pl4_stt_scr_cor_lft_hip_lft' pl4_stt_scr_cor_lft_hip_lft(1,2); ...
%             'pl4_stt_scr_cor_lft_hip_rgh' pl4_stt_scr_cor_lft_hip_rgh(1,2); ...
%             'pl4_stt_scr_cor_rgh_hip_lft' pl4_stt_scr_cor_rgh_hip_lft(1,2); ...
%             'pl4_stt_scr_cor_rgh_hip_rgh' pl4_stt_scr_cor_rgh_hip_rgh(1,2); ...
%             'pl4_stt_scr_cor_lat_hip_lft' pl4_stt_scr_cor_lat_hip_lft(1,2); ...
%             'pl4_stt_scr_cor_lat_hip_rgh' pl4_stt_scr_cor_lat_hip_rgh(1,2); ...
%             '' '' ; ...
%             'pl5_stt_srg_lft' pl5_stt_srg_lft; ...
%             'pl5_stt_srg_rgh' pl5_stt_srg_rgh ; ...
%             '' '' ; ...
%             'pl6_stt_mts_lft' pl6_stt_mts_lft ; ...
%             'pl6_stt_mts_lft' pl6_stt_mts_rgh ; ...
%             '' '' ; ...
%             'pl7_stt_ste_all' pl7_stt_ste_all ; ...
%             'pl7_stt_ste_lft' pl7_stt_ste_lft ; ...
%             'pl7_stt_ste_rgh' pl7_stt_ste_rgh };
%             
% cell2csv([ out_plt_dir '/' 'stat_out.csv'],stt_out);






