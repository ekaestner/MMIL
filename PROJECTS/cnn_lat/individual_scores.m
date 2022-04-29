clear; clc;

dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/performance_final';
dem_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/';

out_plt_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/scratch/individual';

%% Load Data
fcfg = [];
fcfg.dta_loc = [ dem_loc '/' 'demographic_table.csv'];
[ dem_dta, dem_dta_sbj, dem_dta_col ] = ejk_dta_frm(fcfg);

dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')) = strcat('Rush',cellfun(@num2str,dem_dta_sbj(strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Site')),'Rush')),'uni',0));

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

% MTS Status
grp.mts.L.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'yes') ),grp.side.L);
grp.mts.L.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'no') ),grp.side.L);

grp.mts.R.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'yes') ),grp.side.R);
grp.mts.R.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'MTS')),'no') ),grp.side.R);

% Surgical Status
grp.surgery.L.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'yes') ),grp.side.L);
grp.surgery.L.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'no') ),grp.side.L);

grp.surgery.R.yes  = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'yes') ),grp.side.R);
grp.surgery.R.no   = intersect(find( strcmpi(dem_dta(:,strcmpi(dem_dta_col,'Surgery')),'no') ),grp.side.R);

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
    end
end

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
fcfg.out_nme = 'p0_Legends.png';

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
fcfg.out_nme = 'p0_Number.png';

ejk_scatter(fcfg)

[~, pl0_stt_num] = ttest2(scr_num(grp.side.L),scr_num(grp.side.R));

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
fcfg.out_nme = 'p1_Side.png';

ejk_scatter(fcfg)

[~, pl1_stt_sde] = ttest2(scr_hld(grp.side.L),scr_hld(grp.side.R));

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
fcfg.out_nme = 'p2_Side_by_Number.png';

ejk_scatter(fcfg)

[~, pl2_stt_scr_cor_num_lft] = corrcoef(scr_num(grp.side.L),scr_hld(grp.side.L),'Rows','complete');
[~, pl2_stt_scr_cor_num_rgh] = corrcoef(scr_num(grp.side.R),scr_hld(grp.side.R),'Rows','complete');

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
fcfg.out_nme = 'p3_Side_by_automatedQC.png';

ejk_scatter(fcfg)

[~, pl3_stt_scr_cor_qal_lft] = corrcoef(qal_dta(grp.side.L),scr_hld(grp.side.L),'Rows','complete');
[~, pl3_stt_scr_cor_qal_rgh] = corrcoef(qal_dta(grp.side.R),scr_hld(grp.side.R),'Rows','complete');

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
fcfg.out_nme = 'p4_Side_by_hippocampus_left.png';

ejk_scatter(fcfg)

[~, pl4_stt_scr_cor_lft_hip_lft] = corrcoef(hip_dta(grp.side.L,1),scr_hld(grp.side.L),'Rows','complete');
[~, pl4_stt_scr_cor_lft_hip_rgh] = corrcoef(hip_dta(grp.side.R,1),scr_hld(grp.side.R),'Rows','complete');

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
fcfg.out_nme = 'p4_Side_by_hippocampus_right.png';

ejk_scatter(fcfg)

[~, pl4_stt_scr_cor_rgh_hip_lft] = corrcoef(hip_dta(grp.side.L,2),scr_hld(grp.side.L),'Rows','complete');
[~, pl4_stt_scr_cor_rgh_hip_rgh] = corrcoef(hip_dta(grp.side.R,2),scr_hld(grp.side.R),'Rows','complete');

% Right
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
fcfg.out_nme = 'p4_Side_by_hippocampus_LI.png';

ejk_scatter(fcfg)

[~, pl4_stt_scr_cor_lat_hip_lft] = corrcoef(hip_dta(grp.side.L,3),scr_hld(grp.side.L),'Rows','complete');
[~, pl4_stt_scr_cor_lat_hip_rgh] = corrcoef(hip_dta(grp.side.R,3),scr_hld(grp.side.R),'Rows','complete');

%% Plot 5: BarScatter Score -by- SurgicalStatus 
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
fcfg.out_nme = 'p5_SurgicalStatus.png';

ejk_scatter(fcfg)

[~, pl5_stt_srg_lft] = ttest2(scr_hld(grp.surgery.L.yes),scr_hld(grp.surgery.L.no));
[~, pl5_stt_srg_rgh] = ttest2(scr_hld(grp.surgery.R.yes),scr_hld(grp.surgery.R.no));

%% Plot 6: BarScatter Score -by- MTSstatus
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
fcfg.out_nme = 'p6_MTSStatus.png';

ejk_scatter(fcfg)

[~, pl6_stt_mts_lft] = ttest2(scr_hld(grp.mts.L.yes),scr_hld(grp.mts.L.no));
[~, pl6_stt_mts_rgh] = ttest2(scr_hld(grp.mts.R.yes),scr_hld(grp.mts.R.no));

%% Plot 7: BarScatter Score -by- Site
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
fcfg.out_nme = 'p7_Site_L.png';

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
fcfg.out_nme = 'p7_Site_R.png';

ejk_scatter(fcfg)


pl7_stt_ste_all = anova1(scr_hld, dem_dta(:,strcmpi(dem_dta_col,'Site')),'off');
pl7_stt_ste_lft = anova1(scr_hld(grp.side.L), dem_dta(grp.side.L,strcmpi(dem_dta_col,'Site')),'off');
pl7_stt_ste_rgh = anova1(scr_hld(grp.side.R), dem_dta(grp.side.R,strcmpi(dem_dta_col,'Site')),'off');

%% Stat table out
stt_out = { 'pl0_stt_num' pl0_stt_num ; ...
            '' '' ; ...
            'pl1_stt_sde' pl1_stt_sde ; ...
            '' '' ; ...
            'pl2_stt_scr_cor_num_lft' pl2_stt_scr_cor_num_lft(1,2) ; ...
            'pl2_stt_scr_cor_num_rgh' pl2_stt_scr_cor_num_rgh(1,2); ...
            '' '' ; ...
            'pl3_stt_scr_cor_qal_lft' pl3_stt_scr_cor_qal_lft(1,2); ...
            'pl3_stt_scr_cor_qal_rgh' pl3_stt_scr_cor_qal_rgh(1,2); ...
            '' '' ; ...
            'pl4_stt_scr_cor_lft_hip_lft' pl4_stt_scr_cor_lft_hip_lft(1,2); ...
            'pl4_stt_scr_cor_lft_hip_rgh' pl4_stt_scr_cor_lft_hip_rgh(1,2); ...
            'pl4_stt_scr_cor_rgh_hip_lft' pl4_stt_scr_cor_rgh_hip_lft(1,2); ...
            'pl4_stt_scr_cor_rgh_hip_rgh' pl4_stt_scr_cor_rgh_hip_rgh(1,2); ...
            'pl4_stt_scr_cor_lat_hip_lft' pl4_stt_scr_cor_lat_hip_lft(1,2); ...
            'pl4_stt_scr_cor_lat_hip_rgh' pl4_stt_scr_cor_lat_hip_rgh(1,2); ...
            '' '' ; ...
            'pl5_stt_srg_lft' pl5_stt_srg_lft; ...
            'pl5_stt_srg_rgh' pl5_stt_srg_rgh ; ...
            '' '' ; ...
            'pl6_stt_mts_lft' pl6_stt_mts_lft ; ...
            'pl6_stt_mts_lft' pl6_stt_mts_rgh ; ...
            '' '' ; ...
            'pl7_stt_ste_all' pl7_stt_ste_all ; ...
            'pl7_stt_ste_lft' pl7_stt_ste_lft ; ...
            'pl7_stt_ste_rgh' pl7_stt_ste_rgh };
            
cell2csv([ out_plt_dir '/' 'stat_out.csv'],stt_out);






