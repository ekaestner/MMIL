%% Look at DATA
use_sbj = find(~isnan(cell2mat(lod_dta(2:end,32))));

% LHS 313 %%%
lhs_dta = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313_epd006switch_epd096add.mat'); % surf_wmparc_fa_lhs_sm313_epd006switch_epd096add

lhs_roi = lod_dta(2:end,[1 31]);

[col_loc,albl,actbl]=fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/lh.aparc.annot');
ent_roi = find(col_loc==7);

for iS = 1:size(lhs_roi,1)
   dta_out_lhs{iS,1} = lhs_roi{iS,1};
   dta_out_lhs{iS,2} = lhs_roi{iS,2};
   
   if any(iS==use_sbj) && ~isnan(dta_out_lhs{iS,2}) && ~strcmpi(rhs_roi{iS,1},'epd_ucsf048') && ~strcmpi(rhs_roi{iS,1},'epd096') && ~strcmpi(rhs_roi{iS,1},'epd006') && ~strcmpi(rhs_roi{iS,1},'epd041')
       dta_out_lhs{iS,3} = mean(lhs_dta.srf_dta( find(strcmpi(lhs_dta.srf_dta_sbj,dta_out_lhs{iS,1})), ent_roi));
   else
       dta_out_lhs{iS,3} = nan;
   end
   
end

figure()
subplot(2,2,1)
scatter( cell2mat(dta_out_lhs(:,2)), cell2mat(dta_out_lhs(:,3)) )
title('LHS Enotrhinal 22mm^2')

% RHS 313 %%%
rhs_dta = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm313_epd006switch_epd096add.mat'); % surf_wmparc_fa_rhs_sm313_epd006switch_epd096add

rhs_roi = lod_dta(2:end,[1 32]);

[col_loc,albl,actbl]=fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/rh.aparc.annot');
ent_roi = find(col_loc==7);

for iS = 1:size(rhs_roi,1)
   dta_out_rhs{iS,1} = rhs_roi{iS,1};
   dta_out_rhs{iS,2} = rhs_roi{iS,2};
   
   if ~isnan(dta_out_rhs{iS,2}) && ~strcmpi(rhs_roi{iS,1},'epd_ucsf048') && ~strcmpi(rhs_roi{iS,1},'epd096') && ~strcmpi(rhs_roi{iS,1},'epd006') && ~strcmpi(rhs_roi{iS,1},'epd041')
       dta_out_rhs{iS,3} = mean(rhs_dta.srf_dta( find(strcmpi(rhs_dta.srf_dta_sbj,dta_out_rhs{iS,1})), ent_roi));
   else
       dta_out_rhs{iS,3} = nan;
   end
   
end

subplot(2,2,2)
scatter( cell2mat(dta_out_rhs(:,2)), cell2mat(dta_out_rhs(:,3)) )
title('RHS Enotrhinal 22mm^2')

clear dta_out_lhs dta_out_rhs

% LHS 176 %%%
lhs_dta = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm176_epd006switch.mat'); % surf_wmparc_fa_lhs_sm313_epd006switch_epd096add

lhs_roi = lod_dta(2:end,[1 31]);

[col_loc,albl,actbl]=fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/lh.aparc.annot');
ent_roi = find(col_loc==7);

for iS = 1:size(lhs_roi,1)
   dta_out_lhs{iS,1} = lhs_roi{iS,1};
   dta_out_lhs{iS,2} = lhs_roi{iS,2};
   
   if any(iS==use_sbj) && ~isnan(dta_out_lhs{iS,2}) && ~strcmpi(rhs_roi{iS,1},'epd_ucsf048') && ~strcmpi(rhs_roi{iS,1},'epd096') && ~strcmpi(rhs_roi{iS,1},'epd006') && ~strcmpi(rhs_roi{iS,1},'epd041')           
       dta_out_lhs{iS,3} = mean(lhs_dta.srf_dta( find(strcmpi(lhs_dta.srf_dta_sbj,dta_out_lhs{iS,1})), ent_roi));
   else
       dta_out_lhs{iS,3} = nan;
   end
   
end

subplot(2,2,3)
scatter( cell2mat(dta_out_lhs(:,2)), cell2mat(dta_out_lhs(:,3)) )
title('LHS Enotrhinal 16mm^2')

% RHS 176 %%%
rhs_dta = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_rhs_sm176_epd006switch.mat'); % surf_wmparc_fa_rhs_sm313_epd006switch_epd096add

rhs_roi = lod_dta(2:end,[1 32]);

[col_loc,albl,actbl]=fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/rh.aparc.annot');
ent_roi = find(col_loc==7);

for iS = 1:size(rhs_roi,1)
   dta_out_rhs{iS,1} = rhs_roi{iS,1};
   dta_out_rhs{iS,2} = rhs_roi{iS,2};
   
   if ~isnan(dta_out_rhs{iS,2}) && ~strcmpi(rhs_roi{iS,1},'epd_ucsf048') && ~strcmpi(rhs_roi{iS,1},'epd096') && ~strcmpi(rhs_roi{iS,1},'epd006') && ~strcmpi(rhs_roi{iS,1},'epd041')       
       dta_out_rhs{iS,3} = mean(rhs_dta.srf_dta( find(strcmpi(rhs_dta.srf_dta_sbj,dta_out_rhs{iS,1})), ent_roi));
   else
       dta_out_rhs{iS,3} = nan;
   end
   
end

subplot(2,2,4)
scatter( cell2mat(dta_out_rhs(:,2)), cell2mat(dta_out_rhs(:,3)) )
title('RHS Enotrhinal 16mm^2')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at RIGHT HEMISPHERE
rvl_lhs = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February/VPA2_n21_Covariate_none_ATL_only/rvalues_lhs_dep_var_left.mat'); %%% Find LHS Vertex
[ ~, lhs_ent_ind ] = max(rvl_lhs.rvalues);

rvl_rhs = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February/VPA2_n21_Covariate_none_ATL_only/rvalues_rhs_dep_var_left.mat'); %%% Find RHS Vertex
[ ~, rhs_ent_ind ] = max(rvl_rhs.rvalues);

fcfg = []; %%% make plot

fcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test';
fcfg.out_pre_fix = 'rhs_entorhinal_cortex';

fcfg.vtx_cor = { [ 5483          ]        [ 7034               39922               62522              89649] }; % 
fcfg.vtx_col = { { rgb('bright red') }    { rgb('bright green') rgb('bright purple') rgb('bright red') rgb('cyan') } };

fcfg.hme_wrk = 1;

ejk_highlight_vertex(fcfg)

%%
ent_ent_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February/EntorhinalRight_n21_Covariate_none_ATL_only_Raw/exact_data_2.csv');
figure()
subplot(2,2,1)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_ent_chk(2:end,5))); hold on;
title('ENTORHINAL LHS -BY- SURF FA')

ent_lm2_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February/LM2_n21_Covariate_none_ATL_only_Raw/exact_data_2.csv');
subplot(2,2,2)
scatter( cell2mat(ent_lm2_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,5))); hold on;
title('LM2 -BY- SURF FA')

ent_vpa_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February/VPA2_n21_Covariate_none_ATL_only/exact_data_2.csv');
subplot(2,2,3)
scatter( cell2mat(ent_vpa_chk(2:end,3)), cell2mat(ent_vpa_chk(2:end,5))); hold on;
title('VPA -BY- SURF FA')

subplot(2,2,4)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,3))); hold on;
title('ENTORHINAL RHS -BY- LM2')

tightfig();

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at LEFT HEMISPHERE
%% Find LHS ENTORHINAL VERTEX - 5483, r=.4108, p=.08
rvl_lhs = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/wmparc_FA_wm_lh_entorhinal/rvalues_lhs_dep_var_3_left.mat'); %%% Find LHS Vertex
[ ~, lhs_ent_ind ] = max(rvl_lhs.rvalues);

rvl_rhs = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/wmparc_FA_wm_rh_entorhinal/rvalues_rhs_dep_var_3_left.mat'); %%% Find RHS Vertex
[ ~, rhs_ent_ind ] = max(rvl_rhs.rvalues);

fcfg = []; %%% make plot

fcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test';
fcfg.out_pre_fix = 'lhs_entorhinal_cortex';

fcfg.vtx_cor = { [ lhs_ent_ind          ] [ rhs_ent_ind          ] };
fcfg.vtx_col = { { rgb('bright red') } { rgb('bright green') } };

fcfg.hme_wrk = 1;

ejk_highlight_vertex(fcfg)

%% Scatter LHS - LM2 
ent_ent_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/wmparc_FA_wm_lh_entorhinal/exact_data.csv');
figure()
subplot(2,2,1)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_ent_chk(2:end,6))); hold on;
scatter( cell2mat(ent_ent_chk(17,3)), cell2mat(ent_ent_chk(17,6)),'r');
title('ENTORHINAL LHS -BY- SURF FA')

ent_lm2_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/Logical_Memory_II_RCI/exact_data.csv');
subplot(2,2,2)
scatter( cell2mat(ent_lm2_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,6))); hold on;
scatter( cell2mat(ent_lm2_chk(17,3)), cell2mat(ent_lm2_chk(17,6)),'r');
title('LM2 -BY- SURF FA')

subplot(2,2,3)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,3))); hold on;
scatter( cell2mat(ent_ent_chk(17,3)), cell2mat(ent_lm2_chk(17,3)),'r');
title('ENTORHINAL LHS -BY- LM2')

tightfig();

%% Scatter LHS - LM1
ent_ent_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/wmparc_FA_wm_lh_entorhinal/exact_data.csv');
figure()
subplot(2,2,1)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_ent_chk(2:end,6))); hold on;
scatter( cell2mat(ent_ent_chk(17,3)), cell2mat(ent_ent_chk(17,6)),'r');
title('ENTORHINAL LHS -BY- SURF FA')

ent_lm2_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/Logical_Memory_I_RCI/exact_data.csv');
subplot(2,2,2)
scatter( cell2mat(ent_lm2_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,6))); hold on;
scatter( cell2mat(ent_lm2_chk(17,3)), cell2mat(ent_lm2_chk(17,6)),'r');
title('LM1 -BY- SURF FA')

subplot(2,2,3)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,3))); hold on;
scatter( cell2mat(ent_ent_chk(17,3)), cell2mat(ent_lm2_chk(17,3)),'r');
title('ENTORHINAL LHS -BY- LM1')

tightfig();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find LHS ENTORHINAL VERTEX - 5483
ent_ent_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/wmparc_FA_wm_lh_entorhinal/exact_data.csv');
figure()
subplot(2,2,1)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_ent_chk(2:end,6))); hold on;
title('ENTORHINAL LHS -BY- SURF FA')

ent_lm2_chk = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/Logical_Memory_II_RCI/exact_data.csv');
subplot(2,2,2)
scatter( cell2mat(ent_lm2_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,6))); hold on;
title('LM2 -BY- SURF FA')

subplot(2,2,3)
scatter( cell2mat(ent_ent_chk(2:end,3)), cell2mat(ent_lm2_chk(2:end,3))); hold on;
title('ENTORHINAL LHS -BY- LM2')

tightfig();

ttt = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/Logical_Memory_II_RCI/rvalues_lhs_dep_var_3_left.mat');
tt2 =           load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/Logical_Memory_II_RCI/rvalues_lhs_dep_var_3_left.mat');

ppp = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/Logical_Memory_II_RCI/pvalues_lhs_dep_var_3_left.mat');
pp2 =           load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/Logical_Memory_II_RCI/pvalues_lhs_dep_var_3_left.mat');

ttt.rvalues(5483)
tt2.rvalues(5483)

ppp.pvalues(5483)
pp2.pvalues(5483)


ttt = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/Logical_Memory_I_RCI/rvalues_lhs_dep_var_3_left.mat');
tt2 =           load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/Logical_Memory_I_RCI/rvalues_lhs_dep_var_3_left.mat');

ppp = load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_no_epd006/Logical_Memory_I_RCI/pvalues_lhs_dep_var_3_left.mat');
pp2 =           load('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/Logical_Memory_I_RCI/pvalues_lhs_dep_var_3_left.mat');

ttt.rvalues(5483)
tt2.rvalues(5483)

ppp.pvalues(5483)
pp2.pvalues(5483)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttt = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_wmparc_fa_lhs_sm313.mat');

subplot(4,2,1)
hist(ttt.srf_dta(30,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{30})

subplot(4,2,2)
hist(ttt.srf_dta(32,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{32})

subplot(4,2,3)
hist(ttt.srf_dta(33,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{33})

subplot(4,2,4)
hist(ttt.srf_dta(34,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{34})

subplot(4,2,5)
hist(ttt.srf_dta(35,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{35})

subplot(4,2,6)
hist(ttt.srf_dta(36,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{36})

subplot(4,2,7)
hist(ttt.srf_dta(37,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{37})

subplot(4,2,8)
hist(ttt.srf_dta(38,:),1000)
xlim([0 .50])
ylim([0 2000])
title(ttt.srf_dta_sbj{38})

tightfig();
print('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_Cor_Test/FA_check.png','-dpng')
close all;

















