clear; clc;

%% Setup Paths & I
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_graph_sample.csv';

thk_dta = 'MRI_thickness_aparc_DESIKAN.csv';

bin_cut_off = 10:5:50;
bin_dir     = 'pos';

sbj_grp_col = {   'Diagnosis' };
sbj_grp_nme = { { 'HC'       'MCI'         'EPD_Old' } };
sbj_grp_clr = { { [.8 .8 .8] [.51 .29 .47] [.33 .59 .64] } };

gbl_mes_nme = { '' '' '' };
gbl_mes_plt = { '' '' '' };

num_rep = 5000;
pvl_lvl = .05;
pvl_con = { 'HC' };
pvl_cmp = { { 'MCI' 'EPD_Old'} };
lve_num_out = 3;
pvl_lvl = .01;

%% Load Data
sbj_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' sbj_nme ] );


thk_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' thk_dta ] );
    thk_sbj_nme = thk_dta(2:end,1);
    thk_roi_nme = thk_dta(1,2:end);
    thk_dta     = cell2mat(thk_dta(2:end,2:end));

%% Calculate Matrix
clear cor_mtx con_smp dff_smp

fcfg = [];

fcfg.sbj_grp = sbj_nme;

fcfg.sbj_grp_col = sbj_grp_col;
fcfg.sbj_grp_nme = sbj_grp_nme;

fcfg.dta     = thk_dta;
fcfg.roi_nme = thk_roi_nme;

fcfg.clc_con_smp_mtx = 1;
fcfg.con_smp_nme = pvl_con;
fcfg.con_smp_rep = num_rep;
fcfg.lve_num_out = lve_num_out;

fcfg.clc_dff_mtx = 0;

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Graph_Scratch' ];
fcfg.plt = 0;

[ cor_mtx , con_smp , ~ ] = GraphThickMatrix(fcfg);

%% Binarize
fcfg = [];

fcfg.sbj_grp_col = sbj_grp_col;
fcfg.sbj_grp_nme = sbj_grp_nme;

fcfg.bin_dir     = 'pos';
fcfg.bin_cut_off = bin_cut_off;

fcfg.cor_mtx         = cor_mtx;

fcfg.clc_con_smp_mtx = 1;
fcfg.con_smp         = con_smp;

fcfg.clc_dff_mtx = 0;

[ cor_mtx_bin , con_smp_bin , ~ ] = GraphBinzarize(fcfg);

%% Generate Distribution Plots

%% Generate Square Plots

%% Calcuation & Plots
gbl_mes_nme = { 'CharacteristicPathLength' 'PathEfficiency' 'Transitivity' 'Modularity' 'LocalPathEfficiency' 'ClusteringCoefficient' 'Degree' 'BetweennessCentrality'};
gbl_mes_plt = { 'p10'                      'p11'            'p12'          'p13'        'p14'                 'p15'                   'p16'    'p17' };

grp_plt     = [ 1                          1                1              1            0                     0                       1        0     ];
reg_plt     = [ 0                          0                0              0            1                     1                       1        1     ];

for iGL = 1:numel(gbl_mes_nme)
    
    % Calculate
    fcfg = [];
    
    fcfg.mes = gbl_mes_nme{iGL};
    
    fcfg.dta = cor_mtx_bin;
    
    fcfg.con_spr_stt = 1;
    fcfg.con_spr_dta = con_smp_bin;
    
    fcfg.dff_cmp_stt = 0;
    
    plt_dta = GraphCalc(fcfg);
    
    % Plot
    for iC = 1:numel(sbj_grp_col); sve_nme{iC} = [ gbl_mes_plt{iGL} '_' gbl_mes_nme{iGL} '_' sbj_grp_col{iC}]; end
    
    fcfg = [];
    
    fcfg.grp_ovr = sbj_grp_col;
    fcfg.grp_nme = sbj_grp_nme;
    fcfg.grp_col = sbj_grp_clr;
    
    fcfg.leg_loc = 'southwest';
    
    fcfg.bin_cut = bin_cut_off;
    
    fcfg.grp_dta = plt_dta;
    fcfg.thk_nme = thk_roi_nme;
    
    fcfg.con_typ = 'con_shf';
    fcfg.con_nme = pvl_con;
    fcfg.cmp_nme = pvl_cmp;
    
    fcfg.con_spr_plt = 1;
    
    fcfg.con_cmp_stt = 0;
    
    fcfg.pvl_lvl = pvl_lvl;
    
    fcfg.grp_plt     = grp_plt(iGL);
    fcfg.reg_lne_plt = reg_plt(iGL);
    fcfg.reg_srf_plt = 0;
    
    fcfg.ttl = gbl_mes_nme{iGL};
    fcfg.sve_loc = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Graph_Scratch' '/' 'Plots' ];
    fcfg.sve_nme = sve_nme;
    
    GraphLinePlot(fcfg)
    
    %
    clear plt_dta
    
end

%% Default Mode Network
fcfg = [];

fcfg.grp_ovr = sbj_grp_col;
fcfg.grp_nme = sbj_grp_nme;
fcfg.grp_col = sbj_grp_clr;

fcfg.leg_loc = 'southwest';

fcfg.bin_cut = bin_cut_off;

fcfg.grp_dta = plt_dta;
fcfg.thk_nme = thk_roi_nme;
cfg.reg_inc_plt = { 'precuneus' 'superiorparietal' 'medialorbitofrontal' ''  };

fcfg.con_typ = 'con_shf';
fcfg.con_nme = pvl_con;
fcfg.cmp_nme = pvl_cmp;

fcfg.con_spr_plt = 1;

fcfg.con_cmp_stt = 0;

fcfg.pvl_lvl = pvl_lvl;

fcfg.grp_plt     = grp_plt(iGL);
fcfg.reg_lne_plt = reg_plt(iGL);
fcfg.reg_srf_plt = 0;

fcfg.ttl = gbl_mes_nme{iGL};
fcfg.sve_loc = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Graph_Scratch' '/' 'Plots' ];
fcfg.sve_nme{1} = [sve_nme{1} '_' 'DefaultModeNetwork'];

GraphLinePlot(fcfg)


%%


































