%%
clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';

sbj_nme = 'epilepsy_aging_graph_sample_v2.csv';
    sbj_nme = mmil_readtext( [prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' sbj_nme ]);

prc_nme = { '' '.a2009s' '.split' };

gbl_mes_nme = { 'PathEfficiency' 'Transitivity' 'Modularity' 'LocalPathEfficiency' 'ClusteringCoefficient' 'Degree' };
gbl_mes_plt = { 'p01'            'p02'          'p03'        'p04'                 'p05'                   'p06'    };

grp_plt     = [ 1                1              1            0                     0                       1        ];
reg_plt     = [ 0                0              0            1                     1                       1        ];

bin_cut_off = 10:5:50;
bin_dir     = 'pos';
bin_foc     = 4; % 

sbj_grp_col = {   'Laterality'                                                  'MTS'   }; % {   'Cognitive'                                              'Diagnosis'                                'Onset'                                   };
sbj_grp_nme = { { 'HC'       'MCI'         'Left'        'Right'       'Bilateral' }  { 'HC'        'MCI'         'Yes'         'No' } }; % { { 'HC'       'MCI'         'TLE_NI'      'TLE_MND' }     { 'HC'       'MCI'         'EPD_Old' }     { 'HC'       'MCI'         'Early'       'Late' }        };
sbj_grp_clr = { { [.8 .8 .8] [.47 .31 .45] [.63 .77 .80] [.00 .33 .37] [.04 .45 .45] } { [.8 .8 .8] [.47 .31 .45] [.63 .77 .80] [.00 .33 .37] } }; % { { [.8 .8 .8] [.47 .31 .45] [.63 .77 .80] [.00 .33 .37] } { [.8 .8 .8] [.47 .31 .45] [.34 .59 .64] } { [.8 .8 .8] [.47 .31 .45] [.63 .77 .80] [.00 .33 .37] } };

num_rep = 5000;
pvl_lvl = .01;
pvl_con = {    'HC'                                           'HC'              }; %{    'HC'                                   'HC'               'HC' };
pvl_cmp = {  { 'MCI'       'Left'     'Right' 'Bilateral' } { 'MCI' 'Yes' 'No'} }; % {  { 'MCI'       'TLE_NI'     'TLE_MND' } { 'MCI' 'EPD_Old'} { 'MCI'         'Early'       'Late' } };
lve_num_out = 3;

%% Plot Example ROIs
for iPR = 1:numel(prc_nme)
    
    fcfg = [];
    
    fcfg.prj_dir = '/home/ekaestne/PROJECTS/';
    fcfg.prc_nme = prc_nme{iPR};
    
    fcfg.fsr_dir = '/home/mmilmcdRSI/data/fsurf';
    fcfg.fsr_nme = 'fsaverage';
    
    fcfg.out_dir = [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/ROI/' 'roi_plot' '_' mmil_spec_char(prc_nme{iPR},{'.'}) ];
    fcfg.out_nme = prc_nme{iPR};
    
    ejk_roi_plot(fcfg)
    
end

%% Load ROIs
for iPR = 1:numel(prc_nme)
    
    fcfg = [];
    
    fcfg.prj_dir = prj_dir;
    fcfg.prj_nme = prj_nme;
    
    fcfg.sbj_nme = sbj_nme;
    
    fcfg.dta_nme = [ 'MRI' '_' 'thickness' ];
    fcfg.prc_nme = prc_nme{iPR};
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'ROI'];
    
    ejk_project_pull_roi(fcfg)
    
end

%% Run Initial Graph Measures
dta_nme = { 'Thickness_Desikan_noTransverse' ...
            'Thickness_Desikan' ...
            'Thickness_Destrieaux' ...
            'Thickness_Desikan_Modified' ...
            'Thickness_Desikan_Modified_Additional'};

dta_hld = { 'MRI_thickness_aparc__Epilepsy_and_Aging_updated_noTransverse.csv' ...
            'MRI_thickness_aparc__Epilepsy_and_Aging_updated_withTransverse.csv' ...
            'MRI_thickness_aparc_xa2009s_Epilepsy_and_Aging_updated.csv' ...
            'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ...
            'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' };

prc_nme_hld = [ 1 ...
                1 ...
                2 ...
                3 ...
                3 ];

for iDT = 5 %3:numel(dta_nme)
    
    out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Graph_Scratch' '/' dta_nme{iDT} ];
    
    % Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thk_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'ROI' '/' dta_hld{iDT} ] );
    thk_sbj_nme = thk_dta(2:end,1);
    thk_roi_nme = thk_dta(1,2:end);
    thk_dta     = cell2mat(thk_dta(2:end,2:end));

    % Run Graph Theory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fcfg = [];
    
    fcfg.sbj_nme = sbj_nme;
    fcfg.prj_dir = prj_dir;
    
    fcfg.prc_nme = prc_nme{ prc_nme_hld(iDT) };
    
    fcfg.sbj_grp_col = sbj_grp_col;
    fcfg.sbj_grp_nme = sbj_grp_nme;
    fcfg.sbj_grp_clr = sbj_grp_clr;
    
    fcfg.thk_dta     = thk_dta;
    fcfg.thk_roi_nme = thk_roi_nme;
    
    fcfg.pvl_con     = pvl_con;
    fcfg.num_rep     = num_rep;
    fcfg.lve_num_out = lve_num_out;
    fcfg.pvl_lvl     = pvl_lvl;
    fcfg.pvl_cmp     = pvl_cmp;
    
    fcfg.gbl_mes_nme = gbl_mes_nme;
    fcfg.gbl_mes_plt = gbl_mes_plt;
    
    fcfg.grp_plt     = grp_plt;
    fcfg.reg_plt     = reg_plt;
    
    fcfg.bin_cut_off = bin_cut_off;
    fcfg.bin_dir     = bin_dir;
    fcfg.reg_col     = bin_foc; % 25
    
    fcfg.out_dir     = out_dir;
    
    cfg.ovr_wrt = 0;
    
    GraphWrapper(fcfg)
    
end
















