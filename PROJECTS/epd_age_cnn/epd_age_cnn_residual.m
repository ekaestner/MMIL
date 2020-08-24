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

%% Calculate Age Residual
% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
thk_sbj_nme = thk_dta(2:end,1);
thk_roi_nme = thk_dta(1,2:end);
thk_dta     = cell2mat(thk_dta(2:end,2:end));

cov_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/epd_age_covariates_no_nan_site_updated.csv');
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

% Collate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_ind = 1;
for iS = 1:size(thk_sbj_nme, 1)
    
    cov_ind = find(strcmpi( cov_sbj_nme, thk_sbj_nme{iS,1} ));
    
    if ~isempty(cov_ind) %&& ~isempty(cov_dta{cov_ind, 5}) % 5 % 6 % 7
        
        thk_dta_sbj_use{num_ind,1} = thk_sbj_nme{iS,1};
        thk_dta_use(num_ind,:)     = thk_dta(iS,:);
               
        cov_dta_use{num_ind,1}     = cov_dta{cov_ind, 1};
        if cov_dta{cov_ind, 2}==1
            cov_dta_use{num_ind,2}     = 'Male';
        elseif cov_dta{cov_ind, 2}==2
            cov_dta_use{num_ind,2}     = 'Female';
        end
        cov_dta_use{num_ind,3}     = cov_dta{cov_ind, 3};
        if cov_dta{cov_ind, 4}==1
           cov_dta_use{num_ind,4}     = '1.5T';
        elseif cov_dta{cov_ind, 4}==2
           cov_dta_use{num_ind,4}     = '3T';
        end

        ste_dta_use{num_ind,1} = cov_dta{cov_ind, 5};
        
        num_ind = num_ind+1;
       
    else
        error('')
    end
end

% Residual Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = thk_dta_sbj_use;

fcfg.dta     = thk_dta_use;
fcfg.dta_nme = thk_roi_nme;

fcfg.cov     = cov_dta_use(:,1);
fcfg.cov_nme = cov_roi_nme(1);

fcfg.out_dir = '/home/ekaestner/Downloads/Braintest/residuals'; % split % desikan

rsd_mtx = ejk_residual_matrix( fcfg );

%%
% Setup
clear thk_dta thk_sbj_nme thk_roi_nme

dta_nme = { 'Thickness_Desikan_Split' ...
            'Thickness_Desikan_Residuals' };

prc_nme_hld = [ 3 ...
                3 ];

thk_dta_hld = mmil_readtext( [ '/home/ekaestne/PROJECTS/' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
thk_sbj_nme{1} = thk_dta_hld(2:end,1);
thk_roi_nme{1} = thk_dta_hld(1,2:end);
thk_dta{1}     = cell2mat(thk_dta_hld(2:end,2:end));            

thk_dta{2} = rsd_mtx;
thk_sbj_nme{2} = thk_sbj_nme{1};
thk_roi_nme{2} = thk_roi_nme{1};
            
for iDT = 1%:2
    
    out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'TestResiduals' '/' dta_nme{iDT} ];
    
    % Run Graph Theory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fcfg = [];
    
    fcfg.sbj_nme = sbj_nme;
    fcfg.prj_dir = prj_dir;
    
    fcfg.prc_nme = prc_nme{ prc_nme_hld(iDT) };
    
    fcfg.sbj_grp_col = sbj_grp_col;
    fcfg.sbj_grp_nme = sbj_grp_nme;
    fcfg.sbj_grp_clr = sbj_grp_clr;
    
    fcfg.thk_dta     = thk_dta{iDT};
    fcfg.thk_roi_nme = thk_roi_nme{iDT};
    
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

