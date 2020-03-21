%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

sbj_nme = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/sbj_nme_no_nan.csv');
lhs_dta = load('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surf_aMRI_thickness_lhs_no_nan.mat');
rhs_dta = load('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surf_aMRI_thickness_rhs_no_nan.mat');

%
cfg = [];

cfg.sbj_grp     = sbj_nme;
cfg.sbj_grp_col = { 'Diagnosis' };
cfg.sbj_grp_nme = { { 'HC' 'EPD_Old' 'MCI'} };
cfg.sbj_grp_cmp = { { [2 1] [3 1] [2 3] } };

cfg.dta_lbl = { 'corticalthickness' };
cfg.dta     = { lhs_dta rhs_dta }; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

cfg.hms     = {'lhs' 'rhs'};

%
fcfg = [];
fcfg.grp_fle = cfg.sbj_grp;
grp_fle      = ejk_collate_groups(fcfg);

% Get and Collate Data
for iHM = 1:numel(cfg.hms)
    
    hms_nme = cfg.hms{iHM};
    
    for iGR = 1:numel(cfg.sbj_grp_col)
        
        grp_ovr_nme = cfg.sbj_grp_col{iGR};
        
        for iCD = 1:numel(cfg.sbj_grp_nme{iGR})
            
            grp_cde_nme = cfg.sbj_grp_nme{iGR}{iCD};
            grp_cde_ind = find( strcmpi( grp_fle.table.(grp_ovr_nme)(:,1) , cfg.sbj_grp_nme{iGR}{iCD} ) );
            grp_cde_sbj = grp_fle.sbj_nme( grp_fle.(grp_ovr_nme) == grp_cde_ind );            
            [ ~ , grp_dta_ind ] = intersect( cfg.dta{iHM}.srf_dta_sbj , grp_cde_sbj );
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).sbj_nme = cfg.dta{iHM}.srf_dta_sbj(grp_dta_ind);
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta  = cfg.dta{iHM}.srf_dta(grp_dta_ind,:);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum    = numel(grp_dta_ind) - sum(isnan(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta(:,1))); 
                         
            fprintf('N for diagnosis %s, hemi %s: %d\n',grp_cde_nme,hms_nme,grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
            vec_sum = nansum( grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta , 1 );
            vec_sm2 = nansum( grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta.^2 , 1 );
            vec_num = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum ;
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecmean   = vec_sum / vec_num;
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd    = sqrt( ( vec_num * vec_sm2 - vec_sum.^2 ) ./ ( vec_num * (vec_num - 1) ) );
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstderr = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd / sqrt(vec_num);
            
        end
    end
end

% Differences
plt_nme_ind = 1;

for iGR = 1:numel(cfg.sbj_grp_col)
    
    grp_ovr_nme = cfg.sbj_grp_col{iGR};
        
    for iC = 1:numel(cfg.sbj_grp_cmp{iGR})
        
        grp_cde_one_nme = cfg.sbj_grp_nme{iGR}{ cfg.sbj_grp_cmp{iGR}{iC}(1) };
        grp_cde_two_nme = cfg.sbj_grp_nme{iGR}{ cfg.sbj_grp_cmp{iGR}{iC}(2) };
        
        % absolute difference
        vecdiff   = grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
        vecstderr = sqrt(grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
        plt_men_dat{1}{plt_nme_ind}   = vecdiff;
        plt_tvl_dat{1}{plt_nme_ind}   = vecdiff ./ vecstderr;
        plt_pvl_dat{1}{plt_nme_ind}   = (1 - normcdf(abs(plt_tvl_dat{1}{plt_nme_ind})))*2;
        
        vecdiff   = grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
        vecstderr = sqrt(grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
        plt_men_dat{2}{plt_nme_ind}   = vecdiff;
        plt_tvl_dat{2}{plt_nme_ind}   = vecdiff ./ vecstderr;
        plt_pvl_dat{2}{plt_nme_ind}   = (1 - normcdf(abs(plt_tvl_dat{2}{plt_nme_ind})))*2;
                
        srf_typ_num(plt_nme_ind)     = grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).nsum + grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).nsum - 2;
        srf_typ_out_men_pre{plt_nme_ind} = sprintf('Cort_%s_diff_mean_%s_%s.mgh',cfg.dta_lbl{iGR},grp_cde_one_nme,grp_cde_two_nme);
        plt_nme_ind = plt_nme_ind + 1;
                
    end
    
end

%%%%%%%%%%%%%%%% TTEST SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TLE vs HC
dta_ind = 1;
plt_nme = 'ZTEST_TLE_vs_HC';

pvl_lhs     = plt_pvl_dat{1}{dta_ind};
pvl_rhs     = plt_pvl_dat{2}{dta_ind};

men_dff_lhs = plt_men_dat{1}{dta_ind};
men_dff_rhs = plt_men_dat{2}{dta_ind};

lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

%
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/diff/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)
  
%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr men_dff_rhs_fdr };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%% MCI vs HC
dta_ind = 2;
plt_nme = 'ZTEST_MCI_vs_HC';

pvl_lhs     = plt_pvl_dat{1}{dta_ind};
pvl_rhs     = plt_pvl_dat{2}{dta_ind};

men_dff_lhs = plt_men_dat{1}{dta_ind};
men_dff_rhs = plt_men_dat{2}{dta_ind};

lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

%
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/diff/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)
  
%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr men_dff_rhs_fdr };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%% TLE vs MCI
dta_ind = 3;
plt_nme = 'ZTEST_TLE_vs_MCI';

pvl_lhs     = plt_pvl_dat{1}{dta_ind};
pvl_rhs     = plt_pvl_dat{2}{dta_ind};

men_dff_lhs = plt_men_dat{1}{dta_ind};
men_dff_rhs = plt_men_dat{2}{dta_ind};

lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

%
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)
  
%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs men_dff_rhs };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr men_dff_rhs_fdr };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%%%%%%%%%%%%%%%% ANCOVA SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TLE vs HC
plt_nme = 'ANCOVA_TLE_vs_HC';

pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean;
    
lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

% %
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/diff/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%% MCI vs HC
plt_nme = 'ANCOVA_MCI_vs_HC';

pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_MCI_HC.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_MCI_HC.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_MCI_HC.mat');
    men_dff_lhs = men_dff_lhs.emmean*-1;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_MCI_HC.mat');
    men_dff_rhs = men_dff_rhs.emmean*-1;
    
lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

% %
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/diff/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)

%% TLE vs MCI
plt_nme = 'ANCOVA_TLE_MCI';

pvl_lhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueL_TLE_MCI.mat');
    pvl_lhs = pvl_lhs.pValue;
pvl_rhs     = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/pValueR_TLE_MCI.mat');
    pvl_rhs = pvl_rhs.pValue;
    
men_dff_lhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanL_TLE_MCI.mat');
    men_dff_lhs = men_dff_lhs.emmean;
men_dff_rhs = load('/home/jxrao/Desktop/Lab/surface_map/data/R_Output/emmeanR_TLE_MCI.mat');
    men_dff_rhs = men_dff_rhs.emmean;
    
lhs_fdr = .005; %FDR(pvl_lhs,.05);
rhs_fdr = .005; %FDR(pvl_rhs,.05);

% %
lhs_bdd_ind = find(pvl_lhs>lhs_fdr);
rhs_bdd_ind = find(pvl_rhs>rhs_fdr);

men_dff_lhs_fdr = men_dff_lhs;
men_dff_rhs_fdr = men_dff_rhs;

men_dff_lhs_fdr(lhs_bdd_ind) = 0;
men_dff_rhs_fdr(rhs_bdd_ind) = 0;

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/diff/';
pcfg.out_pre_fix = [ plt_nme '_thick01'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.low_rng_num = [ -0.01 0.01 ];
pcfg.hgh_rng_num = [ -0.20 0.20 ];

mmil_anat_surf_plot(pcfg)

%
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_thick05'];

pcfg.plt_dta = { men_dff_lhs' men_dff_rhs' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.05 0.05 ];

mmil_anat_surf_plot(pcfg)

%  
pcfg = [];

pcfg.out_dir     = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/ANCOVA_TEST/';
pcfg.out_pre_fix = [ plt_nme '_p005'];

pcfg.plt_dta = { men_dff_lhs_fdr' men_dff_rhs_fdr' };

pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

pcfg.fmr_rng_num = [ -0.01 0.01 ];

mmil_anat_surf_plot(pcfg)





