clear; clc;

dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData';
out_put = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/HippocampalInvestigation/';

sbj_col = { 'SubjId' };

dem_col = { 'VisitID' 'FieldStrength' 'B0scan' 'SideOfSeizureFocus' 'Engel_Outcome_Simple' };

mem_col = { 'CVLT_Total_RCI' 'CVLT_LDFR_RCI' 'Logical_Memory_I_RCI' 'Logical_Memory_II_RCI' };

mem_grp_col = { 'CVLT_Total_1SD_Decline' 'CVLT_LDFR_1SD_Decline' 'Logical_Memory_I_1SD_Decline' 'Logical_Memory_II_1SD_Decline' ...
                'CVLT_Total_80CI_Decline' 'CVLT_LDFR_80CI_Decline' 'Logical_Memory_I_80CI_Decline' 'Logical_Memory_II_80CI_Decline'};
        
mri_col     = { 'subcort_vol_Left_Hippocampus_Total' 'subcort_vol_Right_Hippocampus_Total' ...
                'subcort_vol_Left_Amygdala_Total'    'subcort_vol_Right_Amygdala_Total'    };
mri_lat_col = { 'LI_hippocampalvolume' ...
                'LI_amygdalavolume' };        
        
dti_col     = { 'aseg_FA_Left_Hippocampus' 'aseg_FA_Right_Hippocampus' ...
                'aseg_FA_Left_Amygdala'    'aseg_FA_Right_Amygdala'    ...
                'aseg_MD_Left_Hippocampus' 'aseg_MD_Right_Hippocampus' ...
                'aseg_MD_Left_Amygdala'    'aseg_MD_Right_Amygdala'    };
dti_lat_col = { 'aseg_FA_LI_Hippocampus' ...
                'aseg_FA_LI_Amygdala' ...
                'aseg_MD_LI_Hippocampus' ...
                'aseg_MD_LI_Amygdala' };

ejk_chk_dir(out_put)

%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri_dta = mmil_readtext([ dta_dir '/' 'T1.csv']);
    mri_dta( cellfun(@isempty,mri_dta) ) = {NaN};
    mri_dta( cellfun(@(x) strcmpi(x,'NA'),mri_dta) ) = {NaN}; 
    sbj_col_ind_mri_dta = find(strcmpi( mri_dta(1,:), sbj_col));
dti_dta = mmil_readtext([ dta_dir '/' 'FA_MD.csv']);
    dti_dta( cellfun(@isempty,dti_dta) ) = {NaN};
    dti_dta( cellfun(@(x) strcmpi(x,'NA'),dti_dta) ) = {NaN}; 
    sbj_col_ind_dti_dta = find(strcmpi( dti_dta(1,:), sbj_col)); 
scr_dta = mmil_readtext([ dta_dir '/' 'RCIs_ClinicalData_group.csv']);
    scr_dta( cellfun(@isempty,scr_dta) ) = {'N/A'};
    scr_dta( cellfun(@(x) strcmpi(x,'NA'),scr_dta) ) = {'N/A'}; 
    sbj_col_ind_scr_dta = find(strcmpi( scr_dta(1,:), sbj_col)); 

% Demographics
dem_col_ind = find(ismember( mri_dta(1,:), dem_col ));

dem_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta dem_col_ind]);
dem_sbj_sbj_nme = dem_sbj_dta(2:end,1);
dem_sbj_roi_nme = dem_sbj_dta(1,2:end);
dem_sbj_dta     = dem_sbj_dta(2:end,2:end); 

% Memory Scores
mem_col_ind = find(ismember( mri_dta(1,:), mem_col ));

mem_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta mem_col_ind]);
mem_sbj_sbj_nme = mem_sbj_dta(2:end,1);
mem_sbj_roi_nme = mem_sbj_dta(1,2:end);
mem_sbj_dta     = cell2mat(mem_sbj_dta(2:end,2:end)); 

% Memory Impairment
mem_grp_col_ind = find(ismember( scr_dta(1,:), mem_grp_col ));

mem_grp_sbj_dta     = scr_dta(:,[sbj_col_ind_mri_dta mem_grp_col_ind]);
mem_grp_sbj_sbj_nme = mem_grp_sbj_dta(2:end,1);
mem_grp_sbj_roi_nme = mem_grp_sbj_dta(1,2:end);
mem_grp_sbj_dta     = mem_grp_sbj_dta(2:end,2:end); 

% Hippocampal T1
mri_col_ind = find(ismember( mri_dta(1,:), mri_col ));

mri_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta mri_col_ind]);
mri_sbj_sbj_nme = mri_sbj_dta(2:end,1);
mri_sbj_roi_nme = mri_sbj_dta(1,2:end);
mri_sbj_dta     = mri_sbj_dta(2:end,2:end); 

% Hippocampal TI - Laterality Index
mri_lat_col_ind = find(ismember( mri_dta(1,:), mri_lat_col ));

mri_lat_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta mri_lat_col_ind]);
mri_lat_sbj_sbj_nme = mri_lat_sbj_dta(2:end,1);
mri_lat_sbj_roi_nme = mri_lat_sbj_dta(1,2:end);
mri_lat_sbj_dta     = mri_lat_sbj_dta(2:end,2:end); 

% Hippocampal DTI
dti_col_ind = find(ismember( dti_dta(1,:), dti_col ));

dti_sbj_dta     = dti_dta(:,[sbj_col_ind_dti_dta dti_col_ind]);
dti_sbj_sbj_nme = dti_sbj_dta(2:end,1);
dti_sbj_roi_nme = dti_sbj_dta(1,2:end);
dti_sbj_dta     = dti_sbj_dta(2:end,2:end); 

% Hippocampal DTI - Laterality Index
dti_lat_col_ind = find(ismember( dti_dta(1,:), dti_lat_col ));

dti_lat_sbj_dta     = dti_dta(:,[sbj_col_ind_dti_dta dti_lat_col_ind]);
dti_lat_sbj_sbj_nme = dti_lat_sbj_dta(2:end,1);
dti_lat_sbj_roi_nme = dti_lat_sbj_dta(1,2:end);
dti_lat_sbj_dta     = dti_lat_sbj_dta(2:end,2:end); 

% Combine NeuroBio
neu_bio_dta     = cell2mat([ mri_sbj_dta     dti_sbj_dta ]);
neu_bio_roi_nem = [ mri_sbj_roi_nme dti_sbj_roi_nme ];

neu_bio_lat_dta     = cell2mat([ mri_lat_sbj_dta     dti_lat_sbj_dta ]);
neu_bio_lat_roi_nem = [ mri_lat_sbj_roi_nme dti_lat_sbj_roi_nme ];

%% Select Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select subjects
scn_str_col = find(strcmpi( dem_sbj_roi_nme, 'FieldStrength'));
sze_sde_col = find(strcmpi( dem_sbj_roi_nme, 'SideOfSeizureFocus'));

ts3_lft = find( cell2mat(dem_sbj_dta(:,scn_str_col))==3 & strcmpi( dem_sbj_dta(:,sze_sde_col), 'left'));
ts3_rgh = find( cell2mat(dem_sbj_dta(:,scn_str_col))==3 & strcmpi( dem_sbj_dta(:,sze_sde_col), 'right'));

lft = find( strcmpi( dem_sbj_dta(:,sze_sde_col), 'left'));
rgh = find( strcmpi( dem_sbj_dta(:,sze_sde_col), 'right'));

% Get the groups


%% Run Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left 3T
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = mem_sbj_dta( ts3_lft, :);
fcfg.lbl_one = mem_sbj_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta( ts3_lft, :);
fcfg.lbl_two = neu_bio_roi_nem;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_3T' '/'];

ejk_cross_cor( fcfg );

% Right 3T
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_rgh, 1);

fcfg.dta_one = mem_sbj_dta( ts3_rgh, :);
fcfg.lbl_one = mem_sbj_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta( ts3_rgh, :);
fcfg.lbl_two = neu_bio_roi_nem;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Right_3T' '/'];

ejk_cross_cor( fcfg );

% Left All
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( lft, 1);

fcfg.dta_one = mem_sbj_dta( lft, :);
fcfg.lbl_one = mem_sbj_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_lat_dta( lft, :);
fcfg.lbl_two = neu_bio_lat_roi_nem;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_all' '/'];

ejk_cross_cor( fcfg );

% Right All
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( rgh, 1);

fcfg.dta_one = mem_sbj_dta( rgh, :);
fcfg.lbl_one = mem_sbj_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_lat_dta( rgh, :);
fcfg.lbl_two = neu_bio_lat_roi_nem;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Right_all' '/'];

ejk_cross_cor( fcfg );

%% Run Decliner/Non-Decliner Comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta     = neu_bio_dta( ts3_lft, :);
fcfg.dta_nme = neu_bio_roi_nem;

fcfg.grp     = mem_grp_sbj_dta( ts3_lft, :);
fcfg.grp_nme = mem_grp_sbj_roi_nme;

fcfg.out_dir = [ out_put '/' 'Left_ttest_3T' '/'];

ejk_ttest2_independent( fcfg );

% Right 3T
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_rgh, 1);

fcfg.dta     = neu_bio_dta( ts3_rgh, :);
fcfg.dta_nme = neu_bio_roi_nem;

fcfg.grp     = mem_grp_sbj_dta( ts3_rgh, :);
fcfg.grp_nme = mem_grp_sbj_roi_nme;

fcfg.out_dir = [ out_put '/' 'Right_ttest_3T' '/'];

ejk_ttest2_independent( fcfg );

% Left All
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( lft, 1);

fcfg.dta     = neu_bio_lat_dta( lft, :);
fcfg.dta_nme = neu_bio_lat_roi_nem;

fcfg.grp     = mem_grp_sbj_dta( lft, :);
fcfg.grp_nme = mem_grp_sbj_roi_nme;

fcfg.out_dir = [ out_put '/' 'Left_ttest_All' '/'];

ejk_ttest2_independent( fcfg );

% Right All
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( rgh, 1);

fcfg.dta     = neu_bio_lat_dta( rgh, :);
fcfg.dta_nme = neu_bio_lat_roi_nem;

fcfg.grp     = mem_grp_sbj_dta( rgh, :);
fcfg.grp_nme = mem_grp_sbj_roi_nme;

fcfg.out_dir = [ out_put '/' 'Right_ttest_All' '/'];

ejk_ttest2_independent( fcfg );

%% Logistic Regression
tst_num_hld = [ 2 4 6 8 ];

for iTN = 1:numel(tst_num_hld)
    
    tst_num = tst_num_hld(iTN);
    
    % L-TLE
    lcfg = [];
    
    lcfg.sbj_grp = [ mem_sbj_sbj_nme( lft, 1)  mem_grp_sbj_dta( lft, tst_num) ];
    lcfg.lbl_ord = { 'Decline' 'NoDecline' };
    
    lcfg.dta     = [ mem_sbj_sbj_nme( lft, 1) num2cell(neu_bio_lat_dta( lft, :)) ];
    lcfg.dta_lbl = { 'sbj_nme' neu_bio_lat_roi_nem{:} };
    
    lcfg.mdl     = { { 'LI_hippocampalvolume' } ...
        { 'aseg_MD_LI_Hippocampus' } ...
        { 'LI_hippocampalvolume' 'aseg_MD_LI_Hippocampus' } };
    lcfg.mdl_nme = { 'Volume' ...
        'MD' ...
        'Volume+MD' };
    
    lcfg.nrm_grp = 0;
    
    lcfg.mdl_cmp_plt = { [1 2 3] };
    lcfg.mld_cmp_col = { { rgb('green') rgb('red') rgb('purple') } };
    lcfg.mdl_cmp_nme = { [ mem_grp_sbj_roi_nme{tst_num} '_' 'Hippocampus'] };
    
    lcfg.out_dir = [ out_put '/' 'Left_All_Logistic' '/' mem_grp_sbj_roi_nme{tst_num}];
    
    mmil_log_reg(lcfg)
    
end
