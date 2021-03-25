clear; clc;

dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData';
out_put = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/HippocampalInvestigation_v2/';

sbj_col = { 'SubjId' };

dem_col = { 'VisitID' 'Field Strength' 'SideOfSeizureFocus' 'Engel_Outcome_Simple' }; % 'B0scan'

mem_col = { 'PRE_CVLT_LDFR_Raw' 'PRE_LM_II_SS' 'PRE_VPA_II_SS' 'PRE_BVMT_Delay_Raw' ... 
            'CVLT_LDFR_RCI'     'LM_II_RCI'    'VPA_II_RCI'    'BVMT_Delayed_Recall_RCI' };

mem_grp_col = { 'CVLT_Total_80CI_Decline' 'CVLT_LDFR_80CI_Decline' 'Logical_Memory_I_80CI_Decline' 'Logical_Memory_II_80CI_Decline'};
        
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
dta_hld = mmil_readtext([ dta_dir '/' 'Memory_FA_Updated.csv']);
    dta_hld( cellfun(@isempty,dta_hld) ) = {NaN};
    dta_hld( cellfun(@(x) strcmpi(x,'NA'),dta_hld) ) = {NaN}; 
    sbj_col_dta_hld = find(strcmpi( dta_hld(1,:), sbj_col));


mri_dta = mmil_readtext([ dta_dir '/' 'T1.csv']);
    mri_dta( cellfun(@isempty,mri_dta) ) = {NaN};
    mri_dta( cellfun(@(x) strcmpi(x,'NA'),mri_dta) ) = {NaN}; 
    sbj_col_ind_mri_dta = find(strcmpi( mri_dta(1,:), sbj_col));
dti_dta = mmil_readtext([ dta_dir '/' 'FA_MD.csv']);
    dti_dta( cellfun(@isempty,dti_dta) ) = {NaN};
    dti_dta( cellfun(@(x) strcmpi(x,'NA'),dti_dta) ) = {NaN}; 
    sbj_col_ind_dti_dta = find(strcmpi( dti_dta(1,:), sbj_col)); 
   

dta_hld_mtc = find(~ismember(dta_hld(:,1), mri_dta(:,1)));  
dta_hld(dta_hld_mtc,:) = [];  

% Demographics %%%%%%%%%%%%%%%%%%%%%
dem_col_ind = find(ismember( dta_hld(1,:), dem_col ));

dem_sbj_dta     = dta_hld(:,[sbj_col_dta_hld dem_col_ind]);
dem_sbj_sbj_nme = dem_sbj_dta(2:end,1);
dem_sbj_roi_nme = dem_sbj_dta(1,2:end);
dem_sbj_dta     = dem_sbj_dta(2:end,2:end); 

% Memory Scores %%%%%%%%%%%%%%%%%%%%%
mem_col_ind = find(ismember( dta_hld(1,:), mem_col ));

mem_sbj_dta     = dta_hld(:,[sbj_col_dta_hld mem_col_ind]);
mem_sbj_sbj_nme = mem_sbj_dta(2:end,1);
mem_sbj_roi_nme = mem_sbj_dta(1,2:end);
mem_sbj_dta     = cell2mat(mem_sbj_dta(2:end,2:end)); 

% Memory Impairment %%%%%%%%%%%%%%%%%%%%%
mem_grp_sbj_sbj_nme = mem_sbj_sbj_nme;
mem_grp_sbj_roi_nme = mem_sbj_roi_nme;
mem_grp_sbj_dta     = num2cell(mem_sbj_dta); 

for iR = 1:size(mem_grp_sbj_dta,1); for iC = 1:size(mem_grp_sbj_dta,2)
        if ~isnan(mem_grp_sbj_dta{iR,iC}) && mem_grp_sbj_dta{iR,iC} <= -1
            mem_grp_sbj_dta{iR,iC} = 'Decline';
        elseif ~isnan(mem_grp_sbj_dta{iR,iC}) && mem_grp_sbj_dta{iR,iC} > -1
            mem_grp_sbj_dta{iR,iC} = 'NoDecline';
        else
            mem_grp_sbj_dta{iR,iC} = 'N/A';
        end
    end;end

% Hippocampal T1 %%%%%%%%%%%%%%%%%%%%%
mri_col_ind = find(ismember( mri_dta(1,:), mri_col ));

mri_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta mri_col_ind]);
mri_sbj_sbj_nme = mri_sbj_dta(2:end,1);
mri_sbj_roi_nme = mri_sbj_dta(1,2:end);
mri_sbj_dta     = mri_sbj_dta(2:end,2:end); 

% Hippocampal TI - Laterality Index %%%%%%%%%%%%%%%%%%%%%
mri_lat_col_ind = find(ismember( mri_dta(1,:), mri_lat_col ));

mri_lat_sbj_dta     = mri_dta(:,[sbj_col_ind_mri_dta mri_lat_col_ind]);
mri_lat_sbj_sbj_nme = mri_lat_sbj_dta(2:end,1);
mri_lat_sbj_roi_nme = mri_lat_sbj_dta(1,2:end);
mri_lat_sbj_dta     = mri_lat_sbj_dta(2:end,2:end); 

% Hippocampal DTI %%%%%%%%%%%%%%%%%%%%%
dti_col_ind = find(ismember( dti_dta(1,:), dti_col ));

dti_sbj_dta     = dti_dta(:,[sbj_col_ind_dti_dta dti_col_ind]);
dti_sbj_sbj_nme = dti_sbj_dta(2:end,1);
dti_sbj_roi_nme = dti_sbj_dta(1,2:end);
dti_sbj_dta     = dti_sbj_dta(2:end,2:end); 

% Hippocampal DTI - Laterality Index %%%%%%%%%%%%%%%%%%%%%
dti_lat_col_ind = find(ismember( dti_dta(1,:), dti_lat_col ));

dti_lat_sbj_dta     = dti_dta(:,[sbj_col_ind_dti_dta dti_lat_col_ind]);
dti_lat_sbj_sbj_nme = dti_lat_sbj_dta(2:end,1);
dti_lat_sbj_roi_nme = dti_lat_sbj_dta(1,2:end);
dti_lat_sbj_dta     = dti_lat_sbj_dta(2:end,2:end); 

% Combine NeuroBio %%%%%%%%%%%%%%%%%%%%%
neu_bio_dta     = cell2mat([ mri_sbj_dta     dti_sbj_dta ]);
neu_bio_roi_nme = [ mri_sbj_roi_nme dti_sbj_roi_nme ];
neu_bio_roi_use = [1:2 5:8]; % 1:numel(neu_bio_roi_nme)

neu_bio_lat_dta     = cell2mat([ mri_lat_sbj_dta     dti_lat_sbj_dta ]);
neu_bio_lat_roi_nme = [ mri_lat_sbj_roi_nme dti_lat_sbj_roi_nme ];
neu_bio_lat_roi_use = [1 3:4]; % 1:numel(neu_bio_lat_roi_nme)

%% Select Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select subjects
scn_str_col = find(strcmpi( dem_sbj_roi_nme, 'Field Strength'));
sze_sde_col = find(strcmpi( dem_sbj_roi_nme, 'SideOfSeizureFocus'));

ts3_lft = find( cell2mat(dem_sbj_dta(:,scn_str_col))==3 & strcmpi( dem_sbj_dta(:,sze_sde_col), 'left'));
ts3_rgh = find( cell2mat(dem_sbj_dta(:,scn_str_col))==3 & strcmpi( dem_sbj_dta(:,sze_sde_col), 'right'));

lft = find( strcmpi( dem_sbj_dta(:,sze_sde_col), 'left'));
rgh = find( strcmpi( dem_sbj_dta(:,sze_sde_col), 'right'));

% Get the groups


%% Run Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neurobio %%%%%%%%%%%%%%%%%%%%%%%
% Sides %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = neu_bio_dta( ts3_lft, neu_bio_roi_use);
fcfg.lbl_one = neu_bio_roi_nme(neu_bio_roi_use);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta( ts3_lft, neu_bio_roi_use);
fcfg.lbl_two = neu_bio_roi_nme(neu_bio_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' '3T_Sides' '/'];

ejk_cross_cor( fcfg );

% LI %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = neu_bio_lat_dta( ts3_lft, neu_bio_lat_roi_use);
fcfg.lbl_one = neu_bio_lat_roi_nme(neu_bio_lat_roi_use);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_lat_dta( ts3_lft, neu_bio_lat_roi_use);
fcfg.lbl_two = neu_bio_lat_roi_nme(neu_bio_lat_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' '3T_LI' '/'];

ejk_cross_cor( fcfg );

%% Left 3T %%%%%%%%%%%%%%%%%%%%%%%
% Pre - Sides %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = mem_sbj_dta( ts3_lft, 1:4);
fcfg.lbl_one = mem_sbj_roi_nme(1:4);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta( ts3_lft, neu_bio_roi_use);
fcfg.lbl_two = neu_bio_roi_nme(neu_bio_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_3T_Pre' '/'];

ejk_cross_cor( fcfg );

% Pre - LI %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = mem_sbj_dta( ts3_lft, 1:4);
fcfg.lbl_one = mem_sbj_roi_nme(1:4);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_lat_dta( ts3_lft, neu_bio_lat_roi_use);
fcfg.lbl_two = neu_bio_lat_roi_nme(neu_bio_lat_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_3T_Pre_LI' '/'];

ejk_cross_cor( fcfg );

% Post - Sides %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = mem_sbj_dta( ts3_lft, 5:8);
fcfg.lbl_one = mem_sbj_roi_nme(5:8);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta( ts3_lft, neu_bio_roi_use);
fcfg.lbl_two = neu_bio_roi_nme(neu_bio_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_3T_Post' '/'];

ejk_cross_cor( fcfg );

% Post - LI %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta_one = mem_sbj_dta( ts3_lft, 5:8);
fcfg.lbl_one = mem_sbj_roi_nme(5:8);

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_lat_dta( ts3_lft, neu_bio_lat_roi_use);
fcfg.lbl_two = neu_bio_lat_roi_nme(neu_bio_lat_roi_use);

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_put '/' 'Left_3T_Post_LI' '/'];

ejk_cross_cor( fcfg );

%% Right 3T %%%%%%%%%%%%%%%%%%%%%%%%

%% Run Decliner/Non-Decliner Comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left Post Sides %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta     = neu_bio_dta( ts3_lft, :);
fcfg.dta_nme = neu_bio_roi_nme;

fcfg.grp     = mem_grp_sbj_dta( ts3_lft, 5:8);
fcfg.grp_nme = mem_grp_sbj_roi_nme(5:8);

fcfg.out_dir = [ out_put '/' 'Left_post_ttest_3T' '/'];

ejk_ttest2_independent( fcfg );

% Left Post LI %%%
fcfg = [];

fcfg.sbj_nme = mem_sbj_sbj_nme( ts3_lft, 1);

fcfg.dta     = neu_bio_lat_dta( ts3_lft, :);
fcfg.dta_nme = neu_bio_lat_roi_nme;

fcfg.grp     = mem_grp_sbj_dta( ts3_lft, 5:8);
fcfg.grp_nme = mem_grp_sbj_roi_nme(5:8);

fcfg.out_dir = [ out_put '/' 'Left_post_ttest_3T_LI' '/'];

ejk_ttest2_independent( fcfg );

[ mem_grp_sbj_dta(ts3_lft,6) num2cell(mem_sbj_dta(ts3_lft,6)) num2cell(neu_bio_lat_dta(ts3_lft,4)) ]

%% Logistic Regression
% Sides %%%
tst_num_hld = [ 5 6 7 8 ];

for iTN = 1:numel(tst_num_hld)
    
    tst_num = tst_num_hld(iTN);
    
    % L-TLE
    lcfg = [];
    
    lcfg.sbj_grp = [ mem_sbj_sbj_nme( ts3_lft, 1)  mem_grp_sbj_dta( ts3_lft, tst_num) ];
    lcfg.lbl_ord = { 'Decline' 'NoDecline' };
    
    lcfg.dta     = [ mem_sbj_sbj_nme( ts3_lft, 1) num2cell(neu_bio_dta( ts3_lft, :)) ];
    lcfg.dta_lbl = { 'sbj_nme' neu_bio_roi_nme{:} };
    
    lcfg.mdl     = { { 'subcort_vol_Left_Hippocampus_Total' } ...
                     { 'subcort_vol_Left_Hippocampus_Total' 'aseg_MD_Left_Hippocampus' 'aseg_FA_Left_Hippocampus' } };
    lcfg.mdl_nme = { 'Left Volume' ...
                     'Left Volume+MD+FA' };
    
    lcfg.nrm_grp = 0;
    
    lcfg.mdl_cmp_plt = { [1 2] };
    lcfg.mld_cmp_col = { { rgb('blue') rgb('purple') } };
    lcfg.mdl_cmp_nme = { [ mem_grp_sbj_roi_nme{tst_num} '_' 'Hippocampus'] };
    
    lcfg.out_dir = [ out_put '/' 'Left_3T_Logistic_Sides' '/' mem_grp_sbj_roi_nme{tst_num}];
    
    mmil_log_reg(lcfg)
    
end

% LI %%%
tst_num_hld = [ 5 6 7 8 ];

for iTN = 1:numel(tst_num_hld)
    
    tst_num = tst_num_hld(iTN);
    
    % L-TLE
    lcfg = [];
    
    lcfg.sbj_grp = [ mem_sbj_sbj_nme( ts3_lft, 1)  mem_grp_sbj_dta( ts3_lft, tst_num) ];
    lcfg.lbl_ord = { 'Decline' 'NoDecline' };
    
    lcfg.dta     = [ mem_sbj_sbj_nme( ts3_lft, 1) num2cell(neu_bio_lat_dta( ts3_lft, :)) ];
    lcfg.dta_lbl = { 'sbj_nme' neu_bio_lat_roi_nme{:} };
    
    lcfg.mdl     = { { 'LI_hippocampalvolume' } ...
                     { 'LI_hippocampalvolume' 'aseg_MD_LI_Hippocampus' } };
    lcfg.mdl_nme = { 'Volume' ...
                     'Volume+MD' };
    
    lcfg.nrm_grp = 0;
    
    lcfg.mdl_cmp_plt = { [1 2 ] };
    lcfg.mld_cmp_col = { { rgb('blue') rgb('purple') } };
    lcfg.mdl_cmp_nme = { [ mem_grp_sbj_roi_nme{tst_num} '_' 'Hippocampus'] };
    
    lcfg.out_dir = [ out_put '/' 'Left_3T_Logistic_LI' '/' mem_grp_sbj_roi_nme{tst_num}];
    
    mmil_log_reg(lcfg)
    
end
