
%%
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%%
pst_var = 'lm2_chg';
dti_var = 'L_Unc_fiber_FA';
mri_var = 'Left_Hippocampus_subcort_vol_ICVcor';
srg_var = 'srg_typ';
ons_sde_var = 'sde_sze_ons';
srg_sde_var = 'srg_sde';
eng_one_var = 'eng_out_I';
eng_two_var = 'eng_out_I_II';
lng_lat_var = 'lng_lat';

%%
pst_var_col = strcmpi(tot_col,pst_var); 
dti_var_col = strcmpi(tot_col,dti_var); 
mri_var_col = strcmpi(tot_col,mri_var); 
srg_var_col = strcmpi(tot_col,srg_var); 
ons_sde_col = strcmpi(tot_col,ons_sde_var); 
srg_sde_col = strcmpi(tot_col,srg_sde_var); 
eng_one_col = strcmpi(tot_col,eng_one_var); 
eng_two_col = strcmpi(tot_col,eng_two_var); 
lng_lat_col = strcmpi(tot_col,lng_lat_var); 

sbj_non_pst = logical(cell2mat(isnan(tot_dta(:,pst_var_col))));
sbj_non_dti = logical(cell2mat(isnan(tot_dta(:,dti_var_col))));
sbj_non_mri = logical(cell2mat(isnan(tot_dta(:,mri_var_col))));

%% Classify missing data patients
smp_sbj     = ~sbj_non_pst & ~sbj_non_dti & ~sbj_non_mri;
pst_dti_sbj = tot_sbj(smp_sbj);

pst_non_img_sbj = tot_sbj(~sbj_non_pst & sbj_non_dti & sbj_non_mri);
pst_non_dti_sbj = tot_sbj(~sbj_non_pst & sbj_non_dti & ~sbj_non_mri);

non_pst_img_sbj = tot_sbj(sbj_non_pst & ~sbj_non_dti & ~sbj_non_mri);

non_pst_non_img_sbj  = tot_sbj(sbj_non_pst & sbj_non_dti & sbj_non_mri);
non_pst_non_dti_sbj  = tot_sbj(sbj_non_pst & sbj_non_dti & ~sbj_non_mri);

% Put together
out_mss_csv = [ pst_dti_sbj repmat({'Included'},numel(pst_dti_sbj),1) ; ...
                pst_non_dti_sbj repmat({'Missing DTI'},numel(pst_non_dti_sbj),1) ; ...
                pst_non_img_sbj repmat({'Missing All Imaging'},numel(pst_non_img_sbj),1) ; ...                
                non_pst_img_sbj repmat({'Missing Post-op Change'},numel(non_pst_img_sbj),1) ; ...
                non_pst_non_dti_sbj repmat({'Missing Post-op Change & DTI'},numel(non_pst_non_dti_sbj),1) ; ...
                non_pst_non_img_sbj repmat({'Missing Post-op Change & All Imaging'},numel(non_pst_non_img_sbj),1) ];
cell2csv(  [ dta_dir '/' 'Subject_Status.csv' ], out_mss_csv);

%% Create types
% Surgical Subjects
slh_sbj = strcmpi(tot_dta(:,srg_var_col),'SLAH');
atl_sbj = strcmpi(tot_dta(:,srg_var_col),'ATL');
srg_sbj = ~sbj_non_pst & ~sbj_non_dti & ~sbj_non_mri & (slh_sbj | atl_sbj);

% Side based on Onset
lft_ons_sbj = strcmpi(tot_dta(:,ons_sde_col),'L');
rgh_ons_sbj = strcmpi(tot_dta(:,ons_sde_col),'R');

% Side based on Surgery
lft_srg_sbj = strcmpi(tot_dta(:,srg_sde_col),'L');
rgh_srg_sbj = strcmpi(tot_dta(:,srg_sde_col),'R');

% Engel I vs all
eng_one_sbj     = strcmpi(tot_dta(:,eng_one_col),'I');
eng_not_one_sbj = strcmpi(tot_dta(:,eng_one_col),'II+');

% Engel I/II vs all
eng_two_sbj     = strcmpi(tot_dta(:,eng_two_col),'I_II');
eng_not_two_sbj = strcmpi(tot_dta(:,eng_two_col),'III_IV');

% Language Laterality
lng_lft_sbj     = strcmpi(tot_dta(:,lng_lat_col),'L');
lng_aty_one_sbj = strcmpi(tot_dta(:,lng_lat_col),'Atypical');

%% Analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Side: Make 2 groups, all
grp.sde_ana.initial.lft = find(lft_srg_sbj & srg_sbj);
grp.sde_ana.initial.rgh = find(rgh_srg_sbj & srg_sbj);

% Side: Make 2 groups, engel I only
grp.sde_ana.eng_one.lft = find(lft_srg_sbj & srg_sbj & eng_one_sbj);
grp.sde_ana.eng_one.rgh = find(rgh_srg_sbj & srg_sbj & eng_one_sbj);

% Side: Make 2 groups, engel I/II only
grp.sde_ana.eng_two.lft = find(lft_srg_sbj & srg_sbj & eng_two_sbj);
grp.sde_ana.eng_two.rgh = find(rgh_srg_sbj & srg_sbj & eng_two_sbj);

% Side: Make 2 groups, Language Laterality Left only
grp.sde_ana.lng_lft.lft = find(lft_srg_sbj & srg_sbj & lng_lft_sbj);
grp.sde_ana.lng_lft.rgh = find(rgh_srg_sbj & srg_sbj & lng_lft_sbj);

%% Make 4 groups
grp.srg_sde_ana.initial.lft_atl = find(atl_sbj & lft_srg_sbj & srg_sbj);
grp.srg_sde_ana.initial.rgh_atl = find(atl_sbj & rgh_srg_sbj & srg_sbj);
grp.srg_sde_ana.initial.lft_slh = find(slh_sbj & lft_srg_sbj & srg_sbj);
grp.srg_sde_ana.initial.rgh_slh = find(slh_sbj & rgh_srg_sbj & srg_sbj);

% Side: Make 2 groups, engel I only
grp.srg_sde_ana.eng_one.lft_atl = find(atl_sbj & lft_srg_sbj & srg_sbj & eng_one_sbj);
grp.srg_sde_ana.eng_one.rgh_atl = find(atl_sbj & rgh_srg_sbj & srg_sbj & eng_one_sbj);
grp.srg_sde_ana.eng_one.lft_slh = find(slh_sbj & lft_srg_sbj & srg_sbj & eng_one_sbj);
grp.srg_sde_ana.eng_one.rgh_slh = find(slh_sbj & rgh_srg_sbj & srg_sbj & eng_one_sbj);

% Side: Make 2 groups, engel I/II only
grp.srg_sde_ana.eng_two.lft_atl = find(atl_sbj & lft_srg_sbj & srg_sbj & eng_two_sbj);
grp.srg_sde_ana.eng_two.rgh_atl = find(atl_sbj & rgh_srg_sbj & srg_sbj & eng_two_sbj);
grp.srg_sde_ana.eng_two.lft_slh = find(slh_sbj & lft_srg_sbj & srg_sbj & eng_two_sbj);
grp.srg_sde_ana.eng_two.rgh_slh = find(slh_sbj & rgh_srg_sbj & srg_sbj & eng_two_sbj);

% Side: Make 2 groups, Language Laterality Left only
grp.srg_sde_ana.lng_lft.lft_atl = find(atl_sbj & lft_srg_sbj & srg_sbj & lng_lft_sbj);
grp.srg_sde_ana.lng_lft.rgh_atl = find(atl_sbj & rgh_srg_sbj & srg_sbj & lng_lft_sbj);
grp.srg_sde_ana.lng_lft.lft_slh = find(slh_sbj & lft_srg_sbj & srg_sbj & lng_lft_sbj);
grp.srg_sde_ana.lng_lft.rgh_slh = find(slh_sbj & rgh_srg_sbj & srg_sbj & lng_lft_sbj);

%% Plot/Table Investigation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% Save
save([ dta_dir '/' 'group_initial.mat' ], 'grp');

%% Save dataset for use
cell2csv([ dta_dir '/' 'Initial_Sample_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ] , [ 'sbj_nme' tot_col ; tot_sbj(srg_sbj,1) tot_dta(srg_sbj,:)]);

  
  
  
  
  
  
