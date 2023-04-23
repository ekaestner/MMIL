
%%
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%%
pst_var = 'lm2_chg';
dti_var = 'L_Unc_fiber_FA';
mri_var = 'Left_Hippocampus_subcort_vol_ICVcor';

%%
pst_var_col = strcmpi(tot_col,pst_var); 
dti_var_col = strcmpi(tot_col,dti_var); 
mri_var_col = strcmpi(tot_col,mri_var); 

sbj_non_pst = logical(cell2mat(isnan(tot_dta(:,pst_var_col))));
sbj_non_dti = logical(cell2mat(isnan(tot_dta(:,dti_var_col))));
sbj_non_mri = logical(cell2mat(isnan(tot_dta(:,mri_var_col))));

%% Classify missing data patients
pst_dti_sbj = tot_sbj(~sbj_non_pst & ~sbj_non_dti & ~sbj_non_mri);

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