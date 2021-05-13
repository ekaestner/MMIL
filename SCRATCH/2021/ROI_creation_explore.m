
% Left Fusiform = 3007
dta_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_fc082_fmri_160923_20160923.152947_1/';

wmp_prc_fsr = ejk_load_vol([ dta_dir '/' 'mri' '/' 'wmparc.mgz']);

%%
dta_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti/DTIPROC_fc082_fmri_160923_20160923.152947_1/';

wmp_prc = ejk_load_vol([ dta_dir '/' 'DTanalysis' '/' 'wmparc_resDTI.mgz']);

wmp_prc_mat = load([ dta_dir '/' 'DTanalysis' '/' 'wmparc_resDTI_info.mat']);

wmp_dta = ejk_load_vol([ dta_dir '/' 'DTcalc' '/' 'DTI1_crev_corr_regT1_DT_FA.mgz']);

fus_hld = nanmean(wmp_dta(wmp_prc==3007));

%%
% Bring over fsaverage, bring over fc082
% Load fsaverage Desikan, add new fusiform, replace everything except STG/Precentral/LateralOccipital
% Convert to wmparc
% Resample to DTI space


% vector -TO- new.annot (fsaverage) %%%%%%%%%%%%%%%%%%%%%%%
% Save annot with new fusiform & old STG/Precentral/LateralOccipital
% fs_read_annotation
% fs_write_annotation

% new.annot (fsaverage) -TO- new.annot (native subject) %%%%%%%%%%%%%%%%%%%%%%%
% fs_annot2annot


% new.annot (native subject) -TO- new+aseg.mgz -TO- new_wmparc.mgz %%%%%%%%%%%%%%%%%%%%%%%
% mmil_fparc2wmparc OR recon-all dev table?




% new_wmparc.mgz -TO- new_wmparc_resDTI.mgz %%%%%%%%%%%%%%%%%%%%%%%
% MMIL_Analyze_DTI_Exam
nvox_ref = parms.volsz_DTI;
M_ref = parms.M_DTI;
M_reg = inv(parms.M_T1_to_DTI);

[vol,M] = fs_load_mgh(fname,[],fnum);
[vol,M] = mmil_resample_vol(vol,M,...
                            'nvox_ref',nvox_ref,'M_ref',M_ref,'M_reg',M_reg);
fs_save_mgh(vol,fname_out,M);







