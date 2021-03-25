
function [sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(cfg)

%% Load
sbj_hld = mmil_readtext(cfg.red_fle);
for iC = 1:size(sbj_hld,2)
    if isstr(sbj_hld{2,iC})
        sbj_hld(cellfun(@isempty,sbj_hld(:,iC)),iC) = repmat({''},sum(cellfun(@isempty,sbj_hld(:,iC))),1);
    else
        sbj_hld(cellfun(@isempty,sbj_hld(:,iC)),iC) = repmat({nan},sum(cellfun(@isempty,sbj_hld(:,iC))),1);
    end    
end

%% Fix Surgery Type
srg_col = ismember( sbj_hld(1,:), [strcat( 'surgery_name___', cellfun(@num2str, num2cell(1:14),'uni',0)) ['surgery_name___998'] ['surgery_name___999'] ]);
for iS = 2:size(sbj_hld,1)
    srg_hld{iS,1} = find(cell2mat(sbj_hld(iS, srg_col))); 
    if isempty(srg_hld{iS,1}); srg_hld{iS,1}=NaN; end
end 
sbj_hld(:,srg_col) = [];
sbj_hld = [ sbj_hld [ 'surgery_name' ; srg_hld(2:end,:) ] ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Subject Demographics
% sbj_dem
sbj_dem.sbj_nme     = cell(size(sbj_hld,1)-1,1); 
    sbj_nme_col = 'study_id'; sbj_nme_col = strcmpi(sbj_hld(1,:),sbj_nme_col);
sbj_dem.sbj_scn_dte = cell(size(sbj_hld,1)-1,1); 
    sbj_scn_dte_col = 't1_date'; sbj_scn_dte_col = strcmpi(sbj_hld(1,:),sbj_scn_dte_col); % 1 = EPD, 3 = HC
sbj_dem.sbj_grp     = nan(size(sbj_hld,1)-1,1);      
    sbj_grp_col =  'Group'; sbj_grp_col = strcmpi(sbj_hld(1,:),sbj_grp_col); % 1 = EPD, 3 = HC
sbj_dem.sbj_inc     = nan(size(sbj_hld,1)-1,1); 
    sbj_inc_col =  'inclusion'; sbj_inc_col = strcmpi(sbj_hld(1,:),sbj_inc_col);
sbj_dem.sbj_sex     = cell(size(sbj_hld,1)-1,1); 
    sbj_sex_col =  'sex'; sbj_sex_col = strcmpi(sbj_hld(1,:),sbj_sex_col); % 0 = Woman, 1 = Man
sbj_dem.sbj_hnd     = cell(size(sbj_hld,1)-1,1);      
    sbj_hnd_col =  'hand'; sbj_hnd_col = strcmpi(sbj_hld(1,:),sbj_hnd_col); % 7 = right, 8 = left
sbj_dem.sbj_age     = cell(size(sbj_hld,1)-1,1); 
    sbj_age_col =  'dob'; sbj_age_col = strcmpi(sbj_hld(1,:),sbj_age_col); % calculation needed
sbj_dem.sbj_edu     = nan(size(sbj_hld,1)-1,1);      
    sbj_edu_col =  'educat'; sbj_edu_col = strcmpi(sbj_hld(1,:),sbj_edu_col);

%% Setup Subject Seizure Characteristics
sbj_sze.sbj_nme     = cell(size(sbj_hld,1)-1,1); 
sbj_sze.sbj_age_ons = nan(size(sbj_hld,1)-1,1);      
    sbj_age_ons_col = 'age_epil_onset'; sbj_age_ons_col = strcmpi(sbj_hld(1,:),sbj_age_ons_col);
sbj_sze.sbj_sze_dur = nan(size(sbj_hld,1)-1,1);
sbj_sze.sbj_sze_frq = nan(size(sbj_hld,1)-1,1);      
    sbj_sze_frq_col =  'szfreq'; sbj_sze_frq_col = strcmpi(sbj_hld(1,:),sbj_sze_frq_col);
sbj_sze.sbj_sde_ons = cell(size(sbj_hld,1)-1,1);      
    sbj_sde_ons_col =  'side_of_epil_focus'; sbj_sde_ons_col = strcmpi(sbj_hld(1,:),sbj_sde_ons_col); % 1 = left, 2 = right 
sbj_sze.sbj_sde_mri = nan(size(sbj_hld,1)-1,1);      
    sbj_sde_mri_col =  'mri_lateralization'; sbj_sde_mri_col = strcmpi(sbj_hld(1,:),sbj_sde_mri_col); % 1 = left, 2 = right 
sbj_sze.sbj_mts     = cell(size(sbj_hld,1)-1,1);      
    sbj_mts_col     =  'imaging_mts'; sbj_mts_col = strcmpi(sbj_hld(1,:),sbj_mts_col); % 1 = no, 2 = left MTS, 3 = right MTS
sbj_sze.sbj_sze_eti = nan(size(sbj_hld,1)-1,1);      
    sbj_sze_eti_col =  'epil_etiology'; sbj_sze_eti_col = strcmpi(sbj_hld(1,:),sbj_sze_eti_col); %  10 = MTS
sbj_sze.sbj_feb_sze = nan(size(sbj_hld,1)-1,1);      
    sbj_feb_sze_col =  'febrileszhx'; sbj_feb_sze_col = strcmpi(sbj_hld(1,:),sbj_feb_sze_col);
sbj_sze.sbj_aed_num = nan(size(sbj_hld,1)-1,1);      
    sbj_aed_num_col =  'number_aeds'; sbj_aed_num_col = strcmpi(sbj_hld(1,:),sbj_aed_num_col);
sbj_sze.sbj_aed_typ = cell(size(sbj_hld,1)-1,1);      
    sbj_aed_typ_col = {'currentaed1' 'currentaed1' 'currentaed2' 'currentaed3' 'currentaed4' 'currentaed5'}; sbj_aed_typ_col = ismember(sbj_hld(1,:),sbj_aed_typ_col);% 5 = carbamazepine, 14 = gabapentin, 15 = lamotrigine, 16 = levetiracetam, 24 = oxcarbazepine, 27 = phenobarbital, 30 = phenytoin, 34 = topiramate, 36 = valproate, 38 = zonisamide

%% Setup Subject Scan Information
sbj_scn.sbj_nme     = cell(size(sbj_hld,1)-1,1);
sbj_scn.sbj_mr1 = nan(size(sbj_hld,1)-1,1);      
   sbj_mr1_col  = 't1_acquired'; sbj_mr1_col = strcmpi(sbj_hld(1,:),sbj_mr1_col);
sbj_scn.sbj_mr2 = nan(size(sbj_hld,1)-1,1);      
   sbj_mr2_col  = 't1_acquired_v2'; sbj_mr2_col = strcmpi(sbj_hld(1,:),sbj_mr2_col); 
sbj_scn.sbj_dt1 = nan(size(sbj_hld,1)-1,1);      
   sbj_dt1_col  = 't1_acquired'; sbj_dt1_col = strcmpi(sbj_hld(1,:),sbj_dt1_col);
sbj_scn.sbj_dt2 = nan(size(sbj_hld,1)-1,1);      
   sbj_dt2_col  = 'dti_acquired'; sbj_dt2_col = strcmpi(sbj_hld(1,:),sbj_dt2_col);   
sbj_scn.sbj_frm_b11 = nan(size(sbj_hld,1)-1,1);
   sbj_frm_b11_col = 'fmri_block1_acquired'; sbj_frm_b11_col = strcmpi(sbj_hld(1,:),sbj_frm_b11_col);
sbj_scn.sbj_frm_b22 = nan(size(sbj_hld,1)-1,1);      
   sbj_frm_b22_col = 'fmri_block2_acquired_v2'; sbj_frm_b22_col = strcmpi(sbj_hld(1,:),sbj_frm_b22_col); 
sbj_scn.sbj_frm_e11 = nan(size(sbj_hld,1)-1,1);
   sbj_frm_e11_col = 'fmri_event1_acquired'; sbj_frm_e11_col = strcmpi(sbj_hld(1,:),sbj_frm_e11_col);
sbj_scn.sbj_frm_e22 = nan(size(sbj_hld,1)-1,1);      
   sbj_frm_e22_col = 'fmri_event2_acquired_v2'; sbj_frm_e22_col = strcmpi(sbj_hld(1,:),sbj_frm_e22_col); 
sbj_scn.sbj_frm_rs1 = nan(size(sbj_hld,1)-1,1);
   sbj_frm_rs1_col = 'fmri_rs_acquired'; sbj_frm_rs1_col = strcmpi(sbj_hld(1,:),sbj_frm_rs1_col);
sbj_scn.sbj_frm_rs2 = nan(size(sbj_hld,1)-1,1);      
   sbj_frm_rs2_col = 'fmri_rs_acquired_v2'; sbj_frm_rs2_col = strcmpi(sbj_hld(1,:),sbj_frm_rs2_col);  
sbj_scn.sbj_frm_t21 = nan(size(sbj_hld,1)-1,1);
   sbj_frm_t21_col = 't2_acquired'; sbj_frm_t21_col = strcmpi(sbj_hld(1,:),sbj_frm_t21_col);
sbj_scn.sbj_frm_t22 = nan(size(sbj_hld,1)-1,1);      
   sbj_frm_t22_col = 't2_acquired'; sbj_frm_t22_col = strcmpi(sbj_hld(1,:),sbj_frm_t22_col);
sbj_scn.sbj_frm_ri1 = nan(size(sbj_hld,1)-1,1);
   sbj_frm_ri1_col = 't2_acquired'; sbj_frm_ri1_col = strcmpi(sbj_hld(1,:),sbj_frm_ri1_col);
sbj_scn.sbj_frm_ri2 = nan(size(sbj_hld,1)-1,1);      
   sbj_frm_ri2_col = 't2_acquired'; sbj_frm_ri2_col = strcmpi(sbj_hld(1,:),sbj_frm_ri2_col);
   
%% Setup Subject Cognitive Scores
sbj_cog.sbj_nme         = cell(size(sbj_hld,1)-1,1); 

sbj_cog.log_mem_raw_scr_one     = nan(size(sbj_hld,1)-1,1); 
    log_mem_scr_one_raw_col     = 'lm1_recall_raw'; log_mem_scr_one_raw_col = strcmpi(sbj_hld(1,:),log_mem_scr_one_raw_col);  
sbj_cog.log_mem_nor_scr_one     = nan(size(sbj_hld,1)-1,1); 
    log_mem_scr_one_nor_col     = 'lm1_recall_ss'; log_mem_scr_one_nor_col = strcmpi(sbj_hld(1,:),log_mem_scr_one_nor_col);     
sbj_cog.log_mem_raw_scr_one_pst = nan(size(sbj_hld,1)-1,1); 
    log_mem_scr_one_pst_raw_col = 'post_lm1_recall_raw'; log_mem_scr_one_pst_raw_col = strcmpi(sbj_hld(1,:),log_mem_scr_one_pst_raw_col);
sbj_cog.log_mem_nor_scr_one_pst = nan(size(sbj_hld,1)-1,1); 
    log_mem_scr_one_pst_nor_col = 'post_lm1_recall_ss'; log_mem_scr_one_pst_nor_col = strcmpi(sbj_hld(1,:),log_mem_scr_one_pst_nor_col);    
sbj_cog.log_mem_raw_scr_two     = nan(size(sbj_hld,1)-1,1); 
    log_mem_raw_scr_two_col     = 'lm2_raw'; log_mem_raw_scr_two_col = strcmpi(sbj_hld(1,:),log_mem_raw_scr_two_col);
sbj_cog.log_mem_nor_scr_two     = nan(size(sbj_hld,1)-1,1); 
    log_mem_nor_scr_two_col     = 'lm2_ss'; log_mem_nor_scr_two_col = strcmpi(sbj_hld(1,:),log_mem_nor_scr_two_col);    
sbj_cog.log_mem_raw_scr_two_pst = nan(size(sbj_hld,1)-1,1); 
    log_mem_raw_scr_two_pst_col = 'post_lm2_raw'; log_mem_raw_scr_two_pst_col = strcmpi(sbj_hld(1,:),log_mem_raw_scr_two_pst_col);  
sbj_cog.log_mem_nor_scr_two_pst = nan(size(sbj_hld,1)-1,1); 
    log_mem_nor_scr_two_pst_col = 'post_lm2_ss'; log_mem_nor_scr_two_pst_col = strcmpi(sbj_hld(1,:),log_mem_nor_scr_two_pst_col);      
sbj_cog.cvl_lfr_raw_scr = nan(size(sbj_hld,1)-1,1); 
    cvl_lfr_raw_scr_col =  'cvlt_ldfr_raw'; cvl_lfr_raw_scr_col = strcmpi(sbj_hld(1,:),cvl_lfr_raw_scr_col);
sbj_cog.cvl_lfr_nor_scr = nan(size(sbj_hld,1)-1,1); 
    cvl_lfr_nor_scr_col =  'cvlt_ldfr_z'; cvl_lfr_nor_scr_col = strcmpi(sbj_hld(1,:),cvl_lfr_nor_scr_col);    
sbj_cog.cvl_lfr_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cvl_lfr_raw_scr_pst_col =  'post_cvlt_ldfr_raw'; cvl_lfr_raw_scr_pst_col = strcmpi(sbj_hld(1,:),cvl_lfr_raw_scr_pst_col);
sbj_cog.cvl_lfr_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cvl_lfr_nor_scr_pst_col =  'post_cvlt_ldfr_z'; cvl_lfr_nor_scr_pst_col = strcmpi(sbj_hld(1,:),cvl_lfr_nor_scr_pst_col);    
sbj_cog.vp1_raw_scr         = nan(size(sbj_hld,1)-1,1); 
    vp1_raw_scr_col =  'vpa1_raw'; vp1_raw_scr_col = strcmpi(sbj_hld(1,:),vp1_raw_scr_col);
sbj_cog.vp1_nor_scr         = nan(size(sbj_hld,1)-1,1); 
    vp1_nor_scr_col =  'vpa1_ss'; vp1_nor_scr_col = strcmpi(sbj_hld(1,:),vp1_nor_scr_col);    
sbj_cog.vp1_raw_scr_pst     = nan(size(sbj_hld,1)-1,1); 
    vp1_raw_scr_pst_col =  'post_vpa1_raw'; vp1_raw_scr_pst_col = strcmpi(sbj_hld(1,:),vp1_raw_scr_pst_col);
sbj_cog.vp1_nor_scr_pst     = nan(size(sbj_hld,1)-1,1); 
    vp1_nor_scr_pst_col =  'post_vpa1_ss'; vp1_nor_scr_pst_col = strcmpi(sbj_hld(1,:),vp1_nor_scr_pst_col);    
sbj_cog.vp2_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    vp2_raw_scr_col     =  'vpa2_raw'; vp2_raw_scr_col = strcmpi(sbj_hld(1,:),vp2_raw_scr_col);
sbj_cog.vp2_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    vp2_nor_scr_col     =  'vpa2_ss'; vp2_nor_scr_col = strcmpi(sbj_hld(1,:),vp2_nor_scr_col);    
sbj_cog.vp2_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    vp2_raw_scr_pst_col =  'post_vpa2_raw'; vp2_raw_scr_pst_col = strcmpi(sbj_hld(1,:),vp2_raw_scr_pst_col);
sbj_cog.vp2_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    vp2_nor_scr_pst_col =  'post_vpa2_ss'; vp2_nor_scr_pst_col = strcmpi(sbj_hld(1,:),vp2_nor_scr_pst_col);    
sbj_cog.bnt_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    bnt_raw_scr_col     = 'bnt_correct'; bnt_raw_scr_col = strcmpi(sbj_hld(1,:),bnt_raw_scr_col);
sbj_cog.bnt_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    bnt_nor_scr_col     = 'bnt_tscore'; bnt_nor_scr_col = strcmpi(sbj_hld(1,:),bnt_nor_scr_col);    
sbj_cog.bnt_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    bnt_raw_scr_pst_col =  'post_bnt_correct'; bnt_raw_scr_pst_col = strcmpi(sbj_hld(1,:),bnt_raw_scr_pst_col);
sbj_cog.bnt_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    bnt_nor_scr_pst_col =  'post_bnt_tscore'; bnt_nor_scr_pst_col = strcmpi(sbj_hld(1,:),bnt_nor_scr_pst_col);    
sbj_cog.ant_mem_raw_scr = nan(size(sbj_hld,1)-1,1); 
    ant_mem_raw_scr_col =  'ant_total_correct'; ant_mem_raw_scr_col = strcmpi(sbj_hld(1,:),ant_mem_raw_scr_col);   
sbj_cog.ant_mem_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    ant_mem_raw_scr_pst_col =  'post_ant_total_correct'; ant_mem_raw_scr_pst_col = strcmpi(sbj_hld(1,:),ant_mem_raw_scr_pst_col);  
sbj_cog.cat_flu_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    cat_flu_raw_scr_col     =  'category_fluency_total'; cat_flu_raw_scr_col = strcmpi(sbj_hld(1,:),cat_flu_raw_scr_col);
sbj_cog.cat_flu_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    cat_flu_nor_scr_col     =  'category_fluency_ss'; cat_flu_nor_scr_col = strcmpi(sbj_hld(1,:),cat_flu_nor_scr_col);    
sbj_cog.cat_flu_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cat_flu_raw_scr_pst_col =  'post_category_fluency_total'; cat_flu_raw_scr_pst_col = strcmpi(sbj_hld(1,:),cat_flu_raw_scr_pst_col);
sbj_cog.cat_flu_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cat_flu_nor_scr_pst_col =  'post_category_fluency_ss'; cat_flu_nor_scr_pst_col = strcmpi(sbj_hld(1,:),cat_flu_nor_scr_pst_col);    
sbj_cog.cvl_tot_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    cvl_tot_raw_scr_col     =  'cvlt_total_raw'; cvl_tot_raw_scr_col = strcmpi(sbj_hld(1,:),cvl_tot_raw_scr_col);
sbj_cog.cvl_tot_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    cvl_tot_nor_scr_col     =  'cvlt_total_tscore'; cvl_tot_nor_scr_col = strcmpi(sbj_hld(1,:),cvl_tot_nor_scr_col);    
sbj_cog.cvl_tot_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cvl_tot_raw_scr_pst_col =  'post_cvlt_total_raw'; cvl_tot_raw_scr_pst_col = strcmpi(sbj_hld(1,:),cvl_tot_raw_scr_pst_col);
sbj_cog.cvl_tot_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    cvl_tot_nor_scr_pst_col =  'post_cvlt_total_tscore'; cvl_tot_nor_scr_pst_col = strcmpi(sbj_hld(1,:),cvl_tot_nor_scr_pst_col);  
sbj_cog.ltr_tot_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    ltr_tot_raw_scr_col     =  'letter_fluency_total'; ltr_tot_raw_scr_col = strcmpi(sbj_hld(1,:),ltr_tot_raw_scr_col);
sbj_cog.ltr_tot_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    ltr_tot_nor_scr_col     =  'letter_fluency_ss'; ltr_tot_nor_scr_col = strcmpi(sbj_hld(1,:),ltr_tot_nor_scr_col);    
sbj_cog.ltr_tot_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    ltr_tot_raw_scr_pst_col =  'post_letter_fluency_total'; ltr_tot_raw_scr_pst_col = strcmpi(sbj_hld(1,:),ltr_tot_raw_scr_pst_col);
sbj_cog.ltr_tot_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    ltr_tot_nor_scr_pst_col =  'post_letter_fluency_ss'; ltr_tot_nor_scr_pst_col = strcmpi(sbj_hld(1,:),ltr_tot_nor_scr_pst_col);  
sbj_cog.swt_cor_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    swt_cor_raw_scr_col     =  'category_switching_correct'; swt_cor_raw_scr_col = strcmpi(sbj_hld(1,:),swt_cor_raw_scr_col);
sbj_cog.swt_cor_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    swt_cor_nor_scr_col     =  'cs_correct_ss'; swt_cor_nor_scr_col = strcmpi(sbj_hld(1,:),swt_cor_nor_scr_col);    
sbj_cog.swt_cor_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    swt_cor_raw_scr_pst_col =  'post_category_switching_correct'; swt_cor_raw_scr_pst_col = strcmpi(sbj_hld(1,:),swt_cor_raw_scr_pst_col);
sbj_cog.swt_cor_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    swt_cor_nor_scr_pst_col =  'post_cs_correct_ss'; swt_cor_nor_scr_pst_col = strcmpi(sbj_hld(1,:),swt_cor_nor_scr_pst_col);  
sbj_cog.swt_acc_raw_scr     = nan(size(sbj_hld,1)-1,1); 
    swt_acc_raw_scr_col     =  'switching_accuracy'; swt_acc_raw_scr_col = strcmpi(sbj_hld(1,:),swt_acc_raw_scr_col);
sbj_cog.swt_acc_nor_scr     = nan(size(sbj_hld,1)-1,1); 
    swt_acc_nor_scr_col     =  'switching_accuracy_ss'; swt_acc_nor_scr_col = strcmpi(sbj_hld(1,:),swt_acc_nor_scr_col);    
sbj_cog.swt_acc_raw_scr_pst = nan(size(sbj_hld,1)-1,1); 
    swt_acc_raw_scr_pst_col =  'post_switching_accuracy'; swt_acc_raw_scr_pst_col = strcmpi(sbj_hld(1,:),swt_acc_raw_scr_pst_col);
sbj_cog.swt_acc_nor_scr_pst = nan(size(sbj_hld,1)-1,1); 
    swt_acc_nor_scr_pst_col =  'post_switching_accuracy_ss'; swt_acc_nor_scr_pst_col = strcmpi(sbj_hld(1,:),swt_acc_nor_scr_pst_col);
   
%% Setup Subject Emotional Scores
sbj_emo.sbj_nme         = cell(size(sbj_hld,1)-1,1); 

sbj_emo.bdi     = nan(size(sbj_hld,1)-1,1); 
    bdi_col     = 'bdi'; bdi_col = strcmpi(sbj_hld(1,:),bdi_col); 
sbj_emo.bai     = nan(size(sbj_hld,1)-1,1); 
    bai_col     = 'beck_anxiety_inventory'; bai_col = strcmpi(sbj_hld(1,:),bai_col); 
sbj_emo.ndi     = nan(size(sbj_hld,1)-1,1); 
    ndi_col     = 'nddi_e'; ndi_col = strcmpi(sbj_hld(1,:),ndi_col); 
    
%% Setup Subject Surgical Scores    
sbj_srg.sbj_nme = cell(size(sbj_hld,1)-1,1); 

sbj_srg.srg_prf     = nan(size(sbj_hld,1)-1,1); 
    srg_prf_col     = 'surgery_performed'; srg_prf_col = strcmpi(sbj_hld(1,:),srg_prf_col); 
sbj_srg.srg_age     = nan(size(sbj_hld,1)-1,1); 
    srg_age_col     = 'age_surgery'; srg_age_col = strcmpi(sbj_hld(1,:),srg_age_col); 
sbj_srg.srg_sde     = cell(size(sbj_hld,1)-1,1); 
    srg_sde_col     = 'resection_side'; srg_sde_col = strcmpi(sbj_hld(1,:),srg_sde_col); 
sbj_srg.srg_typ     = cell(size(sbj_hld,1)-1,1); 
    srg_typ_col     = 'surgery_name'; srg_typ_col = strcmpi(sbj_hld(1,:),srg_typ_col);  
sbj_srg.eng_out     = cell(size(sbj_hld,1)-1,1); 
    eng_out_col     = 'y1_simple_engel'; eng_out_col = strcmpi(sbj_hld(1,:),eng_out_col);  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subject Demographics
for iS = 2:size(sbj_hld,1)
    sbj_dem.sbj_nme{iS-1,1}     = sbj_hld{iS,sbj_nme_col};
    sbj_dem.sbj_scn_dte{iS-1,1} = sbj_hld{iS,sbj_scn_dte_col};
    sbj_dem.sbj_grp(iS-1,1)     = sbj_hld{iS,sbj_grp_col};
    sbj_dem.sbj_inc(iS-1,1)     = sbj_hld{iS,sbj_inc_col};
    sbj_dem.sbj_sex{iS-1,1}     = sbj_hld{iS,sbj_sex_col}; % 0 = Woman, 1 = Man
    sbj_dem.sbj_hnd{iS-1,1}     = sbj_hld{iS,sbj_hnd_col}; % 7 = right, 8 = left
    sbj_dem.sbj_age{iS-1,1}     = sbj_hld{iS,sbj_age_col}; % calculation needed
    sbj_dem.sbj_edu(iS-1,1)     = sbj_hld{iS,sbj_edu_col};
end

%% Subject Seizure Information
for iS = 2:size(sbj_hld,1)
    sbj_sze.sbj_nme{iS-1,1}     = sbj_hld{iS,sbj_nme_col};
    sbj_sze.sbj_age_ons(iS-1,1) = sbj_hld{iS,sbj_age_ons_col};
    % sbj_sze.sbj_sze_dur = calculation;
    sbj_sze.sbj_sze_frq(iS-1,1) = sbj_hld{iS,sbj_sze_frq_col};
    sbj_sze.sbj_sde_ons{iS-1,1} = sbj_hld{iS,sbj_sde_ons_col}; % 1 = left, 2 = right
    sbj_sze.sbj_sde_mri(iS-1,1) = sbj_hld{iS,sbj_sde_mri_col}; % 1 = left, 2 = right
    sbj_sze.sbj_mts{iS-1,1}     = sbj_hld{iS,sbj_mts_col}; % 1 = no, 2 = left MTS, 3 = right MTS
    sbj_sze.sbj_sze_eti(iS-1,1) = sbj_hld{iS,sbj_sze_eti_col}; %  10 = MTS
    sbj_sze.sbj_feb_sze(iS-1,1) = sbj_hld{iS,sbj_feb_sze_col};
    sbj_sze.sbj_aed_num(iS-1,1) = sbj_hld{iS,sbj_aed_num_col};
    sbj_sze.sbj_aed_typ{iS-1,1} = [sbj_hld{iS,sbj_aed_typ_col}];% 5 = carbamazepine, 14 = gabapentin, 15 = lamotrigine, 16 = levetiracetam, 24 = oxcarbazepine, 27 = phenobarbital, 30 = phenytoin, 34 = topiramate, 36 = valproate, 38 = zonisamide
end

%% Subject Scan Information
for iS = 2:size(sbj_hld,1)
    sbj_scn.sbj_nme{iS-1,1}     = sbj_hld{iS,sbj_nme_col};
    sbj_scn.sbj_mr1(iS-1,1)     = sbj_hld{iS,sbj_mr1_col};
    sbj_scn.sbj_mr2(iS-1,1)     = sbj_hld{iS,sbj_mr2_col}; 
    sbj_scn.sbj_dt1(iS-1,1)     = sbj_hld{iS,sbj_dt1_col};
    sbj_scn.sbj_dt2(iS-1,1) = sbj_hld{iS,sbj_dt2_col};
    sbj_scn.sbj_frm_b11(iS-1,1) = sbj_hld{iS,sbj_frm_b11_col};
    sbj_scn.sbj_frm_b22(iS-1,1) = sbj_hld{iS,sbj_frm_b22_col}; 
    sbj_scn.sbj_frm_e11(iS-1,1) = sbj_hld{iS,sbj_frm_e11_col};
    sbj_scn.sbj_frm_e22(iS-1,1) = sbj_hld{iS,sbj_frm_e22_col}; 
    sbj_scn.sbj_frm_rs1(iS-1,1) = sbj_hld{iS,sbj_frm_rs1_col};
    sbj_scn.sbj_frm_rs2(iS-1,1) = sbj_hld{iS,sbj_frm_rs2_col};  
    sbj_scn.sbj_frm_t21(iS-1,1) = sbj_hld{iS,sbj_frm_t21_col};
    sbj_scn.sbj_frm_t22(iS-1,1) = sbj_hld{iS,sbj_frm_t22_col};
    sbj_scn.sbj_frm_ri1(iS-1,1) = sbj_hld{iS,sbj_frm_ri1_col};
    sbj_scn.sbj_frm_ri2(iS-1,1) = sbj_hld{iS,sbj_frm_ri2_col};
end

%% Subject Cognitive Information
for iS = 2:size(sbj_hld,1)
    sbj_cog.sbj_nme{iS-1,1}         = sbj_hld{iS,sbj_nme_col};
    sbj_cog.log_mem_raw_scr_one(iS-1,1)     = sbj_hld{iS,log_mem_scr_one_raw_col};
    sbj_cog.log_mem_nor_scr_one(iS-1,1)     = sbj_hld{iS,log_mem_scr_one_nor_col};
    sbj_cog.log_mem_raw_scr_one_pst(iS-1,1) = sbj_hld{iS,log_mem_scr_one_pst_raw_col};
    sbj_cog.log_mem_nor_scr_one_pst(iS-1,1) = sbj_hld{iS,log_mem_scr_one_pst_nor_col};
    sbj_cog.log_mem_raw_scr_two(iS-1,1)     = sbj_hld{iS,log_mem_raw_scr_two_col};
    sbj_cog.log_mem_nor_scr_two(iS-1,1)     = sbj_hld{iS,log_mem_nor_scr_two_col};
    sbj_cog.log_mem_raw_scr_two_pst(iS-1,1) = sbj_hld{iS,log_mem_raw_scr_two_pst_col};
    sbj_cog.log_mem_nor_scr_two_pst(iS-1,1) = sbj_hld{iS,log_mem_nor_scr_two_pst_col};
    sbj_cog.cvl_lfr_raw_scr(iS-1,1)     = sbj_hld{iS,cvl_lfr_raw_scr_col};
    sbj_cog.cvl_lfr_nor_scr(iS-1,1)     = sbj_hld{iS,cvl_lfr_nor_scr_col};
    sbj_cog.cvl_lfr_raw_scr_pst(iS-1,1) = sbj_hld{iS,cvl_lfr_raw_scr_pst_col};
    sbj_cog.cvl_lfr_nor_scr_pst(iS-1,1) = sbj_hld{iS,cvl_lfr_nor_scr_pst_col};
    sbj_cog.vp1_raw_scr(iS-1,1)         = sbj_hld{iS,vp1_raw_scr_col};
    sbj_cog.vp1_nor_scr(iS-1,1)         = sbj_hld{iS,vp1_nor_scr_col};
    sbj_cog.vp1_raw_scr_pst(iS-1,1)     = sbj_hld{iS,vp1_raw_scr_pst_col};
    sbj_cog.vp1_nor_scr_pst(iS-1,1)     = sbj_hld{iS,vp1_nor_scr_pst_col};
    sbj_cog.vp2_raw_scr(iS-1,1)         = sbj_hld{iS,vp2_raw_scr_col};
    sbj_cog.vp2_nor_scr(iS-1,1)         = sbj_hld{iS,vp2_nor_scr_col};
    sbj_cog.vp2_raw_scr_pst(iS-1,1)     = sbj_hld{iS,vp2_raw_scr_pst_col};
    sbj_cog.vp2_nor_scr_pst(iS-1,1)     = sbj_hld{iS,vp2_nor_scr_pst_col};
    sbj_cog.bnt_raw_scr(iS-1,1)         = sbj_hld{iS,bnt_raw_scr_col};
    sbj_cog.bnt_nor_scr(iS-1,1)         = sbj_hld{iS,bnt_nor_scr_col};
    sbj_cog.bnt_raw_scr_pst(iS-1,1)     = sbj_hld{iS,bnt_raw_scr_pst_col};
    sbj_cog.bnt_nor_scr_pst(iS-1,1)     = sbj_hld{iS,bnt_nor_scr_pst_col};
    sbj_cog.ant_mem_raw_scr(iS-1,1)     = sbj_hld{iS,ant_mem_raw_scr_col};
    sbj_cog.ant_mem_raw_scr_pst(iS-1,1) = sbj_hld{iS,ant_mem_raw_scr_pst_col};
    sbj_cog.cat_flu_raw_scr(iS-1,1)     = sbj_hld{iS,cat_flu_raw_scr_col};
    sbj_cog.cat_flu_nor_scr(iS-1,1)     = sbj_hld{iS,cat_flu_nor_scr_col};
    sbj_cog.cat_flu_raw_scr_pst(iS-1,1) = sbj_hld{iS,cat_flu_raw_scr_pst_col};  
    sbj_cog.cat_flu_nor_scr_pst(iS-1,1) = sbj_hld{iS,cat_flu_nor_scr_pst_col};      
    sbj_cog.cvl_tot_raw_scr(iS-1,1)     = sbj_hld{iS,cvl_tot_raw_scr_col};
    sbj_cog.cvl_tot_nor_scr(iS-1,1)     = sbj_hld{iS,cvl_tot_nor_scr_col};
    sbj_cog.cvl_tot_raw_scr_pst(iS-1,1) = sbj_hld{iS,cvl_tot_raw_scr_pst_col};  
    sbj_cog.cvl_tot_nor_scr_pst(iS-1,1) = sbj_hld{iS,cvl_tot_nor_scr_pst_col};  
    sbj_cog.ltr_tot_raw_scr(iS-1,1)     = sbj_hld{iS,ltr_tot_raw_scr_col};
    sbj_cog.ltr_tot_nor_scr(iS-1,1)     = sbj_hld{iS,ltr_tot_nor_scr_col};
    sbj_cog.ltr_tot_raw_scr_pst(iS-1,1) = sbj_hld{iS,ltr_tot_raw_scr_pst_col};  
    sbj_cog.ltr_tot_nor_scr_pst(iS-1,1) = sbj_hld{iS,ltr_tot_nor_scr_pst_col};  
    sbj_cog.swt_cor_raw_scr(iS-1,1)     = sbj_hld{iS,swt_cor_raw_scr_col};
    sbj_cog.swt_cor_nor_scr(iS-1,1)     = sbj_hld{iS,swt_cor_nor_scr_col};
    sbj_cog.swt_cor_raw_scr_pst(iS-1,1) = sbj_hld{iS,swt_cor_raw_scr_pst_col};  
    sbj_cog.swt_cor_nor_scr_pst(iS-1,1) = sbj_hld{iS,swt_cor_nor_scr_pst_col};  
    sbj_cog.swt_acc_raw_scr(iS-1,1)     = sbj_hld{iS,swt_acc_raw_scr_col};
    sbj_cog.swt_acc_nor_scr(iS-1,1)     = sbj_hld{iS,swt_acc_nor_scr_col};
    sbj_cog.swt_acc_raw_scr_pst(iS-1,1) = sbj_hld{iS,swt_acc_raw_scr_pst_col};  
    sbj_cog.swt_acc_nor_scr_pst(iS-1,1) = sbj_hld{iS,swt_acc_nor_scr_pst_col};  
end

%% Subject Emotional
for iS = 2:size(sbj_hld,1)
    sbj_emo.sbj_nme{iS-1,1} = sbj_hld{iS,sbj_nme_col};
    sbj_emo.bdi(iS-1,1)     = sbj_hld{iS,bdi_col};
    sbj_emo.bai(iS-1,1)     = sbj_hld{iS,bai_col};
    sbj_emo.ndi(iS-1,1)     = sbj_hld{iS,ndi_col};
end

%% Subject Surgery
for iS = 2:size(sbj_hld,1)
    sbj_srg.sbj_nme{iS-1,1} = sbj_hld{iS,sbj_nme_col};
    sbj_srg.srg_prf(iS-1,1) = sbj_hld{iS,srg_prf_col};
    sbj_srg.srg_age(iS-1,1) = sbj_hld{iS,srg_age_col};
    sbj_srg.srg_sde{iS-1,1} = sbj_hld{iS,srg_sde_col};
    sbj_srg.srg_typ{iS-1,1} = sbj_hld{iS,srg_typ_col}; 
    sbj_srg.eng_out{iS-1,1} = sbj_hld{iS,eng_out_col};     
end

%% Additioanl Calculations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iS = 1:size(sbj_dem.sbj_nme,1)
    if ( ~isempty(sbj_dem.sbj_scn_dte{iS,1}) && ~strcmpi(sbj_dem.sbj_scn_dte{iS,1},'') ) && ...
       ( ~isempty(sbj_dem.sbj_age{iS,1})     && ~strcmpi(sbj_dem.sbj_age{iS,1},'') )
        
        scn_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_scn_dte{iS,1},'-','split'));
        scn_dte = (scn_dte(1)*365) + (scn_dte(2)*30) + scn_dte(3);
        
        brt_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'-','split'));
        brt_dte = (brt_dte(1)*365) + (brt_dte(2)*30) + brt_dte(3);
        
        sbj_dem.sbj_age{iS,1} = roundsd((scn_dte - brt_dte) / 365,3);
    else
        sbj_dem.sbj_age{iS,1} = nan;
    end
end
    
sbj_dem.sbj_age = cell2mat(sbj_dem.sbj_age);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_sex_hld = {0 'F' ; 1 'M' };
sbj_hnd_hld = {7 'R' ; 8 'L' };

for iS = 1:size(sbj_dem.sbj_nme,1)
    try sbj_dem.sbj_sex{iS,1} = sbj_sex_hld{cell2mat(sbj_sex_hld(:,1))==sbj_dem.sbj_sex{iS,1},2}; catch sbj_dem.sbj_sex{iS,1} = ''; end % 0 = Woman, 1 = Man
    try sbj_dem.sbj_hnd{iS,1} = sbj_hnd_hld{cell2mat(sbj_hnd_hld(:,1))==sbj_dem.sbj_hnd{iS,1},2}; catch sbj_dem.sbj_hnd{iS,1} = ''; end % 7 = right, 8 = left
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_sde_ons_hld = {1 'L'   ; 2 'R' };
sbj_mts_hld     = {1 'N/A' ; 2 'L' ; 3 'R'};

for iS = 1:size(sbj_dem.sbj_nme,1)
    try sbj_sze.sbj_sde_ons{iS,1} = sbj_sde_ons_hld{cell2mat(sbj_sde_ons_hld(:,1))== sbj_sze.sbj_sde_ons{iS,1},2}; catch sbj_sze.sbj_sde_ons{iS,1} = ''; end% 1 = left, 2 = right
    try sbj_sze.sbj_mts{iS,1}     = sbj_mts_hld{cell2mat(sbj_mts_hld(:,1))        == sbj_sze.sbj_mts{iS,1},2};     catch sbj_sze.sbj_mts{iS,1} = ''; end % 1 = no, 2 = left MTS, 3 = right MTS
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srg_typ_hld = {1 'diagnostic' ; 2 'ATL' ; 3 'ATL +' ; 4 'amygdalohippocampectomy' ; 5 'neuroablation (hippocampus & amygdala)' ; 6 'neuroablation (hippocampus & amygdala +)' ; 7 'lesionectomy' ; 8 'lesionectomy +' ; 9 'extratemporal resection' ; 10 'multi-lobar resection' ; 11 'hemispherectomy' ; 12 'neocortical resection only' ; 13 'corpus callosotomy' ; 14 'multiple subpial transections' ; 998 'other' ; 999 'unknown' };
eng_out_hld = {1 'I' ; 2 'II' ; 3 'III' ; 4 'IV' };
srg_sde_hld = {1 'L' ; 2 'R' }; 

for iS = 1:size(sbj_dem.sbj_nme,1)
    try sbj_srg.srg_typ{iS,1} = srg_typ_hld{cell2mat(srg_typ_hld(:,1)) == sbj_srg.srg_typ{iS,1},2}; catch sbj_srg.srg_typ{iS,1} = ''; end
    try sbj_srg.eng_out{iS,1} = eng_out_hld{cell2mat(eng_out_hld(:,1)) == sbj_srg.eng_out{iS,1},2}; catch sbj_srg.eng_out{iS,1} = ''; end
    try sbj_srg.srg_sde{iS,1} = srg_sde_hld{cell2mat(srg_sde_hld(:,1)) == sbj_srg.srg_sde{iS,1},2}; catch sbj_srg.srg_sde{iS,1} = ''; end
end

%% Subject Seizure Information
if isfield(cfg,'sbj_nme')

    sbj_ind = ismember( sbj_dem.sbj_nme , cfg.sbj_nme );    
    sbj_dem_nme = fieldnames(sbj_dem); 
    for iF = 1:numel(sbj_dem_nme); sbj_dem.(sbj_dem_nme{iF}) = sbj_dem.(sbj_dem_nme{iF})(sbj_ind); end
    
    sbj_ind = ismember( sbj_sze.sbj_nme , cfg.sbj_nme );    
    sbj_sze_nme = fieldnames(sbj_sze); 
    for iF = 1:numel(sbj_sze_nme); sbj_sze.(sbj_sze_nme{iF}) = sbj_sze.(sbj_sze_nme{iF})(sbj_ind); end
    
    sbj_ind = ismember( sbj_scn.sbj_nme , cfg.sbj_nme );
    sbj_scn_nme = fieldnames(sbj_scn);
    for iF = 1:numel(sbj_scn_nme); sbj_scn.(sbj_scn_nme{iF}) = sbj_scn.(sbj_scn_nme{iF})(sbj_ind); end
    
    sbj_ind = ismember( sbj_cog.sbj_nme , cfg.sbj_nme );
    sbj_cog_nme = fieldnames(sbj_cog); 
    for iF = 1:numel(sbj_cog_nme); sbj_cog.(sbj_cog_nme{iF}) = sbj_cog.(sbj_cog_nme{iF})(sbj_ind); end
    
    sbj_ind = ismember( sbj_emo.sbj_nme , cfg.sbj_nme );
    sbj_emo_nme = fieldnames(sbj_emo); 
    for iF = 1:numel(sbj_emo_nme); sbj_emo.(sbj_emo_nme{iF}) = sbj_emo.(sbj_emo_nme{iF})(sbj_ind); end
    
    sbj_ind = ismember( sbj_srg.sbj_nme , cfg.sbj_nme );
    sbj_srg_nme = fieldnames(sbj_srg); 
    for iF = 1:numel(sbj_srg_nme); sbj_srg.(sbj_srg_nme{iF}) = sbj_srg.(sbj_srg_nme{iF})(sbj_ind); end
    
end

end


