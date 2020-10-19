function abcd_pcqc_redcap(varargin)
%function abcd_pcqc_redcap(varargin)
%
% Optional input:
%   'indir': incoming directory
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%   'outdir': output directory
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'infix': input file suffix
%     {default = 'merged_pcqcinfo'}
%
% Created:  11/05/16 by Jose Teruel
% Prev Mod: 03/10/17 by Don Hagler
% Last Mod: 09/19/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin);
%check lock files
fname_lck=sprintf('%s/.pcqc_redcap.lck',parms.outdir);
if exist(fname_lck,'file')
  fprintf('%s','lock files exist!.');
  return;
end
%Place lock file
fclose(fopen(fname_lck, 'w'));
while exist(sprintf('%s/.merged_pcqcinfo.lck',parms.outdir),'file')
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;



% input files

%% todo: use struct array instead of table

  fname_info = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.infix);

  pcqcinfo = abcd_load_csv(fname_info);
  pcqcinfo = set_pguidevents(pcqcinfo);
  warning off;
  mergedTable = struct2table(pcqcinfo);
  warning on;
  mergedTable.Properties.VariableNames = fieldnames(pcqcinfo);

toDelete = ~(strcmp(mergedTable.ABCD_Compliant,'Yes') | strcmp(mergedTable.ABCD_Compliant,'No'));
mergedTable(toDelete,:) = [];
toDelete = (strcmp(mergedTable.SeriesType,'Undefined'));
mergedTable(toDelete,:) = [];
toDelete = (strcmp(mergedTable.SeriesType,'Undefined_fMRI'));
mergedTable(toDelete,:) = [];

% changes Yes by 1 and No by 0
c_table = find(strcmpi(mergedTable.Properties.VariableNames,'ABCD_Compliant'));
toChange_1 = strcmp(mergedTable.ABCD_Compliant,'Yes');
mergedTable{toChange_1,c_table} =  num2cell([1]);
toChange_0 = strcmp(mergedTable.ABCD_Compliant,'No');
mergedTable{toChange_0,c_table} =  num2cell([0]);

mergedTable = sortrows(mergedTable,{'StudyDate','StudyTime','SeriesTime'},{'ascend','ascend','ascend'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build output table
redcapTable  = cell2table(cell(0,710), 'VariableNames', {'pGUIDevent','pGUID', 'EventName', 'SiteName' 'LatestDate',...,
    'iqc_T1_1_qc_score', 'iqc_T1_1_qc_notes','iqc_T1_1_pc_score', 'iqc_T1_1_pc_notes', 'iqc_T1_1_complete',...,
    'iqc_T1_2_qc_score', 'iqc_T1_2_qc_notes','iqc_T1_2_pc_score', 'iqc_T1_2_pc_notes', 'iqc_T1_2_complete',...,
    'iqc_T1_3_qc_score', 'iqc_T1_3_qc_notes','iqc_T1_3_pc_score', 'iqc_T1_3_pc_notes', 'iqc_T1_3_complete',...,
    'iqc_T1_total_ser','iqc_T1_total_passPC','iqc_T1_good_ser','iqc_T1_bad_ser','iqc_T1_ser_PC_issues','iqc_T1_ser_QCs','iqc_T1_ser_incomp',...,
    'iqc_T1_dis_QC','iqc_T1_dco_QC','iqc_T1_fa_QC','iqc_T1_gh_QC',...,
    'iqc_T1_ht_QC','iqc_T1_hb_QC','iqc_T1_mo_QC','iqc_T1_rf_QC',...,
    'iqc_T1_sd_QC','iqc_T1_si_QC','iqc_T1_sus_QC','iqc_T1_vco_QC',...,
    'iqc_T1_wr_QC','iqc_T1_slice_QC','iqc_T1_dark_QC','iqc_T1_Moire_QC',...,
    'iqc_T1_recon_QC','iqc_T1_miss_QC','iqc_T1_Other_QC',...,    
    'iqc_T2_1_qc_score', 'iqc_T2_1_qc_notes','iqc_T2_1_pc_score', 'iqc_T2_1_pc_notes', 'iqc_T2_1_complete',...,
    'iqc_T2_2_qc_score', 'iqc_T2_2_qc_notes','iqc_T2_2_pc_score', 'iqc_T2_2_pc_notes', 'iqc_T2_2_complete',...,
    'iqc_T2_3_qc_score', 'iqc_T2_3_qc_notes','iqc_T2_3_pc_score', 'iqc_T2_3_pc_notes', 'iqc_T2_3_complete',...,
    'iqc_T2_total_ser','iqc_T2_total_passPC','iqc_T2_good_ser','iqc_T2_bad_ser','iqc_T2_ser_PC_issues','iqc_T2_ser_QCs','iqc_T2_ser_incomp',...,
    'iqc_T2_dis_QC','iqc_T2_dco_QC','iqc_T2_fa_QC','iqc_T2_gh_QC',...,
    'iqc_T2_ht_QC','iqc_T2_hb_QC','iqc_T2_mo_QC','iqc_T2_rf_QC',...,
    'iqc_T2_sd_QC','iqc_T2_si_QC','iqc_T2_sus_QC','iqc_T2_vco_QC',...,
    'iqc_T2_wr_QC','iqc_T2_slice_QC','iqc_T2_dark_QC','iqc_T2_Moire_QC',...,
    'iqc_T2_recon_QC','iqc_T2_miss_QC','iqc_T2_Other_QC',...,   
    'iqc_dMRI_1_FM_qc_score', 'iqc_dMRI_1_FM_qc_notes','iqc_dMRI_1_FM_pc_score', 'iqc_dMRI_1_FM_pc_notes', 'iqc_dMRI_1_FM_complete',...,
    'iqc_dMRI_1_qc_score', 'iqc_dMRI_1_qc_notes','iqc_dMRI_1_pc_score', 'iqc_dMRI_1_pc_notes', 'iqc_dMRI_1_complete',...,
    'iqc_dMRI_2_FM_qc_score', 'iqc_dMRI_2_FM_qc_notes','iqc_dMRI_2_FM_pc_score', 'iqc_dMRI_2_FM_pc_notes', 'iqc_dMRI_2_FM_complete',...,
    'iqc_dMRI_2_qc_score', 'iqc_dMRI_2_qc_notes','iqc_dMRI_2_pc_score', 'iqc_dMRI_2_pc_notes', 'iqc_dMRI_2_complete',...,
    'iqc_dMRI_3_FM_qc_score', 'iqc_dMRI_3_FM_qc_notes','iqc_dMRI_3_FM_pc_score', 'iqc_dMRI_3_FM_pc_notes', 'iqc_dMRI_3_FM_complete',...,
    'iqc_dMRI_3_qc_score', 'iqc_dMRI_3_qc_notes','iqc_dMRI_3_pc_score', 'iqc_dMRI_3_pc_notes', 'iqc_dMRI_3_complete',...,
    'iqc_dMRI_total_ser','iqc_dMRI_total_passPC','iqc_dMRI_good_ser','iqc_dMRI_bad_ser',...,
    'iqc_dMRI_ser_PC_issues','iqc_dMRI_ser_QCs','iqc_dMRI_ser_incomp',...,
    'iqc_dMRI_FM_wronglocus','iqc_dMRI_FM_PC_issues', 'iqc_dMRI_FM_QCs', 'iqc_dMRI_FM_incomp',...,    
    'iqc_dMRI_dis_QC','iqc_dMRI_dco_QC','iqc_dMRI_fa_QC','iqc_dMRI_gh_QC',...,
    'iqc_dMRI_ht_QC','iqc_dMRI_hb_QC','iqc_dMRI_mo_QC','iqc_dMRI_rf_QC',...,
    'iqc_dMRI_sd_QC','iqc_dMRI_si_QC','iqc_dMRI_sus_QC','iqc_dMRI_vco_QC',...,
    'iqc_dMRI_wr_QC','iqc_dMRI_slice_QC','iqc_dMRI_dark_QC','iqc_dMRI_Moire_QC',...,
    'iqc_dMRI_recon_QC','iqc_dMRI_miss_QC','iqc_dMRI_Other_QC',...,
    'iqc_dMRI_FM_dis_QC','iqc_dMRI_FM_dco_FM_QC','iqc_dMRI_FM_fa_QC','iqc_dMRI_FM_gh_QC',...,
    'iqc_dMRI_FM_ht_QC','iqc_dMRI_FM_hb_QC','iqc_dMRI_FM_mo_QC','iqc_dMRI_FM_rf_QC',...,
    'iqc_dMRI_FM_sd_QC','iqc_dMRI_FM_si_QC','iqc_dMRI_FM_sus_QC','iqc_dMRI_FM_vco_QC',...,
    'iqc_dMRI_FM_wr_QC','iqc_dMRI_FM_slice_QC','iqc_dMRI_FM_dark_QC','iqc_dMRI_FM_Moire_QC',...,
    'iqc_dMRI_FM_recon_QC','iqc_dMRI_FM_miss_QC','iqc_dMRI_FM_Other_QC',...,  
    'iqc_rsfMRI_1_FM_qc_score', 'iqc_rsfMRI_1_FM_qc_notes','iqc_rsfMRI_1_FM_pc_score', 'iqc_rsfMRI_1_FM_pc_notes', 'iqc_rsfMRI_1_FM_complete',...,
    'iqc_rsfMRI_1_qc_score', 'iqc_rsfMRI_1_qc_notes','iqc_rsfMRI_1_pc_score', 'iqc_rsfMRI_1_pc_notes', 'iqc_rsfMRI_1_complete','iqc_rsfMRI_1_sub_02',...,
    'iqc_rsfMRI_2_FM_qc_score', 'iqc_rsfMRI_2_FM_qc_notes','iqc_rsfMRI_2_FM_pc_score', 'iqc_rsfMRI_2_FM_pc_notes', 'iqc_rsfMRI_2_FM_complete',...,
    'iqc_rsfMRI_2_qc_score', 'iqc_rsfMRI_2_qc_notes','iqc_rsfMRI_2_pc_score', 'iqc_rsfMRI_2_pc_notes', 'iqc_rsfMRI_2_complete', 'iqc_rsfMRI_2_sub_02',...,
    'iqc_rsfMRI_3_FM_qc_score', 'iqc_rsfMRI_3_FM_qc_notes','iqc_rsfMRI_3_FM_pc_score', 'iqc_rsfMRI_3_FM_pc_notes', 'iqc_rsfMRI_3_FM_complete',...,
    'iqc_rsfMRI_3_qc_score', 'iqc_rsfMRI_3_qc_notes','iqc_rsfMRI_3_pc_score', 'iqc_rsfMRI_3_pc_notes', 'iqc_rsfMRI_3_complete','iqc_rsfMRI_3_sub_02',...,
    'iqc_rsfMRI_4_FM_qc_score', 'iqc_rsfMRI_4_FM_qc_notes','iqc_rsfMRI_4_FM_pc_score', 'iqc_rsfMRI_4_FM_pc_notes', 'iqc_rsfMRI_4_FM_complete',...,
    'iqc_rsfMRI_4_qc_score', 'iqc_rsfMRI_4_qc_notes','iqc_rsfMRI_4_pc_score', 'iqc_rsfMRI_4_pc_notes', 'iqc_rsfMRI_4_complete','iqc_rsfMRI_4_sub_02',...,
    'iqc_rsfMRI_5_FM_qc_score', 'iqc_rsfMRI_5_FM_qc_notes','iqc_rsfMRI_5_FM_pc_score', 'iqc_rsfMRI_5_FM_pc_notes', 'iqc_rsfMRI_5_FM_complete',...,
    'iqc_rsfMRI_5_qc_score', 'iqc_rsfMRI_5_qc_notes','iqc_rsfMRI_5_pc_score', 'iqc_rsfMRI_5_pc_notes', 'iqc_rsfMRI_5_complete','iqc_rsfMRI_5_sub_02',...,    
    'iqc_rsfMRI_6_FM_qc_score', 'iqc_rsfMRI_6_FM_qc_notes','iqc_rsfMRI_6_FM_pc_score', 'iqc_rsfMRI_6_FM_pc_notes', 'iqc_rsfMRI_6_FM_complete',...,
    'iqc_rsfMRI_6_qc_score', 'iqc_rsfMRI_6_qc_notes','iqc_rsfMRI_6_pc_score', 'iqc_rsfMRI_6_pc_notes', 'iqc_rsfMRI_6_complete','iqc_rsfMRI_6_sub_02',...,    
    'iqc_rsfMRI_7_FM_qc_score', 'iqc_rsfMRI_7_FM_qc_notes','iqc_rsfMRI_7_FM_pc_score', 'iqc_rsfMRI_7_FM_pc_notes', 'iqc_rsfMRI_7_FM_complete',...,
    'iqc_rsfMRI_7_qc_score', 'iqc_rsfMRI_7_qc_notes','iqc_rsfMRI_7_pc_score', 'iqc_rsfMRI_7_pc_notes', 'iqc_rsfMRI_7_complete','iqc_rsfMRI_7_sub_02',...,   
    'iqc_rsfMRI_8_FM_qc_score', 'iqc_rsfMRI_8_FM_qc_notes','iqc_rsfMRI_8_FM_pc_score', 'iqc_rsfMRI_8_FM_pc_notes', 'iqc_rsfMRI_8_FM_complete',...,
    'iqc_rsfMRI_8_qc_score', 'iqc_rsfMRI_8_qc_notes','iqc_rsfMRI_8_pc_score', 'iqc_rsfMRI_8_pc_notes', 'iqc_rsfMRI_8_complete','iqc_rsfMRI_8_sub_02',...,       
    'iqc_rsfMRI_9_FM_qc_score', 'iqc_rsfMRI_9_FM_qc_notes','iqc_rsfMRI_9_FM_pc_score', 'iqc_rsfMRI_9_FM_pc_notes', 'iqc_rsfMRI_9_FM_complete',...,
    'iqc_rsfMRI_9_qc_score', 'iqc_rsfMRI_9_qc_notes','iqc_rsfMRI_9_pc_score', 'iqc_rsfMRI_9_pc_notes', 'iqc_rsfMRI_9_complete','iqc_rsfMRI_9_sub_02',...,      
    'iqc_rsfMRI_10_FM_qc_score', 'iqc_rsfMRI_10_FM_qc_notes','iqc_rsfMRI_10_FM_pc_score', 'iqc_rsfMRI_10_FM_pc_notes', 'iqc_rsfMRI_10_FM_complete',...,
    'iqc_rsfMRI_10_qc_score', 'iqc_rsfMRI_10_qc_notes','iqc_rsfMRI_10_pc_score', 'iqc_rsfMRI_10_pc_notes','iqc_rsfMRI_10_complete','iqc_rsfMRI_10_sub_02', ...,   
    'iqc_rsfMRI_11_FM_qc_score', 'iqc_rsfMRI_11_FM_qc_notes','iqc_rsfMRI_11_FM_pc_score', 'iqc_rsfMRI_11_FM_pc_notes', 'iqc_rsfMRI_11_FM_complete',...,
    'iqc_rsfMRI_11_qc_score', 'iqc_rsfMRI_11_qc_notes','iqc_rsfMRI_11_pc_score', 'iqc_rsfMRI_11_pc_notes', 'iqc_rsfMRI_11_complete','iqc_rsfMRI_11_sub_02',...,
    'iqc_rsfMRI_12_FM_qc_score', 'iqc_rsfMRI_12_FM_qc_notes','iqc_rsfMRI_12_FM_pc_score', 'iqc_rsfMRI_12_FM_pc_notes', 'iqc_rsfMRI_12_FM_complete',...,
    'iqc_rsfMRI_12_qc_score', 'iqc_rsfMRI_12_qc_notes','iqc_rsfMRI_12_pc_score', 'iqc_rsfMRI_12_pc_notes', 'iqc_rsfMRI_12_complete','iqc_rsfMRI_12_sub_02',...,
    'iqc_rsfMRI_total_ser','iqc_rsfMRI_total_sub_02','iqc_rsfMRI_total_passPC','iqc_rsfMRI_good_ser_woFM','iqc_rsfMRI_woFM_sub_02','iqc_rsfMRI_good_ser','iqc_rsfMRI_good_sub_02','iqc_rsfMRI_bad_ser',...,
    'iqc_rsfMRI_ser_PC_issues','iqc_rsfMRI_ser_QCs','iqc_rsfMRI_ser_QC_incomp',...,
    'iqc_rsfMRI_FM_wronglocus','iqc_rsfMRI_FM_PC_issues', 'iqc_rsfMRI_FM_QCs', 'iqc_rsfMRI_FM_incomp',...,        
    'iqc_rsfMRI_dis_QC','iqc_rsfMRI_dco_QC','iqc_rsfMRI_fa_QC','iqc_rsfMRI_gh_QC',...,
    'iqc_rsfMRI_ht_QC','iqc_rsfMRI_hb_QC','iqc_rsfMRI_mo_QC','iqc_rsfMRI_rf_QC',...,
    'iqc_rsfMRI_sd_QC','iqc_rsfMRI_si_QC','iqc_rsfMRI_sus_QC','iqc_rsfMRI_vco_QC',...,
    'iqc_rsfMRI_wr_QC','iqc_rsfMRI_slice_QC','iqc_rsfMRI_dark_QC','iqc_rsfMRI_Moire_QC',...,
    'iqc_rsfMRI_recon_QC','iqc_rsfMRI_miss_QC','iqc_rsfMRI_Other_QC',...,  
    'iqc_rsfMRI_FM_dis_QC','iqc_rsfMRI_FM_dco_QC','iqc_rsfMRI_FM_fa_QC','iqc_rsfMRI_FM_gh_QC',...,
    'iqc_rsfMRI_FM_ht_QC','iqc_rsfMRI_FM_hb_QC','iqc_rsfMRI_FM_mo_QC','iqc_rsfMRI_FM_rf_QC',...,
    'iqc_rsfMRI_FM_sd_QC','iqc_rsfMRI_FM_si_QC','iqc_rsfMRI_FM_sus_QC','iqc_rsfMRI_FM_vco_QC',...,
    'iqc_rsfMRI_FM_wr_QC','iqc_rsfMRI_FM_slice_QC','iqc_rsfMRI_FM_dark_QC','iqc_rsfMRI_FM_Moire_QC',...,
    'iqc_rsfMRI_FM_recon_QC','iqc_rsfMRI_FM_miss_QC','iqc_rsfMRI_FM_Other_QC',...,    
    'iqc_MID_1_FM_qc_score', 'iqc_MID_1_FM_qc_notes','iqc_MID_1_FM_pc_score', 'iqc_MID_1_FM_pc_notes', 'iqc_MID_1_FM_complete',...,
    'iqc_MID_1_qc_score', 'iqc_MID_1_qc_notes','iqc_MID_1_pc_score', 'iqc_MID_1_pc_notes', 'iqc_MID_1_complete', 'iqc_MID_1_sub_02',...,
    'iqc_MID_2_FM_qc_score', 'iqc_MID_2_FM_qc_notes','iqc_MID_2_FM_pc_score', 'iqc_MID_2_FM_pc_notes', 'iqc_MID_2_FM_complete',...,
    'iqc_MID_2_qc_score', 'iqc_MID_2_qc_notes','iqc_MID_2_pc_score', 'iqc_MID_2_pc_notes','iqc_MID_2_complete','iqc_MID_2_sub_02',...,
    'iqc_MID_3_FM_qc_score', 'iqc_MID_3_FM_qc_notes','iqc_MID_3_FM_pc_score', 'iqc_MID_3_FM_pc_notes', 'iqc_MID_3_FM_complete',...,
    'iqc_MID_3_qc_score', 'iqc_MID_3_qc_notes','iqc_MID_3_pc_score', 'iqc_MID_3_pc_notes','iqc_MID_3_complete','iqc_MID_3_sub_02',...,    
    'iqc_MID_4_FM_qc_score', 'iqc_MID_4_FM_qc_notes','iqc_MID_4_FM_pc_score', 'iqc_MID_4_FM_pc_notes', 'iqc_MID_4_FM_complete',...,
    'iqc_MID_4_qc_score', 'iqc_MID_4_qc_notes','iqc_MID_4_pc_score', 'iqc_MID_4_pc_notes','iqc_MID_4_complete','iqc_MID_4_sub_02',...,
    'iqc_MID_5_FM_qc_score', 'iqc_MID_5_FM_qc_notes','iqc_MID_5_FM_pc_score', 'iqc_MID_5_FM_pc_notes', 'iqc_MID_5_FM_complete',...,
    'iqc_MID_5_qc_score', 'iqc_MID_5_qc_notes','iqc_MID_5_pc_score', 'iqc_MID_5_pc_notes','iqc_MID_5_complete','iqc_MID_5_sub_02',...,
    'iqc_MID_6_FM_qc_score', 'iqc_MID_6_FM_qc_notes','iqc_MID_6_FM_pc_score', 'iqc_MID_6_FM_pc_notes', 'iqc_MID_6_FM_complete',...,
    'iqc_MID_6_qc_score', 'iqc_MID_6_qc_notes','iqc_MID_6_pc_score', 'iqc_MID_6_pc_notes','iqc_MID_6_complete','iqc_MID_6_sub_02',...,  
    'iqc_MID_total_ser','iqc_MID_total_sub_02','iqc_MID_total_passPC','iqc_MID_good_ser_woFM','iqc_MID_woFM_sub_02','iqc_MID_good_ser','iqc_MID_good_sub_02','iqc_MID_bad_ser',...,
    'iqc_MID_ser_PC_issues','iqc_MID_ser_QCs','iqc_MID_ser_incomp',...,
    'iqc_MID_FM_wronglocus','iqc_MID_FM_PC_issues', 'iqc_MID_FM_QCs', 'iqc_MID_FM_incomp',...,   
    'iqc_MID_dis_QC','iqc_MID_dco_QC','iqc_MID_fa_QC','iqc_MID_gh_QC',...,
    'iqc_MID_ht_QC','iqc_MID_hb_QC','iqc_MID_mo_QC','iqc_MID_rf_QC',...,
    'iqc_MID_sd_QC','iqc_MID_si_QC','iqc_MID_sus_QC','iqc_MID_vco_QC',...,
    'iqc_MID_wr_QC','iqc_MID_slice_QC','iqc_MID_dark_QC','iqc_MID_Moire_QC',...,
    'iqc_MID_recon_QC','iqc_MID_miss_QC','iqc_MID_Other_QC',...,
    'iqc_MID_FM_dis_QC','iqc_MID_FM_dco_QC','iqc_MID_FM_fa_QC','iqc_MID_FM_gh_QC',...,
    'iqc_MID_FM_ht_QC','iqc_MID_FM_hb_QC','iqc_MID_FM_mo_QC','iqc_MID_FM_rf_QC',...,
    'iqc_MID_FM_sd_QC','iqc_MID_FM_si_QC','iqc_MID_FM_sus_QC','iqc_MID_FM_vco_QC',...,
    'iqc_MID_FM_wr_QC','iqc_MID_FM_slice_QC','iqc_MID_FM_dark_QC','iqc_MID_FM_Moire_QC',...,
    'iqc_MID_FM_recon_QC','iqc_MID_FM_miss_QC','iqc_MID_FM_Other_QC',...,   
    'iqc_SST_1_FM_qc_score', 'iqc_SST_1_FM_qc_notes','iqc_SST_1_FM_pc_score', 'iqc_SST_1_FM_pc_notes', 'iqc_SST_1_FM_complete',...,
    'iqc_SST_1_qc_score', 'iqc_SST_1_qc_notes','iqc_SST_1_pc_score', 'iqc_SST_1_pc_notes','iqc_SST_1_complete','iqc_SST_1_sub_02',...,
    'iqc_SST_2_FM_qc_score', 'iqc_SST_2_FM_qc_notes','iqc_SST_2_FM_pc_score', 'iqc_SST_2_FM_pc_notes', 'iqc_SST_2_FM_complete',...,
    'iqc_SST_2_qc_score', 'iqc_SST_2_qc_notes','iqc_SST_2_pc_score', 'iqc_SST_2_pc_notes','iqc_SST_2_complete','iqc_SST_2_sub_02',...,
    'iqc_SST_3_FM_qc_score', 'iqc_SST_3_FM_qc_notes','iqc_SST_3_FM_pc_score', 'iqc_SST_3_FM_pc_notes', 'iqc_SST_3_FM_complete',...,
    'iqc_SST_3_qc_score', 'iqc_SST_3_qc_notes','iqc_SST_3_pc_score', 'iqc_SST_3_pc_notes','iqc_SST_3_complete','iqc_SST_3_sub_02',...,    
    'iqc_SST_4_FM_qc_score', 'iqc_SST_4_FM_qc_notes','iqc_SST_4_FM_pc_score', 'iqc_SST_4_FM_pc_notes', 'iqc_SST_4_FM_complete', ...,
    'iqc_SST_4_qc_score', 'iqc_SST_4_qc_notes','iqc_SST_4_pc_score', 'iqc_SST_4_pc_notes','iqc_SST_4_complete','iqc_SST_4_sub_02',...,
    'iqc_SST_5_FM_qc_score', 'iqc_SST_5_FM_qc_notes','iqc_SST_5_FM_pc_score', 'iqc_SST_5_FM_pc_notes', 'iqc_SST_5_FM_complete',...,
    'iqc_SST_5_qc_score', 'iqc_SST_5_qc_notes','iqc_SST_5_pc_score', 'iqc_SST_5_pc_notes','iqc_SST_5_complete','iqc_SST_5_sub_02',...,
    'iqc_SST_6_FM_qc_score', 'iqc_SST_6_FM_qc_notes','iqc_SST_6_FM_pc_score', 'iqc_SST_6_FM_pc_notes', 'iqc_SST_6_FM_complete',...,
    'iqc_SST_6_qc_score', 'iqc_SST_6_qc_notes','iqc_SST_6_pc_score', 'iqc_SST_6_pc_notes','iqc_SST_6_complete','iqc_SST_6_sub_02',...,
    'iqc_SST_total_ser','iqc_SST_total_sub_02','iqc_SST_total_passPC','iqc_SST_good_ser_woFM','iqc_SST_woFM_sub_02','iqc_SST_good_ser','iqc_SST_good_sub_02','iqc_SST_bad_ser',...,
    'iqc_SST_ser_PC_issues','iqc_SST_ser_QCs','iqc_SST_ser_incomp',...,
    'iqc_SST_FM_wronglocus','iqc_SST_FM_PC_issues', 'iqc_SST_FM_QCs', 'iqc_SST_FM_incomp',...,        
    'iqc_SST_dis_QC','iqc_SST_dco_QC','iqc_SST_fa_QC','iqc_SST_gh_QC',...,
    'iqc_SST_ht_QC','iqc_SST_hb_QC','iqc_SST_mo_QC','iqc_SST_rf_QC',...,
    'iqc_SST_sd_QC','iqc_SST_si_QC','iqc_SST_sus_QC','iqc_SST_vco_QC',...,
    'iqc_SST_wr_QC','iqc_SST_slice_QC','iqc_SST_dark_QC','iqc_SST_Moire_QC',...,
    'iqc_SST_recon_QC','iqc_SST_miss_QC','iqc_SST_Other_QC',...,
    'iqc_SST_FM_dis_QC','iqc_SST_FM_dco_QC','iqc_SST_FM_fa_QC','iqc_SST_FM_gh_QC',...,
    'iqc_SST_FM_ht_QC','iqc_SST_FM_hb_QC','iqc_SST_FM_mo_QC','iqc_SST_FM_rf_QC',...,
    'iqc_SST_FM_sd_QC','iqc_SST_FM_si_QC','iqc_SST_FM_sus_QC','iqc_SST_FM_vco_QC',...,
    'iqc_SST_FM_wr_QC','iqc_SST_FM_slice_QC','iqc_SST_FM_dark_QC','iqc_SST_FM_Moire_QC',...,
    'iqc_SST_FM_recon_QC','iqc_SST_FM_miss_QC','iqc_SST_FM_Other_QC',...,    
    'iqc_nBack_1_FM_qc_score', 'iqc_nBack_1_FM_qc_notes','iqc_nBack_1_FM_pc_score', 'iqc_nBack_1_FM_pc_notes', 'iqc_nBack_1_FM_complete',...,
    'iqc_nBack_1_qc_score', 'iqc_nBack_1_qc_notes','iqc_nBack_1_pc_score', 'iqc_nBack_1_pc_notes','iqc_nBack_1_complete','iqc_nBack_1_sub_02',...,
    'iqc_nBack_2_FM_qc_score', 'iqc_nBack_2_FM_qc_notes','iqc_nBack_2_FM_pc_score', 'iqc_nBack_2_FM_pc_notes', 'iqc_nBack_2_FM_complete',...,
    'iqc_nBack_2_qc_score', 'iqc_nBack_2_qc_notes','iqc_nBack_2_pc_score', 'iqc_nBack_2_pc_notes','iqc_nBack_2_complete','iqc_nBack_2_sub_02',...,
    'iqc_nBack_3_FM_qc_score', 'iqc_nBack_3_FM_qc_notes','iqc_nBack_3_FM_pc_score', 'iqc_nBack_3_FM_pc_notes', 'iqc_nBack_3_FM_complete',...,
    'iqc_nBack_3_qc_score', 'iqc_nBack_3_qc_notes','iqc_nBack_3_pc_score', 'iqc_nBack_3_pc_notes','iqc_nBack_3_complete','iqc_nBack_3_sub_02',...,    
    'iqc_nBack_4_FM_qc_score', 'iqc_nBack_4_FM_qc_notes','iqc_nBack_4_FM_pc_score', 'iqc_nBack_4_FM_pc_notes', 'iqc_nBack_4_FM_complete',...,
    'iqc_nBack_4_qc_score', 'iqc_nBack_4_qc_notes','iqc_nBack_4_pc_score', 'iqc_nBack_4_pc_notes','iqc_nBack_4_complete','iqc_nBack_4_sub_02',...,
    'iqc_nBack_5_FM_qc_score', 'iqc_nBack_5_FM_qc_notes','iqc_nBack_5_FM_pc_score', 'iqc_nBack_5_FM_pc_notes', 'iqc_nBack_5_FM_complete',...,
    'iqc_nBack_5_qc_score', 'iqc_nBack_5_qc_notes','iqc_nBack_5_pc_score', 'iqc_nBack_5_pc_notes','iqc_nBack_5_complete','iqc_nBack_5_sub_02',...,
    'iqc_nBack_6_FM_qc_score', 'iqc_nBack_6_FM_qc_notes','iqc_nBack_6_FM_pc_score', 'iqc_nBack_6_FM_pc_notes', 'iqc_nBack_6_FM_complete',...,
    'iqc_nBack_6_qc_score', 'iqc_nBack_6_qc_notes','iqc_nBack_6_pc_score', 'iqc_nBack_6_pc_notes','iqc_nBack_6_complete','iqc_nBack_6_sub_02',...,
    'iqc_nBack_total_ser','iqc_nBack_total_sub_02','iqc_nBack_total_passPC','iqc_nBack_good_ser_woFM','iqc_nBack_woFM_sub_02','iqc_nBack_good_ser','iqc_nBack_good_sub_02','iqc_nBack_bad_ser',...,
    'iqc_nBack_ser_PC_issues','iqc_nBack_ser_QCs','iqc_nBack_ser_incomp',...,
    'iqc_nBack_FM_wronglocus','iqc_nBack_FM_PC_issues', 'iqc_nBack_FM_QCs','iqc_nBack_FM_incomp',...,        
    'iqc_nBack_dis_QC','iqc_nBack_dco_QC','iqc_nBack_fa_QC','iqc_nBack_gh_QC',...,
    'iqc_nBack_ht_QC','iqc_nBack_hb_QC','iqc_nBack_mo_QC','iqc_nBack_rf_QC',...,
    'iqc_nBack_sd_QC','iqc_nBack_si_QC','iqc_nBack_sus_QC','iqc_nBack_vco_QC',...,
    'iqc_nBack_wr_QC','iqc_nBack_slice_QC','iqc_nBack_dark_QC','iqc_nBack_Moire_QC',...,
    'iqc_nBack_recon_QC','iqc_nBack_miss_QC','iqc_nBack_Other_QC',...,
    'iqc_nBack_FM_dis_QC','iqc_nBack_FM_dco_QC','iqc_nBack_FM_fa_QC','iqc_nBack_FM_gh_QC',...,
    'iqc_nBack_FM_ht_QC','iqc_nBack_FM_hb_QC','iqc_nBack_FM_mo_QC','iqc_nBack_FM_rf_QC',...,
    'iqc_nBack_FM_sd_QC','iqc_nBack_FM_si_QC','iqc_nBack_FM_sus_QC','iqc_nBack_FM_vco_QC',...,
    'iqc_nBack_FM_wr_QC','iqc_nBack_FM_slice_QC','iqc_nBack_FM_dark_QC','iqc_nBack_FM_Moire_QC',...,
    'iqc_nBack_FM_recon_QC','iqc_nBack_FM_miss_QC','iqc_nBack_FM_Other_QC',...,
    'iqc_EventComplete', 'iqc_EventComplete_C_PC_QC'...,
    });

unique_p_events = unique(mergedTable.pguidevent);
content = cell(size(unique_p_events,1),size(redcapTable,2));
content(:,1)=unique_p_events;
contentTable = cell2table(content);
contentTable.Properties.VariableNames = redcapTable.Properties.VariableNames;
redcapTable = [redcapTable;contentTable];

% fill in table
while ~isempty(mergedTable)
  p_event = mergedTable.pguidevent(1);
  index_redcap = find(strcmp(redcapTable.pGUIDevent,p_event));
  events = strcmp(mergedTable.pguidevent(:),mergedTable.pguidevent(1));
  indices_event = find(events==1);
  eventTable = mergedTable(indices_event,:);
  mergedTable(indices_event,:) = [];
  visitID = eventTable.VisitID{1};
  if visitID(1)=='G'
    redcapTable = checkOneEvent_GE(redcapTable, eventTable, index_redcap);
  else
    redcapTable = checkOneEvent(redcapTable, eventTable, index_redcap);
  end
end

redcapTable = check_completeness(redcapTable);

redcapTable.Properties.VariableNames{'pGUID'} = 'id_redcap';
redcapTable.Properties.VariableNames{'EventName'} = 'redcap_event_name';
redcapTable.Properties.VariableNames = lower(redcapTable.Properties.VariableNames);
outfile = sprintf('%s/%s_REDCapinfo.csv',parms.outdir,parms.outstem);
redcap_struct=table2struct(redcapTable);
redcap_struct = mmil_sortstruct(redcap_struct,{'sitename','id_redcap','redcap_event_name'});
mmil_struct2csv(redcap_struct,outfile);

%delete lock file
delete(fname_lck);


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'infix','merged_pcqcinfo',[],...
    ...
    'series_types',{'T1','T2','dMRI','rsfMRI','MID','SST','nBack'},[],...
  });
  if parms.outdir(1) ~= '/'
     parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir);
  end
  if parms.indir(1) ~= '/', parms.indir=parms.outdir; end;
  if ~exist(parms.indir)
    error('input directory %s not found',parms.indir);
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redcapTable = check_completeness(redcapTable)
  for i=1:height(redcapTable)
    blocks_acquired = ones(1,7);
    blocks_c_pc_qc = ones(1,7);

    %% todo: loop over parms.series_types
    % T1s
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_acquired(1) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_c_pc_qc(1) = 0;
    end

    % T2s
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_acquired(2) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_c_pc_qc(2) = 0;
    end

    % dMRI
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_acquired(3) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>0
        blocks_c_pc_qc(3) = 0;
    end

    % rsfMRI
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>2
        blocks_acquired(4) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>2
        blocks_c_pc_qc(4) = 0;
    end

    % MID
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
       blocks_acquired(5) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
        blocks_c_pc_qc(5) = 0;
    end

    % SST
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
       blocks_acquired(6) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
        blocks_c_pc_qc(6) = 0;
    end

    % nBack
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_total_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
       blocks_acquired(7) = 0;
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_good_ser'));
    if  cell2mat(redcapTable{i,var_index})>1
        blocks_c_pc_qc(7) = 0;
    end

    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_EventComplete'));
    if any(blocks_acquired)
        redcapTable{i,var_index}=num2cell(0);
    else
        redcapTable{i,var_index}=num2cell(1);
    end
    var_index = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_EventComplete_C_PC_QC'));
    if any(blocks_c_pc_qc)
        redcapTable{i,var_index}=num2cell(0);
    else
        redcapTable{i,var_index}=num2cell(1);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redcapTable = checkOneEvent(redcapTable, eventTable, index_redcap)
  redcapTable.pGUID(index_redcap) = eventTable.pGUID(1);
  redcapTable.EventName(index_redcap) = eventTable.EventName(1);
  redcapTable.SiteName(index_redcap) = eventTable.SiteName(1);

  auxTab = sortrows(eventTable,'StudyDate');
  redcapTable.LatestDate(index_redcap) = auxTab.StudyDate(end);

  % look for rsfMRI series

  rsfMRI_events= strcmp(eventTable.SeriesType, 'rsfMRI');
  rsfMRI_index = find(rsfMRI_events==1);
  if ~isempty(rsfMRI_index)
      [redcapTable]  = classify_rsfMRI(eventTable, redcapTable, rsfMRI_index, index_redcap);
  end

  fMRI_MID_events= strcmp(eventTable.SeriesType, 'fMRI_MID_task');
  fMRI_MID_index = find(fMRI_MID_events==1);
  if ~isempty(fMRI_MID_index)
      [redcapTable] = classify_fMRI_MID(eventTable, redcapTable, fMRI_MID_index, index_redcap);
  end

  fMRI_SST_events= strcmp(eventTable.SeriesType, 'fMRI_SST_task');
  fMRI_SST_index = find(fMRI_SST_events==1);
  if ~isempty(fMRI_SST_index)
      [redcapTable] = classify_fMRI_SST(eventTable, redcapTable, fMRI_SST_index, index_redcap);
  end

  fMRI_nBack_events= strcmp(eventTable.SeriesType, 'fMRI_nBack_task');
  fMRI_nBack_index = find(fMRI_nBack_events==1);
  if ~isempty(fMRI_nBack_index)
      [redcapTable] = classify_fMRI_nBack(eventTable, redcapTable, fMRI_nBack_index, index_redcap);
  end

  dMRI_events= strcmp(eventTable.SeriesType, 'dMRI');
  dMRI_index = find(dMRI_events==1);
  if ~isempty(dMRI_index)
      [redcapTable] = classify_dMRI(eventTable, redcapTable, dMRI_index, index_redcap);
  end

  T1_events= strcmp(eventTable.SeriesType, 'T1');
  T1_index = find(T1_events==1);
  if ~isempty(T1_index)
      [redcapTable] = classify_T1(eventTable, redcapTable, T1_index, index_redcap);
  end

  T2_events= strcmp(eventTable.SeriesType, 'T2');
  T2_index = find(T2_events==1);
  if ~isempty(T2_index)
      [redcapTable] = classify_T2(eventTable, redcapTable, T2_index, index_redcap);
  end       
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function redcapTable = checkOneEvent_GE(redcapTable, eventTable, index_redcap)
     
  redcapTable.pGUID(index_redcap) = eventTable.pGUID(1);
  redcapTable.EventName(index_redcap) = eventTable.EventName(1);
  redcapTable.SiteName(index_redcap) = eventTable.SiteName(1);
  auxTab = sortrows(eventTable,'StudyDate');
  redcapTable.LatestDate(index_redcap) = auxTab.StudyDate(end);
  
  % look for rsfMRI series

  rsfMRI_events= strcmp(eventTable.SeriesType, 'rsfMRI');
  rsfMRI_index = find(rsfMRI_events==1);
  if ~isempty(rsfMRI_index)
      [redcapTable]  = classify_rsfMRI_GE(eventTable, redcapTable, rsfMRI_index, index_redcap);
  end

  fMRI_MID_events= strcmp(eventTable.SeriesType, 'fMRI_MID_task');
  fMRI_MID_index = find(fMRI_MID_events==1);
  if ~isempty(fMRI_MID_index)
      [redcapTable] = classify_fMRI_MID_GE(eventTable, redcapTable, fMRI_MID_index, index_redcap);
  end

  fMRI_SST_events= strcmp(eventTable.SeriesType, 'fMRI_SST_task');
  fMRI_SST_index = find(fMRI_SST_events==1);
  if ~isempty(fMRI_SST_index)
      [redcapTable] = classify_fMRI_SST_GE(eventTable, redcapTable, fMRI_SST_index, index_redcap);
  end

  fMRI_nBack_events= strcmp(eventTable.SeriesType, 'fMRI_nBack_task');
  fMRI_nBack_index = find(fMRI_nBack_events==1);
  if ~isempty(fMRI_nBack_index)
      [redcapTable] = classify_fMRI_nBack_GE(eventTable, redcapTable, fMRI_nBack_index, index_redcap);
  end

  dMRI_events= strcmp(eventTable.SeriesType, 'dMRI');
  dMRI_index = find(dMRI_events==1);
  if ~isempty(dMRI_index)
      [redcapTable] = classify_dMRI_GE(eventTable, redcapTable, dMRI_index, index_redcap);
  end

  T1_events= strcmp(eventTable.SeriesType, 'T1');
  T1_index = find(T1_events==1);
  if ~isempty(T1_index)
      [redcapTable] = classify_T1(eventTable, redcapTable, T1_index, index_redcap);
  end

  T2_events= strcmp(eventTable.SeriesType, 'T2');
  T2_index = find(T2_events==1);
  if ~isempty(T2_index)
      [redcapTable] = classify_T2(eventTable, redcapTable, T2_index, index_redcap);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = classify_T1(eventTable, redcapTable, T1_index, index_redcap)
  i_series = 0;
  i_pc = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  %Allow for 3 runs max
  if length(T1_index)>3
     T1_index = T1_index(1:3); 
  end
  i_total = length(T1_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_1_qc_score'));

  % fill issues with 0=
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_dis_QC'));
  index_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_Other_QC'));
  for i=index_issues:index_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(T1_index));
  for index=1:length(T1_index)
      c_table_run=(c_table + (index-1)*5);
      index_run = T1_index(index);  
      % Fill table with T1_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end

      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0;
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if pc, i_pc = i_pc+1; end
      if pc && qc && ic
          i_series = i_series+1;
      elseif qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));
      end

  end
  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+2} = num2cell(i_series);
      redcapTable{index_redcap,c_series+3} = num2cell(i_total-i_series);
  end
  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T1_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = classify_T2(eventTable, redcapTable, T2_index, index_redcap)
  i_series = 0;
  i_pc = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  %Allow for 3 runs max
  if length(T2_index)>3
     T2_index = T2_index(1:3); 
  end
  i_total = length(T2_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_dis_QC'));
  index_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_Other_QC'));
  for i=index_issues:index_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(T2_index));
  for index=1:length(T2_index)
      c_table_run=(c_table + (index-1)*5);
      index_run = T2_index(index);  
      % Fill table with T2_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end

      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end    

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc) 
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if isempty(pc), pc=0; end
      if pc, i_pc = i_pc+1; end
      if pc && qc && ic
          i_series = i_series+1;
      elseif qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end
  end

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+2} = num2cell(i_series);
      redcapTable{index_redcap,c_series+3} = num2cell(i_total-i_series);
  end
  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_T2_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function [redcapTable]  = classify_dMRI(eventTable, redcapTable, dMRI_index, index_redcap)
  i_series=0;
  i_pc = 0;
  flag_fm=0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;

  %Allow for 3 runs max
  if length(dMRI_index)>3
     dMRI_index = dMRI_index(1:3); 
  end
  i_total = length(dMRI_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(dMRI_index));
  for index=1:length(dMRI_index)
      c_table_run=(c_table + (index-1)*10);
      index_run = dMRI_index(index);  
      % Fill table with dMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);

      % Check for possible field map, only inmediately before allowed 
      dFM_fwd_index = index_run-2; if dFM_fwd_index<1, dFM_fwd_index=1; end
      dFM_rev_index = index_run-1; if dFM_rev_index<1, dFM_rev_index=1; end
      if (strcmp(eventTable(dFM_fwd_index,:).SeriesType, 'dMRI_FM_PA') && strcmp(eventTable(dFM_rev_index,:).SeriesType, 'dMRI_FM_AP'))
          flag_fm = 1;
          if ~isempty(eventTable.QC{dFM_fwd_index}) && ~isempty(eventTable.QC{dFM_rev_index})
              redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{dFM_fwd_index} && eventTable.QC{dFM_rev_index});
          end 
          QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(dFM_fwd_index,index_notes(1)))));
          QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(dFM_rev_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(dFM_fwd_index,index_notes(i_notes)))));
              QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(dFM_rev_index,index_notes(i_notes)))));
          end
          QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(dFM_fwd_index),' - R2FM1:', eventTable.notes_TT(dFM_fwd_index),' - R1FM2:',eventTable.notes_CS(dFM_rev_index),' - R2FM2:', eventTable.notes_TT(dFM_rev_index));
          redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{dFM_fwd_index} && eventTable.ABCD_Compliant{dFM_rev_index});
          redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(dFM_fwd_index),eventTable.AdditionalInfo(dFM_rev_index)); 
          redcapTable{index_redcap,c_table_run-1} = num2cell(eventTable.Completed{dFM_fwd_index} && eventTable.Completed{dFM_rev_index});
      end  

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = strcat(eventTable.notes_CS(dFM_fwd_index),eventTable.notes_CS(dFM_rev_index));
      %    fm_issues_tmp_2 = strcat(eventTable.notes_TT(dFM_fwd_index),eventTable.notes_TT(dFM_rev_index));   
      %end

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end  

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0;
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end    

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0;
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);    
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));    
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
      end
  end

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+2} = num2cell(i_series);
      redcapTable{index_redcap,c_series+3} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_dMRI_GE(eventTable, redcapTable, dMRI_index, index_redcap)
  i_series=0;
  i_pc = 0;
  flag_fm=0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 3 runs max
  if length(dMRI_index)>3
     dMRI_index = dMRI_index(1:3); 
  end
  i_total = length(dMRI_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(dMRI_index));
  for index=1:length(dMRI_index)
      c_table_run=(c_table + (index-1)*10);
      index_run = dMRI_index(index);  
      % Fill table with dMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);

      % Check for possible field map, only inmediately before allowed
      dFM_index = index_run-1; if dFM_index<1, dFM_index=1; end 
      if (strcmp(eventTable(dFM_index,:).SeriesType, 'dMRI_FM'))
          flag_fm=1;
          redcapTable{index_redcap,c_table_run-5} = eventTable.QC(dFM_index);
          QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(dFM_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(dFM_index,index_notes(i_notes)))));
          end
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(dFM_index),' - R2:', eventTable.notes_TT(dFM_index));
          redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(dFM_index);
          redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(dFM_index);
          redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(dFM_index);
      end

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run); 
      %if flag_fm
      %    fm_issues_tmp_1 = eventTable.notes_CS(dFM_index);
      %    fm_issues_tmp_2 = eventTable.notes_TT(dFM_index);   
      %end

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end  

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index) = 0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end   

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0;
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);           
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));           
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
      end
  end
  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+2} = num2cell(i_series);
      redcapTable{index_redcap,c_series+3} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_dMRI_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_fMRI_MID(eventTable, redcapTable, fMRI_MID_index, index_redcap)
  i_series=0;
  i_pc = 0 ;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_MID_index)>6
     fMRI_MID_index = fMRI_MID_index(1:6); 
  end
  i_total = length(fMRI_MID_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
    redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_MID_index));
  for index=1:length(fMRI_MID_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_MID_index(index);  
      % Fill table with fMRI_MID_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      %% todo: use notes_CP, and notes_LH -- make reviewer list an input parameter
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_fwd_index = index_run-2; if fFM_fwd_index<1, fFM_fwd_index=1; end
      fFM_rev_index = index_run-1; if fFM_rev_index<1, fFM_rev_index=1; end
      if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
          flag_fm=1;
          if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
              redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
          end 
          QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
          QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
              QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
          end
          QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
          redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
      else
          fFM_fwd_index = index_run-3; if fFM_fwd_index<1, fFM_fwd_index=1; end
          fFM_rev_index = index_run-2; if fFM_rev_index<1, fFM_rev_index=1; end 
          if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
              flag_fm=1;
              if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
                  redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
              end 
              QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
              QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
                  QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
              end
              QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
              redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
          end   
      end


      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = strcat(eventTable.notes_CS(fFM_fwd_index),eventTable.notes_CS(fFM_rev_index));
      %    fm_issues_tmp_2 = strcat(eventTable.notes_TT(fFM_fwd_index),eventTable.notes_TT(fFM_rev_index));   
      %end    
      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end      

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end


      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc) 
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end 
      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm)
              qc_fm=0;
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);               
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));               
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 


  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = classify_fMRI_MID_GE(eventTable, redcapTable, fMRI_MID_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_MID_index)>6
     fMRI_MID_index = fMRI_MID_index(1:6); 
  end
  i_total = length(fMRI_MID_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_MID_index));
  for index=1:length(fMRI_MID_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_MID_index(index);  
      % Fill table with rsfMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=1:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_index = index_run-1; if fFM_index<1, fFM_index=1; end 
      if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
          flag_fm=1;
          redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
          QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
          end
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
          redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
          redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
          redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
      else
          fFM_index = index_run-2; if fFM_index<1, fFM_index=1; end
          if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
              flag_fm=1;
              redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
              QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
              end
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
              redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
              redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
              redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
          end   
      end

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = eventTable.notes_CS(fFM_index);
      %    fm_issues_tmp_2 = eventTable.notes_TT(fFM_index);   
      %end        

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end    

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end    

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0;
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end 

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm)
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);           
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));           
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_MID_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_fMRI_SST(eventTable, redcapTable, fMRI_SST_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_SST_index)>6
     fMRI_SST_index = fMRI_SST_index(1:6); 
  end
  i_total = length(fMRI_SST_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_SST_index));
  for index=1:length(fMRI_SST_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_SST_index(index);  
      % Fill table with fMRI_SST_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_fwd_index = index_run-2; if fFM_fwd_index<1, fFM_fwd_index=1; end
      fFM_rev_index = index_run-1; if fFM_rev_index<1, fFM_rev_index=1; end
      if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
          flag_fm=1;
          if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
              redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
          end 
          QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
          QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
              QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
          end
          QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
          redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
      else
          fFM_fwd_index = index_run-3; if fFM_fwd_index<1, fFM_fwd_index=1; end
          fFM_rev_index = index_run-2; if fFM_rev_index<1, fFM_rev_index=1; end 
          if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
              flag_fm=1;
              if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
                  redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
              end 
              QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
              QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
                  QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
              end
              QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
              redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
          end   
      end


      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = strcat(eventTable.notes_CS(fFM_fwd_index),eventTable.notes_CS(fFM_rev_index));
      %    fm_issues_tmp_2 = strcat(eventTable.notes_TT(fFM_fwd_index),eventTable.notes_TT(fFM_rev_index));   
      %end    
      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end      

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end


      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end     
      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm)
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end        
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));               
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 


  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = classify_fMRI_SST_GE(eventTable, redcapTable, fMRI_SST_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_SST_index)>6
     fMRI_SST_index = fMRI_SST_index(1:6); 
  end
  i_total = length(fMRI_SST_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_SST_index));
  for index=1:length(fMRI_SST_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_SST_index(index);  
      % Fill table with rsfMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);
      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_index = index_run-1; if fFM_index<1, fFM_index=1; end 
      if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
          flag_fm=1;
          redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
          QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
          end
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
          redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
          redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
          redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
      else
          fFM_index = index_run-2; if fFM_index<1, fFM_index=1; end
          if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
              flag_fm=1;
              redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
              QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
              end
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
              redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
              redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
              redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
          end   
      end

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = eventTable.notes_CS(fFM_index);
      %    fm_issues_tmp_2 = eventTable.notes_TT(fFM_index);   
      %end     
      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1;
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end    

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end    

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end 

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm)
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));    
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_SST_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_fMRI_nBack(eventTable, redcapTable, fMRI_nBack_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_nBack_index)>6
     fMRI_nBack_index = fMRI_nBack_index(1:6); 
  end
  i_total = length(fMRI_nBack_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_nBack_index));
  for index=1:length(fMRI_nBack_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_nBack_index(index);  
      % Fill table with fMRI_nBack_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_fwd_index = index_run-2; if fFM_fwd_index<1, fFM_fwd_index=1; end
      fFM_rev_index = index_run-1; if fFM_rev_index<1, fFM_rev_index=1; end
      if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
          flag_fm=1;
          if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
              redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
          end 
          QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
          QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
              QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
          end
          QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
          redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
      else
          fFM_fwd_index = index_run-3; if fFM_fwd_index<1, fFM_fwd_index=1; end
          fFM_rev_index = index_run-2; if fFM_rev_index<1, fFM_rev_index=1; end 
          if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))
              flag_fm=1;
              if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
                  redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
              end 
              QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
              QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
                  QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
              end
              QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
              redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
          end   
      end


      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = strcat(eventTable.notes_CS(fFM_fwd_index),eventTable.notes_CS(fFM_rev_index));
      %    fm_issues_tmp_2 = strcat(eventTable.notes_TT(fFM_fwd_index),eventTable.notes_TT(fFM_rev_index));   
      %end   
      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end      

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end


      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end     
      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);                  
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));                  
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 


  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = classify_fMRI_nBack_GE(eventTable, redcapTable, fMRI_nBack_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 6 maximum runs
  if length(fMRI_nBack_index)>6
     fMRI_nBack_index = fMRI_nBack_index(1:6); 
  end
  i_total = length(fMRI_nBack_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(fMRI_nBack_index));
  for index=1:length(fMRI_nBack_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = fMRI_nBack_index(index);  
      % Fill table with rsfMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);
      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_index = index_run-1; if fFM_index<1, fFM_index=1; end 
      if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
          flag_fm=1;
          redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
          QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
          end
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
          redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
          redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
          redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
      else
          fFM_index = index_run-2; if fFM_index<1, fFM_index=1; end
          if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
              flag_fm=1;
              redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);
              QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
              end
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
              redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
              redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
              redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
          end   
      end


      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = eventTable.notes_CS(fFM_index);
      %    fm_issues_tmp_2 = eventTable.notes_TT(fFM_index);   
      %end    

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end    

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end    

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end     

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));
          end

          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_nBack_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_rsfMRI(eventTable, redcapTable, rsfMRI_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 12 maximum runs
  if length(rsfMRI_index)>12
     rsfMRI_index = rsfMRI_index(1:12); 
  end
  i_total = length(rsfMRI_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(rsfMRI_index));
  for index=1:length(rsfMRI_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = rsfMRI_index(index);  
      % Fill table with fMRI_rsfMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);

      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_fwd_index = index_run-2; if fFM_fwd_index<1, fFM_fwd_index=1; end
      fFM_rev_index = index_run-1; if fFM_rev_index<1, fFM_rev_index=1; end
      if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))...,
              || (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_AP') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_PA'))
          flag_fm=1;
          if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
              redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
          end 
          QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
          QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
              QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
          end
          QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
          redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
          redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
      else
          fFM_fwd_index = index_run-3; if fFM_fwd_index<1, fFM_fwd_index=1; end
          fFM_rev_index = index_run-2; if fFM_rev_index<1, fFM_rev_index=1; end 
          if (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_PA') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_AP'))...,
                  || (strcmp(eventTable(fFM_fwd_index,:).SeriesType, 'fMRI_FM_AP') && strcmp(eventTable(fFM_rev_index,:).SeriesType, 'fMRI_FM_PA'))
              flag_fm=1;
              if ~isempty(eventTable.QC{fFM_fwd_index}) && ~isempty(eventTable.QC{fFM_rev_index})
                  redcapTable{index_redcap,c_table_run-5} = num2cell(eventTable.QC{fFM_fwd_index} && eventTable.QC{fFM_rev_index});
              end 
              QC_notes_FM1 = strcat('R1FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(1)))));
              QC_notes_FM2 = strcat('R1FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM1 = strcat(QC_notes_FM1,' - R',num2str(i_notes),'FM1:',cell2mat(table2cell(eventTable(fFM_fwd_index,index_notes(i_notes)))));
                  QC_notes_FM2 = strcat(QC_notes_FM2,' - R',num2str(i_notes),'FM2:',cell2mat(table2cell(eventTable(fFM_rev_index,index_notes(i_notes)))));
              end
              QC_notes_FM = strcat(QC_notes_FM1,' - ',QC_notes_FM2);
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1FM1:',eventTable.notes_CS(fFM_fwd_index),' - R2FM1:', eventTable.notes_TT(fFM_fwd_index),' - R1FM2:',eventTable.notes_CS(fFM_rev_index),' - R2FM2:', eventTable.notes_TT(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-3} = num2cell(eventTable.ABCD_Compliant{fFM_fwd_index} && eventTable.ABCD_Compliant{fFM_rev_index});
              redcapTable{index_redcap,c_table_run-2} = strcat(eventTable.AdditionalInfo(fFM_fwd_index),eventTable.AdditionalInfo(fFM_rev_index));
              redcapTable{index_redcap,c_table_run-1} =  num2cell(eventTable.Completed{fFM_fwd_index} && eventTable.Completed{fFM_rev_index});
          end   
      end

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run);
      %if flag_fm
      %    fm_issues_tmp_1 = strcat(eventTable.notes_CS(fFM_fwd_index),eventTable.notes_CS(fFM_rev_index));
      %    fm_issues_tmp_2 = strcat(eventTable.notes_TT(fFM_fwd_index),eventTable.notes_TT(fFM_rev_index));   
      %end 
      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end      

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end


      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end
      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end   
      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0;
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end        
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));               
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end 


  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable]  = classify_rsfMRI_GE(eventTable, redcapTable, rsfMRI_index, index_redcap)
  i_series=0;
  i_pc = 0;
  i_series_noFM = 0;
  flag_fm=0;
  totalvalidtime = 0;
  totaltime_noFM = 0;
  totaltime = 0;
  pc_issue = 0;
  qc_issue = 0;
  incomplete_issue = 0;
  fm_issue = 0;
  fm_pc_issue = 0;
  fm_qc_issue = 0;
  fm_incomplete_issue =0;
  %Allow for 12 maximum runs
  if length(rsfMRI_index)>12
     rsfMRI_index = rsfMRI_index(1:12); 
  end
  i_total = length(rsfMRI_index);
  c_table = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_1_qc_score'));

  % fill issues with 0
  index_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_dis_QC'));
  index_fm_issues = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_FM_dis_QC'));
  index_fm_issues_last = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_FM_Other_QC'));
  for i=index_issues:index_fm_issues_last
      redcapTable{index_redcap,i} = num2cell(0);
  end

  % seek notes indexes
  index_notes = cellfun(@isempty, regexp(eventTable.Properties.VariableNames,'^notes_'));
  index_notes = find(index_notes==0);
  notes_count = length(index_notes);

  qc_flag = ones(1,length(rsfMRI_index));
  for index=1:length(rsfMRI_index)
      c_table_run=(c_table + (index-1)*11);
      index_run = rsfMRI_index(index);  
      % Fill table with rsfMRI_run
      redcapTable{index_redcap,c_table_run} = eventTable.QC(index_run);

      QC_notes = strcat('R1:',cell2mat(table2cell(eventTable(index_run,index_notes(1)))));
      for i_notes=2:notes_count
          QC_notes = strcat(QC_notes,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(index_run,index_notes(i_notes)))));
      end
      redcapTable{index_redcap,c_table_run+1} = {QC_notes};
      %redcapTable{index_redcap,c_table_run+1} = strcat('R1:',eventTable.notes_CS(index_run),' - R2:', eventTable.notes_TT(index_run));
      redcapTable{index_redcap,c_table_run+2} = eventTable.ABCD_Compliant(index_run);
      redcapTable{index_redcap,c_table_run+3} = eventTable.AdditionalInfo(index_run);
      redcapTable{index_redcap,c_table_run+4} = eventTable.Completed(index_run);
      redcapTable{index_redcap,c_table_run+5} = eventTable.subthresh_02_nody(index_run);

      % Check for possible field map, only inmediately before or two runs
      % earlier considered    
      fFM_index = index_run-1; if fFM_index<1, fFM_index=1; end 
      if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
          flag_fm=1;
          redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);

          QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
          for i_notes=2:notes_count
              QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
          end
          redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
          %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
          redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
          redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
          redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
      else
          fFM_index = index_run-2; if fFM_index<1, fFM_index=1; end
          if (strcmp(eventTable(fFM_index,:).SeriesType, 'fMRI_FM'))
              flag_fm=1;
              redcapTable{index_redcap,c_table_run-5} = eventTable.QC(fFM_index);

              QC_notes_FM = strcat('R1:',cell2mat(table2cell(eventTable(fFM_index,index_notes(1)))));
              for i_notes=2:notes_count
                  QC_notes_FM = strcat(QC_notes_FM,' - R',num2str(i_notes),':',cell2mat(table2cell(eventTable(fFM_index,index_notes(i_notes)))));
              end
              redcapTable{index_redcap,c_table_run-4} = {QC_notes_FM};
              %redcapTable{index_redcap,c_table_run-4} = strcat('R1:',eventTable.notes_CS(fFM_index),' - R2:', eventTable.notes_TT(fFM_index));
              redcapTable{index_redcap,c_table_run-3} = eventTable.ABCD_Compliant(fFM_index);
              redcapTable{index_redcap,c_table_run-2} = eventTable.AdditionalInfo(fFM_index);
              redcapTable{index_redcap,c_table_run-1} = eventTable.Completed(fFM_index);
          end   
      end

      %issues_tmp_1 = eventTable.notes_CS(index_run);
      %issues_tmp_2 = eventTable.notes_TT(index_run); 
      %if flag_fm
      %    fm_issues_tmp_1 = eventTable.notes_CS(fFM_index);
      %    fm_issues_tmp_2 = eventTable.notes_TT(fFM_index);   
      %end

      if ~(cell2mat(eventTable.ABCD_Compliant(index_run)))
         pc_issue = pc_issue + 1;
      end
      if ~(cell2mat(eventTable.QC(index_run)))
         qc_issue = qc_issue + 1; 
      end
      if ~(cell2mat(eventTable.Completed(index_run)))
         incomplete_issue = incomplete_issue + 1;
      end

      if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
          totaltime = totaltime + cell2mat(eventTable.subthresh_02_nody(index_run));
      end

      qc = cell2mat(redcapTable{index_redcap,c_table_run});
      if isempty(qc)
          qc=0; 
          qc_flag(index)=0;
      end
      pc = cell2mat(redcapTable{index_redcap,c_table_run+2});
      if isempty(pc), pc=0; end
      ic = cell2mat(redcapTable{index_redcap,c_table_run+4});
      if isempty(ic), ic=0; end

      if qc==0
          [redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, QC_notes);
          %[redcapTable] = issueHandling(redcapTable, index_redcap, index_issues, cell2mat(issues_tmp_1), cell2mat(issues_tmp_2));    
      end 

      if flag_fm
          qc_fm = cell2mat(redcapTable{index_redcap,c_table_run-5});
          if isempty(qc_fm) 
              qc_fm=0; 
              qc_flag(index)=0;
          end
          if ~qc_fm, fm_qc_issue = fm_qc_issue + 1; end
          pc_fm = cell2mat(redcapTable{index_redcap,c_table_run-3});
          if isempty(pc_fm), pc_fm=0; end
          if ~pc_fm, fm_pc_issue = fm_pc_issue + 1; end
          ic_fm = cell2mat(redcapTable{index_redcap,c_table_run-1});
          if isempty(ic_fm), ic_fm=0; end
          if ~ic_fm, fm_incomplete_issue = fm_incomplete_issue + 1; end
          if pc && pc_fm
              i_pc = i_pc+1;
          end
          if pc && qc && ic && pc_fm && qc_fm && ic_fm
              i_series = i_series+1;
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totalvalidtime = totalvalidtime + cell2mat(eventTable.subthresh_02_nody(index_run));
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          elseif qc_fm==0
              [redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, QC_notes_FM);
              %[redcapTable] = issueHandling(redcapTable, index_redcap, index_fm_issues, cell2mat(fm_issues_tmp_1), cell2mat(fm_issues_tmp_2));
          end
          flag_fm=0;
      else
          fm_issue = fm_issue + 1;
          if pc && qc && ic
              i_series_noFM = i_series_noFM+1;
              if ~isempty(cell2mat(eventTable.subthresh_02_nody(index_run)))
                  totaltime_noFM = totaltime_noFM + cell2mat(eventTable.subthresh_02_nody(index_run));
              end
          end
      end
  end   

  c_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_total_ser'));
  redcapTable{index_redcap,c_series} = num2cell(i_total);
  redcapTable{index_redcap,c_series+1} = num2cell(totaltime);
  redcapTable{index_redcap,c_series+2} = num2cell(i_pc);
  if (any(qc_flag)), 
      redcapTable{index_redcap,c_series+3} = num2cell(i_series_noFM);
      redcapTable{index_redcap,c_series+4} = num2cell(totaltime_noFM);
      redcapTable{index_redcap,c_series+5} = num2cell(i_series);
      redcapTable{index_redcap,c_series+6} = num2cell(totalvalidtime);
      redcapTable{index_redcap,c_series+7} = num2cell(i_total-i_series);
  end

  b_series = find(strcmpi(redcapTable.Properties.VariableNames,'iqc_rsfMRI_ser_PC_issues'));
  redcapTable{index_redcap,b_series} = num2cell(pc_issue);
  redcapTable{index_redcap,b_series+1} = num2cell(qc_issue);
  redcapTable{index_redcap,b_series+2} = num2cell(incomplete_issue);
  redcapTable{index_redcap,b_series+3} = num2cell(fm_issue);
  redcapTable{index_redcap,b_series+4} = num2cell(fm_pc_issue);
  redcapTable{index_redcap,b_series+5} = num2cell(fm_qc_issue);
  redcapTable{index_redcap,b_series+6} = num2cell(fm_incomplete_issue);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [redcapTable] = issueHandling(redcapTable, index_redcap, index_table, QC_notes)
 
  if regexp(QC_notes,'dis'), redcapTable{index_redcap,index_table} = num2cell(cell2mat(redcapTable{index_redcap,index_table}) + 1); end

  if regexp(QC_notes,'dco'), redcapTable{index_redcap,index_table+1} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+1}) + 1); end

  if regexp(QC_notes,'fa'), redcapTable{index_redcap,index_table+2} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+2}) + 1); end

  if regexp(QC_notes,'gh'), redcapTable{index_redcap,index_table+3} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+3}) + 1); end

  if regexp(QC_notes,'ht'), redcapTable{index_redcap,index_table+4} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+4}) + 1); end

  if regexp(QC_notes,'hb'), redcapTable{index_redcap,index_table+5} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+5}) + 1); end

  if regexp(QC_notes,'mo'), redcapTable{index_redcap,index_table+6} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+6}) + 1); end

  if regexp(QC_notes,'rf'), redcapTable{index_redcap,index_table+7} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+7}) + 1); end

  if regexp(QC_notes,'sd'), redcapTable{index_redcap,index_table+8} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+8}) + 1); end

  if regexp(QC_notes,'si'), redcapTable{index_redcap,index_table+9} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+9}) + 1); end

  if regexp(QC_notes,'sus'), redcapTable{index_redcap,index_table+10} = num2cell(cell2mat(redcapTable{index_redcap,index_table+10}) + 1); end

  if regexp(QC_notes,'vco'), redcapTable{index_redcap,index_table+11} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+11}) + 1); end

  if regexp(QC_notes,'wr'), redcapTable{index_redcap,index_table+12} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+12}) + 1); end

  if regexp(QC_notes,'SLICE'), redcapTable{index_redcap,index_table+13} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+13}) + 1); end

  if regexp(QC_notes,'line'), redcapTable{index_redcap,index_table+14} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+14}) + 1); end

  if regexp(QC_notes,'moire'), redcapTable{index_redcap,index_table+15} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+15}) + 1); end

  if regexp(QC_notes,'RECON'), redcapTable{index_redcap,index_table+16} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+16}) + 1); end

  if regexp(QC_notes,'MISS'), redcapTable{index_redcap,index_table+17} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+17}) + 1); end

  if regexp(QC_notes,'OTHER'), redcapTable{index_redcap,index_table+18} =  num2cell(cell2mat(redcapTable{index_redcap,index_table+18}) + 1); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function series_info = set_pguidevents(series_info)
  if isfield(series_info,'pGUID')
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.pGUID},{series_info.EventName},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  else
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.id_redcap},{series_info.redcap_event_name},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  end
return;
