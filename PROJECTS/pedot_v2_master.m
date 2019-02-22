clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/Pedot2/';
clr_fld     = [prj_dat_hld '/' 'clerical' '/'];

sbj_nme_hld = mmil_readtext([clr_fld '/' 'subjects']);

for sbj_num = 1%:numel(sbj_nme_hld);
        
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    
    %% load the data
    cfg.tsk_hld = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    cfg.out_pth = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    mmil_pedot_load_v2(cfg)
    
    %% Create Frequency Data
%     cfg.out_pth = [sbj_dat_hld{sbj_num} '/' 'freq_data'];
%     cfg.tsk = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
%     Pedot_frequency_load(cfg)
    
end