clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/MST/';

%%
sbj_nme_hld = mmil_readtext([prj_dat_hld '/' 'clerical' '/' 'subjects']);

for sbj_num = 1%:numel(sbj_nme_hld);
        
    cfg = [];
    
    cfg.prj_dat_hld = prj_dat_hld;
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    
    %% load the data
    cfg.tsk_hld = mmil_load_subj_info([prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    cfg.out_pth = [prj_dat_hld '/' 'epoch_data'];
    ejk_load_mst(cfg)
    
    %% Create Frequency Data
%     cfg.out_pth = [sbj_dat_hld{sbj_num} '/' 'freq_data'];
%     cfg.tsk = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
%     Pedot_frequency_load(cfg)

    %% Stats
    ejk_stt_mst(cfg);
    
end