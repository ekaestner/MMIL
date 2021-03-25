clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

cfg = [];
for sbj_num = 22:numel(sbj_nme_hld); % 1:16 17 18:23 24:37
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% load the data
    cfg.tsk_hld = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    cfg.out_pth = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    iSASZ_initial_load2(cfg)
    
end

cfg = [];
for sbj_num = 2:numel(sbj_nme_hld);
    
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.fle_out_pth = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    cfg.ovr_wrt = 1;
    
    %% Put together initial events
    iSASZ_initial_events(cfg)
    
    %% run the initial stats
    iSASZ_initial_stats2(cfg)
    
    %% Run the initial plots
    iSASZ_initial_plot_v3(cfg)
    
    %% Catalog Effects
%     cfg.fle_loc = sbj_dat_hld{sbj_num};
%     cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/';
%     cfg.sbj_clr_hld = sbj_clr_hld{sbj_num};
%     iSASZ_effects(cfg)
    
    %% Timing Plot
    
    
    %% Location Plot
    
    
    %% Table of Effects
    
end