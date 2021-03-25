clear; clc;

cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/';

sbj_nme_hld = mmil_readtext([cfg.clr_fld 'subjects']);

%% Run Subjects
for sbj_num = 1:numel(sbj_nme_hld);
    
    %% Setting up Variables
    cfg.sbj_nme = sbj_nme_hld{sbj_num};
    cfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
       
    %% load the data
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
    
    SL_initial_load(cfg)
    
    % Make ECOG/DEPTH Split
    cfg.alt_lab = 'label';
    mmil_create_depth(cfg)
    
    %% run the initial stats
    cfg.out_pth =  '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
    
    sl_initial_stats_2(cfg)
    
    %% run the initial plots
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
    cfg.in_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
    cfg.ovr_wrt = 2;
    
    SL_initial_plot_2(cfg)
    
    % Split depths/ecog
    cfg.plt_spl = [cfg.clr_out_pth '/' cfg.sbj_nme '/' 'ovr'];
    cfg.idn_spl = 1;
    mmil_split_plot(cfg);
    mmil_split_stat(cfg)
    
end

%% Check Stats
