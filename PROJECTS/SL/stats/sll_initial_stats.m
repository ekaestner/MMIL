clear; clc;

chn_chs = 1;
mve_plt = 1;

% Setting up Variables
sbj_hld  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');

for sbj_num = 3 %[1:17 19:numel(sbj_hld)];
   
    sbj  = sbj_hld{sbj_num};
    
    fprintf('Starting work on %s \n',sbj)
    
    infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/' sbj '_overall_data.mat'];
    outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/';
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [outpath '/' sbj '_overall_data.mat'];
    sll_dat  = ft_func([],cfg);
    
    %% Effect Stats
    tic
    cfg = [];
    cfg.stt_fnc  = {'sll_vis_nse' 'sll_aud_nse' 'sll_mtc_mss'};
    cfg.loc      = 'local';
    cfg.fld_nme  = sll_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    sll_dat      = ft_func(@mmil_cloud_stat,cfg,sll_dat);
    toc
    
    %% Base Stats
    tic
    cfg = [];
    cfg.stt_fnc  = {'sll_ovr_all'};
    cfg.loc      = 'local';
    cfg.fld_nme  = sll_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    sll_dat      = ft_func(@mmil_cloud_stat,cfg,sll_dat);
    toc
    
    %% Mask stats based on All 4 stats
    cfg     = [];
    cfg.stt     = {'vis_nse_stt' 'aud_nse_stt' 'vis_mtc_stt'};
    cfg.stt_msk = {'vis_ovr_stt' 'vis_ovr_stt' 'vis_ovr_stt'};
    sll_dat = ft_func(@ft_mask_stats,cfg,sll_dat);
    
    %% Differences in Phonemes
    tic
    cfg = [];
    cfg.stt_fnc = {'sl_anova_2x2'};
    cfg.loc     = 'local';
    cfg.fld_nme  = sll_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    sll_dat = ft_func(@mmil_cloud_stat,cfg,sll_dat);
    toc
    
    %% Save Data & Stats
    cfg = [];
    cfg.str_nme  = 'sll_dat';
    cfg.save     = 'yes';
    cfg.sve_app  = 'app_all';
    cfg.filename =[outpath '/' sbj '_overall_data.mat'];
    ft_func([],cfg,sll_dat);  
    
end