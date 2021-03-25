clear; clc;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');

pot_sbj_num = [1:17 19:numel(subj)];

parfor ind_sbj_num = 1:numel(pot_sbj_num);
    
    sbj_num = pot_sbj_num(ind_sbj_num);
    sbj  = subj{sbj_num};
    
    fprintf('Starting work on %s \n',sbj)
    
    infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/' sbj '_overall_data.mat'];
    outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
    
    %% Load
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [outpath '/' sbj '_overall_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    %% Mask stats based on overall stats
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    % Stat Mask
    if has_aud
        cfg     = [];
        cfg.stt     = {'aud_dif_rep_stt' 'aud_dif_sem_stt'};
        cfg.stt_msk = {'aud_ovr_all_stt' 'aud_ovr_all_stt'};
        sem_dat = ft_func(@ft_mask_stats,cfg,sem_dat);
    end
    if has_vis
        cfg     = [];
        cfg.stt     = {'vis_dif_rep_stt' 'vis_dif_sem_stt'};
        cfg.stt_msk = {'vis_ovr_all_stt' 'vis_ovr_all_stt'};
        sem_dat = ft_func(@ft_mask_stats,cfg,sem_dat);
    end
    
    %% Save Data & Stats
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.sve_app  = 'app_all';
    cfg.filename =[outpath '/' sbj '_overall_data.mat'];
    ft_func([],cfg,sem_dat);
    
end