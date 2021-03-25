function iSASZ_initial_stats3(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile  = [fcfg.fle_out_pth '/' 'epoch_data' '/' sbj '_overall_data.mat'];
outpath = [fcfg.fle_out_pth '/' 'epoch_data' '/'];

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Run Stats
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)

    % Overall
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ovr'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Old/New
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_old_new'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ani_obj'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic Priming
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_sem_prm'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Length
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_lng'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Bigram
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_big'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Frequency
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_frq'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Neighborhood
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)

    % Overall
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ovr'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Old/New
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_old_new'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ani_obj'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic Priming
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_sem_prm'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Length
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_lng'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Bigram
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_big'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Frequency
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_frq'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Neighborhood
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

%% Save Data
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

