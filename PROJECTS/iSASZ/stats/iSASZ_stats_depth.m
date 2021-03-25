function iSASZ_stats_depth(fcfg)

%% Load
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data_depth.mat'];
bcc_dat  = ft_func([],cfg);

%% Stats
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
    
    % Overall Stats
    cfg = [];
    cfg.stt_fnc  = { 'sasz_vis_ovr_new_dpt' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Repetition Stats
    cfg = [];
    cfg.stt_fnc  = { 'sasz_vis_old_new' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantics
    cfg = [];
    cfg.stt_fnc  = { 'sasz_vis_ani_obj' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)
    
    % Overall Stats
    cfg = [];
    cfg.stt_fnc  = { 'sasz_aud_ovr_new_dpt' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Repetition Stats
    cfg = [];
    cfg.stt_fnc  = { 'sasz_aud_old_new' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantics
    cfg = [];
    cfg.stt_fnc  = { 'sasz_aud_ani_obj' };
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

%% Save Data
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data_depth.mat'];
ft_func([],cfg,bcc_dat);

end