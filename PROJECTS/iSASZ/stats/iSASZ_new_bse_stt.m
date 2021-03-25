function iSASZ_new_bse_stt(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Base Stats
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
    % Novel Setup
    cfg = [];
    cfg.eve_rul     = { 'new'         };
    cfg.old_events  = { { [1 2 3 4] } };
    cfg.new_events  = { [ 102 ]       };
    cfg.crt_alt_eve = { 'vis_new_bse' };
    cfg.use_alt_eve = { 'trialinfo'   };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ovr_new'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)
    % Novel Setup
    cfg = [];
    cfg.eve_rul     = { 'new'             };
    cfg.old_events  = { { [11 12 13 14] } };
    cfg.new_events  = { [ 202 ]           };
    cfg.crt_alt_eve = { 'aud_new_bse'     };
    cfg.use_alt_eve = { 'trialinfo'       };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    % Novel Baseline
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ovr_new'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

%% Save Data
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end