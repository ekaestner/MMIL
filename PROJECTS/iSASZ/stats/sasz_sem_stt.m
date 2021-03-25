has_vis = any(stt_dat.trialinfo < 9);
has_aud = any(stt_dat.trialinfo > 9);

% Stats
if has_aud
    cfg = [];
    cfg.events   = [123 124];
    cfg.add_stt  = 'aud_dif_sem_stt';
    cfg.alt_eve  = 'sem_all_eve';
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end
if has_vis
    cfg = [];
    cfg.events   = [121 122];
    cfg.add_stt  = 'vis_dif_sem_stt';
    cfg.alt_eve  = 'sem_all_eve';
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end
