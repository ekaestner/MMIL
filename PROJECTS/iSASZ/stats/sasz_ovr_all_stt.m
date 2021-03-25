has_vis = any(stt_dat.trialinfo < 9);
has_aud = any(stt_dat.trialinfo > 9);

% Stats
if has_aud
    cfg = [];
    cfg.events   = [102];
    cfg.add_stt  = 'aud_ovr_all_stt';
    cfg.alt_eve  = 'ovr_all_evt';
    cfg.alpha    = .05;
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end
if has_vis
    cfg = [];
    cfg.events   = [101];
    cfg.add_stt  = 'vis_ovr_all_stt';
    cfg.alt_eve  = 'ovr_all_evt';
    cfg.alpha    = .05;
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end