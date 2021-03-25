has_vis = any(stt_dat.trialinfo < 9);
has_aud = any(stt_dat.trialinfo > 9);

% Stats
if has_aud
    cfg = [];
    cfg.events   = [113 114];
    cfg.add_stt  = 'aud_dif_rep_stt';
    cfg.alt_eve  = 'rep_all_eve';
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end
if has_vis
    cfg = [];
    cfg.events   = [111 112];
    cfg.add_stt  = 'vis_dif_rep_stt';
    cfg.alt_eve  = 'rep_all_eve';
    stt_dat      = ft_fdrstats(cfg,stt_dat);
end
