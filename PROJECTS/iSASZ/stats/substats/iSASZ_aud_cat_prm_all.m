% Stats
cfg = [];
cfg.events   = [2001 2002];
cfg.add_stt  = 'aud_cat_prm_all_stt';
cfg.alt_eve  = 'aud_cat_prm_all';
cfg.alpha    = .05;
stt_dat      = ft_fdrstats(cfg,stt_dat);
