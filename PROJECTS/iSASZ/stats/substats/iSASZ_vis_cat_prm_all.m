% Stats
cfg = [];
cfg.events   = [1001 1002];
cfg.add_stt  = 'vis_cat_prm_all_stt';
cfg.alt_eve  = 'vis_cat_prm_all';
cfg.alpha    = .05;
stt_dat      = ft_fdrstats(cfg,stt_dat);
