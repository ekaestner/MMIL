% Stats
cfg = [];
cfg.events   = [261 262];
cfg.add_stt  = 'vis_Orth_med_stt';
cfg.alt_eve  = 'vis_Orth_med';
cfg.alpha    = .05;
stt_dat      = ft_fdrstats(cfg,stt_dat);
