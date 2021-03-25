% Stats
cfg.eve   = [111];
cfg.add_stt  = 'vis_new_ovr_stt';
cfg.alt_eve  = 'vis_new_old';
cfg.alp    = .05;
cfg.cor_mth  = 'fdr';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);