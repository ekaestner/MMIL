% Stats
cfg.eve   = [102];
cfg.add_stt  = 'vis_new_ovr_stt_avg';
cfg.alt_eve  = 'vis_new_bse';
cfg.alp    = .01
cfg.cor_mth  = 'fdr';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);