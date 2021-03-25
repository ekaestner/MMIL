% Stats
cfg = [];
cfg.eve   = [1011 1012];
cfg.add_stt  = 'cat_prm_new_stt';
cfg.alt_eve  = 'cat_prm_new';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);
