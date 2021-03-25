% Stats
cfg = [];
cfg.eve   = [291 192];
cfg.add_stt  = 'sem_mnp_stt';
cfg.alt_eve  = 'sem_mnp';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);