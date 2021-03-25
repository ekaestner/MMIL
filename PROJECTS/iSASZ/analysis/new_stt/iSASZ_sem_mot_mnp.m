% Stats
cfg = [];
cfg.eve   = [193 194];
cfg.add_stt  = 'sem_mot_stt';
cfg.alt_eve  = 'sem_mot';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);