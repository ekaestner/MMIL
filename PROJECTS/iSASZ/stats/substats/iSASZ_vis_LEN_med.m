% Stats
cfg = [];
cfg.eve   = [131 132];
cfg.add_stt  = 'LEN_med_stt';
cfg.alt_eve  = 'LEN_med';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);
