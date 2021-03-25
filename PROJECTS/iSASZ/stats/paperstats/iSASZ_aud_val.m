% Stats
cfg = [];
cfg.eve   = [261 262];
cfg.add_stt  = 'val_stt';
cfg.alt_eve  = 'val';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);