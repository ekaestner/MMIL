% Stats
cfg = [];
cfg.eve   = [271 272];
cfg.add_stt  = 'ars_stt';
cfg.alt_eve  = 'ars';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);