% Stats
cfg = [];
cfg.eve   = [211 212];
cfg.add_stt  = 'new_old_stt';
cfg.alt_eve  = 'new_old';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);