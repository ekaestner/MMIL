% Stats
cfg = [];
cfg.eve  = [221 222];
cfg.add_stt  = 'ani_obj_stt';
cfg.alt_eve  = 'ani_obj';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);