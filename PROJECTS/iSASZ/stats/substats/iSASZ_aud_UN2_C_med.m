% Stats
cfg = [];
cfg.eve   = [241 242];
cfg.add_stt  = 'UN2_C_med_stt';
cfg.alt_eve  = 'UN2_C_med';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);
