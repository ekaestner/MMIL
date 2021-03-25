% Stats
cfg = [];
cfg.eve   = [151 152];
cfg.add_stt  = 'FREQ_med_stt';
cfg.alt_eve  = 'FREQ_med';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics3(cfg,stt_dat);
