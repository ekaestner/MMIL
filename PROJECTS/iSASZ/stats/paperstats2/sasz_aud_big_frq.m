% Stats
cfg.eve   = [241 242];
cfg.add_stt  = 'aud_big_frq_stt';
cfg.alt_eve  = 'aud_big_med_crr_by_aud_frq';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);