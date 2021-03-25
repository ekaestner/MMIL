% Stats
cfg.eve   = [251 252];
cfg.add_stt  = 'aud_frq_lng_stt';
cfg.alt_eve  = 'aud_frq_med_crr_by_aud_lng';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);