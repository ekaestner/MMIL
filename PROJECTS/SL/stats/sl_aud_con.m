cfg.eve   = unique(dat.cfg.alt_eve.aud_con(~isnan(dat.cfg.alt_eve.aud_con)));
cfg.add_stt  = 'aud_con_stt';
cfg.alt_eve  = 'aud_con';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);