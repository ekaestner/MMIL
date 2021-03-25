cfg.eve   = unique(dat.cfg.alt_eve.vis_con(~isnan(dat.cfg.alt_eve.vis_con)));
cfg.add_stt  = 'vis_con_stt';
cfg.alt_eve  = 'vis_con';
cfg.alp    = .05;
cfg.cor_mth  = 'mnt_car';
stt_dat      = mmil_ieeg_statistics5(cfg,stt_dat);