function bgr_dat = SL_pedot_initial_stats(bgr_dat)
%% Run Stats
tic
cfg = [];
cfg.stt_fnc  = {'sll_ovr_all' 'sll_vis_nse' 'sll_aud_nse' 'sll_mtc_mss'};
cfg.loc      = 'local';
cfg.fld_nme  = bgr_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bgr_dat      = ft_func(@mmil_cloud_stat,cfg,bgr_dat);
toc

%% Mask stats based on All 4 stats
cfg     = [];
cfg.stt     = {'vis_nse_stt' 'aud_nse_stt' 'vis_mtc_stt'};
cfg.stt_msk = {'vis_ovr_stt' 'vis_ovr_stt' 'vis_ovr_stt'};
bgr_dat = ft_func(@ft_mask_stats,cfg,bgr_dat);

end