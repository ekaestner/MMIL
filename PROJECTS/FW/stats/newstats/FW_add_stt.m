function FW_add_stt(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Stats
% Make
cfg = [];
cfg.stt_fnc  = {'fw_vis_wrd_ffn'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

%% Mask
cfg     = [];
cfg.stt     = { 'vis_wrd_ffn' 'vis_wrd_ffn' }; %   
cfg.stt_msk = { 'vis_stm_01'  'vis_stm_01'  }; %
cfg.pst_fix = '_msk_01';
bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
toc
ft_func([],cfg,bcc_dat);

end